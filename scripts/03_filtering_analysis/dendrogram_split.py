"""
dendrogram_split.py
===================
MHC-dendrogram-stratified train/val/test splits with peptide-level
leakage control. Supports four experiment modes:
  hla_a, hla_b, hla_c  — leaf-cluster-based HLA family splits
  cross_species         — train on HLA, test/val on non-human species

Steps
-----
1.  Load combined parquet.
2.  Encode MHC pseudosequences (physicochemical, 49x14=686 features).
3.  Ward hierarchical clustering; cut at --leaf_distance_threshold to
    obtain small leaf-level clusters of 2-5 alleles.
    Save dendrogram + cut line to mhc_clustering/dendrogram.png.
    Save allele→cluster assignments to mhc_clustering/allele_cluster_assignments.csv.
4.  Dispatch to mode-specific pipeline (allele assignments fixed once;
    only peptide-cluster exclusions change per threshold).
5.  For each peptide-cluster threshold file:
      a. Assign peptide cluster IDs (precomputed + singleton fallback).
      b. Run set-cover for test alleles → excluded_test_clusters.
      c. Per fold: run set-cover for val alleles → save train/val parquets.
      d. Save binder/nonbinder counts, peptide cluster summaries,
         fold_summary.csv, split_meta.json.

Strict train mask (always):
    train = ~in_test_cluster & ~in_val_cluster & ~in_test_allele & ~in_val_allele

Usage
-----
python dendrogram_split.py \\
    --input_dir   .../combined.parquet \\
    --pseudo_csv  .../mhc_pseudo_class1_n.csv \\
    --pep_cluster_files .../pep50.tsv .../pep60.tsv .../pep70.tsv .../pep80.tsv \\
    --thresholds 50 60 70 80 \\
    --output_dir  .../output_dendro_split \\
    --mode hla_a \\
    --leaf_distance_threshold 25.0 \\
    --n_folds 5
"""

import argparse
import json
import logging
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

warnings.filterwarnings("ignore")

# ── column names ───────────────────────────────────────────────────────────────
COL_PEPTIDE = "long_mer"
COL_ALLELE  = "allele"

# ── peptide cluster set cover hyper-parameters ─────────────────────────────────
PEPTIDE_COVER_COST_FN         = "log"
PEPTIDE_COVER_MIN_COV         = 0.90
PEPTIDE_COVER_MIN_ROWS        = 100
PEPTIDE_COVER_PER_ALLELE_FRAC = 0.20
PEPTIDE_COVER_PER_ALLELE_CAP  = 1000

# ── known non-HLA species prefixes (cleaned-allele format) ────────────────────
KNOWN_NONHLA_PREFIXES = ["H2", "BOLA", "DLA", "EQCA", "GOGO", "MAMU", "PATR", "SLA"]

log = logging.getLogger(__name__)


# ── clean_key (verbatim from filter_map_resample.py) ─────────────────────────

def clean_key(allele_key: str) -> str:
    """Normalise allele name: 'HLA-A*02:01' -> 'HLAA0201'"""
    if allele_key is None:
        return "None"
    mapping = str.maketrans({"*": "", ":": "", " ": "", "/": "_", "-": ""})
    return allele_key.translate(mapping).upper()


# ── physicochemical properties (verbatim from pepclust_asgary_clust.py) ────────

Physicochemical_Properties = {
    "A": [-0.7115, -1.3796, -0.3237,  1.7644, -0.0249, -1.7046, -1.5572, -0.3873, -0.7860, -0.7889, -0.9536, -0.1956,  0.0770, -0.1703],
    "C": [-1.4936,  0.5087, -0.3881, -1.1136, -0.0977, -0.5661, -0.8431, -0.3873, -0.7860, -0.7889, -1.0949, -1.1785, -0.4755, -0.9901],
    "D": [ 0.9952,  0.3345, -1.7827, -0.2679, -2.0926, -0.1392, -0.6552, -0.3873, -0.7860,  2.8401,  1.0251,  1.6568, -1.8416,  0.8496],
    "E": [ 1.3145, -1.5410,  0.7795,  0.1454, -0.4864,  0.3589,  0.0212, -0.3873, -0.7860,  1.6305,  0.8452,  1.3922, -1.5743,  0.9259],
    "F": [-1.1431, -0.6188,  0.9861, -0.4213,  0.3478,  0.9993,  0.9983,  2.5820, -0.7860, -0.7889, -1.5060, -1.2919, -0.2319, -1.6955],
    "G": [-0.4962,  1.7772,  0.7061,  1.1811,  1.4511, -2.2027, -2.2712, -0.3873, -0.7860, -0.7889,  1.4490,  0.1446,  0.0591,  0.4207],
    "H": [ 0.2526, -0.4339, -0.7929, -1.6181,  0.0205,  0.6435,  0.3595, -0.3873,  0.5531,  0.4208, -1.3904,  0.6739,  1.0214,  0.0394],
    "I": [-1.3855, -0.5728,  1.1059,  0.4566,  0.6176, -0.2104,  0.5850, -0.3873, -0.7860, -0.7889,  0.3955, -1.2919,  0.0888, -1.3905],
    "K": [ 1.8074, -0.5878,  0.3083, -0.2879,  1.1733,  0.3233,  0.9983, -0.3873,  1.8922, -0.7889,  0.4469,  1.0141,  2.2985,  1.1546],
    "L": [-1.1567, -1.0430, -0.7090,  1.4266, -0.5365, -0.2104,  0.5850, -0.3873, -0.7860, -0.7889, -0.8508, -1.4053,  0.0651, -1.2951],
    "M": [-0.7864, -1.6169,  1.1498, -1.0969,  0.8821,  0.4300,  0.5850, -0.3873, -0.7860, -0.7889, -1.7501, -1.1029, -0.0775, -0.9139],
    "N": [ 0.8860,  0.8966,  0.6906, -0.1679,  0.6957, -0.1748, -0.4673, -0.3873,  0.5531,  0.4208,  0.4083,  1.1276, -0.2735,  1.2309],
    "P": [ 0.0997,  2.2356, -0.7704,  0.4877, -0.8570, -0.7796, -0.6928, -0.3873, -0.7860, -0.7889,  0.9994, -0.2334,  0.2552,  0.2205],
    "Q": [ 0.8714, -0.1795, -1.4577, -0.5391, -1.1649,  0.3233,  0.2091, -0.3873,  0.5531,  0.4208,  0.7938,  0.7117, -0.1309,  1.1546],
    "R": [ 1.5027, -0.0470,  0.7919,  0.5088,  2.0074,  1.3194,  1.4869, -0.3873,  3.2313, -0.7889,  1.2563,  0.7117,  2.9044,  1.4692],
    "S": [-0.3340,  1.5068, -2.3337,  0.7644, -1.6952, -1.1354, -1.3317, -0.3873,  0.5531,  0.4208,  0.9737,  0.2203, -0.1131,  0.3063],
    "T": [-0.1301,  0.3601,  1.1468,  1.0288,  0.9495, -0.6373, -0.5801, -0.3873,  0.5531,  0.4208,  0.1642, -0.0066, -0.1606,  0.1156],
    "V": [-1.4874, -0.2864, -0.2293,  1.4000, -0.7702, -0.7084, -0.1291, -0.3873, -0.7860, -0.7889, -0.5810, -1.0273,  0.0532, -0.9806],
    "W": [-0.7157,  0.0214,  0.3776, -2.3448, -0.0503,  2.3868,  2.0506,  2.5820,  0.5531, -0.7889, -1.6217, -1.2163,  0.0116, -1.5716],
    "Y": [ 0.1736,  0.8987,  1.5881, -0.9113,  1.0824,  1.5685,  1.2238,  2.5820,  0.5531,  0.4208, -0.1441, -0.9139, -0.1250, -0.9043],
    "B": [ 0.9411,  0.6155, -0.5463, -0.2179, -0.6988, -0.1570, -0.5613, -0.3873, -0.1164,  1.6305,  0.7167,  1.3922, -1.0576,  1.0403],
    "Z": [ 1.0929, -0.8603, -0.3391, -0.1968, -0.8257,  0.3411,  0.1152, -0.3873, -0.1164,  1.0256,  0.8195,  1.0519, -0.8497,  1.0403],
    "X": [-0.0969,  0.0118,  0.0422,  0.0199,  0.0726, -0.0681, -0.1291, -0.3873, -0.1164, -0.1841, -0.4011, -0.2334,  0.0770, -0.0560],
}

_PROP_MATRIX = np.array(
    [Physicochemical_Properties.get(chr(i), Physicochemical_Properties["X"])
     for i in range(128)],
    dtype=np.float32,
)  # (128, 14)
for _i in range(97, 123):
    _PROP_MATRIX[_i] = _PROP_MATRIX[_i - 32]


# ── verbatim from pepclust_asgary_clust.py ─────────────────────────────────────

def encode_full_vectorized(sequences):
    results = []
    for s in sequences:
        indices = np.frombuffer(s.upper().encode("ascii"), dtype=np.uint8)
        results.append(_PROP_MATRIX[indices].ravel())
    return np.vstack(results)


def load_precomputed_peptide_clusters(tsv_path):
    """
    Load external peptide clusters (Asgary TSV) and return
    [COL_PEPTIDE, "cluster"] with integer cluster IDs.
    """
    log.info("Loading precomputed peptide clusters...")
    t0 = time.time()

    if not tsv_path.exists():
        raise FileNotFoundError(
            f"Missing precomputed cluster file: {tsv_path}. "
            "Expected columns: cluster_id, sequence"
        )

    raw = pd.read_csv(
        tsv_path,
        sep="\t",
        usecols=["cluster_id", "sequence"],
        dtype={"cluster_id": "string", "sequence": "string"},
    )
    raw = raw.dropna(subset=["cluster_id", "sequence"]).copy()
    raw["cluster_id"] = raw["cluster_id"].str.strip()
    raw["sequence"]   = raw["sequence"].str.strip()
    raw = raw[raw["sequence"] != ""].copy()

    n_before = len(raw)
    raw = raw.drop_duplicates(subset=["sequence"], keep="first").copy()
    n_dedup = n_before - len(raw)

    parsed = pd.to_numeric(
        raw["cluster_id"].str.extract(r"(-?\d+)$", expand=False),
        errors="coerce",
    )
    if parsed.notna().all():
        cluster_ids = parsed.astype(np.int32)
    else:
        labels = raw["cluster_id"].astype(str).values
        uniq = sorted(set(labels))
        label_to_id = {lbl: i for i, lbl in enumerate(uniq)}
        cluster_ids = pd.Series(
            [label_to_id[v] for v in labels], index=raw.index, dtype=np.int32
        )

    pep_clusters = pd.DataFrame({
        COL_PEPTIDE: raw["sequence"].astype(str),
        "cluster":   cluster_ids.values,
    })

    n_cluster = pep_clusters["cluster"].nunique()
    log.info(
        f"  Loaded {len(pep_clusters):,} unique peptides, {n_cluster:,} clusters"
        f" from {tsv_path}"
        + (f" (dropped {n_dedup:,} duplicate sequences)" if n_dedup else "")
        + f" [{time.time()-t0:.1f}s]"
    )
    return pep_clusters


def peptide_cluster_set_cover(
    target_alleles,
    pep_clust_to_alleles,
    global_cluster_sizes,
    target_cluster_rows=None,
    cluster_allele_rows=None,
    per_allele_min_rows=None,
    excluded_clusters=frozenset(),
    cost_fn=PEPTIDE_COVER_COST_FN,
    min_coverage=PEPTIDE_COVER_MIN_COV,
    min_target_rows=PEPTIDE_COVER_MIN_ROWS,
):
    """
    Greedy weighted set cover: select peptide clusters that collectively
    cover the most target alleles while minimising training-data loss.
    (Verbatim from pepclust_asgary_clust.py)
    """
    target = set(target_alleles)
    target_cluster_rows = target_cluster_rows or {}

    covers: dict = {
        c: (als & target)
        for c, als in pep_clust_to_alleles.items()
        if c not in excluded_clusters and (als & target)
    }

    cluster_sizes = {c: global_cluster_sizes.get(c, 1) for c in covers}
    cluster_target_rows = {
        c: int(target_cluster_rows.get(c, cluster_sizes.get(c, 0)))
        for c in covers
    }

    reachable = set()
    for als in covers.values():
        reachable |= als
    universe = target & reachable

    if len(universe) == 0:
        return set(), set(), {
            "n_candidates": 0, "n_selected": 0, "n_pruned": 0,
            "coverage_frac": 0.0, "target_rows_selected": 0,
            "target_rows_required": int(min_target_rows),
            "total_rows_held_out": 0,
        }

    uncovered  = set(universe)
    selected   = []
    target_cnt = int(np.ceil(min_coverage * len(universe)))

    # ── Phase 1: greedy coverage ───────────────────────────────────────────────
    while (len(universe) - len(uncovered)) < target_cnt and uncovered:
        best_c, best_score, best_size = None, -1.0, float("inf")
        for c, als in covers.items():
            gain = len(als & uncovered)
            if gain == 0:
                continue
            sz    = cluster_sizes.get(c, 1)
            cost  = np.log(1 + sz) if cost_fn == "log" else sz
            score = gain / cost
            if score > best_score or (score == best_score and sz < best_size):
                best_score, best_c, best_size = score, c, sz
        if best_c is None:
            break
        selected.append(best_c)
        uncovered -= covers[best_c]

    selected_set = set(selected)

    # ── Phase 2b: global row backfill ──────────────────────────────────────────
    target_rows_required = max(0, int(min_target_rows))

    def _rows_in_selected(sel_set):
        return int(sum(cluster_target_rows.get(c, 0) for c in sel_set))

    selected_target_rows = _rows_in_selected(selected_set)
    if target_rows_required > 0 and selected_target_rows < target_rows_required:
        remaining = set(covers.keys()) - selected_set
        while remaining and selected_target_rows < target_rows_required:
            best_c, best_score, best_size = None, -1.0, float("inf")
            for c in remaining:
                add_rows = cluster_target_rows.get(c, 0)
                if add_rows <= 0:
                    continue
                sz    = cluster_sizes.get(c, 1)
                cost  = np.log(1 + sz) if cost_fn == "log" else sz
                score = add_rows / cost
                if score > best_score or (score == best_score and sz < best_size):
                    best_score, best_c, best_size = score, c, sz
            if best_c is None:
                break
            selected_set.add(best_c)
            remaining.remove(best_c)
            selected_target_rows += cluster_target_rows.get(best_c, 0)

    # ── Phase 2c: per-allele row backfill ──────────────────────────────────────
    cluster_allele_rows = cluster_allele_rows or {}
    per_allele_min_rows = per_allele_min_rows or {}

    reachable_per_allele: dict = {a: 0 for a in target}
    for c in covers:
        for a, n in cluster_allele_rows.get(c, {}).items():
            if a in reachable_per_allele:
                reachable_per_allele[a] += int(n)

    effective_per_allele = {
        a: min(int(per_allele_min_rows.get(a, 0)), reachable_per_allele.get(a, 0))
        for a in target
    }

    def _per_allele_rows(sel_set):
        rows = {a: 0 for a in target}
        for cc in sel_set:
            for a, n in cluster_allele_rows.get(cc, {}).items():
                if a in rows:
                    rows[a] += int(n)
        return rows

    sel_per_allele = _per_allele_rows(selected_set)
    deficit = {a: max(0, effective_per_allele[a] - sel_per_allele[a]) for a in target}
    if any(deficit.values()) and cluster_allele_rows:
        remaining = set(covers.keys()) - selected_set
        while remaining and any(deficit.values()):
            best_c, best_score, best_size = None, -1.0, float("inf")
            for c in remaining:
                gain = 0
                for a, n in cluster_allele_rows.get(c, {}).items():
                    d = deficit.get(a, 0)
                    if d > 0:
                        gain += min(int(n), d)
                if gain <= 0:
                    continue
                sz    = cluster_sizes.get(c, 1)
                cost  = np.log(1 + sz) if cost_fn == "log" else sz
                score = gain / cost
                if score > best_score or (score == best_score and sz < best_size):
                    best_score, best_c, best_size = score, c, sz
            if best_c is None:
                break
            selected_set.add(best_c)
            remaining.remove(best_c)
            for a, n in cluster_allele_rows.get(best_c, {}).items():
                if a in deficit and deficit[a] > 0:
                    deficit[a] = max(0, deficit[a] - int(n))

    # ── Phase 3: redundancy pruning (largest-first) ────────────────────────────
    n_before = len(selected_set)

    def _union_covered(sel_set):
        cov = set()
        for cc in sel_set:
            cov |= covers[cc]
        return cov

    for c in sorted(selected_set, key=lambda c: cluster_sizes.get(c, 0), reverse=True):
        trial            = selected_set - {c}
        trial_cov        = _union_covered(trial)
        trial_per_allele = _per_allele_rows(trial)
        if (
            len(trial_cov) >= target_cnt
            and _rows_in_selected(trial) >= target_rows_required
            and all(trial_per_allele[a] >= effective_per_allele[a] for a in target)
        ):
            selected_set.remove(c)

    covered              = _union_covered(selected_set)
    selected_target_rows = _rows_in_selected(selected_set)
    held_out             = sum(cluster_sizes.get(c, 0) for c in selected_set)
    final_per_allele     = _per_allele_rows(selected_set)
    alleles_meeting_floor = sum(
        1 for a in target
        if final_per_allele.get(a, 0) >= effective_per_allele.get(a, 0)
    )
    stats = {
        "n_candidates":                    len(covers),
        "n_selected":                      len(selected_set),
        "n_pruned":                        n_before - len(selected_set),
        "coverage_frac":                   len(covered) / len(target) if target else 0.0,
        "target_rows_selected":            selected_target_rows,
        "target_rows_required":            target_rows_required,
        "total_rows_held_out":             held_out,
        "alleles_meeting_per_allele_floor": alleles_meeting_floor,
        "alleles_with_per_allele_floor":   sum(
            1 for a in target if effective_per_allele.get(a, 0) > 0
        ),
    }
    return selected_set, covered, stats


# ── helpers ────────────────────────────────────────────────────────────────────

def save_parquet(df: pd.DataFrame, path: Path, label: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    n_b = int((df["label"] == 1).sum())
    n_n = int((df["label"] == 0).sum())
    log.info(f"  {label:35s}  {len(df):,} rows  |  {n_b:,} binders  {n_n:,} non-binders")


def setup_logging(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fmt = "%(asctime)s | %(levelname)s | %(message)s"
    logging.basicConfig(
        level=logging.INFO,
        format=fmt,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(output_dir / "dendrogram_split.log", mode="w"),
        ],
    )


def load(input_dir: Path, output_dir: Path) -> pd.DataFrame:
    log.info("Loading parquets...")
    df = pd.read_parquet(input_dir)
    log.info(f"  Binder rows    : {len(df[df['label'] == 1]):,}  ({df[df['label'] == 1][COL_ALLELE].nunique():,} alleles)")
    log.info(f"  Nonbinder rows : {len(df[df['label'] == 0]):,}  ({df[df['label'] == 0][COL_ALLELE].nunique():,} alleles)")
    return df


def print_peptide_cluster_coverage(df: pd.DataFrame, pep_clusters: pd.DataFrame) -> None:
    unique_peps  = set(df[COL_PEPTIDE].unique())
    covered_peps = unique_peps & set(pep_clusters[COL_PEPTIDE].unique())
    frac_pep     = len(covered_peps) / len(unique_peps) if unique_peps else 0.0
    rows_covered = int(df[COL_PEPTIDE].isin(covered_peps).sum())
    frac_rows    = rows_covered / len(df) if len(df) else 0.0
    log.info(
        f"Peptide cluster coverage: {len(covered_peps):,}/{len(unique_peps):,} unique peptides "
        f"({frac_pep:.1%})  |  {rows_covered:,}/{len(df):,} rows ({frac_rows:.1%})"
    )


def assign_cluster_ids(df: pd.DataFrame, pep_clusters: pd.DataFrame) -> pd.DataFrame:
    if "cluster" in df.columns:
        df = df.drop(columns=["cluster"])

    df = df.merge(pep_clusters, on=COL_PEPTIDE, how="left")
    n_matched   = int(df["cluster"].notna().sum())
    n_unmatched = int(df["cluster"].isna().sum())
    log.info(f"  Peptide cluster merge: {n_matched:,} rows matched, {n_unmatched:,} unmatched")

    if n_unmatched > 0:
        existing_ids = df.loc[df["cluster"].notna(), "cluster"]
        max_existing = int(existing_ids.max()) if len(existing_ids) > 0 else -1

        unmatched_mask = df["cluster"].isna()
        anchor_pairs = df.loc[unmatched_mask, COL_PEPTIDE].apply(
            lambda p: (p[1] if len(p) >= 2 else "X") + (p[-1] if len(p) >= 1 else "X")
        )
        unique_pairs = sorted(set(anchor_pairs))
        pair_to_id   = {pair: max_existing + 1 + i for i, pair in enumerate(unique_pairs)}
        df.loc[unmatched_mask, "cluster"] = anchor_pairs.map(pair_to_id)
        log.info(
            f"  Assigned {n_unmatched:,} unmatched peptides to {len(unique_pairs)} "
            f"anchor-pair singleton clusters "
            f"(IDs {max_existing + 1}–{max_existing + len(unique_pairs)})"
        )

    df["cluster"] = df["cluster"].fillna(-1).astype(np.int32)
    n_noise     = int((df["cluster"] == -1).sum())
    n_clustered = int((df["cluster"] != -1).sum())
    n_distinct  = int(df.loc[df["cluster"] != -1, "cluster"].nunique())
    log.info(
        f"  Final: {n_clustered:,} rows with cluster, {n_noise:,} noise (-1), "
        f"{n_distinct} distinct clusters"
    )
    return df


def encode_pseudosequences(pseudo_csv_path: Path, alleles_in_data: set):
    """
    Load, clean-key-normalize, filter, and physicochemically encode pseudosequences.
    Returns (X, allele_list, alleles_missing).
    """
    log.info(f"Loading pseudosequences from {pseudo_csv_path}...")
    ps = pd.read_csv(pseudo_csv_path).drop_duplicates(subset=[COL_ALLELE]).reset_index(drop=True)
    log.info(f"  Pseudo CSV: {len(ps):,} alleles before clean_key()")

    ps["allele_clean"] = ps[COL_ALLELE].apply(clean_key)
    ps = ps.drop_duplicates(subset=["allele_clean"]).reset_index(drop=True)

    ps_filtered     = ps[ps["allele_clean"].isin(alleles_in_data)].reset_index(drop=True)
    alleles_missing = alleles_in_data - set(ps_filtered["allele_clean"])

    log.info(
        f"  Matched {len(ps_filtered):,}/{len(alleles_in_data):,} alleles after clean_key(); "
        f"{len(alleles_missing):,} have no pseudosequence (-> always train)"
    )
    if alleles_missing:
        sample = sorted(alleles_missing)[:10]
        log.info(f"  Missing sample: {sample}" + (" ..." if len(alleles_missing) > 10 else ""))

    X = encode_full_vectorized(ps_filtered["mhc_sequence"].values)
    log.info(f"  Encoded: {X.shape}  (49 positions x 14 props per allele)")
    return X, ps_filtered["allele_clean"].tolist(), alleles_missing


def run_hierarchical_clustering(
    X: np.ndarray, allele_list: list, leaf_distance_threshold: float, output_dir: Path
):
    """Ward linkage + distance-threshold fcluster; save dendrogram.png with cut line."""
    log.info(f"Running Ward hierarchical clustering (leaf_distance_threshold={leaf_distance_threshold})...")
    t0 = time.time()
    Z           = linkage(X, method="ward")
    hier_labels = fcluster(Z, t=leaf_distance_threshold, criterion="distance")

    unique, counts = np.unique(hier_labels, return_counts=True)
    log.info(f"  {len(unique)} clusters produced at threshold {leaf_distance_threshold}")
    for cl, cnt in zip(unique, counts):
        log.info(f"  Cluster {cl:3d}: {cnt:3d} alleles")

    n       = len(allele_list)
    leaf_fs = max(3, min(7, 350 // max(n, 1)))
    fig_w   = max(16, n * 0.12)
    fig, ax = plt.subplots(figsize=(fig_w, 7))
    dendrogram(
        Z,
        labels=allele_list,
        leaf_rotation=90,
        leaf_font_size=leaf_fs,
        ax=ax,
        color_threshold=leaf_distance_threshold,
    )
    ax.axhline(
        y=leaf_distance_threshold, color="red", linestyle="--", linewidth=1.5,
        label=f"Cut threshold = {leaf_distance_threshold}",
    )
    ax.legend(fontsize=10)
    ax.set_title(
        f"MHC Pseudosequence Dendrogram  (Ward linkage, threshold={leaf_distance_threshold})",
        fontsize=14,
    )
    ax.set_ylabel("Distance", fontsize=12)
    plt.tight_layout()
    out_path = output_dir / "dendrogram.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    log.info(f"  Saved: {out_path}  [{time.time()-t0:.1f}s]")
    return Z, hier_labels


def precompute_pep_cluster_maps(df: pd.DataFrame) -> tuple:
    """
    Build the data structures needed by peptide_cluster_set_cover.
    Returns (pep_clust_to_alleles, global_cluster_sizes,
             cluster_allele_rows, allele_cluster_rows, allele_total_rows)
    """
    _nc_df = df[df["cluster"] != -1]
    pep_clust_to_alleles = {
        c: set(grp[COL_ALLELE].unique())
        for c, grp in _nc_df.groupby("cluster")
    }
    global_cluster_sizes = _nc_df.groupby("cluster").size().to_dict()
    cluster_allele_rows: dict = {}
    allele_cluster_rows: dict = {}
    for (c, a), n in _nc_df.groupby(["cluster", COL_ALLELE]).size().items():
        cluster_allele_rows.setdefault(int(c), {})[a] = int(n)
        allele_cluster_rows.setdefault(a, {})[int(c)] = int(n)
    allele_total_rows = {
        a: int(sum(m.values())) for a, m in allele_cluster_rows.items()
    }
    log.info(
        f"Precomputed cluster maps: {len(pep_clust_to_alleles)} non-noise peptide clusters"
    )
    return (
        pep_clust_to_alleles,
        global_cluster_sizes,
        cluster_allele_rows,
        allele_cluster_rows,
        allele_total_rows,
    )


# ── species helpers ────────────────────────────────────────────────────────────

def get_species_prefix(cleaned_allele: str) -> str:
    """Return species prefix from a cleaned allele name."""
    if cleaned_allele.startswith("HLA"):
        return "HLA"
    for prefix in KNOWN_NONHLA_PREFIXES:
        if cleaned_allele.startswith(prefix):
            return prefix
    m = re.match(r"^([A-Z]+)", cleaned_allele)
    return m.group(1) if m else cleaned_allele[:4]


def detect_nonhuman_species(alleles: set) -> dict:
    """Group non-HLA alleles by species prefix. Returns {prefix: [allele, ...]}."""
    species_map: dict = {}
    for a in alleles:
        prefix = get_species_prefix(a)
        if prefix != "HLA":
            species_map.setdefault(prefix, []).append(a)
    return species_map


# ── monitoring helpers ─────────────────────────────────────────────────────────

def save_binder_nonbinder_counts(
    df: pd.DataFrame, out_dir: Path, label: str = "", save_plot: bool = False
) -> None:
    """Save per-allele binder/nonbinder counts CSV; log if ratio > 5:1; optional bar plot."""
    out_dir.mkdir(parents=True, exist_ok=True)
    grp = df.groupby(COL_ALLELE)["label"].value_counts().unstack(fill_value=0)
    grp.columns.name = None
    grp = grp.rename(columns={0: "nonbinders", 1: "binders"}).reset_index()
    for col in ("binders", "nonbinders"):
        if col not in grp.columns:
            grp[col] = 0
    grp["total"] = grp["binders"] + grp["nonbinders"]
    grp = grp.sort_values("total", ascending=False).reset_index(drop=True)
    grp.to_csv(out_dir / "binder_nonbinder_counts.csv", index=False)

    total_b = int(grp["binders"].sum())
    total_n = int(grp["nonbinders"].sum())
    if total_n > 0:
        ratio = total_b / total_n
        msg = f"  [{label}] binder:nonbinder = {total_b:,}:{total_n:,} ({ratio:.2f}:1)"
        if ratio > 5.0:
            log.warning(msg + "  << RATIO > 5:1")
        else:
            log.info(msg)
    else:
        log.info(f"  [{label}] binder:nonbinder = {total_b:,}:0  (no non-binders)")

    if save_plot and len(grp) > 0:
        fig, ax = plt.subplots(figsize=(max(8, len(grp) * 0.4), 5))
        x = np.arange(len(grp))
        ax.bar(x, grp["nonbinders"], label="Non-binders", color="steelblue")
        ax.bar(x, grp["binders"],    label="Binders",     color="coral", bottom=grp["nonbinders"])
        ax.set_xticks(x)
        ax.set_xticklabels(grp[COL_ALLELE], rotation=90, fontsize=6)
        ax.set_ylabel("Count")
        ax.set_title(f"Binder/Non-binder counts per allele  [{label}]")
        ax.legend()
        plt.tight_layout()
        plt.savefig(out_dir / "binder_nonbinder_counts.png", dpi=150, bbox_inches="tight")
        plt.close()


def save_peptide_cluster_summary(
    df: pd.DataFrame,
    sel_clusters: set,
    pep_clust_to_alleles: dict,
    global_cluster_sizes: dict,
    target_alleles: set,
    out_dir: Path,
    label: str = "",
) -> None:
    """Save peptide-cluster overlap summary CSV + bar plot."""
    out_dir.mkdir(parents=True, exist_ok=True)

    total_clusters      = len(pep_clust_to_alleles)
    overlapping         = {c for c, als in pep_clust_to_alleles.items() if als & target_alleles}
    excluded_by_sc      = overlapping & sel_clusters
    remaining_overlap   = overlapping - sel_clusters

    if remaining_overlap:
        log.warning(
            f"  [{label}] {len(remaining_overlap)} peptide clusters overlap target alleles "
            f"but were NOT excluded by set-cover (potential leakage)"
        )

    collateral_mask = df["cluster"].isin(sel_clusters) & ~df[COL_ALLELE].isin(target_alleles)
    collateral_rows = int(collateral_mask.sum())

    summary = {
        "total_clusters":         total_clusters,
        "overlapping_with_target": len(overlapping),
        "excluded_by_setcover":   len(excluded_by_sc),
        "remaining_overlap":      len(remaining_overlap),
        "collateral_rows_lost":   collateral_rows,
    }
    pd.DataFrame([summary]).to_csv(out_dir / "peptide_cluster_summary.csv", index=False)
    log.info(
        f"  [{label}] pep-cluster summary: total={total_clusters}, "
        f"overlap={len(overlapping)}, excluded={len(excluded_by_sc)}, "
        f"remaining={len(remaining_overlap)}, collateral={collateral_rows:,}"
    )

    categories = ["Total clusters", "Overlap w/ target", "Excluded (set-cover)",
                  "Remaining overlap", "Collateral rows lost"]
    values     = [total_clusters, len(overlapping), len(excluded_by_sc),
                  len(remaining_overlap), collateral_rows]
    colors     = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3"]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(categories, values, color=colors)
    ax.set_ylabel("Count")
    ax.set_title(f"Peptide cluster summary  [{label}]")
    for i, v in enumerate(values):
        ax.text(i, v + max(values) * 0.01, str(v), ha="center", fontsize=9)
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(out_dir / "peptide_cluster_summary.png", dpi=150, bbox_inches="tight")
    plt.close()


# ── set-cover call helpers ─────────────────────────────────────────────────────

def _rows_by_cluster_for_alleles(alleles, allele_cluster_rows):
    rows: dict = {}
    for a in alleles:
        for c, n in allele_cluster_rows.get(a, {}).items():
            rows[c] = rows.get(c, 0) + n
    return rows


def _per_allele_floor(alleles, allele_total_rows):
    return {
        a: min(
            int(np.ceil(PEPTIDE_COVER_PER_ALLELE_FRAC * allele_total_rows.get(a, 0))),
            PEPTIDE_COVER_PER_ALLELE_CAP,
        )
        for a in alleles
    }


# ── per-threshold data loader ──────────────────────────────────────────────────

def load_threshold_data(df_base: pd.DataFrame, pep_cluster_file: Path):
    """Load one peptide cluster file, assign IDs, precompute maps."""
    pep_clusters = load_precomputed_peptide_clusters(pep_cluster_file)
    print_peptide_cluster_coverage(df_base, pep_clusters)
    df = assign_cluster_ids(df_base.copy(), pep_clusters)
    maps = precompute_pep_cluster_maps(df)
    return df, maps


# ── mode: hla_a / hla_b / hla_c ───────────────────────────────────────────────

def run_mode_hla_family(
    df_base: pd.DataFrame,
    allele_list: list,
    hier_labels: np.ndarray,
    target_prefix: str,
    pep_cluster_files: list,
    thresholds: list,
    mode_output_dir: Path,
    n_folds: int,
):
    """Run hla_a / hla_b / hla_c mode for all peptide-cluster thresholds."""
    log.info(f"\n{'='*60}")
    log.info(f"Mode: {target_prefix} family split")
    log.info(f"{'='*60}")

    # ── Allele assignment: done once from the dendrogram ──────────────────────
    unique_clusters = sorted(np.unique(hier_labels).tolist())
    clusters_with_target = []
    for cl in unique_clusters:
        cl_alleles = [a for a, c in zip(allele_list, hier_labels) if c == cl]
        n_target   = sum(1 for a in cl_alleles if a.startswith(target_prefix))
        if n_target > 0:
            purity = n_target / len(cl_alleles)
            clusters_with_target.append({
                "label":      cl,
                "purity":     purity,
                "all":        cl_alleles,
                "target":     [a for a in cl_alleles if a.startswith(target_prefix)],
                "n_target":   n_target,
            })

    if not clusters_with_target:
        log.warning(f"No clusters contain {target_prefix} alleles. Skipping mode.")
        return

    clusters_with_target.sort(key=lambda x: -x["purity"])

    test_info            = clusters_with_target[0]
    test_all_cluster     = set(test_info["all"])
    test_alleles         = set(test_info["target"])
    log.info(
        f"Test cluster: label={test_info['label']}, purity={test_info['purity']:.2f}, "
        f"{len(test_alleles)} target alleles, {len(test_all_cluster)} total alleles"
    )
    log.info(f"  Test alleles  : {sorted(test_alleles)}")
    log.info(f"  Discarded     : {sorted(test_all_cluster - test_alleles)}")

    remaining = clusters_with_target[1:]
    max_folds = n_folds if n_folds > 0 else len(remaining)
    fold_groups = [remaining[i:i+3] for i in range(0, min(len(remaining), max_folds * 3), 3)]
    log.info(f"Val fold groups: {len(fold_groups)} folds (3 clusters/fold max, n_folds cap={max_folds})")
    for fi, grp in enumerate(fold_groups):
        target_in_grp = [a for info in grp for a in info["target"]]
        log.info(f"  Fold {fi}: clusters={[g['label'] for g in grp]}, "
                 f"{len(target_in_grp)} target val alleles")

    # ── Per-threshold loop ────────────────────────────────────────────────────
    for pep_file, thr_label in zip(pep_cluster_files, thresholds):
        log.info(f"\n  -- threshold {thr_label}  ({pep_file.name}) --")
        thr_dir = mode_output_dir / f"threshold_{thr_label}"

        df, (pep_clust_to_alleles, global_cluster_sizes,
             cluster_allele_rows, allele_cluster_rows, allele_total_rows) = \
            load_threshold_data(df_base, pep_file)

        # Test set-cover on ALL alleles in test cluster (prevents peptide leakage)
        log.info("  Running set-cover for TEST alleles...")
        sel_test_clust, cov_test_a, test_st = peptide_cluster_set_cover(
            target_alleles      = test_all_cluster,
            pep_clust_to_alleles= pep_clust_to_alleles,
            global_cluster_sizes= global_cluster_sizes,
            target_cluster_rows = _rows_by_cluster_for_alleles(test_all_cluster, allele_cluster_rows),
            cluster_allele_rows = cluster_allele_rows,
            per_allele_min_rows = _per_allele_floor(test_all_cluster, allele_total_rows),
            min_target_rows     = PEPTIDE_COVER_MIN_ROWS,
        )
        excluded_test_clusters = sel_test_clust
        log.info(
            f"  Test set-cover: {len(sel_test_clust)} clusters selected, "
            f"coverage={test_st['coverage_frac']:.1%}"
        )

        test_dir  = thr_dir / "test"
        test_mask = (
            df[COL_ALLELE].isin(test_alleles)
            & df["cluster"].isin(excluded_test_clusters)
        )
        save_parquet(df[test_mask].reset_index(drop=True), test_dir / "test.parquet", "test")
        pd.DataFrame({"allele": sorted(test_alleles)}).to_csv(
            test_dir / "test_alleles.csv", index=False
        )
        save_binder_nonbinder_counts(df[test_mask], test_dir, label=f"test thr={thr_label}", save_plot=True)
        save_peptide_cluster_summary(
            df, excluded_test_clusters, pep_clust_to_alleles, global_cluster_sizes,
            test_alleles, test_dir, label=f"test thr={thr_label}",
        )

        fold_summary = []
        for fold_i, fold_group in enumerate(fold_groups):
            val_all_cluster = set(a for info in fold_group for a in info["all"])
            val_alleles     = set(a for info in fold_group for a in info["target"])

            log.info(f"\n  Fold {fold_i}: {len(val_alleles)} target val alleles, "
                     f"{len(val_all_cluster)} total in val clusters")

            sel_val_clust, cov_val_a, val_st = peptide_cluster_set_cover(
                target_alleles      = val_all_cluster,
                pep_clust_to_alleles= pep_clust_to_alleles,
                global_cluster_sizes= global_cluster_sizes,
                target_cluster_rows = _rows_by_cluster_for_alleles(val_all_cluster, allele_cluster_rows),
                cluster_allele_rows = cluster_allele_rows,
                per_allele_min_rows = _per_allele_floor(val_all_cluster, allele_total_rows),
                excluded_clusters   = excluded_test_clusters,
                min_target_rows     = PEPTIDE_COVER_MIN_ROWS,
            )
            log.info(
                f"  Val set-cover: {len(sel_val_clust)} clusters selected, "
                f"coverage={val_st['coverage_frac']:.1%}"
            )

            in_test_cluster = df["cluster"].isin(excluded_test_clusters)
            in_val_cluster  = df["cluster"].isin(sel_val_clust)
            in_test_allele  = df[COL_ALLELE].isin(test_all_cluster)
            in_val_allele   = df[COL_ALLELE].isin(val_all_cluster)

            val_mask   = df[COL_ALLELE].isin(val_alleles) & in_val_cluster
            train_mask = ~in_test_cluster & ~in_val_cluster & ~in_test_allele & ~in_val_allele

            df_val   = df[val_mask].reset_index(drop=True)
            df_train = df[train_mask].reset_index(drop=True)

            fold_dir = thr_dir / f"fold_{fold_i}"
            save_parquet(df_train, fold_dir / "train.parquet", f"fold_{fold_i} train")
            save_parquet(df_val,   fold_dir / "val.parquet",   f"fold_{fold_i} val")
            save_binder_nonbinder_counts(df_train, fold_dir, label=f"fold_{fold_i}/train thr={thr_label}")
            save_binder_nonbinder_counts(df_val,   fold_dir / "val_counts",
                                         label=f"fold_{fold_i}/val thr={thr_label}", save_plot=True)
            save_peptide_cluster_summary(
                df, sel_val_clust, pep_clust_to_alleles, global_cluster_sizes,
                val_alleles, fold_dir, label=f"fold_{fold_i} thr={thr_label}",
            )

            n_excluded = int(
                ((in_test_cluster & ~in_test_allele) | (in_val_cluster & ~in_val_allele)).sum()
            )
            fold_summary.append({
                "fold":                             fold_i,
                "val_cluster_labels":               str([g["label"] for g in fold_group]),
                "n_val_target_alleles":             len(val_alleles),
                "n_val_all_cluster_alleles":        len(val_all_cluster),
                "n_train":                          len(df_train),
                "n_val":                            len(df_val),
                "n_train_binders":                  int((df_train["label"] == 1).sum()),
                "n_train_nonbinders":               int((df_train["label"] == 0).sum()),
                "n_val_binders":                    int((df_val["label"] == 1).sum()),
                "n_val_nonbinders":                 int((df_val["label"] == 0).sum()),
                "n_pep_clusters_test_excluded":     len(excluded_test_clusters),
                "n_pep_clusters_val_excluded":      len(sel_val_clust),
                "n_rows_excluded_collateral":       n_excluded,
                "val_setcover_coverage":            round(val_st["coverage_frac"], 4),
            })

        pd.DataFrame(fold_summary).to_csv(thr_dir / "fold_summary.csv", index=False)
        log.info(f"  Saved: {thr_dir / 'fold_summary.csv'}")

        meta = {
            "mode":                   target_prefix,
            "threshold":              thr_label,
            "pep_cluster_file":       str(pep_file),
            "test_cluster_label":     test_info["label"],
            "test_cluster_purity":    round(test_info["purity"], 4),
            "test_alleles":           sorted(test_alleles),
            "test_all_cluster_alleles": sorted(test_all_cluster),
            "n_folds":                len(fold_groups),
        }
        with open(thr_dir / "split_meta.json", "w") as fh:
            json.dump(meta, fh, indent=2)


# ── mode: cross_species ────────────────────────────────────────────────────────

def run_mode_cross_species(
    df_base: pd.DataFrame,
    pep_cluster_files: list,
    thresholds: list,
    mode_output_dir: Path,
):
    """Run cross_species mode: train=HLA, test=one species, val=other non-human species."""
    log.info(f"\n{'='*60}")
    log.info("Mode: cross_species")
    log.info(f"{'='*60}")

    all_alleles  = set(df_base[COL_ALLELE].unique())
    species_map  = detect_nonhuman_species(all_alleles)
    hla_alleles  = {a for a in all_alleles if get_species_prefix(a) == "HLA"}

    log.info(f"HLA alleles: {len(hla_alleles)}")
    for sp, als in sorted(species_map.items()):
        log.info(f"  Species {sp}: {len(als)} alleles")

    if not species_map:
        log.warning("No non-HLA species found in data. Skipping cross_species mode.")
        return

    for species_prefix, species_alleles_list in sorted(species_map.items()):
        test_alleles = set(species_alleles_list)
        other_nonhuman = {
            a for sp, als in species_map.items()
            if sp != species_prefix
            for a in als
        }
        log.info(
            f"\nSpecies {species_prefix}: test={len(test_alleles)} alleles, "
            f"val(other non-human)={len(other_nonhuman)} alleles, "
            f"train(HLA)={len(hla_alleles)} alleles"
        )

        for pep_file, thr_label in zip(pep_cluster_files, thresholds):
            log.info(f"  -- threshold {thr_label} --")
            thr_dir = mode_output_dir / species_prefix / f"threshold_{thr_label}"

            df, (pep_clust_to_alleles, global_cluster_sizes,
                 cluster_allele_rows, allele_cluster_rows, allele_total_rows) = \
                load_threshold_data(df_base, pep_file)

            # Test set-cover
            log.info(f"  Running set-cover for TEST ({species_prefix})...")
            sel_test_clust, _, test_st = peptide_cluster_set_cover(
                target_alleles      = test_alleles,
                pep_clust_to_alleles= pep_clust_to_alleles,
                global_cluster_sizes= global_cluster_sizes,
                target_cluster_rows = _rows_by_cluster_for_alleles(test_alleles, allele_cluster_rows),
                cluster_allele_rows = cluster_allele_rows,
                per_allele_min_rows = _per_allele_floor(test_alleles, allele_total_rows),
                min_target_rows     = PEPTIDE_COVER_MIN_ROWS,
            )
            excluded_test_clusters = sel_test_clust
            log.info(
                f"  Test set-cover: {len(sel_test_clust)} clusters selected, "
                f"coverage={test_st['coverage_frac']:.1%}"
            )

            # Val set-cover (other non-human species)
            val_alleles = other_nonhuman
            if val_alleles:
                log.info("  Running set-cover for VAL (other non-human)...")
                sel_val_clust, _, val_st = peptide_cluster_set_cover(
                    target_alleles      = val_alleles,
                    pep_clust_to_alleles= pep_clust_to_alleles,
                    global_cluster_sizes= global_cluster_sizes,
                    target_cluster_rows = _rows_by_cluster_for_alleles(val_alleles, allele_cluster_rows),
                    cluster_allele_rows = cluster_allele_rows,
                    per_allele_min_rows = _per_allele_floor(val_alleles, allele_total_rows),
                    excluded_clusters   = excluded_test_clusters,
                    min_target_rows     = PEPTIDE_COVER_MIN_ROWS,
                )
                log.info(
                    f"  Val set-cover: {len(sel_val_clust)} clusters selected, "
                    f"coverage={val_st['coverage_frac']:.1%}"
                )
            else:
                sel_val_clust = set()
                val_st = {"coverage_frac": 0.0}

            in_test_cluster = df["cluster"].isin(excluded_test_clusters)
            in_val_cluster  = df["cluster"].isin(sel_val_clust)
            in_test_allele  = df[COL_ALLELE].isin(test_alleles)
            in_val_allele   = df[COL_ALLELE].isin(val_alleles)

            test_mask  = in_test_allele & in_test_cluster
            val_mask   = in_val_allele  & in_val_cluster
            train_mask = ~in_test_cluster & ~in_val_cluster & ~in_test_allele & ~in_val_allele

            df_test  = df[test_mask].reset_index(drop=True)
            df_val   = df[val_mask].reset_index(drop=True)
            df_train = df[train_mask].reset_index(drop=True)

            test_dir = thr_dir / "test"
            fold_dir = thr_dir / "fold_0"

            save_parquet(df_test,  test_dir / "test.parquet",    f"{species_prefix} test")
            save_parquet(df_train, fold_dir  / "train.parquet",  f"{species_prefix} fold_0 train")
            save_parquet(df_val,   fold_dir  / "val.parquet",    f"{species_prefix} fold_0 val")

            pd.DataFrame({"allele": sorted(test_alleles)}).to_csv(
                test_dir / "test_alleles.csv", index=False
            )

            save_binder_nonbinder_counts(df_test,  test_dir, label=f"{species_prefix} test thr={thr_label}", save_plot=True)
            save_binder_nonbinder_counts(df_train, fold_dir, label=f"{species_prefix} train thr={thr_label}")
            save_binder_nonbinder_counts(df_val,   fold_dir / "val_counts",
                                         label=f"{species_prefix} val thr={thr_label}", save_plot=True)

            save_peptide_cluster_summary(
                df, excluded_test_clusters, pep_clust_to_alleles, global_cluster_sizes,
                test_alleles, test_dir, label=f"{species_prefix} test thr={thr_label}",
            )
            if val_alleles:
                save_peptide_cluster_summary(
                    df, sel_val_clust, pep_clust_to_alleles, global_cluster_sizes,
                    val_alleles, fold_dir, label=f"{species_prefix} val thr={thr_label}",
                )

            n_excluded = int(
                ((in_test_cluster & ~in_test_allele) | (in_val_cluster & ~in_val_allele)).sum()
            )
            fold_summary = [{
                "fold":                          0,
                "species_test":                  species_prefix,
                "n_test_alleles":                len(test_alleles),
                "n_val_alleles":                 len(val_alleles),
                "n_train_alleles":               len(hla_alleles),
                "n_test":                        len(df_test),
                "n_train":                       len(df_train),
                "n_val":                         len(df_val),
                "n_test_binders":                int((df_test["label"]  == 1).sum()),
                "n_test_nonbinders":             int((df_test["label"]  == 0).sum()),
                "n_train_binders":               int((df_train["label"] == 1).sum()),
                "n_train_nonbinders":            int((df_train["label"] == 0).sum()),
                "n_val_binders":                 int((df_val["label"]   == 1).sum()),
                "n_val_nonbinders":              int((df_val["label"]   == 0).sum()),
                "n_pep_clusters_test_excluded":  len(excluded_test_clusters),
                "n_pep_clusters_val_excluded":   len(sel_val_clust),
                "n_rows_excluded_collateral":    n_excluded,
                "test_setcover_coverage":        round(test_st["coverage_frac"], 4),
                "val_setcover_coverage":         round(val_st["coverage_frac"], 4),
            }]
            pd.DataFrame(fold_summary).to_csv(thr_dir / "fold_summary.csv", index=False)

            meta = {
                "mode":           "cross_species",
                "species":        species_prefix,
                "threshold":      thr_label,
                "pep_cluster_file": str(pep_file),
                "test_alleles":   sorted(test_alleles),
                "val_alleles":    sorted(val_alleles),
                "train_alleles":  sorted(hla_alleles),
            }
            with open(thr_dir / "split_meta.json", "w") as fh:
                json.dump(meta, fh, indent=2)
            log.info(f"  Saved: {thr_dir / 'fold_summary.csv'}")


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Dendrogram-stratified train/val/test splits with peptide-level "
            "leakage control.  Modes: hla_a, hla_b, hla_c, cross_species."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input_dir",  type=Path, required=True,
                   help="Input parquet (file or directory)")
    p.add_argument("--pseudo_csv", type=Path, required=True,
                   help="mhc_pseudo_class1_n.csv with 49-char mhc_sequence column")
    p.add_argument("--pep_cluster_files", type=Path, nargs="+", required=True,
                   help="Asgary TSV files, one per threshold (must match --thresholds order)")
    p.add_argument("--thresholds", type=str, nargs="+", default=["50", "60", "70", "80"],
                   help="Threshold labels matching --pep_cluster_files")
    p.add_argument("--output_dir", type=Path, required=True,
                   help="Root output directory")
    p.add_argument("--mode", required=True,
                   choices=["hla_a", "hla_b", "hla_c", "cross_species"],
                   help="Experiment mode")
    p.add_argument("--leaf_distance_threshold", type=float, default=25.0,
                   help="Distance threshold for dendrogram cut (hla_* modes)")
    p.add_argument("--n_folds", type=int, default=5,
                   help="Max number of validation folds (hla_* modes; 3 leaf clusters per fold)")
    return p.parse_args()


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    if len(args.pep_cluster_files) != len(args.thresholds):
        raise ValueError(
            f"--pep_cluster_files ({len(args.pep_cluster_files)}) and "
            f"--thresholds ({len(args.thresholds)}) must have the same length."
        )

    output_dir = args.output_dir
    setup_logging(output_dir)
    log.info("=" * 60)
    log.info("dendrogram_split.py")
    log.info("=" * 60)
    log.info(f"Arguments: {vars(args)}")

    t_total = time.time()

    # ── Step 1: load data ──────────────────────────────────────────────────────
    df_base = load(args.input_dir, output_dir)

    # ── Step 2: encode pseudosequences ─────────────────────────────────────────
    alleles_in_data = set(df_base[COL_ALLELE].unique())
    X_mhc, allele_list, alleles_missing = encode_pseudosequences(
        args.pseudo_csv, alleles_in_data
    )

    # ── Step 3: Ward clustering + dendrogram (once) ────────────────────────────
    mhc_cluster_dir = output_dir / "mhc_clustering"
    Z, hier_labels = run_hierarchical_clustering(
        X_mhc, allele_list, args.leaf_distance_threshold, mhc_cluster_dir
    )

    # Save allele→cluster assignment (fixed across all thresholds)
    assignments_df = pd.DataFrame({
        "allele":          allele_list,
        "cluster_label":   hier_labels.tolist(),
    })
    assignments_path = mhc_cluster_dir / "allele_cluster_assignments.csv"
    assignments_df.to_csv(assignments_path, index=False)
    log.info(f"Saved: {assignments_path}")

    # ── Step 4: dispatch to mode ───────────────────────────────────────────────
    if args.mode in ("hla_a", "hla_b", "hla_c"):
        prefix_map = {"hla_a": "HLAA", "hla_b": "HLAB", "hla_c": "HLAC"}
        target_prefix = prefix_map[args.mode]
        run_mode_hla_family(
            df_base          = df_base,
            allele_list      = allele_list,
            hier_labels      = hier_labels,
            target_prefix    = target_prefix,
            pep_cluster_files= args.pep_cluster_files,
            thresholds       = args.thresholds,
            mode_output_dir  = output_dir / args.mode,
            n_folds          = args.n_folds,
        )
    else:  # cross_species
        run_mode_cross_species(
            df_base          = df_base,
            pep_cluster_files= args.pep_cluster_files,
            thresholds       = args.thresholds,
            mode_output_dir  = output_dir / "cross_species",
        )

    log.info(f"\nDone.  Total: {time.time()-t_total:.1f}s")


if __name__ == "__main__":
    main()
