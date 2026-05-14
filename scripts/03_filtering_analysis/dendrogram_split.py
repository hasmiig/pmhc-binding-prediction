"""
dendrogram_split.py
===================
Creates MHC-dendrogram-stratified train/val/test splits with peptide-level
leakage control.

Steps
-----
1. Load combined parquet.
2. Load --pep_cluster_file, print peptide-cluster coverage of combined frame.
3. Assign peptide cluster IDs: use precomputed clusters where available;
   for unmatched peptides, create anchor-pair singleton clusters (P2 + last residue).
4. Load --pseudo_csv (49-char mhc_sequence), apply clean_key() to match parquet
   allele names, encode with physicochemical properties (49x14 = 686 features).
5. Ward hierarchical clustering on encoded alleles, cut into --k_clusters groups.
   Save dendrogram.png.
6. Hold out one cluster as the fixed test set (--test_cluster or largest).
7. Run peptide_cluster_set_cover for test alleles; save test.parquet.
   excluded_test_clusters are permanently off-limits for all folds.
8. Distribute remaining clusters across --n_folds folds (merging smallest if needed).
9. For each fold: run set_cover for val alleles (excluding test clusters),
   assign rows, save train/val parquets.
10. Save fold_summary.csv and split_meta.json.

Row assignment per fold
-----------------------
- allele in test_alleles AND cluster in excluded_test_clusters  -> test.parquet
- allele in val_alleles  AND cluster in sel_val_clusters        -> val.parquet
- cluster in excluded_test_clusters AND allele NOT in test      -> excluded (neither split)
- cluster in sel_val_clusters       AND allele NOT in val       -> excluded (neither split)
- everything else                                               -> train.parquet
  (includes test/val allele rows whose peptide cluster was NOT selected by set cover)

Usage
-----
python dendrogram_split.py \\
    --input_dir        .../combined.parquet \\
    --pseudo_csv       .../mhc_pseudo_class1_n.csv \\
    --pep_cluster_file .../pep_clusters.tsv \\
    --output_dir       .../output_dendro_split \\
    --k_clusters 10 \\
    --n_folds 5
"""

import argparse
import json
import logging
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
        log.info("  Dropping existing 'cluster' column from combined frame before merge.")
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

    Parameters
    ----------
    pseudo_csv_path : Path
        CSV with columns [allele, mhc_sequence] (49-char, already aligned).
    alleles_in_data : set
        Cleaned allele keys as stored in the combined parquet.

    Returns
    -------
    X            : np.ndarray (n_alleles, 49*14)
    allele_list  : list[str]  — cleaned allele keys, aligned with X rows
    alleles_missing : set[str] — alleles in data with no pseudosequence (-> always train)
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

    X = encode_full_vectorized(ps_filtered["mhc_sequence"].values)  # (n, 49*14)
    log.info(f"  Encoded: {X.shape}  (49 positions x 14 props per allele)")
    return X, ps_filtered["allele_clean"].tolist(), alleles_missing


def run_hierarchical_clustering(
    X: np.ndarray, allele_list: list, k_clusters: int, output_dir: Path
):
    """Ward linkage + fcluster; save dendrogram.png. Returns (Z, hier_labels)."""
    log.info(f"Running Ward hierarchical clustering (k={k_clusters})...")
    t0 = time.time()
    Z           = linkage(X, method="ward")
    hier_labels = fcluster(Z, t=k_clusters, criterion="maxclust")  # 1-indexed int array

    unique, counts = np.unique(hier_labels, return_counts=True)
    for cl, cnt in zip(unique, counts):
        log.info(f"  Cluster {cl:2d}: {cnt:3d} alleles")

    # Dendrogram
    n          = len(allele_list)
    leaf_fs    = max(3, min(7, 350 // max(n, 1)))
    fig_w      = max(16, n * 0.12)
    color_thr  = float(Z[-(k_clusters - 1), 2]) if k_clusters > 1 else 0.0
    fig, ax    = plt.subplots(figsize=(fig_w, 7))
    dendrogram(
        Z,
        labels=allele_list,
        leaf_rotation=90,
        leaf_font_size=leaf_fs,
        ax=ax,
        color_threshold=color_thr,
    )
    ax.set_title(
        f"MHC Pseudosequence Dendrogram  (Ward linkage, k={k_clusters})", fontsize=14
    )
    ax.set_ylabel("Distance", fontsize=12)
    plt.tight_layout()
    out_path = output_dir / "dendrogram.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    log.info(f"  Saved: {out_path}  [{time.time()-t0:.1f}s]")
    return Z, hier_labels


def select_test_cluster(
    hier_labels: np.ndarray, allele_list: list, test_cluster_arg
) -> int:
    """Return the cluster label to hold out as the fixed test set."""
    unique, counts = np.unique(hier_labels, return_counts=True)
    if test_cluster_arg is not None:
        if test_cluster_arg not in unique:
            raise ValueError(
                f"--test_cluster {test_cluster_arg} is not a valid cluster label. "
                f"Valid labels: {sorted(unique.tolist())}"
            )
        chosen = int(test_cluster_arg)
        cnt    = int(counts[unique == chosen][0])
    else:
        idx    = int(np.argmax(counts))
        chosen = int(unique[idx])
        cnt    = int(counts[idx])
    log.info(f"Test cluster: label={chosen}  ({cnt} alleles)")
    return chosen


def merge_clusters_to_folds(
    remaining_cluster_ids: list,
    hier_labels: np.ndarray,
    allele_list: list,
    n_folds: int,
) -> list:
    """
    Distribute remaining MHC cluster labels into min(n_folds, len(remaining)) groups.
    When more clusters than folds, iteratively merge the two smallest groups.

    Returns
    -------
    list of sets, each set containing one or more cluster label ints.
    """
    remaining_set = set(remaining_cluster_ids)

    # Count alleles per remaining cluster
    cluster_size_map: dict = {}
    for a, lbl in zip(allele_list, hier_labels):
        if lbl in remaining_set:
            cluster_size_map[lbl] = cluster_size_map.get(lbl, 0) + 1

    groups = [{c} for c in sorted(remaining_set)]
    sizes  = [cluster_size_map.get(c, 0) for c in sorted(remaining_set)]

    while len(groups) > n_folds:
        i1        = int(np.argmin(sizes))
        saved     = sizes[i1]
        sizes[i1] = float("inf")
        i2        = int(np.argmin(sizes))
        sizes[i1] = saved
        groups[i2] |= groups[i1]
        sizes[i2]  += sizes[i1]
        groups.pop(i1)
        sizes.pop(i1)

    actual = len(groups)
    log.info(
        f"Fold groups: {actual} folds "
        f"(requested {n_folds}, {len(remaining_cluster_ids)} remaining clusters)"
    )
    for fi, grp in enumerate(groups):
        n_alleles = sum(cluster_size_map.get(c, 0) for c in grp)
        log.info(f"  Fold {fi}: clusters={sorted(grp)}  alleles={n_alleles}")
    return groups


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


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Dendrogram-stratified train/val/test splits with peptide-level "
            "leakage control."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--input_dir", type=Path, required=True,
        help="Input directory containing the combined binder/nonbinder dataset",
    )

    p.add_argument(
        "--pseudo_csv", type=Path, required=True,
        help="mhc_pseudo_class1_n.csv with 49-char mhc_sequence column",
    )
    p.add_argument(
        "--pep_cluster_file", type=Path, required=True,
        help="Asgary TSV with cluster_id and sequence columns",
    )
    p.add_argument(
        "--output_dir", type=Path, required=True,
        help="Root output directory",
    )
    p.add_argument(
        "--k_clusters", type=int, default=10,
        help="Number of MHC hierarchical clusters",
    )
    p.add_argument(
        "--test_cluster", type=int, default=None,
        help="Cluster label (1-indexed) to hold out as test; default = largest cluster",
    )
    p.add_argument(
        "--n_folds", type=int, default=5,
        help="Number of CV folds from remaining (non-test) clusters",
    )
    return p.parse_args()


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    args       = parse_args()
    output_dir = args.output_dir
    splits_dir = output_dir / "splits"

    setup_logging(output_dir)
    log.info("=" * 60)
    log.info("dendrogram_split.py")
    log.info("=" * 60)
    log.info(f"Arguments: {vars(args)}")

    t_total = time.time()

    # ── Step 1: load ───────────────────────────────────────────────
    df = load(args.input_dir, output_dir)

    # ── Step 2: load peptide clusters, print coverage ──────────────────────────
    pep_clusters = load_precomputed_peptide_clusters(args.pep_cluster_file)
    print_peptide_cluster_coverage(df, pep_clusters)

    # ── Step 3: assign cluster IDs (file + singleton fallback) ─────────────────
    df = assign_cluster_ids(df, pep_clusters)

    # ── Step 4: encode pseudosequences ─────────────────────────────────────────
    alleles_in_data = set(df[COL_ALLELE].unique())
    X_mhc, allele_list, alleles_missing = encode_pseudosequences(
        args.pseudo_csv, alleles_in_data
    )

    # ── Step 5: Ward clustering + dendrogram ───────────────────────────────────
    _, hier_labels = run_hierarchical_clustering(
        X_mhc, allele_list, args.k_clusters, output_dir
    )

    # ── Step 6: select test cluster ────────────────────────────────────────────
    test_cluster_label = select_test_cluster(
        hier_labels, allele_list, args.test_cluster
    )
    test_alleles = {
        a for a, c in zip(allele_list, hier_labels) if c == test_cluster_label
    }

    # ── Step 7: precompute peptide cluster maps ─────────────────────────────────
    (
        pep_clust_to_alleles,
        global_cluster_sizes,
        cluster_allele_rows,
        allele_cluster_rows,
        allele_total_rows,
    ) = precompute_pep_cluster_maps(df)

    def rows_by_cluster_for_alleles(alleles):
        rows: dict = {}
        for a in alleles:
            for c, n in allele_cluster_rows.get(a, {}).items():
                rows[c] = rows.get(c, 0) + n
        return rows

    def per_allele_floor(alleles):
        return {
            a: min(
                int(np.ceil(PEPTIDE_COVER_PER_ALLELE_FRAC * allele_total_rows.get(a, 0))),
                PEPTIDE_COVER_PER_ALLELE_CAP,
            )
            for a in alleles
        }

    # ── Step 8: set cover for test alleles ─────────────────────────────────────
    log.info("\nRunning peptide cluster set cover for TEST alleles...")
    sel_test_clust, cov_test_a, test_st = peptide_cluster_set_cover(
        target_alleles      = test_alleles,
        pep_clust_to_alleles= pep_clust_to_alleles,
        global_cluster_sizes= global_cluster_sizes,
        target_cluster_rows = rows_by_cluster_for_alleles(test_alleles),
        cluster_allele_rows = cluster_allele_rows,
        per_allele_min_rows = per_allele_floor(test_alleles),
        min_target_rows     = PEPTIDE_COVER_MIN_ROWS,
    )
    excluded_test_clusters = sel_test_clust
    log.info(
        f"  Test: {len(sel_test_clust)} pep-clusters selected, "
        f"{len(cov_test_a)}/{len(test_alleles)} alleles covered, "
        f"coverage={test_st['coverage_frac']:.1%}, "
        f"target_rows={test_st['target_rows_selected']:,}"
    )

    # Save test split
    test_mask = (
        df[COL_ALLELE].isin(test_alleles)
        & df["cluster"].isin(excluded_test_clusters)
    )
    #test_mask = df[COL_ALLELE].isin(test_alleles)
    save_parquet(
        df[test_mask].reset_index(drop=True),
        splits_dir / "test" / "test.parquet",
        "test",
    )
    pd.DataFrame({"allele": sorted(test_alleles)}).to_csv(
        splits_dir / "test" / "test_alleles.csv", index=False
    )
    log.info(f"  Saved: {splits_dir / 'test' / 'test_alleles.csv'}")

    # ── Step 9: distribute remaining clusters into fold groups ─────────────────
    remaining_cluster_labels = sorted(set(hier_labels.tolist()) - {test_cluster_label})
    fold_groups = merge_clusters_to_folds(
        remaining_cluster_labels, hier_labels, allele_list, args.n_folds
    )

    # ── Step 10: per-fold set cover + splits ───────────────────────────────────
    fold_summary = []
    for fold_i, val_cluster_group in enumerate(fold_groups):
        val_alleles = {
            a for a, c in zip(allele_list, hier_labels) if c in val_cluster_group
        }
        log.info(
            f"\nFold {fold_i}: val cluster group={sorted(val_cluster_group)}, "
            f"{len(val_alleles)} val alleles"
        )

        sel_val_clust, cov_val_a, val_st = peptide_cluster_set_cover(
            target_alleles      = val_alleles,
            pep_clust_to_alleles= pep_clust_to_alleles,
            global_cluster_sizes= global_cluster_sizes,
            target_cluster_rows = rows_by_cluster_for_alleles(val_alleles),
            cluster_allele_rows = cluster_allele_rows,
            per_allele_min_rows = per_allele_floor(val_alleles),
            excluded_clusters   = excluded_test_clusters,
            min_target_rows     = PEPTIDE_COVER_MIN_ROWS,
        )
        log.info(
            f"  Val set cover: {len(sel_val_clust)} pep-clusters selected, "
            f"{len(cov_val_a)}/{len(val_alleles)} alleles covered, "
            f"coverage={val_st['coverage_frac']:.1%}, "
            f"target_rows={val_st['target_rows_selected']:,}"
        )

        # Boolean masks
        in_test_cluster = df["cluster"].isin(excluded_test_clusters)
        in_val_cluster  = df["cluster"].isin(sel_val_clust)
        in_test_allele  = df[COL_ALLELE].isin(test_alleles)
        in_val_allele   = df[COL_ALLELE].isin(val_alleles)

        val_mask   = in_val_allele & in_val_cluster
        # val_mask = in_val_allele
        train_mask = ~in_test_cluster & ~in_val_cluster
        # train_mask = ~in_test_cluster & ~in_val_cluster & ~in_test_allele & ~in_val_allele
        # Excluded rows (neither test, val, nor train):
        #   (in_test_cluster & ~in_test_allele)  — cross-contaminated by test clusters
        #   (in_val_cluster  & ~in_val_allele)   — cross-contaminated by val clusters

        n_excluded = int(
            ((in_test_cluster & ~in_test_allele) | (in_val_cluster & ~in_val_allele)).sum()
        )
        n_recycled = int(((in_test_allele | in_val_allele) & train_mask).sum())

        df_val   = df[val_mask].reset_index(drop=True)
        df_train = df[train_mask].reset_index(drop=True)

        fold_dir = splits_dir / f"fold_{fold_i}"
        save_parquet(df_train, fold_dir / "train.parquet", f"fold_{fold_i} train")
        save_parquet(df_val,   fold_dir / "val.parquet",   f"fold_{fold_i} val")
        log.info(
            f"  fold_{fold_i}: train={len(df_train):,}  val={len(df_val):,}  "
            f"excluded={n_excluded:,}  recycled_to_train={n_recycled:,}"
        )

        fold_summary.append({
            "fold":                             fold_i,
            "n_train":                          len(df_train),
            "n_val":                            len(df_val),
            "n_train_binders":                  int((df_train["label"] == 1).sum()),
            "n_train_nonbinders":               int((df_train["label"] == 0).sum()),
            "n_val_binders":                    int((df_val["label"] == 1).sum()),
            "n_val_nonbinders":                 int((df_val["label"] == 0).sum()),
            "n_val_alleles":                    len(val_alleles),
            "n_peptide_clusters_selected":      len(sel_val_clust),
            "n_rows_excluded_by_peptide_cover": n_excluded,
            "n_rows_recycled_to_train":         n_recycled,
        })

    # ── Save summary outputs ───────────────────────────────────────────────────
    summary_path = splits_dir / "fold_summary.csv"
    pd.DataFrame(fold_summary).to_csv(summary_path, index=False)
    log.info(f"\nSaved: {summary_path}")

    meta = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    meta_path = splits_dir / "split_meta.json"
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=2)
    log.info(f"Saved: {meta_path}")

    log.info(f"\nDone.  Total: {time.time()-t_total:.1f}s")


if __name__ == "__main__":
    main()
