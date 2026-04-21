"""
filter_map_train_prep.py
========================
Pipeline for preparing ProteinMPNN fine-tuning training data from PMGen pMHC
structure predictions. Four sequential stages:

  1. FILTER   — select best structure per (allele, peptide) pair; apply pLDDT threshold
  2. MAP      — join back to raw parquet to recover metadata + MHC sequences
  3. RESAMPLE — re-run iterative median sampling to correct allele bias
  4. SPLIT    — produce train/val/test splits in two modes:
                  hla    : rare alleles (<20% by frequency) held out as independent test;
                           k-fold CV on remainder
                  anchor : no independent test set; k-fold CV stratified by anchor residues

Usage examples
--------------
# Binders, HLA-stratified split
python filter_map_train_prep.py \\
    --plddt_csv      .../plddt_means_binder.csv \\
    --parquet        .../PMDb_2025_11_18_class1.parquet \\
    --mhc_encodings  .../mhc1_encodings.csv \\
    --pdb_base_dir   .../outputs/binder \\
    --output_dir     .../trainprep/binder/ \\
    --mode           binder \\
    --plddt_threshold 80 \\
    --k              5 \\
    --split_mode     hla

# Non-binders, anchor-stratified split, keep both structures
python filter_map_train_prep.py \\
    --plddt_csv      .../plddt_means_nonbinder.csv \\
    --parquet        .../PMDb_2025_11_18_class1.parquet \\
    --mhc_encodings  .../mhc1_encodings.csv \\
    --pdb_base_dir   .../outputs/nonbinder \\
    --output_dir     .../trainprep/nonbinder/ \\
    --mode           nonbinder \\
    --keep_all \\
    --plddt_threshold 0 \\
    --k              5 \\
    --split_mode     anchor
"""

import argparse
import logging
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
_PMHC_SAMPLING_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(_PMHC_SAMPLING_DIR))
try:
    from pmhc_sampling import (
        plot_comparison,
        plot_per_allele_anchor_fractions,
        plot_anchor_kl_divergence,
        compute_anchor_combo_stats,
        compute_stats,
        save_stats,
    )
    _SAMPLING_PLOTS_AVAILABLE = True
except ImportError:
    _SAMPLING_PLOTS_AVAILABLE = False

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
log = logging.getLogger(__name__)


def setup_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))
    logging.getLogger().addHandler(file_handler)


# ─────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────
SAMPLE_CAP      = 1000
MAX_PEPTIDE_LEN = 15
VALID_AA        = set("ACDEFGHIKLMNPQRSTVWY")
AA_ORDER        = list("ACDEFGHIKLMNPQRSTVWY")
SMOOTH_EPS      = 1e-6

COL_PEPTIDE = "long_mer"
COL_MHC     = "allele"
COL_LABEL   = "assigned_label"
MHC_CLASS   = "mhc_class"

HLA_TEST_FRAC = 0.20


# ══════════════════════════════════════════════════════════════════
# UTILITIES
# ══════════════════════════════════════════════════════════════════

def clean_key(allele_key: str) -> str:
    """Normalise allele name: 'HLA-A*02:01' -> 'HLAA0201'"""
    if allele_key is None:
        return "None"
    mapping = str.maketrans({"*": "", ":": "", " ": "", "/": "_", "-": ""})
    return allele_key.translate(mapping).upper()


def add_anchor_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add a1_res (P2) and a2_res (last position) columns."""
    df = df.copy()
    df["a1_res"] = df[COL_PEPTIDE].str[1]
    df["a2_res"] = df[COL_PEPTIDE].str[-1]
    return df


def save_parquet(df: pd.DataFrame, path: Path, label: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    log.info(f"  Saved {label} -> {path}  ({len(df):,} rows)")


# ══════════════════════════════════════════════════════════════════
# STAGE 1 — FILTER
# ══════════════════════════════════════════════════════════════════

def filter_structures(
    plddt_csv: str,
    mode: str,
    keep_all: bool,
    plddt_threshold: float,
) -> pd.DataFrame:
    log.info("=" * 60)
    log.info("STAGE 1 — FILTER")
    log.info("=" * 60)

    df = pd.read_csv(plddt_csv)
    log.info(f"  Loaded pLDDT CSV: {len(df):,} rows")

    df["_score"] = df[["pep_mean_plddt", "anchor_mean_plddt"]].mean(axis=1)

    if mode == "binder" or not keep_all:
        best_idx = df.groupby(["allele", "peptide"])["_score"].idxmax()
        df = df.loc[best_idx].copy()
        log.info(f"  After best-structure selection: {len(df):,} rows")
    else:
        log.info(f"  --keep_all: retaining both structures per pair")

    df_all = df.drop(columns=["_score"], errors="ignore")

    if plddt_threshold > 0:
        df_filtered = df_all[df_all["pep_mean_plddt"] >= plddt_threshold].copy()
        log.info(
            f"  After pLDDT >= {plddt_threshold} filter: "
            f"{len(df_filtered):,} rows  (removed {len(df_all) - len(df_filtered):,})"
        )
    else:
        df_filtered = df_all.copy()
        log.info("  plddt_threshold=0: no filtering applied")

    log.info(f"  Unique (allele, peptide) pairs: "
             f"{df_filtered.groupby(['allele','peptide']).ngroups:,}")
    return df_all, df_filtered


# ══════════════════════════════════════════════════════════════════
# STAGE 2 — MAP
# ══════════════════════════════════════════════════════════════════

def map_to_parquet(
    filtered_df: pd.DataFrame,
    parquet_path: str,
    mhc_encodings_path: str,
    mode: str,
) -> pd.DataFrame:
    log.info("=" * 60)
    log.info("STAGE 2 — MAP")
    log.info("=" * 60)

    label_val = 1.0 if mode == "binder" else 0.0

    log.info(f"  Loading parquet: {parquet_path}")
    raw = pd.read_parquet(parquet_path)
    log.info(f"  Raw parquet rows: {len(raw):,}")

    raw = raw[raw[MHC_CLASS] == 1.0].copy()
    raw = raw[raw[COL_LABEL] == label_val].copy()
    log.info(f"  After class-I + {mode} filter: {len(raw):,} rows")

    raw["allele_clean"] = raw[COL_MHC].apply(clean_key)

    merge_keys = filtered_df[
        ["allele", "peptide", "chunk", "struct_idx", "pep_mean_plddt", "anchor_mean_plddt"]
    ].copy()

    merged = raw.merge(
        merge_keys,
        left_on=["allele_clean", COL_PEPTIDE],
        right_on=["allele", "peptide"],
        how="inner",
    )
    merged = merged.drop(columns=["allele_y", "peptide", "allele_clean"], errors="ignore")
    merged = merged.rename(columns={"allele_x": COL_MHC})
    merged[COL_MHC] = merged[COL_MHC].apply(clean_key)

    log.info(f"  After parquet join: {len(merged):,} rows")
    n_missing = len(filtered_df) - len(merged)
    if n_missing > 0:
        log.warning(f"  {n_missing:,} filtered rows had no match in parquet")

    log.info(f"  Loading MHC encodings: {mhc_encodings_path}")
    enc = pd.read_csv(mhc_encodings_path)
    enc["key_clean"] = enc["key"].apply(clean_key)

    merged["allele_clean"] = merged[COL_MHC].apply(clean_key)
    merged = merged.merge(
        enc[["key_clean", "mhc_sequence"]],
        left_on="allele_clean",
        right_on="key_clean",
        how="left",
    )
    merged = merged.drop(columns=["key_clean", "allele_clean"], errors="ignore")

    n_no_seq = merged["mhc_sequence"].isna().sum()
    if n_no_seq > 0:
        log.warning(f"  {n_no_seq:,} rows have no MHC sequence — will be dropped")
        merged = merged[merged["mhc_sequence"].notna()].copy()

    log.info(f"  Final mapped rows: {len(merged):,}")
    log.info(f"  Unique alleles   : {merged[COL_MHC].nunique():,}")
    log.info(f"  Unique peptides  : {merged[COL_PEPTIDE].nunique():,}")
    return merged


# ══════════════════════════════════════════════════════════════════
# STAGE 2b — ADD PDB PATHS
# ══════════════════════════════════════════════════════════════════

def add_pdb_paths(df: pd.DataFrame, pdb_base_dir: str) -> pd.DataFrame:
    if pdb_base_dir is None:
        log.info("  --pdb_base_dir not provided — skipping PDB path construction")
        return df

    log.info(f"  Constructing PDB paths from base dir: {pdb_base_dir}")

    df = df.copy()
    df["pdb_path"] = (
        pdb_base_dir.rstrip("/") + "/" +
        df["chunk"] + "/alphafold/" +
        df["allele"].apply(clean_key) + "_" +
        df[COL_PEPTIDE] + "_" +
        df["struct_idx"].astype(str)
    )

    sample = df["pdb_path"].head(5)
    missing_sample = [p for p in sample if not Path(p).exists()]
    if missing_sample:
        log.warning(f"  Sample check: {len(missing_sample)}/5 PDB directories not found — check pdb_base_dir")
        for p in missing_sample:
            log.warning(f"    Missing: {p}")
    else:
        log.info(f"  Sample check (5 files): all found ✓ — {len(df):,} paths constructed total")

    return df


# ══════════════════════════════════════════════════════════════════
# STAGE 3 — RESAMPLE
# ══════════════════════════════════════════════════════════════════

def _iterative_median_sample(
    df: pd.DataFrame,
    group_col: str,
    label: str = "",
) -> pd.DataFrame:
    pools   = {g: grp.copy() for g, grp in df.groupby(group_col)}
    sampled = {g: []          for g in pools}
    tallies = {g: 0           for g in pools}
    active  = set(pools.keys())

    iteration = 0
    while active:
        iteration += 1
        med = int(np.median([len(pools[g]) for g in active]))
        log.info(f"  [{label}] Iter {iteration:>3} | active: {len(active):>4} | median remaining: {med}")
        to_retire = set()
        for group in active:
            already  = tallies[group]
            n_needed = min(med, SAMPLE_CAP) - already
            if n_needed <= 0:
                to_retire.add(group); continue
            pool = pools[group]
            if len(pool) == 0:
                to_retire.add(group); continue
            n_draw = min(n_needed, len(pool))
            drawn  = pool.sample(n=n_draw, random_state=iteration)
            sampled[group].append(drawn)
            pools[group]    = pool.drop(drawn.index)
            tallies[group] += n_draw
            if len(pools[group]) == 0 or tallies[group] >= SAMPLE_CAP:
                to_retire.add(group)
        active -= to_retire
        log.info(f"             retired: {len(to_retire):>4}")

    all_chunks = [chunk for chunks in sampled.values() for chunk in chunks]
    return pd.concat(all_chunks).reset_index(drop=True)


def phase1_mhc_sampling(df: pd.DataFrame) -> pd.DataFrame:
    log.info("--- Phase 1: MHC allele balancing ---")
    result = _iterative_median_sample(df, group_col=COL_MHC, label="MHC")
    log.info(f"  Phase 1 done | {len(result):,} rows | {result[COL_MHC].nunique()} alleles")
    return result


def phase2_anchor_sampling(df: pd.DataFrame) -> pd.DataFrame:
    log.info("--- Phase 2: Anchor residue balancing (per allele) ---")
    df      = add_anchor_columns(df)
    alleles = df[COL_MHC].unique()
    parts   = []
    for i, allele in enumerate(alleles):
        sub = df[df[COL_MHC] == allele].copy()
        sub = _iterative_median_sample(sub, group_col="a1_res", label=f"{allele}/a1")
        sub = _iterative_median_sample(sub, group_col="a2_res", label=f"{allele}/a2")
        parts.append(sub)
        if (i + 1) % 50 == 0 or (i + 1) == len(alleles):
            log.info(f"  Processed {i+1}/{len(alleles)} alleles")
    result = pd.concat(parts).reset_index(drop=True)
    result = result.drop(columns=["a1_res", "a2_res"], errors="ignore")
    log.info(f"  Phase 2 done | {len(result):,} rows")
    return result


def resample(df: pd.DataFrame, mode: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    log.info("=" * 60)
    log.info("STAGE 3 — RESAMPLE")
    log.info("=" * 60)

    before = len(df)
    valid_mask = df[COL_PEPTIDE].apply(lambda p: all(c in VALID_AA for c in str(p)))
    df = df[valid_mask].copy()
    len_mask = df[COL_PEPTIDE].str.len() < MAX_PEPTIDE_LEN
    df = df[len_mask].copy()
    df = df.drop_duplicates(subset=[COL_PEPTIDE, COL_MHC]).copy()
    log.info(f"  After cleaning: {len(df):,} rows (removed {before - len(df):,})")

    df_p1 = phase1_mhc_sampling(df)

    if mode == "nonbinder":
        df_p2 = phase2_anchor_sampling(df_p1)
        return df_p1, df_p2
    else:
        return df_p1, None


# ══════════════════════════════════════════════════════════════════
# STAGE 4 — SPLIT
# ══════════════════════════════════════════════════════════════════

def split_hla(df: pd.DataFrame, k: int, out_dir: Path) -> None:
    log.info("--- Split mode: HLA ---")
    split_dir = out_dir / "splits" / "hla"
    split_dir.mkdir(parents=True, exist_ok=True)

    allele_counts  = df[COL_MHC].value_counts().sort_values(ascending=True)
    n_test_alleles = max(1, int(np.floor(len(allele_counts) * HLA_TEST_FRAC)))
    test_alleles   = set(allele_counts.index[:n_test_alleles])
    train_alleles  = [a for a in allele_counts.index if a not in test_alleles]

    log.info(f"  Total alleles      : {len(allele_counts):,}")
    log.info(f"  Test alleles       : {len(test_alleles):,}  (bottom {HLA_TEST_FRAC:.0%} by frequency)")
    log.info(f"  CV alleles         : {len(train_alleles):,}")

    df_test = df[df[COL_MHC].isin(test_alleles)].copy()
    df_cv   = df[df[COL_MHC].isin(train_alleles)].copy().reset_index(drop=True)

    save_parquet(df_test, split_dir / "test.parquet", "test set")
    pd.DataFrame({"allele": sorted(test_alleles)}).to_csv(split_dir / "test_alleles.csv", index=False)

    allele_arr   = np.array(train_alleles)
    rng          = np.random.default_rng(42)
    rng.shuffle(allele_arr)
    allele_folds = np.array_split(allele_arr, k)

    fold_summary = []
    for i in range(k):
        val_alleles        = set(allele_folds[i])
        train_alleles_fold = set(allele_arr) - val_alleles

        df_train_fold = df_cv[df_cv[COL_MHC].isin(train_alleles_fold)].copy()
        df_val_fold   = df_cv[df_cv[COL_MHC].isin(val_alleles)].copy()

        fold_dir = split_dir / f"fold_{i+1}"
        fold_dir.mkdir(parents=True, exist_ok=True)
        save_parquet(df_train_fold, fold_dir / "train.parquet", f"fold {i+1} train")
        save_parquet(df_val_fold,   fold_dir / "val.parquet",   f"fold {i+1} val")

        fold_summary.append({
            "fold":            i + 1,
            "n_train":         len(df_train_fold),
            "n_val":           len(df_val_fold),
            "n_train_alleles": len(train_alleles_fold),
            "n_val_alleles":   len(val_alleles),
        })

    pd.DataFrame(fold_summary).to_csv(split_dir / "fold_summary.csv", index=False)
    log.info(f"  Fold summary -> {split_dir}/fold_summary.csv")

    meta = {
        "split_mode":         "hla",
        "k":                  k,
        "hla_test_frac":      HLA_TEST_FRAC,
        "n_test_alleles":     len(test_alleles),
        "n_cv_alleles":       len(train_alleles),
        "n_test_rows":        len(df_test),
        "n_cv_rows":          len(df_cv),
        "pdb_paths_included": "pdb_path" in df.columns,
    }
    with open(split_dir / "split_meta.json", "w") as f:
        json.dump(meta, f, indent=2)


def split_anchor(df: pd.DataFrame, k: int, out_dir: Path) -> None:
    log.info("--- Split mode: Anchor ---")
    split_dir = out_dir / "splits" / "anchor"
    split_dir.mkdir(parents=True, exist_ok=True)

    df = add_anchor_columns(df).copy()
    df["_anchor_combo"] = df["a1_res"] + df["a2_res"]

    all_combos = df["_anchor_combo"].unique()
    rng        = np.random.default_rng(42)
    rng.shuffle(all_combos)

    log.info(f"  Total unique anchor combos : {len(all_combos)}")

    combo_chunks = np.array_split(all_combos, k)
    log.info(f"  Combos per chunk (~)       : {len(all_combos) // k}")

    drop_cols = ["_anchor_combo", "a1_res", "a2_res"]

    fold_summary = []
    for i in range(k):
        test_chunk   = set(combo_chunks[i])
        val_chunk    = set(combo_chunks[(i + 1) % k])
        train_combos = set(np.concatenate([
            combo_chunks[j] for j in range(k)
            if j != i and j != (i + 1) % k
        ]))

        df_test  = df[df["_anchor_combo"].isin(test_chunk)  ].drop(columns=drop_cols, errors="ignore").copy()
        df_val   = df[df["_anchor_combo"].isin(val_chunk)   ].drop(columns=drop_cols, errors="ignore").copy()
        df_train = df[df["_anchor_combo"].isin(train_combos)].drop(columns=drop_cols, errors="ignore").copy()

        fold_dir = split_dir / f"fold_{i+1}"
        fold_dir.mkdir(parents=True, exist_ok=True)
        save_parquet(df_test,  fold_dir / "test.parquet",  f"fold {i+1} test")
        save_parquet(df_val,   fold_dir / "val.parquet",   f"fold {i+1} val")
        save_parquet(df_train, fold_dir / "train.parquet", f"fold {i+1} train")

        pd.DataFrame({"anchor_combo": sorted(test_chunk)}).to_csv(
            fold_dir / "test_combos.csv", index=False
        )

        fold_summary.append({
            "fold":         i + 1,
            "test_combos":  len(test_chunk),
            "val_combos":   len(val_chunk),
            "train_combos": len(train_combos),
            "n_test":       len(df_test),
            "n_val":        len(df_val),
            "n_train":      len(df_train),
        })
        log.info(
            f"  Fold {i+1}: "
            f"test={len(test_chunk)} combos ({len(df_test):,} rows) | "
            f"val={len(val_chunk)} combos ({len(df_val):,} rows) | "
            f"train={len(train_combos)} combos ({len(df_train):,} rows)"
        )

    pd.DataFrame(fold_summary).to_csv(split_dir / "fold_summary.csv", index=False)
    log.info(f"  Fold summary -> {split_dir}/fold_summary.csv")
    log.info(f"  Total files  : {k * 3} (= {k} folds x 3 splits)")

    meta = {
        "split_mode":    "anchor",
        "k":             k,
        "n_rows":        len(df),
        "n_combos":      len(all_combos),
        "files_per_fold": 3,
        "total_files":   k * 3,
        "note": (
            "Unit of splitting is anchor combo (P2+last), not rows. "
            "Combos shuffled then split into k chunks. "
            "Rotating assignment: test=chunk_i, val=chunk_(i+1)%k, train=rest. "
            "Every combo is held out as test exactly once across k folds."
        ),
    }
    with open(split_dir / "split_meta.json", "w") as f:
        json.dump(meta, f, indent=2)


def run_splits(df: pd.DataFrame, split_mode: str, k: int, out_dir: Path) -> None:
    log.info("=" * 60)
    log.info("STAGE 4 — SPLIT")
    log.info("=" * 60)
    if split_mode in ("hla", "both"):
        split_hla(df, k=k, out_dir=out_dir)
    if split_mode in ("anchor", "both"):
        split_anchor(df, k=k, out_dir=out_dir)


# ══════════════════════════════════════════════════════════════════
# DIAGNOSTIC PLOTS
# ══════════════════════════════════════════════════════════════════

def plot_plddt_distribution(
    df_all: pd.DataFrame,
    df_filtered: pd.DataFrame,
    out_dir: Path,
    mode: str,
    plddt_threshold: float,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharey=False)

    for row, col in enumerate(["pep_mean_plddt", "anchor_mean_plddt"]):
        col_label = "Peptide" if col == "pep_mean_plddt" else "Anchor"
        for column, (df, stage) in enumerate([
            (df_all,      "Before filtering"),
            (df_filtered, "After filtering"),
        ]):
            ax = axes[row, column]
            ax.hist(df[col].dropna(), bins=50, color="steelblue", edgecolor="white", alpha=0.85)
            if plddt_threshold > 0:
                ax.axvline(plddt_threshold, color="red", linestyle="--",
                           linewidth=1.2, label=f"threshold={plddt_threshold}")
            ax.set_title(f"{col_label} mean pLDDT — {stage}\n{mode} | n={len(df):,}")
            ax.set_xlabel("mean pLDDT")
            ax.set_ylabel("Count")
            ax.legend()

    plt.suptitle(f"pLDDT distribution before vs after filtering  ({mode})", fontsize=12, y=1.01)
    plt.tight_layout()
    path = out_dir / "plddt_distribution_before_after.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"  Saved -> {path}")


def plot_allele_distribution(df: pd.DataFrame, out_dir: Path, label: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    counts = df[COL_MHC].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(range(len(counts)), counts.values, width=1.0, color="steelblue")
    ax.set_title(
        f"Peptide count per allele — {label}\n"
        f"n_alleles={len(counts):,}  median={int(counts.median()):,}  max={counts.max():,}"
    )
    ax.set_ylabel("Count")
    ax.set_xlabel("MHC allele (sorted by count)")
    ax.set_xticks([])
    plt.tight_layout()
    safe = label.replace(" ", "_").lower()
    path = out_dir / f"allele_distribution_{safe}.png"
    fig.savefig(path, dpi=150)
    plt.close()
    log.info(f"  Saved -> {path}")


# ══════════════════════════════════════════════════════════════════
# SUMMARY REPORT
# ══════════════════════════════════════════════════════════════════

def save_summary(
    args,
    df_filtered: pd.DataFrame,
    df_mapped: pd.DataFrame,
    df_resampled: pd.DataFrame,
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "pipeline_summary.txt"
    with open(path, "w") as f:
        f.write("filter_map_train_prep.py — Pipeline Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"mode             : {args.mode}\n")
        f.write(f"plddt_threshold  : {args.plddt_threshold}\n")
        f.write(f"keep_all         : {args.keep_all}\n")
        f.write(f"split_mode       : {args.split_mode}\n")
        f.write(f"k                : {args.k}\n\n")
        f.write(f"Input files\n")
        f.write(f"  plddt_csv      : {args.plddt_csv}\n")
        f.write(f"  parquet        : {args.parquet}\n")
        f.write(f"  mhc_encodings  : {args.mhc_encodings}\n")
        f.write(f"  pdb_base_dir   : {args.pdb_base_dir}\n\n")
        f.write(f"Stage results\n")
        f.write(f"  After filter   : {len(df_filtered):,} rows\n")
        f.write(f"  After map      : {len(df_mapped):,} rows\n")
        f.write(f"  After resample : {len(df_resampled):,} rows  "
                f"| {df_resampled[COL_MHC].nunique():,} alleles\n")
    log.info(f"  Pipeline summary -> {path}")


# ══════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Filter, map, resample and split PMGen pMHC structures for ProteinMPNN fine-tuning.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--plddt_csv", required=True)
    parser.add_argument("--parquet", required=True)
    parser.add_argument("--mhc_encodings", required=True)
    parser.add_argument("--pdb_base_dir", default=None,
        help="Base directory containing PMGen PDB outputs. "
             "Expected: {pdb_base_dir}/{chunk}/alphafold/{allele}_{peptide}_{struct_idx}/")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--mode", required=True, choices=["binder", "nonbinder"])
    parser.add_argument("--keep_all", action="store_true",
        help="(nonbinder only) Keep both structures per pair")
    parser.add_argument("--plddt_threshold", type=float, default=80.0,
        help="Minimum pep_mean_plddt. Set 0 to disable.")
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--split_mode", required=True, choices=["hla", "anchor", "both"])

    args = parser.parse_args()

    out_dir  = Path(args.output_dir)
    plot_dir = out_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(out_dir / "pipeline.log")

    # Stage 1: Filter
    df_all_structures, df_filtered = filter_structures(
        plddt_csv=args.plddt_csv, mode=args.mode,
        keep_all=args.keep_all, plddt_threshold=args.plddt_threshold,
    )
    plot_plddt_distribution(df_all_structures, df_filtered, plot_dir, args.mode, args.plddt_threshold)

    # Stage 2: Map
    df_mapped = map_to_parquet(
        filtered_df=df_filtered, parquet_path=args.parquet,
        mhc_encodings_path=args.mhc_encodings, mode=args.mode,
    )
    df_mapped = add_pdb_paths(df_mapped, args.pdb_base_dir)
    plot_allele_distribution(df_mapped, plot_dir, "after mapping")
    save_parquet(df_mapped, out_dir / "mapped.parquet", "mapped dataset")

    # Stage 3: Resample
    df_p1, df_p2 = resample(df_mapped, mode=args.mode)
    df_resampled  = df_p2 if df_p2 is not None else df_p1
    plot_allele_distribution(df_resampled, plot_dir, "after resampling")
    save_parquet(df_resampled, out_dir / "resampled.parquet", "resampled dataset")

    # Diagnostic plots
    if _SAMPLING_PLOTS_AVAILABLE:
        log.info("=" * 60)
        log.info("DIAGNOSTIC PLOTS  (pmhc_sampling style)")
        log.info("=" * 60)
        mode_label = "Binders" if args.mode == "binder" else "Non-binders"

        if df_p2 is not None:
            stages    = [(f"Mapped ({mode_label})", df_mapped), ("Post-Phase 1", df_p1), ("Post-Phase 2", df_p2)]
            kl_stages = [("Post-Phase 1", df_p1), ("Post-Phase 2", df_p2)]
        else:
            stages    = [(f"Mapped ({mode_label})", df_mapped), ("Post-Phase 1", df_p1)]
            kl_stages = [("Post-Phase 1", df_p1)]

        plot_comparison(stages=stages, out_dir=plot_dir)
        plot_per_allele_anchor_fractions(df_p1, "Post-Phase 1", plot_dir)
        if df_p2 is not None:
            plot_per_allele_anchor_fractions(df_p2, "Post-Phase 2", plot_dir)
        plot_anchor_kl_divergence(stages=kl_stages, out_dir=plot_dir)
        compute_anchor_combo_stats(stages=stages, out_dir=plot_dir)
        stats_list = [compute_stats(df, label) for label, df in stages]
        save_stats(stats_list, out_dir=plot_dir)
    else:
        log.warning(f"pmhc_sampling.py not found at {_PMHC_SAMPLING_DIR / 'pmhc_sampling.py'}")

    # Stage 4: Split
    run_splits(df_resampled, split_mode=args.split_mode, k=args.k, out_dir=out_dir)

    # Summary
    save_summary(args, df_filtered, df_mapped, df_resampled, out_dir)

    log.info("=" * 60)
    log.info("DONE — all outputs saved to: %s", out_dir)
    log.info("=" * 60)


if __name__ == "__main__":
    main()