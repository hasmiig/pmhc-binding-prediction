"""
pMHC Non-Binder Sampling Pipeline
===================================
Two-phase iterative median sampling strategy:
  - Phase 1 : balance across MHC alleles          (cap: SAMPLE_CAP per allele)
  - Phase 2 : balance anchor residue diversity     (per allele, independently for a1 and a2)

Usage:
    python pmhc_sampling.py --input data.parquet --inspect
    python pmhc_sampling.py --input data.parquet --explore  --plots /path/to/plots
    python pmhc_sampling.py --input data.parquet --output out.parquet --plots /path/to/plots
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.special import kl_div  # noqa: F401 – imported for reference; KL computed manually

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────
SAMPLE_CAP = 1000

ANCHOR_RULES = {
    8:  (1, 7),
    9:  (1, 8),
    10: (1, 9),
    11: (1, 10),
    12: (1, 11),
    13: (1, 12),
    14: (1, 13),
}

VALID_AA   = set("ACDEFGHIKLMNPQRSTVWY")
AA_ORDER   = list("ACDEFGHIKLMNPQRSTVWY")   # alphabetical, used in new plots
CHUNK_SIZE = 50                              # alleles per page in per-allele plots
SMOOTH_EPS = 1e-6                           # smoothing for KL divergence

def _aa_palette() -> dict:
    """20 visually distinct colours, one per amino acid."""
    cmap_a = plt.get_cmap("tab20")
    cmap_b = plt.get_cmap("tab20b")
    colours = [cmap_a(i / 20) for i in range(20)]
    colours[18] = cmap_b(0 / 20)
    colours[19] = cmap_b(4 / 20)
    return dict(zip(AA_ORDER, colours))

AA_COLOURS = _aa_palette()

MAX_PEPTIDE_LEN = 15   # peptides of length >= this are removed due to low sampling survival

# ─────────────────────────────────────────────────────────────────
# COLUMN NAMES
# ─────────────────────────────────────────────────────────────────
COL_PEPTIDE = "long_mer"
COL_MHC     = "allele"
COL_LABEL   = "assigned_label"
MHC_CLASS   = "mhc_class"


# ══════════════════════════════════════════════════════════════════
# 0.  I/O
# ══════════════════════════════════════════════════════════════════

def inspect_file(path: str) -> None:
    ext = Path(path).suffix.lower()
    log.info(f"=== Inspecting: {path} ===")
    if ext == ".parquet":
        try:
            import pyarrow.parquet as pq
            print("\n── Schema ──")
            for field in pq.read_schema(path):
                print(f"  {field.name:<35} {field.type}")
        except ImportError:
            log.warning("pyarrow not installed, skipping schema.")
    df = _read_file(path, nrows=5)
    print("\n── Dtypes ──")
    print(df.dtypes.to_string())
    print("\n── First 5 rows ──")
    print(df.to_string())
    print("\n── Column names ──")
    print(df.columns.tolist())


def _read_file(path: str, nrows: int = None) -> pd.DataFrame:
    ext = Path(path).suffix.lower()
    if ext == ".parquet":
        df = pd.read_parquet(path)
        return df.head(nrows) if nrows else df
    elif ext in (".tsv", ".txt"):
        return pd.read_csv(path, sep="\t", nrows=nrows)
    elif ext == ".csv":
        return pd.read_csv(path, nrows=nrows)
    raise ValueError(f"Unsupported format: '{ext}'")


def clean_key(allele_key: str) -> str:
    """Normalise allele name: 'HLA-A*02:01' -> 'HLAA0201'"""
    if allele_key is None:
        return "None"
    mapping = str.maketrans({"*": "", ":": "", " ": "", "/": "_", "-": ""})
    return allele_key.translate(mapping).upper()


def get_seq(key: str, seq_dict: dict):
    """Look up MHC sequence from dict by partial key match."""
    return next(
        (seq for seq_key, seq in seq_dict.items()
         if seq_key.upper().startswith(key.upper())),
        None,
    )


def load_data(path: str, mode: str = "nonbinder", phases: str = "both") -> pd.DataFrame:
    """Load file, keep MHC class I binders or non-binders, deduplicate, normalise allele names."""
    assert mode in ("binder", "nonbinder"), f"mode must be 'binder' or 'nonbinder', got '{mode}'"
    assert phases in ("only_phase1", "both"), f"phases must be 'only_phase1' or 'both', got '{phases}'"
    label_val  = 1.0 if mode == "binder" else 0.0
    label_name = "Binder" if mode == "binder" else "Non-binder"

    log.info(f"Loading data from {path} (mode={mode}, phases={phases}) ...")
    df = _read_file(path)
    log.info(f"  Total rows      : {len(df):,}")
    log.info(f"  Columns         : {list(df.columns)}")

    df = df[df[MHC_CLASS] == 1.0].copy()
    log.info(f"  Class I rows    : {len(df):,}")

    df = df[df[COL_LABEL] == label_val].copy()
    log.info(f"  {label_name} rows  : {len(df):,}")

    df = df.drop_duplicates(subset=[COL_PEPTIDE, COL_MHC]).copy()
    log.info(f"  After dedup     : {len(df):,}")

    # remove peptides with non-standard characters
    # catches: accession numbers (e.g. 48761), placeholders (XXXX), rare AAs (U, B, Z)
    valid_mask = df[COL_PEPTIDE].apply(lambda p: all(c in VALID_AA for c in str(p)))
    n_invalid  = (~valid_mask).sum()
    if n_invalid > 0:
        log.warning(f"  Removing {n_invalid:,} peptides with non-standard characters")
        df = df[valid_mask].copy()
    log.info(f"  After AA filter : {len(df):,}")

    # remove peptides of length 15 or longer — too few survive sampling to be useful
    len_mask  = df[COL_PEPTIDE].str.len() < MAX_PEPTIDE_LEN
    n_toolong = (~len_mask).sum()
    if n_toolong > 0:
        log.warning(f"  Removing {n_toolong:,} peptides with length >= {MAX_PEPTIDE_LEN}")
        df = df[len_mask].copy()
    log.info(f"  After len filter: {len(df):,}")

    df[COL_MHC] = df[COL_MHC].apply(clean_key)
    log.info(f"  Allele names normalised")
    return df


def save_data(df: pd.DataFrame, path: str) -> None:
    ext = Path(path).suffix.lower()
    if ext == ".parquet":
        df.to_parquet(path, index=False)
    elif ext in (".tsv", ".txt"):
        df.to_csv(path, sep="\t", index=False)
    elif ext == ".csv":
        df.to_csv(path, index=False)
    else:
        raise ValueError(f"Unsupported format: '{ext}'")
    log.info(f"Saved {len(df):,} rows -> {path}")


# ══════════════════════════════════════════════════════════════════
# 1.  ANCHOR HELPERS
# ══════════════════════════════════════════════════════════════════

def add_anchor_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add a1_res (P2) and a2_res (last position) columns.
    Uses fully vectorised string ops -- fast even on 40M+ rows.
    Anchor rule: position 2 (index 1) and last position for all peptide lengths.
    """
    df = df.copy()
    df["a1_res"] = df[COL_PEPTIDE].str[1]    # P2 — always index 1
    df["a2_res"] = df[COL_PEPTIDE].str[-1]   # last position
    return df


# ══════════════════════════════════════════════════════════════════
# 2.  STATISTICS  (saved to file)
# ══════════════════════════════════════════════════════════════════

def compute_stats(df: pd.DataFrame, label: str) -> dict:
    """Compute summary statistics. Returns a dict."""
    counts = df[COL_MHC].value_counts()
    df     = df.copy()
    df["pep_len"] = df[COL_PEPTIDE].str.len()
    df = add_anchor_columns(df)

    return {
        "label":                label,
        "total_rows":           len(df),
        "unique_alleles":       len(counts),
        "rows_per_allele_min":  int(counts.min()),
        "rows_per_allele_25p":  int(counts.quantile(0.25)),
        "rows_per_allele_med":  int(counts.median()),
        "rows_per_allele_75p":  int(counts.quantile(0.75)),
        "rows_per_allele_max":  int(counts.max()),
        "alleles_at_cap":       int((counts >= SAMPLE_CAP).sum()),
        "alleles_below_cap":    int((counts < SAMPLE_CAP).sum()),
        "peptide_length_dist":  df["pep_len"].value_counts().sort_index().to_dict(),
        "anchor1_freq":         df["a1_res"].value_counts().sort_index().to_dict(),
        "anchor2_freq":         df["a2_res"].value_counts().sort_index().to_dict(),
    }


def save_stats(stats_list: list, out_dir: Path) -> None:
    """Save all stage statistics to a text report and a CSV summary."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # text report
    txt_path = out_dir / "stats_report.txt"
    with open(txt_path, "w") as f:
        for s in stats_list:
            f.write(f"\n{'='*60}\n")
            f.write(f"  {s['label']}\n")
            f.write(f"{'='*60}\n")
            f.write(f"  Total rows              : {s['total_rows']:,}\n")
            f.write(f"  Unique alleles          : {s['unique_alleles']:,}\n")
            f.write(f"  Rows/allele  min        : {s['rows_per_allele_min']:,}\n")
            f.write(f"  Rows/allele  25pct      : {s['rows_per_allele_25p']:,}\n")
            f.write(f"  Rows/allele  median     : {s['rows_per_allele_med']:,}\n")
            f.write(f"  Rows/allele  75pct      : {s['rows_per_allele_75p']:,}\n")
            f.write(f"  Rows/allele  max        : {s['rows_per_allele_max']:,}\n")
            f.write(f"  Alleles >= cap ({SAMPLE_CAP:,})  : {s['alleles_at_cap']:,}\n")
            f.write(f"  Alleles <  cap           : {s['alleles_below_cap']:,}\n")
            f.write(f"  Peptide length dist      : {s['peptide_length_dist']}\n")
            f.write(f"  Anchor 1 (P2) freq       : {s['anchor1_freq']}\n")
            f.write(f"  Anchor 2 (last) freq     : {s['anchor2_freq']}\n")
    log.info(f"Stats report -> {txt_path}")

    # CSV summary (numeric columns only, one row per stage)
    csv_rows = [{
        "stage":             s["label"],
        "total_rows":        s["total_rows"],
        "unique_alleles":    s["unique_alleles"],
        "min_per_allele":    s["rows_per_allele_min"],
        "p25_per_allele":    s["rows_per_allele_25p"],
        "median_per_allele": s["rows_per_allele_med"],
        "p75_per_allele":    s["rows_per_allele_75p"],
        "max_per_allele":    s["rows_per_allele_max"],
        "alleles_at_cap":    s["alleles_at_cap"],
        "alleles_below_cap": s["alleles_below_cap"],
    } for s in stats_list]
    pd.DataFrame(csv_rows).to_csv(out_dir / "stats_summary.csv", index=False)
    log.info(f"Stats CSV     -> {out_dir}/stats_summary.csv")


# ══════════════════════════════════════════════════════════════════
# 3.  EXPLORATION PLOTS  (pre-sampling)
# ══════════════════════════════════════════════════════════════════

def explore(df: pd.DataFrame, out_dir: Path) -> None:
    """Generate pre-sampling exploration plots."""
    out_dir.mkdir(parents=True, exist_ok=True)
    df = df.copy()
    df["pep_len"] = df[COL_PEPTIDE].str.len()
    df = add_anchor_columns(df)   # fast — vectorised str ops

    # 1. allele distribution
    allele_counts = df[COL_MHC].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(range(len(allele_counts)), allele_counts.values, width=1.0, color="steelblue")
    ax.set_title(f"Peptide count per MHC allele\n"
                 f"n_alleles={len(allele_counts):,}  median={int(allele_counts.median()):,}  "
                 f"max={allele_counts.max():,}")
    ax.set_ylabel("Count")
    ax.set_xlabel("MHC allele (sorted by count)")
    ax.set_xticks([])
    plt.tight_layout()
    fig.savefig(out_dir / "00_raw_allele_distribution.png", dpi=150)
    plt.close()

    # 2. top 20 / bottom 20 alleles
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    for ax, subset, title in zip(
        axes,
        [allele_counts.tail(20), allele_counts.head(20)],
        ["Top 20 most frequent alleles", "Top 20 least frequent alleles"],
    ):
        subset.plot(kind="barh", ax=ax, color="steelblue")
        ax.set_title(title)
        ax.set_xlabel("Count")
    plt.tight_layout()
    fig.savefig(out_dir / "00_raw_top_bottom_alleles.png", dpi=150)
    plt.close()

    # 3. peptide length distribution (global)
    len_counts = df["pep_len"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(8, 4))
    len_counts.plot(kind="bar", ax=ax, color="coral", edgecolor="white")
    ax.set_title("Peptide length distribution")
    ax.set_xlabel("Length (aa)")
    ax.set_ylabel("Count")
    for bar in ax.patches:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.01,
                f"{int(bar.get_height()):,}", ha="center", va="bottom", fontsize=7, rotation=45)
    plt.tight_layout()
    fig.savefig(out_dir / "00_raw_peptide_lengths.png", dpi=150)
    plt.close()

    # 4. anchor residues (global, both positions)
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    for ax, col, title in zip(
        axes,
        ["a1_res", "a2_res"],
        ["Anchor 1 (P2) residue frequency",
         "Anchor 2 (last pos) residue frequency"],
    ):
        counts = df[col].value_counts().sort_index()
        counts.plot(kind="bar", ax=ax, color="mediumseagreen", edgecolor="white")
        ax.set_title(title)
        ax.set_xlabel("Amino acid")
        ax.set_ylabel("Count")
    plt.tight_layout()
    fig.savefig(out_dir / "00_raw_anchor_residues.png", dpi=150)
    plt.close()

    # 5. source distribution
    if "source" in df.columns:
        source_counts = df["source"].value_counts()
        fig, ax = plt.subplots(figsize=(8, 4))
        source_counts.plot(kind="barh", ax=ax, color="mediumpurple")
        ax.set_title("Data source distribution")
        ax.set_xlabel("Count")
        plt.tight_layout()
        fig.savefig(out_dir / "00_raw_source_distribution.png", dpi=150)
        plt.close()

    log.info(f"Exploration plots saved to {out_dir}/")


# ══════════════════════════════════════════════════════════════════
# 4.  CORE -- Generic iterative median sampler
# ══════════════════════════════════════════════════════════════════

def _iterative_median_sample(df: pd.DataFrame, group_col: str, label: str = "") -> pd.DataFrame:
    """
    Iterative median sampling over any grouping column.
    Each group has its own independent pool so sampling one group
    never affects another group's available rows.
    """
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
                to_retire.add(group)
                continue

            pool = pools[group]
            if len(pool) == 0:
                to_retire.add(group)
                continue

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


# ══════════════════════════════════════════════════════════════════
# 5.  PHASE 1 -- MHC allele balancing
# ══════════════════════════════════════════════════════════════════

def phase1_mhc_sampling(df: pd.DataFrame) -> pd.DataFrame:
    """Balance dataset across MHC alleles using iterative median sampling."""
    log.info("=== Phase 1: MHC allele sampling ===")
    result = _iterative_median_sample(df, group_col=COL_MHC, label="MHC")
    log.info(f"Phase 1 done | {len(result):,} rows | {result[COL_MHC].nunique()} alleles")
    return result


# ══════════════════════════════════════════════════════════════════
# 6.  PHASE 2 -- Anchor residue diversity (per allele)
# ══════════════════════════════════════════════════════════════════

def phase2_anchor_sampling(df: pd.DataFrame) -> pd.DataFrame:
    """
    Balance anchor residue diversity within each MHC allele independently.
    For each allele: apply iterative median sampling on a1_res, then a2_res.
    """
    log.info("=== Phase 2: Anchor residue sampling (per allele) ===")
    df      = add_anchor_columns(df)
    alleles = df[COL_MHC].unique()
    parts   = []

    for i, allele in enumerate(alleles):
        allele_df = df[df[COL_MHC] == allele].copy()
        allele_df = _iterative_median_sample(allele_df, group_col="a1_res", label=f"{allele}/a1")
        allele_df = _iterative_median_sample(allele_df, group_col="a2_res", label=f"{allele}/a2")
        parts.append(allele_df)
        if (i + 1) % 50 == 0 or (i + 1) == len(alleles):
            log.info(f"  Processed {i+1}/{len(alleles)} alleles")

    result = pd.concat(parts).reset_index(drop=True)
    result = result.drop(columns=["a1_idx", "a2_idx", "a1_res", "a2_res", "pep_len"],
                         errors="ignore")
    log.info(f"Phase 2 done | {len(result):,} rows | avg {len(result)/len(alleles):.0f} rows/allele")
    return result


# ══════════════════════════════════════════════════════════════════
# 7.  COMPARISON PLOTS  (raw vs post-phase1 vs post-phase2)
# ══════════════════════════════════════════════════════════════════

def plot_comparison(stages: list, out_dir: Path) -> None:
    """
    Plot allele distribution, peptide lengths, anchor 1 and anchor 2 residues
    for each stage side by side.
    stages = [("Raw", df_raw), ("Post-Phase 1", df_p1), ("Post-Phase 2", df_p2)]
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    n = len(stages)

    # add helper columns to each stage
    stages_prep = []
    for label, df in stages:
        df = df.copy()
        df["pep_len"] = df[COL_PEPTIDE].str.len()
        df = add_anchor_columns(df)
        stages_prep.append((label, df))

    def _save(fig, name):
        fig.savefig(out_dir / name, dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info(f"  Saved -> {name}")

    # allele distribution
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 4))
    for ax, (label, df) in zip(axes, stages_prep):
        counts = df[COL_MHC].value_counts().sort_values()
        ax.bar(range(len(counts)), counts.values, width=1.0, color="steelblue")
        ax.set_title(f"{label}\nn={len(df):,}  median={int(counts.median()):,}  max={counts.max():,}")
        ax.set_ylabel("Count")
        ax.set_xlabel(f"MHC allele (n={len(counts)})")
        ax.set_xticks([])
    plt.tight_layout()
    _save(fig, "comparison_allele_distribution.png")

    # peptide length distribution (global)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 4))
    for ax, (label, df) in zip(axes, stages_prep):
        df["pep_len"].value_counts().sort_index().plot(
            kind="bar", ax=ax, color="coral", edgecolor="white")
        ax.set_title(f"Peptide lengths -- {label}")
        ax.set_xlabel("Length (aa)")
        ax.set_ylabel("Count")
        for bar in ax.patches:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.01,
                    f"{int(bar.get_height()):,}", ha="center", va="bottom",
                    fontsize=6, rotation=45)
    plt.tight_layout()
    _save(fig, "comparison_peptide_lengths.png")

    # anchor 1 residues (global)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 4))
    for ax, (label, df) in zip(axes, stages_prep):
        df["a1_res"].value_counts().sort_index().plot(
            kind="bar", ax=ax, color="mediumseagreen", edgecolor="white")
        ax.set_title(f"Anchor 1 (P2) -- {label}")
        ax.set_xlabel("Amino acid")
        ax.set_ylabel("Count")
    plt.tight_layout()
    _save(fig, "comparison_anchor1_residues.png")

    # anchor 2 residues (global)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 4))
    for ax, (label, df) in zip(axes, stages_prep):
        df["a2_res"].value_counts().sort_index().plot(
            kind="bar", ax=ax, color="steelblue", edgecolor="white")
        ax.set_title(f"Anchor 2 (last pos) -- {label}")
        ax.set_xlabel("Amino acid")
        ax.set_ylabel("Count")
    plt.tight_layout()
    _save(fig, "comparison_anchor2_residues.png")

    log.info(f"All comparison plots saved to {out_dir}/")


# ══════════════════════════════════════════════════════════════════
# 8.  ANCHOR COMBINATION RATIO
# ══════════════════════════════════════════════════════════════════

def _anchor_combo_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a 20x20 matrix of anchor combination counts.
    Rows = anchor 1 (P2), Columns = anchor 2 (last position).
    All 20 standard amino acids are always present as rows/columns.
    """
    df = add_anchor_columns(df)
    matrix = pd.crosstab(df["a1_res"], df["a2_res"])
    # ensure all 20 AAs present even if some combos are missing
    matrix = matrix.reindex(index=AA_ORDER, columns=AA_ORDER, fill_value=0)
    return matrix

def compute_anchor_combo_stats(stages: list, out_dir: Path) -> None:
    """
    Anchor combination (a1+a2) analysis.
    stages = [("Raw", df_raw), ("Post-Phase 1", df_p1)]         # phase1-only -> 2 panels
    stages = [("Raw", df_raw), ("Post-Phase 1", df_p1), (...)]  # both phases -> 3 panels
    """
    n  = len(stages)
 
    # ── 1. N-panel 20x20 heatmap ──
    fig, axes = plt.subplots(1, n, figsize=(22 * n // 3, 7))
    if n == 1:
        axes = [axes]
 
    for ax, (label, df) in zip(axes, stages):
        matrix = _anchor_combo_matrix(df)
        im = ax.imshow(matrix.values, cmap="YlOrRd", aspect="auto")
        ax.set_xticks(range(20))
        ax.set_yticks(range(20))
        ax.set_xticklabels(AA_ORDER, fontsize=8)
        ax.set_yticklabels(AA_ORDER, fontsize=8)
        ax.set_xlabel("Anchor 2 (last position)")
        ax.set_ylabel("Anchor 1 (P2)")
        n_combos = int((matrix > 0).values.sum())
        ax.set_title(f"Anchor combo counts -- {label}\nn={len(df):,}  unique combos={n_combos}/400")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
 
    plt.suptitle("20x20 Anchor Combination Matrix (P2 x last position)", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(out_dir / "anchor_combo_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"  Saved -> anchor_combo_heatmap.png")
 
    # ── 2. per-allele combo stats CSV (last stage only) ──
    _, df_last = stages[-1]
    df_last_anc = add_anchor_columns(df_last.copy())
    df_last_anc["anchor_combo"] = df_last_anc["a1_res"] + df_last_anc["a2_res"]
 
    rows = []
    for allele, grp in df_last_anc.groupby(COL_MHC):
        combo_counts = grp["anchor_combo"].value_counts()
        top_combo    = combo_counts.index[0]
        top_ratio    = combo_counts.iloc[0] / len(grp)
        rows.append({
            "allele":          allele,
            "n_peptides":      len(grp),
            "n_unique_combos": grp["anchor_combo"].nunique(),
            "top_combo":       top_combo,
            "top_combo_ratio": round(top_ratio, 4),
        })
 
    combo_df = pd.DataFrame(rows).sort_values("top_combo_ratio", ascending=False)
    combo_df.to_csv(out_dir / "anchor_combo_stats.csv", index=False)
    log.info(f"  Anchor combo stats -> {out_dir}/anchor_combo_stats.csv")
 
    # ── 3. log summary ──
    log.info(f"  Avg unique combos/allele         : {combo_df['n_unique_combos'].mean():.1f}")
    log.info(f"  Alleles with only 1 combo        : {(combo_df['n_unique_combos'] == 1).sum()}")
    log.info(f"  Median top-combo dominance ratio : {combo_df['top_combo_ratio'].median():.2f}")
    log.info("  Most dominated alleles (least diversity):")
    log.info(combo_df[["allele", "n_peptides", "n_unique_combos",
                        "top_combo", "top_combo_ratio"]].head(10).to_string(index=False))
    log.info("  Most diverse alleles:")
    log.info(combo_df.sort_values("n_unique_combos", ascending=False)[
        ["allele", "n_peptides", "n_unique_combos", "top_combo", "top_combo_ratio"]
    ].head(10).to_string(index=False))
 
 

# ══════════════════════════════════════════════════════════════════
# 9.  PER-ALLELE ANCHOR FRACTION PLOTS
# ══════════════════════════════════════════════════════════════════

def _anchor_fraction_matrix(df: pd.DataFrame, anchor_col: str) -> pd.DataFrame:
    """
    Returns a (n_alleles × 20) DataFrame of per-allele AA fractions for one anchor.
    Rows = alleles sorted alphabetically, Columns = AA_ORDER.
    """
    df = add_anchor_columns(df)
    grouped = (
        df.groupby([COL_MHC, anchor_col])
        .size()
        .rename("count")
        .reset_index()
    )
    pivot = grouped.pivot(index=COL_MHC, columns=anchor_col, values="count").fillna(0)
    for aa in AA_ORDER:
        if aa not in pivot.columns:
            pivot[aa] = 0.0
    pivot = pivot[AA_ORDER]
    pivot = pivot.div(pivot.sum(axis=1), axis=0)   # row-normalise to fractions
    pivot = pivot.sort_index()                       # alphabetical allele order
    return pivot


def _stacked_bar_page(ax, frac_df: pd.DataFrame, alleles: list, title: str) -> None:
    """Draw one stacked-bar page onto ax."""
    sub    = frac_df.loc[alleles]
    x      = np.arange(len(alleles))
    bottom = np.zeros(len(alleles))
    for aa in AA_ORDER:
        vals = sub[aa].values
        ax.bar(x, vals, bottom=bottom, color=AA_COLOURS[aa], width=0.85, label=aa)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(alleles, rotation=90, fontsize=6)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of peptides")
    ax.set_title(title, fontsize=9)


def plot_per_allele_anchor_fractions(df: pd.DataFrame,
                                      stage_label: str,
                                      out_dir: Path) -> None:
    """
    Produces exactly TWO output files per call (one per anchor position):
        per_allele_anchor1_<safe_stage>.png
        per_allele_anchor2_<safe_stage>.png

    Each file is a single figure with multiple subplots — one subplot per
    chunk of CHUNK_SIZE alleles — so all alleles are visible in one image.
    Call once for Post-Phase 1 and once for Post-Phase 2 to get 4 files total.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    safe_stage = stage_label.replace(" ", "_").replace("-", "").lower()  # e.g. postphase1

    legend_handles = [
        mpatches.Patch(color=AA_COLOURS[aa], label=aa) for aa in AA_ORDER
    ]

    for anchor_col, anchor_name, anchor_tag in [
        ("a1_res", "Anchor 1 (P2)",       "anchor1"),
        ("a2_res", "Anchor 2 (last pos)", "anchor2"),
    ]:
        frac_df = _anchor_fraction_matrix(df, anchor_col)
        alleles = frac_df.index.tolist()   # sorted alphabetically
        chunks  = [alleles[i:i + CHUNK_SIZE]
                   for i in range(0, len(alleles), CHUNK_SIZE)]
        n_chunks = len(chunks)

        # one tall figure: n_chunks rows, 1 column
        fig, axes = plt.subplots(
            n_chunks, 1,
            figsize=(max(16, CHUNK_SIZE * 0.30), 5 * n_chunks),
        )
        if n_chunks == 1:
            axes = [axes]   # always iterable

        for ax, chunk, part_idx in zip(axes, chunks, range(1, n_chunks + 1)):
            title = (
                f"Alleles {(part_idx - 1) * CHUNK_SIZE + 1}–"
                f"{(part_idx - 1) * CHUNK_SIZE + len(chunk)}"
                f" of {len(alleles)}"
            )
            _stacked_bar_page(ax, frac_df, chunk, title)

        # single shared legend on the last subplot
        axes[-1].legend(
            handles=legend_handles, title="AA",
            bbox_to_anchor=(1.01, 1), loc="upper left",
            fontsize=7, ncol=1, frameon=False,
        )

        fig.suptitle(
            f"{anchor_name} residue distribution per MHC allele — {stage_label}\n"
            f"({len(alleles)} alleles total, {CHUNK_SIZE} per panel)",
            fontsize=11, y=1.002,
        )
        plt.tight_layout()
        fname = f"per_allele_{anchor_tag}_{safe_stage}.png"
        fig.savefig(out_dir / fname, dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info(f"  Saved -> {fname}")

    log.info(f"Per-allele anchor fraction plots done ({stage_label})")



# ══════════════════════════════════════════════════════════════════
# 10. ANCHOR KL DIVERGENCE BOXPLOT
# ══════════════════════════════════════════════════════════════════

def _kl_from_uniform(count_vec: np.ndarray) -> float:
    """
    KL( p || Uniform(1/20) ) using Laplace-smoothed count vector.
    Returns KL in nats.
    """
    p = count_vec.astype(float) + SMOOTH_EPS
    p /= p.sum()
    q  = np.full(20, 1.0 / 20)
    return float(np.sum(p * np.log(p / q)))


def _per_allele_kl(df: pd.DataFrame, anchor_col: str) -> pd.Series:
    """Per-allele KL divergence of anchor distribution vs Uniform."""
    df = add_anchor_columns(df)
    kl_vals = {}
    for allele, grp in df.groupby(COL_MHC):
        counts = grp[anchor_col].value_counts()
        vec    = np.array([counts.get(aa, 0) for aa in AA_ORDER], dtype=float)
        kl_vals[allele] = _kl_from_uniform(vec)
    return pd.Series(kl_vals)


def _kl_baselines() -> list:
    """
    Compute KL( p || Uniform(1/20) ) for canonical 'dominance' distributions:
      - all mass on 1 AA  → most dominated possible
      - all mass on 2 AAs → evenly split across 2
      - all mass on 4, 8, 16 AAs
    Returns list of (n_aa, kl_value, label) sorted descending by KL.
    """
    baselines = []
    for k in [1, 2, 4, 8, 16]:
        vec = np.zeros(20)
        vec[:k] = 1.0 / k          # equal mass on k amino acids, zero elsewhere
        kl = _kl_from_uniform(vec)
        baselines.append((k, kl, f"top-{k} AA  ({kl:.2f})"))
    return baselines   # already ascending in k → descending in KL


def _draw_kl_broken_axis(ax_top, ax_bot,
                          anchor_name: str, stage_labels: list,
                          all_vals_list: list) -> None:
    """
    Draw broken-axis boxplots + jittered strips for multiple stages onto
    a (ax_top, ax_bot) pair.  The y-axis break is computed from the
    pooled whisker maximum so both stages share the same break point.
    Horizontal reference lines show canonical KL baselines.

    ax_bot  — zoomed into the IQR/whisker region
    ax_top  — outliers only
    """
    # compute shared break point from pooled data
    all_pooled = np.concatenate(all_vals_list)
    q3_pool    = np.percentile(all_pooled, 75)
    iqr_pool   = q3_pool - np.percentile(all_pooled, 25)
    whi_pool   = all_pooled[all_pooled <= q3_pool + 1.5 * iqr_pool].max()

    bot_max = whi_pool * 1.40
    top_min = whi_pool * 1.65
    top_max = all_pooled.max() * 1.08

    n         = len(all_vals_list)
    colours   = ["#4C9BE8", "#E8834C", "#5DBB63"]
    positions = list(range(1, n + 1))

    box_kw = dict(patch_artist=True, showfliers=False,
                  medianprops=dict(color="black", linewidth=2.5),
                  whiskerprops=dict(linewidth=1.4),
                  capprops=dict(linewidth=1.4),
                  boxprops=dict(linewidth=1.4))

    rng = np.random.default_rng(42)

    # reference line style — distinct dashes per baseline
    baseline_styles = [
        dict(color="#d62728", lw=1.1, ls=(0, (1, 1))),      # top-1  dotted red
        dict(color="#ff7f0e", lw=1.1, ls=(0, (3, 1))),      # top-2  dashed orange
        dict(color="#9467bd", lw=1.1, ls=(0, (5, 2))),      # top-4  dashed purple
        dict(color="#8c564b", lw=1.1, ls=(0, (7, 2))),      # top-8  dashed brown
        dict(color="#2ca02c", lw=1.1, ls=(0, (3, 1, 1, 1))),# top-16 dash-dot green
    ]
    baselines = _kl_baselines()   # [(k, kl_val, label), ...]

    for ax, ylim in [(ax_bot, (0, bot_max)), (ax_top, (top_min, top_max))]:
        bp = ax.boxplot(all_vals_list, positions=positions, widths=0.45, **box_kw)
        for patch, colour in zip(bp["boxes"], colours[:n]):
            patch.set_facecolor(colour)
            patch.set_alpha(0.75)
        ax.set_ylim(*ylim)
        ax.set_xlim(0.3, n + 0.7)
        ax.set_xticks(positions)
        ax.grid(axis="y", linestyle="--", alpha=0.25)

        # draw reference lines on whichever panel they fall into
        for (k, kl_val, lbl), style in zip(baselines, baseline_styles):
            ylo, yhi = ylim
            if ylo <= kl_val <= yhi:
                ax.axhline(kl_val, xmin=0, xmax=1, **style, zorder=1)
                ax.text(n + 0.72, kl_val, lbl, va="center", ha="left",
                        fontsize=7, color=style["color"])

        # jittered strip for each stage
        for pos, v, colour in zip(positions, all_vals_list, colours[:n]):
            jitter = rng.uniform(-0.12, 0.12, size=len(v))
            ylo, yhi = ylim
            mask = (v >= ylo) & (v <= yhi)
            ax.scatter(pos + jitter[mask], v[mask],
                       color=colour, alpha=0.30, s=11, zorder=3, linewidths=0)

    # annotate median + Q1/Q3 on bottom panel for each stage
    for pos, v in zip(positions, all_vals_list):
        med = float(np.median(v))
        q1  = float(np.percentile(v, 25))
        q3  = float(np.percentile(v, 75))
        ax_bot.text(pos + 0.14, med, f"med={med:.3f}", va="center",
                    ha="left", fontsize=7.5)
        ax_bot.text(pos + 0.14, q1,  f"Q1={q1:.3f}",  va="center",
                    ha="left", fontsize=7)
        ax_bot.text(pos + 0.14, q3,  f"Q3={q3:.3f}",  va="center",
                    ha="left", fontsize=7)

    ax_bot.set_xticklabels(stage_labels, fontsize=10)
    ax_top.set_xticklabels([])
    ax_bot.set_ylabel("KL divergence  (nats)", fontsize=9)
    ax_top.set_title(anchor_name, fontsize=10, pad=6)

    # broken axis cosmetics
    ax_top.spines["bottom"].set_visible(False)
    ax_bot.spines["top"].set_visible(False)
    ax_top.tick_params(bottom=False)

    d  = 0.015
    kw = dict(transform=ax_top.transAxes, color="k", clip_on=False, linewidth=1.2)
    ax_top.plot((-d, +d), (-d, +d), **kw)
    ax_top.plot((1 - d, 1 + d), (-d, +d), **kw)
    kw.update(transform=ax_bot.transAxes)
    ax_bot.plot((-d, +d), (1 - d, 1 + d), **kw)
    ax_bot.plot((1 - d, 1 + d), (1 - d, 1 + d), **kw)


def plot_anchor_kl_divergence(stages: list, out_dir: Path) -> None:
    """
    Broken-axis boxplot + jittered strip of per-allele KL divergence
    (anchor distribution vs Uniform(1/20)) for multiple pipeline stages.

    Parameters
    ----------
    stages : list of (label: str, df: pd.DataFrame)
        e.g. [("Post-Phase 1", df_p1), ("Post-Phase 2", df_p2)]

    One figure with two panels (anchor1 left, anchor2 right).
    The y-axis is broken so the tight IQR box and sparse outliers are
    both clearly visible.

    File saved:
        anchor_kl_divergence_boxplot.png
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(
        2, 2,
        figsize=(12, 7),
        gridspec_kw={"height_ratios": [1, 3], "hspace": 0.08},
    )

    stage_labels = [lbl for lbl, _ in stages]

    for col, (anchor_col, anchor_name) in enumerate([
        ("a1_res", "Anchor 1 (P2)"),
        ("a2_res", "Anchor 2 (last pos)"),
    ]):
        all_vals_list = []
        for label, df in stages:
            kl_series = _per_allele_kl(df, anchor_col)
            all_vals_list.append(kl_series.values)
            log.info(
                f"  KL {anchor_name} | {label} | "
                f"median={kl_series.median():.4f}  "
                f"p25={kl_series.quantile(0.25):.4f}  "
                f"p75={kl_series.quantile(0.75):.4f}  "
                f"max={kl_series.max():.4f}"
            )

        _draw_kl_broken_axis(
            ax_top=axes[0, col],
            ax_bot=axes[1, col],
            anchor_name=anchor_name,
            stage_labels=stage_labels,
            all_vals_list=all_vals_list,
        )

    fig.suptitle(
        "Per-allele anchor distribution KL divergence vs Uniform(1/20)\n"
        "Laplace-smoothed  |  each point = one MHC allele",
        fontsize=11, y=1.01,
    )
    plt.tight_layout()
    fname = "anchor_kl_divergence_boxplot.png"
    fig.savefig(out_dir / fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Saved -> {fname}")
    log.info("KL divergence boxplot saved.")


def report_high_kl_alleles(df_raw: pd.DataFrame,
                            df_phase1: pd.DataFrame,
                            df_phase2: pd.DataFrame,
                            out_dir: Path,
                            kl_delta_rel_threshold: float = 0.05) -> None:
    """
    Identify alleles where anchor sampling FAILED — i.e. domination persists
    after Phase 2.

    Two conditions must BOTH be true to flag an allele:
      1. kl_phase2 > kl_abs_threshold
             — still meaningfully dominated after sampling
             — threshold is the top-16 AA KL baseline, computed from
               _kl_baselines() so it is always consistent with the plot
      2. (kl_phase1 - kl_phase2) / kl_phase1 < kl_delta_rel_threshold
             — KL reduced by less than kl_delta_rel_threshold (default 5%)
               relative to its phase-1 value, meaning sampling did not help

    Outputs:
        high_kl_alleles.csv   — flagged alleles with sample counts + KL per stage
        Logs a summary table.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # derive absolute threshold from the top-16 AA baseline — same value
    # shown as the green reference line on the KL plot
    kl_abs_threshold = next(kl for k, kl, _ in _kl_baselines() if k == 16)

    log.info(f"  Flagging criteria:")
    log.info(f"    kl_phase2 > {kl_abs_threshold:.4f}  (top-16 AA baseline from _kl_baselines())")
    log.info(f"    AND relative KL reduction (p1→p2) / kl_p1 < {kl_delta_rel_threshold:.0%}"
             f"  (sampling did not meaningfully help)")

    raw_counts = df_raw[COL_MHC].value_counts()
    p1_counts  = df_phase1[COL_MHC].value_counts()
    p2_counts  = df_phase2[COL_MHC].value_counts()

    rows = []
    for anchor_col, anchor_name in [("a1_res", "anchor1_P2"),
                                     ("a2_res", "anchor2_last")]:
        kl_raw = _per_allele_kl(df_raw,    anchor_col)
        kl_p1  = _per_allele_kl(df_phase1, anchor_col)
        kl_p2  = _per_allele_kl(df_phase2, anchor_col)

        kl_delta_rel = (kl_p1 - kl_p2) / kl_p1.replace(0, np.nan)  # relative reduction

        still_high  = kl_p2 > kl_abs_threshold
        not_reduced = kl_delta_rel < kl_delta_rel_threshold
        flagged     = kl_p2[still_high & not_reduced].index

        for allele in flagged:
            rows.append({
                "allele":           allele,
                "anchor":           anchor_name,
                "n_raw":            int(raw_counts.get(allele, 0)),
                "n_phase1":         int(p1_counts.get(allele, 0)),
                "n_phase2":         int(p2_counts.get(allele, 0)),
                "kl_raw":           round(float(kl_raw.get(allele, np.nan)), 4),
                "kl_phase1":        round(float(kl_p1.get(allele, np.nan)), 4),
                "kl_phase2":        round(float(kl_p2[allele]), 4),
                "kl_reduction_pct": round(float(kl_delta_rel[allele]) * 100, 2),
            })

    if not rows:
        log.info("  No problematic alleles found — all alleles either have low KL "
                 "or were improved by phase 2 sampling.")
        return

    result = (pd.DataFrame(rows)
                .sort_values(["anchor", "kl_phase2"], ascending=[True, False])
                .reset_index(drop=True))

    csv_path = out_dir / "high_kl_alleles.csv"
    result.to_csv(csv_path, index=False)
    log.info(f"  Problematic alleles -> {csv_path}  ({len(result)} rows)")
    log.info(f"\n{result.to_string(index=False)}")


# ══════════════════════════════════════════════════════════════════
# 11. MAIN
# ══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="pMHC sampling pipeline")
    parser.add_argument("--input",   required=True)
    parser.add_argument("--output",  required=False)
    parser.add_argument("--mode",    required=True, choices=["binder", "nonbinder"],
                        help="Whether to process binder (label=1) or non-binder (label=0) peptides")
    parser.add_argument("--inspect", action="store_true", help="Peek at file and exit")
    parser.add_argument("--explore", action="store_true", help="Exploration only, no sampling")
    parser.add_argument("--plots",   default="plots",     help="Directory for plots and stats")
    parser.add_argument("--phases",  default="both", choices=["only_phase1", "both"],
                        help="Whether to run only phase 1 (MHC sampling) or both phases (MHC + anchor sampling).")
    args = parser.parse_args()
 
    if args.inspect:
        inspect_file(args.input)
        return
 
    mode_label  = "Binders" if args.mode == "binder" else "Non-binders"
    phase1_only = args.phases == "only_phase1"
 
    out_dir = Path(args.plots)
    df_raw = load_data(args.input, mode=args.mode, phases=args.phases)
 
    explore(df_raw, out_dir)
 
    if args.explore:
        save_stats([compute_stats(df_raw, f"RAW ({mode_label}, class I)")], out_dir)
        log.info("Exploration-only mode. Exiting.")
        return
 
    if not args.output:
        parser.error("--output is required unless using --inspect or --explore")
 
    df_phase1 = phase1_mhc_sampling(df_raw)
 
    if phase1_only:
        stages = [("Raw", df_raw), ("Post-Phase 1", df_phase1)]
 
        save_stats([
            compute_stats(df_raw,    f"RAW ({mode_label}, class I)"),
            compute_stats(df_phase1, f"POST-PHASE 1 ({mode_label}, MHC balanced)"),
        ], out_dir)
 
        plot_comparison(stages=stages, out_dir=out_dir)
 
        log.info("=== Per-allele anchor fraction plots ===")
        plot_per_allele_anchor_fractions(df_raw, "Raw", out_dir)
        plot_per_allele_anchor_fractions(df_phase1, "Post-Phase 1", out_dir)
 
        log.info("=== Anchor KL divergence boxplot ===")
        plot_anchor_kl_divergence(stages=[("Raw", df_raw), ("Post-Phase 1", df_phase1)], out_dir=out_dir)
 
        log.info("=== Anchor combination diversity ===")
        compute_anchor_combo_stats(stages=stages, out_dir=out_dir)
 
        log.info("Phase 1 only mode. Skipping Phase 2.")
        save_data(df_phase1, args.output)
 
    else:
        df_phase2 = phase2_anchor_sampling(df_phase1)
        stages = [("Raw", df_raw), ("Post-Phase 1", df_phase1), ("Post-Phase 2", df_phase2)]
 
        save_stats([
            compute_stats(df_raw,    f"RAW ({mode_label}, class I)"),
            compute_stats(df_phase1, f"POST-PHASE 1 ({mode_label}, MHC balanced)"),
            compute_stats(df_phase2, f"POST-PHASE 2 ({mode_label}, anchor balanced)"),
        ], out_dir)
 
        plot_comparison(stages=stages, out_dir=out_dir)
 
        log.info("=== Per-allele anchor fraction plots ===")
        plot_per_allele_anchor_fractions(df_phase1, "Post-Phase 1", out_dir)
        plot_per_allele_anchor_fractions(df_phase2, "Post-Phase 2", out_dir)
 
        log.info("=== Anchor KL divergence boxplot ===")
        plot_anchor_kl_divergence(
            stages=[("Post-Phase 1", df_phase1), ("Post-Phase 2", df_phase2)],
            out_dir=out_dir,
        )
 
        log.info("=== High-KL allele report ===")
        report_high_kl_alleles(df_raw, df_phase1, df_phase2, out_dir)
 
        log.info("=== Anchor combination diversity ===")
        compute_anchor_combo_stats(stages=stages, out_dir=out_dir)
 
        save_data(df_phase2, args.output)
 
 
if __name__ == "__main__":
    main()