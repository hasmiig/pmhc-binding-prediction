"""
joint_split.py
==================
Creates HLA-stratified and anchor-stratified train/val/test splits from a
SINGLE combined binder + nonbinder dataset.

The two input parquets are concatenated and given a `label` column
(1 = binder, 0 = nonbinder).  All splits operate on this combined frame so
that binder/nonbinder rows always land in the same fold.

Output layout
-------------
<output_dir>/
  splits/
    hla/
      test.parquet               ← rows whose allele is in the held-out test set
      test_alleles.csv
      fold_summary.csv
      split_meta.json
      fold_1/
        train.parquet
        val.parquet
      fold_2/ …
    anchor/
      fold_summary.csv
      split_meta.json
      fold_1/
        train.parquet
        val.parquet
        test.parquet
        test_combos.csv
      fold_2/ …

HLA test-allele strategy
------------------------
Test alleles are drawn from alleles that appear in BOTH binders and nonbinders
(the intersection), ranked by binder frequency ascending so the rarest binder
alleles are held out first.  `--hla_test_frac` (default 0.20) controls the
fraction of shared alleles that go to test.

Usage
-----
python joint_hla_split.py \\
    --input_dir .../combined.parquet \\
    --output_dir    .../output_joint_splits \\
    --k 5
"""

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger(__name__)

COL_MHC     = "allele"
COL_PEPTIDE = "long_mer"


# ── helpers ────────────────────────────────────────────────────────────────────

def save_parquet(df: pd.DataFrame, path: Path, label: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    n_b = (df["label"] == 1).sum()
    n_n = (df["label"] == 0).sum()
    log.info(f"    {label:35s} {path.name}  ({len(df):,} rows  |  {n_b:,} binders  {n_n:,} nonbinders)")


# ── HLA split ──────────────────────────────────────────────────────────────────

def split_hla(df: pd.DataFrame, k: int, hla_test_frac: float, out_dir: Path) -> None:
    log.info("=" * 60)
    log.info("HLA SPLIT  (combined dataset, shared test alleles)")
    log.info("=" * 60)

    split_dir = out_dir / "splits" / "hla"
    split_dir.mkdir(parents=True, exist_ok=True)

    binder_alleles    = set(df.loc[df["label"] == 1, COL_MHC].unique())
    nonbinder_alleles = set(df.loc[df["label"] == 0, COL_MHC].unique())
    shared            = binder_alleles & nonbinder_alleles

    log.info(f"  Binder alleles      : {len(binder_alleles):,}")
    log.info(f"  Nonbinder alleles   : {len(nonbinder_alleles):,}")
    log.info(f"  Shared alleles      : {len(shared):,}")
    log.info(f"  Binder-only         : {len(binder_alleles - shared):,}  (go to CV folds)")
    log.info(f"  Nonbinder-only      : {len(nonbinder_alleles - shared):,}  (go to CV folds)")

    # Rank shared alleles by binder count ascending; hold out the rarest ones
    binder_counts = (
        df[df["label"] == 1]
        .groupby(COL_MHC)
        .size()
        .reindex(sorted(shared))          # restrict to shared
        .sort_values(ascending=True)
    )
    n_test       = max(1, int(np.floor(len(binder_counts) * hla_test_frac)))
    test_alleles = set(binder_counts.index[:n_test])

    log.info(f"  Test alleles ({hla_test_frac:.0%}) : {len(test_alleles):,}")

    # Global test set — all rows (binder + nonbinder) whose allele is held out
    df_test = df[df[COL_MHC].isin(test_alleles)].copy()
    df_cv   = df[~df[COL_MHC].isin(test_alleles)].copy().reset_index(drop=True)

    n_b = (df_test["label"] == 1).sum()
    n_n = (df_test["label"] == 0).sum()
    total = len(df_test)
    log.info(
        f"  Global test balance : {n_b:,} binders / {n_n:,} nonbinders "
        f"({n_b / total:.1%} binder)" if total else "  Global test set is empty!"
    )

    save_parquet(df_test, split_dir / "test.parquet", "global test")
    pd.DataFrame({"allele": sorted(test_alleles)}).to_csv(
        split_dir / "test_alleles.csv", index=False
    )

    # k-fold CV on the remaining combined data, split by allele
    cv_alleles = np.array(df_cv[COL_MHC].unique())
    rng = np.random.default_rng(42)
    rng.shuffle(cv_alleles)
    allele_folds = np.array_split(cv_alleles, k)

    fold_summary = []
    for i in range(k):
        val_alleles   = set(allele_folds[i])
        train_alleles = set(cv_alleles) - val_alleles

        df_val   = df_cv[df_cv[COL_MHC].isin(val_alleles)].copy()
        df_train = df_cv[df_cv[COL_MHC].isin(train_alleles)].copy()

        fold_dir = split_dir / f"fold_{i + 1}"
        save_parquet(df_train, fold_dir / "train.parquet", f"fold {i + 1} train")
        save_parquet(df_val,   fold_dir / "val.parquet",   f"fold {i + 1} val")

        fold_summary.append({
            "fold":              i + 1,
            "n_train":           len(df_train),
            "n_val":             len(df_val),
            "n_train_alleles":   len(train_alleles),
            "n_val_alleles":     len(val_alleles),
            "n_train_binders":   int((df_train["label"] == 1).sum()),
            "n_train_nonbinders":int((df_train["label"] == 0).sum()),
            "n_val_binders":     int((df_val["label"]   == 1).sum()),
            "n_val_nonbinders":  int((df_val["label"]   == 0).sum()),
        })

    pd.DataFrame(fold_summary).to_csv(split_dir / "fold_summary.csv", index=False)
    with open(split_dir / "split_meta.json", "w") as f:
        json.dump(
            {
                "split_mode":     "hla_joint_combined",
                "k":              k,
                "hla_test_frac":  hla_test_frac,
                "n_test_alleles": len(test_alleles),
                "n_cv_alleles":   int(df_cv[COL_MHC].nunique()),
                "n_test_rows":    len(df_test),
                "n_cv_rows":      len(df_cv),
                "note": (
                    "Single combined binder+nonbinder frame. "
                    "Test alleles drawn from binder∩nonbinder intersection, "
                    "ranked by binder frequency ascending (rarest held out). "
                    "CV folds split the remaining combined data by allele."
                ),
            },
            f,
            indent=2,
        )


# ── Anchor split ───────────────────────────────────────────────────────────────

def _add_anchor_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["a1_res"] = df[COL_PEPTIDE].str[1]
    df["a2_res"] = df[COL_PEPTIDE].str[-1]
    df["_combo"] = df["a1_res"] + df["a2_res"]
    return df


def split_anchor(df: pd.DataFrame, k: int, out_dir: Path) -> None:
    log.info("\n" + "=" * 60)
    log.info("ANCHOR SPLIT  (combined dataset)")
    log.info("=" * 60)

    split_dir = out_dir / "splits" / "anchor"
    split_dir.mkdir(parents=True, exist_ok=True)

    df = _add_anchor_columns(df)
    drop_cols = ["_combo", "a1_res", "a2_res"]

    all_combos = np.array(df["_combo"].unique())
    rng = np.random.default_rng(42)
    rng.shuffle(all_combos)
    combo_chunks = np.array_split(all_combos, k)

    fold_summary = []
    for i in range(k):
        test_combos  = set(combo_chunks[i])
        val_combos   = set(combo_chunks[(i + 1) % k])
        train_combos = set(np.concatenate([
            combo_chunks[j] for j in range(k) if j != i and j != (i + 1) % k
        ]))

        df_test  = df[df["_combo"].isin(test_combos) ].drop(columns=drop_cols, errors="ignore").copy()
        df_val   = df[df["_combo"].isin(val_combos)  ].drop(columns=drop_cols, errors="ignore").copy()
        df_train = df[df["_combo"].isin(train_combos)].drop(columns=drop_cols, errors="ignore").copy()

        fold_dir = split_dir / f"fold_{i + 1}"
        save_parquet(df_test,  fold_dir / "test.parquet",  f"fold {i + 1} test")
        save_parquet(df_val,   fold_dir / "val.parquet",   f"fold {i + 1} val")
        save_parquet(df_train, fold_dir / "train.parquet", f"fold {i + 1} train")
        pd.DataFrame({"anchor_combo": sorted(test_combos)}).to_csv(
            fold_dir / "test_combos.csv", index=False
        )

        fold_summary.append({
            "fold":               i + 1,
            "test_combos":        len(test_combos),
            "val_combos":         len(val_combos),
            "train_combos":       len(train_combos),
            "n_test":             len(df_test),
            "n_val":              len(df_val),
            "n_train":            len(df_train),
            "n_test_binders":     int((df_test["label"]  == 1).sum()),
            "n_test_nonbinders":  int((df_test["label"]  == 0).sum()),
            "n_val_binders":      int((df_val["label"]   == 1).sum()),
            "n_val_nonbinders":   int((df_val["label"]   == 0).sum()),
            "n_train_binders":    int((df_train["label"] == 1).sum()),
            "n_train_nonbinders": int((df_train["label"] == 0).sum()),
        })

    pd.DataFrame(fold_summary).to_csv(split_dir / "fold_summary.csv", index=False)
    with open(split_dir / "split_meta.json", "w") as f:
        json.dump(
            {
                "split_mode": "anchor_combined",
                "k":          k,
                "n_rows":     len(df),
                "n_combos":   len(all_combos),
                "note": (
                    "Single combined binder+nonbinder frame split by anchor combo (P2+last). "
                    "Rotating assignment: test=chunk_i, val=chunk_(i+1)%k, train=rest."
                ),
            },
            f,
            indent=2,
        )


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Joint HLA + anchor splits from a single combined binder/nonbinder dataset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input_dir", required=True,
        help="Input dir containing the combined binder/nonbinder dataset",
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Where to write the combined dataset and splits",
    )
    parser.add_argument("--k",             type=int,   default=5,    help="Number of CV folds")
    parser.add_argument("--hla_test_frac", type=float, default=0.20, help="Fraction of shared alleles held out as test")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    if out_dir.exists():
        log.warning(f"Output dir already exists: {out_dir}  — files may be overwritten")
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── splits ────────────────────────────────────────────────────────────────
    split_hla(df,    k=args.k, hla_test_frac=args.hla_test_frac, out_dir=out_dir)
    split_anchor(df, k=args.k, out_dir=out_dir)

    log.info(f"\nDone — all splits written to: {out_dir}")


if __name__ == "__main__":
    main()