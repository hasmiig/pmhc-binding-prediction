"""
merge_vectors_into_splits.py
============================
Merges per-residue pLDDT and PAE vectors into existing split parquets.

Steps:
1. For each split parquet (train/val/test across all thresholds and folds),
   left merge on allele + long_mer + struct_idx to add plddt_vector and pae_vector
2. Overwrite the split parquet in place

Usage:
    python merge_vectors_into_splits.py \
        --vectors_dir /path/to/plddt_pae_vectors/ \
        --splits_root /path/to/combined_dendrogram_splits/
"""

import argparse
import os
import pandas as pd
from pathlib import Path

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--vectors_dir", type=str, required=True,
                    help="Directory containing combined_with_vectors.parquet file from extraction")
parser.add_argument("--splits_root", type=str,
                    default="/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/data/combined_dendrogram_splits",
                    help="Root directory of dendrogram splits")
parser.add_argument("--thresholds", nargs="+", default=["50", "60", "70", "80"],
                    help="Thresholds to process")
parser.add_argument("--folds", type=int, default=5,
                    help="Number of folds (0-indexed)")
args = parser.parse_args()

vectors_dir = Path(args.vectors_dir)
splits_root = Path(args.splits_root)

vectors_df = pd.read_parquet(vectors_dir / "combined_with_vectors.parquet")
vectors_df = vectors_df[["allele", "long_mer", "struct_idx", "plddt_vector", "pae_vector"]]

# ── step 1: find all split parquets ───────────────────────────────────────────
split_paths = []

for thr in args.thresholds:
    thr_root = splits_root / thr / "splits"

    # fold train/val
    for fold in range(args.folds):
        fold_dir = thr_root / f"fold_{fold}"
        for part in ("train", "val"):
            p = fold_dir / f"{part}.parquet"
            if p.exists():
                split_paths.append(p)
            else:
                print(f"  WARNING: missing {p}")

    # global test
    test_p = thr_root / "test" / "test.parquet"
    if test_p.exists():
        split_paths.append(test_p)
    else:
        print(f"  WARNING: missing {test_p}")

print(f"\nFound {len(split_paths)} split parquets to update")

# ── step 2: merge vectors into each split parquet ─────────────────────────────
n_total_matched = 0
n_total_missing = 0

for split_path in split_paths:
    df = pd.read_parquet(split_path)
    n_before = len(df)

    # merge on allele + long_mer + struct_idx
    df = df.merge(
        vectors_df,
        on=["allele", "long_mer", "struct_idx"],
        how="left"
    )

    n_matched = df["plddt_vector"].notna().sum()
    n_missing = df["plddt_vector"].isna().sum()
    n_total_matched += n_matched
    n_total_missing += n_missing

    # overwrite in place
    df.to_parquet(split_path, index=False)

    print(f"  {split_path.relative_to(splits_root)} — "
          f"{n_before:,} rows, matched={n_matched:,}, missing={n_missing:,}")

print(f"\nDone.")
print(f"Total matched: {n_total_matched:,}")
print(f"Total missing: {n_total_missing:,}")
print(f"If missing > 0, check that allele/long_mer/struct_idx match between vectors and splits.")