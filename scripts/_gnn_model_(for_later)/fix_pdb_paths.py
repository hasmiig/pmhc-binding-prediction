"""
fix_pdb_paths.py
================
One-time script to remap pdb_path column in all split parquets
from cbscratch paths to the actual accessible path on this cluster.

Run once before training:
    python fix_pdb_paths.py \
        --splits_dir /path/to/splits/hla \
        --old_prefix /cbscratch/amir_share/PMTopo/data/binder_pdbs/chunks \
        --new_prefix /user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/data/pmgen_outputs/binder
"""

import argparse
import pandas as pd
from pathlib import Path

def fix_parquet(path: Path, old_prefix: str, new_prefix: str):
    if not path.exists():
        print(f"  skipping (not found): {path}")
        return
    df = pd.read_parquet(path)
    df['pdb_path'] = df['pdb_path'].str.replace(old_prefix, new_prefix, regex=False)
    df.to_parquet(path, index=False)
    print(f"  fixed: {path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--splits_dir",  required=True)
    parser.add_argument("--old_prefix",  required=True)
    parser.add_argument("--new_prefix",  required=True)
    parser.add_argument("--split_mode",  choices=["hla", "anchor"], default="hla")
    parser.add_argument("--k",           type=int, default=5)
    args = parser.parse_args()

    splits_dir = Path(args.splits_dir)

    if args.split_mode == "hla":
        fix_parquet(splits_dir / "test.parquet", args.old_prefix, args.new_prefix)
        for i in range(1, args.k + 1):
            fix_parquet(splits_dir / f"fold_{i}/train.parquet", args.old_prefix, args.new_prefix)
            fix_parquet(splits_dir / f"fold_{i}/val.parquet",   args.old_prefix, args.new_prefix)
    elif args.split_mode == "anchor":
        for i in range(1, args.k + 1):
            fix_parquet(splits_dir / f"fold_{i}/train.parquet", args.old_prefix, args.new_prefix)
            fix_parquet(splits_dir / f"fold_{i}/val.parquet",   args.old_prefix, args.new_prefix)
            fix_parquet(splits_dir / f"fold_{i}/test.parquet",  args.old_prefix, args.new_prefix)

    print("Done.")

if __name__ == "__main__":
    main()