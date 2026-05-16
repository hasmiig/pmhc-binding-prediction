import argparse
import logging
import numpy as np
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

COL_MHC = "allele"

def main():
    parser = argparse.ArgumentParser(
        description="Combine binder and nonbinder datasets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--binder_dir", required=True,
        help="Binder output dir from filter_map_resample.py (must contain resampled.parquet)",
    )
    parser.add_argument(
        "--nonbinder_dir", required=True,
        help="Nonbinder output dir from filter_map_resample.py (must contain resampled.parquet)",
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Where to write the combined dataset",
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    if out_dir.exists():
        log.warning(f"Output dir already exists: {out_dir}  — files may be overwritten")
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── load & label ──────────────────────────────────────────────────────────
    log.info("Loading resampled parquets...")
    # df_binder    = pd.read_parquet(Path(args.binder_dir)    / "resampled.parquet")
    # df_nonbinder = pd.read_parquet(Path(args.nonbinder_dir) / "resampled.parquet")
    df_binder    = pd.read_parquet(Path(args.binder_dir)    / "mapped.parquet")
    df_nonbinder = pd.read_parquet(Path(args.nonbinder_dir) / "mapped.parquet")

    df_binder["label"]    = 1
    df_nonbinder["label"] = 0

    df = pd.concat([df_binder, df_nonbinder], ignore_index=True)

    log.info(f"  Binder rows         : {len(df_binder):,}  ({df_binder[COL_MHC].nunique():,} alleles)")
    log.info(f"  Nonbinder rows      : {len(df_nonbinder):,}  ({df_nonbinder[COL_MHC].nunique():,} alleles)")
    log.info(f"  Combined rows       : {len(df):,}  ({df[COL_MHC].nunique():,} unique alleles)")

    # Save combined for reference
    combined_path = out_dir / "combined.parquet"
    df.to_parquet(combined_path, index=False)
    log.info(f"  Combined parquet saved → {combined_path}")

if __name__ == "__main__":
    main()