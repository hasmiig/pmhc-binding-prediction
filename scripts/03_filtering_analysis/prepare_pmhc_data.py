#!/usr/bin/env python3
"""
prepare_pmhc_data.py
====================
Converts parquet split files from filter_map_train_prep.py into
ProteinMPNN-ready JSONL data directories and fixed_positions.json files.

Each parquet row has a `pdb_path` column pointing to a directory containing
a .pdb file (PMGen output format). This script finds the .pdb, parses it,
writes data.jsonl, and generates fixed_positions.json with chain A (MHC)
fixed and chain P (peptide) designable.

Usage
-----
# HLA split (1 global test + 5 folds x 2 splits)
python prepare_pmhc_data.py \\
    --splits_dir .../output_trainprep/binder_1/splits/hla \\
    --output_dir .../proteinmpnn_input/binder_hla \\
    --split_mode hla

# Anchor split (5 folds x 3 splits)
python prepare_pmhc_data.py \\
    --splits_dir .../output_trainprep/binder_1/splits/anchor \\
    --output_dir .../proteinmpnn_input/binder_anchor \\
    --split_mode anchor

Output structure (HLA mode)
---------------------------
output_dir/
  test/
    data.jsonl
    fixed_positions.json
  fold_1/
    train/
      data.jsonl
      fixed_positions.json
    val/
      data.jsonl
      fixed_positions.json
  fold_2/ ... fold_5/

Output structure (anchor mode)
-------------------------------
output_dir/
  fold_1/
    train/  val/  test/   (each with data.jsonl + fixed_positions.json)
  fold_2/ ... fold_5/
"""

import json
import logging
import argparse
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
log = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────
MHC_CHAIN     = "A"   # fixed — MHC heavy chain in PMGen PDBs
PEPTIDE_CHAIN = "P"   # designable — peptide chain in PMGen PDBs

COL_PDB_PATH  = "pdb_path"
COL_PEPTIDE   = "long_mer"
COL_MHC       = "allele"


# ══════════════════════════════════════════════════════════════════
# PDB PARSING  (adapted from prepare_data.py)
# ══════════════════════════════════════════════════════════════════

THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}


def find_pdb_file(pdb_dir: str) -> Path:
    """Find the .pdb file inside a PMGen output directory."""
    pdb_dir = Path(pdb_dir)
    if not pdb_dir.exists():
        raise FileNotFoundError(f"PDB directory not found: {pdb_dir}")
    pdbs = list(pdb_dir.glob("*.pdb"))
    if not pdbs:
        raise FileNotFoundError(f"No .pdb file found in: {pdb_dir}")
    return pdbs[0]


def parse_pdb(pdb_path: Path) -> dict | None:
    """
    Parse a PDB file into ProteinMPNN-compatible JSON format.
    Returns None if no valid residues found.
    """
    coords_by_chain = {}
    seqs_by_chain   = {}
    residue_order   = {}   # tracks insertion order per chain

    with open(pdb_path, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            resname   = line[17:20].strip()
            chain_id  = line[21]
            res_num   = line[22:27].strip()
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue

            if resname not in THREE_TO_ONE:
                continue
            if atom_name not in ["N", "CA", "C", "O"]:
                continue

            key = (chain_id, res_num, resname)
            if key not in residue_order:
                residue_order[key] = len(residue_order)
            coords_by_chain.setdefault(key, {})[atom_name] = [x, y, z]

    if not coords_by_chain:
        return None

    # group by chain preserving order
    chain_residues: dict[str, list] = {}
    for key in sorted(residue_order, key=lambda k: residue_order[k]):
        chain_id, res_num, resname = key
        atoms = coords_by_chain[key]
        if all(a in atoms for a in ["N", "CA", "C", "O"]):
            chain_residues.setdefault(chain_id, []).append((resname, atoms))

    result = {"name": pdb_path.parent.name}
    for chain_id in sorted(chain_residues):
        residues = chain_residues[chain_id]
        seq      = "".join(THREE_TO_ONE[r] for r, _ in residues)
        N_list   = [a["N"]  for _, a in residues]
        CA_list  = [a["CA"] for _, a in residues]
        C_list   = [a["C"]  for _, a in residues]
        O_list   = [a["O"]  for _, a in residues]

        result[f"seq_chain_{chain_id}"] = seq
        result[f"coords_chain_{chain_id}"] = {
            f"N_chain_{chain_id}":  N_list,
            f"CA_chain_{chain_id}": CA_list,
            f"C_chain_{chain_id}":  C_list,
            f"O_chain_{chain_id}":  O_list,
        }

    return result


# ══════════════════════════════════════════════════════════════════
# JSONL + FIXED POSITIONS WRITERS
# ══════════════════════════════════════════════════════════════════

def write_jsonl(df: pd.DataFrame, out_dir: Path) -> list[str]:
    """
    Parse all PDB directories in df[COL_PDB_PATH] and write data.jsonl.
    Returns list of successfully parsed structure names.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    jsonl_path = out_dir / "data.jsonl"

    names    = []
    n_ok     = 0
    n_miss   = 0
    n_err    = 0

    with open(jsonl_path, "w") as f:
        for _, row in df.iterrows():
            pdb_dir  = Path(row[COL_PDB_PATH])
            pdb_file = pdb_dir / f"{pdb_dir.name}_model_1_model_2_ptm.pdb"
            try:
                if not pdb_file.exists():
                    n_miss += 1
                    continue
                data = parse_pdb(pdb_file)
                if data is None:
                    n_err += 1
                    continue
                f.write(json.dumps(data) + "\n")
                names.append(data["name"])
                n_ok += 1
            except Exception as e:
                log.warning(f"  Error parsing {pdb_file}: {e}")
                n_err += 1

    log.info(
        f"  {out_dir.name}: {n_ok:,} OK | {n_miss:,} missing | {n_err:,} errors "
        f"→ {jsonl_path}"
    )
    return names


def write_fixed_positions(names: list[str], df: pd.DataFrame, out_dir: Path) -> None:
    """
    Write fixed_positions.json fixing MHC residues (chain A positions 1..len(mhc_sequence))
    for all structures. Peptide residues (last len(long_mer) positions) are left designable.

    PMGen outputs a single chain A where:
      residues 1..len(mhc_sequence)  = MHC
      residues len(mhc_sequence)+1.. = peptide (last N residues)
    """
    # build lookup: structure name -> (mhc_length, pep_length)
    # structure name = pdb stem, which matches data["name"] from _parse_pdb
    name_to_row = {}
    for _, row in df.iterrows():
        pdb_dir = Path(row["pdb_path"])
        name_to_row[pdb_dir.name] = row

    config = {}
    for name in names:
        row = name_to_row.get(name)
        if row is None:
            log.warning(f"    No parquet row found for {name} — skipping fixed positions")
            continue
        mhc_length = len(str(row["mhc_sequence"]))
        fixed_pos  = list(range(1, mhc_length + 1))
        config[name] = {"fixed_positions": {"A": fixed_pos}}

    path = out_dir / "fixed_positions.json"
    with open(path, "w") as f:
        json.dump(config, f, indent=2)
    log.info(f"  fixed_positions.json → {path}  ({len(config):,} entries, MHC residues fixed per structure)")


def process_split(parquet_path: Path, out_dir: Path) -> None:
    """Process a single parquet split into data.jsonl + fixed_positions.json."""
    if not parquet_path.exists():
        log.warning(f"  Parquet not found, skipping: {parquet_path}")
        return
    df    = pd.read_parquet(parquet_path)
    names = write_jsonl(df, out_dir)
    write_fixed_positions(names, df, out_dir)


# ══════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Prepare pMHC parquet splits for ProteinMPNN fine-tuning.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--splits_dir", required=True,
        help="Path to splits directory (e.g. .../splits/hla or .../splits/anchor)",
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Output directory for ProteinMPNN-ready data",
    )
    parser.add_argument(
        "--split_mode", required=True, choices=["hla", "anchor"],
        help=(
            "hla    : expects test.parquet + fold_i/train.parquet + fold_i/val.parquet\n"
            "anchor : expects fold_i/train.parquet + fold_i/val.parquet + fold_i/test.parquet"
        ),
    )
    parser.add_argument(
        "--k", type=int, default=5,
        help="Number of folds",
    )
    parser.add_argument(
        "--fold", type=int, default=None,
        help=(
            "If provided, process only this fold (1-indexed). "
            "Use for SLURM array jobs. If not provided, process all folds. "
            "For HLA mode, fold=0 processes the global test set only."
        ),
    )
    args = parser.parse_args()

    splits_dir = Path(args.splits_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info(f"prepare_pmhc_data.py — split_mode={args.split_mode}")
    log.info(f"  splits_dir : {splits_dir}")
    log.info(f"  output_dir : {output_dir}")
    log.info(f"  fold       : {args.fold if args.fold is not None else 'all'}")
    log.info("=" * 60)

    if args.split_mode == "hla":
        # fold=0 → global test set only
        if args.fold is None or args.fold == 0:
            log.info("--- Global test set ---")
            process_split(splits_dir / "test.parquet", output_dir / "test")

        # fold=1..k → that fold's train + val
        folds_to_run = range(1, args.k + 1) if args.fold is None else                        ([args.fold] if args.fold != 0 else [])
        for i in folds_to_run:
            log.info(f"--- Fold {i} ---")
            fold_dir = splits_dir / f"fold_{i}"
            out_fold = output_dir / f"fold_{i}"
            process_split(fold_dir / "train.parquet", out_fold / "train")
            process_split(fold_dir / "val.parquet",   out_fold / "val")

    elif args.split_mode == "anchor":
        folds_to_run = range(1, args.k + 1) if args.fold is None else [args.fold]
        for i in folds_to_run:
            log.info(f"--- Fold {i} ---")
            fold_dir = splits_dir / f"fold_{i}"
            out_fold = output_dir / f"fold_{i}"
            process_split(fold_dir / "train.parquet", out_fold / "train")
            process_split(fold_dir / "val.parquet",   out_fold / "val")
            process_split(fold_dir / "test.parquet",  out_fold / "test")

    log.info("=" * 60)
    log.info(f"DONE — outputs saved to: {output_dir}")
    log.info("=" * 60)


if __name__ == "__main__":
    main()