"""
extract_full_vectors.py
=======================
Adds full per-residue pLDDT and neighbourhood-averaged PAE vectors to a
resampled parquet. Modeled on add_pae_features.py but saves full vectors
instead of scalar means.

For each structure:
  - plddt_vector : full pLDDT for all residues (MHC + peptide), shape (n_res,)
  - pae_vector   : neighbourhood-averaged PAE for all residues (MHC + peptide),
                   shape (n_res,). For each residue i:
                       PAE_sym = (PAE + PAE.T) / 2
                       pae_vector[i] = mean(PAE_sym[i, neighbours])
                   where neighbours are all residues within 8 Å Cα-Cα distance.

The pdb_path column in the parquet points to the structure directory:
    {pdb_path}/{id_name}_model_1_model_2_ptm.pdb
    {pdb_path}/{id_name}_model_1_model_2_ptm_plddt.npy
    {pdb_path}/{id_name}_model_1_model_2_ptm_predicted_aligned_error.npy

Usage:
    # single run
    python extract_full_vectors.py --parquet .../combined.parquet --output .../combined_with_vectors.parquet

    # parallelised: split into N chunks, run each, then merge
    python extract_full_vectors.py --parquet .../combined.parquet --output .../combined_with_vectors.parquet \
        --n_chunks 10 --chunk_idx 0
    ...
    python extract_full_vectors.py --parquet .../combined.parquet --output .../combined_with_vectors.parquet \
        --n_chunks 10 --chunk_idx 9

    # merge chunk outputs
    python extract_full_vectors.py --merge --merge_dir .../chunks_out/ --output .../combined_with_vectors.parquet
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

DIST_THRESHOLD = 8.0
MODEL_SUFFIX   = "_model_1_model_2_ptm"

# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--parquet",    type=str, default=None,
                    help="Path to input parquet (required unless --merge)")
parser.add_argument("--output",     type=str, required=True,
                    help="Output parquet path")
parser.add_argument("--alphafold_tsv_name", type=str,
                    default="alphafold_input_file.tsv",
                    help="Name of the AlphaFold input TSV inside chunk/alphafold/")

# parallelisation args
parser.add_argument("--n_chunks",   type=int, default=1,
                    help="Total number of chunks to split the parquet into (for parallelisation)")
parser.add_argument("--chunk_idx",  type=int, default=0,
                    help="Which chunk to process (0-indexed)")

# merge mode
parser.add_argument("--merge",      action="store_true",
                    help="Merge chunk output parquets into one")
parser.add_argument("--merge_dir",  type=str, default=None,
                    help="Directory containing chunk output parquets to merge")

args = parser.parse_args()


# ── helpers ────────────────────────────────────────────────────────────────────

def parse_ca_coords(pdb_path: Path) -> np.ndarray:
    """
    Extract Cα coordinates from a PDB file.
    Returns array of shape (n_residues, 3), sorted by residue number.
    """
    entries = []
    seen_resids = set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            resid = line[22:26].strip()
            if resid in seen_resids:
                continue
            seen_resids.add(resid)
            entries.append([float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54])])
    return np.array(entries, dtype=np.float32)


def compute_pae_vector(pae: np.ndarray, ca_coords: np.ndarray) -> list:
    """
    Compute neighbourhood-averaged PAE for every residue (MHC + peptide).

    PAE_sym_ij = (PAE_ij + PAE_ji) / 2
    pae_vector[i] = mean over neighbours j within DIST_THRESHOLD of PAE_sym[i, j]

    Returns a list of length n_residues.
    """
    pae_sym = (pae + pae.T) / 2.0
    n_res   = len(ca_coords)
    result  = []
    for i in range(n_res):
        dists      = np.sqrt(((ca_coords - ca_coords[i]) ** 2).sum(axis=1))
        neighbours = np.where(dists <= DIST_THRESHOLD)[0]
        if len(neighbours) == 0:
            neighbours = np.array([i])
        result.append(float(pae_sym[i, neighbours].mean()))
    return result


def build_length_lookup(chunk_alphafold_dirs: list, tsv_name: str) -> dict:
    """Build id_name -> (mhc_len, pep_len) from AlphaFold input TSVs."""
    lookup = {}
    for chunk_alphafold_dir in chunk_alphafold_dirs:
        tsv_path = Path(chunk_alphafold_dir) / tsv_name
        if not tsv_path.exists():
            continue
        tsv = pd.read_csv(tsv_path, sep="\t")
        for _, row in tsv.iterrows():
            id_name = row["targetid"].split("/")[0]
            mhc_seq, pep_seq = row["target_chainseq"].split("/")
            lookup[id_name] = (len(mhc_seq), len(pep_seq))
    return lookup


# ── merge mode ─────────────────────────────────────────────────────────────────

if args.merge:
    assert args.merge_dir, "--merge_dir is required in merge mode"
    merge_dir = Path(args.merge_dir)
    chunk_files = sorted(merge_dir.glob("chunk_*.parquet"))
    print(f"Merging {len(chunk_files)} chunk parquets from {merge_dir}")
    dfs = [pd.read_parquet(f) for f in chunk_files]
    merged = pd.concat(dfs, ignore_index=True)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_parquet(output_path, index=False)
    print(f"Merged {len(merged):,} rows → {output_path}")
    exit(0)


# ── main extraction ────────────────────────────────────────────────────────────

assert args.parquet, "--parquet is required unless --merge"
parquet_path = Path(args.parquet)
output_path  = Path(args.output)

print(f"Loading parquet: {parquet_path}")
df = pd.read_parquet(parquet_path)
print(f"  {len(df):,} rows total")

# split into chunks for parallelisation
if args.n_chunks > 1:
    chunk_size = len(df) // args.n_chunks
    start = args.chunk_idx * chunk_size
    end   = start + chunk_size if args.chunk_idx < args.n_chunks - 1 else len(df)
    df    = df.iloc[start:end].copy().reset_index(drop=True)
    print(f"  Processing chunk {args.chunk_idx}/{args.n_chunks-1}: "
          f"rows {start}–{end} ({len(df):,} rows)")
    # save chunk output to a temp location next to output
    output_path = output_path.parent / f"chunk_{args.chunk_idx:04d}.parquet"

# pre-load length lookups
print("Pre-loading TSV length lookups...")
chunk_alphafold_dirs = set()
for pdb_path_str in df["pdb_path"].unique():
    chunk_alphafold_dirs.add(Path(pdb_path_str).parent)

length_lookup = build_length_lookup(
    sorted(chunk_alphafold_dirs), args.alphafold_tsv_name
)
print(f"  Loaded lengths for {len(length_lookup):,} structures "
      f"from {len(chunk_alphafold_dirs)} chunk dirs")

# extract vectors
plddt_vectors = []
pae_vectors   = []
n_ok = n_missing = n_error = 0

for idx, row in df.iterrows():
    pdb_dir  = Path(row["pdb_path"])
    id_name  = pdb_dir.name
    pdb_file = pdb_dir / f"{id_name}{MODEL_SUFFIX}.pdb"
    plddt_file = pdb_dir / f"{id_name}{MODEL_SUFFIX}_plddt.npy"
    pae_file   = pdb_dir / f"{id_name}{MODEL_SUFFIX}_predicted_aligned_error.npy"

    if not pdb_file.exists() or not plddt_file.exists() or not pae_file.exists():
        plddt_vectors.append(None)
        pae_vectors.append(None)
        n_missing += 1
        if n_missing <= 5:
            print(f"  WARNING: missing files for {id_name}")
        continue

    try:
        plddt     = np.load(plddt_file)           # shape (n_res,)
        pae       = np.load(pae_file)             # shape (n_res, n_res)
        ca_coords = parse_ca_coords(pdb_file)     # shape (n_res, 3)

        # full pLDDT vector — all residues
        plddt_vec = plddt[:mhc_len + pep_len].tolist()

        # full PAE vector — neighbourhood-averaged for all residues
        pae_vec = compute_pae_vector(pae, ca_coords)

        plddt_vectors.append(plddt_vec)
        pae_vectors.append(pae_vec)
        n_ok += 1

    except Exception as e:
        plddt_vectors.append(None)
        pae_vectors.append(None)
        n_error += 1
        if n_error <= 5:
            print(f"  ERROR {id_name}: {e}")

    if (idx + 1) % 1000 == 0:
        print(f"  {idx+1:,}/{len(df):,}  ok={n_ok:,}  "
              f"missing={n_missing:,}  error={n_error:,}")

df["plddt_vector"] = plddt_vectors
df["pae_vector"]   = pae_vectors

print(f"\nDone.  ok={n_ok:,}  missing={n_missing:,}  error={n_error:,}")
print(f"NaN plddt_vector: {df['plddt_vector'].isna().sum():,}")
print(f"NaN pae_vector:   {df['pae_vector'].isna().sum():,}")

output_path.parent.mkdir(parents=True, exist_ok=True)
df.to_parquet(output_path, index=False)
print(f"Saved: {output_path}")