"""
add_pae_features.py
===================
Computes the neighbourhood-averaged PAE feature (PAE_pep) for each row
in a resampled parquet and saves a new parquet with the extra column.

Formula (from MB baseline spec):
    PAE_sym_ij  = (PAE_ij + PAE_ji) / 2
    N_i         = all residues within 8 Å of peptide residue i (Ca-Ca)
    PAE_pep     = (1/P) * sum_i [ (1/|N_i|) * sum_{j in N_i} PAE_sym_ij ]

The pdb_path column (already in the parquet) points to the structure
directory. File names are constructed as:
    {pdb_path}/{id_name}_model_1_model_2_ptm.pdb
    {pdb_path}/{id_name}_model_1_model_2_ptm_predicted_aligned_error.npy

where id_name = basename of pdb_path.

Usage:
    python add_pae_features.py --parquet .../combined.parquet
    python add_pae_features.py --parquet .../combined.parquet --output .../combined_with_pae.parquet
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

DIST_THRESHOLD = 8.0
MODEL_SUFFIX   = "_model_1_model_2_ptm"

parser = argparse.ArgumentParser()
parser.add_argument("--parquet", type=str, required=True,
                    help="Path to combined.parquet (or any parquet with pdb_path column)")
parser.add_argument("--output",  type=str, default=None,
                    help="Output parquet path (default: same dir, combined_with_pae.parquet)")
parser.add_argument("--alphafold_tsv_name", type=str,
                    default="alphafold_input_file.tsv",
                    help="Name of the AlphaFold input TSV inside chunk/alphafold/")
args = parser.parse_args()

parquet_path = Path(args.parquet)
output_path  = Path(args.output) if args.output else \
               parquet_path.parent / "combined_with_pae.parquet"


# ── helpers ────────────────────────────────────────────────────────────────────

def parse_ca_coords(pdb_path: Path) -> np.ndarray:
    """
    Extract Ca coordinates from a PDB file (single chain).
    Returns array of shape (n_residues, 3), one Ca per residue.
    """
    coords      = []
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
            coords.append([float(line[30:38]),
                           float(line[38:46]),
                           float(line[46:54])])

    return np.array(coords, dtype=np.float32)


def compute_pae_pep(pae: np.ndarray, ca_coords: np.ndarray,
                    mhc_len: int, pep_len: int) -> float:
    """
    Compute PAE_pep scalar.

    PAE_sym_ij = (PAE_ij + PAE_ji) / 2
    N_i        = residues within DIST_THRESHOLD of peptide residue i
    PAE_pep    = mean over i of [ mean over N_i of PAE_sym_ij ]
    """
    pae_sym   = (pae + pae.T) / 2.0
    pep_start = mhc_len
    pep_end   = mhc_len + pep_len

    pep_vals = []
    for i in range(pep_start, pep_end):
        if i >= len(ca_coords):
            continue
        dists     = np.sqrt(((ca_coords - ca_coords[i]) ** 2).sum(axis=1))
        neighbors = np.where(dists <= DIST_THRESHOLD)[0]
        if len(neighbors) == 0:
            neighbors = np.array([i])
        pep_vals.append(pae_sym[i, neighbors].mean())

    return float(np.mean(pep_vals)) if pep_vals else float("nan")


def build_length_lookup(chunks: list, tsv_name: str) -> dict:
    """
    Pre-load all TSV files once and build a dict:
        id_name -> (mhc_len, pep_len)
    This avoids re-reading the TSV for every row.
    """
    lookup = {}
    for chunk_alphafold_dir in chunks:
        tsv_path = chunk_alphafold_dir / tsv_name
        if not tsv_path.exists():
            continue
        tsv = pd.read_csv(tsv_path, sep="\t")
        for _, row in tsv.iterrows():
            id_name = row["targetid"].split("/")[0]
            mhc_seq, pep_seq = row["target_chainseq"].split("/")
            lookup[id_name] = (len(mhc_seq), len(pep_seq))
    return lookup


# ── main ───────────────────────────────────────────────────────────────────────

print(f"Loading parquet: {parquet_path}")
df = pd.read_parquet(parquet_path)
print(f"  {len(df):,} rows")

# Pre-load all TSV length lookups — one read per chunk, not per row
print("Pre-loading TSV length lookups...")
chunk_alphafold_dirs = set()
for pdb_path_str in df["pdb_path"].unique():
    # pdb_path = .../chunk_XXX/alphafold/<id_name>
    # parent   = .../chunk_XXX/alphafold/
    chunk_alphafold_dirs.add(Path(pdb_path_str).parent)

length_lookup = build_length_lookup(
    sorted(chunk_alphafold_dirs), args.alphafold_tsv_name
)
print(f"  Loaded lengths for {len(length_lookup):,} structures "
      f"from {len(chunk_alphafold_dirs)} chunk dirs")

pae_pep_values = []
n_ok = n_missing = n_error = 0

for idx, row in df.iterrows():
    pdb_dir  = Path(row["pdb_path"])
    id_name  = pdb_dir.name
    pdb_file = pdb_dir / f"{id_name}{MODEL_SUFFIX}.pdb"
    pae_file = pdb_dir / f"{id_name}{MODEL_SUFFIX}_predicted_aligned_error.npy"

    if not pdb_file.exists() or not pae_file.exists():
        pae_pep_values.append(float("nan"))
        n_missing += 1
        if n_missing <= 5:
            print(f"  WARNING: missing files for {id_name}")
        continue

    try:
        pae       = np.load(pae_file)
        ca_coords = parse_ca_coords(pdb_file)

        if id_name in length_lookup:
            mhc_len, pep_len = length_lookup[id_name]
        else:
            # fallback: infer from PAE shape and known peptide sequence
            pep_len = len(str(row["long_mer"]))
            mhc_len = pae.shape[0] - pep_len

        val = compute_pae_pep(pae, ca_coords, mhc_len, pep_len)
        pae_pep_values.append(round(val, 4))
        n_ok += 1

    except Exception as e:
        pae_pep_values.append(float("nan"))
        n_error += 1
        if n_error <= 5:
            print(f"  ERROR {id_name}: {e}")

    if (idx + 1) % 1000 == 0:
        print(f"  {idx+1:,}/{len(df):,}  ok={n_ok:,}  "
              f"missing={n_missing:,}  error={n_error:,}")

df["pae_pep"] = pae_pep_values

print(f"\nDone.  ok={n_ok:,}  missing={n_missing:,}  error={n_error:,}")
print(f"pae_pep — min={df['pae_pep'].min():.3f}  max={df['pae_pep'].max():.3f}  "
      f"mean={df['pae_pep'].mean():.3f}  NaN={df['pae_pep'].isna().sum():,}")

df.to_parquet(output_path, index=False)
print(f"Saved: {output_path}")