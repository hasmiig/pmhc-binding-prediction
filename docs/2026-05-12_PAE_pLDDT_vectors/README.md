# PAE + pLDDT Per-Residue Vectors

## Date
2026-05-12

## Analyst
Hasmig Aintablian

## What we did

Extracted full per-residue AlphaFold confidence vectors and merged them into the existing split parquets.

- `extract_full_vectors.py` appends `plddt_vector` and `pae_vector` to a combined parquet.
- `merge_vectors_into_splits.py` merges those vectors into train/val/test split parquets by `allele`, `long_mer`, and `struct_idx`.

## Implementation (`extract_full_vectors.py`)

For each structure row in the input parquet:

1. Derive `id_name` from `basename(pdb_path)`.
2. Load structure outputs from:
   - `{pdb_path}/{id_name}_model_1_model_2_ptm.pdb`
   - `{pdb_path}/{id_name}_model_1_model_2_ptm_plddt.npy`
   - `{pdb_path}/{id_name}_model_1_model_2_ptm_predicted_aligned_error.npy`
3. Build a lookup of `id_name → (mhc_len, pep_len)` from chunk-level `alphafold_input_file.tsv` files.
4. Compute two vector features for every residue in the complex:
   - `plddt_vector`: raw pLDDT values for all residues (MHC + peptide)
   - `pae_vector`: neighbourhood-averaged PAE per residue

### PAE vector formula

```text
PAE_sym_ij = (PAE_ij + PAE_ji) / 2
neighbours = residues within 8 Å Ca–Ca of residue i
pae_vector[i] = mean(PAE_sym[i, neighbours])
```

The vector is computed over the full PAE matrix after symmetrisation, so each residue's value reflects its local structural uncertainty relative to nearby residues.

### Parallel extraction mode

The script supports chunked processing with `--n_chunks` and `--chunk_idx`, and can merge chunk outputs back together with `--merge --merge_dir`.

## Implementation (`merge_vectors_into_splits.py`)

1. Load `combined_with_vectors.parquet` from `--vectors_dir`.
2. Select `allele`, `long_mer`, `struct_idx`, `plddt_vector`, `pae_vector`.
3. Iterate over existing split parquets under `--splits_root/<threshold>/splits/`.
4. Left-merge the vector data into each split parquet on `allele`, `long_mer`, `struct_idx`.
5. Overwrite `train`, `val`, and `test` parquets in place.

## Why we did it

- `plddt_vector` preserves residue-level confidence for the full MHC+peptide complex.
- `pae_vector` encodes relative structural uncertainty with a local 8 Å neighbourhood, capturing pairwise placement quality.
- These vectors avoid reducing confidence information to scalars and make full per-residue features available for downstream modelling.
- Merging into the splits ensures the same training/evaluation datasets receive the additional structural features.

## Data flow

| Stage | Notes |
|-------|-------|
| Input combined parquet | contains `pdb_path`, `allele`, `long_mer`, `struct_idx` |
| AlphaFold outputs | PDB + `_plddt.npy` + `_predicted_aligned_error.npy` |
| Output combined parquet | adds `plddt_vector`, `pae_vector` |
| Split parquets | train/val/test parquets updated in place |

## Usage

### Extract vectors

```bash
python extract_full_vectors.py \
  --parquet /path/to/combined.parquet \
  --output /path/to/combined_with_vectors.parquet
```

Parallel extraction:

```bash
python extract_full_vectors.py \
  --parquet /path/to/combined.parquet \
  --output /path/to/combined_with_vectors.parquet \
  --n_chunks 10 --chunk_idx 0
```

Merge chunk outputs:

```bash
python extract_full_vectors.py \
  --merge --merge_dir /path/to/chunk_outputs/ \
  --output /path/to/combined_with_vectors.parquet
```

### Merge into splits

```bash
python merge_vectors_into_splits.py \
  --vectors_dir /path/to/plddt_pae_vectors/ \
  --splits_root /path/to/combined_dendrogram_splits/
```

## Files

- `scripts/05_baseline_dense/extract_full_vectors.py`
- `scripts/05_baseline_dense/merge_vectors_into_splits.py`

## Inputs

- `combined.parquet` — combined dataset with `pdb_path` and structure metadata
- AlphaFold outputs for each structure
- `alphafold_input_file.tsv` chunk-level TSVs for length lookup
- split parquets under `splits_root/<threshold>/splits/`

## Outputs

- `combined_with_vectors.parquet` — combined parquet with `plddt_vector` and `pae_vector`
- updated split parquets containing the same vector columns

## Notes

- Rows with missing PDB or NPY files are saved with `None` for both vectors.
- If merge results have missing rows, check whether `allele`, `long_mer`, and `struct_idx` match between the vector parquet and split parquets.
