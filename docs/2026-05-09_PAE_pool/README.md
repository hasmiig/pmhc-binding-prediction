# PAE Neighbourhood Pooling Feature

## Date
2026-05-09

## Analyst
Hasmig Aintablian

## What we did

Computed a single per-structure scalar feature, **PAE_pep**, from the AlphaFold2/PMGen Predicted
Aligned Error (PAE) matrix and Ca coordinates, and appended it as a new column to the combined
parquet.

### Formula

```
PAE_sym_ij  = (PAE_ij + PAE_ji) / 2
N_i         = all residues within 8 Å of peptide residue i  (Ca–Ca distance)
PAE_pep     = (1/P) * sum_i [ (1/|N_i|) * sum_{j in N_i} PAE_sym_ij ]
```

Where P is the number of peptide residues. The symmetrised PAE removes directionality
(`PAE_ij ≠ PAE_ji` in the raw matrix). Neighbourhood pooling (8 Å shell) captures the
local structural context of each peptide position rather than isolated per-residue error.

### Implementation (`add_pae_features.py`)

For each row in the input parquet:

1. **Locate files** — derive `id_name` from `basename(pdb_path)` and resolve:
   - `{pdb_path}/{id_name}_model_1_model_2_ptm.pdb`
   - `{pdb_path}/{id_name}_model_1_model_2_ptm_predicted_aligned_error.npy`

2. **Load lengths** — pre-load all chunk-level `alphafold_input_file.tsv` files once to build an
   `id_name → (mhc_len, pep_len)` lookup (avoids re-reading per row). Fallback: infer `pep_len`
   from `long_mer` and `mhc_len` from PAE matrix shape.

3. **Parse Ca coordinates** — read all Ca atoms from the PDB, deduplicated by residue ID.

4. **Compute PAE_pep** — symmetrise the full PAE matrix, identify 8 Å neighbours for each peptide
   residue, pool mean symmetrised PAE over the neighbourhood, average across all peptide residues.
   Returns `NaN` for rows where files are missing or parsing fails.

5. **Write output** — append `pae_pep` column (float, rounded to 4 d.p.) and save new parquet.

## Why we did it

pLDDT measures per-residue confidence in absolute coordinate placement. PAE measures AlphaFold2's
confidence in the **relative positioning** between pairs of residues — specifically, the expected
positional error at residue j when aligning on residue i's frame. For pMHC binding prediction,
relative placement of the peptide within the MHC groove is the structurally meaningful quantity.

Since pLDDT is not used as a filter (filtering by pLDDT would introduce selection bias — only
high-confidence structures would survive, artificially inflating downstream model performance), PAE_pep
provides a complementary structural confidence signal that the classifier can use as a feature. A low
PAE_pep indicates the peptide is placed confidently relative to its local neighbourhood; a high value
suggests poor positional certainty, which may correlate with non-binding.

The neighbourhood pooling (8 Å shell, Ca–Ca) was chosen to match the contact radius commonly used
in structure-based feature extraction and is consistent with the MB baseline feature specification.

## Results

### Data Flow

| Stage | Count |
|-------|-------|
| Input combined parquet rows | ~252,000 |
| Successfully computed PAE_pep | ~252,000 (minus missing / error) |
| Missing files (no PDB or PAE .npy) | reported in run log |
| Parse errors | reported in run log |

### PAE_pep Distribution

| Statistic | Value |
|-----------|-------|
| Min | TBD |
| Max | TBD |
| Mean | TBD |
| NaN count | TBD |

### Diagnostic Plots Generated

None — `add_pae_features.py` produces no plots. PAE_pep distribution can be visualised
downstream from the output parquet.

## Conclusion

`add_pae_features.py` augments the combined parquet with a single neighbourhood-pooled PAE scalar
per structure. The output `combined_with_pae.parquet` is a drop-in replacement for `combined.parquet`
with one additional column (`pae_pep`) and is used as input to all downstream split and training steps.

## Files

### Script
- [`scripts/03_filtering_analysis/add_pae_features.py`](../../scripts/03_filtering_analysis/add_pae_features.py) — computes PAE_pep and writes augmented parquet

### Inputs (on SCC)
- `trainprep/combined.parquet` — combined binder + nonbinder dataset with `pdb_path` column
- `{pdb_base_dir}/{chunk}/alphafold/{id_name}/{id_name}_model_1_model_2_ptm.pdb` — PMGen PDB structure
- `{pdb_base_dir}/{chunk}/alphafold/{id_name}/{id_name}_model_1_model_2_ptm_predicted_aligned_error.npy` — PAE matrix
- `{pdb_base_dir}/{chunk}/alphafold/alphafold_input_file.tsv` — per-chunk TSV with MHC/peptide chain lengths

### Output (on SCC)
- `trainprep/combined_with_pae.parquet` — combined parquet with added `pae_pep` column (~252k rows)
