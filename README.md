# pMHC Binding Prediction — Data Preparation Pipeline

End-to-end pipeline for curating and preparing a high-quality pMHC class I dataset for fine-tuning ProteinMPNN as a binding predictor. Starting from raw IEDB exports, the pipeline produces balanced, structure-predicted, QC-filtered train/val/test splits ready for model training.

For full methodology and design decisions, see [docs/METHODOLOGY.md](docs/METHODOLOGY.md).

---

## Overview

```
IEDB export
    │
    ▼
01_data_preparation/      ← filter, deduplicate, balance alleles & anchors
    │
    ▼
02_structure_prediction/  ← prepare PMGen input, run structure prediction (SLURM)
    │
    ▼
03_filtering_analysis/    ← pLDDT QC, re-balance, produce train/val/test splits
    │
    ▼
Training-ready dataset (parquet + split manifests)
```

---

## Repository Structure

```
pmhc-binding-prediction/
├── scripts/
│   ├── 01_data_preparation/
│   │   ├── prepare_pmhc_data.py     # IEDB filtering (MHC-I, binding labels)
│   │   ├── pmhc_sampling.py         # Iterative median sampling (allele + anchor balancing)
│   │   └── chunk_tsv.py             # Split large files for parallel processing
│   ├── 02_structure_prediction/
│   │   ├── prepare_pmgen_input.py   # Build TSV input for PMGen
│   │   └── run_pmgen.sh             # SLURM job script for PMGen structure prediction
│   ├── 03_filtering_analysis/
│   │   ├── compute_plddt_means.py   # Extract mean pLDDT per structure
│   │   ├── analyse_plddt.py         # Visualize pLDDT distributions and apply threshold
│   │   └── filter_map_train_prep.py # Filter, re-sample, split for training
│   └── utils/
│       └── setup_env.sh             # Environment setup for cluster jobs
├── PMGen/                           # PMGen submodule (structure prediction)
├── docs/
│   ├── METHODOLOGY.md               # Full design rationale and statistics
│   └── <dated run logs>/            # Per-run exploration outputs
└── README.md
```

---

## Requirements

```
pandas
numpy
matplotlib
scipy
pyarrow     # for .parquet files
```

```bash
pip install pandas numpy matplotlib scipy pyarrow
```

Structure prediction requires PMGen and a GPU node (see `scripts/utils/setup_env.sh`).

---

## Stage 1 — Data Preparation

**Directory:** `scripts/01_data_preparation/`

### 1a. Filter raw IEDB data

```bash
python prepare_pmhc_data.py \
    --input iedb_export.csv \
    --output iedb_mhc1_filtered.parquet
```

Keeps MHC class I records, removes non-binders (`assigned_label == 0`) and low-affinity peptides (BA > 500 nM).

### 1b. Inspect the filtered data

```bash
python pmhc_sampling.py --input iedb_mhc1_filtered.parquet --inspect
```

Prints schema and first rows without loading the full file. Verify column names match the constants at the top of `pmhc_sampling.py`.

### 1c. Explore allele and anchor distributions (optional but recommended)

```bash
python pmhc_sampling.py --input iedb_mhc1_filtered.parquet --mode binder --explore --plots plots/binders/
python pmhc_sampling.py --input iedb_mhc1_filtered.parquet --mode nonbinder --explore --plots plots/nonbinders/
```

### 1d. Run the sampling pipeline

```bash
# Binders — phase 1 only (anchor dominance is biological, preserve it)
python pmhc_sampling.py \
    --input iedb_mhc1_filtered.parquet \
    --mode binder \
    --phases only_phase1 \
    --output data_binders_sampled.parquet \
    --plots results/binders/

# Non-binders — both phases (anchor dominance is noise, balance it)
python pmhc_sampling.py \
    --input iedb_mhc1_filtered.parquet \
    --mode nonbinder \
    --phases both \
    --output data_nonbinders_sampled.parquet \
    --plots results/nonbinders/
```

**Phase 1** balances peptide counts across MHC alleles (iterative median sampling, cap 1000/allele).  
**Phase 2** balances anchor residue diversity within each allele (same strategy, applied at P2 and C-terminal).

For large files, use `chunk_tsv.py` to split before sampling:

```bash
python chunk_tsv.py --input large_file.tsv --size 10000 --output chunks/
```

### Input column requirements

| Column | Description |
|---|---|
| `long_mer` | Peptide amino acid sequence |
| `allele` | MHC allele name (e.g. `HLA-A*02:01`) |
| `assigned_label` | Binding label: `1.0` = binder, `0.0` = non-binder |
| `mhc_class` | MHC class: `1.0` = class I, `2.0` = class II |

If your file uses different column names, update the constants at the top of `pmhc_sampling.py`:

```python
COL_PEPTIDE = "long_mer"
COL_MHC     = "allele"
COL_LABEL   = "assigned_label"
MHC_CLASS   = "mhc_class"
SAMPLE_CAP  = 1000
```

---

## Stage 2 — Structure Prediction

**Directory:** `scripts/02_structure_prediction/`

### 2a. Prepare PMGen input

```bash
python prepare_pmgen_input.py
```

Merges sampled peptides with MHC sequences from `mhc1_encodings.csv` and writes a TSV for PMGen. Update the hardcoded paths at the top of the script to point to your parquet and MHC encoding files.

### 2b. Run PMGen on the cluster

```bash
# Submit one job per chunk
for chunk in chunks/*.tsv; do
    sbatch run_pmgen.sh $chunk
done
```

Each job generates two structures per pMHC pair (in `--initial_guess --multiple_anchors` mode) and writes PDB outputs to `outputs/<chunk>/`.

---

## Stage 3 — Filtering, Analysis, and Train Prep

**Directory:** `scripts/03_filtering_analysis/`

### 3a. Extract pLDDT means

```bash
python compute_plddt_means.py
```

Reads PMGen PDB outputs and writes `plddt_means_binder.csv` / `plddt_means_nonbinder.csv` with per-structure mean peptide pLDDT and anchor pLDDT.

### 3b. Analyse pLDDT distributions

```bash
python analyse_plddt.py
```

Generates pLDDT histograms, applies the 80-threshold filter, and saves the best structure per (allele, peptide) pair to `plddt_means_{mode}_best.csv`.

### 3c. Filter, re-sample, and split

```bash
# Binders, allele-stratified split
python filter_map_train_prep.py \
    --plddt_csv      outputs/binder/plddt_means_binder.csv \
    --parquet        iedb_mhc1_filtered.parquet \
    --mhc_encodings  data/mhc1_encodings.csv \
    --pdb_base_dir   outputs/binder \
    --output_dir     trainprep/binder/ \
    --mode           binder \
    --plddt_threshold 80 \
    --k              5 \
    --split_mode     hla
```

This script runs four sequential stages:

1. **Filter** — select best structure per (allele, peptide); apply pLDDT threshold
2. **Map** — join back to the original parquet to recover metadata and MHC sequences
3. **Resample** — re-run iterative median sampling to correct any allele bias introduced during prediction
4. **Split** — produce train/val/test splits in one of two modes:
   - `hla` — rare alleles (bottom 20% by frequency) held out as test; k-fold CV on remainder
   - `anchor` — k-fold CV stratified by anchor residue combinations (P2 × C-terminal)

---

## Key Numbers (from full run)

| Stage | Count |
|---|---|
| Initial IEDB MHC-I pairs | ~500,000 |
| After binding filter | ~250,000 |
| After allele balancing | ~100,000 |
| After structure prediction (pLDDT ≥ 80) | 87,187 |
| After final re-balancing | **63,817** |
| Unique MHC-I alleles | 426 |

---

## Documentation

- [docs/METHODOLOGY.md](docs/METHODOLOGY.md) — full rationale for sampling strategy, anchor analysis, pLDDT thresholds, and splitting design
- `docs/<date>/` — per-run exploration plots and stats
