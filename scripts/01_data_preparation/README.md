# Data Preparation — IEDB Sampling and Balancing

This directory contains scripts for downloading, processing, and balancing pMHC data from IEDB before structure prediction.

## Scripts

### `pmhc_sampling.py`

**Main sampling pipeline** — Implements iterative median sampling to balance MHC allele representation and anchor residue diversity.

**Usage:**

```bash
# Inspect input file schema
python pmhc_sampling.py --input data.parquet --inspect

# Explore raw data distribution (generates plots, no sampling)
python pmhc_sampling.py --input data.parquet --mode binder --explore --plots plots/

# Run full sampling pipeline (two phases)
python pmhc_sampling.py \
    --input data.parquet \
    --mode binder \
    --phases both \
    --output sampled_binders.parquet \
    --plots output_plots/
```

**Key Options:**

| Flag | Description |
|------|-------------|
| `--input` | Input file path (`.parquet`, `.tsv`, `.csv`) |
| `--mode` | `binder` or `nonbinder` |
| `--phases` | `both` (default) or `only_phase1` |
| `--output` | Output file path for sampled data |
| `--plots` | Directory to save diagnostic plots and stats |
| `--inspect` | Print schema and first rows, then exit |
| `--explore` | Generate raw data plots without sampling |

**Outputs:**

- **Plots:** Allele distribution, peptide lengths, anchor residues, anchor combinations
- **Stats:** `stats_report.txt`, `stats_summary.csv`, `anchor_combo_stats.csv`
- **Data:** Final sampled dataset in parquet format

### `prepare_pmhc_data.py`

**Data preprocessing** — Downloads or loads IEDB data, filters for MHC-I, removes non-binders and low-affinity peptides.

**Usage:**

```bash
python prepare_pmhc_data.py \
    --input iedb_export.csv \
    --output iedb_mhc1_filtered.parquet
```

### `chunk_tsv.py`

**Utility** — Splits large TSV/CSV files into chunks for parallel processing.

**Usage:**

```bash
python chunk_tsv.py --input large_file.tsv --size 10000 --output chunks/
```

---

## Workflow

1. **Start with raw IEDB export**
   ```bash
   python prepare_pmhc_data.py --input iedb_export.csv --output iedb_mhc1_filtered.parquet
   ```

2. **Inspect the filtered data**
   ```bash
   python pmhc_sampling.py --input iedb_mhc1_filtered.parquet --inspect
   ```

3. **Explore allele and anchor distributions**
   ```bash
   python pmhc_sampling.py --input iedb_mhc1_filtered.parquet --mode binder --explore --plots plots/
   ```

4. **Run the sampling pipeline**
   ```bash
   # For non-binders (both phases recommended)
   python pmhc_sampling.py \
       --input iedb_mhc1_filtered.parquet \
       --mode nonbinder \
       --phases both \
       --output data_nonbinders_sampled.parquet \
       --plots results_nonbinders/
   
   # For binders (phase 1 only, to preserve anchor preferences)
   python pmhc_sampling.py \
       --input iedb_mhc1_filtered.parquet \
       --mode binder \
       --phases only_phase1 \
       --output data_binders_sampled.parquet \
       --plots results_binders/
   ```

5. **Combine for downstream tasks**
   - Pass `data_binders_sampled.parquet` and `data_nonbinders_sampled.parquet` to structure prediction

---

## Column Requirements

Input files must contain these columns:

| Column | Example | Description |
|--------|---------|-------------|
| `long_mer` | `GGKKKRKV` | Peptide amino acid sequence |
| `allele` | `HLA-A*02:01` | MHC allele name |
| `assigned_label` | `1.0` | Binding label (1.0 = binder, 0.0 = non-binder) |
| `mhc_class` | `1.0` | MHC class (1.0 = class I, 2.0 = class II) |

If column names differ, update the constants at the top of `pmhc_sampling.py`:

```python
COL_PEPTIDE = "your_peptide_col"
COL_MHC     = "your_allele_col"
COL_LABEL   = "your_label_col"
MHC_CLASS   = "your_mhc_class_col"
```

---

## Output Interpretation

### Stats Report

After running, check `stats_report.txt` for:
- Row and allele counts at each stage
- Median peptide count before/after sampling
- Peptide length filtering stats
- Anchor residue diversity metrics

### Plots

- **`comparison_allele_distribution.png`** — Verify alleles are balanced
- **`comparison_anchor1_residues.png` / `anchor2`** — Check anchor diversity
- **`anchor_combo_heatmap.png`** — Visual distribution of anchor combinations

### High KL Divergence Alleles

If using `--phases both`, check `high_kl_alleles.csv` for alleles where anchor dominance persists despite balancing. These may reflect true biological binding constraints.

---

## Tips & Troubleshooting

1. **File won't load?** Check with `--inspect` to verify column names
2. **Too few samples after sampling?** Reduce `SAMPLE_CAP` or use `--phases only_phase1`
3. **Anchor residues look odd?** Verify that position P2 and C-terminal are the true anchor positions for your MHC class
4. **Memory issues on large files?** Use `chunk_tsv.py` to split before sampling

---

## Dependencies

```
pandas
numpy
matplotlib
scipy
pyarrow  # for parquet support
```

Install with:
```bash
pip install pandas numpy matplotlib scipy pyarrow
```
