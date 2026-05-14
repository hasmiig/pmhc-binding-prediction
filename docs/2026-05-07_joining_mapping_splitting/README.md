# Joining, Mapping, & Splitting

## Date
2026-05-07

## Analyst
Hasmig Aintablian

## What we did
### 1. Further Processing and Mapping data rows to PDB paths (binder/nonbinder separately)

Run `filter_map_resample.py` independently on binders and nonbinders. Each run performs three sequential stages:

1. **FILTER** — Selected the best structure per (allele, peptide) pair based on the mean of
   peptide and anchor pLDDT scores.
   - **Binders**: best structure only, `--plddt_threshold 0` (no filtering)
   - **Nonbinders**: `--keep_all` (both seeds retained), `--plddt_threshold 0` (no filtering) —
     nonbinder pLDDT scores are systematically low (~5% ≥ 80) so filtering would discard too many
  Since pLDDT is used as a classification feature, filtering by pLDDT threshold would introduce selection bias — only high-confidence structures would be included, artificially inflating downstream model performance. Therefore no pLDDT threshold is applied for either binders or nonbinders.

2. **MAP** — Joined filtered structures back to raw IEDB parquet data to recover:
   - Peptide sequences (`long_mer`) and binding labels (`assigned_label`)
   - MHC allele annotations
   - MHC protein sequences from `mhc1_encodings.csv`
   - Constructed paths to PDB files (`{pdb_base_dir}/{chunk}/alphafold/{allele}_{peptide}_{struct_idx}/`)

3. **RESAMPLE** — Applied Phase 1 iterative median sampling to balance MHC allele representation.
   - `--phase1_only` used for **both** binders and nonbinders.
   - Phase 2 (anchor residue balancing) was deliberately skipped for nonbinders here to avoid
     over-reducing the nonbinder count, which would worsen class imbalance in the combined dataset.
   - For binders, Phase 2 never runs by design (anchor preferences are biologically meaningful).

### 2. Joined the binder/nonbinder resampled datasets

Concatenated the two Phase 1-resampled parquets and added a `label` column:
- `label = 1` → binder
- `label = 0` → nonbinder

### 3. Splitting

Two complementary splitting strategies were applied to the combined labeled dataset:

#### HLA/Anchor Splits (`joint_split.py`)

- **HLA split**: Test alleles drawn from the intersection of binder and nonbinder alleles (shared
  alleles), ranked by binder frequency ascending so the rarest alleles are held out. `--hla_test_frac 0.20`
  holds out the rarest 20% of shared alleles as a global test set. The remaining combined data is
  split into 5 CV folds by allele.
  - Output: `splits/hla/test.parquet` + `splits/hla/fold_{1-5}/train.parquet` + `fold_{1-5}/val.parquet`

- **Anchor split**: Rows are grouped by anchor combination (P2 × C-terminal residue). All unique
  combinations are divided into 5 chunks with a rotating assignment: test = chunk_i, val = chunk_(i+1)%k,
  train = remaining three chunks.
  - Output: `splits/anchor/fold_{1-5}/{train,val,test}.parquet` + `test_combos.csv` per fold

#### Dendrogram Splits (`dendrogram_split.py`)

MHC-similarity-aware splits with peptide-level leakage control. All alleles were encoded using 49-position MHC pseudosequences × 14 physicochemical properties (686 features) and clustered via Ward hierarchical clustering into `--k_clusters 10` groups. The largest cluster was held out as the fixed test set. Precomputed Asgary peptide clusters were then used to run a greedy set-cover that selects test-split rows covering ≥ 90% of test alleles while minimising training-data loss — test peptide clusters are permanently excluded from all folds to prevent leakage. Unmatched peptides were assigned singleton anchor-pair clusters. The remaining clusters were distributed into 5 folds with per-fold set-cover selecting val rows.

Four peptide-cluster similarity thresholds (50, 60, 70, 80) were evaluated, each producing a separate set of 5-fold splits — 20 train/val combinations and 4 global test sets in total. A lower threshold groups more peptides together (less conservative), while a higher threshold produces more, smaller clusters (more conservative).

## Why we did it

Structure prediction introduces new allele imbalances: PMGen fails at different rates across alleles,
so the allele distribution after prediction differs from the distribution after the initial subsampling
step. Re-applying Phase 1 resampling corrects this structural bias before training.

The pLDDT CSV files from the QC step contain only allele/peptide identifiers and confidence scores —
they do not carry the full metadata needed for model training (MHC sequence, PDB paths, binding labels).
The MAP stage recovers this metadata by re-joining to the original IEDB parquet and encodings file.

Two complementary splitting strategies are maintained to probe distinct failure modes:
- **HLA-based splits** test whether the model generalises to structurally novel MHC alleles not seen
  during training.
- **Anchor-based splits** test generalisation to unseen anchor combinations within alleles that were
  seen during training.
- **Dendrogram splits** add MHC-similarity-aware stratification with peptide-level leakage control —
  this is the most stringent evaluation since entire MHC clades are held out and shared peptide clusters
  are excluded to prevent data leakage.

## Results

### Stage-by-Stage Data Flow

| Stage | Binders | Nonbinders |
|-------|---------|------------|
| Input pLDDT CSV (all structures) | ~328,000 rows (2 seeds × ~164k pairs) | ~240,000 rows (2 seeds × ~120k pairs) |
| After best-structure selection | ~164,000 unique pairs | ~240,000 retained (keep_all) |
| After pLDDT filter | pep_mean_pLDDT ≥ 80 applied | no filter (threshold = 0) |
| After MAP (join to IEDB + encodings) | recovered metadata + PDB paths | recovered metadata + PDB paths |
| After Phase 1 resample (cap = 1000) | ~153,000 rows | ~99,000 rows |

| Stage | Combined |
|-------|----------|
| After joining binder + nonbinder | ~252,000 rows |
| Label distribution | ~61% binder / ~39% nonbinder |
| Unique alleles | shared intersection used for HLA test split |


### Diagnostic Plots Generated

**`filter_map_resample.py` (per mode: binder / nonbinder):**
- `plots/plddt_distribution_before_after.png` — peptide and anchor pLDDT histograms before vs. after
  filtering, with threshold line
- `plots/allele_distribution_after_mapping.png` — peptide count per allele after mapping
- `plots/allele_distribution_after_resampling.png` — peptide count per allele after Phase 1 resample
- `plots/comparison_*.png` — stage-by-stage allele count comparison (from `pmhc_sampling.plot_comparison`)
- `plots/per_allele_anchor_fractions_Post-Phase_1.png` — per-allele P2/C-terminal anchor fractions
- `plots/anchor_kl_divergence_*.png` — KL divergence of anchor distributions across stages

**`dendrogram_split.py`:**
- `dendrogram.png` — Ward hierarchical clustering of MHC alleles (49-position pseudosequences,
  14 physicochemical features), coloured by cluster assignment

## Conclusion

This step bridges the gap between raw pLDDT QC outputs and training-ready split parquets. By
re-mapping structures to full IEDB metadata, re-resampling to correct prediction-induced allele bias,
and combining binders with nonbinders under a unified label scheme, we produce a ~252k-row combined
dataset. Three complementary splitting strategies (HLA-allele, anchor-combo, dendrogram-clade) ensure
downstream model evaluation assesses genuinely distinct axes of generalisation. However, for all downstream model training and evaluation, the dendrogram splits were selected as the primary splitting strategy. This is the most stringent evaluation setting since entire MHC clades are held out and shared peptide clusters are excluded to prevent data leakage. The HLA and anchor splits are retained for reference and ablation studies.

## Files

### Scripts
- [`scripts/03_filtering_analysis/filter_map_resample.py`](../../scripts/03_filtering_analysis/filter_map_resample.py) — Filter, map, resample (binder/nonbinder separately)
- [`scripts/03_filtering_analysis/joint_split.py`](../../scripts/03_filtering_analysis/joint_split.py) — HLA and anchor splits on combined dataset
- [`scripts/03_filtering_analysis/dendrogram_split.py`](../../scripts/03_filtering_analysis/dendrogram_split.py) — Dendrogram-stratified splits with peptide-level leakage control

### Outputs (on SCC)
- `trainprep/binder/mapped.parquet` — binder rows post-filter and post-map
- `trainprep/binder/resampled.parquet` — binder rows post-Phase 1 resample (~153k rows)
- `trainprep/binder/plots/` — diagnostic plots (binder)
- `trainprep/nonbinder/mapped.parquet` — nonbinder rows post-filter and post-map
- `trainprep/nonbinder/resampled.parquet` — nonbinder rows post-Phase 1 resample (~99k rows)
- `trainprep/nonbinder/plots/` — diagnostic plots (nonbinder)
- `trainprep/combined.parquet` — joint binder + nonbinder dataset with `label` column (~252k rows)
- `trainprep/splits/hla/` — HLA-stratified splits (`test.parquet`, `fold_{1-5}/train.parquet`, `fold_{1-5}/val.parquet`, `fold_summary.csv`, `split_meta.json`)
- `trainprep/splits/anchor/` — anchor-stratified splits (`fold_{1-5}/{train,val,test}.parquet`, `fold_{1-5}/test_combos.csv`, `fold_summary.csv`, `split_meta.json`)
- `trainprep/splits/dendro/` — dendrogram splits (`test/test.parquet`, `fold_{0-4}/{train,val}.parquet`, `dendrogram.png`, `fold_summary.csv`, `split_meta.json`)