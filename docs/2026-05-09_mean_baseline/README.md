# Mean Feature Threshold Baseline

## Date
2026-05-09

## Analyst
Hasmig Aintablian

## What we did

Evaluated a threshold-based baseline classifier using three mean structural confidence features —
pLDDT (peptide), pLDDT (anchor), and PAE (peptide neighbourhood) — as standalone binary classifiers.
No model is trained; the optimal decision threshold τ is found on the training data and applied
directly to val/test.

### Features

| Feature | Column | Direction |
|---------|--------|-----------|
| pLDDT (pep) | `pep_mean_plddt` | higher = more confident = more likely binder |
| pLDDT (anchor) | `anchor_mean_plddt` | higher = more confident = more likely binder |
| PAE (pep) | `pae_pep` | lower = better (negated before use so higher = more likely binder) |

### Threshold criterion (`plddt_cv_performance.py`)

For each feature and each training fold, the optimal threshold τ is selected as the point on the
ROC curve closest to the top-left corner:

```
τ* = argmin sqrt(FPR² + (1 - TPR)²)
```

**CV evaluation**: τ fit on `train.parquet`, metrics reported on `val.parquet` for each fold.

**Test evaluation**: median of the 5 per-fold τs applied to the global `test/test.parquet`. The
test set is never touched during threshold selection.

**Metrics**: F1, Precision, Recall, AUC, AP (Average Precision), TPR@FPR10% (true positive rate
at a 10% false positive rate budget).

### Diagnostic checks (`plddt_diagnostics.py`)

Two sanity checks run before interpreting results:

- **Q1 — pLDDT distributions per fold**: KDE plots of binder and nonbinder pLDDT distributions
  for each eval fold vs. the full combined dataset, with the per-fold optimal threshold marked.
  Checks whether fold splitting causes distribution shift that could inflate performance.

- **Q2 — Per-allele binder/nonbinder ratio**: Stacked bar chart of sample counts and binder
  fraction per allele, sorted by total count. Visualises class imbalance across alleles.

### Per-allele analysis (`plddt_per_allele_metrics.py` + `plddt_per_allele_plots.py`)

Per-allele metrics were computed across all split types and all four dendrogram peptide-cluster
thresholds (50, 60, 70, 80). Five analytical plots were generated from `per_allele_metrics.csv`:

1. **Plot 1** — Data quantity vs. performance: n_binders (sum across folds) vs. median AUC per
   allele, coloured by dendrogram threshold. Spearman ρ computed on alleles with ≥ 10 binders.

2. **Plot 2** — Cross-threshold stability heatmap: per-allele median AUC across all four
   dendrogram thresholds, coloured by RdYlGn. Alleles sorted by mean AUC descending.

3. **Plot 3** — Feature agreement scatter: pairwise AUC comparisons between pLDDT (pep),
   pLDDT (anchor), and PAE (pep), with Spearman ρ and Pearson r annotated.

4. **Plot 4** — Allele metric heatmap: alleles × (feature × metric) dense summary, median
   across all folds and dendrogram thresholds.

5. **Plot 5** — Evaluability heatmap: per-allele evaluability status (evaluable / single-class /
   not in val) across all four dendrogram thresholds, with a histogram of how many thresholds
   each allele is evaluable in.

## Why we did it

A threshold baseline on raw structural features establishes the minimum performance any downstream
model must surpass to be worth training. pLDDT and PAE are the primary structural confidence signals
output by PMGen/AlphaFold2 — if these alone separate binders from nonbinders, they provide strong
discriminative signal and the dataset is well-structured for classification. If the baseline already
achieves high AUC, it also validates that the combined parquet and splits were constructed correctly.

The "always predict binder" naive baseline gives F1 ≈ 0.77 on the ~61/39 binder/nonbinder
combined dataset — any meaningful classifier must exceed this.

## Results

### CV Performance — Dendrogram Splits (val set, 5 folds, `--feature both`)

Median across 5 folds for each dendrogram peptide-cluster threshold:

| Threshold | Feature | Median F1 | Median AUC | Median AP | Median TPR@FPR10 |
|-----------|---------|-----------|------------|-----------|-----------------|
| 50 | pLDDT (pep) | ~0.933 | ~0.927 | ~0.980 | ~0.779 |
| 50 | pLDDT (anchor) | ~0.935 | ~0.935 | ~0.983 | ~0.811 |
| 50 | PAE (pep) | ~0.914 | ~0.908 | ~0.974 | ~0.727 |
| 60 | pLDDT (pep) | ~0.937 | ~0.936 | ~0.984 | ~0.814 |
| 60 | pLDDT (anchor) | ~0.937 | ~0.941 | ~0.986 | ~0.831 |
| 60 | PAE (pep) | ~0.919 | ~0.919 | ~0.979 | ~0.766 |
| 70 | pLDDT (pep) | ~0.933 | ~0.944 | ~0.986 | ~0.836 |
| 70 | pLDDT (anchor) | ~0.935 | ~0.947 | ~0.988 | ~0.844 |
| 70 | PAE (pep) | ~0.918 | ~0.930 | ~0.983 | ~0.810 |
| 80 | pLDDT (pep) | ~0.927 | ~0.948 | ~0.985 | ~0.855 |
| 80 | pLDDT (anchor) | ~0.923 | ~0.949 | ~0.986 | ~0.857 |
| 80 | PAE (pep) | ~0.909 | ~0.933 | ~0.978 | ~0.820 |

### Test Performance — Dendrogram Splits (global test, median τ from 5 folds)

| Threshold | Feature | F1 | AUC | AP | TPR@FPR10 | τ |
|-----------|---------|-----|-----|-----|-----------|---|
| 50 | pLDDT (pep) | 0.946 | 0.930 | 0.992 | 0.793 | 67.1 |
| 50 | pLDDT (anchor) | 0.948 | 0.934 | 0.993 | 0.799 | 71.5 |
| 50 | PAE (pep) | 0.925 | 0.912 | 0.989 | 0.749 | −3.38 |
| 60 | pLDDT (pep) | 0.947 | 0.930 | 0.992 | 0.791 | 67.4 |
| 60 | pLDDT (anchor) | 0.948 | 0.934 | 0.993 | 0.799 | 71.7 |
| 60 | PAE (pep) | 0.925 | 0.912 | 0.989 | 0.750 | −3.38 |
| **70** | **pLDDT (pep)** | **0.947** | **0.930** | **0.992** | **0.792** | **67.4** |
| **70** | **pLDDT (anchor)** | **0.948** | **0.934** | **0.993** | **0.798** | **71.95** |
| **70** | **PAE (pep)** | **0.925** | **0.912** | **0.989** | **0.749** | **−3.36** |
| 80 | pLDDT (pep) | 0.936 | 0.930 | 0.986 | 0.784 | 67.5 |
| 80 | pLDDT (anchor) | 0.935 | 0.935 | 0.988 | 0.802 | 72.4 |
| 80 | PAE (pep) | 0.913 | 0.910 | 0.982 | 0.725 | −3.35 |

Results are **consistent across all four peptide-cluster thresholds** — performance is stable
regardless of how conservatively peptide clusters were defined for the dendrogram splits.

pLDDT (anchor) marginally outperforms pLDDT (pep) on AUC and TPR@FPR10 across all thresholds.
PAE (pep) is a weaker standalone classifier but still substantially above random (AUC ~0.91).

### Diagnostic Plots Generated

**`plddt_cv_performance.py`** (one PNG per threshold):
- `cv_performance_output_50_both.png` — boxplot + diamond CV summary (threshold 50)
- `cv_performance_output_60_both.png` — boxplot + diamond CV summary (threshold 60)
- `cv_performance_output_70_both.png` — boxplot + diamond CV summary (threshold 70)
- `cv_performance_output_80_both.png` — boxplot + diamond CV summary (threshold 80)

**`plddt_diagnostics.py`:**
- `plddt_diagnostics_q1.png` — pLDDT KDE distributions per fold vs. full dataset with τ marked
- `plddt_diagnostics_q2.png` — per-allele stacked binder/nonbinder counts and binder fraction

**`plddt_per_allele_metrics.py`** (3 features × 4 thresholds = 12 PNGs):
- `per_allele_boxplots_{feature}_{split}.png` — per-allele F1/Precision/Recall/AUC/AP/TPR@FPR10,
  boxplots across CV folds (val) with test diamonds overlaid

**`plddt_per_allele_plots.py`:**
- `plot1_nbinders_vs_auc.png` — data quantity vs. AUC per allele
- `plot2_threshold_stability.png` — per-allele AUC stability across 4 dendrogram thresholds
- `plot3_feature_agreement.png` — pairwise feature AUC scatter with ρ / r
- `plot4_allele_metric_heatmap.png` — alleles × (feature × metric) dense heatmap
- `plot5_evaluability.png` — per-allele evaluability across thresholds

## Conclusion

The mean feature threshold baseline achieves AUC ~0.93 and AP ~0.99 on the global test set across
all four dendrogram peptide-cluster thresholds, confirming that raw structural confidence scores
strongly separate binders from nonbinders in this dataset. Results are stable across thresholds,
indicating that the dendrogram split design does not materially alter performance estimates.
pLDDT (anchor) is the strongest single feature. Any downstream structural model must exceed these
numbers to demonstrate genuine predictive value beyond pLDDT/PAE alone.

## Files

### Scripts
- [`scripts/04_baseline_model/plddt_cv_performance.py`](../../scripts/04_baseline_model/plddt_cv_performance.py) — CV threshold baseline: τ selection, metrics, summary plots
- [`scripts/04_baseline_model/plddt_diagnostics.py`](../../scripts/04_baseline_model/plddt_diagnostics.py) — Q1 (fold distribution) and Q2 (per-allele imbalance) diagnostics
- [`scripts/04_baseline_model/plddt_per_allele_metrics.py`](../../scripts/04_baseline_model/plddt_per_allele_metrics.py) — per-allele metrics across all splits and dendrogram thresholds
- [`scripts/04_baseline_model/plddt_per_allele_plots.py`](../../scripts/04_baseline_model/plddt_per_allele_plots.py) — 5 analytical plots from per_allele_metrics.csv
- [`scripts/04_baseline_model/plddt_threshold.py`](../../scripts/04_baseline_model/plddt_threshold.py) — early exploratory ROC analysis (trial script; superseded by plddt_cv_performance.py)

### Results
- `results/01_baseline/cv_results/results_{50,60,70,80}_both.csv` — per-fold and test metrics for all features and thresholds
- `results/01_baseline/cv_performance_output_{50,60,70,80}_both.png` — CV summary plots
- `results/01_baseline/per_allele_results/per_allele_metrics.csv` — full per-allele metrics table
- `results/01_baseline/per_allele_results/per_allele_boxplots_*.png` — per-allele boxplots (12 files)
- `results/01_baseline/per_allele_results/plot{1-5}_*.png` — analytical summary plots
