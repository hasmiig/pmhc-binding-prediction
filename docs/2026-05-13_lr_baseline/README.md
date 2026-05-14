# LR Baseline

## Date
2026-05-13

## Analyst
Hasmig Aintablian

## What we did

Built and evaluated a logistic regression baseline with attention pooling over per-residue structural confidence vectors.

- `lr_model.py` trains a pMHC binding classifier using concatenated pooled representations from `plddt_vector` and `pae_vector`.
- `lr_plot_results.py` plots the ablation study across regularization and layer-normalization settings.
- Training is performed on dendrogram splits with both CV folds and a final retrained test evaluation.

## Implementation (`lr_model.py`)

The model encodes each structure as two feature sequences:

- `plddt_vector` + MHC/peptide token flags
- `pae_vector` + MHC/peptide token flags

Each feature sequence is projected from `(L, 3)` to `(L, d)`, sinusoidal positional encodings are added, RoPE is applied, and attention pooling reduces each sequence to a fixed-size vector.

The pooled pLDDT and PAE representations are concatenated and passed through a single linear classification layer producing a sigmoid probability.

### Training procedure

- For each threshold (`50`, `60`, `70`, `80`):
  - Train on 5 CV folds using `train.parquet` and `val.parquet` under
    `splits_root/<threshold>/splits/fold_{fold}`.
  - Track validation AUC and save the best model for each fold.
  - Retrain a final model on all combined train+val folds and evaluate on
    `splits_root/<threshold>/splits/test/test.parquet`.

### Model configuration

- base architecture: attention pooling over per-feature sequences
- embedding dimension: `d=8`
- maximum sequence length: `200`
- optimizer: Adam
- learning rate: `1e-3`
- epochs: `50`
- regularization: `none`, `l1`, or `l2`
- optional `LayerNorm` after attention pooling

## Implementation (`lr_plot_results.py`)

- Reads the `all_metrics.csv` produced by `lr_model.py`.
- Generates ablation plots for CV validation folds and test results.
- Configurations shown include:
  - `baseline`
  - `layer_norm`
  - `l2`
  - `l1`
  - `layer_norm_l2`
  - `layer_norm_l1`

## Inputs

- `splits_root` — directory containing dendrogram splits:
  `splits_root/<threshold>/splits/`
- `plddt_vector` and `pae_vector` already merged into split parquets via the
  separate vector extraction workflow
- `train.parquet`, `val.parquet`, `test.parquet` for each threshold/fold

## Outputs

- `results/02_baseline_lr/<run_tag>/fold_{fold}/best_model.pt`
- `results/02_baseline_lr/<run_tag>/final_model/final_model.pt`
- `results/02_baseline_lr/all_metrics.csv`
- plots produced by `lr_plot_results.py`

## Usage

### Run baseline training

```bash
python lr_model.py \
  --splits_root /path/to/combined_dendrogram_splits/ \
  --output_dir /path/to/results/02_baseline_lr/ \
  --thresholds 50 60 70 80 \
  --folds 5
```

### Run with L1 regularization

```bash
python lr_model.py \
  --splits_root /path/to/combined_dendrogram_splits/ \
  --output_dir /path/to/results/02_baseline_lr/ \
  --thresholds 50 60 70 80 \
  --folds 5 \
  --regularization l1 \
  --reg_lambda 1e-4
```

### Run with L2 regularization

```bash
python lr_model.py \
  --splits_root /path/to/combined_dendrogram_splits/ \
  --output_dir /path/to/results/02_baseline_lr/ \
  --thresholds 50 70 80 \
  --folds 5 \
  --regularization l2 \
  --reg_lambda 1e-4
```

### Run with LayerNorm

```bash
python lr_model.py \
  --splits_root /path/to/combined_dendrogram_splits/ \
  --output_dir /path/to/results/02_baseline_lr/ \
  --thresholds 50 60 70 80 \
  --folds 5 \
  --layer_norm
```

### Plot results

```bash
python lr_plot_results.py \
  --metrics_csv /path/to/results/02_baseline_lr/all_metrics.csv \
  --output_png /path/to/results/02_baseline_lr/lr_ablation.png
```

## Files

- `scripts/05_baseline_dense/lr_model.py`
- `scripts/05_baseline_dense/lr_plot_results.py`
- `scripts/05_baseline_dense/run_lr_model.sh`
- `scripts/05_baseline_dense/run_lr_model_l1.sh`
- `scripts/05_baseline_dense/run_lr_model_l2.sh`
- `scripts/05_baseline_dense/run_lr_model_layer_norm.sh`
- `scripts/05_baseline_dense/run_lr_model_layer_norm_l1.sh`
- `scripts/05_baseline_dense/run_lr_model_layer_norm_l2.sh`

## Notes

- This doc does not cover vector extraction; `plddt_vector` and `pae_vector`
  are expected to exist in the split parquets already.
- Missing vector rows are dropped before training.
