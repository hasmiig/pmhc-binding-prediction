#!/bin/bash
#SBATCH --job-name=lr_model
#SBATCH --output=logs/lr_model_layer_norm_l2_%a_%j.out
#SBATCH --error=logs/lr_model_layer_norm_l2_%a_%j.err
#SBATCH --partition=scc-gpu
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --constraint=inet
#SBATCH --array=0-1   # 2 thresholds

THRESHOLDS=(70 80)
THR=${THRESHOLDS[$SLURM_ARRAY_TASK_ID]}

SPLITS_ROOT=/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/data/combined_dendrogram_splits
OUTPUT_DIR=/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/pmhc-binding-prediction/results/02_baseline_lr
SCRIPT=/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/pmhc-binding-prediction/scripts/05_baseline_dense/lr_model.py

mkdir -p logs
mkdir -p $OUTPUT_DIR

source /mnt/ceph-hdd/projects/scc_mgmn_soeding/hasmig/miniforge3/etc/profile.d/conda.sh
conda activate pmtopo

cd /user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/pmhc-binding-prediction/scripts/05_baseline_dense/

python $SCRIPT \
    --splits_root $SPLITS_ROOT \
    --output_dir  $OUTPUT_DIR \
    --thresholds  $THR \
    --folds 5 \
    --layer_norm \
    --regularization l2 \
    --reg_lambda 0.0001