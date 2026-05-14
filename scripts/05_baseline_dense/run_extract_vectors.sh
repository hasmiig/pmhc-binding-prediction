#!/bin/bash
#SBATCH --job-name=extract_vectors
#SBATCH --output=logs/extract_vectors_%a.out
#SBATCH --error=logs/extract_vectors_%a.err
#SBATCH --partition=scc-cpu
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-19   # adjust to n_chunks - 1

N_CHUNKS=20            # adjust to match --array upper bound + 1
PARQUET=/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/data/output_filter_map_resample/combined.parquet
OUTPUT_DIR=/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/data/plddt_pae_vectors
SCRIPT=/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/pmhc-binding-prediction/scripts/05_baseline_dense/extract_full_vectors.py

mkdir -p logs
mkdir -p $OUTPUT_DIR

source $(conda info --base)/etc/profile.d/conda.sh
conda activate pmtopo

python $SCRIPT \
    --parquet $PARQUET \
    --output $OUTPUT_DIR/chunk_${SLURM_ARRAY_TASK_ID}.parquet \
    --n_chunks $N_CHUNKS \
    --chunk_idx $SLURM_ARRAY_TASK_ID