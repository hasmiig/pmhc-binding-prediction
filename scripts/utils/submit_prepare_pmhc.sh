#!/bin/bash
# ================================================================
# submit_prepare_pmhc.sh
# ================================================================
# Parallelizes prepare_pmhc_data.py across folds using SLURM arrays.
#
# Array jobs 1-5 : process fold_1 to fold_5 (HLA train/val + anchor train/val/test)
# Array job  0   : process HLA global test set only
#
# Submit with:
#   sbatch --array=0-5 submit_prepare_pmhc.sh binder
#   sbatch --array=0-5 submit_prepare_pmhc.sh nonbinder
# ================================================================
#SBATCH --job-name=prep_pmhc
#SBATCH --output=logs/prep_pmhc_%a.out
#SBATCH --error=logs/prep_pmhc_%a.err
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=2
#SBATCH --partition=soeding

MODE=${1:-binder}   # binder or nonbinder — passed as first argument
FOLD=$SLURM_ARRAY_TASK_ID

# ── paths — adjust as needed ──────────────────────────────────────
PROJECT=/projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig
SCRIPTS=$PROJECT/scripts
DATA=$PROJECT/data

SPLITS_BASE=$DATA/output_trainprep/${MODE}/splits
OUTPUT_BASE=$DATA/proteinmpnn_input/${MODE}

mkdir -p logs

module load python/3.11.9

echo "================================================"
echo "MODE  : $MODE"
echo "FOLD  : $FOLD"
echo "START : $(date)"
echo "================================================"

# ── HLA split ────────────────────────────────────────────────────
python $SCRIPTS/prepare_pmhc_data.py \
    --splits_dir $SPLITS_BASE/hla \
    --output_dir $OUTPUT_BASE/hla \
    --split_mode hla \
    --fold $FOLD

# ── Anchor split (fold=0 is only used for HLA test set, skip here) ──
if [ "$FOLD" -ne 0 ]; then
    python $SCRIPTS/prepare_pmhc_data.py \
        --splits_dir $SPLITS_BASE/anchor \
        --output_dir $OUTPUT_BASE/anchor \
        --split_mode anchor \
        --fold $FOLD
fi

echo "================================================"
echo "DONE  : $(date)"
echo "================================================"