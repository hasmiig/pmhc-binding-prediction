#!/bin/bash
#SBATCH --job-name=pmgen_binder_new
#SBATCH --output=logs/pmgen_%x_%j.out
#SBATCH --error=logs/pmgen_%x_%j.err
#SBATCH --partition=scc-gpu
#SBATCH --time=48:00:00
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --constraint=inet
 
CHUNK=$1
cd /projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig/PMGen/

source /projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig/pmhc-binding-prediction/scripts/utils/setup_env.sh
 
python /projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig/PMGen/run_PMGen.py \
    --initial_guess \
    --multiple_anchors \
    --df $CHUNK \
    --output_dir "/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/data/pmgen_outputs/binder_new/$(basename $CHUNK .tsv)/" \
    --mode wrapper \
    --run parallel \
    --max_cores 10 \
    --no_netmhcpan \
    --best_structures