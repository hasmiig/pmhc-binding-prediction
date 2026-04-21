#!/bin/bash
#SBATCH --job-name=pmgen_binder
#SBATCH --output=logs/pmgen_%x_%j.out
#SBATCH --error=logs/pmgen_%x_%j.err
#SBATCH --partition=scc-gpu
#SBATCH --time=48:00:00
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --constraint=inet
 
CHUNK=$1
 
source /projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig/scripts/setup_env.sh
 
python /projects/scc/MPG/MGMN/scc_mgmn_soeding/dir.project/hasmig/PMGen/run_PMGen.py \
    --initial_guess \
    --multiple_anchors \
    --df $CHUNK \
    --output_dir outputs/$(basename $CHUNK .tsv)/ \
    --mode wrapper \
    --run parallel \
    --max_cores 10 \
    --no_netmhcpan \
    --only_last_anchor \