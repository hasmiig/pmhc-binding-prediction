module purge
module load gcc/13.2.0-nvptx
module load cuda/12.6.2
module load cudnn/9.8.0.87-12
source $(conda info --base)/etc/profile.d/conda.sh
conda activate PMGen