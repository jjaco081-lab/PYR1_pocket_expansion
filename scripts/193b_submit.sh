#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 1:58:00
#SBATCH -J rfdgen
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/rfdgen_%j.log
set -uo pipefail
module load gcc/12.2.0 2>/dev/null; module load cuda/12.1 2>/dev/null
source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda deactivate 2>/dev/null
conda activate /bigdata/cutlerlab/jjaco081/conda_envs/foundry
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion || exit 1
python scripts/193_rfd3_gen.py --arm "$1" --n "$2" --master-seed "$3"
