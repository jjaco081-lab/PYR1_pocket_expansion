#!/bin/bash
#SBATCH --job-name=desrec
#SBATCH --partition=cutlerlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=0-8
#SBATCH --output=logs/desrec_%a.log
# 3 systems x 3 arms. CPU only -- the GPUs stay with RFd3.
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
SYS=(aba mandi win); ARM=(fastdesign hbnet_terms hbnet_mover)
s=${SYS[$((SLURM_ARRAY_TASK_ID / 3))]}
a=${ARM[$((SLURM_ARRAY_TASK_ID % 3))]}
echo "system=$s arm=$a"
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python -u \
  scripts/210_design_recovery.py --system "$s" --arm "$a" --nrep 1
test -s results/design_recovery/${s}_${a}.json || { echo "NO OUTPUT"; exit 1; }
echo "OK $s $a"
