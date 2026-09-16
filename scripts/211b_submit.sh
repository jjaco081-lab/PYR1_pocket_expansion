#!/bin/bash
#SBATCH --job-name=posscan
#SBATCH --partition=cutlerlab
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=12:00:00
#SBATCH --array=0-23
#SBATCH --output=logs/posscan_%a.log
# 3 systems x 2 arms x 4 chunks = 24 tasks. CPU only.
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
SYS=(aba mandi win); ARM=(ref2015 hbnet_terms)
t=$SLURM_ARRAY_TASK_ID
s=${SYS[$((t / 8))]}; a=${ARM[$(((t % 8) / 4))]}; c=$((t % 4))
echo "system=$s arm=$a chunk=$c"
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python -u \
  scripts/211_position_scan.py --system "$s" --arm "$a" --chunk "$c" --nchunks 4
test -s results/position_scan/${s}_${a}_${c}.json || { echo "NO OUTPUT"; exit 1; }
echo "OK $s $a $c"
