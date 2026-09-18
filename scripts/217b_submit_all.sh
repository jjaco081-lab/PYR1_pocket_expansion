#!/bin/bash
#SBATCH --job-name=trunc_all
#SBATCH --partition=cutlerlab
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/trunc_all_%j.log
# 217 globs results/rfd3/arm* so this now covers the four Stage 3 ABA arms too
# (2,820 earlier designs + 160 Stage 3). casc3 measured them; nothing truncated them.
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/foundry/bin/python -u \
  scripts/217_truncate_cascade.py --procs 64 --out results/rfd3/truncated_all.csv
n=$(( $(wc -l < results/rfd3/truncated_all.csv) - 1 ))
echo "RESULT: $n designs truncated"
test "$n" -gt 2900
