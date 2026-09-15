#!/bin/bash
#SBATCH --job-name=trunc
#SBATCH --partition=cutlerlab
#SBATCH --cpus-per-task=40
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/trunc_%j.log
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/foundry/bin/python -u scripts/203_truncation_rescue.py 40
test -s results/rfd3/truncation_rescue.json || { echo "NO OUTPUT"; exit 1; }
echo "OK"
