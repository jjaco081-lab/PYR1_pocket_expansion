#!/bin/bash
#SBATCH --job-name=trunc2
#SBATCH --partition=cutlerlab
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/truncate_%j.log
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/foundry/bin/python -u scripts/217_truncate_cascade.py --procs 64
test -s results/rfd3/truncated.csv || { echo "NO OUTPUT"; exit 1; }
echo "RESULT: $(( $(wc -l < results/rfd3/truncated.csv) - 1 )) designs"
