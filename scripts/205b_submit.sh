#!/bin/bash
#SBATCH --job-name=predcav
#SBATCH --partition=cutlerlab
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=logs/predcav_%j.log
# Cavity measurement on 1,810 EXISTING Boltz-2 co-folds. No GPU: the structures
# were computed in the tractability run and are on disk.
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/foundry/bin/python -u scripts/205_predicted_cavity.py --procs 64
test -s results/predicted_cavity/predicted_cavity.csv || { echo "NO CSV"; exit 1; }
echo "RESULT: $(( $(wc -l < results/predicted_cavity/predicted_cavity.csv) - 1 )) rows"
