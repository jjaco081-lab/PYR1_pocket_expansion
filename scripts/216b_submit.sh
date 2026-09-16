#!/bin/bash
#SBATCH --job-name=cavv2
#SBATCH --partition=cutlerlab
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --output=logs/cavity_v2_%j.log
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/foundry/bin/python -u scripts/216_cavity_recompute.py 48
test -s results/homolog_cavities/cavity_v2.json || { echo "NO OUTPUT"; exit 1; }
echo OK
