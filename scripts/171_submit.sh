#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 8:00:00
#SBATCH -J agro2
#SBATCH -a 0-29
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/agro2_%a.log
# Three more arms of the 166 scan, all against Park's 475 in the K59R background:
#   mandipropamid  PREDICTED pose -- isolates crystal vs predicted for one ligand
#   azoxystrobin   30 heavy atoms, 0 of 475 responders -- specificity
#   lufenuron      30 heavy atoms, 0 of 475 responders -- specificity
set -uo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion || exit 1
NAMES=(mandipropamid azoxystrobin lufenuron)
C=${NAMES[$((SLURM_ARRAY_TASK_ID / 10))]}
K=$((SLURM_ARRAY_TASK_ID % 10))
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/166_agro_scan.py --compound "$C" --chunk "$K" --nchunks 10 --nrep 3
