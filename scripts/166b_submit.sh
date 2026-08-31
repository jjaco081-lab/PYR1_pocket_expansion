#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 8:00:00
#SBATCH -J agro
#SBATCH -a 0-29
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/agro_%a.log
set -uo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion || exit 1
NAMES=(benzothiadiazole benoxacor fludioxonil)
C=${NAMES[$((SLURM_ARRAY_TASK_ID / 10))]}
K=$((SLURM_ARRAY_TASK_ID % 10))
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/166_agro_scan.py --compound "$C" --chunk "$K" --nchunks 10 --nrep 3
