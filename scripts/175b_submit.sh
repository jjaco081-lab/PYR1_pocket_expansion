#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 10:00:00
#SBATCH -J terms
#SBATCH -a 0-39
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/terms_%a.log
set -uo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion || exit 1
NAMES=(mandipropamid_crystal fludioxonil benzothiadiazole benoxacor)
C=${NAMES[$((SLURM_ARRAY_TASK_ID / 10))]}
K=$((SLURM_ARRAY_TASK_ID % 10))
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/175_term_decomposition.py --compound "$C" --chunk "$K" --nchunks 10 --nrep 3
