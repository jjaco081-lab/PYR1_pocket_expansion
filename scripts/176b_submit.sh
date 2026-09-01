#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 16:00:00
#SBATCH -J dsmterm
#SBATCH -a 0-39
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/dsmterm_%a.log
set -uo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/176_dsm_terms.py --chunk "$SLURM_ARRAY_TASK_ID" --nchunks 40 --nrep 2
