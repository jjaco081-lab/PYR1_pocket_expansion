#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 10:00:00
#SBATCH -J parkcond
#SBATCH -a 0-19
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/parkcond_%a.log
set -uo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/160_conditional_park.py --chunk "$SLURM_ARRAY_TASK_ID" --nchunks 20 --nrep 3
