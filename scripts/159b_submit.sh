#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 14:00:00
#SBATCH -J dsm2
#SBATCH -a 0-39
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/dsm2_%a.log
set -uo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/159_dsm_doubles.py --chunk "$SLURM_ARRAY_TASK_ID" --nchunks 40 --nrep 2
