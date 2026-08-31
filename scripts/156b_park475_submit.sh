#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 8:00:00
#SBATCH -J park475
#SBATCH -a 0-19
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/park475_%a.log
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/156_park475.py --chunk "$SLURM_ARRAY_TASK_ID" --nchunks 20 --nrep 3
