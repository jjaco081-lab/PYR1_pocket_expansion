#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 10:00:00
#SBATCH -J ladder
#SBATCH -a 0-12
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/ladder_%a.log
set -uo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/167_pose_ladder.py --chunk "$SLURM_ARRAY_TASK_ID" --nchunks 13 --nrep 3
