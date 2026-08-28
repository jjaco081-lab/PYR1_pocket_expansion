#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 12:00:00
#SBATCH -J belddg
#SBATCH -a 0-19
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/belddg_%a.log
#
# 154b -- 380 single substitutions at Beltran's 20 design positions x
# (mutant + paired wild-type) x 3 Cartesian relax replicates = 2,280 relaxes.
# CPU only on cutlerlab, so neither GPU allowance is touched and this runs
# alongside the tractability array.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/154_beltran_ddg.py --chunk "$SLURM_ARRAY_TASK_ID" --nchunks 20 --nrep 3
