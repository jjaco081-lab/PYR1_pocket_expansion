#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 8:00:00
#SBATCH -J viab
#SBATCH -a 0-19
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/viab_%a.log
#
# 122b -- protein-only viability scoring of 678 variants (78 real coumarin
# sensors + 300 random members of Tian's library + 300 random combinations from
# the full DSM-Hao menu at the same 11 positions). CPU only, no ligand, so it
# runs on cutlerlab and does not touch either GPU allowance.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/122_combination_viability.py --chunk "$SLURM_ARRAY_TASK_ID" --nchunks 20
