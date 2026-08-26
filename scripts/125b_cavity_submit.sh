#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 12:00:00
#SBATCH -J cavmatch
#SBATCH -a 0-39
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/cav_%x_%a.log
#
# 125b -- the §68e cavity-match test, two arms:
#   coumarin  target 163 A^3, BELOW PYR1's ~166 A^3 -- where ddG already works
#   pfas      target 225 A^3, ABOVE it -- the pocket-EXPANSION regime we care about
# Cartesian relax + cavity volume per variant. CPU only, cutlerlab, no GPU touched.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/125_cavity_match.py --klass "$KLASS" --chunk "$SLURM_ARRAY_TASK_ID" \
    --nchunks 40 --nrep 3
