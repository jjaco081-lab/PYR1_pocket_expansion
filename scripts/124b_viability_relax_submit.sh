#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 10:00:00
#SBATCH -J viabrlx
#SBATCH -a 0-39
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/viabrlx_%a.log
#
# 124b -- the FastRelax version of the §66a viability test. 670 variants x
# (mutant + paired wild-type) x 3 Cartesian relax replicates = ~4,000 relaxes,
# ~21 CPU-hours total. CPU only, so neither GPU allowance is touched.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/124_viability_relax.py --chunk "$SLURM_ARRAY_TASK_ID" --nchunks 40 --nrep 3
