#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 0:50:00
#SBATCH -J dsm2smk
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/dsm2_smoke.log
set -uo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/159_dsm_doubles.py --chunk 0 --nchunks 3000 --nrep 2
