#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 0:40:00
#SBATCH -J p475smk
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/park475_smoke.log
set -uo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/156_park475.py --chunk 0 --nchunks 158 --nrep 3
