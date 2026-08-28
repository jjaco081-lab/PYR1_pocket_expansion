#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 0:40:00
#SBATCH -J belsmoke
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/belddg_smoke.log
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    scripts/154_beltran_ddg.py --chunk 0 --nchunks 190 --nrep 1
