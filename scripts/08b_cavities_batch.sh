#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 08:00:00
#SBATCH -J cavities
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/08_cavities_batch.log
#
# 08b_cavities_batch.sh -- cavity measurement for all 266 helix-grip hits, as a
# BATCH job rather than an interactive one.
#
# WHY BATCH. The measurement itself is light -- 135 MB peak RSS, measured, so the
# script's old "must run under SLURM" warning was wrong about the reason. What it
# does need is to OUTLIVE THE SESSION: an interactive allocation kills every process
# it owns when it ends, nohup included, and this takes ~2 h for 266 structures.
#
# Results now stream to homolog_cavities_partial.csv as each structure is measured,
# so even a killed run keeps everything it had finished. Output goes to
# homolog_cavities_full.csv, NOT the original homolog_cavities.csv, so a partial run
# can never clobber the good 147-row file.
set -uo pipefail
P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
PY=/bigdata/cutlerlab/jjaco081/conda_envs/esmfold2/bin/python
echo "=== $(date) cavity measurement, 266 hits ==="
$PY "$P/scripts/08_homolog_cavities.py" --stage compute --min-tm 0.5 --max-n 300
rc=$?
echo "CAVITIES_DONE rc=$rc  rows=$(tail -n +2 $P/results/homolog_cavities/homolog_cavities_partial.csv 2>/dev/null | wc -l)"
date
exit $rc
