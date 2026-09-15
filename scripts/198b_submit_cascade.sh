#!/bin/bash
#SBATCH --job-name=cascade2
#SBATCH --partition=cutlerlab
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/cascade2_%j.log
# Cascade v2 over every arm generated so far. ~3.7 min per design single-core,
# so 64 procs clears ~1000 designs in well under an hour of wall time; the
# generous walltime is for the case where generation is still finishing.
# The CSV is resumable -- rerunning skips designs already present.
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
PY=/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python

# wild-type PYR1 must survive the cascade before it is applied to any design
$PY -u scripts/198_cascade2.py --calibrate

$PY -u scripts/198_cascade2.py --procs 64 \
  --out results/rfd3/cascade2.csv \
  --arms results/rfd3/batch03 \
         results/rfd3/arm0a_s2026  results/rfd3/arm0b_s2027 \
         results/rfd3/arm1slack_s3001 results/rfd3/arm4rep_s3004 \
         results/rfd3/arm2noanch_s3002 results/rfd3/arm3b3rasa_s3003

# assert on the RESULT
test -s results/rfd3/cascade2.csv || { echo "NO CSV"; exit 1; }
echo "RESULT: $(( $(wc -l < results/rfd3/cascade2.csv) - 1 )) designs measured"
$PY -u scripts/199_cascade_report.py
