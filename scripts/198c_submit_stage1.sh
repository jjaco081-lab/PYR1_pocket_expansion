#!/bin/bash
#SBATCH --job-name=casc1
#SBATCH --partition=cutlerlab
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/casc_stage1_%j.log
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
PY=/bigdata/cutlerlab/jjaco081/conda_envs/foundry/bin/python
$PY -u scripts/198_cascade2.py --procs 64 --out results/rfd3/cascade2.csv \
  --arms results/rfd3/batch03 results/rfd3/arm0a_s2026 results/rfd3/arm0b_s2027 \
         results/rfd3/arm2noanch_s3002 results/rfd3/arm3b3rasa_s3003 \
         results/rfd3/arm1slack_s3001 results/rfd3/arm4rep_s3004
test -s results/rfd3/cascade2.csv || { echo "NO CSV"; exit 1; }
echo "RESULT: $(( $(wc -l < results/rfd3/cascade2.csv) - 1 )) designs measured"
$PY -u scripts/199_cascade_report.py
