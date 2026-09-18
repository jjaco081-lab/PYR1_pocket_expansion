#!/bin/bash
#SBATCH --job-name=casc3
#SBATCH --partition=cutlerlab
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/casc_stage3_%j.log
# Stage 3: the four arms that finally have ABA IN THE POCKET (README 129).
# Every earlier design was scaffolded around an empty pocket.
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
PY=/bigdata/cutlerlab/jjaco081/conda_envs/foundry/bin/python
$PY -u scripts/198_cascade2.py --procs 64 --out results/rfd3/cascade3.csv \
  --arms results/rfd3/arm9aba_s9001 results/rfd3/arm10ext1_s10011 \
         results/rfd3/arm10ext2_s10021 results/rfd3/arm10ext3_s10031
test -s results/rfd3/cascade3.csv || { echo "NO CSV"; exit 1; }
echo "RESULT: $(( $(wc -l < results/rfd3/cascade3.csv) - 1 )) designs measured"
$PY -u scripts/217_truncate_cascade.py --csv results/rfd3/cascade3.csv \
  --out results/rfd3/truncated3.csv 2>/dev/null || \
  $PY -u scripts/220_core_truncate.py --csv results/rfd3/cascade3.csv \
  --out results/rfd3/truncated3.csv
echo "done"
