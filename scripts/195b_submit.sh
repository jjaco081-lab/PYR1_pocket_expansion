#!/bin/bash
#SBATCH --job-name=armcmp
#SBATCH --partition=cutlerlab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=logs/armcmp_%j.log
# Stage-0 arm comparison. ~60 s per design, 30 designs, so ~35 min of compute --
# but it ran on the login node under three-way contention and dragged past an
# hour, and /scratch is node-local so nothing may be written there.
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
PY=/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python
$PY -u scripts/195_arm_compare.py \
    results/rfd3/batch03 results/rfd3/arm0a_s2026 results/rfd3/arm0b_s2027
# assert on the RESULT, not the exit code
test -s results/rfd3/arm_compare.json || { echo "NO SUMMARY WRITTEN"; exit 1; }
echo "OK: $(wc -c < results/rfd3/arm_compare.json) bytes of summary"
