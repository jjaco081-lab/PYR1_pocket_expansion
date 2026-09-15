#!/bin/bash
#SBATCH --job-name=rfd3s1
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu05,gpu11,gpu13,gpu14
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/rfd3s1_%j.log
# Stage-1 RFd3 generation. $1=arm $2=n $3=master_seed
# Partition and time come from the sbatch command line so the same script
# serves both the `gpu` (7-day) and `short_gpu` (2 h) submissions.
# ~47 s per design measured on Stage 0 (10 designs in 7:47).
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
PY=/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python
$PY -u scripts/193_rfd3_gen.py --arm "$1" --n "$2" --master-seed "$3"
# assert on the RESULT, not the exit code: count structures actually written
OUT=results/rfd3/arm${1}_s${3}
N=$(find "$OUT" -name '*.cif.gz' | wc -l)
echo "RESULT: $N structures in $OUT (requested $2)"
test "$N" -gt 0 || { echo "NO STRUCTURES WRITTEN"; exit 1; }
