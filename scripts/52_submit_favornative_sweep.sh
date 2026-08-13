#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 04:00:00
#SBATCH -J pe_s1_fnsweep
#SBATCH -a 0-7
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/stage1_fnsweep_%a.log
#
# 52_submit_favornative_sweep.sh -- calibrate the favor-native bonus on the NULL arm.
#
# WHY
#   Job 27421344 ran with no favor-native term. Its ABA arm -- WT protein with its
#   native ligand, correct answer ZERO mutations -- kept WT at only 4 of 15 designable
#   positions and deleted the K59 carboxylate salt bridge in 100% of trajectories.
#   That makes every position the null also mutates uninterpretable: "no ligand-
#   conditional signal" is not separable from "the protocol cannot hold a native
#   contact". K59R and V81I are exactly those positions.
#
# THE CALIBRATION, AND WHY IT STAYS HONEST
#   The bonus weight is chosen ONLY on the ABA null arm, by a criterion fixed before
#   looking at any output: pick the SMALLEST weight whose null keeps WT at >= 70% of
#   the 15 designable positions AND restores K59 (>= 50% Lys). The mandipropamid arms
#   are not run in this sweep and play no part in choosing the weight, so the weight
#   cannot be tuned toward the answer we are testing for.
#
#   Smallest-that-passes matters: a large bonus trivially freezes the sequence and
#   would produce a perfect null with zero design signal in either arm. The null
#   getting quieter is only meaningful if design can still move. The mandipropamid
#   arms are NOT in this sweep, so 53 cannot check that directly; it instead reports
#   distinct sequences and mean entropy within the null, which detects a bonus large
#   enough to freeze the packer. The real floor check -- that the mandipropamid arm
#   still mutates F108 -- comes after the full six-arm re-run at the chosen weight.
#
# SHAPE
#   4 weights x 2 alphabets x 10 trajectories = 80 trajectories.
#   Units 20 (wt_aba/dsm_hao block 0) and 25 (wt_aba/free block 0) are the ABA arm;
#   see the arm_idx mapping in 50_stage1_rosetta.py (units 20-29 are wt_aba).
#   ~5.1 min/trajectory measured at favor_native=0, so ~51 min/task, ~6.8 CPU-hours.
#   Each weight writes to its own outdir so job 27421344's results are untouched.
#
# NEXT
#   bash scripts/52_submit_favornative_sweep.sh   (or sbatch)
#   python scripts/53_favornative_pick.py         # applies the criterion above
set -euo pipefail

ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
PY=/bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python

WEIGHTS=(0.25 0.5 1.0 1.5)
UNITS=(20 25)                     # wt_aba dsm_hao b0, wt_aba free b0

w=${WEIGHTS[$((SLURM_ARRAY_TASK_ID / 2))]}
u=${UNITS[$((SLURM_ARRAY_TASK_ID % 2))]}

OUTDIR=$ROOT/results/stage1_rosetta_fnsweep/w${w}
mkdir -p "$OUTDIR"

echo "task $SLURM_ARRAY_TASK_ID: favor_native=$w unit=$u -> $OUTDIR"
$PY $ROOT/scripts/50_stage1_rosetta.py \
    --unit "$u" \
    --n-per-unit 10 \
    --favor-native "$w" \
    --outdir "$OUTDIR"
