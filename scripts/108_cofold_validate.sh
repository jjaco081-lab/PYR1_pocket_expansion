#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 8:00:00
#SBATCH -J cofold
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/cofold.log
#
# 108_cofold_validate.sh -- can co-folding reproduce a KNOWN complex?
#
# WHY THIS GATES HALF THE DESIGN STACK
# §49c rated pose-derived methods DO-NOT-ATTEMPT, because §46b measured docking at
# 8.33 A from the crystal pose into a pocket that does not exist yet. But an
# initial HIT changes the problem: it hands you a sequence KNOWN to bind, so the
# task becomes modelling a complex that exists rather than designing one. That is
# what co-folding is for.
#
# Whether it actually works here has never been tested. If Boltz reproduces 4WVO,
# structure becomes usable for the steric positions (159/160/83/87/89) and the
# polar ones (120/163) that §49b left open because H-bond geometry is directional.
# If it does not, the design stack stays pose-free and those positions fall back
# to frequency priors.
#
# FOUR RUNS, and the controls are the point
#   quad_mandi_binary   PYR1^MANDI + mandipropamid          <- the design scenario
#   quad_mandi_ternary  + HAB1                              <- matches the crystal
#   wt_mandi_binary     WT + mandipropamid                  <- NEGATIVE control
#   wt_aba_binary       WT + ABA                            <- POSITIVE control (3QN1)
#
# The negative control is not optional. 4WVO is a TERNARY complex and the
# gate-latch only closes on HAB1 binding, so a binary prediction may not reach the
# closed state at all. And if WT + mandipropamid scores as confidently as the
# quadruple does, the method has no discriminative power regardless of how good
# the quadruple's pose looks -- WT is known NOT to bind mandipropamid, which is
# the entire reason the quadruple was evolved.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
BOLTZ=/bigdata/cutlerlab/jjaco081/conda_envs/boltz2/bin/boltz
OUT=$ROOT/data/cofold/out
CACHE=/bigdata/cutlerlab/jjaco081/conda_envs/boltz2/.boltz_cache
mkdir -p "$OUT" "$CACHE"

echo "=== co-folding validation $(date -Is)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null

for y in quad_mandi_binary quad_mandi_ternary wt_mandi_binary wt_aba_binary; do
    f=$ROOT/data/cofold/$y.yaml
    if [[ -d $OUT/boltz_results_$y ]]; then echo "  $y: exists, skip"; continue; fi
    echo "--- $y  $(date -Is)"
    $BOLTZ predict "$f" \
        --out_dir "$OUT" \
        --cache "$CACHE" \
        --accelerator gpu --devices 1 \
        --recycling_steps 3 --diffusion_samples 5 \
        --output_format mmcif \
        --use_msa_server \
        --override 2>&1 | tail -15
    echo "    exit=$?"
done
echo "=== done $(date -Is)"
ls -R "$OUT" 2>/dev/null | head -40
