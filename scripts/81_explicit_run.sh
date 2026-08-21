#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu13,gpu14
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 1:58:00
#SBATCH -J exp_ref
#SBATCH -a 0-5%4
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/exp_ref_%a.log
#
# 81_explicit_run.sh -- the two MM-GBSA reference systems, in explicit water, x3 seeds.
#
# WHAT THIS IS TESTING (README 40, 41)
#   In implicit solvent the WT+ABA reference never converged: -32.74 at 50 ps,
#   -25.00 at 250 ps, still drifting +2.62 inside the 250 ps window, with ABA
#   sliding 2-4 A out of its crystal pose. Because every ddG in the ABA arm
#   subtracts that one number, the drift propagated into all 22 of them at full
#   strength and inverted the pre-registered K59R verdict.
#
#   Hypothesis under test: with real waters present the crystallographic pose is a
#   minimum of the simulated Hamiltonian too, so the reference stops moving.
#
#   FALSIFIABLE. If ABA still slides 2-4 A and the score still drifts by several
#   kcal/mol across the last half, explicit solvent is NOT the fix and MM-GBSA
#   should be retired for this ligand rather than tuned further.
#
# WHY THREE SEEDS AND NOT ONE LONGER RUN
#   Twice now a tight WITHIN-run error bar has hidden the real problem: first the
#   8 byte-identical repacks (sd 0.00), then the 250 ps block SE, which stayed at
#   0.3-1.3 while the mean marched 7.75. An error bar cannot see a drift it is
#   centred on. Three independent velocity seeds give a BETWEEN-run spread, which
#   is the honest error bar for a quantity this drifty, and it is the check that
#   would have caught this at 50 ps.
#
# WHAT IS DELIBERATELY UNCHANGED
#   Scoring is still MM-GBSA (igb=8) on water-stripped frames -- the standard
#   protocol, and it keeps these numbers commensurable with the 46 implicit runs.
#   Exactly one variable moves: the solvent the SAMPLING happened in.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22_mpi_cuda >/dev/null 2>&1
source scripts/lib_mdinputs.sh

# 10 ns, chosen to FIT the 2 h short_gpu wall rather than to be sufficient a priori
# (user, 2026-08-21). This is not a commitment: progress is derived from the frames
# actually on disk (81b_frames.py), so resubmitting with PROD_NS=20 CONTINUES from
# where the last wave stopped instead of restarting. Pick the floor, extend if the
# reference has not plateaued.
PROD_NS=${PROD_NS:-10}

# short_gpu caps at 2 h, so production runs in CHUNKS under a wall-clock deadline
# rather than as one long leg. The `gpu` partition would have been simpler, but the
# cutlerlab association allows gres/gpu=4 there and the factorial (27547457/27547617)
# already holds 3 -- 6 more would have queued behind it. short_gpu carries a
# SEPARATE gres/gpu=4 allowance, so these run alongside rather than competing.
CHUNK_NS=${CHUNK_NS:-1}          # small enough that even a slow card finishes one
BUDGET_S=${BUDGET_S:-5700}       # 95 min of the 118 min wall; the rest is margin
T0=$SECONDS

UNITS=(aba_WT_s0 aba_WT_s1 aba_WT_s2 mandi_WT_s0 mandi_WT_s1 mandi_WT_s2)
U=${UNITS[${SLURM_ARRAY_TASK_ID:-0}]}
D=$ROOT/data/mmgbsa_explicit/$U
cd "$D" || { echo "no such system: $D"; exit 1; }

echo "=== $U on $(hostname) $(date -Is)"
nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1

TOP=$(topology_for_run "$D") || exit 1
echo "topology: $TOP"

run_stage () {   # $1 = stage name, $2 = input restart
    local nm=$1 prev=$2
    if [[ -f ${nm}.rst7 ]]; then echo "  $nm: done, skip"; return 0; fi
    echo "  $nm ..."
    pmemd.cuda -O -i ${nm}.in -p "$TOP" -c "$prev" \
        -o ${nm}.out -r ${nm}.rst7 -x ${nm}.nc -inf ${nm}.mdinfo -ref "$prev" \
        || { echo "  !! $nm FAILED"; record_seed "$nm"; return 2; }
    record_seed "$nm"
}

write_md_inputs
run_stage min1 system.inpcrd || exit 2
run_stage min2 min1.rst7     || exit 2
run_stage heat min2.rst7     || exit 2
run_stage eq1  heat.rst7     || exit 2
run_stage eq2  eq1.rst7      || exit 2

# ---- production, chunked under a wall-clock deadline ----
ns_done () { python3 "$ROOT/scripts/81b_frames.py" | \
             xargs -I{} python3 -c "print({} * $FRAME_PS / 1000.0)"; }

NS_DONE=$(ns_done)
echo "production so far: $NS_DONE of $PROD_NS ns"
START=eq2.rst7; [[ -f prod.rst7 ]] && START=prod.rst7
LAST_S=0

while :; do
    if python3 -c "import sys; sys.exit(0 if $PROD_NS - $NS_DONE <= 0.001 else 1)"; then
        echo "=== $U production COMPLETE ($NS_DONE ns)"; break
    fi
    ELAPSED=$((SECONDS - T0))
    # only start a chunk if the LAST one's measured duration fits in what is left,
    # so a slow card stops cleanly instead of being killed mid-write at the wall
    if [[ $LAST_S -gt 0 && $((ELAPSED + LAST_S + 120)) -gt $BUDGET_S ]]; then
        echo "=== stopping cleanly: ${ELAPSED}s used, a chunk costs ~${LAST_S}s of"
        echo "    a ${BUDGET_S}s budget. $NS_DONE / $PROD_NS ns done; resubmit to continue."
        break
    fi
    THIS_NS=$(python3 -c "print(min($CHUNK_NS, round($PROD_NS - $NS_DONE, 3)))")
    THIS_STEPS=$(python3 -c "print(int($THIS_NS * 1000 / $DT))")
    write_prod_input "$THIS_STEPS"
    SEG=prod_seg$(date +%s)
    CS=$SECONDS
    pmemd.cuda -O -i prod.in -p "$TOP" -c "$START" \
        -o ${SEG}.out -r prod.rst7 -x ${SEG}.nc -inf prod.mdinfo \
        || { echo "  !! chunk $SEG failed"; exit 2; }
    LAST_S=$((SECONDS - CS)); record_seed "$SEG"; START=prod.rst7
    NS_DONE=$(ns_done)
    echo "  +${THIS_NS} ns in ${LAST_S}s -> $NS_DONE / $PROD_NS ns"
done

echo "=== $U done $(date -Is)"
ls -la prod*.nc 2>/dev/null | tail -3
