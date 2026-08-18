#!/bin/bash
#SBATCH -p preempt_gpu
#SBATCH -A preempt
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 7-00:00:00
#SBATCH -J fact_md
#SBATCH -a 0-11
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/fact_md_%a.log
#
# 70_md_run_factorial.sh -- the conformation x occupancy factorial, all four cells
# on ONE build protocol. README 30.
#
# THE DESIGN
#                     apo                      + ABA
#     open      S1_apo_open   (tasks 0-2)   S10_holo_open  (tasks 9-11)
#     closed    S9_apo_closed (tasks 6-8)   S2_holo_closed (tasks 3-5)
#
#   README 19b set this up and only the diagonal was ever run, so conformation and
#   occupancy have never been separated (24e). The two off-diagonal cells are the
#   informative ones: each is PYR1 in the "wrong" form.
#
# WHY ALL FOUR ARE RE-RUN, INCLUDING THE TWO WE ALREADY HAVE
#   The existing S1/S2 trajectories are NaCl, 178 residues, blanket HIE, crystal
#   chain breaks. The old-tree S9 was built later and is KCl -- so pairing it with
#   the old S2 would confound ligand removal with a cation swap, in exactly the
#   comparison the factorial exists to make, and Na+ binds carboxylates more
#   strongly than K+. A factorial whose cells differ in the build is not a
#   factorial. The md/ trajectories are NOT discarded: they remain the evidence
#   behind README 24 and 29 on their own protocol.
#
# WHY S10 IS THE INTERESTING ONE
#   Gate closure on ligand binding is the physiological direction, so it is the one
#   arm where the barrier this project has never crossed (README 24: no transition
#   in 1.8 us of aggregate sampling) is plausibly downhill. ABA was transplanted
#   from the closed structure by core superposition, NOT docked, so it starts in
#   its known correct position; script 69 froze the open backbone bit-for-bit while
#   doing it, because a relaxation toward closed would pre-bias the answer.
#
# PRE-REGISTERED READING (README 30c), repeated here so it travels with the job
#   * S10 has three outcomes and only two are informative: the gate CLOSES (the
#     ligand drives closure and MD can see it), or ABA LEAVES (the open pocket does
#     not retain it). Nothing happening is the uninformative case -- 300 ns against
#     a barrier -- and must not be reported as "ABA does not close the gate".
#   * S9 vs S2 is the ligand-removal contrast WITHOUT the dimer clamp that
#     confounds S3 protomer B (29c: dimerisation alone changes gate RMSF 2.9x).
#   * the primary observables are those already pre-registered in 19d. Nothing new
#     is added here, so the answer cannot be shopped for.
#
# PREEMPTION: PreemptMode=REQUEUE on this partition, so a preempted task returns to
# the queue and re-runs from the top. Every resume writes a fresh prod_cont_NNN.nc
# and recomputes the remaining steps from the restart file's own clock, so a
# requeue never runs a further full 300 ns. S4 rep0 needed 23 segments to finish
# this way. --gres=gpu:1 rather than a named card: preempt_gpu holds a100, ada6000,
# blackwell6000, h100, p100 and k80, and pinning to one type only lengthens the
# wait. (k80 is too old for this build; if a task lands on one, cancel that task.)
set -uo pipefail
module load amber/22_mpi_cuda >/dev/null 2>&1

P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
MD=$P/data/md191
PROD_NS=300
# sourced by ABSOLUTE path: SLURM copies this script to /var/spool/slurmd, so
# $(dirname $0) does not point into the repo and the source would fail silently
source "$P/scripts/lib_mdinputs.sh"
md_settings_check
TARGET_PS=$(( EQUIL_PS + PROD_NS * 1000 ))

SYSTEMS=(S1_apo_open S2_holo_closed S9_apo_closed S10_holo_open)
SYS=${SYSTEMS[$(( SLURM_ARRAY_TASK_ID / 3 ))]}
REP=$(( SLURM_ARRAY_TASK_ID % 3 ))
TOPDIR=$MD/$SYS
D=$TOPDIR/rep$REP

echo "=== $(date)  task $SLURM_ARRAY_TASK_ID = $SYS rep$REP  job $SLURM_JOB_ID  restarts=${SLURM_RESTART_COUNT:-0} ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

TOP=$(topology_for_run "$TOPDIR") || exit 1
CRD=$TOPDIR/system.inpcrd
[[ -s "$TOP" && -s "$CRD" ]] || { echo "missing topology/coords in $TOPDIR -- run 68 and 69b first"; exit 1; }

rst_time () { ncdump -v time "$1" 2>/dev/null | awk '/^ time =/{gsub(/[^0-9.]/,"",$3); print $3; exit}'; }

run_stage () {
    local nm=$1 prev=$2
    if [[ -s ${nm}.rst7 ]]; then echo "  $nm done, skipping"; return 0; fi
    echo "  running $nm ..."
    local ref=()
    [[ "${3:-}" == "ref" ]] && ref=(-ref "$prev")
    pmemd.cuda -O -i ${nm}.in -p "$TOP" -c "$prev" \
               -o ${nm}.out -r ${nm}.rst7 -x ${nm}.nc -inf ${nm}.info "${ref[@]}" \
        || { echo "  !! $nm FAILED"; record_seed "$nm"; return 2; }
    record_seed "$nm"
}

mkdir -p "$D"; cd "$D" || exit 1
write_md_inputs

run_stage min1 "$CRD" ref || exit 1
run_stage min2 min1.rst7  || exit 1

# The AMBER-minimised coordinates are the real MD starting structure, so the
# ligand-contact check belongs here rather than at the Rosetta stage: script 69
# left S10's tightest contact at 2.72 A (the K59-carboxylate salt bridge) and a
# continuous minimiser is what decides whether that is comfortable. Reported, not
# asserted -- a first-frame geometry is not grounds to kill a job, but it must be
# in the log where anyone reading the result can find it.
if [[ "$SYS" == S2_holo_closed || "$SYS" == S10_holo_open ]] && [[ -s min2.rst7 ]]; then
    cat > _ligchk.in <<CPP
parm $TOP
trajin min2.rst7
strip :WAT,K+,Cl-,Na+
distance k59_aba :59@NZ :A8S&!@H= out _k59_aba.dat
run
quit
CPP
    cpptraj -i _ligchk.in > _ligchk.log 2>&1 || true
    echo "  post-minimisation K59 NZ to ABA centre: $(awk 'NR==2{printf "%.2f A", $2}' _k59_aba.dat 2>/dev/null || echo unknown)"
fi

run_stage heat min2.rst7 ref || exit 1
run_stage eq1  heat.rst7 ref || exit 1
run_stage eq2  eq1.rst7      || exit 1

START=eq2.rst7
NOW=$EQUIL_PS
if [[ -s prod.rst7 ]]; then
    T=$(rst_time prod.rst7)
    if [[ -z "$T" && -s prod_backup.rst7 ]]; then
        T=$(rst_time prod_backup.rst7)
        [[ -n "$T" ]] && { cp prod_backup.rst7 prod.rst7; echo "  recovered from backup at $T ps"; }
    fi
    if [[ -n "$T" ]]; then START=prod_resume.rst7; NOW=$T; cp prod.rst7 "$START"; fi
fi
REMAIN_PS=$(awk -v a=$TARGET_PS -v b=$NOW 'BEGIN{printf "%.0f", a-b}')
if [[ "$REMAIN_PS" -le 0 ]]; then
    echo "FACT_RUN_DONE sys=$SYS rep=$REP rc=0 time=$NOW/$TARGET_PS"
    exit 0
fi
NSTEPS=$(( REMAIN_PS * STEPS_PER_PS ))
echo "  production: at $NOW ps -> $REMAIN_PS ps remaining ($NSTEPS steps)"
write_prod_input "$NSTEPS"

if [[ "$START" == "eq2.rst7" ]]; then
    OUT_NC=prod.nc; OUT_O=prod.out
else
    n=1; while [[ -e $(printf "prod_cont_%03d.nc" $n) ]]; do n=$((n+1)); done
    OUT_NC=$(printf "prod_cont_%03d.nc" $n); OUT_O=$(printf "prod_cont_%03d.out" $n)
    cp prod.rst7 prod_backup.rst7
fi
echo "  writing $OUT_NC"
pmemd.cuda -O -i prod.in -p "$TOP" -c "$START" \
           -o $OUT_O -r prod.rst7 -x $OUT_NC -inf prod.info
rc=$?
s=$(grep -m1 "Setting random seed to" $OUT_O 2>/dev/null | awk '{print $NF}')
[[ -n "$s" ]] && echo "$(date -Is) ${OUT_O%.out} ig=$s" >> seeds.txt
FIN=$(rst_time prod.rst7)
echo "FACT_RUN_DONE sys=$SYS rep=$REP rc=$rc time=${FIN:-unknown}/${TARGET_PS} ps"
echo "=== $(date) task $SLURM_ARRAY_TASK_ID finished ==="
exit $rc
