#!/bin/bash
#SBATCH -p preempt_gpu
#SBATCH -A preempt
#SBATCH --gres=gpu:ada6000:1
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 7-00:00:00
#SBATCH -J pe_md_804
#SBATCH -a 0
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/md_804_%a.log
#
# 57_md_run_804_1.sh -- S5_804_1_ternary (PYR1 + ABA + the 804_1 designed binder).
#
# WHAT S5 IS FOR
#   S4 is the WT signalling complex, PYR1 + ABA + HAB1 + Mn2+, and it now has one
#   complete 300 ns replicate. S5 is that complex with HAB1 replaced by 804_1, the
#   RFdiffusion binder that was actually sequenced (4/4 colonies, 2026-08-13). Run
#   at the SAME 300 ns under the SAME protocol precisely so the two are comparable
#   -- a candidate trajectory means nothing without the WT reference under
#   identical conditions.
#
#   Built from an AF3 model, not a crystal: 79,290 atoms vs S4's 111,680, because
#   804_1 (143 aa) is far smaller than HAB1. NO Mn2+ -- the metal in S4 is HAB1's
#   catalytic phosphatase centre, and 804_1 is not a phosphatase.
#
#   NUMBERING: chain A of the AF3 model starts GAMASEL... where the real Y2H
#   construct starts MPSEL..., so residue indices here are OFFSET BY 2 from native
#   PYR1 numbering (native K59 = K61 here, F108 = F110, latch H115 = H117).
#   Convert before comparing anything against PYR1_pocket_expansion numbering.
#
#   ONE replicate (-a 0). The preempt account's QOS allows gres/gpu=1, so extra
#   array tasks serialise rather than run in parallel; add rep1/rep2 as separate
#   submissions once rep0 shows the system is stable.
#
# THIS FILE IS A COPY of 42_md_run_preempt.sh with the system name changed. 42 is
# left untouched as the exact record of what produced S4.
#
# WHY A SEPARATE SCRIPT
#   preempt_gpu has ada6000/blackwell cards that `gpu` cannot see. The cost is
#   preemption: PreemptMode=REQUEUE, so a preempted job is automatically put back
#   in the queue and this script re-runs from the top.
#
#   CORRECTION 2026-08-14: this header used to say all three A100 nodes in `gpu`
#   had gone to DRAINING. They had not. `sinfo` prints `mixed-`, which is
#   MIXED+PLANNED -- fully allocated with a backfill reservation -- and the
#   trailing `-` was misread as a drain flag. Read node state with
#   `scontrol show node <n>`, which spells it out, not the sinfo suffix. `gpu` is
#   healthy and returns real start estimates, so it stays a viable home for this
#   job; preempt_gpu was chosen for latency, not because `gpu` was broken.
#
#   The real constraint runs the other way: the `preempt` account is capped at
#   gres/gpu=1 on preempt_gpu, so this job queues behind ANY other preempt_gpu job
#   of the same user -- pending reason AssocGrpGRES, not Priority. At submission it
#   sat behind two esmf2_pyr1fix jobs and had still not started a day later.
#
#   That auto-requeue is exactly why 39_md_run.sh could NOT simply be pointed at
#   preempt_gpu. It has two bugs that are harmless on a non-preemptible
#   partition and destructive on a preemptible one:
#
#     1. LOST SEGMENTS. Its resume path always writes -x prod_cont.nc with -O.
#        A second preemption overwrites the trajectory the first one produced.
#        Here every resume writes prod_cont_NNN.nc with a fresh index, so
#        segments accumulate instead of clobbering each other.
#
#     2. RUNAWAY LENGTH. Its prod.in has a fixed nstlim of 300 ns, and restarts
#        use irest=1, which continues the clock. Resuming therefore ran a FURTHER
#        300 ns rather than the remainder. Here the remaining step count is
#        computed from the restart file's own clock each time.
#
#   The clock includes equilibration: heat 200 ps + eq1 500 ps + eq2 1 ns =
#   1700 ps, so a finished replicate reads 301700 ps, not 300000.
#
# CRASH SAFETY
#   prod.rst7 is rewritten every 1 ns (ntwr=500000). A preemption during that
#   write could truncate it, so the previous restart is kept as prod_backup.rst7
#   and validated with ncdump before use; if prod.rst7 is unreadable the backup
#   is used instead. Worst case loss is ~1 ns.
#
# NO COORDINATION NEEDED
#   Unlike S4, nothing else writes data/md/S5_804_1_ternary/rep0. Job 27333712's
#   task _9 targeted S4_ternary/rep0 and was unrelated to this system; it was held
#   for its own reasons and then CANCELLED on 2026-08-14, once S4 rep0 had finished
#   via 42_md_run_preempt.sh, so that nobody could release it into a redundant
#   300 ns on top of a complete trajectory.
#   The general rule stands: never let two pmemd processes write one prod.rst7.
set -uo pipefail
module load amber/22_mpi_cuda >/dev/null 2>&1

P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
MD=$P/data/md
SYS=S5_804_1_ternary
PROD_NS=300
EQUIL_PS=1700                       # heat 200 + eq1 500 + eq2 1000
TARGET_PS=$(( EQUIL_PS + PROD_NS * 1000 ))

rep=$SLURM_ARRAY_TASK_ID
D=$MD/$SYS/rep$rep
mkdir -p "$D"
cd "$D" || exit 1

TOP=$MD/$SYS/system.prmtop
CRD=$MD/$SYS/system.inpcrd
[[ -s "$TOP" ]] || { echo "missing $TOP"; exit 1; }

echo "=== $(date)  $SYS rep$rep  job $SLURM_JOB_ID  restarts=${SLURM_RESTART_COUNT:-0} ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ---------------- equilibration inputs (identical to 39_md_run.sh) ----------------
cat > min1.in <<EOF
minimisation 1, solute restrained
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=1, cut=10.0,
  ntr=1, restraintmask='!:WAT,Na+,Cl- & !@H=', restraint_wt=10.0,
 /
EOF
cat > min2.in <<EOF
minimisation 2, unrestrained
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=1, cut=10.0, ntr=0,
 /
EOF
cat > heat.in <<EOF
NVT heating 0 -> 300 K over 200 ps
 &cntrl
  imin=0, irest=0, ntx=1,
  nstlim=100000, dt=0.002,
  ntc=2, ntf=2, cut=10.0,
  ntb=1, ntp=0,
  ntt=3, gamma_ln=2.0, ig=-1,
  tempi=0.0, temp0=300.0,
  nmropt=1,
  ntr=1, restraintmask='!:WAT,Na+,Cl- & !@H=', restraint_wt=5.0,
  ntpr=5000, ntwx=5000, ntwr=50000,
 /
 &wt type='TEMP0', istep1=0, istep2=80000, value1=0.0, value2=300.0 /
 &wt type='TEMP0', istep1=80001, istep2=100000, value1=300.0, value2=300.0 /
 &wt type='END' /
EOF
cat > eq1.in <<EOF
NPT equilibration 1, 500 ps, light restraints
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=250000, dt=0.002,
  ntc=2, ntf=2, cut=10.0,
  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  ntt=3, gamma_ln=2.0, ig=-1, temp0=300.0,
  ntr=1, restraintmask='!:WAT,Na+,Cl- & !@H=', restraint_wt=1.0,
  ntpr=5000, ntwx=5000, ntwr=50000,
 /
EOF
cat > eq2.in <<EOF
NPT equilibration 2, 1 ns, unrestrained
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=500000, dt=0.002,
  ntc=2, ntf=2, cut=10.0,
  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  ntt=3, gamma_ln=2.0, ig=-1, temp0=300.0,
  ntr=0,
  ntpr=5000, ntwx=5000, ntwr=50000,
 /
EOF

run_stage () {   # name  prev_rst  [ref]
    local nm=$1 prev=$2
    if [[ -s ${nm}.rst7 ]]; then echo "  $nm done, skipping"; return 0; fi
    echo "  running $nm ..."
    local ref=()
    [[ "${3:-}" == "ref" ]] && ref=(-ref "$prev")
    pmemd.cuda -O -i ${nm}.in -p "$TOP" -c "$prev" \
               -o ${nm}.out -r ${nm}.rst7 -x ${nm}.nc -inf ${nm}.info "${ref[@]}" \
        || { echo "  !! $nm FAILED"; exit 2; }
}

run_stage min1 "$CRD" ref
run_stage min2 min1.rst7
run_stage heat min2.rst7 ref
run_stage eq1  heat.rst7 ref
run_stage eq2  eq1.rst7

# ---------------- production, restart-aware ----------------
rst_time () {    # echo simulated ps in $1, or nothing if unreadable
    ncdump -v time "$1" 2>/dev/null | awk '/^ time =/{gsub(/[^0-9.]/,"",$3); print $3; exit}'
}

START=eq2.rst7
NOW=$EQUIL_PS
if [[ -s prod.rst7 ]]; then
    T=$(rst_time prod.rst7)
    if [[ -z "$T" ]]; then
        echo "  !! prod.rst7 unreadable (truncated by preemption?)"
        if [[ -s prod_backup.rst7 ]]; then
            T=$(rst_time prod_backup.rst7)
            [[ -n "$T" ]] && { cp prod_backup.rst7 prod.rst7; echo "  recovered from prod_backup.rst7 at $T ps"; }
        fi
    fi
    if [[ -n "$T" ]]; then START=prod_resume.rst7; NOW=$T; cp prod.rst7 "$START"; fi
fi

REMAIN_PS=$(awk -v a=$TARGET_PS -v b=$NOW 'BEGIN{printf "%.0f", a-b}')
if [[ "$REMAIN_PS" -le 0 ]]; then
    echo "  production already complete ($NOW / $TARGET_PS ps)"
    echo "MD_TASK_DONE sys=$SYS rep=$rep rc=0"
    exit 0
fi
NSTEPS=$(( REMAIN_PS * 500 ))
echo "  production: at $NOW ps, target $TARGET_PS ps -> $REMAIN_PS ps ($NSTEPS steps)"

cat > prod.in <<EOF
NPT production, remaining ${REMAIN_PS} ps, frames every 10 ps
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=${NSTEPS}, dt=0.002,
  ntc=2, ntf=2, cut=10.0,
  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  ntt=3, gamma_ln=2.0, ig=-1, temp0=300.0,
  ntr=0,
  ntpr=25000, ntwx=5000, ntwr=500000,
  iwrap=1,
 /
EOF

# next free segment index so repeated preemption never overwrites a segment
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

FIN=$(rst_time prod.rst7)
echo "MD_TASK_DONE sys=$SYS rep=$rep rc=$rc time=${FIN:-unknown}/${TARGET_PS} ps"
date
