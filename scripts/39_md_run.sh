#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 7-00:00:00
#SBATCH -J pe_md
#SBATCH -a 0-11%4
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/md_%a.log
#
# 39_md_run.sh -- WT PYR1 MD baseline: 4 systems x 3 replicates.
#
# JOB SHAPE
#   The account is limited to 4 concurrent GPUs, so this is a 12-task array
#   throttled with %4: one task per (system, replicate), never more than four
#   running. One replicate per task keeps every job inside the 7-day limit --
#   packing 3 replicates into one job would put S4 (111,680 atoms) at ~10 days
#   and it would be killed.
#
#     task = system_index * 3 + replicate      systems below, replicates 0-2
#
#   --gres=gpu:a100:1 is deliberate. The gpu partition also holds p100 and k80
#   (k80 currently down); A100 is the deep pool (24 cards) and several times
#   faster. The newest cards (ada6000, blackwell6000) exist ONLY under
#   preempt_gpu, which is preemptible -- not used here by request.
#
# PROTOCOL (Amber 22 pmemd.cuda)
#   min1  5000 steps, solute heavy atoms restrained 10 kcal/mol/A^2
#   min2  5000 steps, unrestrained
#   heat  NVT 0->300 K over 200 ps, restraints 5.0
#   eq1   NPT 300 K, 500 ps, restraints 1.0
#   eq2   NPT 300 K, 1 ns, unrestrained
#   prod  NPT 300 K, 300 ns, 2 fs, SHAKE on H, PME, 10 A cutoff,
#         Langevin gamma_ln=2.0, Monte Carlo barostat, frames every 10 ps
#
#   2 fs without hydrogen-mass repartitioning is the conservative choice. HMR at
#   4 fs would roughly double throughput; if this campaign needs to be repeated
#   for many candidates, that is the first optimisation to make.
#
# RESTART
#   Production auto-resumes from prod.rst7 if present, so a job that hits the
#   walltime can simply be resubmitted. Completed stages are skipped.
set -uo pipefail
module load amber/22_mpi_cuda >/dev/null 2>&1

P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
MD=$P/data/md
SYSTEMS=(S1_apo_open S2_holo_closed S3_apo_dimer S4_ternary)
NREP=3
PROD_NS=300

i=$SLURM_ARRAY_TASK_ID
sys=${SYSTEMS[$((i / NREP))]}
rep=$((i % NREP))
D=$MD/$sys/rep$rep
mkdir -p "$D"
cd "$D" || exit 1

TOP=$MD/$sys/system.prmtop
CRD=$MD/$sys/system.inpcrd
[[ -s "$TOP" ]] || { echo "missing $TOP -- run 38b_md_build.py"; exit 1; }

echo "=== task $i : $sys replicate $rep ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
date

# ---------------- input files ----------------
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
NSTEPS=$((PROD_NS * 500000))
cat > prod.in <<EOF
NPT production, ${PROD_NS} ns, frames every 10 ps
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

run_stage () {   # name  input  prev_rst  [-ref]
    local nm=$1 inp=$2 prev=$3
    if [[ -s ${nm}.rst7 ]]; then echo "  $nm already done, skipping"; return 0; fi
    echo "  running $nm ..."
    local ref=()
    [[ "${4:-}" == "ref" ]] && ref=(-ref "$prev")
    pmemd.cuda -O -i ${nm}.in -p "$TOP" -c "$prev" \
               -o ${nm}.out -r ${nm}.rst7 -x ${nm}.nc -inf ${nm}.info "${ref[@]}" \
        || { echo "  !! $nm FAILED (see ${nm}.out)"; exit 2; }
}

run_stage min1 min1.in "$CRD" ref
run_stage min2 min2.in min1.rst7
run_stage heat heat.in min2.rst7 ref
run_stage eq1  eq1.in  heat.rst7 ref
run_stage eq2  eq2.in  eq1.rst7

# production, resumable
if [[ -s prod.rst7 ]]; then
    echo "  resuming production from prod.rst7"
    cp prod.rst7 prod_prev.rst7
    pmemd.cuda -O -i prod.in -p "$TOP" -c prod_prev.rst7 \
        -o prod_cont.out -r prod.rst7 -x prod_cont.nc -inf prod.info
else
    pmemd.cuda -O -i prod.in -p "$TOP" -c eq2.rst7 \
        -o prod.out -r prod.rst7 -x prod.nc -inf prod.info
fi

echo "MD_TASK_DONE task=$i sys=$sys rep=$rep rc=$?"
date
ls -la *.nc 2>/dev/null | awk '{print "  "$9" "$5" bytes"}'
