#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14   # amber22 has no kernels for k80/p100/h100/blackwell
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 12:00:00
#SBATCH -J us
#SBATCH -a 0-61%4
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/us_%a.log
#
# 112_us_run.sh -- umbrella sampling along the gate coordinate, holo vs apo.
#
# WHAT THIS IS FOR (README 59b)
# §24: unbiased MD never crosses the gate barrier -- neither state converts once in
# 1.8 us. So the open/closed FREE ENERGY difference, which is the quantity that
# separates a signalling switch from a constitutively closed protein (§58a), has
# never been measurable in this project. A PMF gives it.
#
# THIS RUN IS THE CALIBRATION, NOT THE RESULT. It must reproduce a known answer:
#     WT + ABA  -> CLOSED favoured
#     WT apo    -> OPEN favoured
# If the two PMFs do not separate in that direction the scheme is wrong, and no
# number it later produces for a designed pocket means anything. Same discipline
# as §13e, applied to a new method.
#
# COORDINATE: P88 CA (atom 1403) - R116 CA (atom 1819). Chosen by measuring every
# gate-latch CA pair against the two crystal references and taking the largest
# separation, then checked against 5 x 300 ns of unbiased MD:
#     closed basin  6.06-6.08 A       open basin  16.3-16.9 A     no overlap
# Atom indices are identical in the holo and apo topologies (protein comes first).
#
# WINDOWS: 5.0-20.0 A in 0.5 A steps, 31 per arm. With k = 10 kcal/mol/A^2 the
# thermal width is sqrt(kT/k) = 0.24 A, so neighbouring windows overlap at ~2 sigma.
#
# THE COORDINATE IS THE RISK, and it is stated up front. Gate closure is a loop
# rearrangement with orthogonal slow modes -- latch, ligand pose, side-chain
# repacking -- so a single distance will show hysteresis. Two mitigations are built
# in: windows are seeded from BOTH basins of real equilibrated trajectories rather
# than from one steered pull, so the low and high arms approach the barrier from
# opposite directions and disagreement in the overlap region is visible; and the
# quantity finally reported is the DIFFERENCE between holo and apo PMFs, in which
# systematic coordinate error largely cancels -- the same reason TI selectivity was
# usable when its absolute ddG was not (§43b, §46a).
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22_mpi_cuda >/dev/null 2>&1

MD=$ROOT/data/md191
US=$ROOT/data/umbrella
IAT1=1403        # P88 CA
IAT2=1819        # R116 CA
K=10.0           # kcal/mol/A^2
EQ_NS=${EQ_NS:-2}
PROD_NS=${PROD_NS:-20}

# 62 tasks = 2 arms x 31 windows
T=${SLURM_ARRAY_TASK_ID:-0}
if (( T < 31 )); then TAG=holo; TOP=$MD/S2_holo_closed/system.prmtop; I=$T
else                  TAG=apo;  TOP=$MD/S9_apo_closed/system.prmtop;  I=$((T-31)); fi
W=$(python3 -c "print(f'{5.0 + 0.5*$I:04.1f}')")
D=$US/$TAG/w$W
[[ -s $D/seed.rst7 ]] || { echo "no seed for $TAG w$W"; exit 1; }
cd "$D" || exit 1

# assert the restrained atoms really are P88 CA and R116 CA in THIS topology --
# the two arms use different prmtops and a silent index shift would bias the PMF
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/python - "$TOP" <<'PYX' || exit 1
import sys, parmed as pmd
t = pmd.load_file(sys.argv[1])
for idx, want_res, want_num in ((1403, "PRO", 88), (1819, "ARG", 116)):
    a = t.atoms[idx - 1]
    if a.name != "CA" or a.residue.name != want_res or a.residue.idx + 1 != want_num:
        sys.exit(f"atom {idx} is {a.residue.name}{a.residue.idx+1}@{a.name}, "
                 f"expected {want_res}{want_num}@CA")
print(f"    restrained atoms verified: PRO88 CA (1403), ARG116 CA (1819)")
PYX

cat > rst.dat <<RST
# umbrella restraint, window $W A
 &rst iat=$IAT1,$IAT2, r1=0.0, r2=$W, r3=$W, r4=99.0, rk2=$K, rk3=$K, /
RST

mk_in () {   # $1 name  $2 nstlim  $3 irest/ntx  $4 dumpfreq
cat > $1.in <<IN
umbrella $1, window $W A
 &cntrl
  imin=0, irest=$3, ntx=$( [[ $3 == 1 ]] && echo 5 || echo 1 ), nstlim=$2, dt=0.002,
  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  ntc=2, ntf=2, cut=10.0, ntt=3, gamma_ln=2.0, ig=-1, temp0=300.0,
  ntpr=5000, ntwx=0, ntwr=50000,
  nmropt=1,
 /
 &wt type='DUMPFREQ', istep1=$4 /
 &wt type='END' /
DISANG=rst.dat
DUMPAVE=$1.rc
IN
}

NEQ=$(python3 -c "print(int($EQ_NS*1000/0.002))")
NPR=$(python3 -c "print(int($PROD_NS*1000/0.002))")

echo "=== $TAG window $W A   k=$K   eq ${EQ_NS} ns, prod ${PROD_NS} ns   $(date -Is)"

run () {   # name prev
    local nm=$1 prev=$2
    [[ -s ${nm}.rst7 ]] && { echo "    $nm done"; return 0; }
    echo "    $nm ..."
    pmemd.cuda -O -i ${nm}.in -p "$TOP" -c "$prev" -o ${nm}.out -r ${nm}.rst7 \
        -inf ${nm}.mdinfo -x /dev/null 2>/dev/null \
        || { echo "    !! $nm FAILED"; tail -20 ${nm}.out 2>/dev/null; return 2; }
}

# equilibration under the restraint -- this is also what pulls the 4 windows whose
# seed sits up to 1.15 A off target onto their target value
mk_in eq  "$NEQ" 1 5000
run eq seed.rst7 || exit 2
mk_in prod "$NPR" 1 100          # dump the coordinate every 100 steps = 0.2 ps
run prod eq.rst7 || exit 2

n=$(grep -vc '^#' prod.rc 2>/dev/null || echo 0)
echo "    prod.rc: $n samples"
echo "=== $TAG w$W done $(date -Is)"
