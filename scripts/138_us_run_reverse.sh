#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14   # amber22 has no kernels for k80/p100/h100/blackwell
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 1:58:00
#SBATCH -J usrev2
#SBATCH -a 0-61%4        # short_gpu carries its own 4-GPU allowance, separate from `gpu`
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/us_rev2_%a.log
#
# 138_us_run_reverse.sh -- rerun the umbrella windows from the REVERSE seeds.
#
# THE HYSTERESIS TEST (README 81b, 60's original design). §70's repair pulled every
# open window OUTWARD from the closed basin, destroying the two-directional seeding
# §60 built in so that "disagreement in the overlap region is visible". §136 restored
# the missing direction by pulling INWARD from the now-equilibrated 20.0 A window --
# a legitimate open-basin structure in the correct topology.
#
# Running the same windows from both seed sets is the comparison §113's docstring
# promised and the code never implemented:
#
#   PMFs agree     -> converged; §81's failure is physics, not seeding
#   PMFs disagree  -> the pull direction sets the answer, and the calibration was
#                     measuring the protocol
#
# Output goes to w<val>/reverse/ so the forward run is untouched and both are
# available to 113 side by side. Only 10.5-19.5 A have reverse seeds (§136 pulled
# inward from 20.0), so tasks outside that range exit 0 immediately.
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
# COORDINATE: P88 CA - R116 CA. Chosen by measuring every gate-latch CA pair
# against the two crystal references and taking the largest separation, then
# checked against 5 x 300 ns of unbiased MD:
#     closed basin  6.06-6.08 A       open basin  16.3-16.9 A     no overlap
# The atom INDICES are computed per topology, never hardcoded: they are 1403/1819
# in the WT arms, but K59R adds atoms before residue 88 so the PYR1^MANDI arms
# shift, and a hardcoded pair would silently restrain the wrong atoms there.
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
K=10.0           # kcal/mol/A^2
EQ_NS=${EQ_NS:-2}
PROD_NS=${PROD_NS:-20}

# 62 tasks = 2 arms x 31 windows
T=${SLURM_ARRAY_TASK_ID:-0}
if (( T < 31 )); then TAG=holo; TOP=$MD/S2_holo_closed/system.prmtop; I=$T
else                  TAG=apo;  TOP=$MD/S9_apo_closed/system.prmtop;  I=$((T-31)); fi
W=$(python3 -c "print(f'{5.0 + 0.5*$I:04.1f}')")
D=$US/$TAG/w$W/reverse
mkdir -p "$D"; [[ -s $US/$TAG/w$W/seed_rev.rst7 ]] || { echo "no REVERSE seed for $TAG w$W (only 10.5-19.5 A have one)"; exit 0; }; cp -n $US/$TAG/w$W/seed_rev.rst7 $D/seed.rst7
cd "$D" || exit 1

# P88 CA and R116 CA, COMPUTED per topology. They are 1403/1819 in the WT arms,
# but K59R adds atoms before residue 88, so the quad arms shift -- hardcoding
# would silently restrain the wrong pair in half the systems.
read IAT1 IAT2 <<< "$(env -u PYTHONPATH \
  /bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/python - "$TOP" <<'PYX'
import sys, parmed as pmd
t = pmd.load_file(sys.argv[1])
out = []
for num, want in ((88, "PRO"), (116, "ARG")):
    r = t.residues[num - 1]
    if r.name != want:
        sys.exit(f"residue {num} is {r.name}, expected {want}")
    ca = [a for a in r.atoms if a.name == "CA"]
    if len(ca) != 1:
        sys.exit(f"residue {num} has {len(ca)} CA atoms")
    out.append(ca[0].idx + 1)
print(*out)
PYX
)"
[[ -n "${IAT1:-}" && -n "${IAT2:-}" ]] || { echo "could not resolve restraint atoms"; exit 1; }
echo "    restraint atoms: PRO88 CA = $IAT1, ARG116 CA = $IAT2"


cat > rst.dat <<RST
# umbrella restraint, window $W A
 &rst iat=$IAT1,$IAT2, r1=0.0, r2=$W, r3=$W, r4=99.0, rk2=$K, rk3=$K, /
RST

mk_in () {   # $1 name  $2 nstlim  $3 irest/ntx  $4 dumpfreq
cat > $1.in <<IN
umbrella $1, window $W A
 &cntrl
  imin=0, irest=$3, ntx=$( [[ $3 == 1 ]] && echo 5 || echo 1 ), nstlim=$2, dt=0.002,
$( [[ $3 == 0 ]] && echo "  tempi=300.0," )
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

# Equilibration under the restraint. irest=0 / ntx=1, NOT irest=1: the seeds were
# extracted from trajectory frames by cpptraj, and a trajectory carries no
# velocities, so pmemd rejects them with "could not find enough velocities in
# seed.rst7". Velocities are reassigned from a 300 K Maxwell distribution instead,
# which is right here anyway -- the frame came from a 300 K trajectory and this
# stage exists to re-equilibrate it under the bias. It is also what pulls the four
# windows whose seed sits up to 1.15 A off target onto their target value.
mk_in eq "$NEQ" 0 5000
run eq seed.rst7 || exit 2

# ---------------------------------------------------------------------------
# DEADLINE-AWARE CHUNKING, so a 20 ns window fits short_gpu's 2 h limit.
#
# 20 ns at the measured 225 ns/day is 2.14 h, which does not fit; 2 ns of
# equilibration plus 20 ns of production certainly does not. So production runs
# in CHUNK_NS pieces and the job stops cleanly while it still has time to write a
# restart. Re-submitting the array continues from wherever each window got to.
#
# This is continuation, not restart-from-scratch: irest=1/ntx=5 carries velocities
# and box across, so the trajectory is unbroken. `ig=-1` does draw a fresh Langevin
# seed each chunk, which is harmless for equilibrium sampling -- the samples remain
# a valid draw from the canonical ensemble, and the seam only matters if one were
# computing a time correlation function across it, which WHAM does not.
#
# Accumulated production is MEASURED from the .rc files each pass, never counted
# from a variable, so a preempted or requeued chunk cannot inflate the total.
# ---------------------------------------------------------------------------
CHUNK_NS=${CHUNK_NS:-5}
BUDGET_S=${BUDGET_S:-6300}       # 105 min of the 118 min limit
START=$SECONDS

have_ns () {
    python3 - "$PWD" <<'PYX'
import glob, os, sys
wd = sys.argv[1]
n = 0
for f in sorted(glob.glob(os.path.join(wd, "prod*.rc"))):
    n += sum(1 for l in open(f) if l.strip() and not l.startswith("#"))
print(f"{n * 0.2 / 1000:.2f}")          # one sample per 100 steps = 0.2 ps
PYX
}

NCH=$(python3 -c "print(int($CHUNK_NS*1000/0.002))")
i=1
while :; do
    CUR=$(have_ns)
    if python3 -c "import sys; sys.exit(0 if $CUR >= $PROD_NS-0.2 else 1)"; then
        echo "    have ${CUR} ns >= ${PROD_NS} ns target -- window complete"; break
    fi
    # stop while there is still time to finish a chunk and write its restart
    LEFT=$((BUDGET_S - (SECONDS - START)))
    NEED=$(python3 -c "print(int($CHUNK_NS/225.0*86400*1.25))")
    if (( LEFT < NEED )); then
        echo "    ${LEFT}s left, a ${CHUNK_NS} ns chunk needs ~${NEED}s -- stopping cleanly at ${CUR} ns"
        echo "    resubmit the array to continue"; break
    fi
    SEG=prod$(printf '%02d' $i)
    if [[ -s ${SEG}.rst7 ]]; then i=$((i+1)); continue; fi
    PREV=eq.rst7
    j=$((i-1)); [[ $j -ge 1 ]] && PREV=prod$(printf '%02d' $j).rst7
    echo "    chunk $i (have ${CUR} ns) from $PREV"
    mk_in "$SEG" "$NCH" 1 100
    run "$SEG" "$PREV" || exit 2
    i=$((i+1))
done
echo "    total production: $(have_ns) ns"
echo "=== $TAG w$W $(date -Is)"
