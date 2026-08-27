#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 1:58:00
#SBATCH -J usrev
#SBATCH -a 0-1
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/us_reverse_%a.log
#
# 136_us_reverse_seed.sh -- restore the two-directional seeding §81b showed was lost.
#
# THE PROBLEM (README 81b). §60's design deliberately seeded the low windows from
# CLOSED-basin trajectories and the high windows from OPEN-basin ones, "so the low
# and high arms approach the barrier from opposite directions and disagreement in
# the overlap region is visible". §70's repair replaced every open seed in BOTH arms
# with an adiabatic pull OUTWARD FROM CLOSED. The windows then ran -- but every open
# window descends from the closed state, the hysteresis control is gone, and a
# one-directional pulling bias would inflate the open free energy in both arms,
# which is exactly what the PMF shows (apo +6.92, ddG -1.90, both wrong).
#
# WHY NOT REBUILD THE TOPOLOGIES. The obvious fix -- re-solvate the S10/S1 open
# frames into the S2/S9 topologies -- founders on water counts (10,867 vs 12,120)
# and on a HID/HIE tautomer difference at H60, which sits next to K59 in the pocket.
# Copying coordinates across would need the waters re-placed and one histidine
# hydrogen rebuilt, and every window would have to be re-equilibrated.
#
# THE CHEAPER AND CLEANER FIX. Window 20.0 A has now run 20 ns of PRODUCTION under
# its own restraint. Whatever it descended from, it is a legitimate, equilibrated
# open-basin structure in the correct topology. So pull INWARD from it, 0.5 A per
# 200 ps, generating a second, independent seed for every window from the OPPOSITE
# direction. Running the windows from both seed sets is the hysteresis test §113
# promised and never implemented:
#
#   agreement  -> the PMF is converged and §81's failure is physics, not seeding
#   divergence -> the pulling direction determines the answer, and the calibration
#                 is measuring the protocol rather than the free energy
#
# Seeds are written to w<val>/seed_rev.rst7 so the outward-pulled seed.rst7 is
# preserved for the comparison.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22_mpi_cuda >/dev/null 2>&1
MD=$ROOT/data/md191
US=$ROOT/data/umbrella
K=10.0
STEP_PS=${STEP_PS:-200}

T=${SLURM_ARRAY_TASK_ID:-0}
if (( T == 0 )); then TAG=holo; TOP=$MD/S2_holo_closed/system.prmtop
else                  TAG=apo;  TOP=$MD/S9_apo_closed/system.prmtop; fi

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
echo "=== $TAG reverse seed: PRO88 CA=$IAT1  ARG116 CA=$IAT2"

WORK=$US/$TAG/_reverse
mkdir -p "$WORK"; cd "$WORK" || exit 1

# start from the fully-equilibrated 20.0 A window
PREV=$(ls -t $US/$TAG/w20.0/prod*.rst7 2>/dev/null | head -1)
[[ -s $PREV ]] || { echo "no completed 20.0 A restart"; exit 1; }
NS=$(env -u PYTHONPATH python3 -c "
import glob,os
n=sum(sum(1 for l in open(f) if l.strip() and not l.startswith('#'))
      for f in glob.glob('$US/$TAG/w20.0/prod*.rc'))
print(f'{n*0.2/1000:.1f}')")
echo "    origin: $PREV  (that window has ${NS} ns of production behind it)"

NST=$(python3 -c "print(int($STEP_PS*1000/2))")
for i in $(seq 0 18); do
    W=$(python3 -c "print(f'{19.5 - 0.5*$i:04.1f}')")
    D=$US/$TAG/w$W
    mkdir -p "$D"
    [[ -s $D/seed_rev.rst7 ]] && { echo "    w$W reverse seed exists, skip"; PREV=$D/seed_rev.rst7; continue; }
    cat > rst_$W.dat <<RST
 &rst iat=$IAT1,$IAT2, r1=0.0, r2=$W, r3=$W, r4=99.0, rk2=$K, rk3=$K, /
RST
    cat > pull_$W.in <<IN
adiabatic pull INWARD to $W A
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=$NST, dt=0.002,
  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  ntc=2, ntf=2, cut=10.0, ntt=3, gamma_ln=2.0, ig=-1, temp0=300.0,
  ntpr=5000, ntwx=0, ntwr=$NST,
  nmropt=1,
 /
 &wt type='DUMPFREQ', istep1=500 /
 &wt type='END' /
DISANG=rst_$W.dat
DUMPAVE=pull_$W.rc
IN
    echo "    pulling inward to $W A from $(basename $PREV)"
    pmemd.cuda -O -i pull_$W.in -p "$TOP" -c "$PREV" -o pull_$W.out \
        -r pull_$W.rst7 -inf pull_$W.mdinfo -x /dev/null 2>/dev/null \
        || { echo "    !! pull to $W FAILED"; tail -20 pull_$W.out; exit 2; }
    cp pull_$W.rst7 "$D/seed_rev.rst7"
    ACH=$(tail -2 pull_$W.rc | head -1 | awk '{print $2}')
    echo "        seeded w$W (reverse), coordinate reached ${ACH:-?} A"
    PREV=$D/seed_rev.rst7
done
echo "=== $TAG reverse seeding complete $(date -Is)"
