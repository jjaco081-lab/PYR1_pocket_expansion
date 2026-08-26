#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 1:58:00
#SBATCH -J usreseed
#SBATCH -a 0-1
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/us_reseed_%a.log
#
# 127_us_reseed_open.sh -- repair the 38 dead umbrella windows.
#
# ⚠ THE BUG (README 70). The four md191 systems were solvated INDEPENDENTLY and so
# have different atom counts:
#     S2_holo_closed 46,570   S10_holo_open 51,588
#     S9_apo_closed  46,483   S1_apo_open   51,585
# Script 111 seeded the open-basin windows (11.0-20.0 A) from the OPEN systems,
# while 112 runs every window in an arm under the CLOSED system's topology. Every
# one of those 38 windows therefore died instantly with
#     ERROR: natom mismatch in inpcrd/restrt and prmtop files!
# and had never produced a single frame. 111 guarded holo-vs-apo ("a seed frame may
# only be taken from a trajectory built on the same topology") and missed
# closed-vs-open within an arm.
#
# Consequence: the completed windows span 5.0-10.5 A -- the CLOSED BASIN ONLY. A
# PMF needs both basins, so no open/closed free-energy difference was ever
# computable from what had run.
#
# THE FIX: adiabatic pulling, in the correct topology.
# Start from window 10.5 A, which is already complete, equilibrated, and in the
# right topology, and step the restraint centre outward in 0.5 A increments of
# 200 ps each, saving the final restart of each step as that window's seed. 19
# steps per arm = 3.8 ns, against ~10 ns for a Jarzynski pull, and it yields one
# seed per window directly rather than needing frames extracted afterwards.
#
# Slow stepping matters: the equilibration in 112 is only 2 ns, which can settle a
# structure already near its target but cannot drag one 10 A. That is what makes
# this a separate stage rather than something the restraint could absorb.
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

# restraint atoms COMPUTED per topology, identity asserted -- same rule as 112
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
echo "=== $TAG reseed: PRO88 CA = $IAT1, ARG116 CA = $IAT2, topology $TOP"

WORK=$US/$TAG/_reseed
mkdir -p "$WORK"; cd "$WORK" || exit 1

# start from the completed 10.5 A window: correct topology, already equilibrated
PREV=$US/$TAG/w10.5/prod04.rst7
[[ -s $PREV ]] || PREV=$(ls -t $US/$TAG/w10.5/prod*.rst7 2>/dev/null | head -1)
[[ -s $PREV ]] || { echo "no completed 10.5 A restart to start from"; exit 1; }
echo "    starting from $PREV"

NST=$(python3 -c "print(int($STEP_PS*1000/2))")   # 2 fs steps

for i in $(seq 0 18); do
    W=$(python3 -c "print(f'{11.0 + 0.5*$i:04.1f}')")
    D=$US/$TAG/w$W
    mkdir -p "$D"
    if [[ -s $D/seed.rst7 ]] && [[ $(sed -n '2p' "$D/seed.rst7" | awk '{print $1}') == \
          $(env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/python -c \
            "import parmed as pmd,sys;print(len(pmd.load_file('$TOP').atoms))") ]]; then
        echo "    w$W already has a topology-matched seed, skipping"
        PREV=$D/seed.rst7
        continue
    fi
    cat > rst_$W.dat <<RST
 &rst iat=$IAT1,$IAT2, r1=0.0, r2=$W, r3=$W, r4=99.0, rk2=$K, rk3=$K, /
RST
    cat > pull_$W.in <<IN
adiabatic pull to $W A
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
    echo "    pulling to $W A ($STEP_PS ps) from $(basename $PREV)"
    pmemd.cuda -O -i pull_$W.in -p "$TOP" -c "$PREV" -o pull_$W.out \
        -r pull_$W.rst7 -inf pull_$W.mdinfo -x /dev/null 2>/dev/null \
        || { echo "    !! pull to $W FAILED"; tail -20 pull_$W.out; exit 2; }
    cp pull_$W.rst7 "$D/seed.rst7"
    ACH=$(tail -2 pull_$W.rc | head -1 | awk '{print $2}')
    echo "        seeded w$W, coordinate reached ${ACH:-?} A"
    PREV=$D/seed.rst7
done
echo "=== $TAG reseed complete $(date -Is)"
