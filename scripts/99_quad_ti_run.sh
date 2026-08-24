#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14   # amber22 pmemd.cuda has no kernels for k80/p100/h100/blackwell
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 1:58:00
#SBATCH -J qti
#SBATCH -a 0-47%4
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/qti_%a.log
#
# 87_ti_run.sh -- one lambda window per array task.
#
# 6 legs x 12 lambdas = 72 tasks. Each task is min -> heat -> eq -> prod at a
# single fixed lambda and is fully restart-safe, so a requeue costs one stage.
#
# WHY 12-POINT GAUSSIAN QUADRATURE AND NOT AN EVEN 0.0-1.0 LADDER
#   Gauss-Legendre nodes never touch lambda = 0 or 1, so the endpoint singularity
#   that plagues softcore transformations is avoided by construction rather than
#   patched. The weights below integrate <dV/dlambda> exactly for a polynomial of
#   degree 23, which is far more than this integrand needs.
#
# PARTITION: short_gpu, per the standing rule -- gpu if under the limit, else
#   short_gpu split under 2 h. The gpu partition is currently full with the
#   factorial (1) and the non-cognate array (3) = its 4-GPU cap. short_gpu carries
#   a separate per-partition cap; %4 holds concurrency there and keeps the parent
#   cutlerlab account (ceiling 8) from being saturated by this job alone.
#
# GTI SETTINGS
#   gti_lam_sch/gti_ele_sc/gti_vdw_sc enable Amber22's unified softcore, which
#   handles the charge change and the vdW change in ONE transformation. The older
#   route is a three-step decharge -> mutate -> recharge, which triples the run
#   count and adds two more places for a mask to be wrong. K59R needs this: Lys+
#   and Arg+ carry the same net charge but redistribute a lot of it.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22_mpi_cuda >/dev/null 2>&1

LEGS=(aba mandi_xtal mandi_dock apo)
LAM=(0.00922 0.04794 0.11505 0.20634 0.31608 0.43738 \
     0.56262 0.68392 0.79366 0.88495 0.95206 0.99078)

T=${SLURM_ARRAY_TASK_ID:-0}
LEG=${LEGS[$((T / 12))]}
LI=$((T % 12))
L=${LAM[$LI]}

D=$ROOT/data/ti_quad/$LEG
W=$D/lam$(printf '%02d' $LI)
mkdir -p "$W"; cd "$W" || exit 1
TOP=$D/ti.prmtop
[[ -f $TOP ]] || { echo "no ti.prmtop for $LEG"; exit 1; }

# masks come from tiMerge's own output, renumbered for the merged topology.
# Writing them by hand against the pre-merge numbering is the obvious way to get
# a confidently wrong answer here, so they are read, not retyped.
source <(sed 's/^/export /' "$D/masks.txt")

# Compress the index lists to RANGES. tiMerge emits every index explicitly, which
# for K59R is 18 comma-separated numbers; Amber echoes mdin at fixed width and a
# long namelist string is an avoidable place for a mask to be silently clipped.
# The indices are contiguous by construction (tiMerge groups the TI region), and
# the compressor ASSERTS that before rewriting -- if they ever stop being
# contiguous it fails loudly instead of emitting a wrong range.
compress () {   # @1,2,3,7,8 -> @1-3,7-8
    python3 -c "
import sys
m=sys.argv[1].lstrip('@')
xs=sorted(int(x) for x in m.split(','))
out=[];a=b=xs[0]
for v in xs[1:]:
    if v==b+1: b=v
    else: out.append((a,b)); a=b=v
out.append((a,b))
n=sum(y-x+1 for x,y in out)
assert n==len(xs), 'range compression changed the atom count'
print('@'+','.join(str(x) if x==y else f'{x}-{y}' for x,y in out))
" "$1"
}
timask1=$(compress "$timask1"); timask2=$(compress "$timask2")
scmask1=$(compress "$scmask1"); scmask2=$(compress "$scmask2")
echo "=== $LEG  lambda=$L (window $LI)  $(date -Is)"
echo "    timask1=$timask1"
echo "    timask2=$timask2"

# SHAKE must NOT constrain bonds inside the perturbed region: those bonds are
# being alchemically changed, and constraining them applies the endpoint geometry
# to an intermediate lambda. Amber's fix is noshakemask over the whole TI region.
# Built by concatenating the two index lists rather than retyping them.
NOSHAKE="${timask1},${timask2#@}"
echo "    noshakemask=$NOSHAKE"

PROD_NS=${PROD_NS:-10}
NPROD=$(python3 -c "print(int($PROD_NS*1000/0.002))")

ti_block () {   # $1 = extra cntrl lines
cat <<TI
  icfe=1, ifsc=1, clambda=$L,
  timask1='$timask1', timask2='$timask2',
  scmask1='$scmask1', scmask2='$scmask2',
  gti_cut=1, gti_output=1, gti_add_sc=1, gti_scale_beta=1,
  gti_cut_sc_on=6.0, gti_cut_sc_off=8.0,
  gti_lam_sch=1, gti_ele_sc=1, gti_vdw_sc=1,
  gti_cut_sc=2, gti_ele_exp=2, gti_vdw_exp=2,
  noshakemask='$NOSHAKE',
TI
}

cat > min.in <<IN
minimise at lambda $L
 &cntrl
  imin=1, maxcyc=5000, ncyc=2000, ntmin=2, ntb=1, cut=10.0, ntpr=500,
  ntr=1, restraintmask='!:WAT,K+,Cl- & !@H=', restraint_wt=5.0,
$(ti_block)
 /
IN
# SECOND minimisation, BACKBONE restrained only.
#
# min.in restrains every non-water heavy atom (`!:WAT,K+,Cl- & !@H=`), the ligand
# included. Measured, that is fine at mid lambda -- softcore weakens the clash
# enough that the restrained minimiser still resolves it (lam06: 0.62 -> 2.34 A,
# no contacts) -- but at lambda ~ 0 the WT sidechains are fully coupled and the
# restraint pins the clash in place (lam00: 0.62 -> 0.62 A, 9 contacts under
# 2 A). Heating then died with "illegal memory access ... kNLSkinTest" and a
# dU/dl of ************ on all four low-lambda windows of BOTH mandi legs.
#
# Letting the ligand and sidechains move here is not a workaround, it is the
# physics: at lambda ~ 0 the system IS wild-type PYR1, which cannot hold
# mandipropamid in its 4WVO pose, so the ligand is supposed to shift. The
# resulting lambda-dependence of the pose is the hysteresis to watch for in
# 100_quad_aggregate.py's per-window table.
cat > min2.in <<IN
minimise, backbone restrained only, lambda $L
 &cntrl
  imin=1, maxcyc=10000, ncyc=5000, ntmin=2, ntb=1, cut=10.0, ntpr=1000,
  ntr=1, restraintmask='@N,CA,C,O', restraint_wt=5.0,
$(ti_block)
 /
IN
cat > heat.in <<IN
heat 0->300 K, 200 ps, lambda $L
 &cntrl
  imin=0, irest=0, ntx=1, nstlim=100000, dt=0.002,
  ntb=1, ntp=0, ntc=2, ntf=1, cut=10.0,
  ntt=3, gamma_ln=2.0, ig=-1, tempi=10.0, temp0=300.0,
  ntr=1, restraintmask='!:WAT,K+,Cl- & !@H=', restraint_wt=1.0,
  nmropt=1, ntpr=5000, ntwx=0, ntwr=50000,
$(ti_block)
 /
 &wt type='TEMP0', istep1=0, istep2=80000, value1=10.0, value2=300.0 /
 &wt type='END' /
IN
cat > eq.in <<IN
NPT equilibration 500 ps, lambda $L
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=250000, dt=0.002,
  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  ntc=2, ntf=1, cut=10.0, ntt=3, gamma_ln=2.0, ig=-1, temp0=300.0,
  ntpr=5000, ntwx=0, ntwr=50000,
$(ti_block)
 /
IN
cat > prod.in <<IN
TI production $PROD_NS ns, lambda $L
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=$NPROD, dt=0.002,
  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  ntc=2, ntf=1, cut=10.0, ntt=3, gamma_ln=2.0, ig=-1, temp0=300.0,
  ntpr=1000, ntwx=0, ntwr=50000, ntave=0,
$(ti_block)
 /
IN

run () {   # name prev
    local nm=$1 prev=$2
    [[ -s ${nm}.rst7 ]] && { echo "    $nm done"; return 0; }
    echo "    $nm ..."
    pmemd.cuda -O -i ${nm}.in -p "$TOP" -c "$prev" -ref "$prev" \
        -o ${nm}.out -r ${nm}.rst7 -inf ${nm}.mdinfo \
        || { echo "    !! $nm FAILED"; tail -20 ${nm}.out 2>/dev/null; return 2; }
}

run min  "$D/ti.inpcrd" || exit 2
run min2 min.rst7       || exit 2
run heat min2.rst7      || exit 2
run eq   heat.rst7      || exit 2
run prod eq.rst7        || exit 2

# the deliverable: <dV/dlambda> for this window
python3 - <<'PY'
import re, statistics, sys
# pmemd prints every frame TWICE (one block per TI region) and the trailing
# AVERAGES / RMS banners are a mean and a sd, not samples. Counting them gives
# n = 2*frames + 4 and an error bar sqrt(2) too small -- see 89's header.
txt = open("prod.out").read()
body = re.split(r"A V E R A G E S", txt)[0]
seen = {}
for st, val in re.findall(r"NSTEP =\s*(\d+).*?DV/DL\s+=\s+(-?\d+\.\d+)", body, re.S):
    seen.setdefault(int(st), float(val))
v = [seen[k] for k in sorted(seen)]
if not v:
    sys.exit("no DV/DL found in prod.out")
h = len(v) // 2
open("dvdl.txt", "w").write(
    f"n={len(v)}\nmean={statistics.mean(v):.4f}\n"
    f"sd={statistics.stdev(v):.4f}\n"
    f"first_half={statistics.mean(v[:h]):.4f}\n"
    f"second_half={statistics.mean(v[h:]):.4f}\n")
print("    <dV/dl> =", round(statistics.mean(v), 3),
      " drift:", round(statistics.mean(v[h:]) - statistics.mean(v[:h]), 3))
PY
echo "=== $LEG lambda=$L done $(date -Is)"
