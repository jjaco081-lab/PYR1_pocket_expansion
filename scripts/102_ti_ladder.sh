#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14   # amber22 pmemd.cuda has no kernels for k80/p100/h100/blackwell
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 12:00:00
#SBATCH -J qtild
#SBATCH -a 0-1
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/qtild_%a.log
#
# 102_ti_ladder.sh -- rescue the low-lambda mandipropamid windows with a
# DESCENDING LAMBDA LADDER.
#
# WHY THE PARALLEL WINDOWS FAILED
# In 99 every window is built independently from ti.inpcrd, which holds the
# CRYSTAL (or docked) pose. That is right at high lambda, where the quad
# sidechains are coupled and the pose fits. It is impossible at low lambda:
# there the WT sidechains are fully coupled and F108 sits 0.62 A from a ligand
# heavy atom. Measured, minimisation cannot escape it --
#
#     lam06  0.62 -> 2.34 A, 0 contacts        (softcore still weakens the clash)
#     lam00  0.62 -> 0.62 A, 9 contacts        (WT fully coupled; stuck)
#     min2 at lam00: E = 4.5e8 kcal/mol, |F|max = 8.4e6, FLAT over 10,000 steps
#
# and heating then died with "illegal memory access ... kNLSkinTest" and a
# dU/dl of ************ on all four low-lambda windows of BOTH mandi arms.
# A second minimisation stage does not fix this and was the wrong diagnosis:
# there is no downhill path out of a 0.62 A contact, so no minimiser finds one.
#
# WHAT THIS DOES INSTEAD
# Walk DOWN in lambda, seeding each window from the equilibrated structure of the
# window above it:
#
#     lam04 (done) -> lam03 -> lam02 -> lam01 -> lam00
#
# Each step is a small change in lambda, so the ligand leaves the clash
# gradually, re-equilibrating at every rung. This is the standard treatment for a
# perturbation with strong ligand-pose coupling, and it is not a trick to make
# the number come out: at lambda ~ 0 the system IS wild-type PYR1, which cannot
# hold mandipropamid in its 4WVO pose, so the ligand is SUPPOSED to move. The
# resulting pose-vs-lambda curve is the hysteresis to read in
# 100_quad_aggregate.py's per-window table.
#
# The cost is that these four windows are sequential rather than parallel, which
# is why this runs on `gpu` (7-day limit) rather than short_gpu (2 h). The gpu
# partition is free now that the non-cognate array has finished.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22_mpi_cuda >/dev/null 2>&1

LEGS=(mandi_xtal mandi_dock)
LAM=(0.00922 0.04794 0.11505 0.20634 0.31608 0.43738 \
     0.56262 0.68392 0.79366 0.88495 0.95206 0.99078)

LEG=${LEGS[${SLURM_ARRAY_TASK_ID:-0}]}
D=$ROOT/data/ti_quad/$LEG
TOP=$D/ti.prmtop
[[ -f $TOP ]] || { echo "no ti.prmtop for $LEG"; exit 1; }

source <(sed 's/^/export /' "$D/masks.txt")

compress () {
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
NOSHAKE="${timask1},${timask2#@}"

echo "=== $LEG descending ladder  $(date -Is)"
echo "    timask1=$timask1"
echo "    timask2=$timask2"

PROD_NS=${PROD_NS:-10}
NPROD=$(python3 -c "print(int($PROD_NS*1000/0.002))")

run () {   # name prev workdir
    local nm=$1 prev=$2 wd=$3
    [[ -s $wd/${nm}.rst7 ]] && { echo "      $nm done"; return 0; }
    echo "      $nm ..."
    ( cd "$wd" && pmemd.cuda -O -i ${nm}.in -p "$TOP" -c "$prev" -ref "$prev" \
        -o ${nm}.out -r ${nm}.rst7 -inf ${nm}.mdinfo ) \
        || { echo "      !! $nm FAILED"; tail -25 "$wd/${nm}.out" 2>/dev/null; return 2; }
}

# rungs to climb down, and the window each is seeded from
for W in 3 2 1 0; do
    L=${LAM[$W]}
    WD=$D/lam$(printf '%02d' $W)
    SEED=$D/lam$(printf '%02d' $((W + 1)))/prod.rst7
    mkdir -p "$WD"
    echo "  --- window $W  lambda=$L   seed <- lam$(printf '%02d' $((W+1)))/prod.rst7"
    if [[ -s $WD/dvdl.txt ]]; then echo "      already complete, skip"; continue; fi
    [[ -s $SEED ]] || { echo "      !! seed $SEED missing -- the rung above must finish first"; exit 3; }

    ti_block () {
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

    # BACKBONE-only restraint throughout: the ligand and the sidechains must be
    # free to follow lambda down, which is the whole point of the ladder.
    cat > "$WD/min.in" <<IN
minimise at lambda $L, backbone restrained
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500, ntmin=2, ntb=1, cut=10.0, ntpr=500,
  ntr=1, restraintmask='@N,CA,C,O', restraint_wt=5.0,
$(ti_block)
 /
IN
    cat > "$WD/heat.in" <<IN
heat to 300 K, 200 ps, lambda $L
 &cntrl
  imin=0, irest=0, ntx=1, nstlim=100000, dt=0.002,
  ntb=1, ntp=0, ntc=2, ntf=1, cut=10.0,
  ntt=3, gamma_ln=2.0, ig=-1, tempi=10.0, temp0=300.0,
  ntr=1, restraintmask='@N,CA,C,O', restraint_wt=1.0,
  nmropt=1, ntpr=5000, ntwx=0, ntwr=50000,
$(ti_block)
 /
 &wt type='TEMP0', istep1=0, istep2=80000, value1=10.0, value2=300.0 /
 &wt type='END' /
IN
    cat > "$WD/eq.in" <<IN
NPT equilibration 500 ps, lambda $L
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=250000, dt=0.002,
  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  ntc=2, ntf=1, cut=10.0, ntt=3, gamma_ln=2.0, ig=-1, temp0=300.0,
  ntpr=5000, ntwx=0, ntwr=50000,
$(ti_block)
 /
IN
    cat > "$WD/prod.in" <<IN
TI production $PROD_NS ns, lambda $L
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=$NPROD, dt=0.002,
  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  ntc=2, ntf=1, cut=10.0, ntt=3, gamma_ln=2.0, ig=-1, temp0=300.0,
  ntpr=1000, ntwx=0, ntwr=50000, ntave=0,
 $(ti_block)
 /
IN

    run min  "$SEED"          "$WD" || exit 2
    run heat "$WD/min.rst7"   "$WD" || exit 2
    run eq   "$WD/heat.rst7"  "$WD" || exit 2
    run prod "$WD/eq.rst7"    "$WD" || exit 2

    # same de-duplicating parse as 89: pmemd prints every frame twice and the
    # trailing AVERAGES/RMS banners are not samples
    python3 - "$WD" <<'PY'
import re, statistics, sys, os
wd = sys.argv[1]
txt = open(os.path.join(wd, "prod.out")).read()
body = re.split(r"A V E R A G E S", txt)[0]
pairs = re.findall(r"NSTEP =\s*(\d+).*?DV/DL\s+=\s+(-?\d+\.\d+)", body, re.S)
seen = {}
for s, v in pairs:
    seen.setdefault(int(s), float(v))
v = [seen[k] for k in sorted(seen)]
h = len(v) // 2
open(os.path.join(wd, "dvdl.txt"), "w").write(
    f"n={len(v)}\nmean={statistics.mean(v):.4f}\nsd={statistics.stdev(v):.4f}\n"
    f"first_half={statistics.mean(v[:h]):.4f}\n"
    f"second_half={statistics.mean(v[h:]):.4f}\n")
print(f"      <dV/dl> = {statistics.mean(v):.3f}  "
      f"drift {statistics.mean(v[h:]) - statistics.mean(v[:h]):+.3f}")
PY
    echo "  --- window $W done $(date -Is)"
done
echo "=== $LEG ladder complete $(date -Is)"
