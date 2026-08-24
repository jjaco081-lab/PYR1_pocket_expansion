#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14   # amber22 has no kernels for k80/p100/h100/blackwell
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 3-00:00:00
#SBATCH -J qtiext
#SBATCH -a 0-23%4
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/qtiext_%a.log
#
# 105_ti_extend.sh -- take the crystal arm to the sampling its own error bar asked for.
#
# WHY
# S46a: the quadruple's crystal arm returned selectivity -1.12 +/- 2.34 -- the
# right sign, but the error is twice the effect, so it is a NO CALL rather than a
# result. The statistical term was stat = 1.74 from n_eff ~ 10 per window
# (tau ~ 245 frames = 0.5 ns over a 10 ns window). stat scales as 1/sqrt(N), so
# reaching 0.5 needs ~12x more sampling:
#
#     10 ns -> 120 ns per window, 24 windows, ~2.9 us, ~8.3 GPU-days
#
# Only the ABA and mandi_xtal legs are extended. The apo leg is NOT: it cancels
# exactly out of the selectivity, which is the quantity being resolved, and it
# was the worst-drifting leg (conv 4.61) so extending it would buy nothing that
# the cycle does not already discard. The docked arm is NOT extended either --
# it is structurally invalid (S46b), not undersampled, and more sampling of a
# ligand sitting inside a ghost phenylalanine would only make the wrong number
# more precise.
#
# WHAT THIS BUYS AND WHAT IT DOES NOT
# More sampling shrinks `stat`. It shrinks `conv` only if the windows genuinely
# equilibrate rather than drifting somewhere new -- which is exactly why 100
# reports the two separately and why the per-window drift table is printed. If
# conv does not come down with 12x the data, the answer is not undersampling.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22_mpi_cuda >/dev/null 2>&1

LEGS=(aba mandi_xtal)
LAM=(0.00922 0.04794 0.11505 0.20634 0.31608 0.43738 \
     0.56262 0.68392 0.79366 0.88495 0.95206 0.99078)

T=${SLURM_ARRAY_TASK_ID:-0}
LEG=${LEGS[$((T / 12))]}
LI=$((T % 12))
L=${LAM[$LI]}
D=$ROOT/data/ti_quad/$LEG
WD=$D/lam$(printf '%02d' $LI)
TOP=$D/ti.prmtop
[[ -f $TOP && -s $WD/prod.rst7 ]] || { echo "missing topology or prod.rst7 for $LEG lam$LI"; exit 1; }

TARGET_NS=${TARGET_NS:-120}
SEG_NS=${SEG_NS:-30}        # chunked so a failure costs one segment, not the run

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

echo "=== $LEG lam$(printf '%02d' $LI) lambda=$L  extend to ${TARGET_NS} ns  $(date -Is)"

# how much production exists already? MEASURED from the outputs, never counted.
# Frames are 2 ps apart (ntpr=1000 x dt=0.002).
have_ns () {
    env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/python - "$WD" <<'PYX'
import glob, os, re, sys
wd = sys.argv[1]
tot = 0
for f in sorted(glob.glob(os.path.join(wd, "prod*.out"))):
    body = re.split(r"A V E R A G E S", open(f).read())[0]
    steps = {int(s) for s, _ in
             re.findall(r"NSTEP =\s*(\d+).*?DV/DL\s+=\s+(-?\d+\.\d+)", body, re.S)}
    tot += len(steps)
print(f"{tot * 2 / 1000:.1f}")
PYX
}

NSEG=$(python3 -c "print(int($SEG_NS*1000/0.002))")
i=1
while :; do
    CUR=$(have_ns)
    ok=$(python3 -c "print(1 if $CUR >= $TARGET_NS - 0.5 else 0)")
    if [[ "$ok" == "1" ]]; then echo "  have ${CUR} ns >= ${TARGET_NS} ns, done"; break; fi
    SEG=$WD/prod$(printf '%02d' $i)
    if [[ -s ${SEG}.rst7 ]]; then i=$((i+1)); continue; fi
    PREV=$WD/prod.rst7
    j=$((i-1)); [[ $j -ge 1 ]] && PREV=$WD/prod$(printf '%02d' $j).rst7
    echo "  segment $i  (have ${CUR} ns)  from $(basename "$PREV")"
    cat > ${SEG}.in <<IN
TI production segment $i, $SEG_NS ns, lambda $L
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=$NSEG, dt=0.002,
  ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  ntc=2, ntf=1, cut=10.0, ntt=3, gamma_ln=2.0, ig=-1, temp0=300.0,
  ntpr=1000, ntwx=0, ntwr=50000, ntave=0,
$(ti_block)
 /
IN
    ( cd "$WD" && pmemd.cuda -O -i ${SEG}.in -p "$TOP" -c "$PREV" \
        -o ${SEG}.out -r ${SEG}.rst7 -inf ${SEG}.mdinfo ) || {
        echo "  !! segment $i FAILED"; tail -20 ${SEG}.out 2>/dev/null; exit 2; }
    i=$((i+1))
done
echo "=== $LEG lam$(printf '%02d' $LI) now $(have_ns) ns  $(date -Is)"
