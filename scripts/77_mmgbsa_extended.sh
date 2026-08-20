#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 16:00:00
#SBATCH -J mmext
#SBATCH -a 0-45
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/mmext_%a.log
#
# 77_mmgbsa_extended.sh -- MM-GBSA over the STRATIFIED pairwise candidate set,
# with ensembles 5x longer than README 37's. One task per variant x ligand arm.
#
# WHAT CHANGED FROM 75
#   candidates  4 known singles -> 23 variants (README 32's ranking, stratified by
#               POSITION PAIR so that F108, which relieves the largest clash and
#               fills all 24 of the raw top pairs, cannot saturate the set)
#   ensemble    50 ps -> 250 ps, still 100 frames but sampled every 2.5 ps, so the
#               frames are less correlated with each other
#   seeds       8 identical repacks -> 1 structure (README 37c: they were byte-
#               identical, and the ensemble comes from MD)
#
# WHY CPU AND NOT GPU
#   pmemd.cuda is available and would be far faster per system, but this account is
#   capped at 1 concurrent GPU on preempt_gpu and 4 on gpu, while cutlerlab has
#   already run 100 concurrent CPU tasks. 46 tasks x ~6 h wall beats 46 tasks queued
#   behind a cap of 4. Throughput here is a concurrency problem, not a speed one.
set -uo pipefail
# ⚠ EVERY Amber input below starts with a TITLE LINE. sander does not treat it as
# optional: without it the parser reports 'Could not find cntrl namelist' and exits,
# which killed all 46 tasks on the previous submission.
module load amber/22
P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
D=$P/data/mmgbsa
PAR=$D/params
PY=/bigdata/cutlerlab/jjaco081/conda_envs/pyr1_docking/bin/python

# The ordered variant list is written by 75_mmgbsa_build.py. It is NOT rebuilt here:
# the previous version regenerated the names inline with an f-string containing a
# regex backslash, which Python 3.10 rejects, and all 46 tasks died before starting.
VLIST=$D/variants.txt
[[ -s $VLIST ]] || { echo "missing $VLIST -- run 75_mmgbsa_build.py first"; exit 1; }
mapfile -t VARIANTS < "$VLIST"
NV=${#VARIANTS[@]}
ARMS=(mandi aba)
ARM=${ARMS[$(( SLURM_ARRAY_TASK_ID / NV ))]}
VAR=${VARIANTS[$(( SLURM_ARRAY_TASK_ID % NV ))]}
case $ARM in mandi) LIG=3UZ ;; aba) LIG=A8S ;; esac
W=$D/${ARM}_${VAR}
echo "=== $(date) task $SLURM_ARRAY_TASK_ID -> $ARM / $VAR (of $NV variants) ==="
[[ -s $W/seed0_protein.pdb ]] || { echo "no build for $W"; exit 1; }
cd "$W" || exit 1

# ---- topology (mbondi3: igb=8 refuses to run without it) ----
if [[ ! -s s0.prmtop ]]; then
  pdb4amber -i seed0_protein.pdb -o s0_clean.pdb --dry --nohyd > p4a.log 2>&1
  cat > leap.in <<LIN
source leaprc.protein.ff14SB
source leaprc.gaff2
set default PBRadii mbondi3
loadamberparams $PAR/${LIG}.frcmod
LIG = loadmol2 $PAR/${LIG}.mol2
prot = loadpdb s0_clean.pdb
com = combine { prot LIG }
saveamberparm com s0.prmtop s0.inpcrd
quit
LIN
  tleap -f leap.in > leap.log 2>&1
fi
[[ -s s0.prmtop ]] || { echo "tleap FAILED"; tail -6 leap.log; exit 1; }

cat > min.in <<CIN
minimise in implicit solvent
 &cntrl
  imin=1, maxcyc=2000, ncyc=1000, ntb=0, igb=8, saltcon=0.15, cut=999.0, ntpr=500,
 /
CIN
cat > heat.in <<CIN
heat to 300 K over 10 ps, implicit solvent
 &cntrl
  imin=0, irest=0, ntx=1, nstlim=5000, dt=0.002, ntb=0, igb=8, saltcon=0.15,
  cut=16.0, ntc=2, ntf=2, ntt=3, gamma_ln=2.0, ig=-1,
  tempi=10.0, temp0=300.0, nmropt=1, ntpr=2500, ntwx=0, ntwr=5000,
 /
 &wt type='TEMP0', istep1=0, istep2=4000, value1=10.0, value2=300.0 /
 &wt type='END' /
CIN
cat > prod.in <<CIN
250 ps, snapshot every 2.5 ps -> 100 less-correlated frames
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=125000, dt=0.002, ntb=0, igb=8, saltcon=0.15,
  cut=16.0, ntc=2, ntf=2, ntt=3, gamma_ln=2.0, ig=-1, temp0=300.0,
  ntpr=12500, ntwx=1250, ntwr=25000,
 /
CIN

[[ -s s0_min.rst7 ]] || sander -O -i min.in  -p s0.prmtop -c s0.inpcrd    -o min.out  -r s0_min.rst7 >/dev/null 2>&1
[[ -s s0_min.rst7 ]] || { echo "minimisation FAILED"; tail -8 min.out; exit 1; }
[[ -s heat.rst7   ]] || sander -O -i heat.in -p s0.prmtop -c s0_min.rst7  -o heat.out -r heat.rst7   >/dev/null 2>&1
[[ -s heat.rst7   ]] || { echo "heating FAILED"; tail -8 heat.out; exit 1; }
[[ -s md.nc       ]] || sander -O -i prod.in -p s0.prmtop -c heat.rst7    -o md.out   -r md.rst7 -x md.nc >/dev/null 2>&1
[[ -s md.nc       ]] || { echo "production FAILED"; tail -8 md.out; exit 1; }

# ⚠ no -c: ante-MMPBSA writes no complex topology from an already-unsolvated input
ante-MMPBSA.py -p s0.prmtop -r rec.prmtop -l lig.prmtop -n ":LIG" > ante.log 2>&1
for f in rec.prmtop lig.prmtop; do
  [[ -s $f ]] || { echo "ante-MMPBSA failed on $f"; tail -5 ante.log; exit 1; }
done
cat > mmgbsa_md.in <<MIN
 &general
  startframe=1, endframe=9999, interval=1, verbose=2, keep_files=0,
 /
 &gb
  igb=8, saltcon=0.15,
 /
MIN
MMPBSA.py -O -i mmgbsa_md.in -o mmgbsa_md.dat -eo mmgbsa_md_frames.csv \
          -cp s0.prmtop -rp rec.prmtop -lp lig.prmtop -y md.nc > mmpbsa_md.log 2>&1
rc=$?
rm -f reference.frc
[[ -s mmgbsa_md.dat ]] && grep -A2 "DELTA TOTAL" mmgbsa_md.dat | head -3 || { echo "MMPBSA FAILED"; tail -10 mmpbsa_md.log; }
echo "MMEXT_DONE arm=$ARM var=$VAR rc=$rc"
date
exit $rc
