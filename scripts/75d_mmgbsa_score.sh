#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -J mmscore
#SBATCH -a 0-13
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/mmscore_%a.log
#
# 75d_mmgbsa_score.sh -- the MM-GBSA step alone, on structures 75b already minimised.
#
# Split out because 75b's first run did ~90 minutes of correct minimisation and then
# fell over on one ante-MMPBSA flag. Redoing the minimisations to fix a scoring call
# would waste all of it, so the expensive stage is left on disk and only the cheap
# stage is repeated. 75b itself is fixed for future runs.
set -uo pipefail
module load amber/22
P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
D=$P/data/mmgbsa
ARMS=(mandi aba)
VARIANTS=(WT K59R K59Q K59N V81I F108A F159L)
ARM=${ARMS[$(( SLURM_ARRAY_TASK_ID / 7 ))]}
VAR=${VARIANTS[$(( SLURM_ARRAY_TASK_ID % 7 ))]}
W=$D/${ARM}_${VAR}
cd "$W" || exit 1
echo "=== $(date) $ARM / $VAR ==="

NSEED=$(ls s*_min.rst7 2>/dev/null | wc -l)
[[ $NSEED -ge 2 ]] || { echo "only $NSEED minimised seeds"; exit 1; }
echo "minimised seeds: $NSEED"

cat > cat.in <<CPP
parm s0.prmtop
$(for s in $(seq 0 $((NSEED-1))); do [[ -s s${s}_min.rst7 ]] && echo "trajin s${s}_min.rst7"; done)
trajout ensemble.nc netcdf
run
quit
CPP
cpptraj -i cat.in > cat.log 2>&1
[[ -s ensemble.nc ]] || { echo "cpptraj failed"; tail -5 cat.log; exit 1; }

ante-MMPBSA.py -p s0.prmtop -r rec.prmtop -l lig.prmtop -n ":LIG" > ante.log 2>&1
for f in rec.prmtop lig.prmtop; do
    [[ -s $f ]] || { echo "ante-MMPBSA failed to write $f"; tail -5 ante.log; exit 1; }
done

cat > mmgbsa.in <<MIN
MM-GBSA over independently repacked, minimised structures
 &general
  startframe=1, endframe=$NSEED, interval=1, verbose=2, keep_files=0,
 /
 &gb
  igb=8, saltcon=0.15,
 /
MIN
MMPBSA.py -O -i mmgbsa.in -o mmgbsa.dat -eo mmgbsa_frames.csv \
          -cp s0.prmtop -rp rec.prmtop -lp lig.prmtop -y ensemble.nc > mmpbsa.log 2>&1
rc=$?
[[ -s mmgbsa.dat ]] && grep -A2 "DELTA TOTAL" mmgbsa.dat | head -3 || { echo "MMPBSA FAILED"; tail -12 mmpbsa.log; }
echo "MMSCORE_DONE arm=$ARM var=$VAR rc=$rc frames=$NSEED"
date
exit $rc
