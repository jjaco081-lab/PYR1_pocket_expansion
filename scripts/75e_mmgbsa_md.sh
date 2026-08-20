#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 06:00:00
#SBATCH -J mmmd
#SBATCH -a 0-13
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/mmmd_%a.log
#
# 75e_mmgbsa_md.sh -- replace the fake ensemble with a real one.
#
# ⚠ WHY THIS EXISTS. The first attempt built 8 "independent" structures by repacking
# with 8 different Rosetta seeds and treated the spread as the error bar. It is not
# one: for the mandipropamid arm all 8 repacks came out BYTE-IDENTICAL, and the ABA
# arm collapsed to 2 distinct structures after minimisation. A 7-residue packing
# problem has one optimum and simulated annealing finds it every time, so the seeds
# sample nothing. Every sd of 0.00 in that table was an artefact of zero diversity,
# not a measurement of precision.
#
# Real conformational sampling instead: heat, then short Langevin MD in implicit
# solvent, snapshots every 5 ps. Those frames are genuinely different structures, so
# the spread across them means something.
#
# CUTOFFS: MD uses cut=16 for speed -- that affects SAMPLING only. MMPBSA computes
# each frame's energy with its own (effectively infinite) treatment, so the reported
# energies are not cut off at 16 A.
#
# ⚠ LENGTH, AND ITS HONEST LIMIT. GBn2 (igb=8) runs at only ~380 steps/min on one
# core for this 2,900-atom system -- measured, and pmemd is no faster than sander
# here. A 500 ps production would need ~11 h and the first submission was cancelled
# because it could not finish inside its wall. This is 50 ps with a snapshot every
# 0.5 ps: 100 GENUINELY DIFFERENT frames in ~80 min.
#
# 50 ps is short for a converged MM-GBSA average and that limit must travel with the
# numbers. What it does buy is the thing that was actually broken: an ensemble whose
# frames differ from each other at all. Eight byte-identical repacks gave sd = 0.00
# and no information; 100 frames over 50 ps give a real, if under-converged, spread.
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
[[ -s s0_min.rst7 ]] || { echo "no minimised start"; exit 1; }

cat > heat.in <<CIN
heat to 300 K over 10 ps, implicit solvent
 &cntrl
  imin=0, irest=0, ntx=1, nstlim=5000, dt=0.002,
  ntb=0, igb=8, saltcon=0.15, cut=16.0,
  ntc=2, ntf=2, ntt=3, gamma_ln=2.0, ig=-1,
  tempi=10.0, temp0=300.0, nmropt=1,
  ntpr=5000, ntwx=0, ntwr=5000,
 /
 &wt type='TEMP0', istep1=0, istep2=4000, value1=10.0, value2=300.0 /
 &wt type='END' /
CIN
cat > prod.in <<CIN
50 ps production, snapshot every 0.5 ps -> 100 frames
 &cntrl
  imin=0, irest=1, ntx=5, nstlim=25000, dt=0.002,
  ntb=0, igb=8, saltcon=0.15, cut=16.0,
  ntc=2, ntf=2, ntt=3, gamma_ln=2.0, ig=-1,
  temp0=300.0,
  ntpr=2500, ntwx=250, ntwr=5000,
 /
CIN

[[ -s heat.rst7 ]] || sander -O -i heat.in -p s0.prmtop -c s0_min.rst7 -o heat.out -r heat.rst7 >/dev/null 2>&1
[[ -s heat.rst7 ]] || { echo "heating FAILED"; tail -12 heat.out; exit 1; }
[[ -s md.nc ]] || sander -O -i prod.in -p s0.prmtop -c heat.rst7 -o md.out -r md.rst7 -x md.nc >/dev/null 2>&1
[[ -s md.nc ]] || { echo "production FAILED"; tail -12 md.out; exit 1; }
NF=$(cpptraj -p s0.prmtop -y md.nc -tl 2>/dev/null | grep -oE "[0-9]+$" | tail -1)
echo "frames: ${NF:-unknown}"

ante-MMPBSA.py -p s0.prmtop -r rec.prmtop -l lig.prmtop -n ":LIG" > ante2.log 2>&1
cat > mmgbsa_md.in <<MIN
MM-GBSA over a short implicit-solvent MD ensemble
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
[[ -s mmgbsa_md.dat ]] && grep -A2 "DELTA TOTAL" mmgbsa_md.dat | head -3 || { echo "MMPBSA FAILED"; tail -12 mmpbsa_md.log; }
echo "MMMD_DONE arm=$ARM var=$VAR rc=$rc"
date
exit $rc
