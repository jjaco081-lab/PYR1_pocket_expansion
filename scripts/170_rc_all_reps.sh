#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 6:00:00
#SBATCH -J rcall
#SBATCH -a 0-11
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/rcall_%a.log
#
# 170 -- measure the umbrella reaction coordinate (P88 CA - R116 CA) on EVERY
# factorial replicate, reading all preemption-continuation segments.
#
# README 95a: the original rc extraction covered 5 of 11 finished trajectories,
# because 110 read only prod.nc. README 62a's distribution-width result is two
# single-trajectory estimates as a consequence. This measures all of them so 62a
# can be recomputed at n=3/n=2 instead of n=1/n=1.
set -uo pipefail
R=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
module load amber/22
cd $R || exit 1
SYS=(S1_apo_open S2_holo_closed S9_apo_closed S10_holo_open)
S=${SYS[$((SLURM_ARRAY_TASK_ID / 3))]}
REP=rep$((SLURM_ARRAY_TASK_ID % 3))
D=$R/data/md191/$S/$REP
TOP=$R/data/md191/$S/system.prmtop
OUT=$R/results/factorial/rc
mkdir -p "$OUT"
[[ -f $TOP && -f $D/prod.nc ]] || { echo "missing inputs for $S/$REP"; exit 1; }
CONT=$(ls -1 "$D"/prod_cont_*.nc 2>/dev/null | sort)
{
  echo "parm $TOP"
  echo "trajin $D/prod.nc"
  for f in $CONT; do echo "trajin $f"; done
  echo "distance rc :88@CA :116@CA out $OUT/${S}_${REP}.dat"
  echo "go"
  echo "quit"
} > "$OUT/in_${S}_${REP}.in"
cpptraj -i "$OUT/in_${S}_${REP}.in" 2>&1 | tail -3
N=$(grep -vc '^#' "$OUT/${S}_${REP}.dat" 2>/dev/null || echo 0)
echo "=== $S/$REP  $(( 1 + $(echo \"$CONT\" | grep -c .) )) segments  $N frames"
