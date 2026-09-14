#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH -J pcx
#SBATCH -a 0-17
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/pcx_%a.log
#
# 184 -- ligand RMSD to the STARTING pose over 20 ns, for the six pose-check
# systems. The question (README 107a): does MD retain a predicted pose, or shed
# it? A forced pose should reorient; a correct one should hold.
#
# Superposition is on the protein CORE ONLY (gate 85-89 and latch 115-117
# excluded) so ligand-driven loop motion cannot absorb ligand displacement into
# the fit. RMSD is then ligand heavy atoms, no re-fitting of the ligand itself,
# so it reports where the ligand SITS in the pocket frame.
#
# ⚠ Reads prod.nc AND every prod_cont_*.nc, in order. README 95a: the factorial
# extraction hardcoded prod.nc, silently analysed first segments, and produced a
# false "n=1" that propagated into four sections and one retracted result.
set -uo pipefail
R=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
module load amber/22 2>/dev/null
cd $R || exit 1
TAGS=(aba_crystal mandi_crystal mandi_predicted win_crystal fludioxonil_pred benoxacor_pred)
LIGS=(ABX MDX MDY WNX FLX BNY)
I=$(( SLURM_ARRAY_TASK_ID / 3 ))
T=${TAGS[$I]}; RN=${LIGS[$I]}
REP=rep$(( SLURM_ARRAY_TASK_ID % 3 ))
D=$R/data/posecheck/$T/$REP
TOP=$R/data/posecheck/$T/system.prmtop
OUT=$R/results/posecheck/${T}_${REP}
mkdir -p "$OUT"
[[ -f $TOP && -f $D/prod.nc ]] || { echo "missing inputs for $T/$REP"; exit 1; }
CONT=$(ls -1 "$D"/prod_cont_*.nc 2>/dev/null | sort)
NSEG=$(( 1 + $(echo "$CONT" | grep -c .) ))
# reference = the minimised/equilibrated start, i.e. the pose as built
REF=$D/eq2.rst7
[[ -s $REF ]] || REF=$D/min2.rst7
{
  echo "parm $TOP"
  echo "reference $REF [start]"
  echo "trajin $D/prod.nc"
  for f in $CONT; do echo "trajin $f"; done
  # ⚠ AUTOIMAGE FIRST. Without it, the moment the solute diffuses across a
  # periodic boundary the ligand is measured against its WRAPPED copy and the
  # RMSD jumps discontinuously to ~70-86 A. The first cut of this script omitted
  # it and reported crystal ABA "leaving the pocket" at 8.8 ns -- a jump of
  # 2.1 -> 77.4 A in 10 ps, which no ligand can do. Anchoring on the protein
  # keeps the complex together across the boundary.
  echo "autoimage anchor :6-181"
  # core = everything except the mobile gate/latch, backbone only
  echo "rms core :6-84,90-114,118-181@CA,C,N,O ref [start] out $OUT/core.dat"
  # ligand RMSD in that frame -- nofit so the ligand is NOT re-superposed
  echo "rms lig :$RN&!@H= ref [start] nofit out $OUT/lig_rmsd.dat"
  echo "rms gate :85-89@CA,C,N,O ref [start] nofit out $OUT/gate.dat"
  echo "distance d :88@CA :116@CA out $OUT/rc.dat"
  echo "go"
  echo "quit"
} > $OUT/extract.in
cpptraj -i $OUT/extract.in > $OUT/extract.log 2>&1
N=$(grep -vc '^#' $OUT/lig_rmsd.dat 2>/dev/null || echo 0)
echo "=== $T/$REP  ligand $RN  $NSEG segment(s)  $N frames"
tail -3 $OUT/extract.log | grep -i error && exit 1
exit 0
