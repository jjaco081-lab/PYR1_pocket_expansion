#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 8:00:00
#SBATCH -J factx
#SBATCH -a 0-11
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/factx_%a.log
#
# 110_factorial_extract.sh -- the conformation x occupancy factorial (README 30).
#
# THE QUESTION THIS ANSWERS, and why it is now the important one
# Pocket expansion has a mechanistic problem that is prior to any scoring method:
# the gate closes ONTO the ligand. Enlarge the cavity without enlarging the
# ligand and the ligand no longer reaches the gate, giving either
#   (a) no closure           -> no signal, or
#   (b) ligand-INDEPENDENT closure -> constitutive activity, which is exactly what
#       the counter-selection removes.
# So the design question is whether an under-filled cavity can hold the closed
# state at all, and at what energetic cost.
#
# The 2x2 that addresses it has been sitting complete and unanalysed:
#
#                   apo                     + ABA
#   open      S1  (resting state)      S10  (does the ligand CLOSE it?)
#   closed    S9  (EMPTY but closed)   S2   (the reference)
#
#   S9 vs S2 is the whole question: what does the ligand contribute to holding
#   the closed state? If S9 holds as well as S2, closure does not need filling --
#   which is permissive for an expanded pocket but raises the constitutive risk.
#   If S9 opens, filling matters, and an under-filled expanded pocket loses signal.
#
# All four cells were rebuilt on one protocol (README 30a) because the old-tree S9
# was KCl while S1/S2 were NaCl, and pairing them would have confounded ligand
# removal with a cation swap inside the one comparison the factorial exists to make.
#
# NUMBERING: md191 is the 191-residue rebuild and uses NATIVE numbering directly --
# verified by sequence here, gate 85-89 = SGLPA and latch 115-117 = HRL. That is
# NOT true of the older data/md tree (178 residues, gate = sequential 82-86), so
# the masks are not interchangeable between trees.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22

SYS=(S1_apo_open S2_holo_closed S9_apo_closed S10_holo_open)
T=${SLURM_ARRAY_TASK_ID:-0}
S=${SYS[$((T / 3))]}
R=rep$((T % 3))
D=$ROOT/data/md191/$S/$R
TOP=$ROOT/data/md191/$S/system.prmtop
TRJ=$D/prod.nc
[[ -f $TOP && -f $TRJ ]] || { echo "missing $TOP or $TRJ"; exit 1; }
# ⚠ PREEMPTION CONTINUATIONS. Every run on preempt_gpu that was interrupted
# resumes into prod_cont_NN.nc. Reading only prod.nc silently analyses the FIRST
# SEGMENT and calls it the run -- which is what produced README 58b's "the
# factorial is n=1". 11 of 12 replicates actually finished 300 ns; 7 of them
# have continuation segments. The sister script 58b_loop_dynamics_run.sh:115
# already globbed these correctly, so this was a regression, not an oversight.
CONT=$(ls -1 "$D"/prod_cont_*.nc 2>/dev/null | sort)
NSEG=$(( 1 + $(echo "$CONT" | grep -c . ) ))
OUT=$ROOT/results/factorial/${S}_${R}
mkdir -p "$OUT"

CLOSED=$ROOT/data/pyr1_closed_A.pdb
OPEN=$ROOT/data/pyr1_open_A.pdb

GATE=85-89
LATCH=115-117
# core excludes the three mobile loops so they cannot drag the superposition
CORE=":1-84,90-114,118-147,157-181@CA,C,N,O"

# the ligand residue name differs by cell (apo cells have none); read it, never assume
LIG=$(env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/python - "$TOP" <<'PYX'
import sys, parmed as pmd
t = pmd.load_file(sys.argv[1])
l = [r.name for r in t.residues if r.name in ("A8S",)]
print(l[0] if l else "")
PYX
)
echo "=== $S/$R   ligand='${LIG:-none}'   ${NSEG} trajectory segment(s)   $(date -Is)"

{
  echo "parm $TOP"
  # each reference PDB needs its OWN parm -- it is a 1409-atom protein-only file
  # and cpptraj otherwise reads it against the 46,570-atom solvated topology
  echo "parm $CLOSED [pclosed]"
  echo "parm $OPEN [popen]"
  echo "reference $CLOSED parm [pclosed] [closed]"
  echo "reference $OPEN parm [popen] [open]"
  echo "trajin $TRJ"
  # continuation segments, in order -- see the NSEG comment above
  for f in $CONT; do echo "trajin $f"; done
  echo "autoimage"
  # validity + equilibration, from the CORE only
  echo "rms core $CORE first mass out $OUT/core_rmsd.dat"
  # the open/closed coordinate: gate and latch against BOTH references.
  # fit on the core each time so the loop RMSD is not absorbed by the fit.
  echo "rms fitC $CORE ref [closed]"
  echo "rms g2c :$GATE@CA,C,N,O ref [closed] nofit out $OUT/gate_to_closed.dat"
  echo "rms l2c :$LATCH@CA,C,N,O ref [closed] nofit out $OUT/latch_to_closed.dat"
  echo "rms fitO $CORE ref [open]"
  echo "rms g2o :$GATE@CA,C,N,O ref [open] nofit out $OUT/gate_to_open.dat"
  echo "rms l2o :$LATCH@CA,C,N,O ref [open] nofit out $OUT/latch_to_open.dat"
  # the gate-latch staple
  echo "nativecontacts \":$GATE&!@H=\" \":$LATCH&!@H=\" distance 4.0 first out $OUT/gatelatch.dat"
  if [[ -n "$LIG" ]]; then
      echo "rms fitL $CORE first"
      echo "rms lig :$LIG&!@H= first nofit out $OUT/lig_rmsd.dat"
      echo "nativecontacts \":$LIG&!@H=\" \":1-191&!@H=\" distance 4.0 first out $OUT/lig_contacts.dat"
  fi
  echo "atomicfluct out $OUT/rmsf.dat :1-191@CA byres"
  echo "go"
  echo "quit"
} > "$OUT/extract.in"

cpptraj -i "$OUT/extract.in" > "$OUT/extract.log" 2>&1
rc=$?
if [[ $rc -ne 0 || ! -s $OUT/gate_to_closed.dat ]]; then
    echo "  !! cpptraj failed"; grep -iE "error" "$OUT/extract.log" | head -5; exit 2
fi
echo "  frames: $(grep -oP 'Read \K[0-9]+(?= frames)' "$OUT/extract.log" | tail -1)"
echo "=== $S/$R done $(date -Is)"
