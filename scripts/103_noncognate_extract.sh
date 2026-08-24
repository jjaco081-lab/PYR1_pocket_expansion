#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 8:00:00
#SBATCH -J ncx
#SBATCH -a 0-11
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/ncx_%a.log
#
# 103_noncognate_extract.sh -- per-frame observables for the non-cognate MD.
#
# THE QUESTION (README 27, 43a)
# S24b licensed a STABILITY filter on the closed state. But S1 is apo AND open
# while S2 is holo AND closed, so conformation was confounded with occupancy. If
# closed PYR1 holds its state just as well around a ligand it was never built
# for, that filter is CONFORMATION-reporting and cannot rank ligands -- which
# disqualifies it as a design filter. Three ABA-sized ligands from three chemical
# classes, three independently docked poses each, against S2 (ABA) as reference.
#
# THE READOUT IS POSE STABILITY, NOT GATE OPENING. S24/S29 showed open and closed
# are both kinetically trapped at 300 ns and apo-closed does not open either, so
# opening is unobservable here for any ligand. Pose stability separated ABA from
# mandipropamid within 10 ns in S41-42.
#
# NUMBERING IS VERIFIED, NOT ASSUMED (58's header, feedback: verify residue identity)
# tleap renumbers sequentially, so native gate 85-89 is sequential 82-86 and
# native latch 115-117 is sequential 112-114 in these 178-residue systems.
# Checked by SEQUENCE in all twelve units before this ran: gate = SGLPA,
# latch = HRL. A wrong mask here returns RMSD for the wrong loop, silently.
#
# FIT BEFORE nofit. Ligand RMSD is measured AFTER superposing on the protein
# core. Measuring it without a prior fit reports whole-box tumbling instead --
# a mistake this project has made twice (S24, and the TI pose check in S44c).
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22

UNITS=(
  "S2_holo_closed/rep0"  "S2_holo_closed/rep1"  "S2_holo_closed/rep2"
  "S6_imperatorin/pose0" "S6_imperatorin/pose1" "S6_imperatorin/pose2"
  "S7_flutamide/pose0"   "S7_flutamide/pose1"   "S7_flutamide/pose2"
  "S8_estradiol/pose0"   "S8_estradiol/pose1"   "S8_estradiol/pose2"
)
U=${UNITS[${SLURM_ARRAY_TASK_ID:-0}]}
D=$ROOT/data/md/$U
OUT=$ROOT/results/noncognate/$(echo "$U" | tr / _)
mkdir -p "$OUT"

# S6-S8 keep a topology PER POSE (each pose was solvated separately, so the
# water count differs); S2 keeps ONE at the system level shared by its three
# seed replicates. Look in both places rather than assuming a layout.
TOP=$D/system.prmtop
[[ -f $TOP ]] || TOP=$(dirname "$D")/system.prmtop
TRJ=$D/prod.nc
[[ -f $TOP && -f $TRJ ]] || { echo "missing topology or trajectory for $U"; \
    echo "  looked for $D/system.prmtop and $(dirname "$D")/system.prmtop"; exit 1; }
echo "  topology: ${TOP#$ROOT/}"

# The ligand residue name DIFFERS per system -- A8S, IMP, FLU, EST -- so a
# hardcoded ":LIG" mask would silently select nothing and cpptraj would happily
# write a file of zeros. Read it out of the topology and assert it is one of the
# four expected, with the expected heavy-atom count.
read LIGRES NLIG <<< "$(env -u PYTHONPATH \
  /bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/python - "$TOP" <<'PYX'
import sys, parmed as pmd
EXPECT = {"A8S": 38, "IMP": 34, "FLU": 30, "EST": 44}
t = pmd.load_file(sys.argv[1])
lig = [r for r in t.residues if r.name in EXPECT]
if len(lig) != 1:
    sys.exit(f"expected exactly one ligand, found {[r.name for r in lig]}")
r = lig[0]
if len(r.atoms) != EXPECT[r.name]:
    sys.exit(f"{r.name} has {len(r.atoms)} atoms, expected {EXPECT[r.name]}")
print(r.name, sum(1 for a in r.atoms if a.atomic_number > 1))
PYX
)"
[[ -n "${LIGRES:-}" ]] || { echo "  !! could not identify the ligand in $TOP"; exit 1; }
echo "  ligand: $LIGRES ($NLIG heavy atoms)"

# S2 ran 300 ns, the non-cognates 150 ns. Compare over the COMMON window only --
# otherwise ABA gets twice the chance to drift and the contrast is rigged.
NFRAME=15000          # 150 ns at one frame per 10 ps

# sequential masks (see header)
GATE=82-86
LATCH=112-114
LB7A5=145-153         # native 148-156 - 3
# core = everything except the three mobile loops, backbone only
CORE=":1-81,87-111,115-144,154-178@CA,C,N,O"

cat > "$OUT/extract.in" <<CPPTRAJ
parm $TOP
trajin $TRJ 1 $NFRAME 1
autoimage

# 1. VALIDITY + equilibration: core backbone must settle
rms core $CORE first mass out $OUT/core_rmsd.dat

# 2. THE READOUT: ligand heavy-atom RMSD to its own starting pose, measured in
#    the core-superposed frame (the rms above already moved the coordinates)
rms lig :$LIGRES&!@H= first nofit out $OUT/lig_rmsd.dat

# 3. is the ligand still engaged with the protein at all?
nativecontacts :$LIGRES&!@H= ":1-178&!@H=" distance 4.0 first \
    out $OUT/lig_contacts.dat
# centre-of-mass drift, same superposed frame
vector ligcom center :$LIGRES out $OUT/lig_com.dat

# 4. the S24b filter itself -- does it move differently per ligand?
rms gate :$GATE@CA,C,N,O first nofit out $OUT/gate_rmsd.dat
rms latch :$LATCH@CA,C,N,O first nofit out $OUT/latch_rmsd.dat
rms lb7a5 :$LB7A5@CA,C,N,O first nofit out $OUT/lb7a5_rmsd.dat

# 5. the gate-latch staple: minimum heavy-atom distance between the two loops
nativecontacts ":$GATE&!@H=" ":$LATCH&!@H=" distance 4.0 first \
    out $OUT/gatelatch.dat

# 6. per-residue fluctuation, for the loops
atomicfluct out $OUT/rmsf.dat :1-178@CA byres

go
quit
CPPTRAJ

echo "=== $U  $(date -Is)"
cpptraj -i "$OUT/extract.in" > "$OUT/extract.log" 2>&1
rc=$?
if [[ $rc -ne 0 ]] || [[ ! -s $OUT/lig_rmsd.dat ]]; then
    echo "  !! cpptraj failed for $U"
    grep -iE "error|Error" "$OUT/extract.log" | head -5
    exit 2
fi
echo "  frames read: $(grep -oP 'Read \K[0-9]+(?= frames)' "$OUT/extract.log" | tail -1)"
echo "  wrote $(ls "$OUT"/*.dat | wc -l) data files"
echo "=== $U done $(date -Is)"
