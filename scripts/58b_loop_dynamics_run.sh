#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH -t 08:00:00
#SBATCH -J loopdyn
#SBATCH -a 0-6
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/loopdyn_%a.log
#
# 58b_loop_dynamics_run.sh -- gate/latch loop dynamics on the finished WT MD.
#
# THE QUESTION
#   Does the open state stay open and the closed state stay closed over 300 ns,
#   and do the two have separable loop dynamics we can later use to filter
#   pocket-expanded designs?
#
# WHAT IS COMPUTED -- ONLY the observables pre-registered in README 19d, which
# were fixed before these runs finished so the answer cannot be shopped for:
#     gate (85-89) and latch (115-117) backbone RMSD and RMSF
#     gate-latch minimum heavy-atom distance
#     Lbeta7-alpha5 (148-156) RMSF          [Dorosh 2013, README 17a]
#
# Plus two things that are NOT filter observables and are labelled as such:
#     - whole-core backbone RMSD vs frame 1, used ONLY to pick the equilibration
#       discard window. Chosen from the core, never from the loops being tested,
#       so the discard cannot be tuned to flatter the answer.
#     - ABA heavy-atom RMSD in S2/S4, a VALIDITY CONTROL. If the ligand leaves the
#       pocket then "closed" no longer means what we think and S2 is void.
#
# NUMBERING -- THE TRAP THIS SCRIPT EXISTS TO AVOID
#   tleap renumbers residues 1..N, and the crystal files have gaps, so the
#   native-to-sequential offset is NOT constant (it takes values 0, 1 and 3).
#   A mask of ':85-89' does NOT select the gate and would fail silently.
#
#   Worse, THE MAP IS NOT THE SAME IN EVERY SYSTEM. S1/S2 come from script 30,
#   which intersected 3K3K and 3QN1 and so dropped residue 2; S4 was built from
#   3QN1 alone and keeps it (as ALA, the P2A). So the gate is sequential 82-86 in
#   S1/S2 but 83-87 in S4. Masks below come from 58_loop_dynamics_prep.py, which
#   rebuilds the map per system and verifies every loop by residue identity
#   (gate = SER GLY LEU PRO ALA, latch = HIS ARG LEU). Do not hand-edit them.
#
# TWO-REFERENCE PROJECTION
#   Each frame is superposed on the RIGID CORE (161 residues, gate/latch/Lb7a5
#   excluded so their motion cannot leak into the fit) and then the loop RMSD is
#   measured to BOTH the open and the closed crystal reference with `nofit`.
#   A frame that is genuinely open sits near 0 on the open axis and near 5 on the
#   closed axis, and vice versa. The `rms` fit moves the coordinates and the
#   following `nofit` rmsd measures in that frame -- this is the standard cpptraj
#   idiom and the order below is deliberate.
set -euo pipefail
module load amber/22

P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
MD=$P/data/md
OUT=$P/data/loop_dynamics
MAP=$OUT/residue_map.json
[[ -s "$MAP" ]] || { echo "run 58_loop_dynamics_prep.py first"; exit 1; }

SYSREP=(S1_apo_open:0 S1_apo_open:1 S1_apo_open:2
        S2_holo_closed:0 S2_holo_closed:1 S2_holo_closed:2
        S4_ternary:0)
entry=${SYSREP[$SLURM_ARRAY_TASK_ID]}
SYS=${entry%%:*}
REP=${entry##*:}

# Masks come from the verified per-system map, never retyped and never shared
# between systems. S4 carries 3QN1's residue 2 (the P2A) that the script-30
# intersection dropped from S1/S2, so its gate is :83-87 where S1/S2 have :82-86.
# The REFERENCE pdbs keep the script-30 numbering, so trajectory masks and
# reference masks are read separately and passed to cpptraj as `<mask> <refmask>`.
q() { python3 -c "import json,sys;m=json.load(open('$MAP'));print(m['systems'][sys.argv[1]][sys.argv[2]] if sys.argv[2]!='m' else m['systems'][sys.argv[1]]['masks_sequential'][sys.argv[3]])" "$@"; }
GATE=$(q "$SYS" m gate)
LATCH=$(q "$SYS" m latch)
LB=$(q "$SYS" m lb7a5)
CORE=$(q "$SYS" core_mask_sequential)
RGATE=$(q _reference m gate)
RLATCH=$(q _reference m latch)
RCORE=$(q _reference core_mask_sequential)
echo "traj masks: gate=:$GATE latch=:$LATCH lb7a5=:$LB"
NRES=$(q "$SYS" n_residues)
echo "ref  masks: gate=:$RGATE latch=:$RLATCH   PYR1 residues=$NRES"
D=$MD/$SYS/rep$REP
O=$OUT/$SYS/rep$REP
mkdir -p "$O"

[[ -s $D/prod.nc ]] || { echo "no production trajectory in $D -- skipping"; exit 0; }

echo "=== $SYS rep$REP  $(date) ==="

# trajectory list: prod.nc then any preemption-continuation segments, in order
TRAJIN="trajin $D/prod.nc parm [sysparm]"
for f in $(ls -1 $D/prod_cont_*.nc 2>/dev/null | sort); do
    TRAJIN="$TRAJIN
trajin $f parm [sysparm]"
done

# ABA validity control, for the systems that have a ligand (S2 and S4).
#
# The ligand is addressed BY RESIDUE NAME, `:A8S`, not by an index. An earlier
# version computed it as n_PYR1_residues+1, which is right for S2 (179) and badly
# wrong for S4, where residue 180 is HAB1's first residue and A8S is 478. That
# would have reported an HAB1 side chain as the ligand without any error. The
# name mask cannot drift when a system gains a chain.
#
# NOTE the reference is `first`, NOT the open/closed reference PDBs: script 30
# deliberately wrote both references APO, so they contain no ligand at all and an
# rmsd against them would be meaningless. Measuring against frame 1 immediately
# after the core-to-frame-1 fit is the right control anyway -- it asks "has ABA
# moved relative to the protein it started in", which is what decides whether the
# system is still holo.
LIGBLOCK=""
if [[ "$SYS" == "S2_holo_closed" || "$SYS" == "S4_ternary" ]]; then
LIGBLOCK="rmsd aba_rmsd :A8S&!@H= first nofit out $O/aba_rmsd.dat
distance aba_gate :A8S&!@H= :$GATE&!@H= out $O/aba_gate_dist.dat"
fi

# Each reference PDB needs its OWN parm. cpptraj otherwise reads the reference
# against the last-loaded topology and dies with
#   "No frames read. atom=1428 expected 46820"
# because the solvated system has 46,820 atoms and the reference has 1,428.
# The two references also differ from each other (1428 vs 1409 atoms: script 30
# matched residues, not every side-chain atom), which is harmless here because
# every mask below is backbone-only, and script 30 guarantees a complete N/CA/C/O
# for all 178 residues in both.
#
# cpptraj renumbers residues sequentially when it reads a PDB, exactly as tleap
# did for the prmtop. That makes the REFERENCE masks (:82-86 etc) fixed, but the
# TRAJECTORY masks are per system -- S4 is shifted by one -- so each rmsd below
# passes both: `rmsd <name> <traj mask> <ref mask> ref [...]`. Verified with
# `resinfo`, whose #Orig column shows reference sequential 82 = native 85 SER and
# 112-114 = 115-117 HIS ARG LEU.
cat > "$O/analysis.in" <<EOF
parm $MD/$SYS/system.prmtop [sysparm]
parm $P/data/pyr1_open_A.pdb   [openparm]
parm $P/data/pyr1_closed_A.pdb [closedparm]
reference $P/data/pyr1_open_A.pdb   parm [openparm]   [open_ref]
reference $P/data/pyr1_closed_A.pdb parm [closedparm] [closed_ref]
$TRAJIN

autoimage

# ---- equilibration diagnostic (NOT a filter observable) ----
# ACTION ORDER MATTERS: each rms re-fits the frame, and every following nofit
# measurement is taken in whatever frame the last fit established. The ABA control
# is placed directly after this fit-to-frame-1 on purpose.
# (NO BACKTICKS anywhere below this line: the heredoc is unquoted so that \$GATE
#  and friends expand, which means bash also COMMAND-SUBSTITUTES backticks. A
#  comment reading 'nofit' in backticks made bash run nofit and silently delete
#  the word from the generated cpptraj input.)
rms core_vs_first :$CORE@N,CA,C,O first out $O/core_rmsd_first.dat mass

$LIGBLOCK

# ---- two-reference projection: fit on core, then measure loops with nofit ----
rms fit_open :$CORE@N,CA,C,O :$RCORE@N,CA,C,O ref [open_ref] out $O/core_fit_open.dat
rmsd gate_to_open  :$GATE@N,CA,C,O  :$RGATE@N,CA,C,O  ref [open_ref] nofit out $O/gate_to_open.dat
rmsd latch_to_open :$LATCH@N,CA,C,O :$RLATCH@N,CA,C,O ref [open_ref] nofit out $O/latch_to_open.dat

rms fit_closed :$CORE@N,CA,C,O :$RCORE@N,CA,C,O ref [closed_ref] out $O/core_fit_closed.dat
rmsd gate_to_closed  :$GATE@N,CA,C,O  :$RGATE@N,CA,C,O  ref [closed_ref] nofit out $O/gate_to_closed.dat
rmsd latch_to_closed :$LATCH@N,CA,C,O :$RLATCH@N,CA,C,O ref [closed_ref] nofit out $O/latch_to_closed.dat

# ---- gate-latch minimum heavy-atom distance (pre-registered) ----
nativecontacts :$GATE&!@H= :$LATCH&!@H= mindist maxdist \
    out $O/gate_latch_contacts.dat distance 4.5

# ---- RMSF: fit to the average structure, not to frame 1 ----
rms fit_for_avg :$CORE@N,CA,C,O first
average crdset AVGSTRUCT
run

rms fit_to_avg :$CORE@N,CA,C,O ref AVGSTRUCT
atomicfluct out $O/rmsf_byres.dat :1-$NRES&!@H= byres
atomicfluct out $O/rmsf_bb_byres.dat :1-$NRES@N,CA,C,O byres
run
quit
EOF

# `set -e` would abort here before rc is captured, so take the exit code inline.
rc=0
cpptraj -i "$O/analysis.in" > "$O/cpptraj.log" 2>&1 || rc=$?

# cpptraj returns 0 on some errors, so verify by artefact
fail=0
EXPECT_FILES="core_rmsd_first gate_to_open gate_to_closed latch_to_open latch_to_closed rmsf_byres rmsf_bb_byres gate_latch_contacts"
# the ligand control is only expected where there IS a ligand
if [[ -n "$LIGBLOCK" ]]; then EXPECT_FILES="$EXPECT_FILES aba_rmsd aba_gate_dist"; fi
for f in $EXPECT_FILES; do
    if [[ ! -s "$O/$f.dat" ]]; then echo "MISSING $f.dat"; fail=1; fi
done
# every .dat must have real rows, not just a header
for f in $EXPECT_FILES; do
    n=$(grep -vc '^#' "$O/$f.dat" 2>/dev/null || echo 0)
    if [[ "$n" -lt 10 ]]; then echo "TOO FEW ROWS in $f.dat ($n)"; fail=1; fi
done
if grep -qiE "^Error|Could not|not found" "$O/cpptraj.log"; then
    echo "--- cpptraj reported errors ---"
    grep -iE "^Error|Could not|not found" "$O/cpptraj.log" | head
    fail=1
fi

NFRAME=$(awk 'END{print NR-1}' "$O/gate_to_open.dat" 2>/dev/null || echo 0)
echo "frames analysed: $NFRAME"
echo "LOOPDYN_DONE sys=$SYS rep=$REP rc=$rc fail=$fail frames=$NFRAME"
date
exit $fail
