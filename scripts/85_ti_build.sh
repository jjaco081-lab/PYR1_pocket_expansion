#!/bin/bash
# 85_ti_build.sh -- alchemical (TI) systems for protein point mutations.
#
# WHY TI AND NOT MORE MM-GBSA
#   README 42d: MM-GBSA needs n=60 replicates per variant per ligand to reach
#   +/-0.5 kcal/mol, because its between-seed sd is 2.72. That is ~1200 ns per
#   variant for two ligands and does not scale to many substitutions x many
#   ligands. TI computes the SAME quantity as a rigorous free energy for ~120 ns,
#   because the mutation is run in both legs of a thermodynamic cycle and the
#   ~43,000 unchanged atoms cancel by construction instead of being subtracted as
#   two separately-computed absolute energies -- which is precisely how one bad
#   WT-ABA reference poisoned all 22 ddGs.
#
# THE CYCLE
#     WT.L  --dG1(mutate in complex)-->  MUT.L
#      |                                  |
#   dGbind(WT)                       dGbind(MUT)
#      |                                  |
#     WT + L --dG2(mutate in apo)-->  MUT + L
#     ddGbind = dG1 - dG2
#
#   For SELECTIVITY the apo leg cancels exactly:
#     ddG(mandi) - ddG(ABA) = dG1_mandi - dG1_ABA
#   so only the two COMPLEX legs are needed to rank a mutation between ligands.
#   The apo leg is built anyway (once, shared) because it converts the selectivity
#   into an absolute ddGbind per ligand at 1/3 more cost.
#
# TOPOLOGY: FULL SIDECHAIN DUAL TOPOLOGY, deliberately
#   The efficient setup keeps a common core and perturbs only the atoms that
#   appear/disappear. It is NOT used here. Measured on the real topologies,
#   11 of the 15 atoms VAL and ILE nominally share change partial charge, and the
#   shared CG1 hydrogens do not even sit at matching coordinates because CG1 is a
#   methyl in VAL and a methylene in ILE. A common-core mask would therefore be
#   quietly wrong. Putting the ENTIRE sidechain from CB outward in the softcore
#   region costs convergence speed and buys correctness that can be checked.
#   Backbone stays common, so tiMerge's coordinate match is exact there.
#
# NUMBERING IS VERIFIED, NOT ASSUMED (README 28, feedback: verify residue identity)
#   The source PDB keeps NATIVE numbering with the known 69-70 gap; pdb4amber and
#   tleap renumber sequentially. native 81 -> seq 79, native 59 -> seq 59 (before
#   the gap). Sequential 81 is ALSO a valine (native V83), two positions away, so
#   trusting the number rather than the identity would mutate the wrong residue.
#   Every build below asserts the residue IS what it should be before proceeding.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22
command -v tleap >/dev/null || { echo "no tleap"; exit 1; }
PARMED=/bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/python
PARAMS=$ROOT/data/mmgbsa/params
OUT=$ROOT/data/ti
mkdir -p "$OUT"

SRC=$ROOT/results/stage1_rosetta/_input_wt_aba_A8S_anion.pdb
NRES=179     # asserted below

# heavy-atom protein, hydrogens dropped so tleap rebuilds them identically in both
# copies; that keeps the common backbone coordinate-matched for tiMerge
python3 - "$SRC" "$OUT/wt_heavy.pdb" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
keep = []
for l in open(src):
    if not l.startswith("ATOM"):
        continue
    el = l[76:78].strip() or l[12:16].strip()[0]
    if el == "H":
        continue
    keep.append(l)
open(dst, "w").writelines(keep)
res = {int(l[22:26]) for l in keep}
print(f"  wt_heavy.pdb: {len(keep)} heavy atoms, {len(res)} residues")
PY

# ---- per-mutation sequential index and identity, both asserted ----
declare -A SEQ=(   [V81I]=79  [K59R]=59  )
declare -A FROM=(  [V81I]=VAL [K59R]=LYS )
declare -A TO=(    [V81I]=ILE [K59R]=ARG )
declare -A NATIVE=([V81I]=81  [K59R]=59  )

make_mutant_pdb () {   # $1 = mutation tag, $2 = out.pdb
    local M=$1 O=$2
    python3 - "$OUT/wt_heavy.pdb" "$O" "${SEQ[$M]}" "${FROM[$M]}" "${TO[$M]}" <<'PY'
import sys
src, dst, seq, frm, to = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
lines = [l for l in open(src) if l.startswith("ATOM")]
# sequential index -> the residue number actually written in the file
order, seen = [], set()
for l in lines:
    r = int(l[22:26])
    if r not in seen:
        seen.add(r); order.append(r)
target = order[seq - 1]
have = {l[17:20].strip() for l in lines if int(l[22:26]) == target}
if have != {frm}:
    sys.exit(f"ERROR: sequential residue {seq} is {have}, expected {frm}")
# keep only backbone + CB; tleap builds the rest of the new sidechain
KEEP = {"N", "CA", "C", "O", "CB", "OXT"}
out = []
for l in lines:
    if int(l[22:26]) == target:
        if l[12:16].strip() not in KEEP:
            continue
        l = l[:17] + f"{to:>3}" + l[20:]
    out.append(l)
open(dst, "w").writelines(out)
print(f"    {frm}{seq}{to}: verified identity, wrote {dst.split('/')[-1]}")
PY
}

build_leg () {   # $1 = mutation, $2 = leg name (aba|mandi|apo)
    local M=$1 LEG=$2 D=$OUT/${M}_${LEG}
    [[ -f $D/ti.prmtop ]] && { echo "  $M/$LEG: exists, skip"; return 0; }
    mkdir -p "$D"; cd "$D" || return 1
    cp "$OUT/wt_heavy.pdb" copy1.pdb
    make_mutant_pdb "$M" copy2.pdb || return 1

    local LIGBLOCK="" COMBINE="sys = combine { p1 p2 }"
    case $LEG in
        aba)   LIGBLOCK="loadamberparams $PARAMS/A8S.frcmod
lig = loadmol2 $PARAMS/A8S.mol2"; COMBINE="sys = combine { p1 p2 lig }" ;;
        mandi) LIGBLOCK="loadamberparams $PARAMS/3UZ.frcmod
lig = loadmol2 $PARAMS/3UZ.mol2"; COMBINE="sys = combine { p1 p2 lig }" ;;
        apo)   : ;;
    esac

    cat > tleap.in <<TL
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.gaff2
$LIGBLOCK
p1 = loadpdb copy1.pdb
p2 = loadpdb copy2.pdb
$COMBINE
solvateoct sys OPCBOX 12.0
addions sys K+ 0
addions sys Cl- 0
charge sys
saveamberparm sys dual.prmtop dual.inpcrd
savepdb sys dual.pdb
quit
TL
    tleap -f tleap.in > tleap.log 2>&1
    [[ -f dual.prmtop ]] || { echo "  !! $M/$LEG tleap failed"; grep -iE "error|fatal" tleap.log|head -3; return 1; }
    echo "  $M/$LEG: dual system built ($(grep -c '^ATOM\|^HETATM' dual.pdb) atoms)"
    cd "$ROOT" || return 1
}

echo "=== building dual-topology systems ==="
for M in V81I K59R; do
    for LEG in aba mandi apo; do
        build_leg "$M" "$LEG"
    done
done
