#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH -J quadbuild
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/quad_build.log
#
# 91_quad_build.sh -- TI systems for the FULL PYR1^MANDI quadruple.
#
# WHY THE QUADRUPLE AND NOT MORE SINGLES
#   README 44b: the 43d test asked whether V81I and K59R each shift selectivity
#   toward mandipropamid ALONE. Nothing in the ground truth requires that --
#   PYR1^MANDI = K59R/V81I/F108A/F159L was selected as a SET. The calibration
#   13e actually specifies, and which TI has never been given, is the quadruple:
#   PYR1^MANDI + mandipropamid must rank above WT + mandipropamid.
#
# THE TWO ARMS -- one variable, deliberately
#   mandi_xtal : mandipropamid in its 4WVO crystallographic pose
#   mandi_dock : mandipropamid placed by docking into a MODEL of the quad pocket,
#                using no 4WVO complex information at all
#   Both arms share the SAME ABA leg and the SAME tleap-built mutant sidechains,
#   so the ONLY thing that differs between them is the ligand pose. Giving the
#   crystal arm crystal rotamers too would confound pose provenance with rotamer
#   quality and neither arm would answer the question.
#
#   This matters for production: novel ligands will never have a complex
#   structure, so the docked arm is the realistic accuracy, and the crystal arm
#   is the ceiling. The gap between them is the cost of not having a structure.
#
# SALT: 0.15 M KCl, two-pass, per lib_solvate.sh. The pilot was NEUTRALISED ONLY
#   (README 44d) -- corrected here.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22
command -v tleap >/dev/null || { echo "no tleap"; exit 1; }
# `module load amber/22` puts amber's own ParmEd 3.4.1 on PYTHONPATH, where it
# shadows the conda ParmEd 4.3 and dies on `parmed.utils.six.moves`. Strip
# PYTHONPATH for every interpreter call rather than relying on import order.
PY="env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/python"
RDPY="env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/esmfold2/bin/python"  # RDKit only lives here
SMINA=/bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/smina
OBABEL=/bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/obabel
PARAMS=$ROOT/data/mmgbsa/params
Q=$ROOT/data/ti_quad
mkdir -p "$Q"
SRC=$ROOT/results/stage1_rosetta/_input_wt_aba_A8S_anion.pdb

echo "=== 1. residue mapping, asserted as a SET ==="
$PY - "$SRC" "$Q" <<'PY' || exit 1
import sys, json
src, out = sys.argv[1], sys.argv[2]
MUT = {59: ("LYS", "ARG"), 81: ("VAL", "ILE"), 108: ("PHE", "ALA"), 159: ("PHE", "LEU")}
lines = [l for l in open(src) if l.startswith("ATOM")]
order, seen, name = [], set(), {}
for l in lines:
    r = int(l[22:26]); name[r] = l[17:20].strip()
    if r not in seen: seen.add(r); order.append(r)
seqmap = {}
for nat, (frm, to) in MUT.items():
    assert name[nat] == frm, f"native {nat} is {name[nat]}, expected {frm}"
    seqmap[nat] = order.index(nat) + 1
    print(f"    native {nat:3d} {frm}->{to}  sequential {seqmap[nat]:3d}")
# the two-valine trap: sequential 81 is native 83, ALSO a VAL. Asserting only
# "is it a VAL" cannot tell them apart, so pin the whole mapping instead.
assert seqmap == {59: 59, 81: 79, 108: 106, 159: 157}, f"mapping changed: {seqmap}"
assert name[order[80]] == "VAL" and order[80] == 83, "two-valine trap check failed"
print("    two-valine trap: sequential 81 = native 83 VAL, correctly NOT a target")
# copy1: WT, heavy atoms only (tleap rebuilds H identically in both copies)
keep = [l for l in lines if (l[76:78].strip() or l[12:16].strip()[0]) != "H"]
open(f"{out}/wt_heavy.pdb", "w").writelines(keep)
# copy2: same backbone, 4 residues reduced to backbone+CB and renamed
KEEP = {"N", "CA", "C", "O", "CB", "OXT"}
o = []
for l in keep:
    r = int(l[22:26])
    if r in MUT:
        if l[12:16].strip() not in KEEP: continue
        l = l[:17] + f"{MUT[r][1]:>3}" + l[20:]
    o.append(l)
open(f"{out}/quad_bbcb.pdb", "w").writelines(o)
json.dump({str(k): v for k, v in seqmap.items()}, open(f"{out}/seqmap.json", "w"))
print(f"    wt_heavy.pdb {len(keep)} atoms | quad_bbcb.pdb {len(o)} atoms")
PY

echo "=== 2. quad receptor model for docking (NO 4WVO coordinates) ==="
cd "$Q" || exit 1
if [[ ! -f quad_model.pdb ]]; then
cat > leap_quad.in <<'EOF'
source leaprc.protein.ff19SB
p = loadpdb quad_bbcb.pdb
saveamberparm p quad_vac.prmtop quad_vac.inpcrd
quit
EOF
tleap -f leap_quad.in > leap_quad.log 2>&1
[[ -f quad_vac.prmtop ]] || { echo "  !! tleap failed"; grep -iE "error|fatal" leap_quad.log|head -3; exit 1; }
cat > min_quad.in <<'EOF'
relax the built sidechains, backbone restrained
 &cntrl
  imin=1, maxcyc=3000, ncyc=1500, ntmin=1,
  ntb=0, igb=1, cut=16.0, ntpr=500,
  ntr=1, restraintmask='@N,CA,C,O', restraint_wt=10.0,
 /
EOF
sander -O -i min_quad.in -p quad_vac.prmtop -c quad_vac.inpcrd -ref quad_vac.inpcrd \
       -o min_quad.out -r quad_min.rst7 || { echo "  !! sander failed"; exit 1; }
grep -q "FINAL RESULTS" min_quad.out || { echo "  !! minimisation did not finish"; exit 1; }
ambpdb -p quad_vac.prmtop -c quad_min.rst7 > quad_model.pdb 2>/dev/null
fi
echo "    quad_model.pdb: $(grep -c '^ATOM' quad_model.pdb) atoms"

echo "=== 3. dock mandipropamid into the MODEL (no complex structure used) ==="
# antechamber's mol2 carries GAFF2 atom types. openbabel cannot translate them
# and smina types the whole molecule as Du (dummy) / Ca (calcium) -- docking then
# scores a molecule made of dummy atoms and returns a plausible number for it.
# Retype to Sybyl first; the BOND block already carries real bond orders, so
# atom order, names, coordinates and charges are untouched.
$PY "$ROOT/scripts/92_lig_sybyl.py" "$PARAMS/3UZ.mol2" lig_sybyl.mol2 || exit 1
# fresh 3D conformer so the search does not start from the crystal answer
[[ -f lig_rand.mol2 ]] || $OBABEL lig_sybyl.mol2 -O lig_rand.mol2 --gen3d >/dev/null 2>&1
$OBABEL lig_rand.mol2 -oreport 2>/dev/null | grep -i "^FORMULA" | sed 's/^/    /'
# box centre from POCKET-LINING PROTEIN ATOMS ONLY -- placing the box on the
# ligand would smuggle the answer back into a test about not having one
read CX CY CZ <<< "$($PY "$ROOT/scripts/95_dock_box.py" quad_model.pdb)"
echo "    box centre $CX $CY $CZ (from protein only)"
if [[ ! -f docked.mol2 ]]; then
$SMINA --receptor quad_model.pdb --ligand lig_rand.mol2 \
       --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
       --size_x 22 --size_y 22 --size_z 22 \
       --exhaustiveness 32 --num_modes 9 --seed 42 \
       --out docked.mol2 > dock.log 2>&1 || { echo "  !! smina failed"; tail -5 dock.log; exit 1; }
fi
grep -E "^[0-9]+ +-[0-9]" dock.log | head -5 | sed 's/^/    pose /'
# smina both REORDERS atoms and discards nonpolar hydrogens, so the pose has to
# be mapped back by graph matching, not by atom order (see 97's header).
$RDPY "$ROOT/scripts/97_pose_map.py" lig_sybyl.mol2 docked.mol2 \
      "$PARAMS/3UZ.mol2" mandi_dock_noH.mol2 || exit 1

# smina discards nonpolar hydrogens, so the 21 it dropped are still sitting at
# their CRYSTAL positions on a molecule whose heavy atoms have moved. Relax them
# with every heavy atom restrained: the ligand alone, its own gaff2 parameters,
# so the docked heavy-atom pose is preserved exactly and only H moves.
if [[ ! -f mandi_dock.mol2 ]]; then
cat > leap_ligh.in <<'TL'
source leaprc.gaff2
loadamberparams PARAMSDIR/3UZ.frcmod
lig = loadmol2 mandi_dock_noH.mol2
saveamberparm lig ligH.prmtop ligH.inpcrd
quit
TL
sed -i "s|PARAMSDIR|$PARAMS|" leap_ligh.in
tleap -f leap_ligh.in > leap_ligh.log 2>&1
[[ -f ligH.prmtop ]] || { echo "  !! ligand tleap failed"; grep -iE "error|fatal" leap_ligh.log|head -3; exit 1; }
cat > min_ligh.in <<'TL'
relax hydrogens only, heavy atoms restrained
 &cntrl
  imin=1, maxcyc=1500, ntmin=2, drms=0.05,
  ntb=0, igb=0, cut=99.0, ntpr=200,
  ntr=1, restraintmask='!@H=', restraint_wt=100.0,
 /
TL
# ntmin=2 (pure steepest descent): the conjugate-gradient default converges to
# RMS ~0.05 and then dies in LINMIN, which looks like a failure but is not.
sander -O -i min_ligh.in -p ligH.prmtop -c ligH.inpcrd -ref ligH.inpcrd \
       -o min_ligh.out -r ligH_min.rst7 >/dev/null 2>&1
# judge the result by the MEASURED gradient, not by whether a banner was printed
[[ -s ligH_min.rst7 ]] || { echo "  !! H relaxation wrote no restart"; exit 1; }
$PY - <<'CHK' || exit 1
import re, sys
rows = re.findall(r"^\s*\d+\s+(-?[\d.]+E[+-]\d+)\s+([\d.]+E[+-]\d+)",
                  open("min_ligh.out").read(), re.M)
if not rows:
    sys.exit("  !! no minimisation steps parsed from min_ligh.out")
g = float(rows[-1][1])
print(f"    H relaxation: {len(rows)} reported steps, final RMS gradient {g:.3f}")
sys.exit("  !! H relaxation did not converge" if g > 0.5 else None)
CHK
$PY "$ROOT/scripts/96_mol2_setcoords.py" mandi_dock_noH.mol2 ligH_min.rst7 mandi_dock.mol2 || exit 1
fi
# the stereocentre could not be enforced during the graph match; now that the
# hydrogens are back it is perceptible again, so check it rather than assume it
$PY "$ROOT/scripts/92_lig_sybyl.py" mandi_dock.mol2 mandi_dock_sybyl.mol2 >/dev/null || exit 1
$RDPY "$ROOT/scripts/98_check_stereo.py" lig_sybyl.mol2 mandi_dock_sybyl.mol2 || exit 1



echo "=== 3e. crystal rotamers for the crystal arm ==="
# The crystal arm takes BOTH the ligand pose and the four mutant sidechains from
# 4WVO. Giving it tleap-built rotamers instead left it clashed at lambda=1 as
# well as lambda=0 (0.64 A), i.e. worse than the docked arm it is meant to bound.
# See 101's header.
[[ -f quad_xtal.pdb ]] || $PY "$ROOT/scripts/101_xtal_rotamers.py" \
    wt_heavy.pdb "$ROOT/data/4WVO.cif" quad_xtal.pdb || exit 1

echo "=== 4. build the four dual-topology legs (0.15 M KCl, two-pass) ==="
build_leg () {   # $1 = leg, $2 = ligand mol2 ("" for apo), $3 = frcmod, $4 = copy2 pdb
    local LEG=$1 MOL2=$2 FRC=$3 CP2=${4:-quad_bbcb.pdb}
    local D=$Q/$LEG          # separate statement: $LEG is not visible inside its own `local`
    [[ -f $D/ti.prmtop ]] && { echo "    $LEG: already merged, skip"; return 0; }
    mkdir -p "$D"; cd "$D" || return 1
    cp "$Q/wt_heavy.pdb" copy1.pdb; cp "$Q/$CP2" copy2.pdb
    local LIGBLOCK="" COMBINE="sys = combine { p1 p2 }"
    if [[ -n "$MOL2" ]]; then
        LIGBLOCK="loadamberparams $FRC
lig = loadmol2 $MOL2"
        COMBINE="sys = combine { p1 p2 lig }"
    fi
    cat > tleap.in.template <<TL
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
addionsrand sys K+ {n} Cl- {n}
charge sys
saveamberparm sys dual.prmtop dual.inpcrd
savepdb sys dual.pdb
quit
TL
    sed 's/{n}/0/g' tleap.in.template > tleap.in
    tleap -f tleap.in > tleap_pass1.log 2>&1 || true
    local NWAT NION
    NWAT=$(grep -oP 'Added \K[0-9]+(?= residues)' tleap_pass1.log | tail -1)
    [[ -n "${NWAT:-}" ]] || { echo "    !! $LEG: no water count"; tail -5 tleap_pass1.log; return 1; }
    NION=$(python3 -c "print(int(round($NWAT*0.15/55.5)))")
    sed "s/{n}/$NION/g" tleap.in.template > tleap.in
    rm -f dual.prmtop dual.inpcrd
    tleap -f tleap.in > tleap_pass2.log 2>&1 || true
    [[ -f dual.prmtop ]] || { echo "    !! $LEG tleap failed"; grep -iE "error|fatal" tleap_pass2.log|head -3; return 1; }
    local CHG
    CHG=$(grep -oP 'Total unperturbed charge:\s*\K-?[0-9.]+' tleap_pass2.log | tail -1)
    echo "    $LEG: $(grep -c -E '^(ATOM|HETATM)' dual.pdb) atoms, $NWAT waters, ${NION}x KCl, net charge $CHG"
    python3 -c "import sys; sys.exit('    !! $LEG net charge $CHG' ) if abs(float('$CHG'))>0.01 else None" || return 1
    cd "$ROOT" || return 1
}
build_leg aba        "$PARAMS/A8S.mol2"  "$PARAMS/A8S.frcmod"  || exit 1
build_leg mandi_xtal "$PARAMS/3UZ.mol2"  "$PARAMS/3UZ.frcmod"  quad_xtal.pdb || exit 1
build_leg mandi_dock "$Q/mandi_dock.mol2" "$PARAMS/3UZ.frcmod" || exit 1
build_leg apo        ""                  ""                    || exit 1

echo "=== 5. tiMerge ==="
$PY "$ROOT/scripts/93_quad_merge.py" || exit 1
echo "=== BUILD COMPLETE $(date -Is) ==="
