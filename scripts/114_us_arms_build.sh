#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH -J usbuild
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/us_build.log
#
# 114_us_arms_build.sh -- three more umbrella arms, turning the calibration into a
# full receptor x ligand cross-over.
#
#                    WT PYR1              PYR1^MANDI (K59R/V81I/F108A/F159L)
#   apo         negative (open)                    --
#   + ABA       POSITIVE (closed)            negative
#   + mandi        negative              POSITIVE (closed)
#
# The two WT arms are already running (112). These three complete the square.
# The method has to get all four right AND the diagonals have to FLIP -- that
# tests SELECTIVITY, not merely whether a gate can close, which is a far stronger
# calibration than the two-arm version.
#
# STARTING GEOMETRY, and the one that carries risk
#   quad + mandi   the cognate pair. 4WVO's own mutant sidechain rotamers are
#                  transplanted (101), so this is as close to a crystal complex as
#                  the 191-residue frame allows.
#   quad + ABA     ABA in a pocket evolved away from it.
#   WT + mandi     mandipropamid in a pocket that cannot hold it -- F108 sits
#                  0.62 A from a ligand heavy atom (§44c), and this is exactly the
#                  configuration whose minimisation was measured to be
#                  inescapable at low lambda in the TI (§46c: E = 4.5e8,
#                  |F|max = 8.4e6, flat over 10,000 steps).
#
#                  That is the physically correct answer -- WT genuinely cannot
#                  accommodate mandipropamid, which is why the quadruple exists --
#                  but it must be RELAXED rather than minimised in place. Here the
#                  ligand and side chains are free while only the backbone is
#                  restrained, which is the fix that worked in 102's ladder.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22
PY="env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/python"
PARAMS=$ROOT/data/mmgbsa/params
OUT=$ROOT/data/umbrella_arms
mkdir -p "$OUT"

echo "=== 1. closed 191-residue frame + ligand poses  $(date -Is)"
$PY "$ROOT/scripts/115_us_arm_prep.py" || exit 1

echo "=== 2. tleap build, 0.15 M KCl (canonical)  $(date -Is)"
build () {   # $1 arm  $2 protein.pdb  $3 ligand mol2 ("" = apo)  $4 frcmod
    local A=$1 P=$2 L=$3 F=$4
    local D=$OUT/$A          # separate statement: $A is not visible inside its own `local`
                             # (same bug as 91's `local D=$Q/$LEG`, belief row 55)
    [[ -s $D/system.prmtop ]] && { echo "    $A: exists, skip"; return 0; }
    mkdir -p "$D"; cd "$D" || return 1
    local LIGB="" COMB="sys = combine { prot }"
    if [[ -n "$L" ]]; then
        LIGB="loadamberparams $F
lig = loadmol2 $L"
        COMB="sys = combine { prot lig }"
    fi
    cat > tleap.in.template <<TL
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.gaff2
$LIGB
prot = loadpdb $P
$COMB
solvateoct sys OPCBOX 12.0
addions sys K+ 0
addions sys Cl- 0
addionsrand sys K+ {n} Cl- {n}
charge sys
saveamberparm sys system.prmtop system.inpcrd
savepdb sys system_solvated.pdb
quit
TL
    sed 's/{n}/0/g' tleap.in.template > tleap.in
    tleap -f tleap.in > tleap_pass1.log 2>&1 || true
    local NW NI
    NW=$(grep -oP 'Added \K[0-9]+(?= residues)' tleap_pass1.log | tail -1)
    [[ -n "${NW:-}" ]] || { echo "    !! $A: no water count"; tail -5 tleap_pass1.log; return 1; }
    NI=$(python3 -c "print(int(round($NW*0.15/55.5)))")
    sed "s/{n}/$NI/g" tleap.in.template > tleap.in
    rm -f system.prmtop system.inpcrd
    tleap -f tleap.in > tleap_pass2.log 2>&1 || true
    # -s not -f: tleap leaves a 0-BYTE prmtop behind when it fails, which -f accepts
    [[ -s system.prmtop ]] || { echo "    !! $A tleap failed"; grep -iE "error|fatal" tleap_pass2.log|head -3; return 1; }
    local CHG; CHG=$(grep -oP 'Total unperturbed charge:\s*\K-?[0-9.]+' tleap_pass2.log | tail -1)
    echo "    $A: $(grep -c -E '^(ATOM|HETATM)' system_solvated.pdb) atoms, $NW waters, ${NI}x KCl, charge $CHG"
    python3 -c "import sys; sys.exit('    !! $A net charge $CHG') if abs(float('$CHG'))>0.01 else None" || return 1
    cd "$ROOT" || return 1
}
build quad_mandi "$OUT/quad_closed.pdb" "$OUT/mandi_in_frame.mol2" "$PARAMS/3UZ.frcmod" || exit 1
build quad_aba   "$OUT/quad_closed.pdb" "$OUT/aba_in_frame.mol2"   "$PARAMS/A8S.frcmod" || exit 1
build wt_mandi   "$OUT/wt_closed.pdb"   "$OUT/mandi_in_frame.mol2" "$PARAMS/3UZ.frcmod" || exit 1

echo "=== 3. verify each arm  $(date -Is)"
$PY - <<'PYX'
import parmed as pmd, numpy as np, os
ROOT="/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion"
EXP={"quad_mandi":("3UZ-like",51,"ARG"),"quad_aba":("A8S",38,"ARG"),"wt_mandi":("3UZ-like",51,"LYS")}
for a,(lg,n,r59) in EXP.items():
    p=f"{ROOT}/data/umbrella_arms/{a}/system.prmtop"
    if not os.path.exists(p): print(f"  {a}: MISSING"); continue
    t=pmd.load_file(p, f"{ROOT}/data/umbrella_arms/{a}/system.inpcrd")
    L=[x for x in t.residues if x.name in ("LIG","A8S","3UZ")]
    na=sum(len(x.atoms) for x in L)
    got59=t.residues[58].name
    prot=[x for x in t.residues if x.name not in ("WAT","HOH","K+","Cl-","LIG","A8S","3UZ")]
    # restraint atoms, computed the same way 112 does
    iat=[[y for y in t.residues[k-1].atoms if y.name=="CA"][0].idx+1 for k in (88,116)]
    x=np.array(t.coordinates)
    lh=[y.idx for y in t.atoms if y.residue.name in ("LIG","A8S","3UZ") and y.atomic_number>1]
    ph=[y.idx for y in t.atoms if y.residue.name not in ("WAT","HOH","K+","Cl-","LIG","A8S","3UZ") and y.atomic_number>1]
    D=np.linalg.norm(x[ph][:,None,:]-x[lh][None,:,:],axis=2)
    rc=np.linalg.norm(x[iat[0]-1]-x[iat[1]-1])
    ok="OK" if (na==n and got59==r59) else "MISMATCH"
    print(f"  {a:<11} {len(prot):3d} res  res59={got59}({r59})  LIG {na}({n})  {ok}  "
          f"P88-R116 {rc:5.2f} A  min prot-lig {D.min():5.2f} A  <2A {int((D<2.0).sum())}")
    print(f"              restraint atoms P88 CA={iat[0]} R116 CA={iat[1]}")
PYX
echo "=== build complete $(date -Is)"
