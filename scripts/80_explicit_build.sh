#!/bin/bash
# 80_explicit_build.sh -- rebuild the two MM-GBSA REFERENCE systems in EXPLICIT water.
#
# WHY
#   README 40 traced the retracted "K59 flip" to one number: WT+ABA moved +7.75
#   kcal/mol between 50 and 250 ps and was still moving +2.62 inside the 250 ps
#   window. Measurement (README 41) showed the mechanism: ABA slides 2-4 A out of
#   its crystallographic pose while its K59 salt bridge stays intact 100% of the
#   time -- the ligand, the lysine side chain and the gate migrate together.
#
#   The cause is that the crystal pose is a minimum of REALITY, which includes
#   ordered bridging waters, not a minimum of the implicit-solvent Hamiltonian we
#   integrated. There were zero explicit waters in that build (verified: the input
#   PDB had none and leap did no solvation). GB replaces every water with a
#   continuum, and a continuum cannot donate a directional hydrogen bond. So the
#   one system positioned by real physics had the furthest to fall.
#
#   The energy penalty landed on EGB (+7.65) and not VDWAALS (+1.51) because ABA
#   carries a localised formal -1 (A8S.mol2 net charge -1.001) and the Born term
#   goes as q^2/R_eff -- exquisitely sensitive to burial depth. Mandipropamid is
#   neutral (3UZ.mol2 net charge -0.003), so its EGB moved only +1.60 and its pose
#   held at 1.7 A.
#
# WHAT THIS CHANGES, AND WHAT IT DELIBERATELY DOES NOT
#   Sampling moves to explicit solvent. SCORING stays MM-GBSA on water-stripped
#   frames, which is the standard protocol and keeps these numbers commensurable
#   with the 46 implicit runs. Exactly one variable changes.
#
# NO CRYSTALLOGRAPHIC WATERS -- a decision, not an omission (user, 2026-08-21)
#   Placing ABA's ordered waters while mandipropamid and any future test ligand get
#   only bulk solvent would hand the cognate ligand a structural advantage no test
#   ligand can have, and every ligand-swap difference would inherit it. That is the
#   same asymmetry that produced the artefact this script exists to fix. Uniform
#   bulk solvation for every ligand is the only version where f(test) - f(cognate)
#   stays interpretable.
#
# PROTEIN COORDINATES ARE IDENTICAL BETWEEN ARMS
#   Verified: both stage-1 inputs carry the same 1414 protein atoms. Only the
#   ligand differs. That is the ligand-swap design from README 23 and it is what
#   makes the receptor and the force field cancel.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1

# tleap/pdb4amber are AmberTools -- `module load amber/22`, NOT amber/22_mpi_cuda
module load amber/22
command -v tleap      >/dev/null || { echo "tleap not found after 'module load amber/22'"; exit 1; }
command -v pdb4amber  >/dev/null || { echo "pdb4amber not found"; exit 1; }

source scripts/lib_solvate.sh

OUT=$ROOT/data/mmgbsa_explicit
PARAMS=$ROOT/data/mmgbsa/params
mkdir -p "$OUT"

# protein-only PDB, taken from the same stage-1 input the implicit runs used so the
# starting coordinates are bit-identical and solvation is the only difference
PROT=$OUT/pyr1_wt.pdb
awk '/^ATOM|^TER/ && substr($0,18,3)!="LIG"' \
    results/stage1_rosetta/_input_wt_aba_A8S_anion.pdb > "$PROT"
NP=$(grep -c '^ATOM' "$PROT")
[[ $NP -eq 1414 ]] || { echo "ERROR: protein is $NP atoms, expected 1414"; exit 1; }
echo "protein: $NP atoms"

# assert the ligands are what their names claim BEFORE 6 GPU-days depend on it.
# README 31 lost a whole stage to params that silently never loaded, and README 33
# to two ligands sharing a name in one process.
python3 - <<'PY' || exit 1
import sys
exp = {"A8S": (-1.0, 38), "3UZ": (0.0, 51)}
for name, (q_exp, n_exp) in exp.items():
    f = f"/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/data/mmgbsa/params/{name}.mol2"
    L = open(f).read().split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0].strip().split("\n")
    q = sum(float(l.split()[-1]) for l in L)
    if len(L) != n_exp or abs(q - q_exp) > 0.05:
        sys.exit(f"ERROR: {name} has {len(L)} atoms / charge {q:+.3f}, "
                 f"expected {n_exp} / {q_exp:+.1f}")
    print(f"  {name}: {len(L)} atoms, net charge {q:+.3f}  OK")
PY

for arm in aba mandi; do
    case $arm in
        aba)   LIG=A8S ;;
        mandi) LIG=3UZ ;;
    esac
    for seed in 0 1 2; do
        D=$OUT/${arm}_WT_s${seed}
        if [[ -f $D/system.prmtop ]]; then echo "  $(basename $D): exists, skip"; continue; fi
        build_system "$D" "$PROT" "$PARAMS/${LIG}.mol2" "$PARAMS/${LIG}.frcmod" \
            || { echo "  BUILD FAILED: $D"; continue; }
        echo "$LIG" > "$D/ligand.txt"
    done
done

echo
echo "=== built ==="
for d in "$OUT"/*/; do
    [[ -f $d/build_settings.txt ]] || continue
    printf "  %-18s %s\n" "$(basename "$d")" \
        "$(grep -E 'atoms=|net_charge=|waters=' "$d/build_settings.txt" | tr '\n' ' ')"
done
