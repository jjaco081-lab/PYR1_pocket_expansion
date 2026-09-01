#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH -J pcbuild
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/posecheck_build.log
#
# 179 -- build six ligand/receptor systems to ask whether MD RETAINS the pose it
# was started from.
#
# Jannis: most of the poses this project now trains on were predicted ONCE, and a
# co-folding model will place a ligand in the pocket whether or not it belongs
# there (README 93a: every one of 362 ligands landed on the ABA site, hits and
# non-hits alike). Short MD is the standard check -- it was used this way in the
# nitazene sensor work. A pose that is forced should shed it quickly; a pose that
# is right should hold.
#
# THE DESIGN IS A 3x3: three CRYSTAL poses that are known-correct as positive
# controls, and three PREDICTED poses. Within it sit two paired comparisons that
# do the real work:
#   mandipropamid crystal vs mandipropamid PREDICTED -- same ligand, same
#     receptor, only the pose source differs
#   fludioxonil (seed spread 0.71 A) vs benoxacor (3.27 A) -- does seed
#     agreement, which README 93c showed understates true error ~3x, predict
#     which poses survive?
#
# Ligand parameterisation follows 64 exactly: antechamber GAFF2/AM1-BCC with
# maxcyc=0, then parmchk2, and ATTN lines in the frcmod are reported not ignored.
# Solvation is lib_solvate.sh so the protocol cannot drift from the rest of the
# project (ff19SB/OPC, 0.15 M KCl, solvate -> neutralise -> salt).
set -uo pipefail
R=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd $R || exit 1
source $R/scripts/lib_solvate.sh
module load amber/22 2>/dev/null
PC=$R/data/posecheck
mkdir -p $PC/params

#   tag                     ligand pose PDB                              SMILES                                              receptor frame
#   tag | ligand .mol (bond orders already correct, at the pose) | net charge | receptor frame | resname
#
# The .mol files were built by 135/165 with AssignBondOrdersFromTemplate and
# verified to sit on the pose to within 0.01 A. Rebuilding them here was tried
# and failed six different ways -- explicit hydrogens defeat the template match,
# sanitize=False defeats the substructure match, and a half-built SDF makes
# antechamber die on the mandipropamid alkyne. Reuse the validated artefact.
SYS=(
"aba_crystal|$R/data/stage1/params/A8S_anion.mol|-1|$R/data/stage1/wt_aba.pdb|ABX"
"mandi_crystal|$R/data/stage1/params/3UZ.mol|0|$R/data/stage1/wt_mandi.pdb|MDX"
"mandi_predicted|$R/data/agro_params/MDP.mol|0|$R/data/agro_params/frame_MDP.pdb|MDY"
"win_crystal|$R/results/win_crossover/WI5.mol|0|$R/results/win_crossover/sensor.pdb|WNX"
"fludioxonil_pred|$R/data/agro_params/FLD.mol|0|$R/data/agro_params/frame_FLD.pdb|FLX"
"benoxacor_pred|$R/data/agro_params/BNX.mol|0|$R/data/agro_params/frame_BNX.pdb|BNY"
)
fail=0
for row in "${SYS[@]}"; do
    IFS='|' read -r TAG LIGMOL NC FRAME RN <<< "$row"
    W=$PC/params/$TAG; mkdir -p $W
    if [[ ! -s $W/$RN.mol2 ]]; then
        # 135/165 wrote these with kekulize=False, so aromatic bonds are MDL
        # type 4 -- antechamber rejects that with "Check Weird Bonds / Fatal
        # Error". Rewrite KEKULIZED (orders 1/2 only) before parameterising.
        env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/dockenv/bin/python - <<PY
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
m = Chem.MolFromMolFile("$LIGMOL", removeHs=False, sanitize=True)
assert m is not None, "cannot read $LIGMOL"
Chem.Kekulize(m, clearAromaticFlags=True)
Chem.MolToMolFile(m, "$W/lig.mol", kekulize=True)
n4 = sum(1 for b in m.GetBonds() if b.GetBondTypeAsDouble() == 1.5)
print("$TAG kekulised:", m.GetNumAtoms(), "atoms,", n4, "aromatic bonds left")
PY
        ( cd $W && antechamber -i lig.mol -fi mdl -o $RN.mol2 -fo mol2 \
            -c bcc -nc $NC -rn $RN -pf y -ek "maxcyc=0" > ac.log 2>&1 )
    fi
    if [[ ! -s $W/$RN.mol2 ]]; then echo "!! $TAG: antechamber failed"; tail -12 $W/ac.log; fail=1; continue; fi
    [[ -s $W/$RN.frcmod ]] || ( cd $W && parmchk2 -i $RN.mol2 -f mol2 -o $RN.frcmod -s gaff2 )
    if grep -q ATTN $W/$RN.frcmod 2>/dev/null; then
        echo "   $TAG: parmchk2 GUESSED parameters (ATTN) -- $(grep -c ATTN $W/$RN.frcmod) lines"
    fi
    # protein-only pdb for tleap
    grep '^ATOM' $FRAME > $W/protein.pdb
    build_system "$PC/$TAG" "$W/protein.pdb" "$W/$RN.mol2" "$W/$RN.frcmod" \
        > $W/build.log 2>&1 || { echo "!! $TAG: build_system failed"; tail -12 $W/build.log; fail=1; continue; }
    N=$(grep -c . $PC/$TAG/system.prmtop 2>/dev/null || echo 0)
    echo "=== $TAG built  ($(grep -m1 'atoms' $W/build.log 2>/dev/null || echo 'see build.log'))"
done
echo "build finished, fail=$fail"
exit $fail
