#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 6:00:00
#SBATCH -J panelbuild
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/panel_build_%j.log
#
# 226 -- parameterise the four control ligands and build the 8 MM-GBSA systems
# (4 ligands x {WT, known hit}).
#
# Follows 179 exactly: antechamber GAFF2/AM1-BCC maxcyc=0, then parmchk2, with
# ATTN lines REPORTED not ignored; solvation via lib_solvate.sh so the protocol
# cannot drift (ff19SB/OPC, 0.15 M KCl, solvate -> neutralise -> salt).
#
# The .mol files were written by 223 from the SMILES that generated the
# conformer, kekulised (0 aromatic bonds left) with the atom order asserted
# element-by-element against 209's pose. MDL bond type 4 kills antechamber.
set -uo pipefail
R=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd $R || exit 1
source $R/scripts/lib_solvate.sh
module load amber/22 2>/dev/null
PM=$R/data/panel_mmgbsa
PY=/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python

mapfile -t ROWS < <($PY - <<'PYEOF'
import json
for f in json.load(open("data/panel_mmgbsa/frames.json")):
    d = f["mol"].rsplit("/",1)[0]
    print(f"{f['ligand']}|{f['resname']}|{f['mol']}|{f['wt']}|{f['hit']}|{d}")
PYEOF
)
echo "${#ROWS[@]} ligands to build"
fail=0
for row in "${ROWS[@]}"; do
    IFS='|' read -r LIG RN MOL WT HIT D <<< "$row"
    W=$D/params; mkdir -p $W
    if [[ ! -s $W/$RN.mol2 ]]; then
        ( cd $W && antechamber -i $MOL -fi mdl -o $RN.mol2 -fo mol2 \
            -c bcc -nc 0 -rn $RN -pf y -ek "maxcyc=0" > ac.log 2>&1 )
    fi
    if [[ ! -s $W/$RN.mol2 ]]; then
        echo "!! $LIG: antechamber FAILED"; tail -15 $W/ac.log; fail=1; continue; fi
    [[ -s $W/$RN.frcmod ]] || ( cd $W && parmchk2 -i $RN.mol2 -f mol2 -o $RN.frcmod -s gaff2 )
    if grep -q ATTN $W/$RN.frcmod 2>/dev/null; then
        echo "   $LIG: parmchk2 GUESSED $(grep -c ATTN $W/$RN.frcmod) parameters (ATTN)"
    fi
    for ARM in WT HIT; do
        FRAME=$WT; [[ $ARM == HIT ]] && FRAME=$HIT
        grep '^ATOM' $FRAME > $W/protein_$ARM.pdb
        OUTD=$PM/sys/${LIG//-/}_$ARM
        build_system "$OUTD" "$W/protein_$ARM.pdb" "$W/$RN.mol2" "$W/$RN.frcmod" \
            > $W/build_$ARM.log 2>&1 || { echo "!! $LIG/$ARM build_system FAILED";
                tail -15 $W/build_$ARM.log; fail=1; continue; }
        # assert on the RESULT, not the exit code
        if [[ ! -s $OUTD/system.prmtop ]]; then
            echo "!! $LIG/$ARM: no prmtop"; fail=1; continue; fi
        NA=$(grep -A2 "FLAG POINTERS" $OUTD/system.prmtop | tail -1 | awk '{print $1}')
        echo "=== $LIG/$ARM built: $NA atoms -> $OUTD"
    done
done
N=$(ls -d $PM/sys/*/ 2>/dev/null | wc -l)
echo "systems built: $N (expect 8), fail=$fail"
test "$N" -eq 8 || exit 1
exit $fail
