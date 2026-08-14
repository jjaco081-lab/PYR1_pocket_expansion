#!/bin/bash
# lib_solvate.sh -- THE canonical solvate + neutralise routine for this project.
#
# WHY THIS EXISTS
#   Every MD system here is built to be compared with the others, so the build
#   protocol has to be identical across all of them. It was not: S5_804_1 and the
#   first cut of S9_apo_closed neutralised BEFORE solvating, while S1-S4 solvate
#   first. That places counter-ions in vacuum against the solute instead of
#   substituting them into water -- a different starting state, discovered only by
#   diffing the tleap.in files after the fact.
#
#   So the settings live here, once. Source this file; do not retype tleap blocks.
#
# CANONICAL SETTINGS (README 19b)
#   protein   ff19SB
#   water     OPC          (the model ff19SB was fitted against)
#   ligand    GAFF2 + AM1-BCC
#   box       truncated octahedron, 12.0 A buffer
#   salt      0.15 M NaCl, ions-per-water = conc / 55.5
#   ORDER     solvate -> neutralise -> bulk salt      <- NOT interchangeable
#
#   leaprc.water.opc already loads frcmod.ionslm_126_opc, which covers Mn2+.
#   Adding frcmod.ions234lm_126_opc makes tleap exit "Could not open file": that is
#   TIP3P-era naming and no such file ships for OPC.
#
# OVERRIDING
#   Every constant reads from the environment, so a deliberate change is possible:
#       CONC=0.05 build_system ...
#   but the function PRINTS A WARNING for any value that differs from canonical, so
#   an accidental divergence cannot pass silently the way the ion ordering did.
#
# USAGE
#   source "$(dirname "$0")/lib_solvate.sh"
#   build_system <outdir> <protein.pdb> [ligand.mol2] [ligand.frcmod]

FF_PROTEIN=${FF_PROTEIN:-leaprc.protein.ff19SB}
FF_WATER=${FF_WATER:-leaprc.water.opc}
FF_LIGAND=${FF_LIGAND:-leaprc.gaff2}
BOX_CMD=${BOX_CMD:-solvateoct}
BOX_TYPE=${BOX_TYPE:-OPCBOX}
BUFFER=${BUFFER:-12.0}
CONC=${CONC:-0.15}

_canon_check () {
    local n=0
    [[ "$FF_PROTEIN" == "leaprc.protein.ff19SB" ]] || { echo "  !! FF_PROTEIN overridden: $FF_PROTEIN"; n=1; }
    [[ "$FF_WATER"   == "leaprc.water.opc"      ]] || { echo "  !! FF_WATER overridden: $FF_WATER"; n=1; }
    [[ "$FF_LIGAND"  == "leaprc.gaff2"          ]] || { echo "  !! FF_LIGAND overridden: $FF_LIGAND"; n=1; }
    [[ "$BOX_CMD"    == "solvateoct"            ]] || { echo "  !! BOX_CMD overridden: $BOX_CMD"; n=1; }
    [[ "$BOX_TYPE"   == "OPCBOX"                ]] || { echo "  !! BOX_TYPE overridden: $BOX_TYPE"; n=1; }
    [[ "$BUFFER"     == "12.0"                  ]] || { echo "  !! BUFFER overridden: $BUFFER"; n=1; }
    [[ "$CONC"       == "0.15"                  ]] || { echo "  !! CONC overridden: $CONC"; n=1; }
    [[ $n -eq 1 ]] && echo "  !! NON-CANONICAL BUILD -- this system is not directly comparable"
    return 0
}

# build_system <outdir> <protein.pdb> [ligand.mol2] [ligand.frcmod]
build_system () {
    local D=$1 PROT=$2 LIGMOL2=${3:-} FRC=${4:-}
    _canon_check
    mkdir -p "$D" || return 1
    cp "$PROT" "$D/protein.pdb" || return 1
    ( cd "$D" || exit 1

    {
      echo "source $FF_PROTEIN"
      echo "source $FF_WATER"
      # gaff2 is sourced even with no ligand so the loaded force-field stack is
      # identical in every system (S1/S3 have no ligand and source it too)
      echo "source $FF_LIGAND"
      if [[ -n "$LIGMOL2" ]]; then
          echo "loadamberparams $FRC"
          echo "lig = loadmol2 $LIGMOL2"
          echo "prot = loadpdb protein.pdb"
          echo "sys = combine { prot lig }"
      else
          echo "prot = loadpdb protein.pdb"
          echo "sys = combine { prot }"
      fi
      # ORDER IS LOAD-BEARING -- see the header
      echo "$BOX_CMD sys $BOX_TYPE $BUFFER"
      echo "addions sys Na+ 0"
      echo "addions sys Cl- 0"
      echo "addionsrand sys Na+ {nna} Cl- {ncl}"
      echo "charge sys"
      echo "check sys"
      echo "saveamberparm sys system.prmtop system.inpcrd"
      echo "savepdb sys system_solvated.pdb"
      echo "quit"
    } > tleap.in.template

    # pass 1: solvate with no bulk salt, purely to count waters
    sed -e "s/{nna}/0/" -e "s/{ncl}/0/" tleap.in.template > tleap.in
    tleap -f tleap.in > tleap_pass1.log 2>&1 || true
    local NWAT
    NWAT=$(grep -oP 'Added \K[0-9]+(?= residues)' tleap_pass1.log | tail -1 || true)
    if [[ -z "${NWAT:-}" ]]; then
        echo "  ERROR: no water count in $D/tleap_pass1.log"; tail -20 tleap_pass1.log; exit 1
    fi
    local NION
    NION=$(python3 -c "print(int(round($NWAT * $CONC / 55.5)))")

    # pass 2: rebuild at the target ionic strength
    sed -e "s/{nna}/$NION/" -e "s/{ncl}/$NION/" tleap.in.template > tleap.in
    tleap -f tleap.in > tleap_pass2.log 2>&1 || true

    # tleap exits non-zero on warnings, so verify by artefact and by charge
    [[ -f system.prmtop && -f system.inpcrd ]] || {
        echo "  ERROR: tleap wrote no prmtop/inpcrd in $D"
        grep -iE "error|fatal" tleap_pass2.log | head; exit 1; }
    local CHG NATOM
    CHG=$(grep -oP 'Total unperturbed charge:\s*\K-?[0-9.]+' tleap_pass2.log | tail -1)
    NATOM=$(grep -c -E '^(ATOM|HETATM)' system_solvated.pdb)
    echo "  $(basename "$D"): $NATOM atoms, $NWAT waters, ${NION}x NaCl (${CONC} M), net charge $CHG"
    python3 -c "
import sys
q=float('$CHG')
sys.exit('  ERROR: net charge %s is not neutral; PME would be wrong' % q) if abs(q)>0.01 else None" || exit 1
    if grep -qiE "Could not find|no parameters|Unknown residue|ATOM NOT FOUND" tleap_pass2.log; then
        echo "  ERROR: unparameterised atoms in $D"
        grep -iE "Could not find|no parameters|Unknown residue" tleap_pass2.log | head
        exit 1
    fi
    # record what was actually used, next to the system it produced
    {
      echo "protein_ff=$FF_PROTEIN"; echo "water_ff=$FF_WATER"; echo "ligand_ff=$FF_LIGAND"
      echo "box=$BOX_CMD/$BOX_TYPE buffer=$BUFFER"; echo "conc_M=$CONC"
      echo "waters=$NWAT"; echo "na=$NION"; echo "cl=$NION"; echo "atoms=$NATOM"
      echo "net_charge=$CHG"; echo "order=solvate,neutralise,bulk_salt"
      echo "built=$(date -Is)"; echo "builder=lib_solvate.sh"
    } > build_settings.txt
    ) || return 1
    return 0
}
