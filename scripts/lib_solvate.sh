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
#   salt      0.15 M KCl, ions-per-water = conc / 55.5
#   ORDER     solvate -> neutralise -> bulk salt      <- NOT interchangeable
#
# WHY KCl AND NOT NaCl
#   0.15 M NaCl is the reflex, but it describes mammalian EXTRACELLULAR fluid. Every
#   host relevant here has a K+-dominated cytosol: E. coli ~150-250 mM K+ (Na+
#   5-15 mM), S. cerevisiae ~200-300 mM K+, Arabidopsis cytosol ~100 mM K+ with Na+
#   actively excluded. PYR1 is a cytosolic plant protein assayed in yeast.
#   There is also a project-specific reason: Na+ binds carboxylates more strongly
#   than K+, and the central interaction in this work is the K59-ABA carboxylate
#   salt bridge. Excess Na+ affinity for that carboxylate is a plausible artefact
#   in precisely the place being measured.
#
# ION MODEL
#   Li-Merz 12-6 (loaded by leaprc.water.opc) is well validated for MONOVALENT ions
#   at this concentration, so it is the default. It is NOT adequate for divalent
#   metals -- a 12-6 model cannot match hydration free energy and ion-oxygen
#   distance simultaneously. Any system carrying Mn2+ (S4/S5, HAB1's catalytic
#   centre) should set ION_FRCMOD=frcmod.ionslm_1264_opc, and pmemd then also needs
#   lj1264=1 in &cntrl. Not enabled by default because the current systems have no
#   metal and it would change monovalent behaviour for no benefit.
#
# HISTIDINE PROTONATION
#   tleap silently makes every HIS an HIE. In this protein that is a live problem:
#   H115 is a LATCH residue AND sits 3.86 A from ABA, so its tautomer sets the
#   H-bonding at the exact interface these simulations measure; H60 (7.1 A) is a
#   dimerisation residue. `reduce` assigns tautomers by optimising the local
#   H-bond network, which is far better than a blanket default.
#   LIMITATION, stated plainly: reduce is geometry-based, NOT a pKa calculation. It
#   will not propose HIP (doubly protonated). propka/pdb2pqr are not installed here;
#   if a His turns out to matter quantitatively, compute its pKa properly.
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
CATION=${CATION:-K+}          # see header: cytosolic, not extracellular
ANION=${ANION:-Cl-}
ION_FRCMOD=${ION_FRCMOD:-}    # set to frcmod.ionslm_1264_opc for divalent metals
HMR=${HMR:-0}                 # OFF: 2 fs canonical, see lib_mdinputs.sh header
PROTONATE=${PROTONATE:-1}     # assign His tautomers with reduce instead of tleap default

_canon_check () {
    local n=0
    [[ "$FF_PROTEIN" == "leaprc.protein.ff19SB" ]] || { echo "  !! FF_PROTEIN overridden: $FF_PROTEIN"; n=1; }
    [[ "$FF_WATER"   == "leaprc.water.opc"      ]] || { echo "  !! FF_WATER overridden: $FF_WATER"; n=1; }
    [[ "$FF_LIGAND"  == "leaprc.gaff2"          ]] || { echo "  !! FF_LIGAND overridden: $FF_LIGAND"; n=1; }
    [[ "$BOX_CMD"    == "solvateoct"            ]] || { echo "  !! BOX_CMD overridden: $BOX_CMD"; n=1; }
    [[ "$BOX_TYPE"   == "OPCBOX"                ]] || { echo "  !! BOX_TYPE overridden: $BOX_TYPE"; n=1; }
    [[ "$BUFFER"     == "12.0"                  ]] || { echo "  !! BUFFER overridden: $BUFFER"; n=1; }
    [[ "$CONC"       == "0.15"                  ]] || { echo "  !! CONC overridden: $CONC"; n=1; }
    [[ "$CATION"     == "K+"                    ]] || { echo "  !! CATION overridden: $CATION"; n=1; }
    [[ "$ANION"      == "Cl-"                   ]] || { echo "  !! ANION overridden: $ANION"; n=1; }
    [[ $n -eq 1 ]] && echo "  !! NON-CANONICAL BUILD -- this system is not directly comparable"
    return 0
}

# build_system <outdir> <protein.pdb> [ligand.mol2] [ligand.frcmod]
build_system () {
    local D=$1 PROT=$2 LIGMOL2=${3:-} FRC=${4:-}
    _canon_check
    mkdir -p "$D" || return 1
    cp "$PROT" "$D/protein.pdb" || return 1
    if [[ "$PROTONATE" == "1" ]]; then
        # reduce -build adds H, flips Asn/Gln/His and picks His tautomers from the
        # H-bond network; -Quiet keeps the log readable. tleap then reads the HID/
        # HIE/HIP names reduce assigned rather than defaulting everything to HIE.
        # `pdb4amber --reduce` runs reduce AND renames HIS to HID/HIE/HIP from the
        # hydrogens it added. Do NOT add --nohyd: stripping the hydrogens destroys
        # exactly the information reduce just produced, pdb4amber writes plain HIS,
        # and tleap silently defaults everything back to HIE. Verified: with
        # --reduce --nohyd all six histidines come back HIS; with --reduce alone
        # they come back 5x HIE + 1x HID. `reduce` on its own does not rename at
        # all -- the renaming is pdb4amber's job.
        ( cd "$D" || exit 0
          pdb4amber -i protein.pdb -o protein_prot.pdb --reduce > pdb4amber.log 2>&1
          if [[ -s protein_prot.pdb ]]; then
              # Strip the hydrogens OURSELVES, after the renaming. Keeping them makes
              # tleap fail with "Atom .R<PHE 68>.A<H3 23> does not have a type": this
              # protein has chain breaks, pdb4amber marks each as a new terminus and
              # reduce puts H1/H2/H3 there, which tleap's internal-residue templates
              # reject. Using pdb4amber's own --nohyd is NOT the fix -- it collapses
              # HID/HIE back to HIS and loses the assignment entirely. Doing it here
              # keeps the residue NAMES and lets tleap rebuild hydrogens from them.
              awk '!(/^(ATOM|HETATM)/ && (substr($0,77,2) ~ /H/ || substr($0,14,1)=="H"))' \
                  protein_prot.pdb > protein_noh.pdb
              mv protein_noh.pdb protein.pdb; rm -f protein_prot.pdb
          else
              echo "  !! pdb4amber --reduce failed; falling back to tleap defaults"
          fi
          awk '/^ATOM/ && $3=="CA" && substr($0,18,3) ~ /HI[DEP]/ \
               {print substr($0,18,3), substr($0,23,4)}' protein.pdb > his_states.txt || true )
    fi
    ( cd "$D" || exit 1

    {
      echo "source $FF_PROTEIN"
      echo "source $FF_WATER"
      # gaff2 is sourced even with no ligand so the loaded force-field stack is
      # identical in every system (S1/S3 have no ligand and source it too)
      echo "source $FF_LIGAND"
      [[ -n "$ION_FRCMOD" ]] && echo "loadamberparams $ION_FRCMOD"
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
      echo "addions sys $CATION 0"
      echo "addions sys $ANION 0"
      echo "addionsrand sys $CATION {nna} $ANION {ncl}"
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
    echo "  $(basename "$D"): $NATOM atoms, $NWAT waters, ${NION}x ${CATION}${ANION} (${CONC} M), net charge $CHG"
    python3 -c "
import sys
q=float('$CHG')
sys.exit('  ERROR: net charge %s is not neutral; PME would be wrong' % q) if abs(q)>0.01 else None" || exit 1
    if grep -qiE "Could not find|no parameters|Unknown residue|ATOM NOT FOUND" tleap_pass2.log; then
        echo "  ERROR: unparameterised atoms in $D"
        grep -iE "Could not find|no parameters|Unknown residue" tleap_pass2.log | head
        exit 1
    fi
    # hydrogen-mass repartitioning: move mass from heavy atoms into the hydrogens
    # they carry so the fastest bond vibrations slow down and 4 fs becomes stable.
    # WATER IS DELIBERATELY EXCLUDED (no `dowater`) -- it is already rigid under
    # SHAKE, so repartitioning it buys nothing and perturbs solvent dynamics.
    if [[ "$HMR" == "1" ]]; then
        printf 'HMassRepartition\noutparm system_hmr.prmtop\nquit\n' > hmr.parmed
        parmed -p system.prmtop -i hmr.parmed > hmr.log 2>&1 || true
        if [[ -s system_hmr.prmtop ]]; then
            echo "    + system_hmr.prmtop (4 fs timestep)"
        else
            echo "    !! HMR failed; see $D/hmr.log"; exit 1
        fi
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
