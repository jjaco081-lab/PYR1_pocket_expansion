#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 08:00:00
#SBATCH -J nc_build
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/nc_build.log
#
# 64_build_noncognate.sh -- parameterise and build the §27 ligand-dependence systems.
#
# BUILDS
#   S6_imperatorin   pose0/1/2   furanocoumarin     20 heavy, MW 270.3
#   S7_flutamide     pose0/1/2   nitroaromatic      19 heavy, MW 276.2
#   S8_estradiol     pose0/1/2   steroid            20 heavy, MW 272.4
#   S9_apo_closed                closed PYR1, NO ligand
#
# S9 is the control that matters most: README 24e flagged that the existing set
# confounds conformation with ligand occupancy, because S1 is apo AND open while S2
# is holo AND closed, and no apo-closed system exists. Without S9, "closed is rigid"
# and "ABA rigidifies" cannot be separated -- so S9 makes the SIX REPLICATES WE
# ALREADY HAVE interpretable, independent of anything the new ligands do.
#
# The protein is `data/md/S2_holo_closed/protein.pdb` for every system, byte for
# byte, so the only difference from S2 is the ligand. Any other choice reintroduces
# the confound this is meant to remove.
#
# ff19SB / OPC / GAFF2 + AM1-BCC, identical to README 19b. All three ligands are
# NEUTRAL (verified in 63), so -nc 0 throughout and no charge correction is needed.
#
# tleap is AmberTools: `module load amber/22`, NOT amber/22_mpi_cuda.
set -uo pipefail
module load amber/22

P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
MD=$P/data/md
NC=$P/data/noncognate
PROT=$MD/S2_holo_closed/protein.pdb
CONC=0.15

[[ -s "$PROT" ]] || { echo "missing receptor $PROT"; exit 1; }

declare -A SYS=( [Imperatorin]=S6_imperatorin [Flutamide]=S7_flutamide
                 [Alpha-Estradiol]=S8_estradiol )
declare -A RES=( [Imperatorin]=IMP [Flutamide]=FLU [Alpha-Estradiol]=EST )

build_one () {   # $1 = system dir, $2 = ligand mol2 (or "" for apo), $3 = frcmod
    local D=$1 LIGMOL2=$2 FRC=$3
    mkdir -p "$D"; cd "$D" || return 1
    cp "$PROT" protein.pdb

    if [[ -n "$LIGMOL2" ]]; then
        cat > tleap.in.template <<EOF
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.gaff2
# NO frcmod.ions234lm_126_opc -- that is TIP3P-era naming and does not exist for
# OPC; leaprc.water.opc has already loaded frcmod.ionslm_126_opc (README 19b).
loadamberparams $FRC
lig = loadmol2 $LIGMOL2
prot = loadpdb protein.pdb
sys = combine { prot lig }
# ORDER IS COPIED FROM S2's tleap.in AND MUST NOT BE REARRANGED: solvate, then
# neutralise, then add bulk salt. Neutralising BEFORE solvation places counter-ions
# in vacuum against the solute instead of substituting them into water, which is a
# different starting state -- and these systems exist only to be compared with S2.
solvateoct sys OPCBOX 12.0
addions sys Na+ 0
addions sys Cl- 0
addionsrand sys Na+ {nna} Cl- {ncl}
charge sys
check sys
saveamberparm sys system.prmtop system.inpcrd
savepdb sys system_solvated.pdb
quit
EOF
    else
        cat > tleap.in.template <<EOF
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.gaff2
# gaff2 is sourced even with no ligand present, so the force-field stack loaded is
# identical across every system in the comparison (S1 and S2 both source it too).
prot = loadpdb protein.pdb
sys = combine { prot }
# same order as S2/S1 -- see the note in the ligand branch
solvateoct sys OPCBOX 12.0
addions sys Na+ 0
addions sys Cl- 0
addionsrand sys Na+ {nna} Cl- {ncl}
charge sys
check sys
saveamberparm sys system.prmtop system.inpcrd
savepdb sys system_solvated.pdb
quit
EOF
    fi

    # two-pass: count waters with no salt, then rebuild at 0.15 M
    sed -e "s/{nna}/0/" -e "s/{ncl}/0/" tleap.in.template > tleap.in
    tleap -f tleap.in > tleap_pass1.log 2>&1 || true
    local NWAT
    NWAT=$(grep -oP 'Added \K[0-9]+(?= residues)' tleap_pass1.log | tail -1 || true)
    if [[ -z "${NWAT:-}" ]]; then
        echo "  ERROR: no water count in $D/tleap_pass1.log"; tail -20 tleap_pass1.log; return 1
    fi
    local NION
    NION=$(python3 -c "print(int(round($NWAT * $CONC / 55.5)))")
    sed -e "s/{nna}/$NION/" -e "s/{ncl}/$NION/" tleap.in.template > tleap.in
    tleap -f tleap.in > tleap_pass2.log 2>&1 || true

    if [[ ! -f system.prmtop ]]; then
        echo "  ERROR: tleap wrote no prmtop in $D"
        grep -iE "error|fatal" tleap_pass2.log | head; return 1
    fi
    local CHG NATOM
    CHG=$(grep -oP 'Total unperturbed charge:\s*\K-?[0-9.]+' tleap_pass2.log | tail -1)
    NATOM=$(grep -c -E '^(ATOM|HETATM)' system_solvated.pdb)
    echo "  $(basename "$D"): $NATOM atoms, $NWAT waters, ${NION} NaCl, net charge $CHG"
    python3 -c "
import sys
q=float('$CHG')
sys.exit('  ERROR: net charge %.3f is not neutral' % q) if abs(q)>0.01 else None"
    if grep -qiE "Could not find|no parameters|Unknown residue|ATOM NOT FOUND" tleap_pass2.log; then
        echo "  ERROR: unparameterised atoms in $D"
        grep -iE "Could not find|no parameters|Unknown residue" tleap_pass2.log | head
        return 1
    fi
    return 0
}

fail=0

echo "=== ligand parameterisation (GAFF2 / AM1-BCC) ==="
for LIG in Imperatorin Flutamide Alpha-Estradiol; do
    RN=${RES[$LIG]}
    for POSE in 0 1 2; do
        SDF=$NC/${LIG}_pose${POSE}.sdf
        [[ -s "$SDF" ]] || { echo "MISSING $SDF -- run 63 first"; fail=1; continue; }
        WORK=$NC/params/${LIG}_pose${POSE}
        mkdir -p "$WORK"; cd "$WORK" || exit 1
        if [[ ! -s ${RN}.mol2 ]]; then
            echo "  antechamber $LIG pose$POSE ..."
            # maxcyc=0 makes sqm a SINGLE POINT. Without it antechamber runs a
            # geometry optimisation and hands back a relaxed conformer, which would
            # silently discard the docked pose these runs are built around.
            antechamber -i "$SDF" -fi sdf -o ${RN}.mol2 -fo mol2 \
                        -c bcc -nc 0 -at gaff2 -rn $RN -s 0 \
                        -ek "maxcyc=0, qm_theory='AM1', scfconv=1.d-10" > ac.log 2>&1
        fi
        if [[ ! -s ${RN}.mol2 ]]; then
            echo "  ERROR: antechamber failed for $LIG pose$POSE"; tail -15 ac.log; fail=1; continue
        fi
        [[ -s ${RN}.frcmod ]] || parmchk2 -i ${RN}.mol2 -f mol2 -o ${RN}.frcmod -s gaff2
        # parmchk2 flags guessed parameters with ATTN; they are not fatal but must be seen
        if grep -q "ATTN" ${RN}.frcmod; then
            echo "  NOTE: $LIG pose$POSE has ATTN (guessed) parameters:"
            grep -c "ATTN" ${RN}.frcmod
        fi
    done
done

echo
echo "=== system builds ==="
for LIG in Imperatorin Flutamide Alpha-Estradiol; do
    RN=${RES[$LIG]}; S=${SYS[$LIG]}
    for POSE in 0 1 2; do
        W=$NC/params/${LIG}_pose${POSE}
        [[ -s $W/${RN}.mol2 ]] || { fail=1; continue; }
        build_one "$MD/${S}/pose${POSE}" "$W/${RN}.mol2" "$W/${RN}.frcmod" || fail=1
    done
done

build_one "$MD/S9_apo_closed" "" "" || fail=1

echo
echo "NC_BUILD_DONE fail=$fail"
exit $fail
