#!/bin/bash
# 56b_md_build_804_1.sh -- two-pass tleap build for S5 (PYR1 + ABA + 804_1).
#
# Pass 1 solvates with no salt just to count waters; pass 2 rebuilds with
# 0.15 M NaCl computed from that count (ions per water = conc / 55.5). Same
# procedure 38b used for S1-S4, so S5 is directly comparable to the S4 WT
# ternary baseline it exists to be compared against.
#
# tleap is AmberTools, so `module load amber/22` (NOT amber/22_mpi_cuda, which is
# the pmemd.cuda build).
set -euo pipefail

CONC=${CONC:-0.15}
SYS=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/data/md/S5_804_1_ternary
cd "$SYS"

module load amber/22

run_pass () {   # $1 = nna, $2 = ncl, $3 = tag
    sed -e "s/{nna}/$1/" -e "s/{ncl}/$2/" tleap.in.template > tleap.in
    tleap -f tleap.in > "tleap_$3.log" 2>&1 || true
}

echo "== pass 1: solvate, count waters =="
run_pass 0 0 pass1
NWAT=$(grep -oP 'Added \K[0-9]+(?= residues)' tleap_pass1.log | tail -1 || true)
if [[ -z "${NWAT:-}" ]]; then
    echo "ERROR: could not read water count from tleap_pass1.log" >&2
    tail -30 tleap_pass1.log >&2
    exit 1
fi
NION=$(python3 -c "print(int(round($NWAT * $CONC / 55.5)))")
echo "   waters: $NWAT  ->  ${NION} Na+ / ${NION} Cl- at ${CONC} M"

echo "== pass 2: rebuild with salt =="
run_pass "$NION" "$NION" pass2

# tleap reports a nonzero exit for warnings, so verify by artefact and by charge
if [[ ! -f system.prmtop || ! -f system.inpcrd ]]; then
    echo "ERROR: tleap did not write system.prmtop/inpcrd" >&2
    grep -iE "error|fatal" tleap_pass2.log | head -20 >&2
    exit 1
fi

CHG=$(grep -oP 'Total unperturbed charge:\s*\K-?[0-9.]+' tleap_pass2.log | tail -1)
NATOM=$(grep -c -E '^(ATOM|HETATM)' system_solvated.pdb)
echo
echo "   net charge : $CHG"
echo "   atoms      : $NATOM"
echo "   waters     : $NWAT"

# a nonzero net charge silently ruins PME electrostatics -- refuse to hand it on
python3 -c "
import sys
q=float('$CHG')
if abs(q) > 0.01:
    sys.exit(f'ERROR: net charge {q} is not neutral; do not run MD on this.')
print('   neutral OK')
"

echo
echo "unparameterised-atom / missing-parameter check:"
if grep -qiE "Could not find|no parameters|Unknown residue|ATOM NOT FOUND" tleap_pass2.log; then
    grep -iE "Could not find|no parameters|Unknown residue|ATOM NOT FOUND" tleap_pass2.log | head
    echo "  ^^ RESOLVE THESE BEFORE RUNNING"
    exit 1
fi
echo "   none"
echo
echo "NEXT: bash scripts/56c_md_stage_804_1.sh   (copies the S4 min/heat/eq/prod"
echo "      protocol into rep0 and submits to preempt_gpu)"
