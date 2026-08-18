#!/bin/bash
# 69b_verify_md191.sh -- check that every built md191 system is in the conformation
# and the occupancy its NAME claims, reading the topology cpptraj will actually use.
#
# WHY
#   The factorial only works if each cell is the cell it says it is, and this
#   project has already been burned once by trusting a name: `S3_apo_dimer` was
#   built from a protomer that is closed and was ABA-bound (README 28g), and four
#   memories recording the truth lost to one confident label. So the claims are
#   re-derived from system.prmtop and system.inpcrd, not from the directory name.
#
# HOW
#   cpptraj strips solvent and writes frame 1 as a PDB; every measurement is then
#   made on that PDB. An earlier version read `nativecontacts ... mindist` column 3
#   directly and got 0.00 A for every system -- a false PASS. Measuring from
#   coordinates we can see beats trusting a column index.
set -euo pipefail
module load amber/22

P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
MD=$P/data/md191
OUT=$MD/_verify
mkdir -p "$OUT"

for name in S1_apo_open S2_holo_closed S9_apo_closed S10_holo_open; do
    D=$MD/$name
    [[ -f $D/system.prmtop ]] || { echo "$name: NOT BUILT"; exit 1; }
    cat > "$OUT/$name.in" <<CPPTRAJ
parm $D/system.prmtop
trajin $D/system.inpcrd
strip :WAT,K+,Cl-,Na+
trajout $OUT/$name.frame1.pdb pdb
run
quit
CPPTRAJ
    cpptraj -i "$OUT/$name.in" > "$OUT/$name.log" 2>&1
    [[ -s $OUT/$name.frame1.pdb ]] || { echo "$name: cpptraj wrote no PDB"; exit 1; }
done

exec /bigdata/cutlerlab/jjaco081/conda_envs/pyr1_docking/bin/python "$P/scripts/69b_verify_md191.py"
