#!/bin/bash
# 68_build_md_systems.sh -- take the rebuilt 191-residue structures through the
# canonical solvate/neutralise path, one table for every system.
#
#   bash scripts/68_build_md_systems.sh [system ...]
#
# WHY A NEW OUTPUT TREE
#   Output goes to data/md191/, NOT data/md/. The old tree holds 13 completed 300 ns
#   replicates (S1-S3 x3, S4 x3, S5 x1) and they are not disposable: they are the
#   evidence behind README 24. Building on top of them would destroy that record, and
#   the two generations are not comparable anyway -- different protein (191 aa vs a
#   gapped 178-179), different salt (KCl vs NaCl), different histidine assignment.
#   Keep both, label them, never pool them.
#
# WHAT IS DIFFERENT FROM THE FIRST GENERATION
#   protein   full 191 aa, no chain breaks, wild-type P2, tail 182-191 modelled
#   salt      0.15 M KCl                     (was NaCl)
#   His       reduce-assigned HID/HIE        (was blanket HIE)
#   order     solvate -> neutralise -> salt  (S5 and the first S9 had this backwards)
#   numbering asserted at build time, residue_map.json emitted beside the topology
#
# THE LIGANDS ARE NOT RE-DOCKED, AND THAT IS CHECKED NOT ASSUMED
#   67c asserts that the rebuilt receptor presents an unchanged pocket to the crystal
#   ABA -- same contacting residues, same distances to within 0.05 A -- because script
#   67 froze every residue it did not rebuild. That is what licenses reusing
#   data/md/A8S.mol2, which carries the CRYSTAL pose, and the non-cognate poses from
#   script 63. If 67c ever fails, this script must not be run until the poses are
#   regenerated.
set -euo pipefail

P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion

# tleap/pdb4amber are AmberTools: amber/22, NOT amber/22_mpi_cuda (which ships pmemd
# but no tleap). Loaded here rather than assumed on PATH -- nothing else in this
# script's environment provides it.
module load amber/22
command -v tleap >/dev/null || { echo "tleap not found after 'module load amber/22'"; exit 1; }
command -v pdb4amber >/dev/null || { echo "pdb4amber not found"; exit 1; }

# absolute, not $(dirname "$0"): SLURM copies the submitted script to
# /var/spool/slurmd/job<N>/slurm_script, so $0 does not point into the repo
source "$P/scripts/lib_solvate.sh"

# Rosetta's dump_pdb appends a pose-energies table as lines beginning '#', and the
# structures here come straight from script 67. Strip everything that is not a
# coordinate or chain-terminator record before handing the file to pdb4amber, rather
# than trusting it to ignore what it does not recognise.
sanitize () {
    local src=$1 dst=$2
    awk '/^(ATOM|HETATM|TER|END)/' "$src" > "$dst"
    local n_in n_out
    n_in=$(grep -c '^\(ATOM\|HETATM\)' "$src")
    n_out=$(grep -c '^\(ATOM\|HETATM\)' "$dst")
    [[ "$n_in" == "$n_out" ]] || { echo "  ERROR: sanitize dropped coordinates ($n_in -> $n_out)"; return 1; }
    return 0
}

STRUCT=$P/data/structures_191
OUT=$P/data/md191
LIGDIR=$P/data/md            # A8S.lib/.mol2/.frcmod and the non-cognate ligands live here

mkdir -p "$OUT"

# name | protein pdb | ligand mol2 (or -) | ligand frcmod (or -)
#
# S1  3K3K chain A -- genuinely apo AND open, so no ligand is correct here.
# S2  3QN1 chain A + the crystal ABA. A8S.mol2 carries the 3QN1 pose and 67c asserts
#     the rebuilt receptor still presents an identical pocket (19/19 contacts, within
#     0.05 A), which is what licenses reusing it rather than re-docking.
# S9  apo-closed BY DESIGN -- the missing cell of the 19b factorial (§27). The empty
#     pocket is the experiment, not an oversight.
#
# ⚠ S3 IS DELIBERATELY ABSENT. 3K3K is a MIXED dimer: chain A is apo-open but chain B
#   is CLOSED WITH ABA BOUND (§28g). Building it protein-only would simulate a
#   ligand-shaped protomer around an empty cavity -- exactly the artefact §28d's
#   ligand_shell() exists to prevent, and exactly what the first-generation
#   S3_apo_dimer did for 3 x 300 ns. It needs an ABA in 3K3K's frame first:
#   data/md/A8S.mol2 holds the 3QN1 pose and is in the WRONG FRAME for this dimer.
SYSTEMS=(
  "S1_apo_open|$STRUCT/pyr1_open_191.pdb|-|-"
  "S2_holo_closed|$STRUCT/pyr1_closed_191.pdb|$LIGDIR/A8S.mol2|$LIGDIR/A8S.frcmod"
  "S9_apo_closed|$STRUCT/pyr1_closed_191.pdb|-|-"
)

want=("$@")
selected () {
    [[ ${#want[@]} -eq 0 ]] && return 0
    local n=$1
    for w in "${want[@]}"; do [[ "$w" == "$n" ]] && return 0; done
    return 1
}

for row in "${SYSTEMS[@]}"; do
    IFS='|' read -r name prot lig frc <<< "$row"
    selected "$name" || continue
    echo "=== $name ==="
    [[ -f "$prot" ]] || { echo "  MISSING $prot -- run 67_build_pyr1.py and 67c first"; exit 1; }
    mkdir -p "$OUT/$name"
    clean=$OUT/$name/input_protein.pdb
    sanitize "$prot" "$clean" || exit 1
    if [[ "$lig" == "-" ]]; then
        build_system "$OUT/$name" "$clean"
    else
        [[ -f "$lig" && -f "$frc" ]] || { echo "  MISSING ligand files for $name"; exit 1; }
        build_system "$OUT/$name" "$clean" "$lig" "$frc"
    fi
    # carry the residue map next to the topology so analysis never re-derives offsets
    src=$(basename "$prot" .pdb)_residue_map.json
    [[ -f "$STRUCT/$src" ]] && cp "$STRUCT/$src" "$OUT/$name/residue_map.json"
done

echo
echo "built systems:"
for row in "${SYSTEMS[@]}"; do
    IFS='|' read -r name _ _ _ <<< "$row"
    selected "$name" || continue
    [[ -f "$OUT/$name/system.prmtop" ]] && \
        printf "  %-18s %s\n" "$name" "$(tr '\n' ' ' < "$OUT/$name/build_settings.txt")"
done
