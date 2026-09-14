#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 1:58:00
#SBATCH -J rfd3b2
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/rfd3b2_%j.log
#
# 183 -- RFdiffusion3 motif scaffolding around PYR1's HAB1 interface.
#
# THE MOTIF is the four discontiguous segments README 9c identifies as the HAB1
# interface, taken as contiguous spans so a contig can express them:
#   A60-63    HFIK   helix-60s (interface 60, 61, 63)
#   A84-89    ISGLPA gate
#   A115-117  HRL    latch (interface 116, 117)
#   A148-166  PEGNSEDDTRMFADTVVKL   C-lobe (interface 148,151,155,156,158,159,162,166)
# All 19 interface residues verified present in data/pyr1_A.pdb with the expected
# identities before this was written -- chain A, residues 1-181, gaps only at
# 69-70 which fall inside a linker RFdiffusion rebuilds.
#
# ⚠ THE TRADE-OFF, stated: fixing A148-166 preserves the HAB1 contacts but also
# freezes F159/V163/V164, which are POCKET wall (README 66). So this design can
# enlarge the pocket everywhere except the C-lobe face. That is inherent to
# keeping the interface, not a choice made here.
#
# ⚠ WHAT THIS CANNOT SHOW (README 7, 85): the gate and latch are a two-state
# switch, and no design -- diffused, grafted or mutated -- can be shown
# computationally to cycle, because 85 closed the only route to the open/closed
# dG. Scaffolds are scored on ENVELOPE (106: depth and d/w, not cavity volume,
# since a bigger cavity from smaller side chains is obtainable in PYR1 far more
# cheaply). Switching remains a Y2H question.
#
# FIRST BATCH IS 10, for inspection before the remaining 90 -- Jannis's standing
# rule that generated structures are looked at before anything is built on them.
set -uo pipefail
module load gcc/12.2.0 2>/dev/null
module load cuda/12.1 2>/dev/null
source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda deactivate 2>/dev/null
conda activate /bigdata/cutlerlab/jjaco081/conda_envs/foundry
R=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
OUT=$R/results/rfd3/batch02
mkdir -p "$OUT"
cd "$R" || exit 1
echo "python: $(which python)   rfd3: $(which rfd3)"

# ---------------------------------------------------------------------------
# BATCH 2 -- three changes, each from an inspection Jannis made on batch 1.
#
# 1. HAB1 IS PRESENT. Jannis: designs 8 and 9 built alpha-helical bundles IN
#    FRONT of the gate and latch, where they would block both gate closure and
#    HAB1 binding. Batch 1 scaffolded PYR1 alone, so the diffusion model had no
#    reason to keep that face clear. 3QN1 carries PYR1 (chain A, 1-181) and HAB1
#    (chain B, 185-505) in the closed ternary complex; fixing all of chain B
#    makes the HAB1 surface a hard steric constraint instead of an afterthought.
#
# 2. MINIMAL LOOP ANCHORS, Jannis's span. Batch 1 fixed 32 residues (the four
#    interface segments only). This fixes 54 -- each loop plus the starts of its
#    flanking strands, including his 81-92 gate span, which captures six strand
#    residues in, the turn, and five out. All 19 HAB1-interface residues are
#    still covered. README 107c: the literal extend-to-structure reading would
#    fix 139 of 181 residues, which freezes the sheet and can only reproduce
#    PYR1 with new loops.
#
# 3. BURIAL CONDITIONING instead of secondary-structure bias. RFd3 dialect 2
#    exposes no sheet/helix conditioning (checked: inference/input_parsing.py
#    offers select_buried / select_partially_buried / select_exposed and
#    select_hotspots, nothing for SS). select_buried is the better lever anyway:
#    marking the POCKET-FACING motif residues buried asks the model to enclose
#    them, which is what a pocket is. 7 residues (60, 61, 87, 88, 89, 117, 159)
#    are BOTH pocket wall and HAB1 interface and are deliberately left
#    unconditioned -- they cannot be asked to be buried and exposed at once.
# ---------------------------------------------------------------------------
# ⚠ TWO NUMBERING TRAPS IN 3QN1.cif, BOTH FOUND THE HARD WAY.
# (a) HAB1 chain B is missing 26 residues, so a single span "B185-505" fails
#     with "Residue B288 not found in atom array". Four blocks are needed.
# (b) RFd3's CIF parser reads LABEL numbering, and 3QN1's chain A carries an
#     auth-label offset of -2. Under label numbering "A58-65" would have
#     silently selected the wrong residues -- the same 2-residue offset already
#     on record for 3QN1-derived CIFs.
# Both are removed by feeding a PDB written with AUTH numbering, where there is
# only one scheme. data/3QN1_complex_auth.pdb is generated with the four motif
# sequences asserted: YKHFIKSC / VIVISGLPANTS / IGGEHRLTNYK / DMPEGNSEDDTRMFADTVVKLNL.
# ⚠ CHAIN BREAKS ARE "/0" AND ARE NOT OPTIONAL. Without them every segment is
# one continuous chain, so the first cut asked RFd3 to build a single ~500-residue
# chain threading through PYR1 AND HAB1. It could not, and said so loudly: motif
# max CA deviation 33-50 A (batch 1: 0.1 A), 5-8 chain breaks per design, Rg
# 25-32 A against PYR1's 15.0. foundry/utils/components.py:104 handles "/0" as a
# pn_unit boundary. HAB1's four blocks are each their own unit so no linker is
# built through its 26 missing residues -- for a steric constraint the chain
# identity does not matter, only the atoms.
CONTIG="50-70,A58-65,15-30,A81-92,20-35,A111-121,25-40,A146-168,10-25,/0,B185-221,/0,B232-270,/0,B283-461,/0,B466-505"
# pocket-facing motif residues that are NOT also HAB1 interface
BURIED="A59,A62,A81,A83,A91,A92,A115,A120,A160,A163,A164,A167"
echo "contig: $CONTIG"
echo "buried: $BURIED"

rfd3 \
  inputs=null \
  out_dir=$OUT \
  n_batches=1 \
  diffusion_batch_size=10 \
  +seed=20260914 \
  +specification.input=$R/data/3QN1_complex_auth.pdb \
  +specification.contig=\'$CONTIG\' \
  +specification.select_buried=\'$BURIED\' \
  +specification.dialect=2
rc=$?
echo "rfd3 rc=$rc  designs=$(ls $OUT/*.pdb 2>/dev/null | wc -l)"
exit $rc
