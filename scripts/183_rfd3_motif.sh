#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 1:58:00
#SBATCH -J rfd3
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/rfd3_%j.log
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
OUT=$R/results/rfd3/batch01
mkdir -p "$OUT"
cd "$R" || exit 1
echo "python: $(which python)   rfd3: $(which rfd3)"
python -c "import torch;print('torch',torch.__version__,'cuda',torch.cuda.is_available())"

# Linker ranges are wide so the sampler can vary scaffold length, which is what
# varies the envelope. Native gaps are 63->84 = 20, 89->115 = 25, 117->148 = 30.
CONTIG="50-70/A60-63/15-30/A84-89/20-35/A115-117/25-40/A148-166/10-25"
echo "contig: $CONTIG"

rfd3 \
  +inputs.pdb=$R/data/pyr1_A.pdb \
  +out_dir=$OUT \
  +n_batches=1 \
  +diffusion_batch_size=10 \
  +seed=20260901 \
  +specification="{
    \"input_specification\": {
      \"dialect\": 2,
      \"contig\": \"$CONTIG\"
    }
  }"
rc=$?
echo "rfd3 rc=$rc  designs=$(ls $OUT/*.pdb 2>/dev/null | wc -l)"
exit $rc
