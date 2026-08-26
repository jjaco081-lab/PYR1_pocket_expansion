#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 12:00:00
#SBATCH -J ligaware
#SBATCH -a 0-39
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/ligaware_%a.log
#
# 131c -- the ligand-aware gate test (README 76). Per variant: relax the shell
# protein-only, smina-dock the cognate coumarin, graph-match the pose onto the
# params molecule, relax protein+ligand jointly, score dG_bind. CPU only.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python -u \
    scripts/131_ligand_aware_score.py --chunk "$SLURM_ARRAY_TASK_ID" --nchunks 40
