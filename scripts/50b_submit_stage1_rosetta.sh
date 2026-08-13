#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 04:00:00
#SBATCH -J pe_s1_rosetta
#SBATCH -a 0-29
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/stage1_rosetta_%a.log
#
# 50b_submit_stage1_rosetta.sh -- the Rosetta FastDesign arm of stage 1 (23i).
#
# SHAPE OF THE JOB
#   6 arms (wt_mandi / polygly_mandi / wt_aba x dsm_hao / free) x 50 trajectories
#   = 300 FastRelax-with-design runs. Split as 30 array tasks of 10 trajectories,
#   which is 5 blocks per arm; 50_stage1_rosetta.py maps --unit onto (arm, block).
#
# SIZING, measured rather than guessed (job 27421259)
#   5.1 min per trajectory at --pack-shell 8 --rounds 1, 887 MB peak RSS. So a
#   task is ~51 min and the campaign is ~26 CPU-hours. cutlerlab had 210 idle
#   CPUs on r11 and 254 on r12 at submission, and the account cap is cpu=384,
#   so all 30 tasks start at once and the whole thing lands in ~1 h wall.
#   4 h walltime is deliberate slack, not an estimate -- a task that gets a
#   slower core must not be truncated, and results are written incrementally
#   after every trajectory so even a kill preserves finished work.
#
# WHY 1 CORE PER TASK
#   FastRelax is single-threaded. Asking for more cores would not speed a task
#   up and would make it harder to schedule. Parallelism comes from the array.
#
# TWO EARLIER SIZINGS WERE WRONG, both worth recording:
#   * Repacking the whole protein (no shell) did not finish ONE trajectory in
#     28 min. The pack shell is what makes this tractable.
#   * NeighborhoodResidueSelector with an 8 A cutoff selected only 3 residues,
#     because it measures between neighbour atoms and molfile_to_params gives a
#     29-atom ligand a single NBR atom. That froze the pocket solid, made every
#     trajectory identical, and produced a falsely tidy dG of -12.29. The shell
#     is now computed all-heavy-atom-to-all-heavy-atom: 43 repackable residues,
#     dG -19.11. See pocket_shell() in 50_stage1_rosetta.py.
#
# Merge with: python 51_stage1_rosetta_merge.py
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
  scripts/50_stage1_rosetta.py --unit "${SLURM_ARRAY_TASK_ID}" --n-per-unit 10
