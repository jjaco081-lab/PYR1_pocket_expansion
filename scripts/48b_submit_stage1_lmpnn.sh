#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J pe_s1_lmpnn
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/stage1_lmpnn_%j.log
#
# 48b_submit_stage1_lmpnn.sh -- re-run stage 1 on LigandMPNN after the HETATM
# column bug.
#
# WHY THIS EXISTS AT ALL. Job 27392339 ran this same test on 2026-08-11 and
# returned a clean-looking negative. It was invalid: 47_stage1_inputs.py wrote
# ligand records one column left of spec, ProDy read the '3' of '3UZ' as an
# altLoc, and its default altloc='A' filter dropped all 29 mandipropamid atoms.
# Both mandipropamid arms ran APO. ABA survived only because 'A8S' happens to
# start with 'A'. So the ligand-swap null was the sole arm containing a ligand,
# and `delta = P_mandi - P_aba` measured apo-minus-holo. Retracted; see 23h.
#
# 48_stage1_run.py now runs preflight() before anything else and exits if ProDy
# cannot see every ligand atom under its default settings.
#
# COST. 6 arms (3 structures x 2 alphabets) x 50 sequences. The bad run took
# ~4 minutes wall on CPU; the ligand adds context but not much. 2 h is slack,
# not an estimate.
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python scripts/48_stage1_run.py "$@"
