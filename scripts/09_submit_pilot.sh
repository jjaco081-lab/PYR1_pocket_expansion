#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 08:00:00
#SBATCH -J pe_pilot
#SBATCH -a 0-40
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/pilot_%a.log
#
# 09_submit_pilot.sh -- Arms 1 and 2 over the 41-variant panel.
# One array task per variant; each runs the apo cavity scan (arm 1) then the
# ternary ratchet scan (arm 2), nstruct=3 replicates each.
#
# NOTE: must run under SLURM. FastRelax exceeds the 1 GB login-node cgroup cap.
set -uo pipefail
PY=/bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python
S=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/scripts
i=$SLURM_ARRAY_TASK_ID
echo "=== variant index $i ==="
$PY $S/06_rosetta_cavity_scan.py --index $i --nstruct 3; echo "arm1 rc=$?"
$PY $S/07_rosetta_ratchet.py     --index $i --nstruct 3; echo "arm2 rc=$?"
echo "PILOT_DONE idx=$i"
