#!/bin/bash -l
#SBATCH -J anchpose
#SBATCH -p cutlerlab
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 04:00:00
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/anchor_place_real_%j.log
set -uo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/dockenv/bin/python -u \
  scripts/209_anchor_place.py --all
# assert on the RESULT, not the exit code: poses must exist for all 4 ligands
n=$(ls results/anchor_place/*_ketone.json 2>/dev/null | wc -l)
echo "ketone-arm pose files written: $n (expect 4)"
test "$n" -eq 4 || { echo "FAIL: poses missing"; exit 1; }
