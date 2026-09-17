#!/bin/bash -l
#SBATCH -J panelthread
#SBATCH -p cutlerlab
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/panel_thread_%j.log
set -uo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python -u \
  scripts/225_panel_thread.py
test -s data/panel_mmgbsa/frames.json || { echo "NO frames.json"; exit 1; }
n=$(/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python -c "import json;print(len(json.load(open('data/panel_mmgbsa/frames.json'))))")
echo "ligands with both frames: $n (expect 4)"
test "$n" -eq 4
