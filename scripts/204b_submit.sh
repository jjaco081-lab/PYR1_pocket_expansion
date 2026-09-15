#!/bin/bash
#SBATCH --job-name=latchadv
#SBATCH --partition=cutlerlab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=logs/latchadv_%j.log
# CPU only, by design: the panel is conformationally rigid (median 2 rotatable
# bonds, and anthrone at 1 uM has exactly one conformer), so RDKit ETKDG already
# enumerates the space exhaustively and no GPU-hours are needed. All four GPUs
# stay with RFd3.
set -euo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/foundry/bin/python -u scripts/204_latch_advance.py
test -s results/latch_advance/latch_advance.json || { echo "NO OUTPUT"; exit 1; }
echo "OK: $(python -c "import json;print(len(json.load(open('results/latch_advance/latch_advance.json'))))" 2>/dev/null || echo '?') rows"
