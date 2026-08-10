#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -J pe_smoke2
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/smoke2.log
PY=/bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python
S=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/scripts
$PY $S/06_rosetta_cavity_scan.py --index 0  --nstruct 1
$PY $S/06_rosetta_cavity_scan.py --index 36 --nstruct 1
echo "SMOKE2_DONE"
