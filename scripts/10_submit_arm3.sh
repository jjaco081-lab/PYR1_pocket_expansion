#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 08:00:00
#SBATCH -J pe_arm3
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/arm3.log
# Arm 3 compute stage. The fetch stage must be run FIRST on the login node
# (compute nodes may lack outbound network):
#   python 08_homolog_cavities.py --stage fetch
PY=/bigdata/cutlerlab/jjaco081/conda_envs/esmfold2/bin/python
S=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/scripts
$PY $S/08_homolog_cavities.py --stage compute --min-tm 0.5 --max-n 150
echo "ARM3_DONE rc=$?"
