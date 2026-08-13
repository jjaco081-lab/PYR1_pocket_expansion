#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 00:30:00
#SBATCH -J pe_k59
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/k59_decomposition.log
#
# One pose, six single substitutions at position 59, ref2015 per-residue decomposition.
# Diagnoses WHY the ABA null arm deletes the K59 salt bridge; see the docstring in
# 55_k59_energy_decomposition.py for what has already been ruled out (pose, ligand
# ionisation, sequence prior).
set -euo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python \
    $ROOT/scripts/55_k59_energy_decomposition.py
