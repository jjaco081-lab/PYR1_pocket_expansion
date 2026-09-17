#!/bin/bash -l
#SBATCH -J linkcav
#SBATCH -p cutlerlab
#SBATCH -c 48
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/linker_cavity_%j.log
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python scripts/222_linker_cavity.py
