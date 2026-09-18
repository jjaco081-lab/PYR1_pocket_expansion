#!/bin/bash -l
#SBATCH -J orange
#SBATCH -p cutlerlab
#SBATCH -c 48
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/orange_%j.log
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python scripts/245_orange_filter.py --max-orange 5 --procs 48
