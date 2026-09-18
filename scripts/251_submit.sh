#!/bin/bash -l
#SBATCH -J relaxcav
#SBATCH -p cutlerlab
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/relaxcav_%j.log
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python scripts/251_relaxed_chamber.py
