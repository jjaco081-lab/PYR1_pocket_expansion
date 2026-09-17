#!/bin/bash -l
#SBATCH -J fpvenc
#SBATCH -p cutlerlab
#SBATCH -c 48
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/fpvenc_%j.log
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python scripts/233_fpocket_vs_enclosed.py 48
