#!/bin/bash -l
#SBATCH -J juncall
#SBATCH -p cutlerlab
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/junctions_all_%j.log
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python scripts/237_junctions_all.py --flank 2
