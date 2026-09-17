#!/bin/bash -l
#SBATCH -J donors
#SBATCH -p cutlerlab
#SBATCH -c 48
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/donors_%j.log
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python scripts/234_donor_pockets.py 48
