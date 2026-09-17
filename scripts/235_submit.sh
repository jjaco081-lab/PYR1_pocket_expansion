#!/bin/bash -l
#SBATCH -J donorlin
#SBATCH -p cutlerlab
#SBATCH -c 48
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/donorlin_%j.log
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python scripts/235_donor_by_lining.py 48
