#!/bin/bash -l
#SBATCH -J bbqual
#SBATCH -p cutlerlab
#SBATCH -c 48
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/bbqual_%j.log
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
CSVS=results/rfd3/truncated_all.csv
echo "scoring from: $CSVS"
/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python scripts/246_backbone_quality.py \
  --csv $CSVS --max-orange 5 --min-ss 0.45 --procs 48
