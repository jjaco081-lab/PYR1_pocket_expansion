#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0:30:00
#SBATCH -J tractsmk
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/tract_smoke.log
set -uo pipefail
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion || exit 1
rm -rf data/tractability/out/boltz_results_lig0000
SLURM_ARRAY_TASK_ID=0 bash scripts/142_tractability_run.sh
