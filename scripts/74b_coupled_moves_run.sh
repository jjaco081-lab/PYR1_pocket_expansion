#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 04:00:00
#SBATCH -J cmoves
#SBATCH -a 0-99
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/cmoves_%a.log
#
# 74b_coupled_moves_run.sh -- stage 1 re-run with coupled moves (README 35).
#
# 100 tasks: 0-49 are the mandipropamid arm, 50-99 the ABA ligand-swap null. Each is
# one independent CoupledMoves trajectory of 1000 trials, ~7 min, seeded from its
# replicate index so the set is reproducible.
#
# 50 trajectories per arm matches §23j exactly, so the only thing that changed
# between the two experiments is the SAMPLER -- fixed backbone there, sequence +
# side chain + backbone + ligand pose here.
set -uo pipefail
P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
PY=/bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python

if [[ $SLURM_ARRAY_TASK_ID -lt 50 ]]; then
    ARM=mandi; REP=$SLURM_ARRAY_TASK_ID
else
    ARM=aba;   REP=$(( SLURM_ARRAY_TASK_ID - 50 ))
fi
echo "=== $(date) task $SLURM_ARRAY_TASK_ID -> arm=$ARM rep=$REP ==="
$PY "$P/scripts/74_coupled_moves.py" "$ARM" "$REP"
rc=$?
echo "CM_TASK_DONE arm=$ARM rep=$REP rc=$rc"
date
exit $rc
