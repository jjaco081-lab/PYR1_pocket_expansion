#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 8
#SBATCH --mem=12G
#SBATCH -t 7-00:00:00
#SBATCH -J pe_ratchet2
#SBATCH -a 0-5
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/ratchet2_%a.log
#
# 22_submit_ratchet_v2.sh -- Arm 2 re-run at nstruct=20 with the gate and latch
# released from coordinate constraints (see 21_rosetta_ratchet_v2.py).
#
# SHAPE OF THE JOB
#   Total work is 41 variants x 20 replicates = 820 FastRelax runs on the
#   474-residue ternary complex. Measured on the v2 smoke test (job 27275274):
#   ~20 min and ~1.03 GB RSS per replicate -- faster than the v1 pilot's ~50 min
#   because the gate/latch constraints are gone -- so ~273 CPU-hours total.
#
#   The account is limited to ~6 concurrently RUNNING jobs, so the work is split
#   across 6 array tasks of 8 cores each (48 workers, 12 GB/job = 8 x 1.03 GB
#   plus headroom).
#
#   WHY 8 CORES AND NOT MORE. Sizing here is set by what can actually be
#   scheduled, not by the account cap (384 CPUs). The partition's two nodes are
#   saturated in a lopsided way -- checked with scontrol show node:
#       r11: 20 CPUs idle but only ~7.6 GB SLURM-allocatable memory
#       r12: ~175 GB allocatable but only 10 CPUs idle
#   A 16-core/32 GB task fits NEITHER, and SLURM estimated its start ~14 h out.
#   An 8-core/12 GB task fits r12 immediately. Because the units are
#   work-stolen (below), a narrower job is never wrong -- it just takes fewer
#   units -- so the shape that STARTS beats the shape that is theoretically
#   faster. 48 workers finishes the set in ~6 h; 96 workers would finish in ~3 h
#   but only after waiting half a day to begin.
#
#   Re-check those two numbers before resubmitting; if the nodes have freed up,
#   raising -c and --mem in step is safe and needs no other change.
#
# WHY WORK-STEALING RATHER THAN A FIXED SPLIT
#   Tasks will not start together on a busy partition. A fixed split (task i owns
#   replicates 4i..4i+3) would leave the whole set waiting on whichever task
#   started last. Instead every (variant, replicate) unit is claimed atomically
#   with mkdir, which is atomic on GPFS: the first worker to create the claim
#   directory runs the unit, everyone else skips it instantly. Tasks that start
#   late simply do fewer units, and the set finishes as early as the available
#   cores allow. A unit whose run fails releases its claim so a later pass or a
#   resubmission retries it.
#
#   Re-running this script is therefore safe and idempotent: completed units are
#   skipped on sight (their JSON exists), so it doubles as the recovery path.
#
# NOTE: must run under SLURM. FastRelax exceeds the 1 GB login-node cgroup cap,
# and /scratch is node-local -- all paths must stay on /bigdata.
set -uo pipefail

PY=/bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python
P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
S=$P/scripts
OUT=$P/results/ratchet_v2
CLAIMS=$OUT/claims

NVAR=41          # panel size, must match variants.py
NREP=20          # replicates per variant
WORKERS=8        # parallel relax processes per array task (= #SBATCH -c)

mkdir -p "$CLAIMS" "$OUT/pdb" "$P/logs/units"
export PY P S OUT CLAIMS

# names in panel order, so a unit can test for its own finished output
mapfile -t NAMES < <($PY -c "
import sys; sys.path.insert(0,'$S')
from variants import PANEL
print('\n'.join(n for n,_ in PANEL))")
export NAMES_STR="${NAMES[*]}"

run_unit() {
    local v=$1 r=$2
    local names=($NAMES_STR)
    local name=${names[$v]}
    local json="$OUT/${name}__c${r}.json"

    [[ -s "$json" ]] && return 0                       # already done
    mkdir "$CLAIMS/v${v}_r${r}" 2>/dev/null || return 0 # claimed by someone else

    if $PY "$S/21_rosetta_ratchet_v2.py" \
           --index "$v" --nstruct 1 --rep-offset "$r" --chunk "$r" \
           >"$P/logs/units/unit_v${v}_r${r}.log" 2>&1; then
        echo "OK   v=$v r=$r $name"
    else
        echo "FAIL v=$v r=$r $name (claim released for retry)"
        rmdir "$CLAIMS/v${v}_r${r}" 2>/dev/null
    fi
}
export -f run_unit

i=$SLURM_ARRAY_TASK_ID
echo "=== array task $i: $WORKERS workers, $((NVAR*NREP)) units total ==="
date

# Emit all units, rotated by task id so the 6 tasks enter the list at different
# points and rarely contend for the same claim.
TOTAL=$((NVAR * NREP))
for ((k = 0; k < TOTAL; k++)); do
    j=$(((k + i * TOTAL / 6) % TOTAL))
    echo "$((j / NREP)) $((j % NREP))"
done | xargs -P "$WORKERS" -n 2 bash -c 'run_unit "$0" "$1"'

echo "RATCHET2_TASK_DONE id=$i rc=$?"
date
