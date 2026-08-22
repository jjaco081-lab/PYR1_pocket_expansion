#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 7-00:00:00
#SBATCH -J nc_md
#SBATCH -a 0-2
#SBATCH --open-mode=append
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/nc_md_%a.log
#
# 65_md_run_noncognate.sh -- the §27 ligand-dependence controls.
#
# THREE JOBS, FOUR 300 ns RUNS EACH, SEQUENTIALLY
#   task 0 : Imperatorin     pose0/1/2  + S9_apo_closed rep0
#   task 1 : Flutamide       pose0/1/2  + S9_apo_closed rep1
#   task 2 : Alpha-Estradiol pose0/1/2  + S9_apo_closed rep2
#
# ~44k atoms each, ~320 ns/day on an A100, so 4 x 300 ns is about 3.8 days inside a
# 7-day wall. Runs are sequential inside one job by request -- there is no time
# pressure, and packing them keeps the job count at three.
#
# WHY THE THREE LIGANDS
#   Matched to ABA in size (19-20 heavy atoms, MW 270-276 vs ABA's 19 / 260) but
#   spanning three chemical classes: furanocoumarin, nitroaromatic anilide, steroid.
#   Size is held constant so any difference in loop dynamics is about chemistry
#   rather than bulk. All three are Tian screen hits that need mutations, so all are
#   genuinely non-cognate for WT.
#
#   Each ligand's three replicates start from three DIFFERENT DOCKED POSES, not
#   three seeds of one pose -- registered in §27c so pose sensitivity is measured
#   instead of assumed.
#
# S9_apo_closed WAS REMOVED FROM THESE JOBS (2026-08-21)
#   It has since been run to completion as THREE 300 ns replicates in the md191
#   tree, and §29 reports the result: apo-closed does not open, and has the most
#   rigid gate of the 15 units measured. Re-running it here would also have used
#   the OLD data/md build (NaCl, 178 residues), so it would not have been
#   comparable to the md191 answer anyway. The original rationale is kept below
#   because it explains why the cell mattered.
#
# (historical) WHY S9_apo_closed WAS IN EVERY JOB
#   It is the missing cell of §19b. S1 is apo AND open, S2 is holo AND closed, so
#   "closed is rigid" and "ABA rigidifies" are confounded (§24e). S9 separates them
#   and thereby makes the six replicates we ALREADY have interpretable. Spreading
#   its three replicates one per job means any single job completing still delivers
#   one apo-closed replicate.
#
# READING THIS LATER -- registered in §27c, repeated here so it travels with the job
#   * absence of ligand release is NOT evidence of binding. Residence times are
#     microseconds to milliseconds; 300 ns cannot sample unbinding. Release IF SEEN
#     is informative; not seeing it says nothing.
#   * if closed+ligand looks like closed+ABA for all three chemotypes, the §24b
#     stability filter is conformation-only and cannot rank ligands.
#
# RESTART SAFETY: identical scheme to 42/57. Every resume writes a fresh
# prod_cont_NNN.nc so segments accumulate rather than overwrite, and the remaining
# step count is recomputed from the restart file's own clock, so a requeue never
# runs a further full 300 ns.
set -uo pipefail
module load amber/22_mpi_cuda >/dev/null 2>&1

P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
MD=$P/data/md
# 150 ns, revised down from 300 (2026-08-21). The original length was sized to
# catch GATE OPENING, but §24 and §29 have since shown open and closed are BOTH
# kinetically trapped at 300 ns and that apo-closed does not open either -- so
# opening is not observable on this timescale for ANY ligand, and buying more of it
# buys nothing. The readout that IS ligand-discriminating is LIGAND POSE STABILITY,
# validated today in §41/§42 where ABA and mandipropamid separated cleanly within
# 10 ns (max excursion 4.88 vs 2.45 A implicit). 150 ns is 15x that margin.
# Resumable, so extending later costs only the extra nanoseconds.
PROD_NS=150
# EQUIL_PS comes from lib_mdinputs.sh (1700 = heat 200 + eq1 500 + eq2 1000)
# NOTE: sourced by ABSOLUTE path, not "$(dirname $0)". SLURM copies the batch
# script to /var/spool/slurmd/job<N>/slurm_script, so $0 points at the spool copy
# and dirname finds no lib_*.sh -- the source fails silently and every call becomes
# "command not found". This must also come AFTER P is set.
source "$P/scripts/lib_mdinputs.sh"
md_settings_check

TARGET_PS=$(( EQUIL_PS + PROD_NS * 1000 ))

# "<dir holding system.prmtop>|<dir to run in>"
case $SLURM_ARRAY_TASK_ID in
  0) RUNS=( "$MD/S6_imperatorin/pose0|$MD/S6_imperatorin/pose0"
            "$MD/S6_imperatorin/pose1|$MD/S6_imperatorin/pose1"
            "$MD/S6_imperatorin/pose2|$MD/S6_imperatorin/pose2" ) ;;
  1) RUNS=( "$MD/S7_flutamide/pose0|$MD/S7_flutamide/pose0"
            "$MD/S7_flutamide/pose1|$MD/S7_flutamide/pose1"
            "$MD/S7_flutamide/pose2|$MD/S7_flutamide/pose2" ) ;;
  2) RUNS=( "$MD/S8_estradiol/pose0|$MD/S8_estradiol/pose0"
            "$MD/S8_estradiol/pose1|$MD/S8_estradiol/pose1"
            "$MD/S8_estradiol/pose2|$MD/S8_estradiol/pose2" ) ;;
  *) echo "bad array index"; exit 1 ;;
esac

echo "=== $(date)  task $SLURM_ARRAY_TASK_ID  job $SLURM_JOB_ID  restarts=${SLURM_RESTART_COUNT:-0} ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

rst_time () {    # echo simulated ps in $1, or nothing if unreadable
    ncdump -v time "$1" 2>/dev/null | awk '/^ time =/{gsub(/[^0-9.]/,"",$3); print $3; exit}'
}

run_stage () {   # name  prev_rst  [ref]
    local nm=$1 prev=$2
    if [[ -s ${nm}.rst7 ]]; then echo "    $nm done, skipping"; return 0; fi
    echo "    running $nm ..."
    local ref=()
    [[ "${3:-}" == "ref" ]] && ref=(-ref "$prev")
    pmemd.cuda -O -i ${nm}.in -p "$TOP" -c "$prev" \
               -o ${nm}.out -r ${nm}.rst7 -x ${nm}.nc -inf ${nm}.info "${ref[@]}" \
        || { echo "    !! $nm FAILED"; record_seed "$nm"; return 2; }
    record_seed "$nm"
}

overall=0
for entry in "${RUNS[@]}"; do
    TOPDIR=${entry%%|*}
    D=${entry##*|}
    TOP=$(topology_for_run "$TOPDIR") || { overall=1; continue; }
    CRD=$TOPDIR/system.inpcrd
    echo
    echo "--- $(date)  $D ---"
    if [[ ! -s "$TOP" ]]; then
        echo "  missing topology for $TOPDIR -- run 64_build_noncognate.sh first; SKIPPING"
        overall=1; continue
    fi
    mkdir -p "$D"; cd "$D" || { overall=1; continue; }
    write_md_inputs

    run_stage min1 "$CRD" ref || { overall=1; continue; }
    run_stage min2 min1.rst7  || { overall=1; continue; }
    run_stage heat min2.rst7 ref || { overall=1; continue; }
    run_stage eq1  heat.rst7 ref || { overall=1; continue; }
    run_stage eq2  eq1.rst7      || { overall=1; continue; }

    START=eq2.rst7
    NOW=$EQUIL_PS
    if [[ -s prod.rst7 ]]; then
        T=$(rst_time prod.rst7)
        if [[ -z "$T" && -s prod_backup.rst7 ]]; then
            T=$(rst_time prod_backup.rst7)
            [[ -n "$T" ]] && { cp prod_backup.rst7 prod.rst7; echo "  recovered from backup at $T ps"; }
        fi
        if [[ -n "$T" ]]; then START=prod_resume.rst7; NOW=$T; cp prod.rst7 "$START"; fi
    fi
    REMAIN_PS=$(awk -v a=$TARGET_PS -v b=$NOW 'BEGIN{printf "%.0f", a-b}')
    if [[ "$REMAIN_PS" -le 0 ]]; then
        echo "  already complete ($NOW / $TARGET_PS ps)"
        echo "NC_RUN_DONE dir=$D rc=0 time=$NOW/$TARGET_PS"
        continue
    fi
    NSTEPS=$(( REMAIN_PS * STEPS_PER_PS ))
    echo "  production: at $NOW ps -> $REMAIN_PS ps remaining ($NSTEPS steps)"

    write_prod_input "$NSTEPS"

    if [[ "$START" == "eq2.rst7" ]]; then
        OUT_NC=prod.nc; OUT_O=prod.out
    else
        n=1; while [[ -e $(printf "prod_cont_%03d.nc" $n) ]]; do n=$((n+1)); done
        OUT_NC=$(printf "prod_cont_%03d.nc" $n); OUT_O=$(printf "prod_cont_%03d.out" $n)
        cp prod.rst7 prod_backup.rst7
    fi
    echo "  writing $OUT_NC"
    pmemd.cuda -O -i prod.in -p "$TOP" -c "$START" \
               -o $OUT_O -r prod.rst7 -x $OUT_NC -inf prod.info
    rc=$?
    s=$(grep -m1 "Setting random seed to" $OUT_O 2>/dev/null | awk '{print $NF}')
    [[ -n "$s" ]] && echo "$(date -Is) ${OUT_O%.out} ig=$s" >> seeds.txt
    FIN=$(rst_time prod.rst7)
    echo "NC_RUN_DONE dir=$D rc=$rc time=${FIN:-unknown}/${TARGET_PS} ps"
    [[ $rc -ne 0 ]] && overall=1
done

echo
echo "=== $(date)  task $SLURM_ARRAY_TASK_ID finished, overall=$overall ==="
exit $overall
