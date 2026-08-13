#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 04:00:00
#SBATCH -J pe_s1_ros_fn
#SBATCH -a 0-29
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/stage1_rosetta_fn_%a.log
#
# 54_submit_stage1_rosetta_fn.sh -- re-run all six stage-1 arms with the favor-native
# weight chosen by 53_favornative_pick.py from the NULL arm alone.
#
# This is 50b with one changed knob. Everything else -- 6 arms x 50 trajectories,
# 8 A pack shell, 1 FastRelax round, seeds 1000+17*traj -- is identical to job
# 27421344, so the two campaigns differ only in the term under test.
#
# USAGE (FN is required; there is deliberately no default, so nobody re-runs this
# at a weight nobody chose):
#   sbatch --export=ALL,FN=0.5 scripts/54_submit_stage1_rosetta_fn.sh
#
# Results land in results/stage1_rosetta_fn${FN}/, leaving job 27421344 intact for
# the side-by-side comparison. Merge with:
#   python scripts/51_stage1_rosetta_merge.py --indir results/stage1_rosetta_fn${FN}
#
# READ THE FLOOR CHECK FIRST (see 53): the mandipropamid arm must still mutate F108
# at high frequency. A bonus that quiets the null AND the mandi arm has not fixed
# the benchmark, it has muted it. If F108A drops toward the null's rate, the weight
# is too high and the correct move is to report the arm inconclusive, not to shop
# for a weight that produces a nicer table.
set -euo pipefail

if [[ -z "${FN:-}" ]]; then
    echo "ERROR: FN unset. Submit with: sbatch --export=ALL,FN=<weight> $0" >&2
    exit 1
fi

ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
PY=/bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python
OUTDIR=$ROOT/results/stage1_rosetta_fn${FN}
mkdir -p "$OUTDIR"

echo "task $SLURM_ARRAY_TASK_ID: favor_native=$FN -> $OUTDIR"
$PY $ROOT/scripts/50_stage1_rosetta.py \
    --unit "$SLURM_ARRAY_TASK_ID" \
    --n-per-unit 10 \
    --favor-native "$FN" \
    --outdir "$OUTDIR"
