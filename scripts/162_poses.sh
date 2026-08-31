#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0:40:00
#SBATCH -J poses
#SBATCH -a 0-29%4
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/pose_%a.log
#
# 162 -- Boltz-2 poses for the three agrochemicals that produced responders in
# Park's screen but have no crystal structure in PYR1: benzothiadiazole,
# benoxacor, fludioxonil. Receptor is PYR1(K59R), the actual library background.
#
# TEN SEEDS EACH, and the seeds are the point. README 91 scored 475 variants
# against a CRYSTAL mandipropamid pose. These three have no crystal pose, so the
# question is whether a PREDICTED one is stable enough to anchor the same rigid
# ligand protocol. If the pose moves between seeds the method needs
# crystallography and that is the answer; if it does not, the method transfers
# to novel ligands, which is Goal 2/3.
#
# The PYR1 MSA is reused from the tractability run (README 84b), so each job is
# ~71 s median.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
BOLTZ=/bigdata/cutlerlab/jjaco081/conda_envs/boltz2/bin/boltz
CACHE=/bigdata/cutlerlab/jjaco081/conda_envs/boltz2/.boltz_cache
A3M=$ROOT/data/tractability/pyr1_wt.a3m
[[ -s $A3M ]] || { echo "no MSA at $A3M"; exit 1; }
if LC_ALL=C grep -qaP "\x00" "$A3M"; then echo "MSA has NUL bytes -- refusing"; exit 1; fi
I=${SLURM_ARRAY_TASK_ID:-0}
NAME=$(python3 -c "import json;print(json.load(open('$ROOT/data/poses/jobs.json'))[$I]['compound'])")
SEED=$(python3 -c "import json;print(json.load(open('$ROOT/data/poses/jobs.json'))[$I]['seed'])")
OUTD=$ROOT/data/poses/out/${NAME}_s${SEED}
if compgen -G "$OUTD/boltz_results_*/predictions/*/*.cif" > /dev/null; then
    echo "$NAME seed $SEED: already predicted, skip"; exit 0
fi
rm -rf "$OUTD"; mkdir -p "$OUTD" "$ROOT/data/poses/yaml_msa"
Y=$ROOT/data/poses/yaml_msa/${NAME}.yaml
sed "s|^      sequence: \(.*\)$|      sequence: \1\n      msa: $A3M|" \
    "$ROOT/data/poses/yaml/${NAME}.yaml" > "$Y"
env -u PYTHONPATH $BOLTZ predict "$Y" --out_dir "$OUTD" --cache "$CACHE" \
    --seed "$SEED" --diffusion_samples 1 --output_format mmcif --use_msa_server 2>&1 | tail -3
ls "$OUTD"/boltz_results_*/predictions/*/*.cif >/dev/null 2>&1 \
    && echo "=== $NAME seed $SEED OK" || echo "!! $NAME seed $SEED produced no structure"
