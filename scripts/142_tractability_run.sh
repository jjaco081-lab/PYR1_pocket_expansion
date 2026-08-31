#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 1:58:00
#SBATCH -J tract
#SBATCH -a 0-361%4
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/tract_%a.log
#
# 142 -- Boltz-2 with the affinity head over the PROPERTY-MATCHED tractability set.
#
# 181 ligands that yielded a PYR1 sensor vs 181 that did not, matched on heavy
# atoms, cLogP, TPSA, HBD, HBA, charge and rotatable bonds so that no single
# descriptor separates them (§139: worst |AUC-0.5| = 0.021, against 0.202 for TPSA
# in a random draw). The benchmark therefore cannot be won on "too big for the
# pocket", which is exactly what Jannis asked to rule out.
#
# Both classes are co-folded against the SAME wild-type PYR1 construct -- a hit is a
# ligand some VARIANT bound, so giving hits an evolved pocket and non-hits wild-type
# would be a fatal asymmetry. The question is therefore whether WT co-folding
# predicts LIBRARY tractability: a proxy, one step from the label, but symmetric.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
BOLTZ=/bigdata/cutlerlab/jjaco081/conda_envs/boltz2/bin/boltz
CACHE=/bigdata/cutlerlab/jjaco081/conda_envs/boltz2/.boltz_cache
T=$ROOT/data/tractability
A3M=$T/pyr1_wt.a3m
[[ -s $A3M ]] || { echo "no MSA at $A3M -- run 141 first"; exit 1; }
# A single trailing NUL byte in the shared a3m silently failed all 362 runs
# once: boltz reported "Error: '\x00'" and every task exited 0 with no
# affinity output. ASSERT, do not repair -- a reader that tolerates the
# corruption would hide the next one (README 88e).
if LC_ALL=C grep -qaP "\x00" "$A3M"; then
    echo "MSA $A3M contains NUL bytes -- refusing to run"; exit 1
fi
I=$(printf "%04d" "${SLURM_ARRAY_TASK_ID:-0}")
SRC=$T/yaml_matched/lig$I.yaml
[[ -s $SRC ]] || { echo "no yaml for lig$I"; exit 0; }
OUTD=$T/out
mkdir -p "$OUTD" "$T/yaml_msa"
# Resume on the ARTEFACT, not the directory. The first (NUL-poisoned) run left
# 362 output directories containing only processed/manifest.json, and a
# directory-existence guard then skipped every task of the repaired rerun.
if compgen -G "$OUTD/boltz_results_lig$I/predictions/*/affinity*.json" > /dev/null; then
    echo "lig$I: affinity output present, skip"; exit 0
fi
rm -rf "$OUTD/boltz_results_lig$I"   # stale/partial: start clean
# point this ligand's YAML at the shared MSA
sed "s|^      sequence: \(.*\)$|      sequence: \1\n      msa: $A3M|" "$SRC" > "$T/yaml_msa/lig$I.yaml"
env -u PYTHONPATH $BOLTZ predict "$T/yaml_msa/lig$I.yaml" \
    --out_dir "$OUTD" --cache "$CACHE" \
    --recycling_steps 3 --diffusion_samples 5 \
    --diffusion_samples_affinity 5 --override 2>&1 | tail -3
J=$(find "$OUTD" -name "affinity_lig$I.json" | head -1)
[[ -s $J ]] && echo "=== lig$I affinity: $(cat "$J" | tr -d '\n ')" || echo "!! lig$I no affinity output"
