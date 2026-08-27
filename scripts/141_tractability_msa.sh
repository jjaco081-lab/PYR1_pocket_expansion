#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 0:40:00
#SBATCH -J tractmsa
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/tract_msa.log
#
# 141 -- generate the wild-type PYR1 MSA ONCE.
#
# All 362 tractability runs use the SAME receptor, so paying the public MSA server
# 362 times would dominate the cost, hammer a shared service, and introduce
# run-to-run variation in the one input that should be identical across the
# benchmark. One warm-up run produces the a3m; 142 reuses it for every ligand.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
BOLTZ=/bigdata/cutlerlab/jjaco081/conda_envs/boltz2/bin/boltz
CACHE=/bigdata/cutlerlab/jjaco081/conda_envs/boltz2/.boltz_cache
T=$ROOT/data/tractability
mkdir -p "$T/msa_warmup"
env -u PYTHONPATH $BOLTZ predict "$T/yaml_matched/lig0000.yaml" \
    --out_dir "$T/msa_warmup" --cache "$CACHE" \
    --use_msa_server --recycling_steps 3 --diffusion_samples 1 \
    --override 2>&1 | tail -5
A3M=$(find "$T/msa_warmup" -name "*.a3m" | head -1)
if [[ -z "$A3M" ]]; then echo "!! no a3m produced"; exit 1; fi
cp "$A3M" "$T/pyr1_wt.a3m"
echo "=== MSA saved: $T/pyr1_wt.a3m ($(wc -l < "$T/pyr1_wt.a3m") lines)"
echo "    query line: $(sed -n '2p' "$T/pyr1_wt.a3m" | cut -c1-60)..."
