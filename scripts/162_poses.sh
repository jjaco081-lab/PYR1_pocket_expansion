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
# The MSA query line MUST equal the input sequence. With the WT a3m against the
# K59R input, boltz silently abandons the alignment and folds single-sequence:
# median pLDDT 49 instead of 96, verified by running the identical job with the
# WT sequence (96.0-96.7). One residue, no error message (README 92d).
A3M=$ROOT/data/tractability/pyr1_k59r.a3m
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
python3 - "$A3M" "$ROOT/data/poses/yaml/${NAME}.yaml" <<'PY' || exit 1
import sys, re
q=open(sys.argv[1]).read().split("\n")[1].strip()
s=[l.split("sequence: ")[1].strip() for l in open(sys.argv[2]) if "sequence:" in l][0]
if q!=s:
    d=[i for i,(a,b) in enumerate(zip(q,s)) if a!=b]
    print(f"MSA query != input sequence at {d[:5]} (len {len(q)} vs {len(s)}) -- refusing")
    sys.exit(1)
print("MSA query matches the input sequence")
PY
sed "s|^      sequence: \(.*\)$|      sequence: \1\n      msa: $A3M|" \
    "$ROOT/data/poses/yaml/${NAME}.yaml" > "$Y"
# NO --use_msa_server: it makes boltz ignore the `msa:` path supplied in the
# YAML, and the first run of this script collapsed to median pLDDT 49.4 with a
# 14.0 A core backbone RMSD -- a misfolded receptor, not a pose result. The
# tractability array, identical but without the flag, gets pLDDT 94.4.
env -u PYTHONPATH $BOLTZ predict "$Y" --out_dir "$OUTD" --cache "$CACHE" \
    --seed "$SEED" --diffusion_samples 1 --output_format mmcif 2>&1 | tail -3
C=$(ls "$OUTD"/boltz_results_*/predictions/*/*.cif 2>/dev/null | head -1)
if [[ -z $C ]]; then echo "!! $NAME seed $SEED produced no structure"; exit 1; fi
# ASSERT THE FOLD before anything downstream trusts the ligand position.
python3 - "$C" "$NAME" "$SEED" <<'PY'
import sys
p,name,seed=sys.argv[1],sys.argv[2],sys.argv[3]
cols,rows,inl={},[],False
for line in open(p):
    if line.startswith('_atom_site.'): cols[line.strip().split('.')[1]]=len(cols); inl=True; continue
    if inl:
        if line.startswith(('#','loop_','_')):
            if rows: break
            continue
        f=line.split()
        if len(f)>=len(cols): rows.append(f)
v=[float(f[cols['B_iso_or_equiv']]) for f in rows
   if f[cols['group_PDB']]=='ATOM' and f[cols['label_atom_id']]=='CA']
v.sort(); med=v[len(v)//2] if v else 0
ok = med >= 80
print(f"=== {name} seed {seed} median pLDDT {med:.1f} -- {'OK' if ok else 'REJECTED (fold is unreliable)'}")
sys.exit(0 if ok else 1)
PY
