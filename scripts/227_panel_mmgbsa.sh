#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14   # amber22 pmemd.cuda has no kernels for k80/p100/h100/blackwell
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 1:58:00
#SBATCH -J panelmm
#SBATCH -a 0-63%4
#SBATCH --open-mode=append
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/panelmm_%a.log
#
# 227 -- MM-GBSA on the four control ligands, WT vs known hit, n = 8.
#
# 8 systems (4 ligands x {WT, HIT}) x 8 replicates = 64 runs. Array index maps
# idx/8 -> system, idx%8 -> replicate.
#
# WHY n = 8 AND WHAT IT CAN AND CANNOT RESOLVE (pre-registered).
# §42 retracted MM-GBSA here on precision grounds: ddG carries +/-3.85 kcal/mol
# against a target spacing of 0.01-1.8, giving n = 60/variant. Those grounds do
# not apply at THIS effect size. The control ladder spans 100x in min_conc =
# 2.73 kcal/mol, and n = (3.85/(delta/2))^2 = 8.
# ⚠ So the registered claim is ORDERING OF THE EXTREMES ONLY. n = 8 cannot
# separate adjacent 10x rungs; that needs n ~ 32. Do not read rung-to-rung.
#
# ⚠ EXPLICIT SOLVENT IS NON-NEGOTIABLE (§41/§42): zero explicit waters let ABA
# slide 4 A, and the 250 ps "convergence" was drift AWAY from the answer.
# Solvation came from lib_solvate.sh (ff19SB/OPC, 0.15 M KCl), same as every
# other MD in this project.
#
# ⚠ POSE DEGENERACY IS NOT IN THE ERROR BAR. 209 found eugenol accepts 117 of
# 180 anchored placements at zero clashes, versus anthrone's 6 of 36. Each run
# starts from ONE of those poses, so the n=8 replicate spread measures sampling
# around a pose, not the choice of pose. Same class as §93c, where seed spread
# understated true error ~3x. Report both.
#
# GO/NO-GO: if MM-GBSA cannot order anthrone (1 uM) above eugenol (100 uM),
# the arm closes.
set -uo pipefail
module load amber/22_mpi_cuda >/dev/null 2>&1

P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
PM=$P/data/panel_mmgbsa
source "$P/scripts/lib_mdinputs.sh"
md_settings_check

PROD_NS=5
TARGET_PS=$(( EQUIL_PS + PROD_NS * 1000 ))

mapfile -t SYS < <(ls -d $PM/sys/*/ | sed 's#/$##' | sort)
NS=${#SYS[@]}
if [[ $NS -ne 8 ]]; then echo "expected 8 systems, found $NS"; exit 1; fi
S=${SYS[$(( SLURM_ARRAY_TASK_ID / 8 ))]}
REP=rep$(( SLURM_ARRAY_TASK_ID % 8 ))
TAG=$(basename "$S")
TOP=$S/system.prmtop
CRD=$S/system.inpcrd
[[ -s $TOP ]] || { echo "no topology for $TAG"; exit 1; }

D=$S/$REP
mkdir -p "$D"; cd "$D" || exit 1
echo "=== $TAG $REP  ${PROD_NS} ns  $(date -Is)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

write_md_inputs

run_stage () {
    local nm=$1 prev=$2
    if [[ -s ${nm}.rst7 ]]; then echo "    $nm done, skipping"; return 0; fi
    echo "    running $nm ..."
    local ref=()
    [[ "${3:-}" == "ref" ]] && ref=(-ref "$prev")
    pmemd.cuda -O -i ${nm}.in -p "$TOP" -c "$prev" \
               -o ${nm}.out -r ${nm}.rst7 -x ${nm}.nc -inf ${nm}.info "${ref[@]}" \
        || { echo "    !! $nm FAILED"; return 2; }
}

run_stage min1 "$CRD" ref || exit 2
run_stage min2 min1.rst7  || exit 2
run_stage heat min2.rst7 ref || exit 2
run_stage eq1  heat.rst7 ref || exit 2
run_stage eq2  eq1.rst7      || exit 2

if [[ ! -s prod.rst7 ]]; then
    # 2 fs timestep (HMR is OFF per lib_mdinputs.sh), so ns -> steps is x500000
    write_prod_input $(( PROD_NS * 500000 ))
    pmemd.cuda -O -i prod.in -p "$TOP" -c eq2.rst7 \
        -o prod.out -r prod.rst7 -x prod.nc -inf prod.info \
        || { echo "!! prod FAILED"; exit 2; }
fi

# ---- MM-GBSA on the production trajectory -------------------------------
# strip water/ions once, then MMPBSA.py with igb=8 (GBneck2).
if [[ ! -s mmgbsa.dat ]]; then
    STRIP=$S/stripped.prmtop
    if [[ ! -s $STRIP ]]; then
        cpptraj -p "$TOP" <<CPP > strip.log 2>&1
parmstrip :WAT,K+,Cl-
parmwrite out $STRIP
go
CPP
    fi
    LIG=$(grep -oE '^(ANT|PPH|CXL|EUG)' <<< "$(basename $S)" || true)
    RES=$(/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python -c "
import json,sys
f=[x for x in json.load(open('$PM/frames.json')) if '$TAG'.startswith(x['ligand'].replace('-',''))]
print(f[0]['resname'] if f else '')")
    [[ -n "$RES" ]] || { echo "!! cannot resolve ligand resname for $TAG"; exit 3; }
    cat > mmpbsa.in <<MM
MM-GBSA on $TAG $REP
&general
   startframe=1, endframe=999999, interval=10, verbose=2,
   keep_files=0,
/
&gb
   igb=8, saltcon=0.150,
/
MM
    ante-MMPBSA.py -p "$TOP" -c complex.prmtop -r receptor.prmtop -l ligand.prmtop \
        -s ':WAT,K+,Cl-' -n ":$RES" --radii=mbondi3 > ante.log 2>&1
    MMPBSA.py -O -i mmpbsa.in -o mmgbsa.dat -sp "$TOP" \
        -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop \
        -y prod.nc > mmpbsa.log 2>&1 || { echo "!! MMPBSA failed"; tail -20 mmpbsa.log; exit 4; }
fi

# assert on the RESULT, not the exit code
G=$(grep -A3 "DELTA TOTAL" mmgbsa.dat 2>/dev/null | head -2 | tail -1 | awk '{print $2}')
if [[ -z "$G" ]]; then echo "!! no DELTA TOTAL in mmgbsa.dat"; exit 5; fi
echo "RESULT $TAG $REP  dG_bind = $G kcal/mol"
