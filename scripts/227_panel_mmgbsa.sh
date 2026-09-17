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
# FOUR chained defects had to be cleared here; all four are load-bearing.
#
# 1. ante-MMPBSA.py shells out to amber/22's ParmEd 3.4.1, broken against numpy
#    2.x: `np.array(box, copy=False)` raises the moment it reads a topology that
#    HAS A BOX. Conda ParmEd 4.3 loads fine but SLICING with it writes an LJ
#    table 3.4.1 rejects ("LENNARD_JONES_ACOEF has 378 elements; expected 276"),
#    so swapping interpreters does not fix it. The same ParmEd/numpy trap is
#    already recorded for the TI pipeline; it cost all 30 completed replicates of
#    the first submission.  -> strip with cpptraj (pure C++) + `parmbox nobox`.
# 2. -sp is the documented way to point MMPBSA at a solvated trajectory, but it
#    hands over the BOXED topology and re-triggers (1).
#    -> strip the TRAJECTORY too, so nothing boxed ever reaches MMPBSA.
# 3. igb=8 (GBneck2) REQUIRES mbondi3; tleap gives plain mbondi and
#    mmpbsa_py_energy dies with no useful message. ante-MMPBSA.py had been doing
#    this silently via --radii=mbondi3.  -> ParmEd 4.3 changeRadii, which is safe
#    precisely because it does NOT slice.
# 4. /tmp AND /scratch ARE NODE-LOCAL. Building these on the login node left the
#    compute node unable to see them. Everything here writes to /bigdata.
RES=$(/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python -c "
import json
f=[x for x in json.load(open('$PM/frames.json'))
   if '$TAG'.startswith(x['ligand'].replace('-',''))]
print(f[0]['resname'] if f else '')")
[[ -n "$RES" ]] || { echo "!! cannot resolve ligand resname for $TAG"; exit 3; }

for spec in "complex:!(:WAT,K+,Cl-)" "receptor:!(:WAT,K+,Cl-,$RES)" "ligand::$RES"; do
    NM=${spec%%:*}; SEL=${spec#*:}
    [[ -s $S/$NM.prmtop ]] && continue
    cpptraj > $S/$NM.cpptraj.log 2>&1 <<CPP
parm $TOP
parmstrip !($SEL)
parmbox nobox
parmwrite out $S/$NM.prmtop
go
CPP
    [[ -s $S/$NM.prmtop ]] || { echo "!! cpptraj failed on $NM"; exit 3; }
done
for f in complex receptor ligand; do
  [[ -s $S/${f}_r3.prmtop ]] && continue
  env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/python - <<RADII
import parmed as pmd
from parmed.tools import changeRadii
q = pmd.load_file("$S/$f.prmtop")
changeRadii(q, "mbondi3").execute()
q.save("$S/${f}_r3.prmtop", overwrite=True)
RADII
done
/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/bin/python -c "
from parmed.amber.readparm import LoadParm
n={}
for f in ('complex','receptor','ligand'):
    q=LoadParm('$S/%s_r3.prmtop'%f); n[f]=len(q.atoms)
    assert 'mbondi3' in q.parm_data['RADIUS_SET'][0], (f, q.parm_data['RADIUS_SET'])
assert n['receptor']+n['ligand']==n['complex'], n
print('topologies OK:', n)
" || { echo "!! topology check FAILED for $TAG"; exit 3; }

if [[ ! -s prod_dry.nc ]]; then
cpptraj > strip_traj.log 2>&1 <<CPP
parm $TOP
trajin prod.nc
strip :WAT,K+,Cl-
box nobox
trajout prod_dry.nc netcdf
go
CPP
fi
[[ -s prod_dry.nc ]] || { echo "!! no prod_dry.nc for $TAG $REP"; exit 3; }

if [[ ! -s mmgbsa.dat ]]; then
    cat > mmpbsa.in <<MM
MM-GBSA on $TAG $REP
&general
   startframe=1, endframe=999999, interval=10, verbose=2, keep_files=0,
/
&gb
   igb=8, saltcon=0.150,
/
MM
    MMPBSA.py -O -i mmpbsa.in -o mmgbsa.dat \
        -cp $S/complex_r3.prmtop -rp $S/receptor_r3.prmtop -lp $S/ligand_r3.prmtop \
        -y prod_dry.nc > mmpbsa.log 2>&1 || { echo "!! MMPBSA failed";
            tail -25 mmpbsa.log; exit 4; }
fi

# assert on the RESULT, not the exit code. The value sits on the SAME line as
# the label: "DELTA TOTAL  -23.3875  2.4524  0.7755" -> field 3.
G=$(awk '/^DELTA TOTAL/{print $3; exit}' mmgbsa.dat 2>/dev/null)
if [[ -z "$G" ]]; then echo "!! no DELTA TOTAL in mmgbsa.dat"; exit 5; fi
echo "RESULT $TAG $REP  dG_bind = $G kcal/mol"
