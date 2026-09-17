#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -J mmsmoke2
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/mmgbsa_smoke2_%j.log
# ⚠ EVERYTHING WRITES TO /bigdata. /tmp and /scratch are NODE-LOCAL -- the first
# attempt built the mbondi3 topologies on the login node and the compute node
# could not see them.
set -uo pipefail
module load amber/22_mpi_cuda >/dev/null 2>&1
P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
S=$P/data/panel_mmgbsa/sys/3phenylphenol_HIT
RES=PPH
TOP=$S/system.prmtop

# 1. strip solvent AND box with cpptraj (no python, no numpy)
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
done

# 2. igb=8 (GBneck2) REQUIRES mbondi3; tleap gave us plain mbondi. ante-MMPBSA.py
#    used to do this via --radii=mbondi3. ParmEd 4.3 can set radii safely because
#    it does NOT slice here -- slicing is what produced the incompatible LJ table.
for f in complex receptor ligand; do
  [[ -s $S/${f}_r3.prmtop ]] && continue
  env -u PYTHONPATH /bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/python - <<PY
import parmed as pmd
from parmed.tools import changeRadii
p = pmd.load_file("$S/$f.prmtop")
changeRadii(p, "mbondi3").execute()
p.save("$S/${f}_r3.prmtop", overwrite=True)
print("$f", len(p.atoms), "atoms ->", p.parm_data["RADIUS_SET"][0])
PY
done

# 3. assert amber's own ParmEd can read them, before spending the calculation
/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/bin/python -c "
from parmed.amber.readparm import LoadParm
n={}
for f in ('complex','receptor','ligand'):
    p=LoadParm('$S/%s_r3.prmtop'%f); n[f]=len(p.atoms)
    print(f'{f:9} {n[f]:>5} atoms  radii={p.parm_data[\"RADIUS_SET\"][0]}')
assert n['receptor']+n['ligand']==n['complex'], n
print('ATOM COUNTS CONSISTENT')
" || { echo '!! amber ParmEd cannot read the mbondi3 files'; exit 3; }

# 4. strip the trajectory too -- ParmEd 3.4.1 cannot read ANY boxed topology
cd $S/rep0 || exit 1
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
[[ -s prod_dry.nc ]] || { echo "!! no prod_dry.nc"; exit 3; }

cat > mmpbsa.in <<MM
smoke2
&general
   startframe=1, endframe=999999, interval=50, verbose=2, keep_files=1,
/
&gb
   igb=8, saltcon=0.150,
/
MM
rm -f mmgbsa.dat
MMPBSA.py -O -i mmpbsa.in -o mmgbsa.dat \
  -cp $S/complex_r3.prmtop -rp $S/receptor_r3.prmtop -lp $S/ligand_r3.prmtop \
  -y prod_dry.nc > mmpbsa.log 2>&1 || { echo "!! MMPBSA FAILED"; tail -12 mmpbsa.log;
     echo "--- gb mdout ---"; tail -15 _MMPBSA_complex_gb.mdout.0 2>/dev/null; exit 4; }
echo "=== RESULT ==="
sed -n '/Differences/,+14p' mmgbsa.dat | head -18
