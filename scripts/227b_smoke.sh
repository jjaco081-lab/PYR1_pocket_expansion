#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -J mmsmoke
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/mmgbsa_smoke_%j.log
set -uo pipefail
module load amber/22_mpi_cuda >/dev/null 2>&1
P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
S=$P/data/panel_mmgbsa/sys/3phenylphenol_HIT
RES=PPH
TOP=$S/system.prmtop
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
for f in complex receptor ligand; do
  [[ -s $S/$f.prmtop ]] || { echo "!! missing $f.prmtop"; exit 3; }
done
cd $S/rep0 || exit 1
cat > mmpbsa.in <<MM
smoke
&general
   startframe=1, endframe=999999, interval=50, verbose=2, keep_files=0,
/
&gb
   igb=8, saltcon=0.150,
/
MM
rm -f mmgbsa.dat
MMPBSA.py -O -i mmpbsa.in -o mmgbsa.dat \
  -cp $S/complex.prmtop -rp $S/receptor.prmtop -lp $S/ligand.prmtop \
  -y prod_dry.nc > mmpbsa.log 2>&1 || { echo "!! MMPBSA FAILED"; tail -25 mmpbsa.log; exit 4; }
echo "=== mmgbsa.dat DELTA block ==="
grep -A6 "DELTA TOTAL\|Differences (Complex - Receptor - Ligand)" mmgbsa.dat | head -20
G=$(grep -A3 "DELTA TOTAL" mmgbsa.dat | head -2 | tail -1 | awk '{print $2}')
echo "RESULT dG_bind = ${G:-NONE} kcal/mol"
test -n "$G"
