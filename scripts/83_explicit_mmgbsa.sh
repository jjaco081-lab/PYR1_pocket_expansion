#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 8:00:00
#SBATCH -J exp_gb
#SBATCH -a 0-5
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/exp_gb_%a.log
#
# 83_explicit_mmgbsa.sh -- score the explicit-solvent trajectories with the SAME
# MM-GBSA settings used on the 46 implicit runs.
#
# Exactly one variable differs from README 37/40: the solvent the sampling happened
# in. Scoring is still igb=8, saltcon=0.15, mbondi3, on water-stripped frames --
# the standard MM-PBSA protocol, and the only way these numbers stay commensurable
# with the implicit set. Changing the scorer at the same time would confound the
# comparison this whole exercise exists to make.
#
# Waters and ions are stripped by ante-MMPBSA rather than kept: MM-GBSA has no
# consistent way to assign an explicit water to "receptor" or "ligand", and keeping
# a shell of them would reintroduce, per-frame, the same ligand-specific structural
# advantage that the no-crystal-waters decision in script 80 exists to avoid.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22 >/dev/null 2>&1

UNITS=(aba_WT_s0 aba_WT_s1 aba_WT_s2 mandi_WT_s0 mandi_WT_s1 mandi_WT_s2)
U=${UNITS[${SLURM_ARRAY_TASK_ID:-0}]}
D=$ROOT/data/mmgbsa_explicit/$U
cd "$D" || exit 1
echo "=== $U $(date -Is)"

# dry complex + receptor + ligand topologies, derived from the solvated one
if [[ ! -f com.prmtop ]]; then
    ante-MMPBSA.py -p system.prmtop -c com.prmtop -s ':WAT,K+,Cl-' \
                   -r rec.prmtop -l lig.prmtop -n ':LIG' --radii=mbondi3 \
                   > ante.log 2>&1 || { echo "ante-MMPBSA FAILED"; tail -5 ante.log; exit 2; }
fi
for f in com.prmtop rec.prmtop lig.prmtop; do
    [[ -s $f ]] || { echo "missing $f"; exit 2; }
done

# assert the ligand really is in lig.prmtop and is the expected size, because a
# silently-empty ligand topology would still produce a plausible-looking number
NLIG=$(python3 -c "
import re
t=open('lig.prmtop').read()
m=re.search(r'%FLAG POINTERS.*?%FORMAT\(\S+\)\s+(.*?)%FLAG', t, re.S)
print(int(m.group(1).split()[0]))")
EXP=$([[ $U == aba* ]] && echo 38 || echo 51)
[[ $NLIG -eq $EXP ]] || { echo "ERROR: lig.prmtop has $NLIG atoms, expected $EXP"; exit 2; }
echo "  ligand topology: $NLIG atoms  OK"

cat > mmgbsa.in <<IN
explicit-solvent frames, water-stripped, scored exactly as the implicit runs were
 &general
  startframe=1, endframe=99999, interval=5, verbose=2, keep_files=0,
 /
 &gb
  igb=8, saltcon=0.15,
 /
IN

MMPBSA.py -O -i mmgbsa.in -sp system.prmtop -cp com.prmtop -rp rec.prmtop -lp lig.prmtop \
          -y prod_seg*.nc -o mmgbsa_exp.dat -eo mmgbsa_exp_frames.csv \
          > mmpbsa.log 2>&1 || { echo "MMPBSA FAILED"; tail -20 mmpbsa.log; exit 2; }
grep -A3 "DELTA TOTAL" mmgbsa_exp.dat | head -4
echo "=== $U done $(date -Is)"
