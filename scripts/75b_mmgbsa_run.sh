#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 06:00:00
#SBATCH -J mmgbsa
#SBATCH -a 0-13
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/mmgbsa_%a.log
#
# 75b_mmgbsa_run.sh -- MM-GBSA rescoring of the stage-1 variants (README 37).
#
# 14 tasks = 7 variants x 2 ligand arms. Each builds 8 independently repacked
# structures into one topology, minimises them in implicit solvent, and runs
# MMPBSA.py over the 8 frames. The seed-to-seed spread IS the error bar -- a single
# MM-GBSA number for a charged residue is not interpretable alone.
#
# WHY ff14SB AND NOT ff19SB
#   Our explicit-solvent MD uses ff19SB, which was parameterised with OPC water.
#   ff19SB + generalised Born is not a validated combination; ff14SB + igb=8 is the
#   well-trodden one for MM-GBSA. This is a rescoring, deliberately separate from the
#   MD protocol, so it uses the force field its own method was built for. Recorded
#   here rather than silently mixed.
#
# RADII
#   igb=8 requires mbondi3; the tleap default is not compatible and sander will not
#   run at all with it. Set in every leap script below.
#
# HISTIDINE
#   tleap's default (HIE) is used for every variant identically. The variants differ
#   only at one side chain, so the histidine network cancels in the comparison; what
#   would NOT cancel is assigning it per structure and letting it vary between seeds.
set -uo pipefail
module load amber/22

P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
D=$P/data/mmgbsa
PAR=$D/params
ARMS=(mandi aba)
VARIANTS=(WT K59R K59Q K59N V81I F108A F159L)

ARM=${ARMS[$(( SLURM_ARRAY_TASK_ID / 7 ))]}
VAR=${VARIANTS[$(( SLURM_ARRAY_TASK_ID % 7 ))]}
case $ARM in
  mandi) LIGCODE=3UZ ;;
  aba)   LIGCODE=A8S ;;
esac
W=$D/${ARM}_${VAR}
echo "=== $(date) task $SLURM_ARRAY_TASK_ID -> $ARM / $VAR (ligand $LIGCODE) ==="
[[ -d $W ]] || { echo "no build dir $W -- run 75_mmgbsa_build.py first"; exit 1; }
cd "$W" || exit 1

# ---- one topology per variant, then one minimisation per seed ----
NSEED=$(ls seed*_protein.pdb 2>/dev/null | wc -l)
[[ $NSEED -ge 2 ]] || { echo "only $NSEED seeds; need >=2 for an error bar"; exit 1; }
echo "seeds: $NSEED"

cat > min.in <<CIN
minimise in implicit solvent
 &cntrl
  imin=1, maxcyc=2000, ncyc=1000,
  ntb=0, igb=8, saltcon=0.15, cut=999.0,
  ntpr=200,
 /
CIN

NAT=""
for s in $(seq 0 $((NSEED-1))); do
    pdb4amber -i seed${s}_protein.pdb -o s${s}_clean.pdb --dry --nohyd > s${s}_p4a.log 2>&1
    cat > s${s}_leap.in <<LIN
source leaprc.protein.ff14SB
source leaprc.gaff2
# igb=8 (and 7) require mbondi3 radii. tleap's default assigns hydrogens 0.8 A,
# outside the 1.0-2.0 A range those models accept, and sander refuses to start:
#   "Atom 38 has radius 0.8 outside of allowed range ... Regenerate prmtop with
#    bondi radii." Set here rather than after the fact, because the radii are baked
# into the prmtop and MMPBSA reads the same ones.
set default PBRadii mbondi3
loadamberparams $PAR/${LIGCODE}.frcmod
LIG = loadmol2 $PAR/${LIGCODE}.mol2
prot = loadpdb s${s}_clean.pdb
com = combine { prot LIG }
saveamberparm com s${s}.prmtop s${s}.inpcrd
quit
LIN
    tleap -f s${s}_leap.in > s${s}_leap.log 2>&1
    [[ -s s${s}.prmtop ]] || { echo "  tleap FAILED for seed $s"; tail -5 s${s}_leap.log; exit 1; }
    n=$(awk '/%FLAG POINTERS/{getline;getline;print $1;exit}' s${s}.prmtop)
    if [[ -z "$NAT" ]]; then NAT=$n; else
        [[ "$n" == "$NAT" ]] || { echo "  ATOM COUNT DIFFERS between seeds ($n vs $NAT) -- topologies are not interchangeable"; exit 1; }
    fi
    sander -O -i min.in -p s${s}.prmtop -c s${s}.inpcrd -o s${s}_min.out -r s${s}_min.rst7 >/dev/null 2>&1
    [[ -s s${s}_min.rst7 ]] || { echo "  minimisation FAILED for seed $s"; exit 1; }
done
echo "topology: $NAT atoms, identical across all $NSEED seeds"

# ---- the minimised structures become one trajectory ----
cat > cat.in <<CPP
parm s0.prmtop
$(for s in $(seq 0 $((NSEED-1))); do echo "trajin s${s}_min.rst7"; done)
trajout ensemble.nc netcdf
run
quit
CPP
cpptraj -i cat.in > cat.log 2>&1
[[ -s ensemble.nc ]] || { echo "  cpptraj failed to build the ensemble"; tail -5 cat.log; exit 1; }

# ---- receptor / ligand topologies, then MM-GBSA ----
# ⚠ NO -c HERE. ante-MMPBSA.py only builds a complex topology when it is given a
# SOLVATED one plus a strip mask; ours is already the unsolvated complex, so asking
# for -c com.prmtop silently produces rec/lig and no complex. s0.prmtop IS the
# complex topology and is passed to MMPBSA as -cp directly.
ante-MMPBSA.py -p s0.prmtop -r rec.prmtop -l lig.prmtop -n ":LIG" > ante.log 2>&1
for f in rec.prmtop lig.prmtop; do
    [[ -s $f ]] || { echo "  ante-MMPBSA failed to write $f"; tail -5 ante.log; exit 1; }
done

cat > mmgbsa.in <<MIN
MM-GBSA over independently repacked, minimised structures
 &general
  startframe=1, endframe=$NSEED, interval=1, verbose=2, keep_files=0,
 /
 &gb
  igb=8, saltcon=0.15,
 /
MIN
MMPBSA.py -O -i mmgbsa.in -o mmgbsa.dat -eo mmgbsa_frames.csv \
          -cp s0.prmtop -rp rec.prmtop -lp lig.prmtop -y ensemble.nc > mmpbsa.log 2>&1
rc=$?
if [[ -s mmgbsa.dat ]]; then
    echo "--- DELTA TOTAL ---"
    grep -A3 "DELTA TOTAL" mmgbsa.dat | head -4
else
    echo "  MMPBSA FAILED"; tail -15 mmpbsa.log
fi
echo "MMGBSA_DONE arm=$ARM var=$VAR rc=$rc"
date
exit $rc
