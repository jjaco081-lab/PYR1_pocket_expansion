#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -J pe_mdbench
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/md_bench.log
# Measure real pmemd.cuda throughput on the smallest and largest systems so the
# production campaign can be sized from data instead of estimates.
set -uo pipefail
module load amber/22_mpi_cuda >/dev/null 2>&1
P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
MD=$P/data/md
nvidia-smi --query-gpu=name --format=csv,noheader
for sys in S1_apo_open S4_ternary; do
  D=$MD/$sys/bench; mkdir -p "$D"; cd "$D" || exit 1
  TOP=$MD/$sys/system.prmtop
  cat > bmin.in <<IN
 &cntrl
  imin=1, maxcyc=500, ncyc=250, ntb=1, cut=10.0, ntr=0,
 /
IN
  cat > bmd.in <<IN
 &cntrl
  imin=0, irest=0, ntx=1, nstlim=10000, dt=0.002,
  ntc=2, ntf=2, cut=10.0, ntb=2, ntp=1, barostat=2,
  ntt=3, gamma_ln=2.0, ig=-1, tempi=300.0, temp0=300.0,
  ntpr=1000, ntwx=0, ntwr=10000,
 /
IN
  pmemd.cuda -O -i bmin.in -p "$TOP" -c $MD/$sys/system.inpcrd -o bmin.out -r bmin.rst7 >/dev/null 2>&1
  pmemd.cuda -O -i bmd.in  -p "$TOP" -c bmin.rst7 -o bmd.out -r bmd.rst7 >/dev/null 2>&1
  natom=$(awk '/NATOM/{getline; print $1; exit}' bmd.out 2>/dev/null)
  ns=$(grep -A3 "Average timings for all steps" bmd.out 2>/dev/null | grep -oP 'ns/day\s*=\s*\K[0-9.]+' | head -1)
  [ -z "$ns" ] && ns=$(grep -oP 'ns/day\s*=\s*\K[0-9.]+' bmd.out | tail -1)
  echo "BENCH $sys atoms=${natom:-?} ns_per_day=${ns:-FAILED}"
done
echo BENCH_DONE
