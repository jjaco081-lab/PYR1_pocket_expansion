#!/bin/bash -l
#SBATCH -J cavrender
#SBATCH -p cutlerlab
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/cavity_render_%j.log
# Ray-traced PNGs, no GLX needed at all (DISPLAY unset). Jannis's X server does
# not advertise a GLX FBConfig Qt can use, so interactive PyMOL over SSH cannot
# initialise a context -- that, not the PyMOL version, is the black-viewport bug.
set -uo pipefail
unset DISPLAY
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/results/cavity_check
PM=/rhome/jjaco081/bigdata/envs/pymol_env3/bin/pymol
n=0
for f in *.pml; do
  b=${f%.pml}
  sed "/^zoom /a set ray_opaque_background, 1" "$f" > /tmp/r_$b.pml
  cat >> /tmp/r_$b.pml <<PML
ray 1100, 850
png $PWD/${b}_render.png, dpi=130
PML
  timeout 600 $PM -cq /tmp/r_$b.pml > /tmp/r_$b.log 2>&1
  if [[ -s ${b}_render.png ]]; then echo "OK   $b"; n=$((n+1));
  else echo "FAIL $b"; tail -4 /tmp/r_$b.log; fi
done
echo "rendered $n images"
test "$n" -ge 5
