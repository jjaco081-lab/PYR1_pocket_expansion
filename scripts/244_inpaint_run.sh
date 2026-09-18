#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu05,gpu11,gpu13,gpu14
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 3:58:00
#SBATCH -J inpaint
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/inpaint_%j.log
# Inpaint the 5 seams of the PYR1 / Api g 1 chimera with RFd3, ABA IN THE POCKET.
# The chimera is one chain of 145 residues (102 PYR1, 43 donor) whose two parents'
# backbones do not meet at 5 places; spans are 4.6-10.4 A, needing 2-7 residues each.
set -uo pipefail
module load gcc/12.2.0 2>/dev/null; module load cuda/12.1 2>/dev/null
source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda deactivate 2>/dev/null
conda activate /bigdata/cutlerlab/jjaco081/conda_envs/foundry
R=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd $R || exit 1
IN=$R/results/chimera/inpaint_2bk0A00.pdb
CONTIG="A1-19,2-6,A20-41,2-6,A42-78,2-6,A79-96,3-7,A97-104,2-6,A105-145"
OUTD=$R/results/chimera/inpainted_2bk0A00
mkdir -p $OUTD
n=0
for i in $(seq 0 23); do
  d=$OUTD/d$(printf "%03d" $i)
  [[ -s $d/_0_model_0.cif.gz || -s $d/_0_model_0.cif ]] && { n=$((n+1)); continue; }
  mkdir -p $d
  /bigdata/cutlerlab/jjaco081/conda_envs/foundry/bin/rfd3 \
    inputs=null out_dir=$d n_batches=1 \
    +specification.input=$IN \
    +specification.contig=\'"$CONTIG"\' \
    +specification.dialect=2 \
    +specification.ligand=\'A8S\' \
    seed=$((770000 + i)) > $d/rfd3.log 2>&1
  if compgen -G "$d/*.cif*" > /dev/null; then n=$((n+1)); echo "  d$i OK";
  else echo "  d$i NO OUTPUT"; tail -4 $d/rfd3.log; fi
done
echo "RESULT: $n of 24 inpainted designs produced"
test "$n" -ge 12
