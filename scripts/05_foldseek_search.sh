#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -J fs_pyr1
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/05_foldseek.log
#
# 05_foldseek_search.sh -- structural homolog search for PYR1 chain A.
#
# Searches PYR1 (3QN1 chain A, protein only, from 00_setup_inputs.py) against
# the foldseek CATH50 database. CATH50 is a 50%-sequence-identity-clustered
# representative set of CATH domains including AlphaFold models; the
# helix-grip/SRPBCC superfamily is CATH 3.30.530.20.
#
# MEMORY: peak RSS ~1.15 GB. The HPCC login node enforces a 1 GB per-user
# cgroup cap, so running this interactively is silently OOM-killed with
# "Kmer matching step died". It must go through SLURM.
#
# STORAGE: /scratch is node-local xfs and is NOT visible to compute nodes;
# all paths here are under /bigdata (GPFS, shared).
#
# foldseek commit 941cd33ff0771cd2e3f144e3293e22a2b87e9fda
# CATH50 database downloaded 2024-01-25
set -euo pipefail
P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
FS=/bigdata/cutlerlab/jjaco081/tools/foldseek/foldseek/bin/foldseek
DB=/bigdata/cutlerlab/jjaco081/tools/foldseek/databases/CATH50

mkdir -p "$P/results/foldseek"
rm -rf "$P/results/foldseek/tmp"

"$FS" easy-search \
  "$P/data/pyr1_A.pdb" \
  "$DB" \
  "$P/results/foldseek/fs_hits.tsv" \
  "$P/results/foldseek/tmp" \
  --format-output "query,target,fident,alntmscore,evalue,qstart,qend,tstart,tend,qaln,taln" \
  -e 10 --max-seqs 4000 --threads 8

echo "DONE rc=$? hits=$(wc -l < "$P/results/foldseek/fs_hits.tsv")"
