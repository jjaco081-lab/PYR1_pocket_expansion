#!/bin/bash
#SBATCH -p cutlerlab
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -J pe_view
#SBATCH -o /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/logs/md_view_%j.log
#
# 41_md_view.sh -- build lightweight, visualisable copies of the MD trajectories.
#
# WHY THIS EXISTS
#   The raw prod.nc files are ~15 GB each: ~28,000 frames of a fully solvated
#   box. They are correct but hostile to a viewer. Worse, production ran with
#   iwrap=1, so molecules are wrapped into the primary box -- the protein
#   appears to teleport across the periodic boundary partway through, and a
#   dimer or complex can appear to fly apart. That is a display artefact of
#   wrapping, not a simulation failure, but it makes visual QC useless.
#
#   This produces, per replicate, a stripped + imaged + aligned trajectory:
#     * waters, Na+ and Cl- removed          (~46,800 atoms -> ~2,900)
#     * autoimage                             (undoes the PBC jumps)
#     * rms fit on CA to the first frame      (removes global tumbling)
#     * every 10th frame kept, i.e. 100 ps    (~2,800 frames, tens of MB)
#   Ligand (ABA) and Mn2+ are deliberately KEPT -- they are the point.
#
# RUNTIME
#   Minutes, not hours. NetCDF supports direct frame access, so a stride of 10
#   reads ~1/10 of each file rather than streaming all 15 GB. The 2 h walltime
#   is slack, not an estimate. The job is idempotent, so a walltime kill costs
#   nothing -- just resubmit.
#
# RUN IT TWICE
#   Safe to run while jobs are still in flight: cpptraj reads the frame count
#   from the NetCDF header, so it sees however much has been flushed, and
#   systems with no prod.nc yet are skipped. Run it now for visual QC on the
#   in-flight replicates -- catching a bad setup before S3/S4 burn three more
#   days of GPU is the whole point -- then run it again at the end to rebuild
#   the real artefacts from the completed trajectories.
#
# OUTPUT
#   data/md/view/<system>_<rep>.prmtop   topology matching the stripped traj
#   data/md/view/<system>_<rep>.nc       the trajectory to load
#   data/md/view/<system>_<rep>_frame0.pdb  first frame, for a quick look
#   data/md/view/<system>_<rep>_rmsd.dat    CA RMSD vs frame 0, ns in col 1
#   data/md/view/<system>_<rep>_rgyr.dat    radius of gyration
#
# HOW TO VIEW (VMD)
#   vmd -parm7 S1_apo_open_rep0.prmtop -netcdf S1_apo_open_rep0.nc
# PyMOL
#   load S1_apo_open_rep0.prmtop
#   load_traj S1_apo_open_rep0.nc
set -uo pipefail
module load amber/22 >/dev/null 2>&1

P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
MD=$P/data/md
OUT=$MD/view
mkdir -p "$OUT"

STRIDE=10          # keep every 10th frame -> 100 ps spacing
SOLVENT=":WAT,Na+,Cl-"

for sys in S1_apo_open S2_holo_closed S3_apo_dimer S4_ternary; do
    TOP=$MD/$sys/system.prmtop
    [[ -s "$TOP" ]] || continue
    for rep in 0 1 2; do
        NC=$MD/$sys/rep$rep/prod.nc
        [[ -s "$NC" ]] || continue
        tag=${sys}_rep${rep}
        echo "=== $tag ==="

        # one-time stripped topology
        if [[ ! -s $OUT/$tag.prmtop ]]; then
            cpptraj -p "$TOP" <<EOF
parmstrip $SOLVENT
parmwrite out $OUT/$tag.prmtop
EOF
        fi

        cpptraj -p "$TOP" <<EOF
trajin $NC 1 last $STRIDE
autoimage
strip $SOLVENT
rms first @CA out $OUT/${tag}_rmsd.dat time 0.1
radgyr @CA out $OUT/${tag}_rgyr.dat time 0.1
trajout $OUT/$tag.nc netcdf
trajout $OUT/${tag}_frame0.pdb pdb onlyframes 1
go
EOF
        ls -la $OUT/$tag.nc 2>/dev/null | awk '{print "  wrote "$9" "$5" bytes"}'
    done
done

echo
echo "=== summary ==="
ls -la "$OUT" | awk '{printf "  %-42s %s\n", $9, $5}'
