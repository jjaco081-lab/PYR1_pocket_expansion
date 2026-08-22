#!/bin/bash
# 82_explicit_rmsd.sh -- did the ligand hold its crystallographic pose in real water?
#
# THE TEST, stated before the numbers (README 41): in implicit solvent ABA slid
# 1.09 -> 3.18 A between halves and ended 4.13 A from its start, while its K59 salt
# bridge stayed intact 100% of the time. Mandipropamid, being neutral, held at 1.7 A.
# If explicit water is the fix, ABA should now look like mandipropamid did.
# If ABA still slides 2-4 A, explicit solvent is NOT the fix and MM-GBSA should be
# retired for this ligand rather than tuned further.
#
# The fit is done FIRST and the ligand measured `nofit` afterwards. Skipping the fit
# was an error made once already: with ntb=0 there is no box to anchor the molecule,
# so a bare `nofit` RMSD silently measures global tumbling instead of internal
# motion. Here ntb=2, but the fit is kept explicit so the number means the same
# thing in both solvent models.
set -uo pipefail
ROOT=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
cd "$ROOT" || exit 1
module load amber/22 >/dev/null 2>&1

for D in data/mmgbsa_explicit/*/; do
    U=$(basename "$D")
    cpptraj -p "$D/system.prmtop" >/dev/null 2>&1 <<CPT
$(for f in "$D"/prod_seg*.nc; do echo "trajin $f"; done)
reference $D/eq2.rst7 name REF
autoimage
rms fitprot ref REF :1-190@CA,C,N out /dev/null
rms lig   ref REF :LIG&!@H= nofit out $D/rmsd_lig.dat
rms gate  ref REF :85-89&!@H= nofit out $D/rmsd_gate.dat
$(if [[ $U == aba* ]]; then
    # ONLY meaningful for ABA. O3/O4 are its deprotonated carboxylate; the same
    # atom names in mandipropamid (3UZ) are unrelated oxygens, so measuring a
    # "salt bridge" there reports a number for a bond that does not exist.
    echo "distance sb1 :59@NZ :LIG@O3 out $D/sb.dat"
    echo "distance sb2 :59@NZ :LIG@O4 out $D/sb.dat"
  fi)
run
CPT
done
python3 - <<'PY'
import glob, os
print("="*78)
print("EXPLICIT SOLVENT, 10 ns x 3 seeds -- did the pose hold?")
print("="*78)
print(f"  {'system':<14}{'lig 1st':>9}{'lig 2nd':>9}{'lig max':>9}{'gate 2nd':>10}{'SB<4A':>8}")
for d in sorted(glob.glob("data/mmgbsa_explicit/*/")):
    u = os.path.basename(d.rstrip("/"))
    v = [float(l.split()[1]) for l in open(d+"rmsd_lig.dat") if not l.startswith("#")]
    g = [float(l.split()[1]) for l in open(d+"rmsd_gate.dat") if not l.startswith("#")]
    h = len(v)//2
    sb = "n/a" if "mandi" in u else ""
    if os.path.exists(d+"sb.dat") and "mandi" not in u:
        rows = [l.split() for l in open(d+"sb.dat") if not l.startswith("#")]
        try:
            m = [min(float(r[1]), float(r[2])) for r in rows]
            sb = f"{100.0*sum(1 for x in m if x<4.0)/len(m):.0f}%"
        except Exception:
            sb = "-"
    print(f"  {u:<14}{sum(v[:h])/h:>9.2f}{sum(v[h:])/(len(v)-h):>9.2f}{max(v):>9.2f}"
          f"{sum(g[h:])/(len(g)-h):>10.2f}{sb:>8}")
print()
print("  implicit-solvent reference (README 41):")
print("    aba_WT       1.09     3.18     4.88      (gate 2.07)   100%")
print("    mandi_WT     1.45     1.70     2.45")
PY
