#!/usr/bin/env python
"""
02_rigid_truncation_scan.py -- rank pocket residues by the volume they occlude.

For each first-shell residue, delete its side chain beyond CB (an in-silico
alanine) WITHOUT moving anything else, and re-measure the enclosed cavity.
The delta is an UPPER BOUND on the volume that residue is withholding: no
neighbour is allowed to relax into the new space. Script 06 supplies the
matching post-relax number; the gap between the two is the infilling.

Reading the output
------------------
POSITIVE delta -> the side chain occludes cavity volume (a restrictor).
NEGATIVE delta -> deleting it breaches the cavity to bulk solvent, so the
                  connected component leaks and buriedness drops below the
                  cut-off. These are SEALING residues. They must be left alone:
                  removing them opens the pocket to solvent rather than
                  enlarging it. F61, F159 and L87 behave this way in PYR1.

Run with the esmfold2 env python.
"""
import os, sys, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_cavity import cavity_volume, atoms_from_pose_pdb  # noqa: E402

PDB = os.path.join(ROOT, "data", "pyr1_A.pdb")
LIG = os.path.join(ROOT, "data", "aba_xtal.pdb")

# first shell from 01_pocket_contacts.py (PYR1 chain A, within 5.0 A of ABA)
SHELL = [59, 159, 89, 108, 61, 110, 163, 120, 83, 92, 117, 167, 88, 115, 87,
         141, 94, 164, 91, 160, 122]
SKIP_NATIVE = {"ALA", "GLY", "PRO"}   # nothing to truncate

lig_xyz, _ = atoms_from_pose_pdb(LIG, exclude_resnames=())
REF = lig_xyz.mean(0)

resname = {}
for line in open(PDB):
    if line.startswith("ATOM"):
        resname[int(line[22:26])] = line[17:20].strip()

xyz, el = atoms_from_pose_pdb(PDB, chain="A")
wt, _ = cavity_volume(xyz, el, REF)
print(f"WT closed-state cavity: {wt:.1f} A^3\n")

rows = []
for rn in SHELL:
    if resname.get(rn) in SKIP_NATIVE:
        continue
    x2, e2 = atoms_from_pose_pdb(PDB, chain="A", truncate_to_ala={rn})
    v, _ = cavity_volume(x2, e2, REF)
    rows.append(dict(resnum=rn, wt_aa=resname[rn], cavity_A3=round(v, 1),
                     delta_A3=round(v - wt, 1)))

rows.sort(key=lambda r: -r["delta_A3"])
print(f"{'residue':<10}{'cavity':>9}{'delta':>9}")
for r in rows:
    print(f"{r['wt_aa']}{r['resnum']:<7}{r['cavity_A3']:>9.1f}{r['delta_A3']:>+9.1f}")

print("\n=== combinations ===")
combos = [(59, 108), (59, 159), (108, 159), (59, 108, 159), (59, 108, 159, 61),
          (79, 94), (108, 79), (108, 79, 94), (59, 108, 79, 94)]
crows = []
for c in combos:
    x2, e2 = atoms_from_pose_pdb(PDB, chain="A", truncate_to_ala=set(c))
    v, _ = cavity_volume(x2, e2, REF)
    label = "/".join(f"{resname[i]}{i}A" for i in c)
    crows.append(dict(combo=label, cavity_A3=round(v, 1), delta_A3=round(v - wt, 1)))
    print(f"  {label:<40}{v:>8.1f}{v-wt:>+9.1f}")

outdir = os.path.join(ROOT, "results")
os.makedirs(outdir, exist_ok=True)
with open(os.path.join(outdir, "02_rigid_truncation.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["wt_cavity_A3", round(wt, 1)])
    w.writerow([])
    w.writerow(["resnum", "wt_aa", "cavity_A3", "delta_A3"])
    for r in rows:
        w.writerow([r["resnum"], r["wt_aa"], r["cavity_A3"], r["delta_A3"]])
    w.writerow([])
    w.writerow(["combo", "cavity_A3", "delta_A3"])
    for r in crows:
        w.writerow([r["combo"], r["cavity_A3"], r["delta_A3"]])
print(f"\nwrote {outdir}/02_rigid_truncation.csv")
