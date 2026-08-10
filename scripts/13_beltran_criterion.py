#!/usr/bin/env python
"""
13_beltran_criterion.py -- reproduce the Beltran et al. 2022 library-position
filter and locate R79 within it.

The Beltran supplementary states, verbatim:

    We examined PYR1 bound to ABA (PDB ID 3K90, chain A) and selected all
    residues either:
      - with any atom within 5 A of ABA, or
      - with any atom within 6 A of ABA and a Calpha-Cbeta vector not directed
        away from ABA.
    After manual curation of the resulting list of positions, we removed the
    following residues:
      - H60: Faces away from ABA and involved in PYR1 dimerization
      - F61: Outward facing and involved in PYR1 dimerization
      - R79: Appears to have only second-shell effects on ligand binding
      - P88: Conserved in the gate loop
      - T91: Outward facing
      - H115: Conserved in the latch loop

This script reimplements that geometric filter on 3K90 chain A and asks:

  1. Does R79 PASS the automated geometric filter? (i.e. was it excluded by the
     rule, or by the manual curation step?)
  2. For each of the six manually removed residues, is the stated reason
     supported by an independent measurement -- Calpha-Cbeta orientation,
     cavity lining, and line of sight to the ligand?

The "Calpha-Cbeta vector not directed away from ABA" test is implemented as
cos(theta) > 0, where theta is the angle between the Calpha->Cbeta vector and
the Calpha->(nearest ABA atom) vector. Glycine has no Cbeta and is excluded
from that branch.

Run with the esmfold2 env python.
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CIF = os.path.join(ROOT, "data", "3K90.cif")

REMOVED = {
    60:  "Faces away from ABA and involved in PYR1 dimerization",
    61:  "Outward facing and involved in PYR1 dimerization",
    79:  "Appears to have only second-shell effects on ligand binding",
    88:  "Conserved in the gate loop",
    91:  "Outward facing",
    115: "Conserved in the latch loop",
}
# the 19 positions Beltran actually randomised
DSM19 = [59, 62, 81, 83, 87, 89, 92, 94, 108, 110, 117, 120, 122, 141,
         159, 160, 163, 164, 167]

model = MMCIFParser(QUIET=True).get_structure("3k90", CIF)[0]
chain = model["A"]
res = {r.id[1]: r for r in chain if is_aa(r)}
lig = [r for r in chain if r.get_resname() == "A8S"]
if not lig:
    lig = [r for ch in model for r in ch if r.get_resname() == "A8S"]
lig = lig[0]
L = np.array([a.coord for a in lig if a.element != "H"])
print(f"3K90 chain A: {len(res)} residues; ABA (A8S) with {len(L)} heavy atoms\n")


def min_dist(r):
    return min(float(np.linalg.norm(L - a.coord, axis=1).min())
               for a in r if a.element != "H")


def cb_cos(r):
    """cos(angle) between CA->CB and CA->nearest-ABA-atom. None for Gly."""
    if "CA" not in r or "CB" not in r:
        return None
    ca, cb = r["CA"].coord, r["CB"].coord
    j = int(np.linalg.norm(L - ca, axis=1).argmin())
    v1, v2 = cb - ca, L[j] - ca
    return float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))


sel_5, sel_6 = [], []
for rn, r in sorted(res.items()):
    d = min_dist(r)
    if d <= 5.0:
        sel_5.append(rn)
    elif d <= 6.0:
        c = cb_cos(r)
        if c is not None and c > 0:
            sel_6.append(rn)

selected = sorted(sel_5 + sel_6)
print("=== reproducing the Beltran geometric filter on 3K90 chain A ===")
print(f"  rule 1 (any atom <= 5.0 A of ABA):            {len(sel_5)} residues")
print(f"    {sel_5}")
print(f"  rule 2 (<= 6.0 A AND CA-CB not away from ABA): {len(sel_6)} residues")
print(f"    {sel_6}")
print(f"\n  TOTAL selected by the automated filter: {len(selected)}")
print(f"    {selected}")

print(f"\n  Beltran report 25 'close contact' residues; "
      f"this reproduction gives {len(selected)}.")
print(f"  They randomised 19 after removing 6: {DSM19}")
print(f"  filter minus the 6 removed = "
      f"{sorted(set(selected) - set(REMOVED))}")
missing = sorted(set(DSM19) - set(selected))
extra = sorted(set(selected) - set(REMOVED) - set(DSM19))
print(f"  in their 19 but not reproduced here: {missing}")
print(f"  reproduced here but not in their 19: {extra}")

print("\n" + "=" * 74)
print("R79: was it excluded BY THE RULE, or by manual curation?")
print("=" * 74)
r79 = res[79]
d79, c79 = min_dist(r79), cb_cos(r79)
print(f"  min heavy-atom distance to ABA : {d79:.2f} A")
print(f"  <= 5.0 A (rule 1)              : {d79 <= 5.0}")
print(f"  <= 6.0 A (rule 2 distance)     : {d79 <= 6.0}")
print(f"  CA->CB . CA->ABA cosine        : {c79:+.3f} "
      f"({'TOWARD ABA' if c79 > 0 else 'away from ABA'})")
print(f"  rule 2 satisfied               : {d79 <= 6.0 and c79 > 0}")
print(f"\n  --> R79 {'PASSES' if 79 in selected else 'FAILS'} the automated "
      f"geometric filter.")
if 79 in selected:
    print("      It was therefore removed at the MANUAL CURATION step, on the")
    print("      stated grounds: 'Appears to have only second-shell effects on")
    print("      ligand binding.' That is a judgement, not a measurement --")
    print("      and it is the assumption this project tests directly.")

print("\n" + "=" * 74)
print("the six manually removed residues, independently measured")
print("=" * 74)
print(f"{'res':<8}{'minABA':>8}{'CA-CB cos':>11}{'orientation':>14}   stated reason")
for rn in sorted(REMOVED):
    r = res[rn]
    d = min_dist(r)
    c = cb_cos(r)
    orient = "n/a (Gly/Pro)" if c is None else (
        "TOWARD ABA" if c > 0 else "away from ABA")
    cs = "  n/a" if c is None else f"{c:+.3f}"
    print(f"{r.get_resname()}{rn:<5}{d:>8.2f}{cs:>11}{orient:>14}   "
          f"{REMOVED[rn][:44]}")

print("\nNotes on agreement with this project's independent measurements")
print("-" * 74)
print("  T91  : rigid truncation gave dV = +0.0 A^3 -- fully consistent with")
print("         'outward facing'.")
print("  F61  : rigid truncation gave dV = -16.5 A^3, i.e. a SEALING residue")
print("         whose removal breaches the cavity to solvent. Consistent with")
print("         'outward facing', reached by a different route.")
print("  P88  : gate loop; excluded by this project as well.")
print("  H115 : latch loop; excluded by this project as well.")
print("  H60  : dimerisation; not a cavity residue here either.")
print("  R79  : the ONLY one of the six whose stated reason this project's")
print("         data contradicts -- see scripts/04_cavity_lining.py:")
print("         R79 NH2 and NE have CLEAR line of sight to ABA, R79 borders")
print("         the WT cavity, and R79A/E94A unlocks +84.8 A^3.")
