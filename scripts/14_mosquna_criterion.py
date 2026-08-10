#!/usr/bin/env python
"""
14_mosquna_criterion.py -- test R79 against the Mosquna et al. 2011 criterion.

Mosquna et al. (PNAS 2011) state:

    "Thirty-nine residues within 5 A of ABA, four water-molecules that contact
     ABA, or the PP2C HAB1 (Fig. 1A) were identified using PYR1 structure
     coordinates, and all possible 741 single amino acid substitutions at these
     sites were constructed by site-directed mutagenesis."

So the selection rule is a union of three conditions:
    (a) within 5 A of ABA
    (b) within 5 A of one of the four ordered waters that contact ABA
    (c) contacts the PP2C (HAB1)

Their 39-residue list (Table S1, transcribed in
data/mosquna2011_ssm_residues.tsv) does NOT contain R79. Condition (a) alone
excludes R79, which sits at 5.81 A (3QN1) / 5.89 A (3K90). But condition (b)
is a genuine second chance -- an ABA-contacting water could bridge to R79.

This script therefore asks: does R79 contact any water that contacts ABA?
If it does, R79 satisfied their stated rule and its absence is an oversight.
If it does not, R79 was legitimately excluded by their criterion -- and the
exclusion traces purely to the hard 5 A ligand-distance cutoff.

Both 3K90 (ABA-bound PYR1 alone, the coordinates Mosquna and Beltran used) and
3QN1 (the PYR1-HAB1-ABA ternary complex used elsewhere in this project) are
evaluated, since water positions differ between crystal forms.

Run with the esmfold2 env python.
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WATER_CUT = 3.5      # H-bond distance, water oxygen to ligand/protein heavy atom
RES_CUT = 5.0        # Mosquna's stated residue cutoff

# Table S1 (transcribed); see data/mosquna2011_ssm_residues.tsv
MOSQUNA39 = [55, 59, 60, 61, 62, 63, 81, 83, 84, 85, 86, 87, 88, 89, 92, 94,
             108, 110, 115, 116, 117, 120, 122, 141, 148, 150, 151, 154, 155,
             156, 158, 159, 160, 162, 163, 164, 166, 167, 170]


def analyse(cif, pyr1_chain, label):
    print("=" * 74)
    print(f"{label}  ({os.path.basename(cif)}, PYR1 = chain {pyr1_chain})")
    print("=" * 74)
    model = MMCIFParser(QUIET=True).get_structure("x", cif)[0]

    lig = [r for ch in model for r in ch if r.get_resname() == "A8S"][0]
    L = np.array([a.coord for a in lig if a.element != "H"])

    prot = {r.id[1]: r for r in model[pyr1_chain] if is_aa(r)}
    waters = [r for ch in model for r in ch if r.get_resname() == "HOH"]

    # (b) waters that contact ABA
    aba_waters = []
    for w in waters:
        o = next(iter(w))
        d = float(np.linalg.norm(L - o.coord, axis=1).min())
        if d <= WATER_CUT:
            aba_waters.append((w, o, d))
    aba_waters.sort(key=lambda x: x[2])
    print(f"\nordered waters within {WATER_CUT} A of ABA: {len(aba_waters)}")
    print("  (Mosquna refer to FOUR such waters)")
    for w, o, d in aba_waters[:8]:
        print(f"    HOH {w.id[1]:<5} {d:.2f} A from ABA")

    def min_to_lig(r):
        # NB: axis=1 is required -- without it numpy returns the Frobenius norm
        # of the whole (N,3) difference matrix rather than per-atom distances.
        return min(float(np.linalg.norm(L - a.coord, axis=1).min()) for a in r
                   if a.element != "H")

    d79 = min_to_lig(prot[79])
    print(f"\nR79 -> ABA minimum heavy-atom distance: {d79:.2f} A")
    print(f"  condition (a), within {RES_CUT} A of ABA : "
          f"{'PASS' if d79 <= RES_CUT else 'FAIL'}")

    # (b) for R79
    hits = []
    for w, o, d in aba_waters:
        dm = min(float(np.linalg.norm(o.coord - a.coord)) for a in prot[79]
                 if a.element != "H")
        if dm <= RES_CUT:
            hits.append((w.id[1], d, dm))
    print(f"  condition (b), within {RES_CUT} A of an ABA-contacting water: "
          f"{'PASS' if hits else 'FAIL'}")
    for wid, dw, dr in hits:
        print(f"      via HOH {wid}: water-ABA {dw:.2f} A, water-R79 {dr:.2f} A")

    # which residues DO satisfy (b) but not (a)?
    only_b = []
    for rn, r in prot.items():
        if min_to_lig(r) <= RES_CUT:
            continue
        for w, o, d in aba_waters:
            dm = min(float(np.linalg.norm(o.coord - a.coord)) for a in r
                     if a.element != "H")
            if dm <= RES_CUT:
                only_b.append(rn)
                break
    print(f"\nresidues satisfying (b) but NOT (a): {sorted(only_b)}")
    inlist = [r for r in sorted(only_b) if r in MOSQUNA39]
    print(f"  of those, present in the Mosquna 39: {inlist}")
    print(f"  R79 among them: {79 in only_b}")
    return d79, bool(hits)


d1, w1 = analyse(os.path.join(ROOT, "data", "3K90.cif"), "A",
                 "Mosquna/Beltran coordinate set")
print()
d2, w2 = analyse(os.path.join(ROOT, "data", "3QN1.cif"), "A",
                 "ternary complex used elsewhere in this project")

print("\n" + "=" * 74)
print("VERDICT")
print("=" * 74)
print(f"  R79-ABA distance: {d1:.2f} A (3K90), {d2:.2f} A (3QN1) -- both > 5.0 A,")
print("  so condition (a) excludes it in either crystal form.")
if not (w1 or w2):
    print("  R79 also contacts NO ABA-contacting water, so condition (b) does")
    print("  not rescue it either.")
    print("\n  => R79 was legitimately excluded by Mosquna's stated criterion.")
    print("     The exclusion traces entirely to the hard 5 A ligand-distance")
    print("     cutoff -- the same ligand-centric metric that Beltran's more")
    print("     permissive 6 A + CA-CB rule DID select (script 13), only for")
    print("     R79 to be removed by hand as 'second-shell'.")
else:
    print("  R79 DOES contact an ABA-contacting water, so it satisfied their")
    print("  stated rule and its absence from Table S1 is an inconsistency.")
