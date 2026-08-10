#!/usr/bin/env python
"""
01_pocket_contacts.py -- define the pocket and map the buttress network.

Three outputs:
  (a) every residue with a heavy atom within 5.0 A of ABA (the first shell),
      reported per chain so the PYR1 vs HAB1 split is explicit;
  (b) the K59 / E94 / R79 polar network;
  (c) the residues packed against the F108 ring centroid.

The point of (b) and (c) is that a 5 A contact-radius scan CANNOT see R79 --
it never approaches the ligand -- yet R79 pins the F108 rotamer and is
salt-bridged by E94. The wall is a buttressed cluster, not a single residue.

Run with the esmfold2 env python.
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
model = MMCIFParser(QUIET=True).get_structure("3qn1",
            os.path.join(ROOT, "data", "3QN1.cif"))[0]

lig = [r for ch in model for r in ch if r.get_resname() == "A8S"][0]
prot = [a for ch in model for r in ch if is_aa(r) for a in r]
ns = NeighborSearch(prot)

print("=== (a) residues within 5.0 A of ABA (A8S) ===")
near = {}
for a in lig:
    for b in ns.search(a.coord, 5.0):
        r = b.get_parent()
        key = (r.get_parent().id, r.id[1], r.get_resname())
        near[key] = min(near.get(key, 99), float(np.linalg.norm(a.coord - b.coord)))
for (ch, rn, nm), d in sorted(near.items(), key=lambda x: x[1]):
    tag = "PYR1" if ch == "A" else "HAB1"
    print(f"  {tag} {nm}{rn:<5} {d:.2f}")
n_hab1 = sum(1 for k in near if k[0] == "B")
print(f"\n  --> {len(near)-n_hab1} PYR1 residues, {n_hab1} HAB1 residue(s).")
print("      HAB1 contributes essentially nothing to the pocket: it reads the")
print("      closed-state SURFACE. This is the ratchet, not a shared site.")

A = {r.id[1]: r for r in model["A"] if is_aa(r)}

print("\n=== (b) K59 NZ polar contacts (<4.5 A) ===")
nz = A[59]["NZ"].coord
for rn, r in A.items():
    if rn == 59:
        continue
    for a in r:
        if a.element in ("O", "N"):
            d = np.linalg.norm(nz - a.coord)
            if d < 4.5:
                print(f"  K59 NZ -- {r.get_resname()}{rn} {a.get_id():<4} {d:.2f}")
for a in lig:
    d = np.linalg.norm(nz - a.coord)
    if d < 4.5:
        print(f"  K59 NZ -- ABA {a.get_id():<4} {d:.2f}")

print("\n=== (b') E94 / E141 carboxylate -> N contacts (<4.5 A) ===")
for e in (94, 141):
    for on in ("OE1", "OE2"):
        if on not in A[e]:
            continue
        c = A[e][on].coord
        for rn, r in A.items():
            if rn == e:
                continue
            for a in r:
                if a.element == "N":
                    d = np.linalg.norm(c - a.coord)
                    if d < 4.5:
                        print(f"  E{e} {on} -- {r.get_resname()}{rn} {a.get_id():<4} {d:.2f}")

print("\n=== (c) residues near the F108 ring centroid (<7.5 A) ===")
ring = np.array([A[108][n].coord for n in
                 ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")]).mean(0)
out = []
for rn, r in A.items():
    if rn == 108:
        continue
    d = min(float(np.linalg.norm(ring - a.coord)) for a in r)
    if d < 7.5:
        out.append((d, r.get_resname(), rn))
for d, nm, rn in sorted(out)[:12]:
    print(f"  {nm}{rn:<5} {d:.2f}")
print(f"  ABA        {min(float(np.linalg.norm(ring-a.coord)) for a in lig):.2f}")
print("\n  --> R79 guanidinium stacks on the F108 ring (cation-pi) and is")
print("      locked by the E94 salt bridge. F108 is buttressed.")
