#!/usr/bin/env python
"""
17_ligand_ladder.py -- match candidate ligands to the measured cavity spans.

Cavity spans from 16_cavity_shape.py (rigid-backbone, closed state):
    WT                    10.6 A
    K59A/F108A            16.0 A
    K59A/F108A/R79A/E94A  17.9 A   (quad; long axis 16.7 A, volume 391 A3)
    quad + Y120A/E141A    17.9 A   (volume 479 A3)

Ligands are generated from SMILES, MMFF-minimised, and measured the same way as
the cavity: max pairwise heavy-atom distance = length; max perpendicular
distance from the principal axis = half-width. Rigid rod-like ligands must fit
the span; flexible ones can curl and are less constrained.
"""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

LIGANDS = [
    ("ABA (reference)",        "CC1=CC(=O)CC(C)(C)C1(O)/C=C/C(C)=C/C(=O)O"),
    ("beta-ionone",            "CC1=C(C(C)(C)CCC1)/C=C/C(C)=O"),
    ("retinal",                "CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/C=O"),
    ("retinoic acid",          "CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/C(=O)O"),
    ("GR24 (strigolactone)",   "CC1=C(C(=O)O[C@@H]1O[C@H]2C[C@@H]3CCC(=O)C3=C2)C"),
    ("mandipropamid",          "CC(C)(COC)C(=O)N[C@@H](CC1=CC=C(OCC#C)C=C1)C(=O)NC2=CC=C(Cl)C=C2"),
    ("bixin",                  "COC(=O)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C(=O)O"),
    ("crocetin",               "OC(=O)/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C(=O)O"),
    ("beta-apo-8'-carotenal",  "CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=O"),
    ("zeaxanthin (C40)",       "CC1=C(C(CC(C1)O)(C)C)/C=C/C(=C/C=C(\\C)/C=C/C=C(\\C)/C=C/C=C(\\C)/C=C/C2=C(CC(CC2(C)C)O)C)/C"),
    ("beta-carotene (C40)",    "CC1=C(C(C)(C)CCC1)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C2=C(C)CCCC2(C)C"),
]

def measure(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    mh = Chem.AddHs(m)
    if AllChem.EmbedMolecule(mh, randomSeed=0xC0FFEE) != 0: return None
    AllChem.MMFFOptimizeMolecule(mh, maxIters=2000)
    hv = Chem.RemoveHs(mh)
    P = hv.GetConformer().GetPositions()
    d = np.linalg.norm(P[:,None]-P[None], axis=-1)
    Pc = P - P.mean(0)
    u = np.linalg.eigh(np.cov(Pc.T))[1][:,-1]
    t = Pc @ u
    perp = np.linalg.norm(Pc - np.outer(t,u), axis=1)
    return (float(d.max()), float(perp.max()), P.shape[0],
            Descriptors.MolWt(m), Descriptors.NumRotatableBonds(m))

SPANS = [("WT",10.6),("K59A/F108A",16.0),("quad",17.9)]
print(f"{'ligand':<26}{'MW':>7}{'nHV':>5}{'len A':>8}{'halfW':>7}{'rotB':>6}   fits (rigid, closed-state span)")
print("-"*104)
for name, smi in LIGANDS:
    r = measure(smi)
    if r is None:
        print(f"{name:<26}  (embedding failed)"); continue
    L, w, n, mw, rb = r
    verdict = ", ".join(f"{lbl}:{'Y' if L<=sp else 'N'}" for lbl,sp in SPANS)
    print(f"{name:<26}{mw:>7.1f}{n:>5}{L:>8.1f}{w:>7.2f}{rb:>6}   {verdict}")
print("""
  'fits' compares ligand LENGTH against the cavity's longest internal straight
  segment in the CLOSED state with a rigid backbone. It is a necessary, not
  sufficient, condition: width and the narrowest cross-section still apply, and
  flexible ligands (high rotB) can curl to fit a shorter span.

  Note BmCBP orders only HALF of its bound zeaxanthin (7ZVR), so a C40
  carotenoid does NOT require a fully enclosed 28 A tunnel -- a buried ionone
  ring plus polyene with a solvent-exposed distal tail is the physiological
  binding mode in the one START-fold carotenoid complex that has been solved.""")
