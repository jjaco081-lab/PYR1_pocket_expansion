#!/usr/bin/env python
"""
15_carotenoid_start.py -- how big is a START-fold pocket that actually holds a
carotenoid, measured on the same scale as PYR1?

Background
----------
ABA is itself an APOcarotenoid -- the product of 9-cis-epoxycarotenoid
dioxygenase (NCED) cleavage of violaxanthin/neoxanthin. So PYR1 already binds a
carotenoid-derived fragment. The question is whether the same fold ever binds an
INTACT C40 carotenoid, which would make an ancestral carotenoid-binding PYR1
plausible rather than speculative.

It does. Bombyx mori carotenoid-binding protein (BmCBP) is a bona fide
START-domain protein (STARD3-like fold, orthologous to MLN64) that binds lutein
in vivo and confers the yellow-cocoon phenotype. Crystal structures:

  7ZVR  BmCBP + ZEAXANTHIN (ligand code ZEX)   <-- the holo reference
  7ZTQ  BmCBP apo
  7ZTR / 7ZTU / 7ZVQ  apo point mutants (W232F, D162L, S206V)
  8AAQ  "CRT-416" form

Reference: Structural basis for the carotenoid binding and transport function
of a START domain, Structure (2022), PMID 36356587.

Note the reported architecture: an Omega-1 loop pinned onto the alpha-4 helix by
an R173-D279 SALT BRIDGE -- i.e. this fold also uses an arginine salt bridge to
position a loop over the cavity, structurally analogous to PYR1's R79-E94
(section 2 of the README).

This script measures, on the SAME grid criterion used for PYR1 throughout this
project:
  * the enclosed cavity of holo BmCBP (7ZVR), seeded at the zeaxanthin centroid
  * the enclosed cavity of apo BmCBP (7ZTQ), unbiased largest component
  * ligand heavy-atom counts and end-to-end extents for ABA vs zeaxanthin
so the PYR1 pocket and a working carotenoid pocket are on one scale.

Run with the esmfold2 env python.
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_cavity import cavity_volume, VDW, DIRS  # noqa: E402
from scipy import ndimage                        # noqa: E402
from scipy.spatial import cKDTree                # noqa: E402

DATA = os.path.join(ROOT, "data")


def load(cif):
    return MMCIFParser(QUIET=True).get_structure("x", os.path.join(DATA, cif))[0]


def prot_atoms(model, chain=None):
    xyz, el = [], []
    for ch in model:
        if chain and ch.id != chain:
            continue
        for r in ch:
            if not is_aa(r):
                continue
            for a in r:
                if a.element == "H":
                    continue
                xyz.append(a.coord); el.append(a.element)
    return np.array(xyz), el


def lig_atoms(model, code):
    out = [a for ch in model for r in ch if r.get_resname() == code
           for a in r if a.element != "H"]
    return np.array([a.coord for a in out])


def extent(P):
    """max pairwise distance (end-to-end length) and radius of gyration."""
    d = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
    rg = float(np.sqrt(((P - P.mean(0)) ** 2).sum(1).mean()))
    return float(d.max()), rg


def largest_cavity(xyz, elem, spacing=0.5, probe=1.4, bur_cut=0.88):
    xyz = np.asarray(xyz, float)
    rad = np.array([VDW.get(e, 1.70) for e in elem])
    lo, hi = xyz.min(0) - 3, xyz.max(0) + 3
    axes = [np.arange(lo[i], hi[i], spacing) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    shape = grid.shape[:3]
    pts = grid.reshape(-1, 3)
    tree = cKDTree(xyz)
    d, i = tree.query(pts, k=1)
    clear = d - rad[i]
    occ = (clear < 0).reshape(shape)
    cand = np.where(clear > probe)[0]
    if cand.size == 0:
        return 0.0
    cp = pts[cand]
    steps = np.arange(1.0, 15.0, 0.75)
    hits = np.zeros(len(cp))
    shp = np.array(shape)
    for dv in DIRS:
        blocked = np.zeros(len(cp), bool)
        for s in steps:
            q = cp + dv * s
            idx = ((q - lo) / spacing).astype(int)
            ok = np.all((idx >= 0) & (idx < shp), axis=1)
            b = np.zeros(len(cp), bool)
            v = idx[ok]
            b[ok] = occ[v[:, 0], v[:, 1], v[:, 2]]
            blocked |= b
        hits += blocked
    mask = np.zeros(len(pts), bool)
    mask[cand[(hits / len(DIRS)) >= bur_cut]] = True
    lab, n = ndimage.label(mask.reshape(shape))
    if n == 0:
        return 0.0
    sizes = ndimage.sum(mask.reshape(shape), lab, range(1, n + 1))
    return float(sizes.max() * spacing ** 3)


print("=" * 74)
print("LIGAND GEOMETRY -- what has to fit")
print("=" * 74)
pyr1 = load("3QN1.cif")
aba = lig_atoms(pyr1, "A8S")
zx_m = load("7ZVR.cif")
zx = lig_atoms(zx_m, "ZEX")
for name, P in (("ABA (A8S)", aba), ("zeaxanthin (ZEX, modelled)", zx)):
    L, rg = extent(P)
    print(f"  {name:<26} {len(P):>3} heavy atoms   end-to-end {L:6.2f} A   Rg {rg:5.2f} A")

print("""
  !! ZEAXANTHIN IS ONLY HALF-MODELLED IN 7ZVR.
     Zeaxanthin is C40H56O2 = 42 heavy atoms, but only 20 are deposited
     (atoms C1-C20 plus O3, i.e. one beta-ionone ring and the polyene as far
     as C14/C20), at occupancies of 0.60-0.65. The distal half (C21-C40 and the
     second ring) is disordered and unmodelled.

     Consequences, which must not be glossed over:
       * the 14.06 A "end-to-end" above is for the ordered HALF; intact
         zeaxanthin is roughly 30 A long.
       * any cavity volume seeded on the modelled fragment UNDER-reports the
         true binding channel, because the disordered half's channel is not
         defined by the coordinates and, being disordered, is likely solvent
         exposed -- which this project's buriedness criterion (0.88) excludes
         by construction.
     Treat the BmCBP numbers below as a LOWER BOUND and as qualitative
     evidence about pocket SHAPE, not as a calibrated volume.""")

print("\n" + "=" * 74)
print("CAVITY VOLUMES on this project's grid criterion")
print("=" * 74)

# PYR1 holo, seeded at ABA centroid (matches script 02 / README section 3)
px, pe = prot_atoms(pyr1, "A")
v_pyr1, _ = cavity_volume(px, pe, aba.mean(0))
print(f"  PYR1 (3QN1 chain A), seeded at ABA        {v_pyr1:8.1f} A^3")

# BmCBP holo, seeded at zeaxanthin centroid
bx, be = prot_atoms(zx_m)
v_holo, _ = cavity_volume(bx, be, zx.mean(0), box=22.0)
print(f"  BmCBP holo (7ZVR), seeded at zeaxanthin   {v_holo:8.1f} A^3")

# BmCBP apo, unbiased largest cavity
apo = load("7ZTQ.cif")
ax, ae = prot_atoms(apo)
v_apo = largest_cavity(ax, ae)
print(f"  BmCBP apo  (7ZTQ), largest cavity         {v_apo:8.1f} A^3")

print(f"\n  ratio BmCBP holo : PYR1 = {v_holo / max(v_pyr1, 1):.1f}x")
print(f"  PYR1 expanded quad target (README section 3) was 375.5 A^3 "
      f"= {375.5 / max(v_pyr1,1):.1f}x WT")

print("\n" + "=" * 74)
print("INTERPRETATION")
print("=" * 74)
print("""  ABA is an APOcarotenoid (NCED cleavage product of violaxanthin/neoxanthin),
  so PYR1 already binds a carotenoid FRAGMENT. BmCBP shows the same START fold
  binding an INTACT C40 xanthophyll. An ancestral carotenoid-binding member of
  this fold is therefore not speculative -- the fold demonstrably supports it.

  Caveat for ligand selection: the carotenoid pocket is not merely bigger, it is
  a long hydrophobic TUNNEL matched to an extended polyene. PYR1's cavity is
  globular. Expanding PYR1 toward full C40 carotenoids would require elongation,
  not just volume, and would almost certainly break the gate/latch closure that
  the HAB1 readout depends on. The realistic targets are mid-size apocarotenoids
  and comparably sized rigid hydrophobes -- see the ligand panel note in the
  README.""")
