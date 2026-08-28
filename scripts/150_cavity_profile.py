#!/usr/bin/env python
r"""
150_cavity_profile.py -- is the pocket one chamber or two voids glued together?

RETRACTION FIRST. Script 149 audited lining calls by distance to the
crystallographic ABA and concluded 3OQU has "fewer genuine contacts". That test
was wrong in principle and its verdict is withdrawn. ABA is 19 heavy atoms in a
cavity we are explicitly trying to make bigger; scoring pocket membership by
proximity to it defines the pocket as the volume already occupied and rules out
every expansion by construction. The same mistake sits in 148, whose 9 A
adjacency cut removed 43 of 52 candidate positions on that basis; that filter is
withdrawn too and replaced in 151 by contiguous volume gain.

THE OBJECTION IS STILL LIVE, THOUGH, and this is the honest form of it. The
cavity is taken as the largest enclosed connected component, so if two separate
voids touch anywhere along a thin channel, the labelling reports them as one
pocket and the lining count adds up both. That is what Jannis is seeing in
PyMOL. It is testable without reference to any ligand:

  PROFILE   slice the cavity point cloud along its own long axis and measure the
            cross-sectional area of each slab. One chamber gives a single
            envelope; two chambers joined by a channel give a WAIST.
  WAIST     the minimum cross-section between the two ends, in A^2, and the
            radius of the equivalent circle. A waist under ~1.4 A radius cannot
            pass a water, let alone a ligand, and means the far lobe is not
            usable pocket however contiguous the grid says it is.
  BURIEDNESS
            the ray-occlusion fraction per grid point. Points that only just
            clear the 0.88 cut are surface groove, not interior. Reported for
            each slab so a lobe that is really a groove shows up as a
            buriedness dip rather than an area dip.

The ligand is used ONLY to annotate where along the axis it sits. It is not used
to decide what counts as pocket.
"""
import json
import os
import sys

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_cavity import VDW, DIRS                                    # noqa: E402
from importlib import import_module                                 # noqa: E402
m143 = import_module("143_pocket_backbone_shape")                   # noqa: E402
m149 = import_module("149_pocket_selections")                       # noqa: E402
read_pdb, read_cif_resnames, THREE = m143.read_pdb, m143.read_cif_resnames, m143.THREE

HC = os.path.join(ROOT, "results", "homolog_cavities")
OUT = os.path.join(ROOT, "results", "pocket_shape")
SPACING = 0.6
SLAB = 1.5          # A along the axis


def cavity_points(xyz, elem, spacing=SPACING, probe=1.4, bur_cut=0.88,
                  ray_max=15.0, ray_step=0.75):
    """Largest enclosed component, returning points AND their buriedness."""
    xyz = np.asarray(xyz, float)
    rad = np.array([VDW.get(e, 1.70) for e in elem])
    lo = xyz.min(0) - 3
    axes = [np.arange(lo[i], xyz.max(0)[i] + 3, spacing) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    shape = grid.shape[:3]
    pts = grid.reshape(-1, 3)
    d, i = cKDTree(xyz).query(pts, k=1)
    clear = d - rad[i]
    occ = (clear < 0).reshape(shape)
    cand = np.where(clear > probe)[0]
    cp = pts[cand]
    hits = np.zeros(len(cp))
    shp = np.array(shape)
    for dv in DIRS:
        blocked = np.zeros(len(cp), bool)
        for s in np.arange(1.0, ray_max, ray_step):
            idx = ((cp + dv * s - lo) / spacing).astype(int)
            ok = np.all((idx >= 0) & (idx < shp), axis=1)
            v = idx[ok]
            b = np.zeros(len(cp), bool)
            b[ok] = occ[v[:, 0], v[:, 1], v[:, 2]]
            blocked |= b
        hits += blocked
    bur = hits / len(DIRS)
    mask = np.zeros(len(pts), bool)
    burf = np.zeros(len(pts))
    sel = cand[bur >= bur_cut]
    mask[sel] = True
    burf[cand] = bur
    lab, n = ndimage.label(mask.reshape(shape))
    sizes = ndimage.sum(mask.reshape(shape), lab, range(1, n + 1))
    pick = int(np.argmax(sizes)) + 1
    keep = lab.reshape(-1) == pick
    return pts[keep], burf[keep], float(keep.sum() * spacing ** 3)


def profile(name, P, bur, lig=None):
    Q = P - P.mean(0)
    V = np.linalg.svd(Q, full_matrices=False)[2]
    t = Q @ V[0]
    lo, hi = t.min(), t.max()
    edges = np.arange(lo, hi + SLAB, SLAB)
    print(f"\n{name}: {len(P)} grid points, "
          f"{len(P)*SPACING**3:.0f} A^3, long axis {hi-lo:.1f} A")
    if lig is not None:
        lt = (lig - P.mean(0)) @ V[0]
        print(f"   ligand occupies axis {lt.min():+.1f} to {lt.max():+.1f} "
              f"(cavity runs {lo:+.1f} to {hi:+.1f})")
    print(f"   {'axis':>7}{'area':>8}{'r_eq':>7}{'buried':>8}")
    areas = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (t >= a) & (t < b)
        area = m.sum() * SPACING ** 3 / SLAB
        req = np.sqrt(area / np.pi) if area > 0 else 0.0
        areas.append((0.5 * (a + b), area, req, bur[m].mean() if m.any() else 0))
        bar = "#" * int(min(area, 60) / 2)
        print(f"   {0.5*(a+b):>+7.1f}{area:>8.1f}{req:>7.2f}"
              f"{(bur[m].mean() if m.any() else 0):>8.2f}  {bar}")
    return areas, (lt.min(), lt.max()) if lig is not None else None


def waist(areas, ligspan):
    """Minimum cross-section between the ligand slabs and the far end."""
    if ligspan is None:
        return None
    inside = [i for i, (c, *_ ) in enumerate(areas)
              if ligspan[0] <= c <= ligspan[1]]
    if not inside:
        return None
    out = []
    for side, rng in (("distal", range(max(inside) + 1, len(areas))),
                      ("proximal", range(0, min(inside)))):
        seg = [areas[i] for i in rng]
        if len(seg) < 2:
            continue
        w = min(seg, key=lambda x: x[1])
        far = sum(x[1] * SLAB for x in seg)
        out.append((side, w[1], w[2], w[0], far))
    return out


def main():
    res = {}
    at = read_pdb(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb"), want_atom=False)
    prot = [a for a in at if a[1] in THREE]
    lig1 = np.array([a[4] for a in at if a[1] == "A8S"])
    P, B, v = cavity_points([a[4] for a in prot], [a[3] for a in prot])
    print("=" * 74)
    print("CAVITY CROSS-SECTION ALONG THE LONG AXIS   (the ligand only annotates)")
    print("=" * 74)
    a1, s1 = profile("PYR1 (3QN1 chain A)", P, B, lig1)

    tgt = "3oquB00"
    cif = os.path.join(HC, "raw", tgt + ".cif")
    dom = read_pdb(os.path.join(HC, "domains", tgt + ".pdb"))
    het = m149.cif_atoms(cif, "B", comp="A8S")
    lig2 = np.array([h[2] for h in het])
    P2, B2, v2 = cavity_points([a[4] for a in dom], [a[3] for a in dom])
    a2, s2 = profile("3OQU (chain B)", P2, B2, lig2)

    print("\n" + "=" * 74)
    print("WAIST TEST -- is the far lobe a chamber of the same pocket?")
    print("=" * 74)
    for nm, a, s in (("PYR1", a1, s1), ("3OQU", a2, s2)):
        w = waist(a, s)
        if not w:
            print(f"   {nm}: cavity does not extend past the ligand")
            continue
        for side, area, req, ax, far in w:
            call = ("CONTINUOUS" if req >= 1.4 else
                    "PINCHED -- a separate void, not usable pocket")
            print(f"   {nm:<6}{side:<9} lobe {far:>6.0f} A^3 beyond the ligand; "
                  f"narrowest cross-section {area:>5.1f} A^2 "
                  f"(r_eq {req:.2f} A) at axis {ax:+.1f}   {call}")
    print("\n   r_eq is the radius of the equal-area circle. A probe of 1.4 A")
    print("   (one water) needs r_eq >= 1.4 to pass; anything less means the")
    print("   grid labelled two voids as one component.")
    json.dump({"pyr1": a1, "3oqu": a2}, open(os.path.join(OUT, "cavity_profile.json"), "w"),
              indent=1)
    print(f"\nwritten to {OUT}/cavity_profile.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
