#!/usr/bin/env python
r"""
151_cavity_bottleneck.py -- decompose each cavity into chambers by the size of
the sphere that can pass between them. No ligand is used anywhere.

TWO CORRECTIONS THIS SUPERSEDES
-------------------------------
149 scored pocket membership by distance to the crystallographic ABA. That is
invalid: ABA is 19 heavy atoms in a cavity we are trying to enlarge, so scoring
membership by proximity to it defines the pocket as the volume already occupied
and excludes every expansion by construction. Its verdict is withdrawn, and so
is 148's 9 A adjacency filter, which cut 43 of 52 candidate positions on the
same basis.

150 asked the right question -- one chamber or two? -- but measured the answer
with the wrong quantity. It reported r_eq, the radius of the equal-area circle
of a cross-section of GRID POINTS. Those points are probe CENTRES: every one of
them already satisfies clear > 1.4 A, so a water fits at each by construction
and a narrow ribbon of them does not mean a narrow channel. r_eq is a property
of how many probe positions exist there, not of how wide the passage is.

THE RIGHT QUANTITY. For each grid point, `clear` = distance to the nearest van
der Waals surface, so the largest sphere centred there has radius `clear`. The
bottleneck between two chambers is then the MAXIMIN of `clear` over all paths
joining them -- the biggest sphere that can be walked from one to the other.
Computed by thresholding at increasing r and testing connectivity, which is
exact on the grid and, unlike 150's slabs, does not assume the pocket is
straight.

READING THE OUTPUT
  r_bottle >= 1.4   a water passes -- the chambers are hydraulically one
  r_bottle >= 2.4   roughly a methyl/methylene can pass
  r_bottle >= 3.0   a ligand ring can pass; this is one pocket for design
Chambers are reported with their own volumes so an inflated total is visible.
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
read_pdb, read_cif_resnames, THREE = m143.read_pdb, m143.read_cif_resnames, m143.THREE

HC = os.path.join(ROOT, "results", "homolog_cavities")
OUT = os.path.join(ROOT, "results", "pocket_shape")
SPACING = 0.6
PROBE = 1.4
TARGETS = ["2pcsA00", "2ns9B01", "2bk0A00", "6awvC00", "3qrzA00", "3oquB00",
           "2flhA00", "4dsbB00", "3tfzE00", "1z94E00"]


def clear_map(xyz, elem, spacing=SPACING):
    xyz = np.asarray(xyz, float)
    rad = np.array([VDW.get(e, 1.70) for e in elem])
    lo = xyz.min(0) - 3
    axes = [np.arange(lo[i], xyz.max(0)[i] + 3, spacing) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    shape = grid.shape[:3]
    pts = grid.reshape(-1, 3)
    d, i = cKDTree(xyz).query(pts, k=1)
    return (d - rad[i]).reshape(shape), pts.reshape(shape + (3,)), lo, shape


def enclosed(clear, lo, shape, bur_cut=0.88, ray_max=15.0, ray_step=0.75):
    occ = clear < 0
    cand = np.argwhere(clear > PROBE)
    cp = lo + cand * SPACING
    hits = np.zeros(len(cp))
    shp = np.array(shape)
    for dv in DIRS:
        blocked = np.zeros(len(cp), bool)
        for s in np.arange(1.0, ray_max, ray_step):
            idx = ((cp + dv * s - lo) / SPACING).astype(int)
            ok = np.all((idx >= 0) & (idx < shp), axis=1)
            v = idx[ok]
            b = np.zeros(len(cp), bool)
            b[ok] = occ[v[:, 0], v[:, 1], v[:, 2]]
            blocked |= b
        hits += blocked
    keep = cand[(hits / len(DIRS)) >= bur_cut]
    m = np.zeros(shape, bool)
    m[keep[:, 0], keep[:, 1], keep[:, 2]] = True
    lab, n = ndimage.label(m)
    if n == 0:
        return np.zeros(shape, bool)
    sizes = ndimage.sum(m, lab, range(1, n + 1))
    return lab == int(np.argmax(sizes)) + 1


def bottleneck(clear, mask, a, b, hi):
    """Largest r such that a and b stay connected through {clear > r} & mask."""
    lo_r = PROBE
    for _ in range(24):
        mid = 0.5 * (lo_r + hi)
        sub = mask & (clear > mid)
        lab, n = ndimage.label(sub)
        if n and lab[tuple(a)] and lab[tuple(a)] == lab[tuple(b)]:
            lo_r = mid
        else:
            hi = mid
    return lo_r


def analyse(name, xyz, elem):
    clear, _, lo, shape = clear_map(xyz, elem)
    mask = enclosed(clear, lo, shape)
    tot = float(mask.sum() * SPACING ** 3)
    if tot == 0:
        return None
    # chambers = components that survive a 2.4 A sphere; each keeps its own seed
    core = mask & (clear > 2.4)
    lab, n = ndimage.label(core)
    if n == 0:
        return dict(name=name, total=round(tot, 1), chambers=[], note="no 2.4 A chamber")
    sizes = ndimage.sum(core, lab, range(1, n + 1))
    order = np.argsort(sizes)[::-1]
    seeds = []
    for k in order[:4]:
        idx = np.argwhere(lab == k + 1)
        seeds.append(idx[np.argmax(clear[tuple(idx.T)])])
    # assign every cavity voxel to its nearest chamber seed (grid distance)
    vol = []
    if len(seeds) == 1:
        vol = [tot]
    else:
        cvox = np.argwhere(mask)
        S = np.array(seeds)
        d = np.linalg.norm(cvox[:, None] - S[None], axis=-1)
        own = d.argmin(1)
        vol = [float((own == j).sum() * SPACING ** 3) for j in range(len(seeds))]
    rows = []
    for j in range(1, len(seeds)):
        r = bottleneck(clear, mask, seeds[0], seeds[j],
                       float(min(clear[tuple(seeds[0])], clear[tuple(seeds[j])])))
        rows.append(dict(vol=round(vol[j], 1), r_bottle=round(r, 2)))
    return dict(name=name, total=round(tot, 1),
                main=round(vol[0], 1), satellites=rows,
                r_max=round(float(clear[mask].max()), 2))


def main():
    out = []
    at = read_pdb(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb"), want_atom=False)
    prot = [a for a in at if a[1] in THREE]
    out.append(analyse("PYR1", [a[4] for a in prot], [a[3] for a in prot]))
    for t in TARGETS:
        dom = read_pdb(os.path.join(HC, "domains", t + ".pdb"))
        r = analyse(t.upper()[:4] + "_" + t[4], [a[4] for a in dom], [a[3] for a in dom])
        if r:
            out.append(r)
        print(f"  {t} done", flush=True)

    print("\n" + "=" * 80)
    print("CAVITY DECOMPOSED INTO CHAMBERS -- no ligand used")
    print("=" * 80)
    print(f"{'structure':<10}{'total':>8}{'main':>8}{'r_max':>8}   satellites (volume @ bottleneck radius)")
    for r in out:
        if r is None:
            continue
        sat = ("  ".join(f"{s['vol']:.0f} A^3 @ r={s['r_bottle']:.2f}"
                         for s in r.get("satellites", [])) or "none")
        print(f"{r['name']:<10}{r['total']:>8.0f}{r.get('main', 0):>8.0f}"
              f"{r.get('r_max', 0):>8.2f}   {sat}")
    print("\n   r_max = radius of the largest sphere that fits anywhere in the cavity.")
    print("   A satellite at r < 1.4 is not even water-connected; at r < 2.4 no")
    print("   ligand atom passes and its volume should NOT be counted as pocket.")
    json.dump(out, open(os.path.join(OUT, "cavity_chambers.json"), "w"), indent=1)
    print(f"\nwritten to {OUT}/cavity_chambers.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
