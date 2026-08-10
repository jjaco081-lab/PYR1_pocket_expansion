#!/usr/bin/env python
"""
04_cavity_lining.py -- which residues LINE the cavity (vs merely sit near ABA).

Motivation
----------
Ranking residues by distance-to-ligand answers the wrong question. ABA does not
fill the PYR1 cavity, so a residue can face the cavity surface directly while
still being >5 A from the ligand -- it simply lines a part of the pocket the
ligand does not occupy. R79 is exactly this case, as seen in the ChimeraX
cavity-finder surface.

This script measures cavity lining directly:

  (a) LINING SCORE -- for each residue, the number of cavity grid points within
      `cut` A of any of its heavy atoms, and the minimum distance from the
      residue to the cavity surface. Computed for both the WT cavity and the
      expanded (K59A/F108A/R79A/E94A) cavity, so residues that line only the
      expanded pocket are visible too.

  (b) OCCLUSION / LINE-OF-SIGHT -- for a chosen atom (default R79 NH1 and NH2),
      march along the straight segment toward the nearest ABA atom and report
      whether any OTHER PYR1 heavy atom's van der Waals sphere intersects that
      segment. This tests the specific claim "there is nothing between the
      arginine nitrogen and the ligand", independent of the 5 A cut-off.

Run with the esmfold2 env python.
"""
import os, sys, csv
import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_cavity import (VDW, DIRS, atoms_from_pose_pdb)  # noqa: E402
from scipy import ndimage  # noqa: E402

PDB = os.path.join(ROOT, "data", "pyr1_A.pdb")
LIG = os.path.join(ROOT, "data", "aba_xtal.pdb")
SPACING, PROBE, BUR_CUT, BOX = 0.5, 1.4, 0.88, 14.0
LINING_CUT = 4.5


def cavity_points(truncate=()):
    """Return the (N,3) coordinates of the cavity grid points."""
    xyz, el = atoms_from_pose_pdb(PDB, chain="A", truncate_to_ala=set(truncate))
    rad = np.array([VDW.get(e, 1.70) for e in el])
    lo, hi = REF - BOX, REF + BOX
    axes = [np.arange(lo[i], hi[i], SPACING) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    shape = grid.shape[:3]
    pts = grid.reshape(-1, 3)
    tree = cKDTree(xyz)
    d, i = tree.query(pts, k=1)
    clear = d - rad[i]
    free = clear > PROBE
    occ = (clear < 0).reshape(shape)
    cand = np.where(free)[0]
    cp = pts[cand]
    steps = np.arange(1.0, 15.0, 0.75)
    hits = np.zeros(len(cp))
    shp = np.array(shape)
    for dv in DIRS:
        blocked = np.zeros(len(cp), bool)
        for s in steps:
            q = cp + dv * s
            idx = ((q - lo) / SPACING).astype(int)
            ok = np.all((idx >= 0) & (idx < shp), axis=1)
            b = np.zeros(len(cp), bool)
            v = idx[ok]
            b[ok] = occ[v[:, 0], v[:, 1], v[:, 2]]
            blocked |= b
        hits += blocked
    mask = np.zeros(len(pts), bool)
    mask[cand[(hits / len(DIRS)) >= BUR_CUT]] = True
    lab, n = ndimage.label(mask.reshape(shape))
    ci = np.clip(((REF - lo) / SPACING).astype(int), 0, shp - 1)
    l = lab[ci[0], ci[1], ci[2]]
    if l == 0:
        nz = np.array(np.nonzero(mask.reshape(shape))).T
        l = lab[tuple(nz[np.argmin(np.linalg.norm(nz - ci, axis=1))])]
    sel = np.array(np.nonzero(lab == l)).T
    return sel * SPACING + lo


# residue -> heavy atoms
res_atoms, resname = {}, {}
for line in open(PDB):
    if not line.startswith("ATOM"):
        continue
    e = (line[76:78].strip() or line[12:16].strip()[0])
    if e == "H":
        continue
    rn = int(line[22:26])
    resname[rn] = line[17:20].strip()
    res_atoms.setdefault(rn, []).append(
        (line[12:16].strip(), np.array([float(line[30:38]), float(line[38:46]),
                                        float(line[46:54])]), e))

lig_xyz, _ = atoms_from_pose_pdb(LIG, exclude_resnames=())
REF = lig_xyz.mean(0)

QUAD = (59, 108, 79, 94)
wt_pts = cavity_points()
ex_pts = cavity_points(QUAD)
print(f"WT cavity      : {len(wt_pts)} grid points = {len(wt_pts)*SPACING**3:.1f} A^3")
print(f"expanded cavity: {len(ex_pts)} grid points = {len(ex_pts)*SPACING**3:.1f} A^3")
print(f"   (expanded = K59A/F108A/R79A/E94A in-silico truncation)\n")

wt_tree, ex_tree = cKDTree(wt_pts), cKDTree(ex_pts)
lig_min = {rn: min(float(np.linalg.norm(lig_xyz - a[1], axis=1).min())
                   for a in ats) for rn, ats in res_atoms.items()}

rows = []
for rn, ats in res_atoms.items():
    xyz = np.array([a[1] for a in ats])
    dwt, _ = wt_tree.query(xyz, k=1)
    dex, _ = ex_tree.query(xyz, k=1)
    n_wt = int(sum(len(wt_tree.query_ball_point(p, LINING_CUT)) for p in xyz))
    n_ex = int(sum(len(ex_tree.query_ball_point(p, LINING_CUT)) for p in xyz))
    rows.append(dict(resnum=rn, aa=resname[rn],
                     min_dist_cavity_wt=round(float(dwt.min()), 2),
                     min_dist_cavity_expanded=round(float(dex.min()), 2),
                     lining_pts_wt=n_wt, lining_pts_expanded=n_ex,
                     min_dist_ABA=round(lig_min[rn], 2)))

lining = [r for r in rows if r["lining_pts_wt"] > 0 or r["lining_pts_expanded"] > 0]
lining.sort(key=lambda r: -r["lining_pts_wt"])

print("=== cavity-lining residues (grid points within 4.5 A of a heavy atom) ===")
print(f"{'res':<9}{'lining_WT':>10}{'lining_EXP':>11}{'d_cav_WT':>10}"
      f"{'d_cav_EXP':>11}{'d_ABA':>8}")
for r in lining[:30]:
    print(f"{r['aa']}{r['resnum']:<6}{r['lining_pts_wt']:>10}"
          f"{r['lining_pts_expanded']:>11}{r['min_dist_cavity_wt']:>10.2f}"
          f"{r['min_dist_cavity_expanded']:>11.2f}{r['min_dist_ABA']:>8.2f}")

print("\n=== the six panel positions ===")
print(f"{'res':<9}{'lining_WT':>10}{'lining_EXP':>11}{'d_cav_WT':>10}"
      f"{'d_cav_EXP':>11}{'d_ABA':>8}")
for p in (59, 79, 94, 108, 120, 141):
    r = next(x for x in rows if x["resnum"] == p)
    print(f"{r['aa']}{r['resnum']:<6}{r['lining_pts_wt']:>10}"
          f"{r['lining_pts_expanded']:>11}{r['min_dist_cavity_wt']:>10.2f}"
          f"{r['min_dist_cavity_expanded']:>11.2f}{r['min_dist_ABA']:>8.2f}")

# ---- line-of-sight test ----
print("\n=== line-of-sight: is anything between R79 NH1/NH2 and ABA? ===")
all_xyz = np.array([a[1] for ats in res_atoms.values() for a in ats])
all_rad = np.array([VDW.get(a[2], 1.70) for ats in res_atoms.values() for a in ats])
all_lab = [(rn, a[0]) for rn, ats in res_atoms.items() for a in ats]

for probe_res, probe_atom in ((79, "NH1"), (79, "NH2"), (79, "NE"), (94, "OE2")):
    src = next(a[1] for a in res_atoms[probe_res] if a[0] == probe_atom)
    j = int(np.linalg.norm(lig_xyz - src, axis=1).argmin())
    dst = lig_xyz[j]
    seg = dst - src
    L = float(np.linalg.norm(seg))
    u = seg / L
    ts = np.arange(0.3, L - 0.3, 0.1)
    pts = src + np.outer(ts, u)
    blockers = {}
    for k, (c, rad) in enumerate(zip(all_xyz, all_rad)):
        rnum, aname = all_lab[k]
        if rnum == probe_res:
            continue
        d = np.linalg.norm(pts - c, axis=1)
        if (d < rad).any():
            blockers[(rnum, aname)] = float(d.min())
    tag = f"{resname[probe_res]}{probe_res} {probe_atom}"
    if blockers:
        b = sorted(blockers.items(), key=lambda x: x[1])[:4]
        print(f"  {tag:<14} -> ABA ({L:.2f} A): BLOCKED by "
              + ", ".join(f"{resname[r]}{r}:{a}" for (r, a), _ in b))
    else:
        print(f"  {tag:<14} -> ABA ({L:.2f} A): CLEAR -- no intervening "
              f"PYR1 atom; open cavity space")

out = os.path.join(ROOT, "results", "04_cavity_lining.csv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(sorted(rows, key=lambda r: -r["lining_pts_wt"]))
print(f"\nwrote {out}")
