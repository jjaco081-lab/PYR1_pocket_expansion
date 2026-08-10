#!/usr/bin/env python
"""
19_interlobe_neck.py -- the constriction BETWEEN the two lobes, and how far a
C40 carotenoid still is from fitting.

Script 18 profiled the narrowest slices of the merged channel and found them at
one END of the span axis, lined by gate/latch residues. That answered "where
would a protruding tail exit" but NOT the question that matters for full
enclosure: where is the neck BETWEEN the ABA lobe and the newly opened satellite
lobe, what lines it, and can it be widened without touching the readout?

This script:
  1. profiles cross-sectional radius along the full span axis in fine bins, so
     every local constriction is visible, not just the global minimum;
  2. locates the INTER-LOBE neck as the constriction lying between the two
     volume maxima (the two lobe centres), rather than at either end;
  3. lists the residues lining that neck, flagging which are gate/latch/sealing
     (untouchable) versus free to mutate;
  4. tests whether truncating the free ones widens the neck and lengthens the
     span -- i.e. whether full enclosure of a longer ligand is reachable without
     touching the readout;
  5. reports the residual gap to an intact C40 carotenoid.

Run with the esmfold2 env python.
"""
import os, sys, warnings, collections
warnings.filterwarnings("ignore")
import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_cavity import VDW, DIRS, atoms_from_pose_pdb  # noqa: E402

PDB = os.path.join(ROOT, "data", "pyr1_A.pdb")
LIG = os.path.join(ROOT, "data", "aba_xtal.pdb")
SPACING, PROBE, BUR_CUT, BOX = 0.5, 1.4, 0.88, 17.0
QUAD = (59, 108, 79, 94)
GATE, LATCH = set(range(85, 90)), set(range(115, 118))
SEALING = {61, 87, 159}          # negative-dV residues from script 02
UNTOUCHABLE = GATE | LATCH | SEALING

lig_xyz, _ = atoms_from_pose_pdb(LIG, exclude_resnames=())
REF = lig_xyz.mean(0)

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
        np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))


def cavity(truncate=()):
    xyz, el = atoms_from_pose_pdb(PDB, chain="A", truncate_to_ala=set(truncate))
    rad = np.array([VDW.get(e, 1.70) for e in el])
    lo = REF - BOX
    axes = [np.arange(lo[i], REF[i] + BOX, SPACING) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    shape = grid.shape[:3]
    pts = grid.reshape(-1, 3)
    tree = cKDTree(xyz)
    d, i = tree.query(pts, k=1)
    clear = d - rad[i]
    occ = (clear < 0).reshape(shape)
    cand = np.where(clear > PROBE)[0]
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
    if n == 0:
        return np.empty((0, 3)), tree, rad
    ci = np.clip(((REF - lo) / SPACING).astype(int), 0, shp - 1)
    l = lab[ci[0], ci[1], ci[2]]
    if l == 0:
        nz = np.array(np.nonzero(mask.reshape(shape))).T
        l = lab[tuple(nz[np.argmin(np.linalg.norm(nz - ci, axis=1))])]
    return np.array(np.nonzero(lab == l)).T * SPACING + lo, tree, rad


def span_axis(P, tree, rad, npr=900):
    def clr(X):
        d, i = tree.query(X, k=1)
        return d - rad[i]
    c = P.mean(0)
    S = P[np.argsort(-np.linalg.norm(P - c, axis=1))[:npr]]
    best, pair = 0.0, None
    for a in range(len(S)):
        d = np.linalg.norm(S - S[a], axis=1)
        for b in np.where(d > best)[0]:
            L = d[b]
            t = np.linspace(0, 1, max(int(L / 0.4), 2))[:, None]
            if clr(S[a] * (1 - t) + S[b] * t).min() >= 1.6:
                best, pair = float(L), (S[a], S[b])
    return best, pair


def profile(P, pair, nbin=36):
    a, b = pair
    u = (b - a) / np.linalg.norm(b - a)
    t = (P - a) @ u
    perp = np.linalg.norm((P - a) - np.outer(t, u), axis=1)
    edges = np.linspace(t.min(), t.max(), nbin + 1)
    out = []
    for k in range(nbin):
        m = (t >= edges[k]) & (t < edges[k + 1])
        if m.sum() >= 3:
            out.append((0.5 * (edges[k] + edges[k + 1]),
                        float(perp[m].max()), int(m.sum()), m))
    return out, u, t


def lining(pts, cut=5.5):
    c = collections.Counter()
    for rn, ats in res_atoms.items():
        A = np.array(ats)
        d = np.linalg.norm(pts[:, None, :] - A[None, :, :], axis=-1).min(1)
        n = int((d <= cut).sum())
        if n:
            c[rn] = n
    return c


print("=" * 86)
print("FULL CROSS-SECTIONAL PROFILE of the merged (quad) cavity along its span")
print("=" * 86)
P, tree, rad = cavity(QUAD)
span, pair = span_axis(P, tree, rad)
prof, u, t = profile(P, pair)
print(f"volume {len(P)*SPACING**3:.1f} A^3,  span {span:.1f} A\n")
print(f"{'t (A)':>8}{'radius':>9}{'n pts':>7}   profile")
for tc, r, n, m in prof:
    print(f"{tc:>8.1f}{r:>9.2f}{n:>7}   " + "#" * int(r * 7))

# lobe centres = the two largest-area maxima; neck = minimum BETWEEN them
areas = np.array([p[2] for p in prof])
radii = np.array([p[1] for p in prof])
tc = np.array([p[0] for p in prof])
i_max1 = int(np.argmax(areas))
# find the other lobe: largest area at least 6 A away
far = np.where(np.abs(tc - tc[i_max1]) > 6.0)[0]
i_max2 = far[int(np.argmax(areas[far]))] if len(far) else i_max1
lo_i, hi_i = sorted((i_max1, i_max2))
mid = range(lo_i + 1, hi_i)
if len(list(mid)) == 0:
    print("\n(no interior bins between the two lobe centres)")
    sys.exit(0)
i_neck = min(mid, key=lambda k: radii[k])

print(f"\nlobe centre A at t={tc[lo_i]:+.1f} A (r={radii[lo_i]:.2f}, {areas[lo_i]} pts)")
print(f"lobe centre B at t={tc[hi_i]:+.1f} A (r={radii[hi_i]:.2f}, {areas[hi_i]} pts)")
print(f"INTER-LOBE NECK at t={tc[i_neck]:+.1f} A, radius {radii[i_neck]:.2f} A")

neck_pts = P[prof[i_neck][3]]
print("\nresidues lining the INTER-LOBE neck:")
free_targets = []
for rn, n in lining(neck_pts).most_common(12):
    if rn in GATE:
        tag = "GATE -- untouchable"
    elif rn in LATCH:
        tag = "LATCH -- untouchable"
    elif rn in SEALING:
        tag = "sealing -- do not remove"
    elif rn in QUAD:
        tag = "already truncated in quad"
    else:
        tag = "FREE -- widening target"
        free_targets.append(rn)
    print(f"   {resname[rn]}{rn:<5} weight {n:<4} {tag}")

print("\n" + "=" * 86)
print("CAN THE NECK BE WIDENED WITHOUT TOUCHING THE READOUT?")
print("=" * 86)
tests = [("quad (reference)", QUAD)]
if free_targets:
    tests.append((f"quad + neck({','.join(map(str,free_targets[:4]))})",
                  tuple(QUAD) + tuple(free_targets[:4])))
tests.append(("quad + Y120A/E141A", tuple(QUAD) + (120, 141)))
if free_targets:
    tests.append(("quad + Y120A/E141A + neck",
                  tuple(QUAD) + (120, 141) + tuple(free_targets[:4])))

C40_LEN, C40_HALFW = 27.8, 3.35     # intact zeaxanthin, MMFF (script 17)
print(f"{'variant':<34}{'vol A3':>9}{'span A':>9}{'neck r':>9}   gap to C40")
for name, tr in tests:
    Q, tr_tree, tr_rad = cavity(tr)
    if len(Q) == 0:
        print(f"{name:<34}   none"); continue
    s, pr = span_axis(Q, tr_tree, tr_rad)
    pf, _, _ = profile(Q, pr)
    rr = np.array([p[1] for p in pf]); aa = np.array([p[2] for p in pf])
    tt = np.array([p[0] for p in pf])
    j1 = int(np.argmax(aa))
    fr = np.where(np.abs(tt - tt[j1]) > 6.0)[0]
    j2 = fr[int(np.argmax(aa[fr]))] if len(fr) else j1
    a_, b_ = sorted((j1, j2))
    inner = list(range(a_ + 1, b_))
    nr = float(rr[min(inner, key=lambda k: rr[k])]) if inner else float("nan")
    print(f"{name:<34}{len(Q)*SPACING**3:>9.1f}{s:>9.1f}{nr:>9.2f}"
          f"   {C40_LEN - s:+6.1f} A")

print(f"""
Reference lengths (script 17, MMFF geometries):
   ABA                     11.4 A
   retinoic acid           14.5 A
   crocetin                20.7 A
   beta-apo-8'-carotenal   24.0 A
   zeaxanthin (C40)        {C40_LEN:.1f} A, half-width {C40_HALFW:.2f} A
   beta-carotene (C40)     27.3 A

A C40 needs BOTH span >= ~28 A and neck radius >= ~3.4 A along its whole length.
Note the span metric requires 1.6 A clearance the whole way, so it is a lower
bound on what a snug ligand could thread; but the neck radius is the hard
constraint -- a polyene cannot squeeze through a hole narrower than itself.""")
