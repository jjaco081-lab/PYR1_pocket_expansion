#!/usr/bin/env python
"""
18_exits_and_bottleneck.py -- two design questions for oversized ligands.

Q1. WHERE DOES A PROTRUDING TAIL EXIT, AND DOES IT HIT THE READOUT?
    A C40 carotenoid cannot be fully enclosed (README section 4c; BmCBP orders
    only half of its zeaxanthin). So the tail must leave the protein somewhere.
    If it exits through the gate/latch groove -- the surface HAB1 reads, where
    W385 wedges in -- the sensor is dead regardless of binding affinity. If it
    exits on the opposite face, the readout may survive.

    This script finds the cavity's MOUTHS (points bordering the enclosed cavity
    that are solvent-open rather than buried), clusters them, and reports for
    each mouth:
      - its distance to HAB1 (chain B of 3QN1) and specifically to W385
      - its distance to the gate (85-89) and latch (115-117)
      - the residues lining it
    A mouth far from both HAB1 and the gate/latch is a viable tail exit.

Q2. WHAT LINES THE BOTTLENECK, AND WHAT DO CAROTENOID-BINDING RELATIVES PUT
    THERE?
    The merged cavity is pinched (r_min ~2.1 A, README section 4c). This
    identifies the residues forming that constriction, so they can be targeted
    directly, and reports what the superfamily has at the aligned positions.

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
HAB1 = os.path.join(ROOT, "data", "hab1_B.pdb")
SPACING, PROBE, BUR_CUT, BOX = 0.5, 1.4, 0.88, 16.0
QUAD = (59, 108, 79, 94)
GATE, LATCH = list(range(85, 90)), list(range(115, 118))

lig_xyz, _ = atoms_from_pose_pdb(LIG, exclude_resnames=())
REF = lig_xyz.mean(0)

# protein residues (for lining lookups)
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

# HAB1
hab_xyz, _ = atoms_from_pose_pdb(HAB1, chain="B")
w385 = []
for line in open(HAB1):
    if line.startswith("ATOM") and int(line[22:26]) == 385:
        e = (line[76:78].strip() or line[12:16].strip()[0])
        if e != "H":
            w385.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
w385 = np.array(w385)


def build(truncate):
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
    free = clear > PROBE
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
    bur = np.zeros(len(pts)); bur[cand] = hits / len(DIRS)
    mask = np.zeros(len(pts), bool); mask[cand[(hits / len(DIRS)) >= BUR_CUT]] = True
    m3 = mask.reshape(shape)
    lab, n = ndimage.label(m3)
    ci = np.clip(((REF - lo) / SPACING).astype(int), 0, shp - 1)
    l = lab[ci[0], ci[1], ci[2]]
    if l == 0:
        nz = np.array(np.nonzero(m3)).T
        l = lab[tuple(nz[np.argmin(np.linalg.norm(nz - ci, axis=1))])]
    cav = (lab == l).reshape(-1)
    return pts, cav, free, bur, shape, lo, tree, rad


def lining(points, cut=5.0):
    out = collections.Counter()
    for rn, ats in res_atoms.items():
        A = np.array(ats)
        d = np.linalg.norm(points[:, None, :] - A[None, :, :], axis=-1).min(1)
        n = int((d <= cut).sum())
        if n:
            out[rn] = n
    return out


print("=" * 84)
print("Q1. MOUTHS OF THE EXPANDED CAVITY -- where would a protruding tail exit?")
print("=" * 84)
pts, cav, free, bur, shape, lo, tree, rad = build(QUAD)
cavpts = pts[cav]
print(f"expanded cavity (K59A/F108A/R79A/E94A): {len(cavpts)} pts = "
      f"{len(cavpts)*SPACING**3:.1f} A^3")

# mouth = free, low-buriedness points adjacent to the cavity component
cav3 = cav.reshape(shape)
dil = ndimage.binary_dilation(cav3, iterations=3).reshape(-1)
mouth = dil & free & (bur < 0.75) & (~cav)
mpts = pts[mouth]
print(f"mouth points (free, buriedness<0.75, within 1.5 A of cavity): {len(mpts)}")

if len(mpts):
    # cluster mouths
    mlab, mn = ndimage.label(mouth.reshape(shape))
    sizes = ndimage.sum(mouth.reshape(shape), mlab, range(1, mn + 1))
    order = np.argsort(-sizes)[:5]
    hab_tree = cKDTree(hab_xyz)
    w_tree = cKDTree(w385)
    gate_xyz = np.array([a for rn in GATE for a in res_atoms.get(rn, [])])
    latch_xyz = np.array([a for rn in LATCH for a in res_atoms.get(rn, [])])
    g_tree, l_tree = cKDTree(gate_xyz), cKDTree(latch_xyz)
    print(f"\n{'mouth':<7}{'pts':>6}{'vol A3':>9}{'->HAB1':>9}{'->W385':>9}"
          f"{'->gate':>9}{'->latch':>9}   verdict")
    for k in order:
        sel = np.array(np.nonzero(mlab == k + 1)).T * SPACING + lo
        if len(sel) < 8:
            continue
        dh = float(hab_tree.query(sel)[0].min())
        dw = float(w_tree.query(sel)[0].min())
        dg = float(g_tree.query(sel)[0].min())
        dl = float(l_tree.query(sel)[0].min())
        if dh < 5 or dw < 6:
            v = "BLOCKED -- opens into the HAB1 face"
        elif dg < 5 or dl < 5:
            v = "RISKY -- adjacent to gate/latch"
        else:
            v = "VIABLE tail exit (far from readout)"
        print(f"{k+1:<7}{len(sel):>6}{len(sel)*SPACING**3:>9.1f}"
              f"{dh:>9.2f}{dw:>9.2f}{dg:>9.2f}{dl:>9.2f}   {v}")
        top = lining(sel).most_common(6)
        print("        lined by: " + ", ".join(
            f"{resname[r]}{r}" for r, _ in top))

print("\n" + "=" * 84)
print("Q2. THE BOTTLENECK -- which residues pinch the merged channel?")
print("=" * 84)
# span axis from 16_cavity_shape logic, recomputed cheaply on extremes
c = cavpts.mean(0)
order = np.argsort(-np.linalg.norm(cavpts - c, axis=1))[:900]
S = cavpts[order]


def clearance(P):
    d, i = tree.query(P, k=1)
    return d - rad[i]


best, pair = 0.0, None
for a in range(len(S)):
    d = np.linalg.norm(S - S[a], axis=1)
    for b in np.where(d > best)[0]:
        L = d[b]
        t = np.linspace(0, 1, max(int(L / 0.4), 2))[:, None]
        if clearance(S[a] * (1 - t) + S[b] * t).min() >= 1.6:
            best, pair = float(L), (S[a], S[b])
print(f"span = {best:.1f} A")
a, b = pair
u = (b - a) / np.linalg.norm(b - a)
t = (cavpts - a) @ u
perp = np.linalg.norm((cavpts - a) - np.outer(t, u), axis=1)
bins = np.linspace(t.min(), t.max(), 16)
prof = []
for k in range(len(bins) - 1):
    m = (t >= bins[k]) & (t < bins[k + 1])
    if m.sum() > 3:
        prof.append((float(perp[m].max()), 0.5 * (bins[k] + bins[k + 1]), m))
prof.sort(key=lambda x: x[0])
print("\nnarrowest three slices along the span axis:")
allb = collections.Counter()
for r, tc, m in prof[:3]:
    sl = cavpts[m]
    top = lining(sl, cut=5.5).most_common(8)
    for rn, n in top:
        allb[rn] += n
    print(f"  radius {r:.2f} A at t={tc:+.1f} A -- lined by: "
          + ", ".join(f"{resname[rn]}{rn}" for rn, _ in top))
print("\nBOTTLENECK RESIDUES (union, ranked):")
for rn, n in allb.most_common(12):
    tag = "  <- already in panel" if rn in (59, 79, 94, 108, 120, 141) else ""
    print(f"   {resname[rn]}{rn:<5} weight {n}{tag}")
print("\nThese are the positions to target for WIDENING, as distinct from the")
print("volume-opening positions already in the panel.")
