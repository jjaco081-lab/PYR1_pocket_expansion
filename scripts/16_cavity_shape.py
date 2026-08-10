#!/usr/bin/env python
"""
16_cavity_shape.py -- cavity SHAPE, not just volume: can a merged pocket host an
extended polyene?

Rationale
---------
Volume alone does not decide whether a carotenoid fits. An extended C40
xanthophyll such as zeaxanthin is ~30 A end to end but slim (~5-6 A across), so
what matters is the cavity's LONGEST INTERNAL DIMENSION and its cross-section
along that axis -- i.e. is the pocket a tunnel or a ball?

The WT PYR1 cavity is globular, which is why intact carotenoids look
implausible for WT. But this project's whole premise (README section 2) is that
truncating the K59/F108/R79/E94 wall MERGES the ABA pocket with a sealed
satellite lobe, recruiting I48, V49, R50, V81, I62 and S122 as new lining
residues. Merging two lobes produces ELONGATION, not just added volume -- so
the tunnel-vs-ball objection has to be re-tested on the expanded cavity rather
than inherited from WT.

What this computes, for WT and for a series of truncation sets:
  * volume (as elsewhere in this project)
  * principal axes of the cavity point cloud (PCA eigenvalues -> extents)
  * MAX INTERNAL SPAN: the largest pairwise distance between cavity grid points
    that stays inside the cavity (a straight-line path test), which is the
    honest upper bound on the length of a RIGID rod-like ligand
  * the largest inscribed sphere radius (widest point)
  * the cross-sectional radius along the principal axis (how wide the tunnel is
    where the ligand would lie)

Reference lengths printed for comparison: ABA, the modelled half of zeaxanthin
in 7ZVR, and intact zeaxanthin geometry generated from SMILES.

Run with the esmfold2 env python.
"""
import os, sys, warnings
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
SPACING, PROBE, BUR_CUT, BOX = 0.5, 1.4, 0.88, 16.0

lig_xyz, _ = atoms_from_pose_pdb(LIG, exclude_resnames=())
REF = lig_xyz.mean(0)


def cavity_cloud(truncate=()):
    """Return (points, occupancy_grid, lo) for the cavity nearest the ABA seed."""
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
    if cand.size == 0:
        return np.empty((0, 3)), occ, lo, tree, rad
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
    m3 = mask.reshape(shape)
    lab, n = ndimage.label(m3)
    if n == 0:
        return np.empty((0, 3)), occ, lo, tree, rad
    ci = np.clip(((REF - lo) / SPACING).astype(int), 0, shp - 1)
    l = lab[ci[0], ci[1], ci[2]]
    if l == 0:
        nz = np.array(np.nonzero(m3)).T
        l = lab[tuple(nz[np.argmin(np.linalg.norm(nz - ci, axis=1))])]
    sel = np.array(np.nonzero(lab == l)).T * SPACING + lo
    return sel, occ, lo, tree, rad


def clearance(pts, tree, rad):
    d, i = tree.query(pts, k=1)
    return d - rad[i]


def max_internal_span(P, tree, rad, n_probe=1400, need=1.6):
    """Longest straight segment between cavity points whose whole path keeps at
    least `need` A of clearance from protein atoms -- i.e. the longest rigid rod
    that could lie inside the cavity. Sampled over the convex-hull-ish extremes
    for tractability."""
    if len(P) < 2:
        return 0.0, None
    c = P.mean(0)
    # bias sampling toward extremal points, which is where the max span lives
    order = np.argsort(-np.linalg.norm(P - c, axis=1))
    S = P[order[:min(n_probe, len(P))]]
    best, pair = 0.0, None
    for a in range(len(S)):
        d = np.linalg.norm(S - S[a], axis=1)
        cand = np.where(d > best)[0]
        for b in cand:
            L = d[b]
            steps = max(int(L / 0.4), 2)
            t = np.linspace(0, 1, steps)[:, None]
            path = S[a] * (1 - t) + S[b] * t
            if clearance(path, tree, rad).min() >= need:
                best, pair = float(L), (S[a], S[b])
    return best, pair


def axis_profile(P, pair, tree, rad):
    """Cross-sectional radius along the max-span axis."""
    if pair is None:
        return 0.0, 0.0
    a, b = pair
    u = (b - a) / np.linalg.norm(b - a)
    t = (P - a) @ u
    perp = np.linalg.norm((P - a) - np.outer(t, u), axis=1)
    # median and min of the per-slice max perpendicular reach
    bins = np.linspace(t.min(), t.max(), 14)
    radii = []
    for k in range(len(bins) - 1):
        m = (t >= bins[k]) & (t < bins[k + 1])
        if m.sum() > 3:
            radii.append(perp[m].max())
    if not radii:
        return 0.0, 0.0
    return float(np.median(radii)), float(np.min(radii))


SETS = [
    ("WT", ()),
    ("K59A", (59,)),
    ("F108A", (108,)),
    ("K59A/F108A", (59, 108)),
    ("R79A/E94A", (79, 94)),
    ("F108A/R79A/E94A", (108, 79, 94)),
    ("K59A/F108A/R79A/E94A", (59, 108, 79, 94)),
    ("quad + Y120A", (59, 108, 79, 94, 120)),
    ("quad + Y120A/E141A", (59, 108, 79, 94, 120, 141)),
]

print("=" * 92)
print("CAVITY SHAPE: volume, elongation, and the longest rigid ligand that fits")
print("=" * 92)
print(f"{'variant':<24}{'vol A3':>9}{'span A':>9}{'PCA extents (A)':>26}"
      f"{'r_med':>8}{'r_min':>8}{'max R':>8}")
rows = []
for name, tr in SETS:
    P, occ, lo, tree, rad = cavity_cloud(tr)
    if len(P) == 0:
        print(f"{name:<24}      none")
        continue
    vol = len(P) * SPACING ** 3
    cl = clearance(P, tree, rad)
    maxR = float(cl.max())
    Pc = P - P.mean(0)
    ev = np.linalg.eigvalsh(np.cov(Pc.T))[::-1]
    ext = 2 * np.sqrt(np.maximum(ev, 0)) * 2      # ~2 sigma each side
    span, pair = max_internal_span(P, tree, rad)
    rmed, rmin = axis_profile(P, pair, tree, rad)
    rows.append((name, vol, span, ext, rmed, rmin, maxR))
    print(f"{name:<24}{vol:>9.1f}{span:>9.1f}"
          f"{ext[0]:>9.1f}{ext[1]:>8.1f}{ext[2]:>8.1f}"
          f"{rmed:>8.2f}{rmin:>8.2f}{maxR:>8.2f}")

print("""
  span   = longest straight segment inside the cavity keeping >=1.6 A clearance
           (the honest upper bound on a RIGID rod-like ligand's length)
  PCA    = cavity extents along its three principal axes (elongation indicator)
  r_med  = median cross-sectional radius along the span axis
  r_min  = NARROWEST cross-section along that axis -- the true bottleneck
  max R  = largest inscribed sphere radius (widest single point)""")

# ---- reference ligand dimensions ----
print("\n" + "=" * 92)
print("REFERENCE LIGAND DIMENSIONS")
print("=" * 92)


def dims(P, label):
    d = np.linalg.norm(P[:, None] - P[None], axis=-1)
    L = float(d.max())
    c = P.mean(0)
    Pc = P - c
    ev = np.linalg.eigvalsh(np.cov(Pc.T))[::-1]
    u = np.linalg.eigh(np.cov(Pc.T))[1][:, -1]
    t = Pc @ u
    perp = np.linalg.norm(Pc - np.outer(t, u), axis=1)
    print(f"  {label:<34} {len(P):>3} heavy atoms  length {L:6.2f} A  "
          f"max half-width {perp.max():5.2f} A")
    return L, float(perp.max())


aba_L, aba_w = dims(lig_xyz, "ABA (A8S, 3QN1)")

zx_path = os.path.join(ROOT, "data", "7ZVR.cif")
if os.path.exists(zx_path):
    from Bio.PDB import MMCIFParser
    m = MMCIFParser(QUIET=True).get_structure("x", zx_path)[0]
    Z = np.array([a.coord for ch in m for r in ch if r.get_resname() == "ZEX"
                  for a in r if a.element != "H"])
    dims(Z, "zeaxanthin, MODELLED HALF (7ZVR)")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    SM = ("CC1=C(C(CC(C1)O)(C)C)/C=C/C(=C/C=C(\\C)/C=C/C=C(\\C)/C=C/C=C(\\C)"
          "/C=C/C2=C(CC(CC2(C)C)O)C)/C")
    mol = Chem.AddHs(Chem.MolFromSmiles(SM))
    AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
    AllChem.MMFFOptimizeMolecule(mol)
    mol = Chem.RemoveHs(mol)
    Z2 = mol.GetConformer().GetPositions()
    zx_L, zx_w = dims(Z2, "zeaxanthin, INTACT (MMFF, all-trans)")
except Exception as e:
    print(f"  (intact zeaxanthin geometry unavailable: {e})")
    zx_L, zx_w = 30.0, 3.0

print("\n" + "=" * 92)
print("VERDICT")
print("=" * 92)
best = max(rows, key=lambda r: r[2])
print(f"  longest span achieved: {best[2]:.1f} A  ({best[0]})")
print(f"  intact zeaxanthin needs ~{zx_L:.1f} A length and ~{zx_w:.1f} A half-width")
if best[2] >= zx_L:
    print("  => a C40 xanthophyll is geometrically ACCOMMODATED by the merged cavity.")
else:
    print(f"  => still {zx_L - best[2]:.1f} A short of an intact C40 xanthophyll in")
    print("     the CLOSED state with a rigid backbone. Note this protocol cannot")
    print("     open the backbone; carotenoid binding would additionally require")
    print("     backbone/loop remodelling or a partly solvent-exposed tail, as is")
    print("     in fact seen in BmCBP where half the ligand is disordered.")
