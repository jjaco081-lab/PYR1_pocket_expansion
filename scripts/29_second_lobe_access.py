#!/usr/bin/env python
"""
29_second_lobe_access.py -- how much does each variant actually open the SECOND
LOBE, and what is the maximally open pocket we could design into?

Two questions, both from the user (2026-08-07).

Q1. SECOND-LOBE ACCESS.
    "Bigger cavity" is not the quantity of interest. The satellite lobe is a
    separate, sealed void in WT PYR1; merging it with the ABA site is what
    recruits lining residues (I48, V49, R50, V81, I62, S122) that were never in
    the Tian or Mosquna libraries. So the decision metric is not dV, it is: does
    this variant actually reach the second lobe, and how much of it?

    Method: locate the satellite lobe in WT as a connected cavity component that
    is NOT the ABA-seeded one, and cache its voxels. Then for each relaxed
    variant structure, take the ABA-seeded component and ask
      (a) MERGED?      does it now contain the WT satellite lobe's centre
      (b) V_lobe2      how much of its volume sits inside the WT satellite region
      (c) frac_lobe2   V_lobe2 as a fraction of the variant's total cavity
    A variant with a large dV but frac_lobe2 ~ 0 has just hollowed out the
    existing chamber -- exactly the outcome that adds no new lining residues.

Q2. THE MAXIMALLY OPEN POCKET.
    Even though shrinking Tian-library residues is not a deliverable in itself
    (those substitutions are already inside the existing oligo pools), a maximally
    open backbone is useful as a DESIGN STARTING POINT: build the biggest
    plausible cavity, then let a sequence designer shrink it back to complement a
    specific ligand. This computes cavity volume for progressively more aggressive
    truncation sets, always excluding gate (85-89) and latch (115-117).

    Note that truncating SEALING residues (61, 87, 159; script 02) can REDUCE
    enclosed volume by opening the cavity to bulk solvent, so sets are reported
    with and without them.

Relaxed inputs are results/cavity_scan/pdb/<variant>_r0.pdb from the Arm 1 pilot
(FastRelax, CA-constrained, so the frame matches the WT reference).

Run with the esmfold2 env python.
"""
import os, sys, glob, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_cavity import VDW, DIRS, atoms_from_pose_pdb  # noqa: E402

WT_PDB = os.path.join(ROOT, "data", "pyr1_A.pdb")
LIG = os.path.join(ROOT, "data", "aba_xtal.pdb")
PDB_DIR = os.path.join(ROOT, "results", "cavity_scan", "pdb")
OUT = os.path.join(ROOT, "results", "29_second_lobe_access.csv")

SP, PROBE, BUR, BOX = 0.5, 1.4, 0.88, 16.0
GATE, LATCH = set(range(85, 90)), set(range(115, 118))
SEALING = {61, 87, 159}
# Tian's 18 randomised positions (script 11)
TIAN = [59, 81, 83, 87, 89, 92, 94, 108, 110, 117, 120, 122, 141, 159, 160,
        163, 164, 167]

lig_xyz, _ = atoms_from_pose_pdb(LIG, exclude_resnames=())
REF = lig_xyz.mean(0)
LO = REF - BOX
AXES = [np.arange(LO[i], REF[i] + BOX, SP) for i in range(3)]
GRID = np.stack(np.meshgrid(*AXES, indexing="ij"), -1)
SHAPE = GRID.shape[:3]
PTS = GRID.reshape(-1, 3)
SHP = np.array(SHAPE)
STEPS = np.arange(1.0, 15.0, 0.75)


def components(xyz, elem):
    """Label all enclosed-cavity components on the shared grid."""
    rad = np.array([VDW.get(e, 1.70) for e in elem])
    tree = cKDTree(xyz)
    d, i = tree.query(PTS, k=1)
    clear = d - rad[i]
    occ = (clear < 0).reshape(SHAPE)
    cand = np.where(clear > PROBE)[0]
    if cand.size == 0:
        return np.zeros(SHAPE, int), 0
    cp = PTS[cand]
    hits = np.zeros(len(cp))
    for dv in DIRS:
        blocked = np.zeros(len(cp), bool)
        for s in STEPS:
            q = cp + dv * s
            idx = ((q - LO) / SP).astype(int)
            ok = np.all((idx >= 0) & (idx < SHP), axis=1)
            b = np.zeros(len(cp), bool)
            v = idx[ok]
            b[ok] = occ[v[:, 0], v[:, 1], v[:, 2]]
            blocked |= b
        hits += blocked
    mask = np.zeros(len(PTS), bool)
    mask[cand[(hits / len(DIRS)) >= BUR]] = True
    return ndimage.label(mask.reshape(SHAPE))


def seeded(lab, n):
    """Index of the component containing (or nearest to) the ABA centroid."""
    ci = np.clip(((REF - LO) / SP).astype(int), 0, SHP - 1)
    k = lab[ci[0], ci[1], ci[2]]
    if k == 0:
        nz = np.array(np.nonzero(lab)).T
        if not len(nz):
            return 0
        k = lab[tuple(nz[np.argmin(np.linalg.norm(nz - ci, axis=1))])]
    return int(k)


# ---------- define the WT satellite lobe ----------
xyz, el = atoms_from_pose_pdb(WT_PDB, chain="A")
lab, n = components(xyz, el)
k_aba = seeded(lab, n)
sizes = ndimage.sum(lab > 0, lab, range(1, n + 1))
print("=" * 78)
print("WT CAVITY COMPONENTS (grid seeded at the ABA centroid)")
print("=" * 78)
wall = []
for rn in (79, 94, 108, 59):
    for line in open(WT_PDB):
        if line.startswith("ATOM") and int(line[22:26]) == rn:
            e = (line[76:78].strip() or line[12:16].strip()[0])
            if e != "H":
                wall.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
wall = np.array(wall)
wall_tree = cKDTree(wall)

cands = []
for k in range(1, n + 1):
    v = float(sizes[k - 1] * SP ** 3)
    if v < 8.0:
        continue
    pts = np.array(np.nonzero(lab == k)).T * SP + LO
    dwall = float(wall_tree.query(pts)[0].min())
    tag = "ABA-seeded (lobe 1)" if k == k_aba else ""
    cands.append((k, v, dwall, pts.mean(0), tag))
    print(f"  component {k:>2}: {v:7.1f} A^3   min dist to K59/R79/E94/F108 = "
          f"{dwall:5.2f} A   {tag}")

# satellite = the non-ABA component closest to the wall cluster
sat = [c for c in cands if c[0] != k_aba]
if not sat:
    raise SystemExit("no satellite lobe found -- adjust BOX or buriedness cut")
sat.sort(key=lambda c: c[2])
K_SAT, V_SAT, D_SAT, C_SAT, _ = sat[0]
SAT_MASK = ndimage.binary_dilation(lab == K_SAT, iterations=2)
print(f"\nSATELLITE LOBE = component {K_SAT}: {V_SAT:.1f} A^3, centre "
      f"{np.round(C_SAT,1)}, {D_SAT:.2f} A from the wall cluster")
c_idx = np.clip(((C_SAT - LO) / SP).astype(int), 0, SHP - 1)


def lobe_metrics(path):
    x, e = atoms_from_pose_pdb(path, chain="A")
    lb, nn = components(x, e)
    if nn == 0:
        return 0.0, 0.0, 0.0, False
    k = seeded(lb, nn)
    comp = (lb == k)
    tot = float(comp.sum() * SP ** 3)
    v2 = float((comp & SAT_MASK).sum() * SP ** 3)
    merged = bool(comp[c_idx[0], c_idx[1], c_idx[2]])
    return tot, v2, (v2 / tot if tot else 0.0), merged


# ---------- Q1: per-variant second-lobe access ----------
print("\n" + "=" * 78)
print("Q1. SECOND-LOBE ACCESS, per relaxed variant (rep 0)")
print("=" * 78)
rows = []
for f in sorted(glob.glob(os.path.join(PDB_DIR, "*_r0.pdb"))):
    name = os.path.basename(f)[:-len("_r0.pdb")]
    tot, v2, frac, merged = lobe_metrics(f)
    rows.append(dict(variant=name, cavity_A3=tot, lobe2_A3=v2,
                     frac_lobe2=frac, merged=merged))

wt = next((r for r in rows if r["variant"] == "WT"), None)
rows.sort(key=lambda r: -r["lobe2_A3"])
print(f"{'variant':<26}{'cavity':>9}{'lobe2':>9}{'frac':>7}{'merged':>8}")
for r in rows:
    print(f"{r['variant']:<26}{r['cavity_A3']:>9.1f}{r['lobe2_A3']:>9.1f}"
          f"{r['frac_lobe2']:>7.2f}{'YES' if r['merged'] else '-':>8}")

with open(OUT, "w") as fh:
    fh.write("variant,cavity_A3,lobe2_A3,frac_lobe2,merged\n")
    for r in rows:
        fh.write(f"{r['variant']},{r['cavity_A3']:.2f},{r['lobe2_A3']:.2f},"
                 f"{r['frac_lobe2']:.4f},{int(r['merged'])}\n")
print(f"\nwrote {OUT}")
if wt:
    print(f"WT reference: cavity {wt['cavity_A3']:.1f}, lobe2 {wt['lobe2_A3']:.1f}, "
          f"merged={wt['merged']}")

# ---------- Q2: maximally open backbone ----------
print("\n" + "=" * 78)
print("Q2. MAXIMALLY OPEN POCKET (rigid truncation, gate/latch never touched)")
print("=" * 78)
tian_ok = [p for p in TIAN if p not in GATE and p not in LATCH]
tian_nosealing = [p for p in tian_ok if p not in SEALING]
SETS = [
    ("WT", []),
    ("quad K59/F108/R79/E94", [59, 108, 79, 94]),
    ("quad + Y120 + E141", [59, 108, 79, 94, 120, 141]),
    ("Tian pocket set (no gate/latch)", tian_ok),
    ("Tian set, sealing kept", tian_nosealing),
    ("Tian set + R79 (never sampled)", tian_nosealing + [79]),
    ("Tian set + R79/E94 + lobe2 lining",
     sorted(set(tian_nosealing + [79, 94, 48, 49, 50, 62, 122]))),
]
print(f"{'truncation set':<38}{'n':>4}{'cavity':>9}{'lobe2':>9}{'frac':>7}{'merged':>8}")
for label, trunc in SETS:
    x, e = atoms_from_pose_pdb(WT_PDB, chain="A", truncate_to_ala=set(trunc))
    lb, nn = components(x, e)
    if nn == 0:
        print(f"{label:<38}{len(trunc):>4}     none")
        continue
    k = seeded(lb, nn)
    comp = (lb == k)
    tot = float(comp.sum() * SP ** 3)
    v2 = float((comp & SAT_MASK).sum() * SP ** 3)
    mg = bool(comp[c_idx[0], c_idx[1], c_idx[2]])
    print(f"{label:<38}{len(trunc):>4}{tot:>9.1f}{v2:>9.1f}"
          f"{(v2/tot if tot else 0):>7.2f}{'YES' if mg else '-':>8}")

print("""
READING THIS
  lobe2 / frac  = how much of the cavity lies in the WT satellite lobe. This is
                  the quantity that decides whether NEW lining residues are
                  recruited. A variant with a big cavity and frac ~ 0 has only
                  hollowed the existing chamber and adds nothing to the library.
  merged        = the satellite lobe centre is inside the ABA-connected cavity,
                  i.e. one continuous pocket rather than two voids.
  Q2 rows are RIGID truncations (upper bounds; no repacking), for choosing a
  design starting point -- not predictions of what survives relax.""")
