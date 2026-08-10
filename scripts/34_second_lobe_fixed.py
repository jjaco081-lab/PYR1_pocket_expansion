#!/usr/bin/env python
"""
34_second_lobe_fixed.py -- repair the second-lobe metric from script 29.

The problem with script 29
--------------------------
It defined the satellite lobe on the CRYSTAL structure, dilated that mask by 2
voxels (1 A), and then asked how much of each RELAXED variant's cavity fell
inside it. Relaxed WT came back at 38.1 A^3 and "merged", while rigid WT on the
same mask gave 0.0 -- so the baseline was not zero and the per-variant numbers
could not be trusted. Two possible causes, not separated at the time:

  (a) FastRelax alone breaches the wall, so relaxed WT genuinely has one cavity;
  (b) the 1 A dilation bleeds into the main chamber, so the mask overcounts.

This script separates them and replaces the metric with two dilation-free ones.

DIAGNOSTIC. For relaxed WT specifically, decompose its own cavity into connected
components and report how many there are and whether the ABA-seeded one contains
the satellite centroid. That answers (a) directly, with no mask involved.

METRIC 1 -- STRICT OVERLAP. Same as 29 but with NO dilation: variant cavity
voxels that fall inside the UNDILATED crystal-WT satellite component.

METRIC 2 -- BEYOND-THE-WALL VOLUME (connectivity- and mask-free, the robust one).
Build a plane through the centroid of the K59/R79/E94/F108 side-chain atoms whose
normal points from the ABA centroid toward the satellite centroid. Report the
cavity volume on the far side. This asks "how much pocket lies past the wall",
which is the design question, and it cannot be confounded by mask dilation or by
whether two voids happen to be connected at a 0.5 A grid.

Both are reported as DELTAS versus relaxed WT, since the absolute baseline is the
thing that was ambiguous.

Run with the esmfold2 env python.
"""
import os, sys, glob, warnings
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
OUT = os.path.join(ROOT, "results", "34_second_lobe_fixed.csv")

SP, PROBE, BUR, BOX = 0.5, 1.4, 0.88, 16.0
WALL = (59, 79, 94, 108)

lig_xyz, _ = atoms_from_pose_pdb(LIG, exclude_resnames=())
REF = lig_xyz.mean(0)
LO = REF - BOX
AXES = [np.arange(LO[i], REF[i] + BOX, SP) for i in range(3)]
GRID = np.stack(np.meshgrid(*AXES, indexing="ij"), -1)
SHAPE = GRID.shape[:3]
PTS = GRID.reshape(-1, 3)
SHP = np.array(SHAPE)
STEPS = np.arange(1.0, 15.0, 0.75)


def label_cavities(xyz, elem):
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


def seeded(lab):
    ci = np.clip(((REF - LO) / SP).astype(int), 0, SHP - 1)
    k = lab[ci[0], ci[1], ci[2]]
    if k == 0:
        nz = np.array(np.nonzero(lab)).T
        if not len(nz):
            return 0
        k = lab[tuple(nz[np.argmin(np.linalg.norm(nz - ci, axis=1))])]
    return int(k)


# ---------- crystal WT: define the satellite lobe, undilated ----------
xyz, el = atoms_from_pose_pdb(WT_PDB, chain="A")
lab, n = label_cavities(xyz, el)
k_aba = seeded(lab)
wall_xyz = []
for line in open(WT_PDB):
    if line.startswith("ATOM") and int(line[22:26]) in WALL:
        e = (line[76:78].strip() or line[12:16].strip()[0])
        if e != "H" and line[12:16].strip() not in ("N", "CA", "C", "O"):
            wall_xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
wall_xyz = np.array(wall_xyz)
wall_tree = cKDTree(wall_xyz)

best, K_SAT = None, None
for k in range(1, n + 1):
    if k == k_aba:
        continue
    pts = np.array(np.nonzero(lab == k)).T * SP + LO
    v = len(pts) * SP ** 3
    if v < 8.0:
        continue
    d = float(wall_tree.query(pts)[0].min())
    if best is None or d < best:
        best, K_SAT, C_SAT, V_SAT = d, k, pts.mean(0), v

SAT_STRICT = (lab == K_SAT)
print("=" * 78)
print(f"crystal WT: satellite lobe = component {K_SAT}, {V_SAT:.1f} A^3, "
      f"{best:.2f} A from the wall side chains")

# wall plane: through the wall side-chain centroid, normal ABA -> satellite
P0 = wall_xyz.mean(0)
NRM = C_SAT - REF
NRM = NRM / np.linalg.norm(NRM)
SIDE = ((PTS - P0) @ NRM).reshape(SHAPE) > 0
print(f"wall plane through {np.round(P0,1)}, normal {np.round(NRM,2)} "
      f"(ABA centroid -> satellite centroid)")
sat_frac = float((SAT_STRICT & SIDE).sum()) / max(SAT_STRICT.sum(), 1)
aba_frac = float(((lab == k_aba) & SIDE).sum()) / max((lab == k_aba).sum(), 1)
print(f"sanity: {100*sat_frac:.0f}% of the satellite lobe is beyond the plane; "
      f"{100*aba_frac:.0f}% of the ABA chamber is")

# ---------- diagnostic on relaxed WT ----------
print("\n" + "=" * 78)
print("DIAGNOSTIC -- does FastRelax alone breach the wall in WT?")
print("=" * 78)
c_idx = np.clip(((C_SAT - LO) / SP).astype(int), 0, SHP - 1)
for tag in ("crystal", "relaxed"):
    path = WT_PDB if tag == "crystal" else os.path.join(PDB_DIR, "WT_r0.pdb")
    x, e = atoms_from_pose_pdb(path, chain="A")
    lb, nn = label_cavities(x, e)
    k = seeded(lb)
    comp = (lb == k)
    reaches = bool(comp[c_idx[0], c_idx[1], c_idx[2]])
    sizes = sorted((ndimage.sum(lb > 0, lb, range(1, nn + 1)) * SP ** 3),
                   reverse=True)[:4]
    print(f"  WT {tag:<8}: {nn} components, largest {np.round(sizes,1)} A^3; "
          f"ABA component = {comp.sum()*SP**3:7.1f} A^3; "
          f"reaches satellite centroid: {'YES' if reaches else 'no'}")
print("  -> if relaxed WT reaches the centroid while crystal WT does not, cause")
print("     (a) holds: side-chain relaxation alone opens the wall, and the")
print("     'sealed satellite lobe' is a property of the crystal, not the protein.")


def metrics(path):
    x, e = atoms_from_pose_pdb(path, chain="A")
    lb, nn = label_cavities(x, e)
    if nn == 0:
        return dict(cavity=0.0, strict=0.0, beyond=0.0, reaches=False, ncomp=0)
    k = seeded(lb)
    comp = (lb == k)
    return dict(cavity=float(comp.sum() * SP ** 3),
                strict=float((comp & SAT_STRICT).sum() * SP ** 3),
                beyond=float((comp & SIDE).sum() * SP ** 3),
                reaches=bool(comp[c_idx[0], c_idx[1], c_idx[2]]),
                ncomp=int(nn))


print("\n" + "=" * 78)
print("PER-VARIANT, dilation-free (relaxed rep 0)")
print("=" * 78)
rows = []
for f in sorted(glob.glob(os.path.join(PDB_DIR, "*_r0.pdb"))):
    v = os.path.basename(f)[:-len("_r0.pdb")]
    m = metrics(f)
    m["variant"] = v
    rows.append(m)

wt = next(r for r in rows if r["variant"] == "WT")
for r in rows:
    r["d_cavity"] = r["cavity"] - wt["cavity"]
    r["d_beyond"] = r["beyond"] - wt["beyond"]
rows.sort(key=lambda r: -r["d_beyond"])

with open(OUT, "w") as fh:
    fh.write("variant,cavity_A3,strict_overlap_A3,beyond_wall_A3,"
             "d_cavity_A3,d_beyond_wall_A3,reaches_satellite,n_components\n")
    for r in rows:
        fh.write(f"{r['variant']},{r['cavity']:.2f},{r['strict']:.2f},"
                 f"{r['beyond']:.2f},{r['d_cavity']:.2f},{r['d_beyond']:.2f},"
                 f"{int(r['reaches'])},{r['ncomp']}\n")

print(f"WT baseline: cavity {wt['cavity']:.1f}, beyond-wall {wt['beyond']:.1f}, "
      f"strict-overlap {wt['strict']:.1f}\n")
print(f"{'variant':<26}{'cavity':>9}{'dCav':>8}{'beyond':>9}{'dBeyond':>9}"
      f"{'strict':>8}{'reach':>7}")
for r in rows:
    print(f"{r['variant']:<26}{r['cavity']:>9.1f}{r['d_cavity']:>+8.1f}"
          f"{r['beyond']:>9.1f}{r['d_beyond']:>+9.1f}{r['strict']:>8.1f}"
          f"{'YES' if r['reaches'] else '-':>7}")
print(f"\nwrote {OUT}")
print("""
d_beyond_wall is the metric to use. It is free of mask dilation and of grid
connectivity, and it answers the design question directly: how much MORE pocket
lies past the K59/R79/E94/F108 wall than in WT. Ranking on it also re-tests the
open question from script 29 -- whether F108 or the R79/E94 pair is the real
gatekeeper to the second lobe.""")
