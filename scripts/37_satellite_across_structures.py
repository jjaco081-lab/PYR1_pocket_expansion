#!/usr/bin/env python
"""
37_satellite_across_structures.py -- is the satellite lobe sealed in EVERY PYR1
crystal structure, or only in 3QN1?

Why
---
Script 34 found that FastRelax alone merges the satellite lobe with the ABA
chamber in WT PYR1: the crystal (3QN1 chain A) has them as two separate voids,
the relaxed structure has one. That was read as "the seal is a crystal artefact
-- the wall is made of side-chain rotamers, not backbone." But a single relax
protocol could equally be producing a Rosetta rotamer artefact.

The decisive check needs no simulation at all: look at every independently
determined PYR1 chain available and ask whether the satellite lobe is connected
in any of them. Deposited coordinates from different crystals, space groups and
ligand states are independent samples of the rotamer arrangement.

  3QN1 A  PYR1 + ABA + HAB1        closed, ternary
  3K3K A  PYR1 homodimer           APO protomer (open)
  3K3K B  PYR1 homodimer           ABA-bound protomer
  3K90 A  PYR1                     ABA-bound
  3K90 B  PYR1                     apo (ACY/GOL only)
  3K90 C  PYR1                     apo (GOL only)
  3K90 D  PYR1                     ABA-bound

If the lobe is connected in some deposited chains, the seal is a property of
particular crystal snapshots and the script-34 conclusion holds. If 3QN1 is the
only sealed one, or all are sealed, the interpretation changes accordingly.

Every chain is superposed onto 3QN1 chain A over core CA (excluding gate/latch
and flanks) so the fixed satellite centroid and wall plane apply to all.

Run with the esmfold2 env python.
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_cavity import VDW, DIRS, atoms_from_pose_pdb  # noqa: E402

DATA = os.path.join(ROOT, "data")
SP, PROBE, BUR, BOX = 0.5, 1.4, 0.88, 16.0
WALL = (59, 79, 94, 108)
EXCL = set(range(82, 93)) | set(range(112, 121))

CHAINS = [("3QN1", "A", "ABA + HAB1 (closed, ternary)"),
          ("3K3K", "A", "homodimer, APO protomer (open)"),
          ("3K3K", "B", "homodimer, ABA-bound protomer"),
          ("3K90", "A", "ABA-bound"),
          ("3K90", "B", "apo (ACY/GOL)"),
          ("3K90", "C", "apo (GOL)"),
          ("3K90", "D", "ABA-bound")]

lig_xyz, _ = atoms_from_pose_pdb(os.path.join(DATA, "aba_xtal.pdb"),
                                 exclude_resnames=())
REF = lig_xyz.mean(0)
LO = REF - BOX
AXES = [np.arange(LO[i], REF[i] + BOX, SP) for i in range(3)]
GRID = np.stack(np.meshgrid(*AXES, indexing="ij"), -1)
SHAPE = GRID.shape[:3]
PTS = GRID.reshape(-1, 3)
SHP = np.array(SHAPE)
STEPS = np.arange(1.0, 15.0, 0.75)

p = MMCIFParser(QUIET=True)


def chain_atoms(pdb, ch):
    m = p.get_structure("x", os.path.join(DATA, f"{pdb}.cif"))[0]
    res = {}
    for r in m[ch]:
        if r.id[0] != " " or not is_aa(r, standard=True):
            continue
        d = {}
        for a in r:
            if a.element == "H" or a.get_altloc() not in (" ", "A"):
                continue
            d[a.get_id()] = np.array(a.coord, float)
        if d:
            res[r.id[1]] = (r.get_resname(), d)
    return res


REFRES = chain_atoms("3QN1", "A")


def superpose(res):
    core = [rn for rn in sorted(set(res) & set(REFRES))
            if rn not in EXCL and "CA" in res[rn][1] and "CA" in REFRES[rn][1]]
    P = np.array([res[rn][1]["CA"] for rn in core])
    Q = np.array([REFRES[rn][1]["CA"] for rn in core])
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    out = {}
    for rn, (nm, atoms) in res.items():
        out[rn] = (nm, {k: R @ (v - pc) + qc for k, v in atoms.items()})
    return out, len(core)


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


def flat(res):
    xyz, el = [], []
    for rn, (nm, atoms) in res.items():
        for k, v in atoms.items():
            xyz.append(v)
            el.append("C" if k[0] not in "NOS" else k[0])
    return np.array(xyz), el


# satellite centroid + wall plane, defined once on 3QN1 A
x0, e0 = flat(REFRES)
lab0, n0 = label_cavities(x0, e0)
ci = np.clip(((REF - LO) / SP).astype(int), 0, SHP - 1)
k_aba0 = lab0[ci[0], ci[1], ci[2]]
wall0 = np.array([v for rn in WALL if rn in REFRES
                  for k, v in REFRES[rn][1].items()
                  if k not in ("N", "CA", "C", "O")])
wtree = cKDTree(wall0)
best = None
for k in range(1, n0 + 1):
    if k == k_aba0:
        continue
    pts = np.array(np.nonzero(lab0 == k)).T * SP + LO
    if len(pts) * SP ** 3 < 8.0:
        continue
    dd = float(wtree.query(pts)[0].min())
    if best is None or dd < best[0]:
        best = (dd, pts.mean(0), len(pts) * SP ** 3)
D_SAT, C_SAT, V_SAT = best
c_idx = np.clip(((C_SAT - LO) / SP).astype(int), 0, SHP - 1)
P0 = wall0.mean(0)
NRM = (C_SAT - REF) / np.linalg.norm(C_SAT - REF)
SIDE = ((PTS - P0) @ NRM).reshape(SHAPE) > 0

print("=" * 96)
print(f"satellite lobe defined on 3QN1 A: {V_SAT:.1f} A^3, centre "
      f"{np.round(C_SAT,1)}, {D_SAT:.2f} A from the wall side chains")
print("=" * 96)
print(f"{'structure':<10}{'state':<34}{'fit':>5}{'ABA cav':>9}{'beyond':>8}"
      f"{'sat comp':>10}{'CONNECTED?':>12}")

for pdb, ch, note in CHAINS:
    try:
        res = chain_atoms(pdb, ch)
    except KeyError:
        print(f"{pdb+' '+ch:<10}{note:<34}   chain absent")
        continue
    res, nfit = superpose(res)
    xyz, el = flat(res)
    lab, n = label_cavities(xyz, el)
    if n == 0:
        print(f"{pdb+' '+ch:<10}{note:<34}{nfit:>5}   no cavity")
        continue
    k = lab[ci[0], ci[1], ci[2]]
    if k == 0:
        nz = np.array(np.nonzero(lab)).T
        k = lab[tuple(nz[np.argmin(np.linalg.norm(nz - ci, axis=1))])]
    comp = (lab == k)
    connected = bool(comp[c_idx[0], c_idx[1], c_idx[2]])
    # if not connected, is there a separate void sitting at the satellite site?
    ksat = lab[c_idx[0], c_idx[1], c_idx[2]]
    satvol = float((lab == ksat).sum() * SP ** 3) if ksat > 0 else 0.0
    print(f"{pdb+' '+ch:<10}{note:<34}{nfit:>5}{comp.sum()*SP**3:>9.1f}"
          f"{(comp & SIDE).sum()*SP**3:>8.1f}{satvol:>10.1f}"
          f"{('YES' if connected else 'no'):>12}")

print("""
  ABA cav   = volume of the cavity component containing the ABA site
  beyond    = how much of it lies past the K59/R79/E94/F108 wall plane
  sat comp  = volume of whatever void occupies the satellite centroid (0 = none)
  CONNECTED = the ABA-site cavity reaches the satellite centroid

If CONNECTED is YES in any deposited chain, the seal seen in 3QN1 is a property
of that snapshot rather than of the protein, and the script-34 reading stands.
If every deposited chain is sealed while the relaxed model is not, the merge is
more likely a Rosetta rotamer artefact and should be treated with suspicion.""")
