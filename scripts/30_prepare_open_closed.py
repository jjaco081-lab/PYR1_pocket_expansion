#!/usr/bin/env python
"""
30_prepare_open_closed.py -- build the matched OPEN and CLOSED PYR1 monomers that
Arm 2b (script 31) compares.

Why this script exists
----------------------
Arm 2b asks whether a variant still prefers the CLOSED state, by scoring the same
sequence on two different backbone conformers and taking the difference. That
subtraction is only meaningful if the two poses contain **exactly the same
residues**: Rosetta's total score is extensive, so a single extra residue in one
conformer shifts the difference by far more than the effect being measured.

  closed : 3QN1 chain A -- PYR1 in the HAB1-bound ternary complex (gate/latch shut)
  open   : 3K3K chain A -- the APO protomer of the ABA-bound homodimer

3K3K chain B carries A8S and is already closed-like (script 24: gate backbone
0.93 A from the ternary state). Chain A is ligand-free and is the genuinely open
conformer -- gate backbone 5.23 A, latch 5.34 A from closed. That 5 A stroke is
the conformational change the sensor depends on.

Both are written APO (no ABA, no HAB1, no waters, no metals): Arm 2b measures the
protein's INTRINSIC conformational preference, which is what decides whether a
variant closes without ligand (a constitutive, ligand-independent binder -- the
expensive false positive) or fails to close at all (a silent dead receptor).

What this does
--------------
  1. reads both chains, keeping only standard amino acids
  2. takes the INTERSECTION of residue numbers present in both
  3. drops any residue lacking a complete backbone (N, CA, C, O) in either
  4. keeps altloc A only, strips hydrogens and OXT
  5. verifies the two files have identical residue numbering AND identical
     residue identities, then reports the gate/latch displacement between them
     as a sanity check against script 24

Outputs data/pyr1_closed_A.pdb and data/pyr1_open_A.pdb.

Run with the esmfold2 env python.
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
BB = ("N", "CA", "C", "O")
GATE, LATCH = list(range(85, 90)), list(range(115, 118))

p = MMCIFParser(QUIET=True)


def load(cif, chain):
    m = p.get_structure("x", os.path.join(DATA, cif))[0]
    out = {}
    for r in m[chain]:
        if r.id[0] != " " or not is_aa(r, standard=True):
            continue
        if not all(a in r for a in BB):
            continue
        out[r.id[1]] = r
    return out


closed = load("3QN1.cif", "A")
openp = load("3K3K.cif", "A")
print(f"3QN1 chain A (closed): {len(closed)} residues with complete backbone")
print(f"3K3K chain A (open)  : {len(openp)} residues with complete backbone")

common = sorted(set(closed) & set(openp))
mismatch = [rn for rn in common
            if closed[rn].get_resname() != openp[rn].get_resname()]
if mismatch:
    print(f"  !! {len(mismatch)} residue-identity mismatches, dropping: {mismatch[:10]}")
    common = [rn for rn in common if rn not in mismatch]

print(f"common residue set: {len(common)}  ({common[0]}-{common[-1]})")
missing_c = sorted(set(closed) - set(common))
missing_o = sorted(set(openp) - set(common))
print(f"  dropped from closed: {len(missing_c)} {missing_c[:12]}")
print(f"  dropped from open  : {len(missing_o)} {missing_o[:12]}")

gaps = [(a, b) for a, b in zip(common, common[1:]) if b - a > 1]
if gaps:
    print(f"  NOTE chain breaks in the common set: {gaps}")
    print("       (FastRelax handles these; they are identical in both poses so"
          " they cancel in the difference)")


class Keep(Select):
    def __init__(self, keep):
        self.keep = set(keep)

    def accept_residue(self, r):
        return r.id[0] == " " and r.id[1] in self.keep

    def accept_atom(self, a):
        if a.element == "H":
            return False
        if a.get_id() == "OXT":
            return False
        alt = a.get_altloc()
        return alt in (" ", "A")


def write(cif, chain, path):
    m = p.get_structure("x", os.path.join(DATA, cif))[0]
    for ch in list(m):
        if ch.id != chain:
            m.detach_child(ch.id)
    io = PDBIO()
    io.set_structure(m)
    io.save(path, Keep(common))
    return path


f_closed = write("3QN1.cif", "A", os.path.join(DATA, "pyr1_closed_A.pdb"))
f_open = write("3K3K.cif", "A", os.path.join(DATA, "pyr1_open_A.pdb"))


def superpose_file(path, ref_path, exclude):
    """Rigid-body place `path` into `ref_path`'s frame using core CA only.

    Rosetta scores are frame-independent, so this changes no energy. It is done
    so that the open and closed poses share a coordinate frame, which lets the
    same ABA-derived cavity seed and the same reference geometry be used for
    both states downstream.
    """
    def rd(pth):
        d = {}
        for ln in open(pth):
            if ln.startswith("ATOM"):
                d.setdefault(int(ln[22:26]), {})[ln[12:16].strip()] = np.array(
                    [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        return d
    mob, ref = rd(path), rd(ref_path)
    core = [r for r in sorted(set(mob) & set(ref)) if r not in exclude
            and "CA" in mob[r] and "CA" in ref[r]]
    P = np.array([mob[r]["CA"] for r in core])
    Q = np.array([ref[r]["CA"] for r in core])
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    dsign = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, dsign]) @ U.T
    out = []
    for ln in open(path):
        if ln.startswith("ATOM"):
            v = R @ (np.array([float(ln[30:38]), float(ln[38:46]),
                               float(ln[46:54])]) - pc) + qc
            ln = f"{ln[:30]}{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}{ln[54:]}"
        out.append(ln)
    open(path, "w").writelines(out)
    return len(core)


EXCL = set(range(82, 93)) | set(range(112, 121))
n_fit = superpose_file(f_open, f_closed, EXCL)
print(f"open pose superposed onto the closed frame over {n_fit} core CA "
      "(rigid body; does not change any Rosetta score)")


def readpdb(path):
    res = {}
    for line in open(path):
        if line.startswith("ATOM"):
            rn = int(line[22:26])
            res.setdefault(rn, {})[line[12:16].strip()] = (
                line[17:20].strip(),
                np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    return res


rc, ro = readpdb(f_closed), readpdb(f_open)
assert set(rc) == set(ro), "residue sets differ after writing"
for rn in rc:
    assert rc[rn]["CA"][0] == ro[rn]["CA"][0], f"identity mismatch at {rn}"
print(f"\nwrote {f_closed}  ({len(rc)} residues)")
print(f"wrote {f_open}    ({len(ro)} residues)")
print("residue numbering and identities are IDENTICAL in both -- scores comparable")

# sanity check against script 24: superpose on core, measure gate/latch stroke
core = [rn for rn in rc if rn not in range(82, 93) and rn not in range(112, 121)]
P = np.array([ro[rn]["CA"][1] for rn in core])
Q = np.array([rc[rn]["CA"][1] for rn in core])
pc, qc = P.mean(0), Q.mean(0)
H = (P - pc).T @ (Q - qc)
U, S, Vt = np.linalg.svd(H)
d = np.sign(np.linalg.det(Vt.T @ U.T))
R = Vt.T @ np.diag([1, 1, d]) @ U.T


def rms(rns, names=BB):
    v = []
    for rn in rns:
        for n in names:
            if n in rc[rn] and n in ro[rn]:
                v.append(np.linalg.norm((R @ (ro[rn][n][1] - pc) + qc) - rc[rn][n][1]))
    return float(np.sqrt(np.mean(np.square(v))))


print(f"\nopen -> closed stroke (superposed on {len(core)} core CA):")
print(f"   core CA          {rms(core, ('CA',)):.2f} A")
print(f"   gate 85-89 bb    {rms([r for r in GATE if r in rc]):.2f} A")
print(f"   latch 115-117 bb {rms([r for r in LATCH if r in rc]):.2f} A")
print("   (script 24 reported 5.23 / 5.34 A on the untrimmed structures)")
