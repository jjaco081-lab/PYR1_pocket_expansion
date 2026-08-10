#!/usr/bin/env python
"""
24_state_comparison.py -- how far does PYR1 actually move between states, and is
Arm 2's free set big enough to see it?

Two questions, both raised as design objections to Arm 2.

Q1 (calibration). Arm 2 releases only the gate (85-89) and latch (115-117) from
coordinate constraints, and the v2 smoke test produced gate backbone RMSD of
0.358 A. Is that a real response or a still-hobbled one? The honest yardstick is
the size of the conformational difference between two EXPERIMENTAL PYR1 states:

    3K3K  PYR1 + ABA homodimer      (chains A and B, 183 res each)
    3QN1  PYR1 + ABA + HAB1 ternary (chain A closed, HAB1 chain B)

If the real gate/latch displacement between these is much larger than what the
relax protocol produces, the free set is too small and the loops are being held
by their flanking anchors -- they can bulge but not hinge.

Q2 (interface). The PYR1 surface that contacts HAB1 looks subtly different
between the dimer and the ternary complex. This tests whether that is a real,
localised backbone change or scattered coordinate noise, by measuring deviation
of the HAB1-contacting residues specifically against the deviation of the
protein as a whole.

Method: superpose on CORE Calpha only -- every protein residue EXCEPT the gate,
the latch and their +-3 flanks -- so the loops and the interface are measured as
displacements, not absorbed into the fit. Reporting Calpha and backbone
(N, CA, C, O) separately.

Run with the esmfold2 env python.
"""
import os, sys, warnings, collections
warnings.filterwarnings("ignore")
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

GATE = list(range(85, 90))      # S85-A89
LATCH = list(range(115, 118))   # H115-L117
FLANK = 3
GATE_EXT = list(range(GATE[0] - FLANK, GATE[-1] + FLANK + 1))
LATCH_EXT = list(range(LATCH[0] - FLANK, LATCH[-1] + FLANK + 1))
BB = ("N", "CA", "C", "O")

p = MMCIFParser(QUIET=True)


def load(fn, chain):
    m = p.get_structure("x", os.path.join(DATA, fn))[0]
    return {r.id[1]: r for r in m[chain] if is_aa(r) and r.id[0] == " "}, m


def kabsch(P, Q):
    """rotation+translation taking P onto Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def atoms(res, names):
    return {n: res[n].coord for n in names if n in res}


ternary, m_tern = load("3QN1.cif", "A")     # closed, HAB1-bound
dimerA, _ = load("3K3K.cif", "A")           # ABA-bound homodimer
dimerB, _ = load("3K3K.cif", "B")

print("=" * 78)
print("PYR1 CONFORMATIONAL STATES -- 3K3K (ABA dimer) vs 3QN1 (ternary, closed)")
print("=" * 78)

# HAB1-contacting PYR1 residues, measured in the ternary structure
hab1 = [r for r in m_tern["B"] if is_aa(r) and r.id[0] == " "]
hab_xyz = np.array([a.coord for r in hab1 for a in r if a.element != "H"])
iface = []
for rn, r in ternary.items():
    d = min(np.linalg.norm(hab_xyz - a.coord, axis=1).min()
            for a in r if a.element != "H")
    if d <= 5.0:
        iface.append(rn)
print(f"\nPYR1 residues within 5.0 A of HAB1 in 3QN1: {len(iface)}")
print("   " + ", ".join(f"{ternary[rn].get_resname()}{rn}" for rn in sorted(iface)))

# 3K3K chain A carries NO ligand (open/apo protomer); chain B holds A8S. So the
# A-vs-ternary comparison is the full open->closed switch stroke, while
# B-vs-ternary isolates what HAB1 binding adds ON TOP of ligand-induced closure.
iface_loop = [rn for rn in iface if rn in GATE or rn in LATCH]
iface_only = [rn for rn in iface if rn not in GATE and rn not in LATCH]
print(f"   of these, {len(iface_loop)} are gate/latch "
      f"({', '.join(str(r) for r in sorted(iface_loop))}) and "
      f"{len(iface_only)} are non-loop interface")

for label, ref in (("3K3K chain A  [APO, open protomer]", dimerA),
                   ("3K3K chain B  [ABA-bound protomer]", dimerB)):
    common = sorted(set(ref) & set(ternary))
    core = [rn for rn in common if rn not in GATE_EXT and rn not in LATCH_EXT]
    P = np.array([ref[rn]["CA"].coord for rn in core if "CA" in ref[rn] and "CA" in ternary[rn]])
    Q = np.array([ternary[rn]["CA"].coord for rn in core if "CA" in ref[rn] and "CA" in ternary[rn]])
    R, t = kabsch(P, Q)

    def dev(rns, names=("CA",)):
        d = []
        for rn in rns:
            if rn not in ref or rn not in ternary:
                continue
            a, b = atoms(ref[rn], names), atoms(ternary[rn], names)
            for n in a:
                if n in b:
                    d.append(np.linalg.norm((R @ a[n] + t) - b[n]))
        return np.array(d)

    print(f"\n--- {label} superposed on 3QN1 chain A over {len(P)} core CA ---")
    print(f"{'region':<26}{'n atoms':>9}{'RMSD':>9}{'max dev':>10}")
    for nm, rns, at in (
            ("core (fitted)", core, ("CA",)),
            ("ALL residues", common, ("CA",)),
            ("gate 85-89", GATE, ("CA",)),
            ("gate 85-89 backbone", GATE, BB),
            ("gate +-3 flanks", GATE_EXT, ("CA",)),
            ("latch 115-117", LATCH, ("CA",)),
            ("latch 115-117 backbone", LATCH, BB),
            ("latch +-3 flanks", LATCH_EXT, ("CA",)),
            ("HAB1 interface (all)", iface, ("CA",)),
            ("HAB1 iface NON-loop", iface_only, ("CA",)),
            ("HAB1 iface NON-loop bb", iface_only, BB)):
        d = dev(rns, at)
        if len(d):
            print(f"{nm:<26}{len(d):>9}{np.sqrt((d**2).mean()):>9.2f}{d.max():>10.2f}")

    # biggest movers overall
    per = []
    for rn in common:
        d = dev([rn], ("CA",))
        if len(d):
            per.append((float(d[0]), rn))
    per.sort(reverse=True)
    print("   largest CA displacements: " + ", ".join(
        f"{ternary[rn].get_resname()}{rn} {v:.1f}A" for v, rn in per[:8]))

print("""
=============================================================================
HOW TO READ THIS
=============================================================================
The "gate backbone" and "latch backbone" rows are the magnitude Arm 2 must be
able to reproduce. Compare them against the v2 smoke test (WT, gate_bb_rmsd
0.358 A, latch_bb_rmsd 0.636 A) and against v1 (0.081 / 0.086, constrained).

If the experimental difference is several Angstrom while the relax protocol
moves the loops by a fraction of one, then releasing residues 85-89 and 115-117
alone is NOT enough: their flanking residues are still pinned by coordinate
constraints, so the loops can only bulge locally, not hinge as rigid bodies.
The fix is to extend the free set to the flanks (the GATE_EXT / LATCH_EXT sets
used here) rather than to sample harder.

The "HAB1 interface" rows separate a genuine binding-induced surface change from
crystallographic noise: if interface deviation is at the level of the ALL-residue
row, it is noise; if it stands above it, the surface really does differ between
the dimer and the ternary complex.""")
