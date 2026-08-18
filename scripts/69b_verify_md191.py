#!/usr/bin/env python
"""Measurements behind 69b_verify_md191.sh -- run from that script, not directly."""
import os
import sys

import numpy as np
from Bio.PDB import PDBParser

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V = os.path.join(ROOT, "data", "md191", "_verify")
S = os.path.join(ROOT, "data", "structures_191")

GATE, LATCH = list(range(85, 90)), list(range(115, 118))
LB7A5, TAIL = list(range(148, 157)), list(range(182, 192))
BB = ("N", "CA", "C", "O")
CLAIMS = {"S1_apo_open": ("open", False), "S2_holo_closed": ("closed", True),
          "S9_apo_closed": ("closed", False), "S10_holo_open": ("open", True)}
p = PDBParser(QUIET=True)


def bb(model, nats):
    g = {}
    for ch in model:
        for r in ch:
            if r.id[0] == " " and r.id[1] in nats:
                for a in r:
                    if a.get_id() in BB:
                        g[(r.id[1], a.get_id())] = a.coord
    return g


def loop_rmsd(mob, ref, core, loop):
    X, Y = bb(mob, core), bb(ref, core)
    k = sorted(set(X) & set(Y))
    A, B = np.array([X[i] for i in k]), np.array([Y[i] for i in k])
    ma, mb = A.mean(0), B.mean(0)
    U, _, Vt = np.linalg.svd((A - ma).T @ (B - mb))
    R = U @ np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]) @ Vt
    Xl, Yl = bb(mob, loop), bb(ref, loop)
    kl = sorted(set(Xl) & set(Yl))
    a = (np.array([Xl[i] for i in kl]) - ma) @ R
    b = np.array([Yl[i] for i in kl]) - mb
    return float(np.sqrt(((a - b) ** 2).sum() / len(kl))), len(k)


ref_o = p.get_structure("o", os.path.join(S, "pyr1_open_191.pdb"))[0]
ref_c = p.get_structure("c", os.path.join(S, "pyr1_closed_191.pdb"))[0]
mobile = set(GATE) | set(LATCH) | set(LB7A5) | set(TAIL)
core = [n for n in range(1, 192) if n not in mobile]

fail = False
for name, (want_conf, want_lig) in CLAIMS.items():
    m = p.get_structure("m", os.path.join(V, f"{name}.frame1.pdb"))[0]
    # cpptraj writes the ligand as ATOM, not HETATM, so a hetflag test does NOT
    # separate them -- select by residue name or A8S is counted as a 192nd protein
    # residue and then measured against itself (closest contact 0.00 A).
    lig = [r for ch in m for r in ch if r.get_resname().strip() == "A8S"]
    prot = [r for ch in m for r in ch
            if r.id[0] == " " and r.get_resname().strip() != "A8S"]
    nums = sorted(r.id[1] for r in prot)
    print(f"=== {name}   claims: {want_conf}, ligand={'yes' if want_lig else 'no'} ===")
    print(f"  protein residues {len(prot)} ({nums[0]}..{nums[-1]}) | A8S copies {len(lig)}")

    go, ncore = loop_rmsd(m, ref_o, core, GATE)
    gc, _ = loop_rmsd(m, ref_c, core, GATE)
    lo, _ = loop_rmsd(m, ref_o, core, LATCH)
    lc, _ = loop_rmsd(m, ref_c, core, LATCH)
    got_conf = "open" if go < gc else "closed"
    print(f"  gate  to open {go:5.2f} A | to closed {gc:5.2f} A   -> {got_conf}")
    print(f"  latch to open {lo:5.2f} A | to closed {lc:5.2f} A   "
          f"(core fit {ncore} atoms)")

    if len(prot) != 191 or nums != list(range(1, 192)):
        print("  !! protein is not 191 residues numbered 1..191"); fail = True
    if got_conf != want_conf:
        print(f"  !! CONFORMATION MISMATCH: name says {want_conf}, structure says {got_conf}")
        fail = True
    if bool(lig) != want_lig:
        print(f"  !! LIGAND MISMATCH: name says {want_lig}, topology has {len(lig)}")
        fail = True

    # A ligand built in the WRONG FRAME passes every check above and sits in bulk
    # solvent. Only a real distance to a known pocket residue catches that.
    if lig:
        L = np.array([a.coord for a in lig[0] if a.element != "H"])
        P = np.array([a.coord for r in prot for a in r if a.element != "H"])
        allmin = float(np.linalg.norm(L[:, None] - P[None], axis=2).min())
        k59 = [r for r in prot if r.id[1] == 59][0]
        assert k59.get_resname() == "LYS", f"residue 59 is {k59.get_resname()}"
        nz = k59["NZ"].coord
        d59 = float(np.linalg.norm(L - nz, axis=1).min())
        shell = sorted({r.id[1] for r in prot
                        if min(np.linalg.norm(L - a.coord, axis=1).min()
                               for a in r if a.element != "H") < 4.5})
        print(f"  ABA: closest protein contact {allmin:.2f} A | K59 NZ to nearest "
              f"ABA heavy atom {d59:.2f} A | {len(shell)} lining residues")
        if d59 > 5.0:
            print(f"  !! ABA IS NOT IN THE POCKET (K59 NZ {d59:.2f} A) -- wrong frame?")
            fail = True
        if allmin < 1.8:
            print(f"  !! hard clash at {allmin:.2f} A"); fail = True
    print()

if fail:
    sys.exit("VERIFICATION FAILED -- do not submit")
print("ALL FOUR CELLS VERIFIED: each system is the conformation and occupancy its "
      "name claims, and both ligands are in the pocket.")
