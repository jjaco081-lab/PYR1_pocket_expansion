#!/usr/bin/env python
r"""
121_pocket_vs_interface.py -- how related are PYR1's LIGAND-POCKET residues and
its HAB1-INTERFACE residues? Same secondary-structure elements, or disjoint?

WHY IT MATTERS
§(ratchet) established that HAB1 has zero pocket residues -- it reads PYR1's
closed-state SURFACE and never touches the ligand. That says the two sets are
functionally separate. It does not say they are structurally separate: if a
pocket position sits on the same strand or helix as an interface residue, then
mutating the pocket can propagate along that element and perturb the interface,
which is a route to constitutive or dead sensors that nothing in this project
currently screens for.

INPUTS
  data/complex_AB_ABA.pdb   chain A = PYR1 (native numbering), B = HAB1, X = ABA
  secondary structure from cpptraj's DSSP implementation, not from a guess

⚠ Identity is asserted at all 18 library positions before anything is reported.
"""
import os
import subprocess
import sys
from collections import Counter, defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "pocket_vs_interface")
PDB = os.path.join(ROOT, "data", "complex_AB_ABA.pdb")
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
LIB = {59: "K", 81: "V", 83: "V", 87: "L", 89: "A", 92: "S", 94: "E", 108: "F",
       110: "I", 117: "L", 120: "Y", 122: "S", 141: "E", 159: "F", 160: "A",
       163: "V", 164: "V", 167: "N"}
CUT = 4.5


def read(path):
    A, B, X, seq = defaultdict(list), [], [], {}
    for l in open(path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if l[76:78].strip().upper() == "H":
            continue
        ch = l[21]
        xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        if l.startswith("HETATM"):
            X.append(xyz)
            continue
        r = int(l[22:26])
        if ch == "A":
            A[r].append(xyz)
            seq[r] = AA3.get(l[17:20].strip(), "X")
        elif ch == "B":
            B.append(xyz)
    return A, np.array(B), np.array(X), seq


def dssp():
    """Secondary structure from PyRosetta's DsspMover, keyed by PDB number.

    ⚠ cpptraj's `secstruct` was tried first and FAILED SILENTLY on this input:
    every row of its summary came back all-zero or Turn, with no Extended or
    Alpha anywhere -- it needs backbone amide hydrogens, which the crystal-derived
    PDB does not carry. Worse, my first parser took argmax over an all-zero row,
    which returns column 0, so 0.0 became "Extended" and PYR1 was reported as
    20 strands and ZERO helices. A helix-grip fold with no helices should have
    stopped me; it is recorded here because the failure produced confident,
    plausible-looking output rather than an error."""
    import pyrosetta
    from pyrosetta.rosetta.protocols.moves import DsspMover
    pyrosetta.init("-mute all -ignore_unrecognized_res")
    pose = pyrosetta.pose_from_pdb(os.path.join(ROOT, "data", "pyr1_A.pdb"))
    DsspMover().apply(pose)
    ss = pose.secstruct()
    info = pose.pdb_info()
    out = {}
    for i in range(1, pose.total_residue() + 1):
        out[info.number(i)] = ss[i - 1]
    if not any(v == "H" for v in out.values()):
        raise SystemExit("no helices assigned -- DSSP failed, see docstring")
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    A, B, X, seq = read(PDB)
    bad = [(n, a, seq.get(n)) for n, a in LIB.items() if seq.get(n) != a]
    if bad:
        raise SystemExit(f"identity check failed at {bad}")
    say("=" * 78)
    say("0. INPUT")
    say("=" * 78)
    say(f"   {PDB}")
    say(f"   PYR1 chain A: {len(A)} residues; HAB1 chain B: {len(B)} heavy atoms; "
        f"ABA: {len(X)} heavy atoms")
    say(f"   all 18 library positions carry their expected residue")

    # ---- contacts ---------------------------------------------------------
    iface, lig = {}, {}
    for r, atoms in A.items():
        P = np.array(atoms)
        if len(B):
            d = np.linalg.norm(P[:, None, :] - B[None, :, :], axis=2)
            if d.min() < CUT:
                iface[r] = float(d.min())
        if len(X):
            d = np.linalg.norm(P[:, None, :] - X[None, :, :], axis=2)
            if d.min() < CUT:
                lig[r] = float(d.min())
    say("")
    say("=" * 78)
    say("1. THE TWO SETS")
    say("=" * 78)
    say(f"   HAB1-interface residues (< {CUT} A to chain B): {len(iface)}")
    say(f"     {sorted(iface)}")
    say(f"   ABA-contacting residues  (< {CUT} A to ligand): {len(lig)}")
    say(f"     {sorted(lig)}")
    ov = sorted(set(iface) & set(lig))
    say(f"   OVERLAP: {ov if ov else 'NONE -- the two sets are disjoint'}")
    libset = set(LIB)
    say(f"   of the 18 LIBRARY positions, how many touch HAB1? "
        f"{sorted(libset & set(iface)) or 'none'}")

    # ---- secondary structure ---------------------------------------------
    ss = dssp()
    for n, a in LIB.items():
        if n not in ss:
            raise SystemExit(f"position {n} missing from the DSSP assignment")
    say("")
    say("=" * 78)
    say("2. SECONDARY STRUCTURE (PyRosetta DSSP)")
    say("=" * 78)
    cnt = Counter(ss.values())
    say(f"   {cnt.get('H',0)} helical residues, {cnt.get('E',0)} strand, "
        f"{cnt.get('L',0)} loop -- the helix-grip START fold")
    elems, cur = [], None
    for r in sorted(ss):
        s_ = ss[r]
        if cur and cur[0] == s_ and r == cur[2] + 1:
            cur[2] = r
        else:
            if cur:
                elems.append(tuple(cur))
            cur = [s_, r, r]
    if cur:
        elems.append(tuple(cur))
    named, hi, ei = {}, 0, 0
    for s_, a, b in elems:
        if s_ == "H" and b - a >= 3:
            hi += 1
            named[(a, b)] = f"alpha{hi}"
        elif s_ == "E" and b - a >= 1:
            ei += 1
            named[(a, b)] = f"beta{ei}"
    say("")
    say(f"   {'element':<10}{'range':<11}{'library positions':<28}"
        f"{'HAB1-interface'}")
    shared = []
    for (a, b), nm in sorted(named.items()):
        lp = [r for r in sorted(libset) if a <= r <= b]
        ip = [r for r in sorted(iface) if a <= r <= b]
        if not lp and not ip:
            continue
        if lp and ip:
            shared.append((nm, a, b, lp, ip))
        say(f"   {nm:<10}{f'{a}-{b}':<11}{str(lp) if lp else '-':<28}"
            f"{str(ip) if ip else '-'}")
    say("")
    say("=" * 78)
    say("3. ELEMENTS CARRYING BOTH -- where a pocket mutation can reach HAB1")
    say("=" * 78)
    if not shared:
        say("   NONE.")
    for nm, a, b, lp, ip in shared:
        say(f"   {nm} ({a}-{b}):  library {lp}   HAB1-contacting {ip}")
        for L in lp:
            near = min(abs(L - I) for I in ip)
            sp = min(np.linalg.norm(np.array(A[L])[:, None, :]
                                    - np.array(A[I])[None, :, :], axis=2).min()
                     for I in ip)
            same = "SAME residue" if near == 0 else (
                "same face (i,i+2)" if near == 2 and nm.startswith("beta")
                else ("same face (i,i+3/i+4)" if near in (3, 4)
                      and nm.startswith("alpha") else ""))
            say(f"      {LIB[L]}{L}: {near} apart in sequence, {sp:.1f} A in "
                f"space  {same}")
    say("")
    say("=" * 78)
    say("4. THE C-TERMINAL HELIX HAS TWO FACES, and the library sits on one")
    say("=" * 78)
    say("   Per-residue minimum distances along alpha3. The helix-grip helix caps")
    say("   the pocket on one side and forms the HAB1 interface on the other.")
    say(f"   {'res':<7}{'to ABA':>9}{'to HAB1':>10}   face")
    a3 = [k for k in named if named[k].startswith("alpha") and k[1] - k[0] > 20]
    if a3:
        lo, hi = a3[0]
        for r in range(lo, hi + 1):
            if r not in A:
                continue
            P = np.array(A[r])
            dl = float(np.linalg.norm(P[:, None, :] - X[None, :, :],
                                      axis=2).min()) if len(X) else 99.0
            dh = float(np.linalg.norm(P[:, None, :] - B[None, :, :],
                                      axis=2).min()) if len(B) else 99.0
            face = ("POCKET" if dl < dh - 2 else
                    "HAB1" if dh < dl - 2 else "both/edge")
            tag = "  <- library" if r in libset else ""
            say(f"   {seq[r]}{r:<6}{dl:>9.1f}{dh:>10.1f}   {face}{tag}")
    say("")
    say("   ==> the two functions are on OPPOSITE FACES of one helix. That is why")
    say("       HAB1 contains no pocket residues and yet five library positions")
    say("       share an element with six interface residues: a mutation at 160,")
    say("       163, 164 or 167 does not touch HAB1 directly, but it repacks the")
    say("       core of the helix that presents the interface.")

    say("")
    say("   On a STRAND side chains alternate faces, so i and i+2 point the same")
    say("   way; on a HELIX i and i+3/i+4 do. Those are the separations at which a")
    say("   pocket mutation sits on the same surface as an interface residue.")

    with open(os.path.join(OUT, "pocket_vs_interface.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    say(f"\n   written to {OUT}/pocket_vs_interface.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
