#!/usr/bin/env python
r"""
109_graft_pocket_landscape.py -- what do the DONOR pockets actually look like?

79 measured the graft candidates from PYR1's point of view: how many of PYR1's 19
machinery residues have a structurally equivalent position in the donor. That
answers "can the machinery be carried over", and it is silent on the question
asked here:

    how many residues line the DONOR's own pocket, and where are they?

That matters because the whole premise of the inversion (§9g) is that the donor
brings a bigger cavity. A 570 A^3 cavity lined by the same ~19 residues is a
different proposition from one lined by 30 -- the first is a longer version of
PYR1's pocket, the second is a different architecture with more positions to
design at, and more that can go wrong.

METHOD. The donor's pocket is defined from its own bound ligand where one exists
(residues with any side-chain heavy atom within 4.5 A), which is the same
criterion used for PYR1 throughout this project. Apo donors are reported as such
rather than guessed at -- a cavity-detection proxy would not be comparable to a
ligand-contact count, and mixing the two would make the comparison meaningless.

Backbone-only contacts are excluded: the question is how many SIDE CHAINS line
the pocket, since those are the positions a library would randomise.
"""
import os
import sys
from collections import defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "graft_cache")
CUT = 4.5
BACKBONE = {"N", "CA", "C", "O", "OXT"}
AA3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
SKIP_HET = {"HOH", "GOL", "EDO", "SO4", "PO4", "CL", "NA", "MG", "K", "CA",
            "MN", "ZN", "ACT", "DMS", "PEG", "TRS", "IOD", "BR", "NO3", "FMT"}

# candidate -> chain to use, from results/graft_candidates/graft_geometry.json
CANDS = [("2PCS", "A", 569.6), ("2NS9", "B", 455.3), ("2BK0", "A", 344.5),
         ("6AWV", "C", 319.0), ("3TFZ", "E", 207.1), ("3OQU", "B", 203.7),
         ("4DSB", "B", 161.8), ("3QRZ", "A", 210.0)]


def parse(path):
    """auth-chain-aware mmCIF atom records."""
    prot, het = defaultdict(list), defaultdict(list)
    for l in open(path):
        p = l.split()
        if not p or p[0] not in ("ATOM", "HETATM") or len(p) < 21:
            continue
        comp, atom, ch = p[5], p[3], p[18]
        try:
            xyz = np.array([float(p[10]), float(p[11]), float(p[12])])
        except ValueError:
            continue
        if p[0] == "ATOM" and comp in AA3:
            prot[ch].append((int(p[16]) if p[16].lstrip("-").isdigit() else 0,
                             comp, atom, xyz))
        elif p[0] == "HETATM" and comp not in SKIP_HET:
            het[(ch, comp)].append((atom, xyz))
    return prot, het


def main():
    print(f"{'PDB':<6}{'ch':>3}{'cavity':>9}{'ligand':>14}{'heavy':>7}"
          f"{'pocket residues':>17}{'sidechain-contacting':>22}")
    print("-" * 80)
    rows = []
    for pdb, chain, cav in CANDS:
        f = os.path.join(CACHE, f"{pdb}.cif")
        if not os.path.exists(f):
            print(f"{pdb:<6} missing")
            continue
        prot, het = parse(f)
        # biggest non-solvent hetero group in or near this chain
        best, bn = None, 0
        for (ch, comp), atoms in het.items():
            hv = [a for a in atoms if not a[0].startswith("H")]
            if len(hv) > bn:
                best, bn = (ch, comp, hv), len(hv)
        if best is None or bn < 8:
            print(f"{pdb:<6}{chain:>3}{cav:>9.0f}{'(apo)':>14}{'-':>7}"
                  f"{'not measurable':>17}{'-':>22}")
            rows.append((pdb, cav, None, None))
            continue
        lch, lcomp, lig = best
        L = np.array([a[1] for a in lig])
        # The ligand's auth chain need not match the chain the graft analysis
        # indexed (3OQU/4DSB/3TFZ all returned zero contacts against the recorded
        # chain). Pick the protein chain that actually contacts this ligand
        # rather than trusting the label.
        def contacts(ch):
            r = defaultdict(lambda: [False, None])
            for num, comp, atom, xyz in prot[ch]:
                if np.linalg.norm(L - xyz, axis=1).min() <= CUT:
                    r[num][1] = comp
                    if atom not in BACKBONE:
                        r[num][0] = True
            return r
        cand = {ch: contacts(ch) for ch in prot}
        chain = max(cand, key=lambda c: len(cand[c]))
        res = cand[chain]
        allres = len(res)
        sc = sum(1 for v in res.values() if v[0])
        print(f"{pdb:<6}{chain:>3}{cav:>9.0f}{lcomp:>14}{bn:>7}{allres:>17}{sc:>22}")
        rows.append((pdb, cav, allres, sc))

    print("\nPYR1 reference: cavity 174 A^3, 19-20 side-chain-lining residues "
          "(README §2), 18 of which the Tian libraries randomise.")
    good = [(c, s) for _, c, _, s in rows if s]
    if len(good) >= 3:
        c = np.array([x[0] for x in good], float)
        s = np.array([x[1] for x in good], float)
        r = np.corrcoef(c, s)[0, 1]
        print(f"\ncavity volume vs number of side-chain-lining residues: "
              f"r = {r:+.2f} over n = {len(good)} donors with a ligand")
        print("  cavity per lining residue (A^3):  " +
              ", ".join(f"{p}={cv/sv:.0f}" for (p, cv, _, sv) in rows if sv))


if __name__ == "__main__":
    main()
