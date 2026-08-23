#!/usr/bin/env python
r"""
101_xtal_rotamers.py -- build the quad copy using 4WVO's OWN sidechains.

    101_xtal_rotamers.py <wt_heavy.pdb> <4WVO.cif> <out.pdb>

WHY THIS EXISTS
The first version of the crystal arm gave copy2 tleap-built rotamers, on the
theory that holding rotamer provenance fixed across both arms would isolate the
ligand pose as the single variable. Measured, that backfired:

    mandi_xtal   WT sidechains  min 0.62 A (8 contacts <2 A)
                 QUAD sidechains min 0.64 A (10 contacts <2 A)   <-- wrong
    mandi_dock   WT sidechains  min 0.53 A (17 contacts <2 A)
                 QUAD sidechains min 2.82 A (0 contacts)         <-- clean

The docked arm fits its own lambda=1 endpoint by construction, because the pose
was docked into that model. The crystal arm did not fit EITHER endpoint, because
a crystallographic ligand pose only makes sense against the sidechains it was
solved with. Starting the ceiling arm clashed at both ends makes it worse than
the arm it is supposed to bound.

So the crystal arm now takes both the ligand pose AND the four mutant sidechains
from 4WVO. The arms then differ in pose and rotamer provenance together, which is
the honest contrast: you do not get crystal rotamers without a crystal structure,
so "with a structure" vs "without one" is the comparison that means anything.

The lambda=0 (WT) endpoint stays clashed in both arms, and that is real physics
rather than a setup defect -- wild-type PYR1 genuinely cannot accommodate
mandipropamid, which is the entire reason the quadruple exists.

Backbone comes from the WT frame in both copies so tiMerge can coordinate-match
it; only CB-and-outward is transplanted, and CB is inside the softcore region.
"""
import sys

import numpy as np

AA3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
MUT = {59: ("LYS", "ARG"), 81: ("VAL", "ILE"),
       108: ("PHE", "ALA"), 159: ("PHE", "LEU")}


def read_cif_chainA(path):
    """auth chain A of a PDBx/mmCIF: {(resseq, atomname): (xyz, resname)}."""
    out = {}
    for l in open(path):
        p = l.split()
        if not p or p[0] != "ATOM" or len(p) < 21:
            continue
        if p[18] != "A" or p[17] not in AA3:
            continue
        out[(int(p[16]), p[19])] = (
            np.array([float(p[10]), float(p[11]), float(p[12])]), p[17])
    return out


def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    D = np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))])
    return V @ D @ Wt, P.mean(0), Q.mean(0)


def main():
    wt_f, cif_f, out_f = sys.argv[1], sys.argv[2], sys.argv[3]

    wt = [l for l in open(wt_f) if l.startswith("ATOM")]
    wt_xyz, wt_name = {}, {}
    for l in wt:
        k = (int(l[22:26]), l[12:16].strip())
        wt_xyz[k] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        wt_name[int(l[22:26])] = l[17:20].strip()

    xt = read_cif_chainA(cif_f)

    # superpose 4WVO onto the WT frame on shared BACKBONE atoms
    common = [k for k in xt if k in wt_xyz and k[1] in ("N", "CA", "C", "O")]
    if len(common) < 400:
        raise SystemExit(f"only {len(common)} shared backbone atoms -- numbering suspect")
    P = np.array([xt[k][0] for k in common])
    Q = np.array([wt_xyz[k] for k in common])
    R, mP, mQ = kabsch(P, Q)
    rms = np.sqrt((((P - mP) @ R + mQ - Q) ** 2).sum(1).mean())
    print(f"    superposed 4WVO -> WT frame on {len(common)} backbone atoms, "
          f"RMSD {rms:.2f} A")
    if rms > 1.0:
        raise SystemExit("superposition RMSD too large to transplant rotamers")

    # verify identity at the four positions in BOTH structures before using them
    for nat, (frm, to) in MUT.items():
        got_wt = wt_name.get(nat)
        got_xt = next((v[1] for (r, a), v in xt.items() if r == nat and a == "CA"), None)
        if got_wt != frm or got_xt != to:
            raise SystemExit(f"position {nat}: WT has {got_wt} (want {frm}), "
                             f"4WVO has {got_xt} (want {to})")
        print(f"    position {nat:3d}: WT {frm} / 4WVO {to}  verified")

    out, n_tx = [], 0
    for l in wt:
        r = int(l[22:26])
        atom = l[12:16].strip()
        if r not in MUT:
            out.append(l)
            continue
        if atom not in BACKBONE:
            continue                       # drop the WT sidechain entirely
        out.append(l[:17] + f"{MUT[r][1]:>3}" + l[20:])
    # Append the crystal sidechain (CB outward) using the GLOBAL superposition.
    #
    # A per-residue LOCAL fit was tried and is wrong here. It places CB correctly
    # against the WT backbone (local fit RMSD 0.008-0.034 A) but destroys the
    # ligand-sidechain geometry, because the ligand is positioned against the
    # whole pocket by the global fit: measured, the crystal ligand then clashed
    # with its own crystal rotamers at 1.38 A. The global fit keeps the 4WVO
    # complex internally consistent (2.64 A, no contacts under 2 A), at the cost
    # of moving CB 0.43-1.08 A relative to the WT backbone -- which tleap then
    # propagates into the rebuilt amide H and HA.
    #
    # That is handled in 93 by putting H and HA of the mutated residues INSIDE
    # the softcore region, so they are perturbed atoms rather than common ones
    # and never need coordinate matching. Only N, CA, C, O stay common, and those
    # are copied verbatim from the WT frame, so they match exactly.
    for nat, (frm, to) in MUT.items():
        for (r, a), (xyz, rn) in sorted(xt.items()):
            if r != nat or a in BACKBONE or a.startswith("H"):
                continue
            v = (xyz - mP) @ R + mQ
            # Build the line by SPLICING into a real ATOM line of the same
            # residue rather than formatting one from scratch. Hand-formatting
            # put the chain ID in column 21 instead of 22, which shifted resSeq
            # left, made tleap see a chain break and turned residue 60 into an
            # N-terminus (H1/H2/H3 instead of H). Same shifted-column failure as
            # the stage-1 apo bug -- splicing cannot get the columns wrong.
            tmpl = next(l for l in wt if int(l[22:26]) == nat)
            nm = f" {a:<3s}" if len(a) <= 3 else f"{a:<4s}"
            el = "N" if a.startswith("N") else ("O" if a.startswith("O") else "C")
            out.append(tmpl[:12] + nm + tmpl[16] + f"{to:>3s}" + tmpl[20:30]
                       + f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}"
                       + "  1.00  0.00          " + f"{el:>2s}" + "\n")
            n_tx += 1
    # PDB must stay residue-ordered for tleap
    out.sort(key=lambda l: (int(l[22:26]), l[12:16].strip() not in ("N", "CA", "C", "O")))
    open(out_f, "w").writelines(out)

    # Parse the result back with a real PDB reader at default settings and
    # assert on counts -- a shifted column produces a file that still "looks"
    # fine in a text editor.
    import parmed as pmd
    chk = pmd.load_file(out_f)
    nres = len(chk.residues)
    if nres != len({int(l[22:26]) for l in wt}):
        raise SystemExit(f"{out_f}: reader sees {nres} residues, expected "
                         f"{len({int(l[22:26]) for l in wt})}")
    for nat, (frm, to) in MUT.items():
        seq = sorted({int(l[22:26]) for l in wt}).index(nat)
        got = chk.residues[seq].name
        if got != to:
            raise SystemExit(f"{out_f}: residue {nat} reads back as {got}, not {to}")
    print(f"    transplanted {n_tx} crystal sidechain heavy atoms; "
          f"wrote {out_f.split('/')[-1]} ({len(out)} atoms)")
    print(f"    read back with ParmEd: {nres} residues, all four mutations present")


if __name__ == "__main__":
    main()
