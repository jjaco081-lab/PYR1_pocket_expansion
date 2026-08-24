#!/usr/bin/env python
r"""
115_us_arm_prep.py -- closed 191-residue frames and ligand poses for the three
extra umbrella arms (quad+mandi, quad+ABA, WT+mandi).

Everything is placed in the md191 frame, because that is where the umbrella
windows and the reaction coordinate live. Three things are produced:

  wt_closed.pdb        WT protein, closed, from an equilibrated S2 frame
  quad_closed.pdb      the same backbone carrying K59R/V81I/F108A/F159L, with
                       4WVO's own rotamers for those four side chains
  aba_in_frame.mol2    ABA from that same S2 frame (already in the right place)
  mandi_in_frame.mol2  mandipropamid superposed into the frame from 4WVO

WHY AN EQUILIBRATED FRAME AND NOT THE CRYSTAL
The umbrella windows are seeded from md191 trajectories, so the reference frame
has to be md191's. Taking a solution-equilibrated closed frame also avoids
handing tleap a crystal geometry that then relaxes differently from every window
it will be compared against.

NUMBERING. md191 is the 191-residue rebuild and uses NATIVE numbering directly
(gate 85-89 = SGLPA, latch 115-117 = HRL), unlike the 178-residue data/md tree
where the gate is at sequential 82-86. Identity is asserted at each of the four
mutated positions before anything is written -- native 81 and 83 are BOTH valine,
so "is it a VAL" is not sufficient (§44d).
"""
import os
import subprocess
import sys

import numpy as np
import parmed as pmd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "data", "md191")
OUT = os.path.join(ROOT, "data", "umbrella_arms")
PARAMS = os.path.join(ROOT, "data", "mmgbsa", "params")
MUT = {59: ("LYS", "ARG"), 81: ("VAL", "ILE"),
       108: ("PHE", "ALA"), 159: ("PHE", "LEU")}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
AA3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    D = np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))])
    return V @ D @ Wt, P.mean(0), Q.mean(0)


def cif_chainA(path):
    out = {}
    for l in open(path):
        p = l.split()
        if not p or p[0] != "ATOM" or len(p) < 21 or p[18] != "A" or p[17] not in AA3:
            continue
        out[(int(p[16]), p[19])] = (np.array([float(p[10]), float(p[11]), float(p[12])]),
                                    p[17])
    return out


def cif_ligand(path, code):
    out = []
    for l in open(path):
        p = l.split()
        if p and p[0] == "HETATM" and len(p) > 20 and p[5] == code:
            out.append((p[3], np.array([float(p[10]), float(p[11]), float(p[12])])))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    # ---- 1. a closed, equilibrated WT frame in the md191 frame -------------
    top = os.path.join(MD, "S2_holo_closed", "system.prmtop")
    rst = os.path.join(MD, "S2_holo_closed", "rep0", "prod.rst7")
    p = pmd.load_file(top, rst)
    xyz = np.array(p.coordinates)
    prot = [r for r in p.residues if r.name not in ("WAT", "HOH", "K+", "Cl-", "A8S")]
    if len(prot) != 191:
        raise SystemExit(f"expected 191 protein residues, got {len(prot)}")
    for nat, (frm, _) in MUT.items():
        got = prot[nat - 1].name
        if got != frm:
            raise SystemExit(f"md191 residue {nat} is {got}, expected {frm}")
    if prot[82].name != "VAL":
        raise SystemExit("native 83 should also be VAL -- the two-valine check")
    print(f"    identities verified at 59/81/108/159; native 83 is VAL as expected")

    def write_pdb(path, residues, skip_sidechain=(), extra=()):
        """Atoms MUST come out in residue order. Appending the transplanted side
        chains after residue 191 made tleap read them as continuing the chain and
        it died with 'Atom .R<THR 191>.A<OXT 15> does not have a type'. Same class
        of defect as 101's column shift: the file looks fine and the parser does
        not."""
        by_res = {}
        for nm, resn, resi, v, el in extra:
            by_res.setdefault(resi, []).append((nm, resn, v, el))
        n = 0
        with open(path, "w") as fh:
            for i, r in enumerate(residues, start=1):
                nm = r.name
                if i in MUT and i in skip_sidechain:
                    nm = MUT[i][1]
                for a in r.atoms:
                    if i in skip_sidechain and a.name not in BACKBONE:
                        continue
                    if a.atomic_number == 1:
                        continue
                    n += 1
                    x, y, z = xyz[a.idx]
                    el = a.element_name if hasattr(a, "element_name") else a.name[0]
                    fh.write(f"ATOM  {n:5d} {a.name:<4s}{nm:>4s} A{i:4d}    "
                             f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          "
                             f"{el:>2s}\n")
                for anm, resn, v, el in by_res.get(i, []):
                    n += 1
                    fh.write(f"ATOM  {n:5d} {anm:<4s}{resn:>4s} A{i:4d}    "
                             f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}  1.00  0.00          "
                             f"{el:>2s}\n")
        return n

    nwt = write_pdb(os.path.join(OUT, "wt_closed.pdb"), prot)
    print(f"    wt_closed.pdb: {nwt} heavy atoms, 191 residues")

    # ---- 2. quad: same backbone, 4WVO's own rotamers ----------------------
    xt = cif_chainA(os.path.join(ROOT, "data", "4WVO.cif"))
    wtca = {}
    for i, r in enumerate(prot, start=1):
        for a in r.atoms:
            if a.name in ("N", "CA", "C", "O"):
                wtca[(i, a.name)] = xyz[a.idx]
    common = [k for k in xt if k in wtca]
    if len(common) < 400:
        raise SystemExit(f"only {len(common)} shared backbone atoms")
    R, mP, mQ = kabsch(np.array([xt[k][0] for k in common]),
                       np.array([wtca[k] for k in common]))
    rms = np.sqrt((((np.array([xt[k][0] for k in common]) - mP) @ R + mQ
                    - np.array([wtca[k] for k in common])) ** 2).sum(1).mean())
    print(f"    4WVO -> md191 closed frame: {len(common)} atoms, RMSD {rms:.2f} A")
    extra = []
    for nat, (frm, to) in MUT.items():
        got = next((v[1] for (r, a), v in xt.items() if r == nat and a == "CA"), None)
        if got != to:
            raise SystemExit(f"4WVO residue {nat} is {got}, expected {to}")
        for (r, a), (v, _) in sorted(xt.items()):
            if r != nat or a in BACKBONE or a.startswith("H"):
                continue
            el = "N" if a.startswith("N") else ("O" if a.startswith("O") else "C")
            extra.append((a, to, nat, (v - mP) @ R + mQ, el))
    nq = write_pdb(os.path.join(OUT, "quad_closed.pdb"), prot,
                   skip_sidechain=set(MUT), extra=extra)
    print(f"    quad_closed.pdb: {nq} heavy atoms, {len(extra)} transplanted rotamer atoms")

    # ---- 3. ligands in the same frame -------------------------------------
    aba = [r for r in p.residues if r.name == "A8S"][0]
    def mol2_with(src, coords, dst):
        lines = open(src).read().split("\n")
        fl, i = False, 0
        for k, l in enumerate(lines):
            if l.startswith("@<TRIPOS>ATOM"):
                fl = True; continue
            if l.startswith("@<TRIPOS>BOND"):
                fl = False
            if fl and l.strip():
                q = l.split()
                x, y, z = coords[i]
                lines[k] = (f"{int(q[0]):>7} {q[1]:<8} {x:9.4f} {y:9.4f} {z:9.4f} "
                            f"{q[5]:<8} {q[6]:>4} {q[7]:<8} {float(q[8]):9.6f}")
                i += 1
        if i != len(coords):
            raise SystemExit(f"{src}: wrote {i} of {len(coords)} atoms")
        open(dst, "w").write("\n".join(lines))
        return i

    n = mol2_with(os.path.join(PARAMS, "A8S.mol2"),
                  [xyz[a.idx] for a in aba.atoms],
                  os.path.join(OUT, "aba_in_frame.mol2"))
    print(f"    aba_in_frame.mol2: {n} atoms, taken from the S2 frame directly")

    # mandipropamid: 3UZ.mol2 lives in the 179-residue stage1 frame, so bring it
    # across with the SAME transform used for the rotamers
    src = os.path.join(PARAMS, "3UZ.mol2")
    raw, fl = [], False
    for l in open(src).read().split("\n"):
        if l.startswith("@<TRIPOS>ATOM"):
            fl = True; continue
        if l.startswith("@<TRIPOS>BOND"):
            break
        if fl and l.strip():
            q = l.split()
            raw.append(np.array([float(q[2]), float(q[3]), float(q[4])]))
    # 3UZ is already superposed into the stage1 WT frame; map stage1 -> md191 on
    # the shared backbone rather than assuming the two frames coincide
    s1 = os.path.join(ROOT, "results", "stage1_rosetta", "_input_wt_aba_A8S_anion.pdb")
    s1ca = {}
    for l in open(s1):
        if l.startswith("ATOM") and l[12:16].strip() in ("N", "CA", "C", "O"):
            s1ca[(int(l[22:26]), l[12:16].strip())] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    com = [k for k in s1ca if k in wtca]
    R2, mP2, mQ2 = kabsch(np.array([s1ca[k] for k in com]),
                          np.array([wtca[k] for k in com]))
    r2 = np.sqrt((((np.array([s1ca[k] for k in com]) - mP2) @ R2 + mQ2
                   - np.array([wtca[k] for k in com])) ** 2).sum(1).mean())
    print(f"    stage1 -> md191: {len(com)} atoms, RMSD {r2:.2f} A")
    n = mol2_with(src, [(v - mP2) @ R2 + mQ2 for v in raw],
                  os.path.join(OUT, "mandi_in_frame.mol2"))
    print(f"    mandi_in_frame.mol2: {n} atoms")


if __name__ == "__main__":
    main()
