#!/usr/bin/env python
r"""
218_build_input_with_aba.py -- rebuild the RFd3 input WITH ABA in the pocket.

⚠ EVERY RFd3 DESIGN SO FAR WAS SCAFFOLDED AROUND AN EMPTY POCKET. The input
`data/3QN1_complex_auth.pdb` contains 3,674 ATOM records and ZERO HETATM -- ABA
was stripped when it was generated, and nobody noticed. The model therefore had
no physical reason to leave a cavity anywhere, which is a plausible reason the
cavity rate has been so variable between arms (0-60 %).

Jannis: "is it possible to motif scaffold it with ABA? I still would prefer to
use ABA as the test molecule."

Yes. RFd3 exposes `specification.ligand` -- "Ligand name or index to include in
design" -- and rfd3/utils/inference.py:extract_ligand_array pulls the named
residue out of the INPUT atom array and sets it as a FULLY FIXED MOTIF. So ABA
is held at its crystal coordinates while the protein is built around it. That
requires the ligand to be present in the input file, which is what this fixes.

Three reasons this should be better than the apo scaffolding:
  * the model has a physical object to build a pocket around, instead of being
    asked to leave a void for no reason
  * every cavity number stops being a poly-glycine number measured against
    nothing -- the ligand defines what the pocket is FOR
  * designs become directly testable with ABA, the molecule the lab already has

⚠ ABA's CCD code is A8S (confirmed from the AF3 CIF, not assumed). Residue
identities of the PYR1 motif are re-asserted after writing, because this file is
the input every downstream arm depends on.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "data", "3QN1.cif")
OLD = os.path.join(ROOT, "data", "3QN1_complex_auth.pdb")
OUT = os.path.join(ROOT, "data", "3QN1_complex_auth_aba.pdb")
MOTIF = [("A", 34, 40, "HAQRIHA"), ("A", 58, 65, "YKHFIKSC"),
         ("A", 81, 92, "VIVISGLPANTS"), ("A", 111, 121, "IGGEHRLTNYK"),
         ("A", 146, 168, "DMPEGNSEDDTRMFADTVVKLNL")]
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def main():
    cols, rows, inl = {}, [], False
    for line in open(SRC):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".")[1]] = len(cols); inl = True; continue
        if inl:
            if line.startswith(("#", "loop_", "_")):
                if rows:
                    break
                continue
            x = line.split()
            if len(x) >= len(cols):
                rows.append(x)
    g = lambda r, k: r[cols[k]]                                  # noqa: E731
    aba = [r for r in rows if g(r, "label_comp_id") == "A8S"
           and g(r, "auth_asym_id") == "A" and g(r, "type_symbol") != "H"]
    assert aba, "A8S (ABA) not found in 3QN1 chain A"
    print(f"ABA (A8S): {len(aba)} heavy atoms in 3QN1 chain A")

    old = [l for l in open(OLD) if l.startswith("ATOM")]
    print(f"existing input: {len(old)} ATOM records, "
          f"{sum(1 for l in open(OLD) if l.startswith('HETATM'))} HETATM")

    n = max(int(l[6:11]) for l in old)
    with open(OUT, "w") as fh:
        for l in old:
            fh.write(l)
        for r in aba:
            n += 1
            nm = g(r, "label_atom_id")
            fh.write(f"HETATM{n:>5} {nm:<4}{'':1}A8S L{1:>4}{'':1}   "
                     f"{float(g(r,'Cartn_x')):>8.3f}{float(g(r,'Cartn_y')):>8.3f}"
                     f"{float(g(r,'Cartn_z')):>8.3f}  1.00  0.00          "
                     f"{g(r,'type_symbol'):>2}\n")
        fh.write("END\n")

    # ---- assert the file we just wrote ------------------------------------
    seq, het = {}, 0
    for l in open(OUT):
        if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == "A":
            seq[int(l[22:26])] = AA3.get(l[17:20].strip(), "X")
        if l.startswith("HETATM"):
            het += 1
    bad = []
    for ch, lo, hi, want in MOTIF:
        got = "".join(seq.get(i, "-") for i in range(lo, hi + 1))
        if got != want:
            bad.append(f"{ch}{lo}-{hi}: {got} != {want}")
    assert not bad, f"MOTIF IDENTITY BROKEN: {bad}"
    assert het == len(aba), f"wrote {het} HETATM, expected {len(aba)}"
    print(f"wrote {OUT}")
    print(f"  {het} ABA atoms as chain L, residue A8S 1")
    print(f"  all {len(MOTIF)} PYR1 motif segments re-asserted OK")
    # where does ABA sit relative to the motif?
    ap = np.array([[float(g(r, "Cartn_x")), float(g(r, "Cartn_y")),
                    float(g(r, "Cartn_z"))] for r in aba])
    ca = {}
    for l in open(OUT):
        if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == "A":
            ca[int(l[22:26])] = np.array([float(l[30:38]), float(l[38:46]),
                                          float(l[46:54])])
    near = sorted((float(np.linalg.norm(ap - ca[k], axis=1).min()), k)
                  for k in ca)[:6]
    print(f"  nearest PYR1 CA atoms to ABA: "
          f"{', '.join(f'{k} ({d:.1f} Å)' for d, k in near)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
