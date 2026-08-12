#!/usr/bin/env python
"""
47_stage1_inputs.py -- build the stage-1 retrospective test: can a sequence
designer recover PYR1^MANDI from wild-type PYR1 plus a perfectly placed ligand?

GROUND TRUTH
------------
Comparing 4WVO (PYR1^MANDI + mandipropamid) to 3QN1 (WT + ABA) over the 174
residues resolved in both gives exactly four substitutions:

    K59R   V81I   F108A   F159L

All four lie inside Tian's randomised positions, and -- critically -- all four are
reachable ONLY in the DSM-Hao and TSM libraries. Under the Coumarin alphabet
F108A is not offered at all (position 108 allows W only), so constraining this
test to Coumarin would make it unwinnable by construction. DSM-Hao is used here.

That same fact reframes the coumarin wetlab data: F108W appears in 11/11 sequences
partly because W was the only option the library provided, not purely because it
was preferred.

WHAT THIS BUILDS
----------------
  wt_mandi.pdb      WT PYR1 backbone + mandipropamid in its TRUE 4WVO pose,
                    superposed into the 3QN1 frame. The ligand clashes with WT
                    side chains -- that is the point; relieving those clashes is
                    the design problem, and it is why the four mutations exist.
  polygly_mandi.pdb same, with the 15 designable positions truncated to glycine.
                    No clash signal, so the pocket must be rebuilt rather than
                    merely trimmed. The harder and more honest arm.
  wt_aba.pdb        WT PYR1 + ABA, the NEGATIVE CONTROL. The correct answer here
                    is ZERO mutations, because WT is the evolved solution. Any
                    method proposing substitutions for ABA in WT is revealing a
                    bias toward mutating regardless of ligand, which would make
                    its mandipropamid "recovery" uninterpretable.

DESIGNABLE POSITIONS
--------------------
Tian randomises 18 positions, of which 87, 89 (gate) and 117 (latch) are excluded
here to preserve this project's core constraint of not altering the switch loops.
None of the four ground-truth mutations sits at those positions, so excluding them
costs nothing on this test. Set INCLUDE_GATE_LATCH = True to run the alternative
arm, which is worth doing as a diagnostic: a method that mutates the gate is
proposing to break transduction, and we would want to see that.

TRIVIAL BASELINE
----------------
Also reports per-residue steric overlap with the ligand. If simply ranking
residues by clash recovers F108A, then a learned model recovering it demonstrates
nothing -- the same trap the ligand-blind oracle exposed in section 23b. Any
method must be scored against this, not against zero.

Usage:  python 47_stage1_inputs.py
"""
import os
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Superimposer, Select
import warnings; warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(DATA, "stage1")
os.makedirs(OUT, exist_ok=True)

TRUE = {59: ("K", "R"), 81: ("V", "I"), 108: ("F", "A"), 159: ("F", "L")}
TIAN_POS = [59, 81, 83, 87, 89, 92, 94, 108, 110, 117, 120, 122, 141, 159, 160,
            163, 164, 167]
GATE_LATCH = {87, 89, 117}
INCLUDE_GATE_LATCH = False
DESIGN = [p for p in TIAN_POS if INCLUDE_GATE_LATCH or p not in GATE_LATCH]

# DSM-Hao allowed identities (Tian sd04), designable positions only
DSM_HAO = {
    59: "ADEFGHILMNQRSTVWY", 81: "ILMRTY", 83: "AFGILMSTVWY",
    92: "ADEFGHIKLMNQRSTVWY", 94: "ADFGHIKLMNQRSTVWY",
    108: "ADEGHIKLMNQRSTVWY", 110: "AFGILMSTVWY", 120: "ADEFGHIKLMNQRSTVW",
    122: "ADEFGHIKLMNQRTVWY", 141: "ADEFGHIKLMNQRSTVW",
    159: "ADEGHIKLMNQRSTWY", 160: "DEFGHIKLMNQRSTVWY",
    163: "ADEFGHIKLMNQRSTWY", 164: "ADEFGHIKLMNQRSTWY",
    167: "ADEFGHIKLMQRSTVWY",
}

P = MMCIFParser(QUIET=True)
wt = P.get_structure("wt", os.path.join(DATA, "3QN1.cif"))[0]
mn = P.get_structure("mn", os.path.join(DATA, "4WVO.cif"))[0]


def ca(model, ch):
    return {r.id[1]: r["CA"] for r in model[ch] if r.id[0] == " " and "CA" in r}


a, b = ca(wt, "A"), ca(mn, "A")
common = sorted(set(a) & set(b))
si = Superimposer()
si.set_atoms([a[i] for i in common], [b[i] for i in common])
print(f"superposed 4WVO onto 3QN1 on {len(common)} CA, RMSD {si.rms:.2f} A")

lig = [r for ch in mn for r in ch if r.get_resname() == "3UZ"][0]
si.apply(list(lig.get_atoms()))          # mandipropamid now in the 3QN1 frame

aba = [r for ch in wt for r in ch if r.get_resname() == "A8S"][0]


class Keep(Select):
    """WT chain A protein; optionally truncate DESIGN positions to glycine."""
    def __init__(self, polygly=False):
        self.polygly = polygly

    def accept_chain(self, c):
        return c.id == "A"

    def accept_residue(self, r):
        return r.id[0] == " "

    def accept_atom(self, at):
        r = at.get_parent()
        if self.polygly and r.id[1] in DESIGN:
            return at.get_id() in ("N", "CA", "C", "O")
        return at.element != "H"


io = PDBIO()


def write(path, polygly, extra):
    io.set_structure(wt)
    io.save(path, Keep(polygly))
    with open(path) as fh:
        body = [l for l in fh if l.startswith("ATOM")]
    n = int(body[-1][6:11])
    with open(path, "w") as fh:
        fh.writelines(body)
        for at in extra:
            n += 1
            x, y, z = at.coord
            fh.write(f"HETATM{n:5d} {at.get_id():<4s}{extra.get_resname():>3s} A 900    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          "
                     f"{at.element:>2s}\n")
        fh.write("END\n")
    print(f"  wrote {os.path.basename(path)}")


write(os.path.join(OUT, "wt_mandi.pdb"), False, lig)
write(os.path.join(OUT, "polygly_mandi.pdb"), True, lig)
write(os.path.join(OUT, "wt_aba.pdb"), False, aba)

# ---------------- trivial clash baseline ----------------
lx = np.array([at.coord for at in lig if at.element != "H"])
print(f"\nTRIVIAL BASELINE -- steric overlap of each WT side chain with mandipropamid")
print(f"{'pos':>5} {'wt':>3} {'min dist':>9} {'atoms<4A':>9} {'truth':>7}")
rows = []
for r in wt["A"]:
    if r.id[0] != " " or r.id[1] not in DESIGN:
        continue
    sc = [at for at in r if at.get_id() not in ("N", "CA", "C", "O")
          and at.element != "H"]
    if not sc:
        continue
    d = np.linalg.norm(np.array([at.coord for at in sc])[:, None, :] - lx[None, :, :],
                       axis=-1)
    rows.append((r.id[1], r.get_resname(), float(d.min()), int((d < 4.0).sum())))
for pos, rn, dmin, n4 in sorted(rows, key=lambda x: x[2]):
    t = f"->{TRUE[pos][1]}" if pos in TRUE else ""
    print(f"{pos:>5} {rn:>3} {dmin:>9.2f} {n4:>9} {t:>7}")

hits = [p for p, _, dmin, _ in sorted(rows, key=lambda x: x[2])[:4]]
print(f"\n  4 most-clashing positions: {sorted(hits)}")
print(f"  ground truth positions   : {sorted(TRUE)}")
print(f"  overlap: {len(set(hits) & set(TRUE))}/4  <-- any method must beat this")

# ---------------- chance level under DSM-Hao ----------------
print(f"\nCHANCE LEVEL -- uniform draws from the DSM-Hao alphabet, {len(DESIGN)} positions")
print(f"{'N seqs':>7} {'E[true subs recovered]':>24} {'P(all 4)':>10}")
for N in (1, 5, 10, 25, 50, 100):
    ps = []
    for p, (w, t) in TRUE.items():
        alpha = set(DSM_HAO.get(p, "")) | {w}
        ps.append(1 - (1 - 1.0 / len(alpha)) ** N if t in alpha else 0.0)
    print(f"{N:>7} {sum(ps):>24.2f} {np.prod(ps):>10.3f}")
print("\n  recovery must be scored against these, not against zero.")
