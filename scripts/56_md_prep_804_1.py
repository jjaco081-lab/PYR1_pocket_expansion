#!/usr/bin/env python
"""
56_md_prep_804_1.py -- build the S5 MD system: PYR1 + ABA + the 804_1 designed binder.

WHAT THIS IS
------------
S4 is the WT signalling complex: PYR1 + ABA + HAB1 + Mn2+. S5 is the same complex
with HAB1 replaced by 804_1, an RFdiffusion binder from the HAB1_shrinkinator
project -- the design that was actually sequenced (4/4 colonies, 2026-08-13). The
pair S4 vs S5 is the comparison this system exists to support: does a designed
binder hold the PYR1 closed state the way HAB1 does?

Source structure: AF3 prediction with ABA, chains A (PYR1, 193 aa), B (804_1,
143 aa), L (A8S). pLDDT 90.5 overall, 96.1 on the ligand. ABA makes 19 contacts
with chain A and ZERO with chain B, which is the expected topology -- HAB1 and its
mimics read PYR1's closed-state surface rather than touching the ligand.

NO Mn2+, unlike S4. The metal in S4 is HAB1's catalytic phosphatase centre. 804_1
is not a phosphatase and the AF3 model has no metal; adding one would be inventing
chemistry the design does not have.

TWO INPUT HAZARDS, BOTH HANDLED HERE
------------------------------------
1. ABA ATOM NAMING. `A8S.lib` uses antechamber's sequential names (O1..O4) while
   both the crystal ligand and this AF3 model use CCD names (O7, O10, O11, O12).
   S4 dodged this by loading the ligand from `A8S.mol2`, which carries the CRYSTAL
   pose -- correct there, wrong here, because we need ABA where AF3 put it in this
   complex. So the ligand is loaded from a PDB whose atoms are renamed to the lib
   convention. The mapping below was derived by nearest-neighbour matching of
   crystal to mol2 coordinates and verified one-to-one at 0.0000 A, not assumed
   from the ordering.

2. N-TERMINAL CONSTRUCT MISMATCH. Chain A begins `GAMASEL...`; the real Y2H
   construct and WT PYR1 begin `MPSEL...` -- no `GA` prefix, and A2P. So residue
   numbering in this model is OFFSET BY 2 from native PYR1 numbering (native K59
   appears here as K61, F108 as F110, latch H115 as H117 -- all confirmed among the
   ABA contacts). The run proceeds as-is per the user's call on 2026-08-13: the
   N-terminus is disordered and far from both the pocket and the interface. Any
   residue index quoted from this trajectory must be converted before it is
   compared with anything in the PYR1_pocket_expansion numbering.

Run with the esmfold2 env python, then 56b to build with tleap.
"""
import os
import shutil
import sys
import warnings

warnings.filterwarnings("ignore")
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "data", "md")
SYS = os.path.join(MD, "S5_804_1_ternary")
CIF = ("/bigdata/cutlerlab/jjaco081/HAB1_shrinkinator/AF3/outputs_fixed_10_13_25/"
       "with_ABA/hab1_50_150_1000_af3_chimera_804_dldesign_1_with_aba/"
       "hab1_50_150_1000_af3_chimera_804_dldesign_1_with_aba_model.cif")

# CCD/crystal name -> A8S.lib name. Verified one-to-one at 0.0000 A against
# A8S.mol2; carbons are already identical, only the oxygens move.
ABA_RENAME = {"O7": "O1", "O10": "O2", "O11": "O3", "O12": "O4"}

BB = ("N", "CA", "C", "O")
PROT_CHAINS = ("A", "B")
LIG_CHAIN = "L"


class CleanProtein(Select):
    """standard residues with a complete backbone; drop H, OXT, altloc != A."""

    def accept_chain(self, c):
        return c.id in PROT_CHAINS

    def accept_residue(self, r):
        return r.id[0] == " " and is_aa(r, standard=True) and all(a in r for a in BB)

    def accept_atom(self, a):
        if a.element == "H" or a.get_id() == "OXT":
            return False
        alt = a.get_altloc()
        return alt in (" ", "A")


os.makedirs(SYS, exist_ok=True)
st = MMCIFParser(QUIET=True).get_structure("af3", CIF)[0]

# ---- protein ----
io = PDBIO()
io.set_structure(st)
prot_pdb = os.path.join(SYS, "protein.pdb")
io.save(prot_pdb, CleanProtein())

n_res = {c: sum(1 for r in st[c] if r.id[0] == " " and is_aa(r, standard=True))
         for c in PROT_CHAINS}
print(f"protein.pdb: chain A {n_res['A']} res (PYR1), chain B {n_res['B']} res (804_1)")

# ---- ligand, renamed to the lib convention ----
lig = [r for ch in st for r in ch if r.get_resname() == "A8S"]
assert len(lig) == 1, f"expected exactly 1 A8S, found {len(lig)}"
lig = lig[0]
heavy = [a for a in lig if a.element != "H"]
assert len(heavy) == 19, f"ABA has {len(heavy)} heavy atoms, expected 19"

def het_line(serial, name, element, resname, chain, resseq, xyz):
    """
    PDB v3.3 HETATM with every field placed by explicit column index.

    Written this way because the hand-rolled f-string version of this exact record
    has now failed twice in this project. In 23h it put `3UZ` in cols 17-19, which
    left `3` in altLoc and made ProDy drop the whole ligand -- silently, invalidating
    a finished result. Here the same omission produced resname `8S` and 19 tleap
    "does not have a type" errors, which at least failed loudly. Column 17 is altLoc
    and must be blank; resName is 18-20.
    """
    x, y, z = xyz
    el = element.strip().upper()
    # atom names: 1-char elements start in col 14, 2-char in col 13
    nm = f"{name:<4s}" if len(el) == 2 or len(name) >= 4 else f" {name:<3s}"
    line = (f"HETATM{serial:5d} {nm}{' '}{resname:>3s} {chain}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}"
            f"          {el:>2s}")
    assert line[0:6] == "HETATM", line
    assert line[16] == " ", f"col 17 (altLoc) must be blank: {line!r}"
    assert line[17:20].strip() == resname, f"resName misplaced: {line[17:20]!r}"
    assert line[21] == chain, f"chainID misplaced: {line[21]!r}"
    assert line[12:16].strip() == name, f"atom name misplaced: {line[12:16]!r}"
    assert line[76:78].strip() == el, f"element misplaced: {line[76:78]!r}"
    return line + "\n"


aba_pdb = os.path.join(SYS, "aba.pdb")
with open(aba_pdb, "w") as fh:
    for i, a in enumerate(heavy, start=1):
        nm = ABA_RENAME.get(a.get_id(), a.get_id())
        fh.write(het_line(i, nm, a.element, "A8S", "L", 1, a.coord))
    fh.write("END\n")

# validate the rename against what tleap will actually look for
lib_names = set()
with open(os.path.join(MD, "A8S.lib")) as fh:
    inblk = False
    for l in fh:
        if l.startswith("!entry.A8S.unit.atoms table"):
            inblk = True
            continue
        if inblk:
            if l.startswith("!"):
                break
            lib_names.add(l.split()[0].strip('"'))
written = {l[12:16].strip() for l in open(aba_pdb) if l.startswith("HETATM")}
missing = written - lib_names
assert not missing, (
    f"aba.pdb has atoms absent from A8S.lib: {sorted(missing)}. tleap would "
    f"silently build a broken ligand.")

# Re-parse with a real PDB parser at DEFAULT settings and assert on counts, rather
# than trusting that the file we just wrote says what we meant. This is the check
# that would have caught the 23h apo bug at write time.
from Bio.PDB import PDBParser  # noqa: E402

back = PDBParser(QUIET=True).get_structure("aba", aba_pdb)[0]
bres = [r for ch in back for r in ch]
assert len(bres) == 1, f"aba.pdb re-parses as {len(bres)} residues, expected 1"
assert bres[0].get_resname() == "A8S", (
    f"aba.pdb re-parses with resname {bres[0].get_resname()!r}, not 'A8S' -- "
    f"tleap will report 'does not have a type' for every atom")
assert len(list(bres[0])) == 19, f"re-parsed {len(list(bres[0]))} atoms, expected 19"
assert {a.get_id() for a in bres[0]} == written, "atom names changed on round-trip"

print(f"aba.pdb: {len(written)} heavy atoms, all present in A8S.lib "
      f"(lib also supplies {len(lib_names)-len(written)} hydrogens)")
print(f"  round-trip: resname {bres[0].get_resname()}, "
      f"{len(list(bres[0]))} atoms, names match")

# ---- disulfides: tleap needs explicit bonds, and will silently not make them ----
cys = [r for c in PROT_CHAINS for r in st[c]
       if r.get_resname() == "CYS" and "SG" in r]
ss = []
for i in range(len(cys)):
    for j in range(i + 1, len(cys)):
        d = np.linalg.norm(cys[i]["SG"].coord - cys[j]["SG"].coord)
        if d < 2.5:
            ss.append((cys[i], cys[j], float(d)))
print(f"cysteines: {len(cys)}; disulfides within 2.5 A: {len(ss)}")
for a, b, d in ss:
    print(f"  {a.get_parent().id}{a.id[1]} - {b.get_parent().id}{b.id[1]}  {d:.2f} A")

if ss:
    print("  NOTE: these need explicit `bond` lines and CYX residue names in tleap.")

# ---- tleap ----
# Ligand comes from aba.pdb (the AF3 pose in THIS complex), not from A8S.mol2
# (the crystal pose), which is the whole reason for the rename above.
tleap = f"""# S5: PYR1 + ABA + 804_1 designed binder (AF3 model, no Mn2+)
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.gaff2
# NO `loadamberparams frcmod.ions234lm_126_opc` here. That is TIP3P-era naming, no
# such file ships for OPC, and tleap exits "Could not open file". leaprc.water.opc
# has already loaded frcmod.ionslm_126_opc. The template inside 38_md_prep.py still
# carries the bad line; S4's working tleap.in does not, because 38b rewrote it.
loadamberparams ../A8S.frcmod
loadoff ../A8S.lib
lig = loadpdb aba.pdb
prot = loadpdb protein.pdb
sys = combine {{ prot lig }}
addions sys Na+ 0
addions sys Cl- 0
solvateoct sys OPCBOX 12.0
addionsrand sys Na+ {{nna}} Cl- {{ncl}}
charge sys
check sys
saveamberparm sys system.prmtop system.inpcrd
savepdb sys system_solvated.pdb
quit
"""
open(os.path.join(SYS, "tleap.in.template"), "w").write(tleap)
print(f"\nwrote {SYS}/")
print("  protein.pdb  aba.pdb  tleap.in.template")
print("\nNEXT: bash scripts/56b_md_build_804_1.sh   (two-pass tleap, sets NaCl from")
print("      the actual water count, then reports net charge and atom totals)")
