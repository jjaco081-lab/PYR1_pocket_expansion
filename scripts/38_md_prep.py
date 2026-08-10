#!/usr/bin/env python
"""
38_md_prep.py -- build the four WT MD systems (structures + tleap inputs).

Purpose
-------
Establish the WT baseline that every later candidate simulation is compared
against. Candidate MD is uninterpretable without a reference run under identical
conditions, and two further questions need dynamics rather than static models:

  1. Is the satellite lobe really sealed? Seven deposited chains say yes
     (section 16); one FastRelax says no. That is a question about whether a
     side-chain rearrangement is thermally accessible -- only MD answers it.
  2. Every cavity number in this project is a POINT ESTIMATE from one structure.
     MD turns cavity volume into a distribution, so variant differences can
     finally be tested against thermal fluctuation.

Systems
-------
  S1 apo_open      3K3K chain A          ligand-free protomer, OPEN conformer
  S2 holo_closed   3QN1 chain A + ABA    CLOSED conformer, ligand bound
  S3 apo_dimer     3K3K chains A+B       physiological apo homodimer, ABA removed
  S4 ternary       3QN1 A + ABA + HAB1 B the signalling complex (W385 water bridge)

S1 and S2 reuse data/pyr1_open_A.pdb and data/pyr1_closed_A.pdb from script 30 --
identical residue sets, already superposed -- so MD and Arm 2b share inputs.

Force field
-----------
  protein  ff19SB          water  OPC (matched pairing; ff19SB is parameterised
                                  against OPC, not TIP3P)
  ligand   GAFF2 + AM1-BCC on the DEPROTONATED carboxylate, net charge -1,
                           consistent with scripts/00_setup_inputs.py
  ions     Li/Merz 12-6 for OPC; neutralise then 0.15 M NaCl
  Mn2+     present in HAB1 (3QN1 chain B) and retained in S4

Emits data/md/<system>/system.pdb plus a tleap script per system, and the shared
ABA parameters. Run 38b_md_build.sh afterwards to execute antechamber/tleap.

Run with the esmfold2 env python.
"""
import os, sys, warnings, shutil
warnings.filterwarnings("ignore")
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
MD = os.path.join(DATA, "md")
os.makedirs(MD, exist_ok=True)
p = MMCIFParser(QUIET=True)
BB = ("N", "CA", "C", "O")


class Clean(Select):
    """standard residues only, altloc A, no hydrogens, no OXT, no waters."""
    def __init__(self, chains, keep_het=()):
        self.chains = set(chains)
        self.keep_het = set(keep_het)

    def accept_chain(self, c):
        return c.id in self.chains

    def accept_residue(self, r):
        if r.id[0] == " ":
            return is_aa(r, standard=True) and all(a in r for a in BB)
        return r.get_resname() in self.keep_het

    def accept_atom(self, a):
        return (a.element != "H" and a.get_id() != "OXT"
                and a.get_altloc() in (" ", "A"))


def write(cif, chains, out, keep_het=()):
    m = p.get_structure("x", os.path.join(DATA, cif))[0]
    io = PDBIO()
    io.set_structure(m)
    io.save(out, Clean(chains, keep_het))
    return out


def count(path):
    res, het = set(), {}
    for l in open(path):
        if l.startswith("ATOM"):
            res.add((l[21], int(l[22:26])))
        elif l.startswith("HETATM"):
            het[l[17:20].strip()] = het.get(l[17:20].strip(), 0) + 1
    return len(res), het


SYSTEMS = {
    "S1_apo_open":    dict(desc="3K3K chain A, ligand-free protomer (OPEN)",
                           src=("prepared", "pyr1_open_A.pdb"), lig=False, mn=False),
    "S2_holo_closed": dict(desc="3QN1 chain A + ABA (CLOSED)",
                           src=("prepared", "pyr1_closed_A.pdb"), lig=True, mn=False),
    "S3_apo_dimer":   dict(desc="3K3K chains A+B, APO homodimer (ABA removed)",
                           src=("cif", "3K3K.cif", ["A", "B"], ()), lig=False, mn=False),
    "S4_ternary":     dict(desc="3QN1 A + ABA + HAB1 B (+Mn2+)",
                           src=("cif", "3QN1.cif", ["A", "B"], ("MN",)), lig=True, mn=True),
}

print("=" * 78)
print("BUILDING MD SYSTEMS")
print("=" * 78)
for name, cfg in SYSTEMS.items():
    d = os.path.join(MD, name)
    os.makedirs(d, exist_ok=True)
    out = os.path.join(d, "protein.pdb")
    if cfg["src"][0] == "prepared":
        shutil.copy(os.path.join(DATA, cfg["src"][1]), out)
    else:
        _, cif, chains, het = cfg["src"]
        write(cif, chains, out, het)
    nres, hets = count(out)
    print(f"{name:<16}{cfg['desc']:<44}{nres:>4} res  het={hets or '-'}")

# ---- ABA placed in the crystal frame, shared by S2 and S4 ----
aba_src = os.path.join(DATA, "aba_xtal.pdb")
aba_lines = [l for l in open(aba_src)
             if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() != "H"]
aba_out = os.path.join(MD, "aba_xtal_A8S.pdb")
with open(aba_out, "w") as fh:
    for i, l in enumerate(aba_lines, 1):
        fh.write(f"HETATM{i:>5} {l[12:16]} A8S L{1:>4}    {l[30:54]}"
                 f"  1.00  0.00          {l[76:78]}\n")
    fh.write("END\n")
print(f"\nABA (crystal pose, {len(aba_lines)} heavy atoms) -> {aba_out}")
print("  NOTE: S2/S4 use the CLOSED-state crystal pose. S1/S3 are apo by design.")

# ---- tleap inputs ----
TLEAP = """# {desc}
source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.gaff2
loadamberparams frcmod.ions234lm_126_opc
{ligload}
prot = loadpdb protein.pdb
{combine}
{bonds}
# neutralise, then 0.15 M NaCl
addions sys Na+ 0
addions sys Cl- 0
solvateoct sys OPCBOX 12.0
addionsrand sys Na+ {nna} Cl- {ncl}
charge sys
check sys
saveamberparm sys system.prmtop system.inpcrd
savepdb sys system_solvated.pdb
quit
"""

for name, cfg in SYSTEMS.items():
    d = os.path.join(MD, name)
    ligload, combine = "", "sys = prot"
    if cfg["lig"]:
        shutil.copy(aba_out, os.path.join(d, "aba.pdb"))
        ligload = ("loadamberparams ../A8S.frcmod\n"
                   "loadoff ../A8S.lib\n"
                   "lig = loadpdb aba.pdb")
        combine = "sys = combine { prot lig }"
    # 0.15 M NaCl: ~0.0274 ions per water at 0.15 M; estimated from volume later,
    # so use a placeholder count refined by 38b after the first solvate pass.
    open(os.path.join(d, "tleap.in"), "w").write(TLEAP.format(
        desc=cfg["desc"], ligload=ligload, combine=combine,
        bonds="", nna=0, ncl=0))

print(f"\nwrote tleap inputs for {len(SYSTEMS)} systems under {MD}")
print("""
NEXT: bash scripts/38b_md_build.sh
      - parameterises ABA once (antechamber AM1-BCC, GAFF2, net charge -1)
      - runs tleap per system and reports net charge / atom counts
      - the NaCl count is set there, from the actual solvated box volume""")
