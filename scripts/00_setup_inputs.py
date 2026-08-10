#!/usr/bin/env python
"""
00_setup_inputs.py -- prepare all structural inputs from PDB 3QN1.

Outputs (into ../data):
  pyr1_A.pdb          PYR1 chain A, protein only, apo (closed conformation)
  hab1_B.pdb          HAB1 chain B, protein only
  complex_AB.pdb      PYR1 + HAB1, protein only
  aba_xtal.pdb        ABA (HET code A8S) in its crystallographic pose
  aba_params/A8S.params + A8S_0001.pdb   Rosetta ligand definition, crystal pose

Ligand protonation
------------------
The RCSB "ideal" SDF for A8S is the NEUTRAL carboxylic acid. In 3QN1 the ABA
carboxylate forms a salt bridge with PYR1 K59 NZ at 2.85 A, so the bound state
is the DEPROTONATED carboxylate (-1). Since K59 is one of the positions this
study mutates, getting this charge right is necessary for the K59 ddG to mean
anything. We therefore:
  - read the crystal heavy-atom coordinates from 3QN1,
  - assign bond orders from the ideal-SDF template (RDKit),
  - deprotonate the carboxylic acid (SMARTS [CX3](=O)[OX2H1]),
  - add hydrogens with coordinates,
  - hand the result to molfile_to_params.py.
molfile_to_params emits A8S_0001.pdb with Rosetta-consistent atom names AND the
input coordinates, so that file is used directly as the ligand in the complex --
this sidesteps any atom-name mismatch between PDB and Rosetta conventions.

Run with the esmfold2 env python (has rdkit + biopython).
"""
import os, subprocess, sys, warnings
warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
CIF = os.path.join(DATA, "3QN1.cif")
MOLFILE_TO_PARAMS = ("/opt/linux/rocky/8.x/x86_64/pkgs/rosetta/2023.45/"
                     "scripts/python/public/molfile_to_params.py")

model = MMCIFParser(QUIET=True).get_structure("3qn1", CIF)[0]
io = PDBIO()


class ProteinChain(Select):
    def __init__(self, chains): self.chains = set(chains)
    def accept_residue(self, r):
        return is_aa(r) and r.get_parent().id in self.chains
    def accept_atom(self, a):
        return a.element != "H" and a.get_altloc() in (" ", "A")


io.set_structure(model)
io.save(os.path.join(DATA, "pyr1_A.pdb"), ProteinChain("A"))
io.save(os.path.join(DATA, "hab1_B.pdb"), ProteinChain("B"))
io.save(os.path.join(DATA, "complex_AB.pdb"), ProteinChain("AB"))


class LigOnly(Select):
    def accept_residue(self, r): return r.get_resname() == "A8S"
    def accept_atom(self, a): return a.element != "H"


io.save(os.path.join(DATA, "aba_xtal.pdb"), LigOnly())
print("wrote pyr1_A.pdb, hab1_B.pdb, complex_AB.pdb, aba_xtal.pdb")

# ---- ligand: crystal pose, correct bond orders, deprotonated carboxylate ----
template = Chem.MolFromMolFile(os.path.join(DATA, "A8S_ideal.sdf"), removeHs=True)
lig = Chem.MolFromPDBFile(os.path.join(DATA, "aba_xtal.pdb"), removeHs=True,
                          sanitize=False)
lig = AllChem.AssignBondOrdersFromTemplate(template, lig)

patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1H0-]")
m = Chem.RWMol(lig)
hits = m.GetSubstructMatches(patt)
assert len(hits) == 1, f"expected one carboxylic acid, found {len(hits)}"
o_idx = hits[0][2]
o = m.GetAtomWithIdx(o_idx)
o.SetFormalCharge(-1)
o.SetNoImplicit(True)
o.SetNumExplicitHs(0)
Chem.SanitizeMol(m)
lig = m.GetMol()
print(f"deprotonated carboxylate at atom idx {o_idx}; "
      f"net formal charge {Chem.GetFormalCharge(lig)}")

ligH = Chem.AddHs(lig, addCoords=True)
sdf = os.path.join(DATA, "aba_deprot.sdf")
w = Chem.SDWriter(sdf); w.write(ligH); w.close()
print(f"wrote {sdf} ({ligH.GetNumAtoms()} atoms incl. H)")

# ---- Rosetta params ----
pdir = os.path.join(DATA, "aba_params")
os.makedirs(pdir, exist_ok=True)
r = subprocess.run([sys.executable, MOLFILE_TO_PARAMS, "-n", "A8S",
                    "--clobber", sdf], cwd=pdir, capture_output=True, text=True)
print(r.stdout.strip()[-600:])
if r.returncode:
    print("STDERR:", r.stderr[-800:]); sys.exit(1)

# ---- assemble complex + ligand for the ratchet check ----
lig_lines = [l for l in open(os.path.join(pdir, "A8S_0001.pdb"))
             if l.startswith(("ATOM", "HETATM"))]
lig_lines = [l[:21] + "X" + l[22:] for l in lig_lines]   # ligand on chain X
out = os.path.join(DATA, "complex_AB_ABA.pdb")
with open(out, "w") as fh:
    for l in open(os.path.join(DATA, "complex_AB.pdb")):
        if l.startswith("ATOM"):
            fh.write(l)
    fh.write("TER\n")
    fh.writelines(lig_lines)
    fh.write("TER\nEND\n")
print(f"wrote {out} ({len(lig_lines)} ligand atoms)")
