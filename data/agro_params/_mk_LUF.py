
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
raw = Chem.MolFromPDBFile('/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/data/agro_params/LUF_pose.pdb', removeHs=False, sanitize=False)
tmpl = Chem.MolFromSmiles('O=C(Nc1cc(OC(F)(F)C(F)(F)F)cc(Cl)c1Cl)NC(=O)c1c(F)cccc1F')
if raw is None or tmpl is None: sys.exit("read failed")
m = AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
m = Chem.AddHs(m, addCoords=True)
Chem.MolToMolFile(m, '/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/data/agro_params/LUF.mol', kekulize=False)
print(m.GetNumAtoms(), sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1))
