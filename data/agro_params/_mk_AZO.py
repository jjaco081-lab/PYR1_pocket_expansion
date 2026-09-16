
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
raw = Chem.MolFromPDBFile('/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/data/agro_params/AZO_pose.pdb', removeHs=False, sanitize=False)
tmpl = Chem.MolFromSmiles('CO/C=C(\\C(=O)OC)c1ccccc1Oc1cc(Oc2ccccc2C#N)ncn1')
if raw is None or tmpl is None: sys.exit("read failed")
m = AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
m = Chem.AddHs(m, addCoords=True)
Chem.MolToMolFile(m, '/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/data/agro_params/AZO.mol', kekulize=False)
print(m.GetNumAtoms(), sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1))
