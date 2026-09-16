
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
raw = Chem.MolFromPDBFile('/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/data/agro_params/FLD_pose.pdb', removeHs=False, sanitize=False)
tmpl = Chem.MolFromSmiles('FC1(F)Oc2cccc(c2O1)-c1c[nH]cc1C#N')
if raw is None or tmpl is None: sys.exit("read failed")
m = AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
m = Chem.AddHs(m, addCoords=True)
Chem.MolToMolFile(m, '/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/data/agro_params/FLD.mol', kekulize=False)
print(m.GetNumAtoms(), sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1))
