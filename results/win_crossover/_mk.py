
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
raw = Chem.MolFromPDBFile('/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/results/win_crossover/WI5_xtal.pdb', removeHs=False, sanitize=False)
tmpl = Chem.MolFromSmiles('Cc1c(c2cccc3c2n1[C@@H](CO3)CN4CCOCC4)C(=O)c5cccc6c5cccc6')
if raw is None or tmpl is None:
    sys.exit("read failed")
m = AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
m = Chem.AddHs(m, addCoords=True)
Chem.MolToMolFile(m, '/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/results/win_crossover/WI5.mol', kekulize=False)
print(m.GetNumAtoms(), sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1))
