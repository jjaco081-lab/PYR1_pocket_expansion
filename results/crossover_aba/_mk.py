
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
raw = Chem.MolFromPDBFile('/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/results/crossover_aba/A8S_xtal.pdb', removeHs=False, sanitize=False)
tmpl = Chem.MolFromSmiles('CC1=CC(=O)CC([C@]1(/C=C/C(=CC(=O)O)C)O)(C)C')
if raw is None or tmpl is None:
    sys.exit("read failed")
m = AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
m = Chem.AddHs(m, addCoords=True)
Chem.MolToMolFile(m, '/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion/results/crossover_aba/A8S.mol', kekulize=False)
print(m.GetNumAtoms(), sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1))
