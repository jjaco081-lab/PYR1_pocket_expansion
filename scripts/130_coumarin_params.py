#!/usr/bin/env python
r"""
130_coumarin_params.py -- Rosetta .params for the 11 coumarin targets.

WHY
§67/§69 scored variants PROTEIN-ONLY, deliberately: the ligand pose is the thing
this project cannot get for free (§46b measured 8.33 A docking error into a pocket
that did not yet exist), so a protein-only score sidesteps the circularity. Jannis
asks the obvious next question -- can we relax WITH the ligand, and iterate
dock/relax/mutate? That needs a ligand residue type per target, which is this.

⚠ DIFFERENT FROM §49, AND THE DIFFERENCE MATTERS. §49 built params for ligands in
their CRYSTAL poses (4WVO's mandipropamid, ABA), so the conformer was the bound
one. No coumarin sensor has a structure. These params therefore carry a GENERATED
conformer, and the pose has to come from docking rather than from the file. Any
result built on them inherits that, and it is the single largest reason to treat
the ligand-aware experiment as a test rather than an upgrade.

Pipeline (same as §49 steps 2-5, minus the crystal-coordinate marriage):
  SMILES -> RDKit 3D embed + MMFF -> molfile (kekulize=False so aromatic bonds
  come out as type 4) -> molfile_to_params.py
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "data", "coumarin_params")
M2P = ("/opt/linux/rocky/8.x/x86_64/pkgs/rosetta/2023.45/scripts/python/public/"
       "molfile_to_params.py")

#: three-letter codes must be unique and not collide with standard residues
COUMARINS = {
    "CM1": ("osthole", "CC(C)=CCc1c(OC)ccc2ccc(=O)oc12"),
    "CM2": ("imperatorin", "CC(C)=CCOc1c2ccoc2cc2ccc(=O)oc12"),
    "CM3": ("psoralen", "O=c1ccc2cc3ccoc3cc2o1"),
    "CM4": ("isopsoralen", "O=c1ccc2ccc3occc3c2o1"),
    "CM5": ("bergapten", "COc1c2ccoc2cc2oc(=O)ccc12"),
    "CM6": ("methoxsalen", "COc1c2occc2cc2ccc(=O)oc12"),
    "CM7": ("citropten", "COc1cc(OC)c2ccc(=O)oc2c1"),
    "CM8": ("scopoletin", "COc1cc2ccc(=O)oc2cc1O"),
    "CM9": ("4-methylumbelliferone", "CC1=CC(=O)Oc2cc(O)ccc12"),
    "CMA": ("5,7-dihydroxy-4-methylcoumarin", "CC1=CC(=O)Oc2cc(O)cc(O)c12"),
    "CMB": ("7-methoxycoumarin", "COc1ccc2ccc(=O)oc2c1"),
}


def main():
    os.makedirs(OUT, exist_ok=True)
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    if not os.path.exists(M2P):
        raise SystemExit(f"molfile_to_params not found at {M2P}")
    ok, bad = [], []
    for code, (name, smi) in COUMARINS.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            bad.append((name, "unparseable SMILES"))
            continue
        m = Chem.AddHs(m)
        if AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE) != 0:
            bad.append((name, "embed failed"))
            continue
        AllChem.MMFFOptimizeMolecule(m)
        mol = os.path.join(OUT, f"{code}.mol")
        # kekulize=False keeps aromatic bonds as type 4; molfile_to_params
        # mistypes aromatic rings if handed a Kekule structure (§49 step 4)
        Chem.MolToMolFile(m, mol, kekulize=False)
        r = subprocess.run([sys.executable, M2P, "-n", code, "-p",
                            os.path.join(OUT, code), "--keep-names", mol],
                           capture_output=True, text=True, cwd=OUT)
        pfile = os.path.join(OUT, f"{code}.params")
        if not os.path.exists(pfile):
            bad.append((name, (r.stdout + r.stderr)[-200:]))
            continue
        nheavy = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1)
        # sanity: the params PDB must contain every heavy atom
        pdb = os.path.join(OUT, f"{code}_0001.pdb")
        npdb = sum(1 for l in open(pdb)
                   if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() != "H") \
            if os.path.exists(pdb) else -1
        flag = "OK" if npdb == nheavy else f"MISMATCH pdb={npdb}"
        ok.append((code, name, nheavy, flag))
        print(f"   {code}  {name:<32} {nheavy:>3} heavy atoms   {flag}")
    print(f"\n   {len(ok)} params written to {OUT}")
    if bad:
        print("   FAILED:")
        for n, e in bad:
            print(f"     {n}: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
