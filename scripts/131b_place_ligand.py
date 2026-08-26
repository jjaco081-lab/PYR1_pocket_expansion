#!/usr/bin/env python
r"""
131b_place_ligand.py -- graph-match a docked pose onto its parameterised molecule.

    131b_place_ligand.py <PARAMS_DIR> <CODE> <docked.sdf> <out.pdb>

Runs under the RDKit environment; 131 runs under PyRosetta and the two are not
installed together, so this is invoked as a subprocess rather than imported.

⚠ Why graph matching and not a coordinate copy (§97): smina renames the ligand to
UNL, strips nonpolar hydrogens and REORDERS atoms -- template C1 C2 C3 C4 C5 O1
comes back as C C C C O C. Copying by order silently puts an oxygen's coordinates
on a carbon, and the system would build, run and return a confident number for the
wrong molecule.
"""
import os
import sys


def main():
    params, code, sdf, out = sys.argv[1:5]
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    tmpl = Chem.MolFromMolFile(os.path.join(params, f"{code}.mol"), removeHs=False)
    dock = Chem.MolFromMolFile(sdf, removeHs=False)
    if tmpl is None or dock is None:
        sys.exit("could not read template or docked molecule")
    th = Chem.RemoveHs(Chem.Mol(tmpl))
    dh = Chem.RemoveHs(Chem.Mol(dock))
    if th.GetNumAtoms() != dh.GetNumAtoms():
        sys.exit(f"heavy-atom count differs: template {th.GetNumAtoms()}, "
                 f"docked {dh.GetNumAtoms()}")
    match = th.GetSubstructMatch(dh)
    if not match or len(match) != dh.GetNumAtoms():
        sys.exit("no one-to-one graph match")
    names = [l[12:16] for l in open(os.path.join(params, f"{code}_0001.pdb"))
             if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() != "H"]
    if len(names) != th.GetNumAtoms():
        sys.exit(f"{len(names)} names vs {th.GetNumAtoms()} heavy atoms")
    conf = dh.GetConformer()
    pos = {}
    for j in range(dh.GetNumAtoms()):
        p = conf.GetAtomPosition(j)
        pos[match[j]] = (p.x, p.y, p.z)
    if len(pos) != len(names):
        sys.exit("match is not one-to-one")
    with open(out, "w") as fh:
        for k, nm in enumerate(names):
            x, y, z = pos[k]
            el = nm.strip()[0]
            fh.write(f"HETATM{k+1:5d} {nm}{code:>4} X   1    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2}\n")
        fh.write("END\n")
    print(len(names))
    return 0


if __name__ == "__main__":
    sys.exit(main())
