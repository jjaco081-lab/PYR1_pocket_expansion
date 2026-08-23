#!/usr/bin/env python
r"""98_check_stereo.py -- assert a rebuilt ligand kept its stereochemistry.

    98_check_stereo.py <reference_sybyl.mol2> <candidate_sybyl.mol2>

97_pose_map.py has to match the docked pose on topology alone, because smina
discards the hydrogen on the stereocentre and RDKit then perceives no chiral
centre to match against. That is safe by construction (Vina rotates bonds, it
never inverts a centre) but "safe by construction" is exactly the kind of claim
this project has learned to check rather than assert. Once the hydrogens are
relaxed back on, the centre is perceptible again -- so compare it here.

Needs RDKit: run under conda_envs/esmfold2/bin/python.
"""
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def describe(path):
    m = Chem.MolFromMol2File(path, removeHs=False, sanitize=True)
    if m is None:
        raise SystemExit(f"RDKit could not parse {path}")
    Chem.AssignStereochemistryFrom3D(m)
    centres = Chem.FindMolChiralCenters(m, useLegacyImplementation=False,
                                        includeUnassigned=True)
    return centres, Chem.MolToSmiles(Chem.RemoveHs(m))


def main():
    ref, cand = sys.argv[1], sys.argv[2]
    c_ref, s_ref = describe(ref)
    c_can, s_can = describe(cand)
    print(f"    reference stereocentres {c_ref}")
    print(f"    rebuilt   stereocentres {c_can}")
    if c_ref != c_can or s_ref != s_can:
        raise SystemExit("    STEREOCHEMISTRY CHANGED -- the docked arm would be "
                         f"the wrong molecule\n      ref  {s_ref}\n      cand {s_can}")
    print("    stereochemistry and canonical SMILES both preserved")


if __name__ == "__main__":
    main()
