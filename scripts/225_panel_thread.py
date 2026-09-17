#!/usr/bin/env python
r"""
225_panel_thread.py -- receptor frames for the MM-GBSA panel: WT and the known
hit, for each of the four control ligands.

The comparison MM-GBSA has to make is ddG(hit) - ddG(WT) for the same ligand in
the same pose. That needs two receptor frames per ligand, differing only by the
hit substitutions, with the ligand held fixed.

  anthrone        F108M F159M
  3-phenylphenol  A89M
  chloroxylenol   V83I L87M N167Q
  eugenol         V164L N167V          [HELD OUT -- built, not looked at]

⚠ ASSERT THE RESIDUE IDENTITY, NEVER THE NUMBER. Every substitution names the
wild-type amino acid it replaces; this script checks the residue AT that number
IS that amino acid before mutating, and refuses the whole system otherwise. That
check is what caught the -2 auth/label offset elsewhere in this project.

⚠ Repacking is restricted with a TaskFactory via lib_rosetta.restrict_packing.
FastRelax.set_movemap() restricts MINIMISATION ONLY and has caused two bugs here
that backbone RMSD could not see.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = os.path.join(ROOT, "data", "panel_mmgbsa")
THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
         "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
         "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
         "W": "TRP", "Y": "TYR"}


def main():
    import pyrosetta
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    import lib_rosetta as LR
    pyrosetta.init("-mute all -ex1 -ex2aro -use_input_sc", silent=True)

    systems = json.load(open(os.path.join(OUT, "systems.json")))
    frame = os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")
    if not os.path.exists(frame):
        print(f"no receptor frame at {frame}"); return 1
    pose0 = pyrosetta.pose_from_pdb(frame)
    pdbi = pose0.pdb_info()
    # auth number -> pose index, for chain A only
    rmap = {pdbi.number(i): i for i in range(1, pose0.total_residue() + 1)
            if pdbi.chain(i) == "A"}
    print(f"receptor frame {frame}: {pose0.total_residue()} residues, "
          f"{len(rmap)} in chain A")

    sf = pyrosetta.create_score_function("ref2015")
    made = []
    for s in systems:
        lig = s["ligand"]
        d = os.path.join(OUT, lig.replace("-", ""))
        wt_path = os.path.join(d, "receptor_WT.pdb")
        hit_path = os.path.join(d, "receptor_HIT.pdb")
        pose0.dump_pdb(wt_path)

        p = pose0.clone()
        touched, bad = [], []
        for sub in s["hit"]:
            wt, num, mt = sub[0], int(sub[1:-1]), sub[-1]
            if num not in rmap:
                bad.append(f"{sub}: residue {num} absent from chain A"); continue
            i = rmap[num]
            got = p.residue(i).name3()
            # ⚠ IDENTITY, not number
            if got != THREE[wt]:
                bad.append(f"{sub}: residue {num} IS {got}, expected {THREE[wt]}")
                continue
            MutateResidue(i, THREE[mt]).apply(p)
            touched.append(i)
        if bad:
            print(f"!! {lig}: REFUSED -- " + "; ".join(bad))
            continue
        # repack the shell around the mutations, design nothing
        from pyrosetta.rosetta.core.select.residue_selector import (
            ResidueIndexSelector, NeighborhoodResidueSelector)
        sel = ResidueIndexSelector(",".join(str(i) for i in touched))
        sh = NeighborhoodResidueSelector(sel, 8.0, True).apply(p)
        allowed = [i for i in range(1, p.total_residue() + 1) if sh[i]]
        tf, packable = LR.restrict_packing(p, allowed)
        from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
        pk = PackRotamersMover(sf)
        pk.task_factory(tf)
        pk.apply(p)
        p.dump_pdb(hit_path)
        # assert on the RESULT: the mutations are actually present
        pi = p.pdb_info()
        rm2 = {pi.number(i): i for i in range(1, p.total_residue() + 1)
               if pi.chain(i) == "A"}
        for sub in s["hit"]:
            mt, num = sub[-1], int(sub[1:-1])
            got = p.residue(rm2[num]).name3()
            assert got == THREE[mt], f"{lig} {sub}: after mutation residue is {got}"
        print(f"{lig:<16} {'+'.join(s['hit']):<22} "
              f"{len(touched)} mutations, {len(packable)} repacked  -> {hit_path}")
        made.append(dict(ligand=lig, resname=s["resname"], mol=s["mol"],
                         wt=wt_path, hit=hit_path, subs=s["hit"],
                         uM=s["uM"], holdout=s["holdout"]))
    json.dump(made, open(os.path.join(OUT, "frames.json"), "w"), indent=1)
    print(f"\n{len(made)}/{len(systems)} ligands have both frames "
          f"-> {2*len(made)} systems, x n=8 = {16*len(made)} runs")
    return 0 if len(made) == len(systems) else 1


if __name__ == "__main__":
    sys.exit(main())
