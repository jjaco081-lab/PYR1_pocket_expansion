#!/usr/bin/env python
r"""
167_pose_ladder.py -- how much pose error does the §91 result tolerate?

Jannis: mandipropamid has a crystal structure, so a Boltz pose for it is
contaminated -- 4WVO is in the training set, and "the predicted pose works"
would mostly report memorisation. Correct, and it makes that test ASYMMETRIC:
failure would be conclusive, success would not.

This is the version with no model in the loop at all. Take the CRYSTAL
mandipropamid pose and displace it by a controlled amount -- rigid-body rotation
plus translation, scaled to hit a target heavy-atom RMSD -- then rerun §91's
scan at each level. That measures the quantity actually in question:

    HOW ACCURATE DOES A POSE HAVE TO BE FOR THIS METHOD TO WORK?

and it answers the agrochemical result cleanly:
  * if AUC survives to 2-3 A of pose error, then fludioxonil's 0.71 A seed
    spread was easily good enough, and its failure is CHEMISTRY (no clash) --
    §94a's reading is confirmed and predicted poses are exonerated;
  * if AUC collapses by 1 A, pose accuracy is the binding constraint and §91
    depends on crystallography whatever the chemistry does.

No training-set contamination is possible because no predictor is used.

LEVELS: 0 (the crystal control, must reproduce AUC 0.868), then 0.5, 1.0, 2.0
and 3.0 A, three independent random perturbations each. Perturbations are rigid
-- the ligand's internal geometry is untouched, so the params file stays valid
and only its placement changes.

⚠ A displaced pose creates clashes that the WT pocket did not have. Background
strain is therefore reported at every level: if it explodes, the level is
measuring "ligand rammed into the protein" rather than "pose slightly wrong".
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import lib_rosetta as LR                                            # noqa: E402

OUT = os.path.join(ROOT, "results", "pose_ladder")
FRAME = os.path.join(ROOT, "data", "stage1", "wt_mandi.pdb")
LIGPDB = os.path.join(ROOT, "data", "stage1", "params", "3UZ_0001.pdb")
LIGPRM = os.path.join(ROOT, "data", "stage1", "params", "3UZ.params")
PARK25 = {55: "P", 61: "F", 62: "I", 81: "V", 83: "V", 87: "L", 88: "P", 89: "A",
          92: "S", 94: "E", 141: "E", 108: "F", 110: "I", 115: "H", 116: "R",
          117: "L", 120: "Y", 122: "S", 158: "M", 159: "F", 160: "A", 162: "T",
          163: "V", 164: "V", 167: "N"}
THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
         "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
         "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
         "W": "TRP", "Y": "TYR"}
AA20 = "ACDEFGHIKLMNPQRSTVWY"
LEVELS = [(0.0, 1), (0.5, 3), (1.0, 3), (2.0, 3), (3.0, 3)]


def perturb(X, target, rng):
    """Rigid rotation+translation of X giving heavy-atom RMSD ~= target."""
    if target == 0:
        return X.copy()
    c = X.mean(0)
    for _ in range(200):
        ax = rng.normal(size=3); ax /= np.linalg.norm(ax)
        th = rng.normal() * 0.25
        K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
        R = np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K
        t = rng.normal(size=3)
        t /= np.linalg.norm(t)
        Y = (X - c) @ R.T + c
        # scale the translation so the total RMSD lands on target
        base = float(np.sqrt(((Y - X) ** 2).sum(1).mean()))
        if base >= target:
            Y = (X - c) @ R.T * (target / max(base, 1e-9)) + c
            return Y + t * 0.0
        s = np.sqrt(max(target ** 2 - base ** 2, 0))
        Z = Y + t * s
        got = float(np.sqrt(((Z - X) ** 2).sum(1).mean()))
        if abs(got - target) < 0.05:
            return Z
    return Z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--nchunks", type=int, default=1)
    ap.add_argument("--nrep", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    import pyrosetta
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    from pyrosetta.rosetta.numeric import xyzVector_double_t
    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ex1 -ex2aro "
                   f"-extra_res_fa {LIGPRM}")
    sfxn = pyrosetta.create_score_function("ref2015_cart")
    base0 = pyrosetta.pose_from_pdb(FRAME)
    info = base0.pdb_info()
    idx = {info.number(i): i for i in range(1, base0.total_residue() + 1)}
    lig = base0.total_residue()
    bad = {n: base0.residue(idx[n]).name1() for n, w in PARK25.items()
           if base0.residue(idx[n]).name1() != w}
    assert not bad, f"Park identities disagree: {bad}"
    MutateResidue(idx[59], "ARG").apply(base0)
    L0 = np.array([[base0.residue(lig).xyz(i).x, base0.residue(lig).xyz(i).y,
                    base0.residue(lig).xyz(i).z]
                   for i in range(1, base0.residue(lig).natoms() + 1)])
    e_lig = float(sfxn(pyrosetta.pose_from_pdb(LIGPDB)))

    def dg(p):
        e = float(sfxn(p)); q = p.clone(); q.delete_residue_slow(q.total_residue())
        return e - float(sfxn(q)) - e_lig

    def repack(p0, allowed, nrep):
        v = []
        for _ in range(nrep):
            p = p0.clone()
            tf, _ = LR.restrict_packing(p, allowed)
            pk = PackRotamersMover(sfxn); pk.task_factory(tf); pk.apply(p)
            v.append(dg(p))
        return min(v)

    cells = [(lvl, k) for lvl, n in LEVELS for k in range(n)]
    jobs = [(n, mu) for n in sorted(PARK25) for mu in AA20 if mu != PARK25[n]]
    out = []
    for ci, (lvl, k) in enumerate(cells):
        if ci % a.nchunks != a.chunk:
            continue
        rng = np.random.default_rng(1000 * int(lvl * 10) + k)
        Lp = perturb(L0, lvl, rng)
        got = float(np.sqrt(((Lp - L0) ** 2).sum(1).mean()))
        base = base0.clone()
        for i in range(1, base.residue(lig).natoms() + 1):
            base.residue(lig).set_xyz(i, xyzVector_double_t(*Lp[i - 1]))
        shell, bgc, rows = {}, {}, []
        for num, mu in jobs:
            if num not in shell:
                sel = ResidueIndexSelector(str(idx[num]))
                nb = NeighborhoodResidueSelector(sel, 6.0, True)
                shell[num] = sorted({i + 1 for i, b in enumerate(nb.apply(base)) if b}
                                    | {lig})
                bgc[num] = repack(base.clone(), shell[num], a.nrep)
            p = base.clone()
            MutateResidue(idx[num], THREE[mu]).apply(p)
            s = repack(p, shell[num], a.nrep)
            rows.append({"sub": f"{PARK25[num]}{num}{mu}", "pos": num,
                         "ddG": s - bgc[num]})
        rec = dict(level=lvl, rep=k, actual_rmsd=round(got, 3),
                   bg_median=float(np.median(list(bgc.values()))), rows=rows)
        json.dump(rec, open(os.path.join(OUT, f"L{int(lvl*10):03d}_{k}.json"), "w"))
        print(f"  level {lvl:.1f} A rep {k}: actual RMSD {got:.2f}, "
              f"background dG_bind median {rec['bg_median']:+.1f}, {len(rows)} scored",
              flush=True)
        out.append(rec)
    return 0


if __name__ == "__main__":
    sys.exit(main())
