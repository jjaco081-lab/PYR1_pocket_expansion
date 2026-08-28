#!/usr/bin/env python
r"""
154_beltran_ddg.py -- protein-only Cartesian FastRelax ddG for every single
substitution at Beltran's 20 design positions.

WHY THESE POSITIONS. §153 established the bar on the Beltran-45 benchmark:
Tian's substitution frequency transfers, reaching 64 % recall at a library of
1e5 against a random null of 15 % +- 13 %. But it SATURATES at 67 %, which is
exactly the ceiling imposed by Tian's vocabulary -- 15 of the 45 sensors use
positions Tian never varied (all of CBDA, all of THC, 4F-MDMB-BUTINACA). No
amount of re-ranking Tian's 144 substitutions reaches them. The headroom is in
POSITIONS, not identities, which is the same conclusion the PFAS arm reached.

So the question this run answers is whether the one scorer that beats chance
(§69, protein-only Cartesian FastRelax ddG, AUC 0.250) can nominate
substitutions at positions no library has tried -- where, by construction, there
is no frequency prior at all. 10 of Beltran's 20 positions are outside Tian's
set and this project has never scored any of them.

WHAT IS COMPUTED. 20 positions x 19 substitutions = 380 single mutants, each
against a WILD-TYPE CONTROL RELAXED IN THE IDENTICAL SHELL. The paired control
is not optional: §66a's first version scored against a raw wild type and
produced -130 to -165 REU of pure repack artefact.

Method is 124's exactly -- Cartesian FastRelax (torsion space leaks 0.98 A into
movemap-frozen residues), MoveMap for minimisation plus TaskFactory for the
packer, 8 A shell, 3 replicates, minimum score taken. Nothing here is retuned;
the point is to apply an already-characterised filter to new positions.

⚠ This is a STABILITY filter, not an affinity predictor (§69, §77). It can say a
substitution is tolerated; it cannot say it creates a sensor. Its role in §153's
benchmark is to supply a ranking where frequency is silent.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_rosetta as LR                                        # noqa: E402
import importlib.util                                           # noqa: E402
_s = importlib.util.spec_from_file_location(
    "v122", os.path.join(HERE, "122_combination_viability.py"))
V = importlib.util.module_from_spec(_s)
_s.loader.exec_module(V)

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "beltran_ddg")
THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
         "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
         "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
         "W": "TRP", "Y": "TYR"}
AA20 = "ACDEFGHIKLMNPQRSTVWY"


def heavy_coords(pose, idxs):
    out = {}
    for i in idxs:
        r = pose.residue(i)
        out[i] = np.array([[r.xyz(j).x, r.xyz(j).y, r.xyz(j).z]
                           for j in range(1, r.nheavyatoms() + 1)])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--nchunks", type=int, default=1)
    ap.add_argument("--nrep", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    bel = json.load(open(os.path.join(ROOT, "data", "beltran", "win_sensors.json")))
    positions = bel["positions"]                     # e.g. "Y120" -> (Y, 120)

    import pyrosetta
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    pyrosetta.init("-mute all -ignore_unrecognized_res -ex1 -ex2aro")
    base = pyrosetta.pose_from_pdb(V.PDB)
    sfxn = pyrosetta.create_score_function("ref2015_cart")

    # ---- assert every position IS the residue Beltran says it is -------
    todo = []
    for p in positions:
        wt, num = p[0], int(p[1:])
        one, i = V.wt_at(base, num)
        assert one == wt, (f"position {num} is {one} in {V.PDB}, Beltran's table "
                           f"says {wt} -- numbering frames disagree, stopping")
        for mu in AA20:
            if mu != wt:
                todo.append((num, wt, mu, i))
    print(f"all {len(positions)} Beltran positions confirmed by identity; "
          f"{len(todo)} single substitutions", flush=True)
    todo = [t for k, t in enumerate(todo) if k % a.nchunks == a.chunk]
    print(f"chunk {a.chunk}/{a.nchunks}: {len(todo)} to score, "
          f"{a.nrep} replicates each + paired wild-type", flush=True)

    def relax_in_shell(pose, allowed, nrep):
        outside = [i for i in range(1, pose.total_residue() + 1)
                   if i not in set(allowed)]
        before = heavy_coords(pose, outside)
        scores, worst = [], 0.0
        for _ in range(nrep):
            p = pose.clone()
            mm = MoveMap()
            mm.set_bb(False)
            mm.set_chi(False)
            for i in allowed:
                mm.set_bb(i, True)
                mm.set_chi(i, True)
            tf, _ = LR.restrict_packing(p, allowed)
            fr = FastRelax(sfxn, 1)
            fr.cartesian(True)
            fr.min_type("lbfgs_armijo_nonmonotone")
            fr.set_movemap(mm)
            fr.set_task_factory(tf)
            fr.apply(p)
            scores.append(float(sfxn(p)))
            after = heavy_coords(p, outside)
            for i in outside:
                if i in before and i in after and len(before[i]) == len(after[i]):
                    worst = max(worst, float(np.abs(before[i] - after[i]).max()))
        return min(scores), scores, worst

    rows = []
    for num, wt, mu, i in todo:
        pose = base.clone()
        MutateResidue(i, THREE[mu]).apply(pose)
        sel = ResidueIndexSelector(str(i))
        nb = NeighborhoodResidueSelector(sel, 8.0, True)
        allowed = [k + 1 for k, b in enumerate(nb.apply(pose)) if b]
        s_m, all_m, dm = relax_in_shell(pose, allowed, a.nrep)
        s_w, all_w, dw = relax_in_shell(base.clone(), allowed, a.nrep)
        rows.append({"sub": f"{wt}{num}{mu}", "pos": num, "wt": wt, "mut": mu,
                     "n_shell": len(allowed), "score": s_m, "wt_same_shell": s_w,
                     "ddG": s_m - s_w, "mut_reps": all_m, "wt_reps": all_w,
                     "mut_spread": float(max(all_m) - min(all_m)),
                     "wt_spread": float(max(all_w) - min(all_w)),
                     "nonshell_drift": float(max(dm, dw))})
        print(f"  {wt}{num}{mu:<3} shell {len(allowed):>3}  "
              f"ddG {s_m - s_w:+8.2f}  spread {max(all_m)-min(all_m):5.2f}  "
              f"drift {max(dm, dw):.2f}", flush=True)
        json.dump({"nrep": a.nrep, "rows": rows},
                  open(os.path.join(OUT, f"chunk_{a.chunk:03d}.json"), "w"), indent=1)
    print(f"\nwrote {len(rows)} rows to {OUT}/chunk_{a.chunk:03d}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
