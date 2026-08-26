#!/usr/bin/env python
r"""
124_viability_relax.py -- the FastRelax version of §66a.

WHY REDO IT
§66a scored 670 variants with fixed-backbone repacking and found nothing: REAL vs
LIBRARY AUC 0.542, REAL vs WILD 0.443. Two objections survive that result and both
point the same way -- the protocol was too rigid to be fair to a real sensor.

  1. A sensor works by REMODELLING the cavity. With the backbone frozen, any
     combination that needs a small backbone shift to accommodate its new side
     chains is scored in a conformation it never adopts, and the strain it is
     charged for is an artefact of the constraint.
  2. The score distribution was dominated by an unrelieved-clash tail (23-32 % of
     every set above +200 REU, sd ~350-400 REU). That variance swamps whatever
     signal exists; relaxation is exactly what removes it.

So this repeats the identical experiment -- same three sets, same variants, same
paired control -- with FastRelax over a movemap that allows BACKBONE and side-chain
motion inside the mutated shell.

WHAT IS HELD IDENTICAL TO §66a, so the two are comparable
  * the same 8 A neighbourhood shell around the mutated positions
  * the same paired wild-type control, relaxed in the SAME shell with the SAME
    protocol and the same number of replicates -- §66a's first version scored the
    repacked mutant against a RAW wild-type and made every variant look 130-165
    REU better than wild-type, which was input strain being relieved
  * the same three sets: REAL (sd07 sensors), LIBRARY (Tian's menu), WILD
    (DSM-Hao menu at the same positions)

WHAT IS NEW, AND THE PART THAT NEEDS CARE
FastRelax is STOCHASTIC. A single trajectory's score is not the variant's score, so
each structure is relaxed N_REP times and the MINIMUM is taken (standard practice:
the relax is a search, and a bad trajectory is a search failure, not evidence of
instability -- [[feedback_score_generative_predictions_bestofn]]). The spread across
replicates is recorded for every variant, because if the within-variant spread is
comparable to the between-set difference then the experiment cannot answer the
question and must say so rather than quoting a mean.

⚠ Both movemap AND task factory are set. `FastRelax.set_movemap()` restricts
minimisation only -- it does NOT stop the packer redesigning or repacking the whole
pose (§28d cost two silent bugs to this). A non-shell drift check runs afterwards.

⚠ CARTESIAN, not torsion-space, and the drift check is why. Torsion-space FastRelax
with the backbone free only inside the shell still moves everything downstream of
it: a phi/psi change at residue i rotates the entire chain after i, so residues the
movemap "froze" were measured drifting **0.98 A**. That is a lever arm, not
packing, and it would have differed systematically between sets because the shells
differ. In Cartesian space atoms outside the movemap genuinely do not move, which
the same check now confirms at <0.01 A. Scored with `ref2015_cart` throughout;
mutant and paired wild-type use the identical function, so the difference is
comparable even though the absolute numbers are not ref2015's.

Usage: 124_viability_relax.py --chunk i --nchunks n [--nrep 3]
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_rosetta as LR                                       # noqa: E402
_v = __import__("importlib").import_module("importlib.util")
import importlib.util                                          # noqa: E402
_s = importlib.util.spec_from_file_location(
    "v122", os.path.join(HERE, "122_combination_viability.py"))
V = importlib.util.module_from_spec(_s)
_s.loader.exec_module(V)

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "viability_relax")
THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
         "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
         "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
         "W": "TRP", "Y": "TYR"}


def heavy_coords(pose, idxs):
    out = {}
    for i in idxs:
        r = pose.residue(i)
        out[i] = np.array([[r.xyz(j).x, r.xyz(j).y, r.xyz(j).z]
                           for j in range(1, r.natoms() + 1)
                           if r.atom_type(j).element() != "H"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--nchunks", type=int, default=1)
    ap.add_argument("--nrep", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    import pyrosetta
    from pyrosetta.rosetta.core.scoring import get_score_function
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    pyrosetta.init("-mute all -ignore_unrecognized_res -ex1 -ex2aro")
    base = pyrosetta.pose_from_pdb(V.PDB)
    sfxn = pyrosetta.create_score_function("ref2015_cart")

    variants, n_sub = V.build_sets(base)
    variants = [v for i, v in enumerate(variants) if i % a.nchunks == a.chunk]
    print(f"chunk {a.chunk}/{a.nchunks}: {len(variants)} variants, "
          f"{a.nrep} relax replicates each, paired wild-type control")

    def relax_in_shell(pose, allowed, nrep):
        """FastRelax with bb+chi free inside `allowed` and nothing else movable.
        Returns (min score, [all scores], max non-shell heavy-atom drift)."""
        outside = [i for i in range(1, pose.total_residue() + 1)
                   if i not in set(allowed)]
        before = heavy_coords(pose, outside)
        scores, worst = [], 0.0
        for rep in range(nrep):
            p = pose.clone()
            mm = MoveMap()
            mm.set_bb(False)
            mm.set_chi(False)
            for i in allowed:
                mm.set_bb(i, True)
                mm.set_chi(i, True)
            tf, _ = LR.restrict_packing(p, allowed)
            fr = FastRelax(sfxn, 1)          # 1 = one relax cycle, standard "fast"
            fr.cartesian(True)               # see header: torsion space leaks 0.98 A
            fr.min_type("lbfgs_armijo_nonmonotone")
            fr.set_movemap(mm)               # minimisation only -- see header
            fr.set_task_factory(tf)          # this is what restricts the PACKER
            fr.apply(p)
            scores.append(float(sfxn(p)))
            after = heavy_coords(p, outside)
            for i in outside:
                if i in before and i in after and len(before[i]) == len(after[i]):
                    d = float(np.abs(before[i] - after[i]).max())
                    worst = max(worst, d)
        return min(scores), scores, worst

    rows = []
    for v in variants:
        pose = base.clone()
        idx, ok = [], True
        for num, wt, mu in v["subs"]:
            one, i = V.wt_at(pose, num)
            if one != wt:
                print(f"  !! {v['id']}: residue {num} is {one}, sheet says {wt}")
                ok = False
                break
            MutateResidue(i, THREE[mu]).apply(pose)
            idx.append(i)
        if not ok:
            continue
        sel = ResidueIndexSelector(",".join(str(i) for i in idx))
        nb = NeighborhoodResidueSelector(sel, 8.0, True)
        allowed = [i + 1 for i, b in enumerate(nb.apply(pose)) if b]

        s_mut, all_mut, drift_m = relax_in_shell(pose, allowed, a.nrep)
        s_wt, all_wt, drift_w = relax_in_shell(base.clone(), allowed, a.nrep)

        rows.append({
            "set": v["set"], "id": v["id"], "n_sub": len(v["subs"]),
            "n_shell": len(allowed),
            "score": s_mut, "wt_same_shell": s_wt, "ddG": s_mut - s_wt,
            "mut_reps": all_mut, "wt_reps": all_wt,
            "mut_spread": float(max(all_mut) - min(all_mut)),
            "wt_spread": float(max(all_wt) - min(all_wt)),
            "nonshell_drift": float(max(drift_m, drift_w)),
            "subs": [f"{w}{n}{m}" for n, w, m in v["subs"]]})
        print(f"  {v['set']:<8}{v['id']:<14}{len(v['subs'])} subs shell "
              f"{len(allowed):3d}  mut {s_mut:9.2f} (spread {max(all_mut)-min(all_mut):5.1f})"
              f"  wt {s_wt:9.2f}  ddG {s_mut-s_wt:+8.2f}  drift {max(drift_m,drift_w):.2f}")
    with open(os.path.join(OUT, f"chunk_{a.chunk:03d}.json"), "w") as fh:
        json.dump({"nrep": a.nrep, "rows": rows}, fh, indent=1)
    print(f"wrote {OUT}/chunk_{a.chunk:03d}.json  ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
