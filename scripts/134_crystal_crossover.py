#!/usr/bin/env python
r"""
134_crystal_crossover.py -- score the CROSS-OVER on crystal poses, with and
without relaxation.

THE QUESTION (Jannis, 2026-08-26). §77 found ligand-aware scoring at chance, with
two demonstrated causes it could not apportion: docking noise, and §76's scoring
defect. This removes BOTH sources of doubt at once -- use the crystal pose, and
score without relaxing -- and asks whether the score can do the one thing it must:

    PYR1^MANDI + mandipropamid   must beat   WT PYR1 + mandipropamid
    WT PYR1    + ABA             must beat   PYR1^MANDI + ABA

A receptor x ligand cross-over, so the DIAGONAL has to win in both directions.
Getting one right is easy for the wrong reason (a mutant pocket is emptier, an
unmutated one is tighter); getting both right means the score is tracking
complementarity rather than a size artefact.

THREE PROTOCOLS, because Jannis's point is that relax may be the problem:
    raw     score the coordinates as they are
    repack  side chains only, backbone fixed
    relax   Cartesian FastRelax of the pocket shell + ligand    (what §77 used)
If the answer is right at `raw` and wrong at `relax`, relaxation is destroying it.
If it is wrong at `raw` too, the score itself cannot see the difference.

⚠ FAIRNESS. All four cells are built the SAME way in the SAME stage-1 frame, so
none is a crystal while another is a model. The quad cells carry 4WVO's own side
chains at 59/81/108/159 (transplanted, §101/§133); the wild-type cells carry the
frame's own. No cell gets a packing advantage from provenance.

⚠ NO WIN CELL. §38 pre-registered WIN 55,212-2 as a held-out test, but there is no
PYR1-WIN structure in `data/`, so "give it the right pose" is not available for
WIN. Running it would mean docking, which is the noise §77 already measured.

Usage: 134_crystal_crossover.py --ligand {mandi,aba}
  (the two params files both declare NAME LIG, so they cannot co-load; the
   comparisons are within-ligand anyway)
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
    "r133", os.path.join(HERE, "133_relax_vs_crystal.py"))
R = importlib.util.module_from_spec(_s)
_s.loader.exec_module(R)

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "crossover")
CFG = {
    "mandi": (os.path.join(ROOT, "data", "stage1", "wt_mandi.pdb"),
              os.path.join(ROOT, "data", "stage1", "params", "3UZ_0001.pdb"),
              os.path.join(ROOT, "data", "stage1", "params", "3UZ.params")),
    "aba":   (os.path.join(ROOT, "data", "stage1", "wt_aba.pdb"),
              os.path.join(ROOT, "data", "stage1", "params", "A8S_anion_0001.pdb"),
              os.path.join(ROOT, "data", "stage1", "params", "A8S_anion.params")),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ligand", choices=("mandi", "aba"), required=True)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    frame_pdb, ligpdb, ligparams = CFG[a.ligand]
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    frame = [l for l in open(frame_pdb) if l.startswith("ATOM")]
    ligl = [l for l in open(ligpdb) if l.startswith("HETATM")]
    het = [l for l in open(frame_pdb) if l.startswith("HETATM")]
    d0 = np.linalg.norm(
        np.array([float(ligl[0][30:38]), float(ligl[0][38:46]), float(ligl[0][46:54])])
        - np.array([float(het[0][30:38]), float(het[0][38:46]), float(het[0][46:54])]))
    say("=" * 76)
    say(f"CROSS-OVER on crystal poses -- ligand = {a.ligand}")
    say("=" * 76)
    say(f"   frame {os.path.basename(frame_pdb)}, params ligand matches its HETATM "
        f"to {d0:.4f} A")
    if d0 > 0.01:
        raise SystemExit("params ligand is not the frame's crystal pose")

    xt, xseq = R.read_cif(R.CIF)
    fbb = {}
    for l in frame:
        nm = l[12:16].strip()
        if nm in ("N", "CA", "C", "O"):
            fbb[(int(l[22:26]), nm)] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    com = sorted(set(fbb) & set(xt))
    P = np.array([xt[k] for k in com])
    Q = np.array([fbb[k] for k in com])
    Rm, mP, mQ = R.kabsch(P, Q)
    say(f"   4WVO -> frame on {len(com)} backbone atoms, RMSD "
        f"{float(np.sqrt((((P-mP)@Rm+mQ-Q)**2).sum(1).mean())):.2f} A")

    def build(quad):
        out = []
        for l in frame:
            num = int(l[22:26])
            nm = l[12:16].strip()
            if quad and num in R.MUT:
                if nm not in R.BB:
                    continue
                out.append(l[:17] + f"{R.MUT[num]:>3s}" + l[20:])
            else:
                out.append(l)
        if quad:
            for num in sorted(R.MUT):
                tmpl = next(l for l in frame if int(l[22:26]) == num
                            and l[12:16].strip() not in R.BB)
                for at, v in sorted({k[1]: (vv - mP) @ Rm + mQ
                                     for k, vv in xt.items()
                                     if k[0] == num and k[1] not in R.BB}.items()):
                    nm4 = f" {at:<3s}" if len(at) < 4 else at
                    out.append(tmpl[:12] + nm4 + tmpl[16:17] + f"{R.MUT[num]:>3s}"
                               + tmpl[20:30]
                               + f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}" + tmpl[54:])
        out.sort(key=lambda l: (int(l[22:26]), l[12:16]))
        return out

    import pyrosetta
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ex1 -ex2aro "
                   f"-extra_res_fa {ligparams}")
    sfxn = pyrosetta.create_score_function("ref2015_cart")

    res = {}
    for quad in (False, True):
        tag = "PYR1^MANDI" if quad else "WT PYR1"
        f = os.path.join(OUT, f"{a.ligand}_{'quad' if quad else 'wt'}.pdb")
        with open(f, "w") as fh:
            fh.writelines(build(quad))
            fh.write("TER\n")
            fh.writelines(ligl)
            fh.write("END\n")
        pose = pyrosetta.pose_from_pdb(f)
        info = pose.pdb_info()
        idx = {info.number(i): i for i in range(1, pose.total_residue() + 1)}
        lig = pose.total_residue()
        for num, want in R.MUT.items():
            got = pose.residue(idx[num]).name3()
            exp = want if quad else None
            if quad and got != want:
                raise SystemExit(f"{tag}: position {num} is {got}, expected {want}")
        sel = ResidueIndexSelector(",".join(str(idx[p]) for p in R.POCKET if p in idx))
        nb = NeighborhoodResidueSelector(sel, 6.0, True)
        allowed = sorted({i + 1 for i, b in enumerate(nb.apply(pose)) if b} | {lig})

        def split_score(p):
            """dG_bind = E(complex) - E(protein) - E(ligand), same protocol each."""
            e_c = float(sfxn(p))
            pp = p.clone()
            pp.delete_residue_slow(pp.total_residue())
            e_p = float(sfxn(pp))
            lp = pyrosetta.pose_from_pdb(ligpdb)
            e_l = float(sfxn(lp))
            return e_c - e_p - e_l, e_c

        out = {}
        out["raw"] = split_score(pose.clone())
        pr = pose.clone()
        tf, _ = LR.restrict_packing(pr, allowed)
        pk = PackRotamersMover(sfxn)
        pk.task_factory(tf)
        pk.apply(pr)
        out["repack"] = split_score(pr)
        px = pose.clone()
        mm = MoveMap()
        mm.set_bb(False)
        mm.set_chi(False)
        for i in allowed:
            mm.set_bb(i, True)
            mm.set_chi(i, True)
        tf2, _ = LR.restrict_packing(px, allowed)
        fr = FastRelax(sfxn, 1)
        fr.cartesian(True)
        fr.min_type("lbfgs_armijo_nonmonotone")
        fr.set_movemap(mm)
        fr.set_task_factory(tf2)
        fr.apply(px)
        out["relax"] = split_score(px)
        res[tag] = out
        say(f"   {tag:<12} dG_bind  raw {out['raw'][0]:9.2f}   "
            f"repack {out['repack'][0]:9.2f}   relax {out['relax'][0]:9.2f}")

    say("")
    say(f"   SELECTIVITY for {a.ligand}:  dG_bind(PYR1^MANDI) - dG_bind(WT)")
    for proto in ("raw", "repack", "relax"):
        d = res["PYR1^MANDI"][proto][0] - res["WT PYR1"][proto][0]
        want = "quad should WIN (negative)" if a.ligand == "mandi" \
            else "WT should win (positive)"
        ok = (d < 0) if a.ligand == "mandi" else (d > 0)
        say(f"     {proto:<8}{d:+9.2f} REU   {want:<28}{'CORRECT' if ok else 'WRONG'}")
    with open(os.path.join(OUT, f"{a.ligand}.json"), "w") as fh:
        json.dump({k: {p: list(v) for p, v in d.items()} for k, d in res.items()},
                  fh, indent=1)
    with open(os.path.join(OUT, f"{a.ligand}.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
