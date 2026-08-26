#!/usr/bin/env python
r"""
125_cavity_match.py -- test §68e: is CAVITY MATCH a better objective than ddG?

THE ARGUMENT BEING TESTED (README §68)
§67 found that Cartesian FastRelax ddG separates real coumarin sensors from random
members of the same library (AUC 0.284). §68 then found the mechanism and a
problem: ref2015 penalises cavities, so it prefers GROW substitutions (mean beta
+0.43) over SHRINK (+1.75). Coumarins are 11-15 heavy atoms against ABA's 19, so
their sensors grow the lining -- the same direction ref2015 wants. Across 125 sd03
ligands, Spearman(ligand size, mean dVolume used) = -0.30, and ligands of >= 29
heavy atoms average -15.1 A^3, i.e. they SHRINK the lining.

So ddG is predicted to have worked for a reason that will not generalise, and the
proposed fix is to change the OBJECTIVE rather than the method:

    score by |cavity volume - target ligand volume|, not by ddG

which stays POSE-FREE, because a ligand's volume is a property of its SMILES and
needs no docked pose (§49b's circularity is not re-imported).

THE TWO-ARM DESIGN, and why PFAS is the arm that matters
    COUMARIN  target volume ~163 A^3, BELOW PYR1's ~174 A^3 cavity -- the pocket
              must close slightly. This is where ddG already works.
    PFAS      target volume ~225 A^3, WELL ABOVE it -- the pocket must OPEN.
              This is the pocket-expansion regime, and the one we actually care
              about. sd09's `mut_lib` splits it into round 1 (DSM-Hao, 89 clones)
              and round 2 (PFOS_GenWT, 154 clones); round 2 is the sensor set.

PRE-REGISTERED PREDICTIONS, written before the run
    1. ddG on COUMARIN reproduces §67:                 AUC < 0.40
    2. ddG on PFAS FAILS:                              AUC >= 0.45
    3. cavity match works on BOTH:                     AUC < 0.40 in each arm
  Outcomes and what they mean:
    - 1,2,3 all hold -> the diagnosis is right and cavity match is the fix
    - 3 fails on coumarin -> the idea is dead, ddG stays the only filter
    - 2 fails (ddG works on PFAS too) -> §68d's volume argument is wrong and ddG
      is more general than its mechanism suggested
    - 3 works only where ddG works -> cavity match adds nothing

⚠ CAVITY CAVEAT carried from lib_cavity's own docstring: truncating a residue that
SEALS the cavity lowers the reported volume, because the component leaks and
buriedness drops. A shrinking number is therefore ambiguous between "smaller
cavity" and "opened to solvent". Both arms are treated identically so the
comparison is still fair, but no absolute volume here should be quoted elsewhere.

Usage: 125_cavity_match.py --klass {coumarin,pfas} --chunk i --nchunks n
"""
import argparse
import json
import os
import random
import re
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_xlsx import table                                     # noqa: E402
import lib_rosetta as LR                                        # noqa: E402
import lib_cavity as LC                                         # noqa: E402

ROOT = os.path.dirname(HERE)
SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "cavity_match")
PDB = os.path.join(ROOT, "data", "structures_191", "pyr1_closed_191.pdb")
ABA_PDB = os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")
AA = "ACDEFGHIKLMNPQRSTVWY"
THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
         "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
         "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
         "W": "TRP", "Y": "TYR"}
N_RANDOM = 300
SEED = 0

#: measured with RDKit (MMFF-optimised conformer, ComputeMolVolume) in this file's
#: companion analysis; class MEDIAN is used because a library is designed per class,
#: not per ligand, and LIBRARY/WILD variants have no ligand of their own
TARGET_VOL = {"coumarin": 162.7, "pfas": 224.6}
LIBNAME = {"coumarin": "Coumarin", "pfas": "PFAS"}


def sd04_menu(name):
    _, recs = table(f"{SD}/pnas.2519924122.sd04.xlsx")
    m = defaultdict(set)
    for r in recs:
        if (r.get("Library") or "").strip() == name:
            m[int(float(r["Position"]))] |= set(
                r["Amino Acids allowed for mutation"].strip())
    return dict(m)


def real_sensors(klass):
    """Round-2 sensors for the class. PFAS lives in sd09 and must be filtered on
    mut_lib -- DSM-Hao rows there are ROUND 1 and are not sensors of this library."""
    if klass == "coumarin":
        path, ligcol, want = f"{SD}/pnas.2519924122.sd07(1).xlsx", "compound", None
    else:
        path, ligcol, want = f"{SD}/pnas.2519924122.sd09.xlsx", "chem_name", "PFOS_GenWT"
    hdr, recs = table(path)
    cols = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    out = []
    for r in recs:
        if want and (r.get("mut_lib") or "").strip() != want:
            continue
        lig = (r.get(ligcol) or "").strip()
        if not lig or lig.startswith("Bold:") or lig.startswith("*"):
            continue
        subs = []
        for c in cols:
            v = (r.get(c) or "").strip().upper()
            if len(v) == 1 and v in AA and v != c[0]:
                subs.append((int(c[1:]), c[0], v))
        if subs:
            out.append({"lig": lig, "subs": subs})
    return out


def wt_at(pose, num):
    info = pose.pdb_info()
    for i in range(1, pose.total_residue() + 1):
        if info.number(i) == num:
            return pose.residue(i).name1(), i
    raise SystemExit(f"native residue {num} not found in {PDB}")


def build_sets(pose, klass):
    menu = sd04_menu(LIBNAME[klass])
    dsm = sd04_menu("DSM-Hao")
    positions = sorted(menu)
    for p in positions:
        one, _ = wt_at(pose, p)
        _, recs = table(f"{SD}/pnas.2519924122.sd04.xlsx")
        want = next((r["WT"].strip() for r in recs
                     if (r.get("Library") or "").strip() == LIBNAME[klass]
                     and int(float(r["Position"])) == p), None)
        if want and one != want:
            raise SystemExit(f"position {p}: structure {one}, sd04 {want}")
    rng = random.Random(SEED)
    real = [{"set": "REAL", "id": f"R{i}", "subs": s["subs"]}
            for i, s in enumerate(real_sensors(klass))
            if all(n in menu for n, _w, _m in s["subs"])]
    if not real:
        raise SystemExit(f"no {klass} sensors fall inside the sd04 menu")
    ns = [len(r["subs"]) for r in real]
    lo, hi = min(ns), max(ns)

    def draw(m, tag):
        out = []
        for i in range(N_RANDOM):
            k = rng.randint(lo, hi)
            subs = []
            for p in rng.sample(positions, min(k, len(positions))):
                one, _ = wt_at(pose, p)
                ch = sorted(m[p] - {one})
                if ch:
                    subs.append((p, one, rng.choice(ch)))
            if subs:
                out.append({"set": tag, "id": f"{tag}_{i}", "subs": subs})
        return out

    return real + draw(menu, "LIBRARY") + draw(dsm, "WILD"), (lo, hi)


def pocket_seed(pose):
    """Cavity seed point: the ABA centroid from data/stage1/wt_aba.pdb, carried
    across by superposing shared backbone atoms. Asserted, not assumed."""
    ref, lig = {}, []
    for l in open(ABA_PDB):
        if l.startswith("ATOM") and l[12:16].strip() in ("N", "CA", "C", "O"):
            ref[(int(l[22:26]), l[12:16].strip())] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
        elif l.startswith("HETATM") and l[17:20].strip() == "A8S":
            lig.append(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    if len(lig) < 15:
        raise SystemExit(f"only {len(lig)} A8S atoms in {ABA_PDB}")
    info = pose.pdb_info()
    tgt = {}
    for i in range(1, pose.total_residue() + 1):
        r = pose.residue(i)
        for nm in ("N", "CA", "C", "O"):
            if r.has(nm):
                v = r.xyz(nm)
                tgt[(info.number(i), nm)] = np.array([v.x, v.y, v.z])
    com = sorted(set(ref) & set(tgt))
    if len(com) < 400:
        raise SystemExit(f"only {len(com)} shared backbone atoms for superposition")
    P = np.array([ref[k] for k in com])
    Q = np.array([tgt[k] for k in com])
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    D = np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))])
    R = V @ D @ Wt
    rms = float(np.sqrt((((P - P.mean(0)) @ R + Q.mean(0) - Q) ** 2).sum(1).mean()))
    cen = np.mean(lig, axis=0)
    return (cen - P.mean(0)) @ R + Q.mean(0), rms, len(com)


def pose_atoms(pose):
    xyz, el = [], []
    for i in range(1, pose.total_residue() + 1):
        r = pose.residue(i)
        for j in range(1, r.natoms() + 1):
            e = r.atom_type(j).element().strip().upper()
            if e == "H":
                continue
            v = r.xyz(j)
            xyz.append([v.x, v.y, v.z])
            el.append(e)
    return np.array(xyz), el


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--klass", choices=("coumarin", "pfas"), required=True)
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--nchunks", type=int, default=1)
    ap.add_argument("--nrep", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    import pyrosetta
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    pyrosetta.init("-mute all -ignore_unrecognized_res -ex1 -ex2aro")
    base = pyrosetta.pose_from_pdb(PDB)
    sfxn = pyrosetta.create_score_function("ref2015_cart")

    seed, rms, ncom = pocket_seed(base)
    variants, (lo, hi) = build_sets(base, a.klass)
    variants = [v for i, v in enumerate(variants) if i % a.nchunks == a.chunk]
    wt_xyz, wt_el = pose_atoms(base)
    wt_cav, _ = LC.cavity_volume(wt_xyz, wt_el, seed)
    print(f"{a.klass}: chunk {a.chunk}/{a.nchunks}, {len(variants)} variants, "
          f"{lo}-{hi} substitutions, target volume {TARGET_VOL[a.klass]} A^3")
    print(f"  pocket seed from {ncom} shared backbone atoms, superposition RMSD "
          f"{rms:.2f} A; unmutated closed PYR1 cavity = {wt_cav:.1f} A^3")

    def relax(pose, allowed, nrep):
        best, scores = None, []
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
            s = float(sfxn(p))
            scores.append(s)
            if best is None or s < best[0]:
                best = (s, p)
        return best[0], best[1], scores

    rows = []
    for v in variants:
        pose = base.clone()
        idx, ok = [], True
        for num, wt, mu in v["subs"]:
            one, i = wt_at(pose, num)
            if one != wt:
                print(f"  !! {v['id']}: residue {num} is {one}, expected {wt}")
                ok = False
                break
            MutateResidue(i, THREE[mu]).apply(pose)
            idx.append(i)
        if not ok:
            continue
        sel = ResidueIndexSelector(",".join(str(i) for i in idx))
        nb = NeighborhoodResidueSelector(sel, 8.0, True)
        allowed = [i + 1 for i, b in enumerate(nb.apply(pose)) if b]
        s_mut, p_mut, reps_m = relax(pose, allowed, a.nrep)
        s_wt, p_wt, reps_w = relax(base.clone(), allowed, a.nrep)
        cav_m, _ = LC.cavity_volume(*pose_atoms(p_mut), seed)
        cav_w, _ = LC.cavity_volume(*pose_atoms(p_wt), seed)
        rows.append({
            "set": v["set"], "id": v["id"], "n_sub": len(v["subs"]),
            "n_shell": len(allowed), "ddG": s_mut - s_wt,
            "cav": cav_m, "cav_wt": cav_w,
            "cav_err": abs(cav_m - TARGET_VOL[a.klass]),
            "mut_spread": float(max(reps_m) - min(reps_m)),
            "subs": [f"{w}{n}{m}" for n, w, m in v["subs"]]})
        print(f"  {v['set']:<8}{v['id']:<12}{len(v['subs'])} subs  "
              f"ddG {s_mut-s_wt:+8.2f}  cav {cav_m:6.1f} (wt {cav_w:6.1f})  "
              f"|cav-target| {abs(cav_m-TARGET_VOL[a.klass]):6.1f}")
    with open(os.path.join(OUT, f"{a.klass}_{a.chunk:03d}.json"), "w") as fh:
        json.dump({"klass": a.klass, "target_vol": TARGET_VOL[a.klass],
                   "wt_cavity": wt_cav, "rows": rows}, fh, indent=1)
    print(f"wrote {OUT}/{a.klass}_{a.chunk:03d}.json ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
