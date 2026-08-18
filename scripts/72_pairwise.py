#!/usr/bin/env python
"""
72_pairwise.py -- enumerate mutation PAIRS against a jointly re-evaluated pocket.

WHY PAIRS AND NOT MENUS
-----------------------
README 31 killed the per-position step: steric admissibility admits 14.7 of 20
residues, ABA and mandipropamid menus overlap at Jaccard 0.84, and every miss was a
GROW mutation (V83W, V164W, V164H, V83F, A89V, E141F). The reason is structural --
the real sensors are compensating shrink/grow PAIRS (F159A+A160I in PYR1-WIN,
Y120G+A160G in 27 and 23 of Beltran's 45 sensors), and a per-position filter judges
each mutation in a context where its partner has not happened yet.

So evaluate both mutations simultaneously, and ask the question the empirical data
actually poses: does the PAIR beat the sum of its parts?

THE SIGNAL THIS EXPLOITS, MEASURED FIRST (README 31f)
-----------------------------------------------------
    ABA   in WT: max ligand-protein overlap 0.22 A, 0 of 19 atoms clashing
    mandi in WT: max 2.78 A, 7 of 29 atoms clashing, worst at PHE108 and PHE159

The cognate ligand fits its own pocket; the non-cognate one does not. That makes
ABA a built-in NEGATIVE CONTROL -- no mutation should improve it, and a method that
claims large gains for ABA is measuring noise. And the clash already names F108 and
F159 unprompted, which are 2 of the 4 mandipropamid ground-truth mutations. The
open question is the other 2: K59R and V81I do NOT clash (2.86 and 3.28 A away), so
they cannot be reached by clash relief and are the half stage 1 also missed.

THE OBJECTIVE, AND WHY IT IS DECOMPOSABLE
------------------------------------------
    lig_ov = SUM over (ligand atom, protein atom) pairs of max(0, r1 + r2 - d)

Summed rather than maxed, deliberately: a sum is ADDITIVE over protein atoms, so
the ligand overlap of a double mutant splits exactly into

    lig_ov(fixed environment) + lig_ov(rotamer at i) + lig_ov(rotamer at j)

The two mutated side chains then interact with each other and with nothing else
that is not already counted. That turns 26x20 x 26x20 into a precompute plus a tiny
per-pair search, with no approximation.

Packing is tracked alongside, because relief bought by hollowing out the pocket is
not a design -- that is exactly what the shrink half of a shrink/grow pair does on
its own.

⚠ AND PACKING MUST NOT COUNT CLASHES. The first version counted every protein heavy
atom within 4.5 A of the ligand. That rewards stuffing: `83W+163Y` came out top of
the Pareto front with 335 "contacts" and a relief of **-71.7 A**, i.e. it made the
ligand overlap catastrophically worse and was scored as the best-packed pair on the
board. A contact is only favourable if the atoms are TOUCHING and NOT overlapping,
so `contacts` now counts pairs whose van der Waals overlap lies in
[-CONTACT_FAR, +CONTACT_NEAR] -- in the attractive well, not inside each other.

STAGE A ONLY -- THIS IS A PREFILTER
------------------------------------
No chi minimisation here; rotamers are used as placed. README 31 showed that makes
the test too strict, so the tolerance is set to the loosest end of the measured
range rather than the tightest: TOL = 0.75 A is the smallest value at which the
wild-type positive control passed 25 of 26 positions. A prefilter should err toward
recall. Minimisation belongs in stage B, on the survivors.

Run with the tier1_analysis env python (PyRosetta).
"""
import itertools
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
STAGE1 = os.path.join(ROOT, "data", "stage1")
OUT = os.path.join(ROOT, "data", "pairwise")
os.makedirs(OUT, exist_ok=True)

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "F": 1.47, "CL": 1.75,
       "BR": 1.85, "I": 1.98, "P": 1.80}
AA20 = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
        "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
BB = ("N", "CA", "C", "O")
#: Loosest end of the measured control range: 0.75 A is the smallest tolerance at
#: which the WT positive control passes 25/26 (README 31e sweep). A PREFILTER should
#: favour recall; the tight value belongs in stage B.
TOL = 0.75
#: A favourable contact is one where the atoms touch without interpenetrating.
#: Counting everything within a flat 4.5 A instead makes a clash look like packing.
CONTACT_FAR = 1.0          # up to 1.0 A short of touching still counts as contact
CONTACT_NEAR = 0.25        # more overlap than this is a clash, not a contact
N_ROT = 6                  # rotamers kept per (position, residue) in the pair search
POCKET_CUTOFF = 5.0
COMPLEXES = {
    "ABA":   (os.path.join(STAGE1, "wt_aba.pdb"),
              os.path.join(STAGE1, "params", "A8S_anion.params")),
    "mandi": (os.path.join(STAGE1, "wt_mandi.pdb"),
              os.path.join(STAGE1, "params", "3UZ.params")),
}


def log(m):
    print(m, flush=True)


def good_contacts(xyz, rad, exyz, erad):
    """pairs in van der Waals CONTACT: touching, but not interpenetrating"""
    if len(xyz) == 0 or len(exyz) == 0:
        return 0
    d = np.linalg.norm(xyz[:, None, :] - exyz[None, :, :], axis=2)
    ov = (rad[:, None] + erad[None, :]) - d
    return int(((ov >= -CONTACT_FAR) & (ov <= CONTACT_NEAR)).sum())


def sum_overlap(xyz, rad, exyz, erad):
    if len(xyz) == 0 or len(exyz) == 0:
        return 0.0
    d = np.linalg.norm(xyz[:, None, :] - exyz[None, :, :], axis=2)
    return float(np.clip((rad[:, None] + erad[None, :]) - d, 0, None).sum())


def max_overlap(xyz, rad, exyz, erad):
    if len(xyz) == 0 or len(exyz) == 0:
        return 0.0
    d = np.linalg.norm(xyz[:, None, :] - exyz[None, :, :], axis=2)
    return float(((rad[:, None] + erad[None, :]) - d).max())


def analyse(tag, pdb, params):
    import pyrosetta
    from pyrosetta.rosetta.core.chemical import ChemicalManager
    from pyrosetta.rosetta.core.pack.rotamer_set import bb_independent_rotamers
    pyrosetta.init(f"-mute all -extra_res_fa {params} -detect_disulf false "
                   f"-run:constant_seed -run:jran 20260818")
    pose = pyrosetta.pose_from_pdb(pdb)
    pi = pose.pdb_info()
    rts = ChemicalManager.get_instance().residue_type_set("fa_standard")

    # ---- decompose the structure ----
    prot, lig = {}, None
    for i in range(1, pose.total_residue() + 1):
        r = pose.residue(i)
        bbx, bbr, scx, scr = [], [], [], []
        for a in range(1, r.natoms() + 1):
            e = r.atom_type(a).element().strip().upper()
            if e == "H":
                continue
            v = r.xyz(a)
            nm = r.atom_name(a).strip()
            (bbx if nm in BB else scx).append([v.x, v.y, v.z])
            (bbr if nm in BB else scr).append(VDW.get(e, 1.7))
        d = dict(idx=i, bbx=np.array(bbx).reshape(-1, 3), bbr=np.array(bbr),
                 scx=np.array(scx).reshape(-1, 3), scr=np.array(scr),
                 name3=r.name3().strip())
        if r.is_protein():
            prot[pi.number(i)] = d
        else:
            lig = d
    LX = np.vstack([lig["bbx"], lig["scx"]]) if len(lig["bbx"]) else lig["scx"]
    LR = np.concatenate([lig["bbr"], lig["scr"]]) if len(lig["bbr"]) else lig["scr"]

    pocket = sorted(n for n, d in prot.items()
                    if len(d["scx"]) and np.linalg.norm(
                        d["scx"][:, None, :] - LX[None, :, :], axis=2).min() < POCKET_CUTOFF
                    or np.linalg.norm(d["bbx"][:, None, :] - LX[None, :, :],
                                      axis=2).min() < POCKET_CUTOFF)
    log(f"\n=== {tag}: {len(LX)} ligand heavy atoms, {len(pocket)} pocket positions ===")

    # WT baseline
    allx = np.vstack([np.vstack([d["bbx"], d["scx"]]) if len(d["scx"]) else d["bbx"]
                      for d in prot.values()])
    allr = np.concatenate([np.concatenate([d["bbr"], d["scr"]]) if len(d["scr"])
                           else d["bbr"] for d in prot.values()])
    wt_ov = sum_overlap(LX, LR, allx, allr)
    wt_ct = good_contacts(LX, LR, allx, allr)
    log(f"  WT: lig_ov {wt_ov:.2f} A, contacts {wt_ct}")

    # ---- precompute rotamers ----
    # For each (position, aa, rotamer): its ligand-overlap contribution, its clash
    # against everything that is NOT a pocket side chain, and its clash against each
    # pocket side chain separately -- so a pair can drop exactly the two side chains
    # it replaces without recomputing anything.
    nonpocket_x, nonpocket_r = [], []
    for n, d in prot.items():
        nonpocket_x.append(d["bbx"]); nonpocket_r.append(d["bbr"])
        if n not in pocket and len(d["scx"]):
            nonpocket_x.append(d["scx"]); nonpocket_r.append(d["scr"])
    NPX, NPR = np.vstack(nonpocket_x), np.concatenate(nonpocket_r)

    rot = defaultdict(dict)
    for n in pocket:
        target = pose.residue(prot[n]["idx"])
        # own backbone and immediate neighbours' backbone are bonded context
        own = {prot[n]["idx"], prot[n]["idx"] - 1, prot[n]["idx"] + 1}
        keep = np.array([True] * len(NPX))
        off = 0
        for m, d in prot.items():
            k = len(d["bbx"])
            if d["idx"] in own:
                keep[off:off + k] = False
            off += k
            if m not in pocket and len(d["scx"]):
                off += len(d["scx"])
        EX, ER = NPX[keep], NPR[keep]
        for aa in AA20:
            cands = []
            try:
                rots = bb_independent_rotamers(rts.name_map(aa))
            except Exception:
                continue
            for k in range(1, len(rots) + 1):
                r = rots[k].clone()
                r.orient_onto_residue(target)
                xs, rs = [], []
                for a in range(1, r.natoms() + 1):
                    e = r.atom_type(a).element().strip().upper()
                    nm = r.atom_name(a).strip()
                    if e == "H" or nm in BB:
                        continue
                    v = r.xyz(a)
                    xs.append([v.x, v.y, v.z]); rs.append(VDW.get(e, 1.7))
                X, R = np.array(xs).reshape(-1, 3), np.array(rs)
                if max_overlap(X, R, EX, ER) > TOL:
                    continue                      # infeasible against fixed context
                vs = {m: max_overlap(X, R, prot[m]["scx"], prot[m]["scr"])
                      for m in pocket if m != n and len(prot[m]["scx"])}
                cands.append(dict(X=X, R=R,
                                  lig=sum_overlap(X, R, LX, LR),
                                  ct=good_contacts(LX, LR, X, R),
                                  vs=vs))
            cands.sort(key=lambda c: c["lig"])
            rot[n][aa] = cands
    nrot = sum(len(v) for d in rot.values() for v in d.values())
    log(f"  precomputed {nrot} feasible rotamers over {len(pocket)}x{len(AA20)} slots")
    return dict(pose=pose, prot=prot, pocket=pocket, LX=LX, LR=LR, rot=rot,
                wt_ov=wt_ov, wt_ct=wt_ct, NPX=NPX, NPR=NPR)


def baseline(ctx, i, j):
    """ligand overlap and contacts from everything EXCEPT positions i and j"""
    prot, LX, LR = ctx["prot"], ctx["LX"], ctx["LR"]
    xs, rs = [], []
    for n, d in prot.items():
        xs.append(d["bbx"]); rs.append(d["bbr"])
        if n not in (i, j) and len(d["scx"]):
            xs.append(d["scx"]); rs.append(d["scr"])
    X, R = np.vstack(xs), np.concatenate(rs)
    return sum_overlap(LX, LR, X, R), good_contacts(LX, LR, X, R)


def best_single(ctx, n, aa):
    c = ctx["rot"][n].get(aa)
    return c[0] if c else None


def main():
    ctxs = {t: analyse(t, p, q) for t, (p, q) in COMPLEXES.items()}

    results = {}
    for tag, ctx in ctxs.items():
        pocket, rotset = ctx["pocket"], ctx["rot"]
        # ---- singles ----
        singles = {}
        for n in pocket:
            b0, c0 = baseline(ctx, n, None)
            for aa in AA20:
                c = best_single(ctx, n, aa)
                if c is None:
                    continue
                singles[(n, aa)] = dict(ov=b0 + c["lig"], ct=c0 + c["ct"])
        # ---- pairs ----
        # Only position pairs whose side chains could plausibly interact: a shrink
        # at i can only make room for a grow at j if they share volume. CB-CB 12 A
        # is generous for the longest side chains (Arg/Trp reach ~7 A from CB).
        prot = ctx["prot"]
        def cb(n):
            d = prot[n]
            return d["scx"][0] if len(d["scx"]) else d["bbx"][1]
        cand_pairs = [(i, j) for i, j in itertools.combinations(pocket, 2)
                      if np.linalg.norm(cb(i) - cb(j)) <= 12.0]
        log(f"  {tag}: {len(cand_pairs)} of {len(pocket)*(len(pocket)-1)//2} "
            f"position pairs within 12 A CB-CB")

        pairs = {}
        for i, j in cand_pairs:
            b0, c0 = baseline(ctx, i, j)
            # feasibility against the pocket side chains that are NOT being replaced
            def usable(n, aa, other):
                out = []
                for c in rotset[n].get(aa, [])[:N_ROT]:
                    vs = [v for m, v in c["vs"].items() if m != other]
                    if not vs or max(vs) <= TOL:
                        out.append(c)
                return out
            Ai = {aa: usable(i, aa, j) for aa in AA20}
            Aj = {aa: usable(j, aa, i) for aa in AA20}
            # one distance matrix for the whole position pair, not 400 small ones
            fi = [(aa, k, c) for aa in AA20 for k, c in enumerate(Ai[aa])]
            fj = [(aa, k, c) for aa in AA20 for k, c in enumerate(Aj[aa])]
            if not fi or not fj:
                continue
            Xi = np.vstack([c["X"] for _, _, c in fi]) if any(len(c["X"]) for _,_,c in fi) else None
            starts_i, off = [], 0
            for _, _, c in fi:
                starts_i.append(off); off += len(c["X"])
            starts_j, off = [], 0
            for _, _, c in fj:
                starts_j.append(off); off += len(c["X"])
            Xj = np.vstack([c["X"] for _, _, c in fj]) if any(len(c["X"]) for _,_,c in fj) else None
            if Xi is None or Xj is None or len(Xi) == 0 or len(Xj) == 0:
                OVmax = np.zeros((len(fi), len(fj)))
            else:
                Ri = np.concatenate([c["R"] for _, _, c in fi])
                Rj = np.concatenate([c["R"] for _, _, c in fj])
                D = np.linalg.norm(Xi[:, None, :] - Xj[None, :, :], axis=2)
                OV = (Ri[:, None] + Rj[None, :]) - D
                # collapse atoms -> rotamers with grouped maxima (GLY has no atoms,
                # so its group is empty and must default to "no clash", not -inf)
                red = np.full((len(fi), len(fj)), -np.inf)
                si = [k for k in starts_i] + [len(Xi)]
                sj = [k for k in starts_j] + [len(Xj)]
                for a in range(len(fi)):
                    if si[a] == si[a + 1]:
                        red[a, :] = -np.inf
                        continue
                    blk = OV[si[a]:si[a + 1], :]
                    for b in range(len(fj)):
                        if sj[b] == sj[b + 1]:
                            continue
                        red[a, b] = blk[:, sj[b]:sj[b + 1]].max()
                OVmax = np.where(np.isfinite(red), red, -np.inf)
            for a, (ai, _, ca) in enumerate(fi):
                for b, (aj, _, cb_) in enumerate(fj):
                    if OVmax[a, b] > TOL:
                        continue
                    key = (i, ai, j, aj)
                    ov = b0 + ca["lig"] + cb_["lig"]
                    prev = pairs.get(key)
                    if prev is None or ov < prev["ov"]:
                        pairs[key] = dict(ov=ov, ct=c0 + ca["ct"] + cb_["ct"])

        results[tag] = dict(singles=singles, pairs=pairs,
                            wt_ov=ctx["wt_ov"], wt_ct=ctx["wt_ct"])
        log(f"  {tag}: {len(singles)} singles, {len(pairs)} feasible pairs")

    json.dump({t: dict(wt_ov=r["wt_ov"], wt_ct=r["wt_ct"],
                       singles={f"{k[0]}{k[1]}": v for k, v in r["singles"].items()},
                       pairs={f"{k[0]}{k[1]}|{k[2]}{k[3]}": v
                              for k, v in r["pairs"].items()})
               for t, r in results.items()},
              open(os.path.join(OUT, "pairwise.json"), "w"), indent=0)
    log(f"\nwrote {OUT}/pairwise.json")


if __name__ == "__main__":
    main()
