#!/usr/bin/env python
"""
71_admissibility.py -- HEADROOM CHECK: how much does pure steric admissibility
actually narrow the choice of residue at each pocket position?

THE IDEA BEING TESTED
---------------------
Stage 1 (README 23f-j) failed because we asked a score function to RANK residues,
and ref2015's ranking is wrong in this pocket for a quantified reason: it pays
+10.1 REU of Lazaridis-Karplus desolvation to bury K59's ammonium, so no bonus
below the freezing regime keeps Lys there.

Geometry does not rank. It EXCLUDES. So ask a different question: at each position,
which residues can physically be placed at all, with the ligand present? That is a
hard yes/no, it is ligand-conditional by construction, and it needs no null arm.

WHY THIS IS A HEADROOM CHECK AND NOT YET A METHOD
-------------------------------------------------
The stage-1 post-mortem showed the trivial clash baseline already handed us 3 of 4
ground-truth POSITIONS for free, and the coumarin benchmark (README 25) showed
pocket positions are near ligand-independent while substitutions are not. So the
failure mode here is obvious in advance: if admissibility admits 15 of 20 residues
everywhere, it is saturated and useless no matter how elegant it is.

Therefore this script measures PERMISSIVENESS FIRST and reports it before anything
is built on top:

  1. |admissible| per position, per ligand
  2. how much the ABA menu and the mandipropamid menu DIFFER -- if they are the
     same, the method is reading pocket shape, not the ligand, and fails the
     discrimination check that stage 1's null failed silently
  3. sensitivity to the clash tolerance, which is the one free parameter

WHAT ROSETTA IS AND IS NOT USED FOR
-----------------------------------
Rotamers come from Rosetta's backbone-INDEPENDENT library, used as an observational
catalogue of side-chain conformations -- geometry, not energy. The accept/reject
test is a hard-sphere overlap against explicit van der Waals radii. **ref2015 is
never called**, and neither is any electrostatic, solvation or hydrogen-bond term:
those are precisely the terms the K59 decomposition showed to be wrong here.

⚠ ONE RELAXATION STEP IS NEEDED, AND HERE IS WHY
   The first version placed idealised rotamers on the crystal backbone and tested
   them as-is. Its positive control FAILED: the wild-type residue was not
   admissible at its OWN position for 6 of 26 positions (P88 by 1.96 A, F108 1.62,
   H115 1.45, R79 1.10), and raising the tolerance until they passed also admitted
   ~15 of 20 residues everywhere -- i.e. the only tolerance that saves the control
   destroys the measurement.

   The cause is not the pocket: it is that a library rotamer has idealised chi
   angles while the crystal side chain does not, so the nearest catalogue entry
   sits a few tenths of an angstrom off. Every residue type suffers this equally
   only if every type gets the same chance to relieve it. So the best rotamer per
   (position, residue) is chi-minimised against a **repulsive-only** score function
   -- `fa_rep` alone, no attraction, no electrostatics -- with the backbone and the
   whole environment fixed. That is still a steric calculation; it just lets the
   side chain settle into the same non-ideal geometry the crystal already has.

   Tuning the tolerance until the control passed would have been the other option.
   It is the wrong one: it changes the number being measured instead of fixing the
   thing that is broken.

⚠ TWO MORE CONTROL FAILURES, DIAGNOSED RATHER THAN ABSORBED INTO THE TOLERANCE
   Relaxation alone still left 5 of 26 positions failing, and printing the actual
   worst contacts showed two unrelated causes:

   1. BONDED ATOMS COUNTED AS CLASHES. The worst offender was
      `PRO88:NV -- LEU87:C at d = 1.35 A`, an overlap of 2.05 A. That is a peptide
      BOND, not a clash: proline's ring nitrogen is covalently attached to the
      preceding carbonyl carbon. Any side chain also packs legitimately against its
      own and its neighbours' backbone. So backbone atoms of residues i-1, i and
      i+1 are excluded from the environment -- the standard treatment in a rotamer
      bump check, and the reason this was caught is that PRO made it enormous.

   2. THE NEAREST LIBRARY ROTAMER IS IN THE WRONG WELL. R79, V83, E94 and H115 make
      NO contact worse than 0.4 A in the crystal, so nothing about the pocket
      excludes them -- minimising from a single starting rotamer simply landed in a
      different chi basin. Fixed by minimising the best FIVE rotamers rather than
      the best one. The crystal conformation itself is deliberately NOT added as a
      candidate: that would make the positive control pass by construction and stop
      it being a control.

CONTROLS, so a broken placement cannot look like a result
---------------------------------------------------------
  * the WILD-TYPE residue must be admissible at its own position -- it is in the
    crystal. If it is not, the placement or the tolerance is wrong.
  * GLY must be admissible everywhere (no side chain past CB).
  * TRP should be excluded somewhere, or the test is not discriminating at all.

Run with the tier1_analysis env python (PyRosetta).
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_sterics as st            # noqa: E402  one clash rule, shared
ROOT = os.path.dirname(HERE)
STAGE1 = os.path.join(ROOT, "data", "stage1")
OUT = os.path.join(ROOT, "data", "admissibility")
os.makedirs(OUT, exist_ok=True)

# Radii and the clash rule both come from lib_sterics, so this script, 72 and 73
# cannot drift apart -- and so the hydrogen-bond exemption applies everywhere at once.
VDW = st.VDW
AA20 = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
        "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
POCKET_CUTOFF = 5.0        # a position is "in the pocket" at this range from the ligand
TOL_MAIN = 0.5             # allowed vdW overlap, A -- swept below, not tuned
TOL_SWEEP = (0.0, 0.25, 0.5, 0.75, 1.0)
N_MIN_STARTS = 5           # chi-minimisation restarts per (position, residue)

COMPLEXES = {
    "ABA":  (os.path.join(STAGE1, "wt_aba.pdb"),
             os.path.join(STAGE1, "params", "A8S_anion.params"), "LIG"),
    "mandi": (os.path.join(STAGE1, "wt_mandi.pdb"),
              os.path.join(STAGE1, "params", "3UZ.params"), "3UZ"),
}


def log(m):
    print(m, flush=True)


def elem(name, fallback="C"):
    n = name.strip()
    for e in ("CL", "BR"):
        if n.startswith(e):
            return e
    return n[0] if n and n[0] in VDW else fallback


def load(pdb, params):
    import pyrosetta
    pyrosetta.init(f"-mute all -extra_res_fa {params} -detect_disulf false "
                   f"-run:constant_seed -run:jran 20260818")
    return pyrosetta.pose_from_pdb(pdb)


def pose_arrays(pose):
    """per-residue heavy-atom coords, radii and donor/acceptor flags, plus the ligand's"""
    res, lig = {}, None
    for i in range(1, pose.total_residue() + 1):
        r = pose.residue(i)
        xyz, rad, flag = st.atom_arrays(r)
        names = [r.atom_name(a).strip() for a in range(1, r.natoms() + 1)
                 if r.atom_type(a).element().strip().upper() != "H"]
        d = dict(xyz=xyz, rad=rad, flag=flag, names=names,
                 name3=r.name3().strip(), is_protein=r.is_protein())
        if r.is_protein():
            res[pose.pdb_info().number(i)] = (i, d)
        else:
            lig = d
    return res, lig


def rotamer_objects(restype, target_res):
    """the oriented rotamer residues themselves, for the minimisation step"""
    from pyrosetta.rosetta.core.pack.rotamer_set import bb_independent_rotamers
    out = []
    try:
        rots = bb_independent_rotamers(restype)
    except Exception:
        return out
    for k in range(1, len(rots) + 1):
        r = rots[k].clone()
        r.orient_onto_residue(target_res)
        out.append(r)
    return out


def relaxed_overlap(pose, idx, rotamer_res, env, sf_rep):
    """Place one rotamer, chi-minimise it sterically, and re-measure the overlap.

    Backbone and every other residue are fixed; only this side chain's chi angles
    move, and only against fa_rep. Returns the overlap after settling.
    """
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.protocols.minimization_packing import MinMover
    work = pose.clone()
    work.replace_residue(idx, rotamer_res, True)     # True = orient onto backbone
    mm = MoveMap()
    mm.set_bb(False)
    mm.set_jump(False)
    mm.set_chi(False)
    mm.set_chi(idx, True)
    MinMover(mm, sf_rep, "lbfgs_armijo_nonmonotone", 1e-3, True).apply(work)
    x, rr, f = st.atom_arrays(work.residue(idx), side_chain_only=True)
    return st.max_overlap(x, rr, f, *env)


def main():
    import pyrosetta
    from pyrosetta.rosetta.core.chemical import ChemicalManager

    results, poses = {}, {}
    for tag, (pdb, params, ligname) in COMPLEXES.items():
        assert os.path.exists(pdb), pdb
        assert os.path.exists(params), params
        poses[tag] = load(pdb, params)

    rts = ChemicalManager.get_instance().residue_type_set("fa_standard")
    # repulsive ONLY -- this is a steric relaxation, not an energy evaluation
    from pyrosetta.rosetta.core.scoring import ScoreFunction, fa_rep
    sf_rep = ScoreFunction()
    sf_rep.set_weight(fa_rep, 1.0)

    # ---- pocket definition: union over both ligands, so neither is privileged ----
    pocket = set()
    per_lig = {}
    for tag, pose in poses.items():
        res, lig = pose_arrays(pose)
        assert lig is not None, f"{tag}: no ligand found in the pose"
        per_lig[tag] = (res, lig)
        for num, (_, d) in res.items():
            if len(d["xyz"]) and np.linalg.norm(
                    d["xyz"][:, None, :] - lig["xyz"][None, :, :], axis=2).min() < POCKET_CUTOFF:
                pocket.add(num)
    pocket = sorted(pocket)
    log("=" * 78)
    log("STERIC ADMISSIBILITY -- HEADROOM CHECK")
    log("=" * 78)
    for tag in COMPLEXES:
        log(f"  {tag:<6} ligand heavy atoms: {len(per_lig[tag][1]['xyz'])}")
    log(f"  pocket positions within {POCKET_CUTOFF} A of either ligand: {len(pocket)}")
    log(f"  {pocket}")

    # ---- the measurement ----
    adm = {tag: defaultdict(dict) for tag in COMPLEXES}
    wt_at = {}
    for tag, pose in poses.items():
        res, lig = per_lig[tag]
        for num in pocket:
            idx, d = res[num]
            wt_at[num] = d["name3"]
            target = pose.residue(idx)
            # Environment: every other protein residue, plus the ligand -- EXCEPT
            # the backbone of the sequence neighbours i-1 and i+1. Those contacts
            # are bonded or near-bonded geometry, not steric exclusion: proline's
            # ring N sits 1.35 A from the preceding carbonyl C because they are
            # covalently joined, which the first version scored as a 2.05 A clash.
            env_xyz = [lig["xyz"]]
            env_rad = [lig["rad"]]
            env_flag = [lig["flag"]]
            for onum, (oidx, od) in res.items():
                if oidx == idx or len(od["xyz"]) == 0:
                    continue
                if abs(oidx - idx) == 1:
                    keep = [k for k, n in enumerate(od["names"])
                            if n not in ("N", "CA", "C", "O", "CB")]
                    if not keep:
                        continue
                    env_xyz.append(od["xyz"][keep])
                    env_rad.append(od["rad"][keep])
                    env_flag.append(od["flag"][keep])
                else:
                    env_xyz.append(od["xyz"])
                    env_rad.append(od["rad"])
                    env_flag.append(od["flag"])
            E = (np.vstack(env_xyz), np.concatenate(env_rad), np.concatenate(env_flag))
            for aa in AA20:
                objs = rotamer_objects(rts.name_map(aa), target)
                if not objs:
                    adm[tag][num][aa] = (None, 0)
                    continue
                scored = []
                for r in objs:
                    x, rr, f = st.atom_arrays(r, side_chain_only=True)
                    scored.append((st.max_overlap(x, rr, f, *E), r))
                scored.sort(key=lambda t: t[0])
                # Minimise the best FIVE rotamers, not just the best one. Chi
                # minimisation is local, so a single start can settle into the wrong
                # basin -- that is what made R79/V83/E94/H115 fail their own control.
                # Five starts is still ~100x cheaper than minimising every rotamer.
                best = scored[0][0]
                for _, r in scored[:N_MIN_STARTS]:
                    best = min(best, relaxed_overlap(pose, idx, r, E, sf_rep))
                adm[tag][num][aa] = (best, len(objs))

    # ---- controls ----
    log("\n" + "=" * 78)
    log("CONTROLS")
    log("=" * 78)
    bad = []
    for tag in COMPLEXES:
        for num in pocket:
            wt = wt_at[num]
            if wt not in AA20:
                continue
            ov = adm[tag][num][wt][0]
            if ov is not None and ov > TOL_MAIN:
                bad.append((tag, num, wt, round(ov, 2)))
    log(f"  wild-type admissible at its own position, tol={TOL_MAIN} A: "
        f"{'PASS' if not bad else 'FAIL ' + str(bad[:6])}")
    gly = [(t, n) for t in COMPLEXES for n in pocket
           if adm[t][n]["GLY"][0] is not None and adm[t][n]["GLY"][0] > TOL_MAIN]
    log(f"  GLY admissible everywhere: {'PASS' if not gly else 'FAIL ' + str(gly[:6])}")
    trp_excl = sum(1 for t in COMPLEXES for n in pocket
                   if adm[t][n]["TRP"][0] is not None and adm[t][n]["TRP"][0] > TOL_MAIN)
    log(f"  TRP excluded somewhere: {'PASS' if trp_excl else 'FAIL -- test is not discriminating'}"
        f"  ({trp_excl} of {2*len(pocket)} position-ligand pairs)")

    # ---- headroom ----
    log("\n" + "=" * 78)
    log(f"HOW PERMISSIVE IS IT?  (|admissible| of {len(AA20)}, tol = {TOL_MAIN} A)")
    log("=" * 78)
    log(f"  {'pos':>5} {'WT':<5}{'ABA':>6}{'mandi':>7}   {'shared':>7}{'Jaccard':>9}   "
        f"admissible for mandi but NOT ABA")
    rows = []
    for num in pocket:
        sets = {}
        for tag in COMPLEXES:
            sets[tag] = {aa for aa in AA20
                         if adm[tag][num][aa][0] is not None
                         and adm[tag][num][aa][0] <= TOL_MAIN}
        a, m = sets["ABA"], sets["mandi"]
        inter, union = a & m, a | m
        j = len(inter) / len(union) if union else 1.0
        extra = sorted(m - a)
        rows.append(dict(pos=num, wt=wt_at[num], n_aba=len(a), n_mandi=len(m),
                         jaccard=j, aba=sorted(a), mandi=sorted(m),
                         mandi_only=extra, aba_only=sorted(a - m)))
        log(f"  {num:>5} {wt_at[num]:<5}{len(a):>6}{len(m):>7}   {len(inter):>7}{j:>9.2f}   "
            f"{' '.join(extra) if extra else '-'}")

    na = np.mean([r["n_aba"] for r in rows])
    nm = np.mean([r["n_mandi"] for r in rows])
    jj = np.mean([r["jaccard"] for r in rows])
    log(f"\n  mean |admissible|: ABA {na:.1f}, mandi {nm:.1f} of 20")
    log(f"  mean per-position Jaccard(ABA, mandi): {jj:.2f}")
    log(f"  -> {'SATURATED, little headroom' if na > 14 else 'narrows the choice'}; "
        f"{'ligand-BLIND' if jj > 0.9 else 'ligand-discriminating'}")

    # ---- tolerance sensitivity: the one free parameter ----
    log("\n" + "=" * 78)
    log("SENSITIVITY TO THE CLASH TOLERANCE (the only knob)")
    log("=" * 78)
    log(f"  {'tol A':>7}{'|adm| ABA':>11}{'|adm| mandi':>13}{'Jaccard':>9}{'WT ok':>8}")
    sweep = []
    for tol in TOL_SWEEP:
        sa, sm, js, ok = [], [], [], 0
        for num in pocket:
            a = {aa for aa in AA20 if adm["ABA"][num][aa][0] is not None
                 and adm["ABA"][num][aa][0] <= tol}
            m = {aa for aa in AA20 if adm["mandi"][num][aa][0] is not None
                 and adm["mandi"][num][aa][0] <= tol}
            sa.append(len(a)); sm.append(len(m))
            js.append(len(a & m) / len(a | m) if (a | m) else 1.0)
            if wt_at[num] in a:
                ok += 1
        sweep.append(dict(tol=tol, aba=float(np.mean(sa)), mandi=float(np.mean(sm)),
                          jaccard=float(np.mean(js)), wt_ok=ok, n_pos=len(pocket)))
        log(f"  {tol:>7.2f}{np.mean(sa):>11.1f}{np.mean(sm):>13.1f}"
            f"{np.mean(js):>9.2f}{ok:>5}/{len(pocket)}")

    json.dump(dict(pocket=pocket, wt=wt_at, tol_main=TOL_MAIN,
                   per_position=rows, sweep=sweep,
                   overlaps={t: {str(n): {aa: adm[t][n][aa][0] for aa in AA20}
                                 for n in pocket} for t in COMPLEXES}),
              open(os.path.join(OUT, "admissibility.json"), "w"), indent=1)
    log(f"\nwrote {OUT}/admissibility.json")


if __name__ == "__main__":
    main()
