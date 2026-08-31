#!/usr/bin/env python
r"""
160_conditional_park.py -- score Park's library in the background the SCREEN
would have reached, with a fixed shell.

TWO DEFECTS IN 156, BOTH VISIBLE IN ITS OWN OUTPUT.

1. PER-VARIANT SHELLS. 156 repacked a 6 A neighbourhood of the mutated position,
   so each of the 475 was scored in a different shell and cross-position ranking
   was not meaningful. Worse, K59R sits outside almost every shell, so it is
   frozen and its energy CANCELS in the paired difference: K59R changed ddG for
   22 of 475 variants, and only at positions 55, 108, 141, 158, 164 and 167. The
   A/B Jannis asked for was therefore mostly vacuous. Fixed here by a SINGLE
   shell -- the union of the 6 A neighbourhoods of all 25 Park positions, plus
   the ligand -- identical for every variant and every background.

2. ONE CLASH SWAMPS EVERYTHING. The mandipropamid pose comes from the quadruple
   mutant crystal, so wild-type F108 clashes with it enormously. In 156 the top
   12 of 475 were F108S/A/C/G/P/T/D/Q/N/E/V/M at -1358 to -1484 REU, while V81I
   ranked 260 and F159L 146. That is not the score failing to find them: it is
   the score correctly reporting that nothing else matters until F108 is
   relieved.

So the question becomes conditional, and it is the one a screen actually faces:
GIVEN the dominant clash relieved, can the score find the rest? Three
backgrounds, all 456 non-F108 substitutions scored in each:

    K59R              the library as built
    K59R + F108A      the state after round 1 selects the obvious winner
    K59R + F108A + F159L   after two, leaving V81I as the only target

If V81I and F159L rise in the conditional backgrounds, the search is sequential
rather than impossible, and Goal 3 needs a greedy protocol rather than a better
score. If they do not, the score has no information about them at all.

⚠ Cross-reactivity is not modelled (Jannis), so a high rank is not proof of
selection. And the volume correction of README 88d is NOT applied here: it
failed to transfer to this dataset (see 88f) and applying it would be fitting
the answer.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import lib_rosetta as LR                                          # noqa: E402

OUT = os.path.join(ROOT, "results", "park_conditional")
FRAME = os.path.join(ROOT, "data", "stage1", "wt_mandi.pdb")
LIGPDB = os.path.join(ROOT, "data", "stage1", "params", "3UZ_0001.pdb")
LIGPRM = os.path.join(ROOT, "data", "stage1", "params", "3UZ.params")
PARK25 = {55: "P", 61: "F", 62: "I", 81: "V", 83: "V", 87: "L", 88: "P", 89: "A",
          92: "S", 94: "E", 141: "E", 108: "F", 110: "I", 115: "H", 116: "R",
          117: "L", 120: "Y", 122: "S", 158: "M", 159: "F", 160: "A", 162: "T",
          163: "V", 164: "V", 167: "N"}
BACKGROUNDS = {"K59R": [(59, "ARG")],
               "K59R_F108A": [(59, "ARG"), (108, "ALA")],
               "K59R_F108A_F159L": [(59, "ARG"), (108, "ALA"), (159, "LEU")]}
THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
         "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
         "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
         "W": "TRP", "Y": "TYR"}
AA20 = "ACDEFGHIKLMNPQRSTVWY"


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
    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ex1 -ex2aro "
                   f"-extra_res_fa {LIGPRM}")
    sfxn = pyrosetta.create_score_function("ref2015_cart")
    base = pyrosetta.pose_from_pdb(FRAME)
    info = base.pdb_info()
    idx = {info.number(i): i for i in range(1, base.total_residue() + 1)}
    lig = base.total_residue()
    bad = {n: base.residue(idx[n]).name1() for n, w in PARK25.items()
           if idx.get(n) is None or base.residue(idx[n]).name1() != w}
    assert not bad, f"Park identities disagree: {bad}"
    assert base.residue(idx[59]).name1() == "K"

    # ---- ONE shell for everything: union over all 25 positions ---------
    sel = ResidueIndexSelector(",".join(str(idx[p]) for p in sorted(PARK25)))
    nb = NeighborhoodResidueSelector(sel, 6.0, True)
    SHELL = sorted({i + 1 for i, b in enumerate(nb.apply(base)) if b} | {lig})
    print(f"fixed shell: {len(SHELL)} residues (union of 25 pocket "
          f"neighbourhoods + ligand); identical for every variant", flush=True)

    e_lig = float(sfxn(pyrosetta.pose_from_pdb(LIGPDB)))

    def dg(p):
        e_c = float(sfxn(p))
        q = p.clone()
        q.delete_residue_slow(q.total_residue())
        return e_c - float(sfxn(q)) - e_lig

    def repack(p0, nrep):
        v = []
        for _ in range(nrep):
            p = p0.clone()
            tf, _ = LR.restrict_packing(p, SHELL)
            pk = PackRotamersMover(sfxn)
            pk.task_factory(tf)
            pk.apply(p)
            v.append(dg(p))
        return min(v), float(max(v) - min(v))

    bg_pose, bg_score = {}, {}
    for name, muts in BACKGROUNDS.items():
        p = base.clone()
        for num, three in muts:
            MutateResidue(idx[num], three).apply(p)
        bg_pose[name] = p
        bg_score[name], sp = repack(p.clone(), a.nrep)
        print(f"  background {name:<18} dG_bind {bg_score[name]:+10.2f} "
              f"(spread {sp:.2f})", flush=True)

    fixed = {n for _, ms in BACKGROUNDS.items() for n, _ in ms}
    jobs = [(n, mu) for n in sorted(PARK25) if n not in fixed
            for mu in AA20 if mu != PARK25[n]]
    jobs = [j for i, j in enumerate(jobs) if i % a.nchunks == a.chunk]
    print(f"chunk {a.chunk}/{a.nchunks}: {len(jobs)} variants x "
          f"{len(BACKGROUNDS)} backgrounds", flush=True)

    rows = []
    for num, mu in jobs:
        rec = {"sub": f"{PARK25[num]}{num}{mu}", "pos": num, "mut": mu}
        for name, p0 in bg_pose.items():
            p = p0.clone()
            MutateResidue(idx[num], THREE[mu]).apply(p)
            s, sp = repack(p, a.nrep)
            rec[name] = {"dG": s, "ddG": s - bg_score[name], "spread": sp}
        rows.append(rec)
        print("  " + rec["sub"].ljust(8) + "  ".join(
            f"{k} {rec[k]['ddG']:+9.2f}" for k in BACKGROUNDS), flush=True)
        json.dump({"backgrounds": bg_score, "shell": len(SHELL), "rows": rows},
                  open(os.path.join(OUT, f"chunk_{a.chunk:03d}.json"), "w"))
    print(f"wrote {len(rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
