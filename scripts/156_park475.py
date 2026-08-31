#!/usr/bin/env python
r"""
156_park475.py -- Park's ACTUAL mandipropamid library, and the K59R A/B.

WHAT JANNIS CORRECTED. Every previous statement in this project that "K59R is
invisible to every method" (§23j, §33, §36, §47b) was scoring a prediction
nobody ever had to make. From the mandipropamid paper's methods: the library was
built in the **ABA non-responsive PYR1(K59R) backbone**, because K59R had been
isolated separately in error-prone PCR screens against structurally dissimilar
agrochemicals and removes ABA sensitivity. K59R was FORCED, not selected.

So the benchmark was wrong, not only the methods. The real task is:

    given PYR1(K59R) + mandipropamid, recover V81I, F108A and F159L
    from site-saturation at 25 pocket-lining residues

and their library is stated exactly: 475 variants = **25 positions x 19
substitutions**, at P55 F61 I62 V81 V83 L87 P88 A89 S92 E94 E141 F108 I110 H115
R116 L117 Y120 S122 M158 F159 A160 T162 V163 V164 N167. Note 59 is NOT among
them -- confirmation that the position list is right, since 25 x 19 = 475 on the
nose. This script reproduces that library in silico.

THE A/B, which is Jannis's question. Each of the 475 is scored in TWO
backgrounds:
  WT      the benchmark as this project has been running it
  K59R    the benchmark as the experiment was actually done
If forcing K59R improves the rank of V81I/F108A/F159L, then part of what has
been recorded as method failure was benchmark error, and the recorded
conclusion "stage 1 does not clear" has to be revisited.

PROTOCOL. §78 is the licence for this: on crystal poses the score gets the
cross-over right after REPACK and relax destroys it. So repack only, ligand pose
fixed at its crystal coordinates, dG_bind by the same three-term split as 134,
and a PAIRED background repacked in the IDENTICAL shell (§66a's first version
scored against a raw wild type and produced pure repack artefact).

REPLICATE SPREAD IS REPORTED, NOT HIDDEN. PackRotamersMover is stochastic and
nothing in this project has ever measured its noise floor. The mandi decomposition
says V81I+F108A+F159L are worth ~2.7 REU once F108A's clash is relieved; if the
replicate spread is that size, no search over this score can resolve them and
that is the finding.

⚠ Cross-reactivity is not modelled. Jannis's note: a variant with better affinity
may still lose because it responds to the wrong ligands. So a low rank here is
evidence against the SCORE, but a high rank is not proof a variant would have
been selected.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import lib_rosetta as LR                                          # noqa: E402

OUT = os.path.join(ROOT, "results", "park475")
FRAME = os.path.join(ROOT, "data", "stage1", "wt_mandi.pdb")
LIGPDB = os.path.join(ROOT, "data", "stage1", "params", "3UZ_0001.pdb")
LIGPRM = os.path.join(ROOT, "data", "stage1", "params", "3UZ.params")

#: verbatim from the paper's methods, 25 sites; 25 x 19 = 475 = their stated size
PARK25 = {55: "P", 61: "F", 62: "I", 81: "V", 83: "V", 87: "L", 88: "P", 89: "A",
          92: "S", 94: "E", 141: "E", 108: "F", 110: "I", 115: "H", 116: "R",
          117: "L", 120: "Y", 122: "S", 158: "M", 159: "F", 160: "A", 162: "T",
          163: "V", 164: "V", 167: "N"}
TARGETS = {"V81I", "F108A", "F159L"}
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
    ligres = base.total_residue()

    # ---- assert every Park position IS the residue the paper names ------
    bad = {}
    for num, wt in PARK25.items():
        if num not in idx:
            bad[num] = "absent"
            continue
        got = base.residue(idx[num]).name1()
        if got != wt:
            bad[num] = f"{got} not {wt}"
    assert not bad, f"Park position identities disagree with the frame: {bad}"
    assert base.residue(idx[59]).name1() == "K", "position 59 is not K in the frame"
    print(f"all 25 Park positions confirmed; frame {os.path.basename(FRAME)}, "
          f"ligand {base.residue(ligres).name3()}", flush=True)

    e_lig = float(sfxn(pyrosetta.pose_from_pdb(LIGPDB)))

    def dg_bind(p):
        e_c = float(sfxn(p))
        pp = p.clone()
        pp.delete_residue_slow(pp.total_residue())
        return e_c - float(sfxn(pp)) - e_lig

    def repack_score(pose, allowed, nrep):
        vals = []
        for _ in range(nrep):
            p = pose.clone()
            tf, _ = LR.restrict_packing(p, allowed)
            pk = PackRotamersMover(sfxn)
            pk.task_factory(tf)
            pk.apply(p)
            vals.append(dg_bind(p))
        return min(vals), vals

    backgrounds = {}
    backgrounds["WT"] = base.clone()
    k = base.clone()
    MutateResidue(idx[59], "ARG").apply(k)
    backgrounds["K59R"] = k

    jobs = [(num, mu) for num in sorted(PARK25)
            for mu in AA20 if mu != PARK25[num]]
    assert len(jobs) == 475, f"library is {len(jobs)}, paper says 475"
    jobs = [j for i, j in enumerate(jobs) if i % a.nchunks == a.chunk]
    print(f"chunk {a.chunk}/{a.nchunks}: {len(jobs)} variants x "
          f"{len(backgrounds)} backgrounds x {a.nrep} replicates", flush=True)

    rows = []
    for num, mu in jobs:
        rec = {"sub": f"{PARK25[num]}{num}{mu}", "pos": num, "mut": mu}
        for bg, bpose in backgrounds.items():
            pose = bpose.clone()
            MutateResidue(idx[num], THREE[mu]).apply(pose)
            sel = ResidueIndexSelector(str(idx[num]))
            nb = NeighborhoodResidueSelector(sel, 6.0, True)
            allowed = sorted({i + 1 for i, b in enumerate(nb.apply(pose)) if b}
                             | {ligres})
            s_m, all_m = repack_score(pose, allowed, a.nrep)
            s_b, all_b = repack_score(bpose.clone(), allowed, a.nrep)
            rec[bg] = {"dG": s_m, "bg": s_b, "ddG": s_m - s_b,
                       "n_shell": len(allowed),
                       "spread_mut": float(max(all_m) - min(all_m)),
                       "spread_bg": float(max(all_b) - min(all_b))}
        rows.append(rec)
        print(f"  {rec['sub']:<8} WT {rec['WT']['ddG']:+8.2f} "
              f"(spread {rec['WT']['spread_mut']:4.2f})   "
              f"K59R {rec['K59R']['ddG']:+8.2f} "
              f"(spread {rec['K59R']['spread_mut']:4.2f})", flush=True)
        json.dump({"nrep": a.nrep, "rows": rows},
                  open(os.path.join(OUT, f"chunk_{a.chunk:03d}.json"), "w"), indent=1)
    print(f"\nwrote {len(rows)} to {OUT}/chunk_{a.chunk:03d}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
