#!/usr/bin/env python
r"""
211_position_scan.py -- does HBNet scoring rank the TRUE sensor residue above
wild type, one position at a time?

WHY THIS REPLACES 210's WHOLE-POCKET DESIGN. 210 shaved all 23 lining positions
to alanine and redesigned them together. That reproduced the literature baseline
exactly -- native sequence recovery 0.35 on the ABA arm, which is what Rosetta
gives -- but it is UNDERPOWERED BY CONSTRUCTION for this question: at 0.35 per
position, recovering 4 specific positions by chance is 0.35^4 = 1.5 %, so
`hbnet_terms` scoring 1/4 against FastDesign's 0/4 is inside the noise. Shaving
23 positions at once also opens a cavity large enough that design fills it with
bulk, which is why specificity was 5-8 of 23 in every arm.

A SINGLE-POSITION SCAN asks the question directly. For each lining position,
mutate to all 19 non-proline residues in the otherwise-NATIVE background with the
ligand present, score, and record where the TRUE residue ranks. That isolates the
scoring pathology from the design search, which matters because §35 established
that sampling is NOT the limit -- coupled moves retained K59 0 % of the time.

THE PRE-REGISTERED TARGET is K59R on mandipropamid. §23j measured the cause:
ref2015 pays +10.1 REU to desolvate the salt bridge it correctly rewards, so R
ranks far below K. `hbnet` and `buried_unsatisfied_penalty` exist to make buried
polar networks payable, so if they work, R's rank at 59 improves. If it does not
move, HBNet does not fix this failure mode and the arm closes.

PREDICTIONS, FIXED BEFORE THE RUN:
  1. CONTROL. On the ABA arm the true residue is WILD TYPE at every position, and
     ref2015 must rank it top-3 at most positions. If it does not, the scan is
     mis-set-up and nothing below counts. ⚠ This control is WEAK on its own --
     ref2015's reference energies were fitted to reproduce native sequences.
  2. BASELINE. On mandipropamid, plain ref2015_cart ranks F108A and F159L well
     (they are clash relief, §32) and K59R and V81I poorly. This must reproduce
     the known 2-of-4 result or the scan disagrees with the project record.
  3. THE TEST. hbnet_terms improves the RANK OF R AT POSITION 59 relative to
     baseline. Rank, not a sign, because §40 showed sign-based verdicts are the
     fragile ones.

⚠ Scored as RANK of the true residue among 19, per position -- not as a binary
"recovered", so a near miss is visible instead of being rounded to failure.
"""
import argparse, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = os.path.join(ROOT, "results", "position_scan")
AA19 = "ACDEFGHIKLMNQRSTVWY"          # no proline
LINING = [59, 60, 61, 62, 81, 83, 87, 89, 91, 92, 94, 108, 110, 115, 117,
          120, 122, 141, 159, 160, 163, 164, 167]        # P88 excluded
WT_PYR1 = {59: "K", 60: "H", 61: "F", 62: "I", 81: "V", 83: "V", 87: "L",
           89: "A", 91: "T", 92: "S", 94: "E", 108: "F", 110: "I", 115: "H",
           117: "L", 120: "Y", 122: "S", 141: "E", 159: "F", 160: "A",
           163: "V", 164: "V", 167: "N"}
SYSTEMS = {
    "aba":   dict(pdb="3QN1", lig="A8S", params=f"{ROOT}/data/aba_params/A8S.params"),
    "mandi": dict(pdb="4WVO", lig="3UZ", params=f"{ROOT}/data/stage1/params/3UZ.params"),
    "win":   dict(pdb="7MWN", lig="WI5", params=f"{ROOT}/results/win_crossover/WI5.params"),
}
THREE = dict(A="ALA", C="CYS", D="ASP", E="GLU", F="PHE", G="GLY", H="HIS",
             I="ILE", K="LYS", L="LEU", M="MET", N="ASN", Q="GLN", R="ARG",
             S="SER", T="THR", V="VAL", W="TRP", Y="TYR")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", default="mandi", choices=sorted(SYSTEMS))
    ap.add_argument("--arm", default="ref2015", choices=("ref2015", "hbnet_terms"))
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--nchunks", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    S = SYSTEMS[a.system]
    import pyrosetta
    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ex1 -ex2aro "
                   f"-extra_res_fa {S['params']}")
    from pyrosetta import pose_from_file
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    import importlib
    LR = importlib.import_module("lib_rosetta")

    base = pose_from_file(os.path.join(ROOT, "data", f"{S['pdb']}.cif"))
    pi = base.pdb_info()
    rm = {pi.number(i): i for i in range(1, base.total_residue() + 1)
          if pi.chain(i) == "A"}
    lig_idx = [i for i in range(1, base.total_residue() + 1)
               if base.residue(i).name3().strip() == S["lig"]]
    assert lig_idx, f"ligand {S['lig']} absent from the pose"
    truth = {n: base.residue(rm[n]).name1() for n in LINING if n in rm}
    print(f"{a.system}: {S['pdb']}, ligand {S['lig']} at pose index {lig_idx[0]}, "
          f"{len(truth)}/{len(LINING)} lining positions")
    if a.system != "win":
        d = [f"{WT_PYR1[n]}{n}{truth[n]}" for n in truth if truth[n] != WT_PYR1.get(n)]
        print(f"   truth vs WT PYR1: {d if d else 'wild type'}")

    sf = pyrosetta.create_score_function("ref2015_cart")
    if a.arm == "hbnet_terms":
        from pyrosetta.rosetta.core.scoring import ScoreType
        sf.set_weight(ScoreType.hbnet, 1.0)
        sf.set_weight(ScoreType.buried_unsatisfied_penalty, 1.0)
    print(f"   arm {a.arm}")

    pos = [n for n in LINING if n in rm]
    pos = [p for i, p in enumerate(pos) if i % a.nchunks == a.chunk]
    rows = []
    for n in pos:
        sel = ResidueIndexSelector(str(rm[n]))
        sh = NeighborhoodResidueSelector(sel, 6.0, True).apply(base)
        allowed = [i for i in range(1, base.total_residue() + 1) if sh[i]]
        allowed += [i for i in lig_idx if i not in allowed]
        scores = {}
        for aa in AA19:
            p = base.clone()
            MutateResidue(rm[n], THREE[aa]).apply(p)
            tf = LR.restrict_packing(p, allowed)
            pk = PackRotamersMover(sf)
            pk.task_factory(tf)
            pk.apply(p)
            scores[aa] = float(sf(p))
        order = sorted(scores, key=lambda x: scores[x])
        t = truth[n]
        rows.append(dict(pos=n, truth=t, wt=WT_PYR1.get(n),
                         rank_true=order.index(t) + 1,
                         rank_wt=order.index(WT_PYR1[n]) + 1 if WT_PYR1.get(n) in order else None,
                         best=order[0], order="".join(order),
                         d_true_best=round(scores[t] - scores[order[0]], 2)))
        r = rows[-1]
        mark = "  <-- TRUE != WT" if r["wt"] and r["truth"] != r["wt"] else ""
        print(f"   {n:>4} truth {t}  rank {r['rank_true']:>2}/19  "
              f"(wt {r['wt']} rank {r['rank_wt']})  best {r['best']}  "
              f"dE {r['d_true_best']:+.2f}{mark}", flush=True)
    out = os.path.join(OUT, f"{a.system}_{a.arm}_{a.chunk}.json")
    json.dump(rows, open(out, "w"), indent=1)
    print(f"   wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
