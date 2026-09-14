#!/usr/bin/env python
r"""
194_size_conditioned_prior.py -- does conditioning the prior on LIGAND SIZE beat
the ligand-blind bar on Beltran-45?

§121c found that PYR1's productive positions depend on ligand size. Splitting the
194 sd03 ligands at ABA's 19 heavy atoms:

    F159   48.1 % of small-ligand clones (rank 1)  vs  44.5 % large (rank 1)
    Y120    9.2 % (rank 10)                        vs  32.6 % (rank 2)
    K59    24.8 % (rank 2)                         vs  18.2 % (rank 6)
    V163   12.2 % (rank 8)                         vs   6.3 % (rank 15)

and the direction of substitution flips: small ligands FILL the pocket
(mean +1.7 A^3, 51 % grow) while large ligands OPEN it (-11.3 A^3, 44 % grow).

WHY THIS MIGHT SUCCEED WHERE SIMILARITY FAILED. §103/§104 closed three
ligand-conditioning attempts, all for the same reason: 194 ligands at median
pairwise Tanimoto 0.11 is too sparse a chemical space to condition on. Heavy-atom
count is not sparse -- it is a dense, reliably-measured scalar with 79 ligands
below ABA and 115 at or above. The failure mode that killed similarity does not
obviously apply.

⚠ PRE-REGISTERED, so this cannot be read as post-hoc. The bar is the SAME one
§103/§109 used and failed to beat out of sample:

    ligand-blind frequency prior, recall@20 = 0.429 and recall@40 = 0.660 on
    Beltran-45 (§109, measured on this exact test set)

The size-conditioned prior must beat BOTH, on the same 14 Beltran ligands, with
nothing fitted. The split point is ABA's 19 heavy atoms, chosen in §121c before
this test existed and NOT tuned here. A sensitivity sweep over the split point is
reported afterwards and is explicitly labelled exploratory.

⚠ §87b's ceiling still applies: Tian's whole vocabulary captures only 30/45 of
Beltran's sensors, so re-weighting cannot exceed 67 %.
"""
import json, os, re, sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                          # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "ligand_prior")
AA = set("ACDEFGHIKLMNPQRSTVWY")
SPLIT = 19          # ABA's heavy-atom count, fixed in §121c

#: Beltran's 14 ligands, same SMILES table as §109 so the two are comparable
sys.path.insert(0, HERE)
BEL = import_bel = None


def main():
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m185 = __import__("importlib").import_module("185_gated_prior_beltran")
    BEL_SMILES = m185.BEL_SMILES

    hdr, recs = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    pc = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    clones = []
    for r in recs:
        lg = (r.get("library_name") or "").strip()
        sm = (r.get("canonical_smiles") or "").strip()
        s = [c + (r.get(c) or "").strip().upper() for c in pc
             if len((r.get(c) or "").strip()) == 1
             and (r.get(c) or "").strip().upper() in AA
             and (r.get(c) or "").strip().upper() != c[0]]
        if lg and s and sm:
            m = Chem.MolFromSmiles(sm)
            if m:
                clones.append((m.GetNumHeavyAtoms(), s))
    print(f"training: {len(clones)} sd03 clones with size and substitutions")

    def prior(sel):
        c = Counter(x for _, s in sel for x in s)
        return [x for x, _ in c.most_common()]

    blind = prior(clones)
    small = prior([c for c in clones if c[0] < SPLIT])
    large = prior([c for c in clones if c[0] >= SPLIT])
    print(f"  blind {len(blind)} subs | small(<{SPLIT}) {len(small)} from "
          f"{sum(1 for c in clones if c[0]<SPLIT)} clones | "
          f"large {len(large)} from {sum(1 for c in clones if c[0]>=SPLIT)}")

    bel = json.load(open(os.path.join(ROOT, "data", "beltran", "win_sensors.json")))
    sensors = [s for s in bel["sensors"] if s["muts"]]
    truth = defaultdict(set)
    for s in sensors:
        truth[s["ligand"]] |= set(s["muts"])

    def rec(rank, t, n):
        return len(set(rank[:n]) & t) / max(len(t), 1)

    rows = []
    for lg in sorted(truth):
        m = Chem.MolFromSmiles(BEL_SMILES[lg])
        na = m.GetNumHeavyAtoms()
        cond = small if na < SPLIT else large
        rows.append(dict(lig=lg, na=na, arm="small" if na < SPLIT else "large",
                         n=len(truth[lg]),
                         b20=rec(blind, truth[lg], 20), c20=rec(cond, truth[lg], 20),
                         b40=rec(blind, truth[lg], 40), c40=rec(cond, truth[lg], 40)))
    print(f"\n{'ligand':<20}{'atoms':>6}{'arm':>7}{'n':>4}"
          f"{'blind@20':>10}{'size@20':>9}{'blind@40':>10}{'size@40':>9}")
    for r in sorted(rows, key=lambda x: x["na"]):
        print(f"{r['lig']:<20}{r['na']:>6}{r['arm']:>7}{r['n']:>4}"
              f"{r['b20']:>10.2f}{r['c20']:>9.2f}{r['b40']:>10.2f}{r['c40']:>9.2f}")

    rng = np.random.default_rng(0)
    print(f"\n{'':<24}{'blind (bar)':>13}{'size-cond':>12}{'delta':>9}{'p':>9}")
    for k, lab in (("20", "recall@20"), ("40", "recall@40")):
        b = np.array([r["b" + k] for r in rows])
        c = np.array([r["c" + k] for r in rows])
        d = c - b
        null = np.array([np.mean(d * rng.choice([-1, 1], len(d)))
                         for _ in range(20000)])
        p = float((np.abs(null) >= abs(d.mean())).mean())
        print(f"  {lab:<22}{b.mean():>13.3f}{c.mean():>12.3f}"
              f"{d.mean():>+9.3f}{p:>9.4f}")
        print(f"  {'':22}win/loss/tie "
              f"{int((d>0).sum())}/{int((d<0).sum())}/{int((d==0).sum())}")

    print(f"\n  ⚠ EXPLORATORY split-point sweep (the pre-registered test is "
          f"SPLIT={SPLIT} above):")
    for sp in (14, 16, 19, 22, 25):
        sm = prior([c for c in clones if c[0] < sp])
        lg_ = prior([c for c in clones if c[0] >= sp])
        v = []
        for r in rows:
            m = Chem.MolFromSmiles(BEL_SMILES[r["lig"]])
            cond = sm if m.GetNumHeavyAtoms() < sp else lg_
            v.append(rec(cond, truth[r["lig"]], 20) - r["b20"])
        nsm = sum(1 for r in rows if Chem.MolFromSmiles(BEL_SMILES[r['lig']]).GetNumHeavyAtoms() < sp)
        print(f"     split {sp:>2}: delta@20 {np.mean(v):+.3f}   "
              f"({nsm}/{len(rows)} Beltran ligands fall in the small arm)")
    json.dump(rows, open(os.path.join(OUT, "size_conditioned.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
