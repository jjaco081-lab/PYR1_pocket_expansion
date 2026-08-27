#!/usr/bin/env python
r"""
139_tractability_set.py -- build a PROPERTY-MATCHED hit / non-hit ligand set.

THE QUESTION. §82 showed Boltz-2 cannot rank pocket VARIANTS for a fixed ligand --
a direction its training never constrains. The ligand-level question is different
and runs WITH the training direction: given a molecule, will PYR1 yield a sensor
for it at all? sd01 has the labels: **194 Hit? = Yes, 3,172 Hit? = No**, all with
SMILES. That is the first genuine negative set this project has had.

⚠ JANNIS'S CONSTRAINT, AND IT IS THE WHOLE DESIGN. A random draw of non-hits is a
rigged benchmark: PYR1's pocket is ~174 A^3, so much of the Selleck/LATCA deck is
simply too large, and a model could "succeed" for a reason `heavy_atoms` alone
would also deliver. Any such result would say nothing about chemistry.

So the negatives are MATCHED to the hits on the trivial axes -- heavy atoms, cLogP,
TPSA, H-bond donors and acceptors, formal charge, rotatable bonds -- by greedy
nearest-neighbour on z-scored descriptors, without replacement. Two negative sets
are written:

    matched   the real test; trivial descriptors should NOT separate it
    random    the naive comparator, kept precisely to SHOW how easy it looks

and the script reports, for both, how well each single descriptor separates hits
from non-hits. **If the matched set still separates on size, the matching failed
and the benchmark is void** -- that check is printed before anything is exported.

RECEPTOR CHOICE, and its honest limitation. A "hit" is a ligand for which SOME
library variant worked -- not one that binds wild-type PYR1. Giving hits their
evolved sensor sequence while non-hits get wild-type would be a fatal asymmetry, so
**both classes are co-folded against wild-type PYR1**. The question therefore
becomes "does WT-PYR1 co-folding predict LIBRARY tractability", which is a proxy,
one step removed from the label. It is symmetric, which is what matters.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                      # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "data", "tractability")
DESC = ["heavy", "clogp", "tpsa", "hbd", "hba", "charge", "rotb"]


def props(smiles):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None or m.GetNumHeavyAtoms() < 5:
        return None
    return {"heavy": m.GetNumHeavyAtoms(), "clogp": Crippen.MolLogP(m),
            "tpsa": rdMolDescriptors.CalcTPSA(m),
            "hbd": rdMolDescriptors.CalcNumHBD(m),
            "hba": rdMolDescriptors.CalcNumHBA(m),
            "charge": Chem.GetFormalCharge(m),
            "rotb": rdMolDescriptors.CalcNumRotatableBonds(m)}


def auc(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.mean([(x > y) + 0.5 * (x == y) for x in a for y in b]))


def main():
    os.makedirs(OUT, exist_ok=True)
    _hdr, recs = table(f"{SD}/pnas.2519924122.sd01.xlsx")
    hits, nons = [], []
    seen = set()
    for r in recs:
        s = (r.get("canonical_smiles") or "").strip()
        h = (r.get("Hit?") or "").strip()
        if not s or h not in ("Yes", "No") or s in seen:
            continue
        p = props(s)
        if p is None:
            continue
        seen.add(s)
        rec = {"name": (r.get("Library_name") or "?").strip(), "smiles": s, **p}
        (hits if h == "Yes" else nons).append(rec)
    print(f"unique parseable molecules: {len(hits)} hits, {len(nons)} non-hits")

    # z-score on the pooled distribution so no descriptor dominates the distance
    allr = hits + nons
    mu = {k: float(np.mean([r[k] for r in allr])) for k in DESC}
    sd = {k: float(np.std([r[k] for r in allr])) or 1.0 for k in DESC}

    def z(r):
        return np.array([(r[k] - mu[k]) / sd[k] for k in DESC])

    H = np.array([z(r) for r in hits])
    N = np.array([z(r) for r in nons])
    # greedy nearest-neighbour matching without replacement
    taken = set()
    matched = []
    order = np.argsort([np.linalg.norm(h) for h in H])[::-1]   # hardest first
    for i in order:
        d = np.linalg.norm(N - H[i], axis=1)
        d[list(taken)] = np.inf
        j = int(np.argmin(d))
        taken.add(j)
        matched.append(nons[j])
    rng = np.random.default_rng(0)
    rand = [nons[i] for i in rng.choice(len(nons), len(hits), replace=False)]

    print(f"\nCAN A TRIVIAL DESCRIPTOR ALONE SEPARATE HITS FROM NON-HITS?")
    print(f"   AUC 0.5 = no separation. The MATCHED set must be near 0.5 or the")
    print(f"   benchmark is rigged and any model result would be uninterpretable.")
    print(f"   {'descriptor':<12}{'random negatives':>18}{'MATCHED negatives':>20}")
    worst = 0.0
    for k in DESC:
        a = [r[k] for r in hits]
        ar = auc(a, [r[k] for r in rand])
        am = auc(a, [r[k] for r in matched])
        worst = max(worst, abs(am - 0.5))
        print(f"   {k:<12}{ar:>18.3f}{am:>20.3f}")
    print(f"\n   worst |AUC-0.5| in the matched set: {worst:.3f}")
    if worst > 0.10:
        print("   ⚠ MATCHING INSUFFICIENT -- a descriptor still separates. Do not run.")
    else:
        print("   ✅ no single descriptor separates the matched set; the benchmark")
        print("      cannot be won on size or polarity alone.")
    for k in ("heavy", "clogp"):
        print(f"   {k}: hits median {np.median([r[k] for r in hits]):.1f}, "
              f"matched {np.median([r[k] for r in matched]):.1f}, "
              f"random {np.median([r[k] for r in rand]):.1f}")

    for tag, neg in (("matched", matched), ("random", rand)):
        recs_out = ([{**r, "label": 1} for r in hits] + [{**r, "label": 0} for r in neg])
        with open(os.path.join(OUT, f"set_{tag}.json"), "w") as fh:
            json.dump(recs_out, fh, indent=1)
        print(f"   wrote {len(recs_out)} molecules to {OUT}/set_{tag}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
