#!/usr/bin/env python
r"""
201_specificity_ceiling.py -- OPTION C. Is a clone's substitution set predictable
from its ligand AT ALL?

WHY THIS IS THE RIGHT TEST NOW. Five ligand-conditioning attempts have failed:
similarity lookup (52), the MSA-style prior (103, = lookup), the gated prior
(109, +0.031 n.s.), the size-conditioned prior (115d, RETRACTED on a subset
control) and fill-matching (200, net NEGATIVE, 0 wins/6 losses even in its
parameter-free form). Rather than try a sixth, measure the CEILING: if two clones
selected for the SAME ligand are no more similar to each other than two clones
selected for DIFFERENT ligands, then no ligand-conditioned re-ranking can work,
and all five failures have one explanation.

The eugenol case is the motivating anecdote: its hit V164L+N167V also appears for
flufenamic acid, 2-hydroxy-9-fluorenone, efavirenz, chloroxylenol and
5-bromo-2-hydroxybenzophenone -- six chemically unrelated ligands, one solution.

TEST. Jaccard similarity between substitution sets. Compare the mean within-ligand
pair similarity against the mean between-ligand pair similarity, with a label
permutation null (ligand labels shuffled, clone sets untouched, n = 2000).

⚠ CONFOUND, handled: clones of one ligand often come from the SAME mut_lib (dsm
or tsm), which constrains both the positions available and the substitution
depth. Two clones from one library are similar for reasons that have nothing to
do with the ligand. Every comparison is therefore ALSO run WITHIN a single
mut_lib, exactly as 49a did.

PREDICTIONS, FIXED BEFORE THE RUN:
  1. Within-ligand Jaccard EXCEEDS between-ligand, permutation p < 0.05. If this
     fails, the ceiling is zero and every ligand-conditioned method is dead.
  2. The effect SHRINKS but survives when mut_lib is held fixed. If it vanishes,
     the apparent ligand signal was library structure.
  3. The six V164L+N167V ligands show no more mutual similarity than chance --
     they are a generic solution, not a chemical family.
"""
import json, os, re, sys
from collections import defaultdict
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                          # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "specificity_ceiling")
AA = set("ACDEFGHIKLMNPQRSTVWY")
FAMILY = ["eugenol", "flufenamic", "fluorenone", "efavirenz",
          "chloroxylenol", "bromo-2-hydroxybenzophenone"]


def load():
    hdr, recs = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    assert "library_name" in hdr and "mut_lib" in hdr
    pc = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    out = []
    for r in recs:
        s = frozenset(c + (r.get(c) or "").strip().upper() for c in pc
                      if len((r.get(c) or "").strip()) == 1
                      and (r.get(c) or "").strip().upper() in AA
                      and (r.get(c) or "").strip().upper() != c[0])
        lg = (r.get("library_name") or "").strip()
        ml = (r.get("mut_lib") or "").strip().lower()
        if s and lg:
            out.append((lg, ml, s))
    return out


def jac(a, b):
    return len(a & b) / max(len(a | b), 1)


def within_between(clones, rng, nperm=2000):
    """(within_mean, between_mean, p) with a ligand-label permutation null."""
    by = defaultdict(list)
    for lg, _, s in clones:
        by[lg].append(s)
    multi = {k: v for k, v in by.items() if len(v) >= 2}
    if len(multi) < 3:
        return None
    win = [jac(a, b) for v in multi.values() for a, b in combinations(v, 2)]
    keys = list(by)
    allpairs = [(i, j) for i in range(len(keys)) for j in range(i + 1, len(keys))]
    rng.shuffle(allpairs)
    btw = []
    for i, j in allpairs[:4000]:
        for a in by[keys[i]][:3]:
            for b in by[keys[j]][:3]:
                btw.append(jac(a, b))
    obs = float(np.mean(win)) - float(np.mean(btw))
    # permutation: shuffle ligand labels across clones, keep group sizes
    sets = [s for _, _, s in clones]
    labels = [lg for lg, _, _ in clones]
    null = []
    for _ in range(nperm):
        perm = rng.permutation(len(labels))
        pb = defaultdict(list)
        for k, idx in enumerate(perm):
            pb[labels[k]].append(sets[idx])
        pm = {k: v for k, v in pb.items() if len(v) >= 2}
        w = [jac(a, b) for v in pm.values() for a, b in combinations(v, 2)]
        null.append(float(np.mean(w)) if w else 0.0)
    null = np.array(null)
    p = float((null >= float(np.mean(win))).mean())
    return dict(within=float(np.mean(win)), between=float(np.mean(btw)),
                delta=obs, n_lig=len(multi), n_pairs=len(win),
                null_mean=float(null.mean()), p=p)


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(0)
    clones = load()
    print(f"sd03: {len(clones)} clones with >=1 substitution")
    res = {}

    print("\nP1 ALL clones, ligand labels permuted (n=2000)")
    r = within_between(clones, rng)
    res["all"] = r
    print(f"   within-ligand Jaccard {r['within']:.4f} vs between-ligand "
          f"{r['between']:.4f}   (null within {r['null_mean']:.4f})")
    print(f"   {r['n_lig']} ligands with >=2 clones, {r['n_pairs']} within pairs, "
          f"p = {r['p']:.4f}   "
          f"{'HOLDS' if r['p'] < 0.05 else 'FAILS -- ceiling is ZERO'}")

    print("\nP2 held WITHIN each mut_lib (the 49a control)")
    for ml in sorted({m for _, m, _ in clones if m}):
        sub = [c for c in clones if c[1] == ml]
        rr = within_between(sub, rng, nperm=1000)
        if rr is None:
            print(f"   {ml}: too few"); continue
        res[ml] = rr
        print(f"   {ml:<5} n={len(sub):>4}  within {rr['within']:.4f} vs "
              f"between {rr['between']:.4f}  null {rr['null_mean']:.4f}  "
              f"p = {rr['p']:.4f}  "
              f"{'holds' if rr['p'] < 0.05 else 'FAILS'}")

    print("\nP3 the six V164L+N167V ligands")
    fam = [c for c in clones
           if any(f in c[0].lower() for f in FAMILY)]
    names = sorted({c[0] for c in fam})
    print(f"   matched {len(fam)} clones over {len(names)} names: {names}")
    tgt = frozenset({"V164L", "N167V"})
    carriers = [c for c in clones if tgt <= c[2]]
    print(f"   clones carrying BOTH V164L and N167V: {len(carriers)} over "
          f"{len({c[0] for c in carriers})} distinct ligands")
    if fam:
        fj = [jac(a[2], b[2]) for a, b in combinations(fam, 2)]
        allj = []
        sets = [c[2] for c in clones]
        for _ in range(4000):
            i, j = rng.integers(0, len(sets), 2)
            if i != j:
                allj.append(jac(sets[i], sets[j]))
        print(f"   family mutual Jaccard {np.mean(fj):.4f} vs random clone pairs "
              f"{np.mean(allj):.4f}")
        res["family"] = dict(family=float(np.mean(fj)),
                             random=float(np.mean(allj)), n=len(fj))

    json.dump(res, open(os.path.join(OUT, "ceiling.json"), "w"), indent=1)
    print(f"\nwrote {OUT}/ceiling.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
