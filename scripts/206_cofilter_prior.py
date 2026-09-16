#!/usr/bin/env python
r"""
206_cofilter_prior.py -- condition on ligands that SHARE A HIT, not on ligands
that are chemically similar.

WHY THIS IS NOT ATTEMPT SIX OF THE SAME THING. Five conditioning attempts have
failed and every one used a CHEMICAL neighbourhood: Tanimoto lookup (52, dead --
only 0.4 % of pairs reach T >= 0.5), similarity-weighted aggregation (103, = 
lookup), the gated version (109), ligand SIZE (115d, retracted on a subset
control) and volume deficit (200, net negative). Meanwhile 201 measured that
ligand identity really does structure which substitutions get selected -- 4.57x
within- over between-ligand Jaccard after deduping -- so the information exists
and the conditioning VARIABLE is what is wrong.

This uses an empirical neighbourhood instead: two ligands are neighbours if their
round-1 hits SHARE a substitution. Content-based vs collaborative filtering. It
needs no chemistry at all, which is the point, since 114a showed chemistry
matching is null and 114e showed the real chemical signal runs ANTI-complementary.

THE SETUP IS THE REAL EXPERIMENT, not a benchmark convenience. For a held-out
ligand we reveal ONE of its clones -- the weakest-responding one, which is what a
primary screen actually hands you -- and predict the substitutions in its OTHER
clones. That is the secondary-library problem exactly: 60 says use a hit to
WEIGHT, not to restrict, and pooled round-1 covers 75-100 % of round-2 while
own-ligand carryover is only 0-64 %.

  seed        = the revealed clone's substitutions (eugenol: V164L+N167V)
  neighbours  = every OTHER ligand with a clone sharing >= 1 seed substitution
  prior       = substitution frequency over neighbour clones
  truth       = the held-out ligand's remaining substitutions, seed REMOVED
  baseline    = ligand-blind frequency over all training clones, same truth

⚠ Scored recall-only. Clones are HITS FOUND, positive-unlabeled (51).
⚠ The held-out ligand's own clones are excluded from both prior and neighbours.
⚠ The seed substitutions are removed from the truth so the method cannot score
   by echoing what it was given.

PREDICTIONS, FIXED BEFORE THE RUN:
  1. The co-occurrence prior beats the ligand-blind prior at recall@20, paired
     across ligands, sign test p < 0.05.
  2. It survives a NEIGHBOURHOOD-PERMUTATION control: reassigning each ligand a
     random neighbourhood of equal size must not reproduce the gain. This is the
     115d lesson -- control against a random set of the SAME SIZE.
  3. It does not merely reproduce library structure: the gain survives holding
     mut_lib fixed.
"""
import json, os, re, sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                          # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "cofilter_prior")
AA = set("ACDEFGHIKLMNPQRSTVWY")


def load():
    hdr, recs = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    pc = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    out = []
    for r in recs:
        lg = (r.get("library_name") or "").strip()
        s = frozenset(c + (r.get(c) or "").strip().upper() for c in pc
                      if len((r.get(c) or "").strip()) == 1
                      and (r.get(c) or "").strip().upper() in AA
                      and (r.get(c) or "").strip().upper() != c[0])
        mc = (r.get("min_conc") or "").strip()
        try:
            v = float(mc)
        except ValueError:
            v = None
        if lg and s:
            out.append(dict(lig=lg, subs=s, mc=v,
                            ml=(r.get("mut_lib") or "").strip().lower()))
    return out


def rank(clones):
    return [x for x, _ in Counter(s for c in clones for s in c["subs"]).most_common()]


def recall(r, truth, n):
    return len(set(r[:n]) & truth) / max(len(truth), 1)


def main():
    os.makedirs(OUT, exist_ok=True)
    clones = load()
    by = defaultdict(list)
    for c in clones:
        by[c["lig"]].append(c)
    multi = {k: v for k, v in by.items() if len(v) >= 2}
    print(f"{len(clones)} clones, {len(by)} ligands, {len(multi)} with >=2 clones\n")

    rng = np.random.default_rng(0)
    rows = []
    for lig, cl in multi.items():
        train = [c for c in clones if c["lig"] != lig]
        # the revealed clone: weakest responder, which is what a screen hands you
        seed_c = max(cl, key=lambda c: (c["mc"] if c["mc"] is not None else -1))
        seed = set(seed_c["subs"])
        truth = {s for c in cl for s in c["subs"]} - seed
        if not truth:
            continue
        nb = [c for c in train if seed & c["subs"]]
        if len(nb) < 3:
            continue
        base = rank(train)
        co = rank(nb)
        # control: a RANDOM neighbourhood of the SAME SIZE (the 115d lesson)
        idx = rng.choice(len(train), size=min(len(nb), len(train)), replace=False)
        ctl = rank([train[i] for i in idx])
        rows.append(dict(lig=lig, n_seed=len(seed), n_truth=len(truth),
                         n_nb=len(nb), nb_ligs=len({c["lig"] for c in nb}),
                         b20=recall(base, truth, 20), c20=recall(co, truth, 20),
                         r20=recall(ctl, truth, 20),
                         b40=recall(base, truth, 40), c40=recall(co, truth, 40),
                         r40=recall(ctl, truth, 40), ml=seed_c["ml"]))
    json.dump(rows, open(os.path.join(OUT, "rows.json"), "w"), indent=1)
    from scipy.stats import wilcoxon
    print(f"{len(rows)} ligands evaluated "
          f"(median neighbourhood {np.median([r['n_nb'] for r in rows]):.0f} clones "
          f"from {np.median([r['nb_ligs'] for r in rows]):.0f} ligands)\n")

    def rep(sel, lab):
        if len(sel) < 5:
            print(f"  {lab}: n={len(sel)}, too few"); return
        for k in ("20", "40"):
            b = np.array([x["b" + k] for x in sel])
            c = np.array([x["c" + k] for x in sel])
            r = np.array([x["r" + k] for x in sel])
            d = c - b
            w, l = int((d > 0).sum()), int((d < 0).sum())
            p = wilcoxon(c, b).pvalue if np.any(d) else 1.0
            pc_ = wilcoxon(c, r).pvalue if np.any(c - r) else 1.0
            print(f"  {lab:<22}@{k}  blind {b.mean():.3f}  CO-OCC {c.mean():.3f}  "
                  f"rand-nb {r.mean():.3f}   delta {d.mean():+.3f}  "
                  f"W/L {w}/{l}  p={p:.3g}  vs-control p={pc_:.3g}")

    print("P1/P2 -- co-occurrence vs blind, and vs an equal-size random neighbourhood")
    rep(rows, "ALL ligands")
    print("\nP3 -- held within mut_lib")
    for ml in sorted({r["ml"] for r in rows if r["ml"]}):
        rep([r for r in rows if r["ml"] == ml], f"mut_lib={ml}")
    eu = [r for r in rows if "eugenol" in r["lig"].lower()]
    if eu:
        e = eu[0]
        print(f"\nEUGENOL: seed {e['n_seed']} subs, {e['n_truth']} to predict, "
              f"neighbourhood {e['n_nb']} clones / {e['nb_ligs']} ligands")
        print(f"  blind@20 {e['b20']:.2f}  co-occurrence@20 {e['c20']:.2f}  "
              f"random-nb@20 {e['r20']:.2f}")
    print(f"\nwrote {OUT}/rows.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
