#!/usr/bin/env python
r"""
157_ddg_per_residue.py -- the per-RESIDUE false-negative rate of the ddG filter.

§89b could only answer Jannis's question at the level of whole variants, because
the scored set was all 3-10 substitution combinations. Job 27925009 supplies what
was missing: all 380 single substitutions at Beltran's 20 design positions, so
each residue can be asked directly -- if we had used ddG to decide which
substitutions go INTO the library, would we have thrown this one away?

40 of the 380 are used by at least one of the 45 real Beltran sensors. Those are
the positives. The other 340 are UNLABELLED, not negatives (§51): Beltran's
libraries did not test every substitution, and an unused one may be perfectly
good. So recall is reported and precision is not.

REPLICATE NOISE, measured for the first time. Median spread over 3 repacks is
0.00 REU and the 90th percentile is 0.21. The signal that has to be resolved --
V81I + F108A + F159L are worth ~2.7 REU once F108A's clash is relieved -- is an
order of magnitude above that. Whatever is wrong with this score, it is not
imprecision.
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "library_recall")


def main():
    rows = []
    for f in glob.glob(os.path.join(ROOT, "results", "beltran_ddg", "chunk_*.json")):
        rows += json.load(open(f))["rows"]
    bel = json.load(open(os.path.join(ROOT, "data", "beltran", "win_sensors.json")))
    sensors = [s for s in bel["sensors"] if s["muts"]]
    used = {m for s in sensors for m in s["muts"]}
    for r in rows:
        r["used"] = r["sub"] in used
    rows.sort(key=lambda r: r["ddG"])
    n, npos = len(rows), sum(r["used"] for r in rows)
    print(f"{n} single substitutions at Beltran's 20 positions; "
          f"{npos} are used by a real sensor")
    sp = np.array([r["mut_spread"] for r in rows])
    print(f"repack replicate spread: median {np.median(sp):.2f}, "
          f"90th pct {np.percentile(sp, 90):.2f} REU  "
          f"-- the signal to resolve is ~2.7 REU")

    print("\n" + "=" * 74)
    print("IF ddG CHOSE THE LIBRARY, WHICH REAL SUBSTITUTIONS WOULD BE LOST?")
    print("=" * 74)
    print(f"{'keep top':>10}{'of 380':>9}{'sensor subs kept':>19}{'recall':>9}"
          f"{'expected':>10}")
    res = {}
    for frac in (0.10, 0.25, 0.50, 0.75, 0.90):
        k = int(round(frac * n))
        kept = sum(r["used"] for r in rows[:k])
        print(f"{frac:>9.0%}{k:>9}{kept:>15}/{npos:<3}"
              f"{kept/npos:>9.0%}{frac*npos:>10.1f}")
        res[f"{frac:.2f}"] = dict(k=k, kept=kept, recall=kept / npos)

    # AUC: probability a used substitution scores better than an unused one
    ranks = {r["sub"]: i for i, r in enumerate(rows)}
    pos = [ranks[r["sub"]] for r in rows if r["used"]]
    neg = [ranks[r["sub"]] for r in rows if not r["used"]]
    auc = float(np.mean([[1.0 if p < q else 0.5 if p == q else 0.0
                          for q in neg] for p in pos]))
    print(f"\n   AUC (used vs unlabelled, lower ddG = better) = {auc:.3f}")
    print(f"   0.5 is chance; >0.5 means ddG favours real substitutions")

    print("\n" + "=" * 74)
    print("THE REAL SUBSTITUTIONS ddG SCORES WORST -- the false negatives")
    print("=" * 74)
    print(f"   {'sub':>8}{'ddG':>9}{'rank':>8}   used by")
    bad = sorted((r for r in rows if r["used"]), key=lambda r: -r["ddG"])[:10]
    for r in bad:
        who = sorted({s["ligand"] for s in sensors if r["sub"] in s["muts"]})
        print(f"   {r['sub']:>8}{r['ddG']:>+9.2f}{ranks[r['sub']] + 1:>5}/380"
              f"   {', '.join(who)[:44]}")
    print(f"\n   best-scoring real substitutions:")
    for r in sorted((r for r in rows if r["used"]), key=lambda r: r["ddG"])[:5]:
        print(f"   {r['sub']:>8}{r['ddG']:>+9.2f}{ranks[r['sub']] + 1:>5}/380")

    # per position: does ddG at least pick the right POSITIONS?
    print("\n" + "=" * 74)
    print("POSITIONS: does ddG rank Beltran's productive positions highly?")
    print("=" * 74)
    bypos = {}
    for r in rows:
        bypos.setdefault(r["pos"], []).append(r)
    order = sorted(bypos, key=lambda p: np.median([x["ddG"] for x in bypos[p]]))
    prod = {int(m[1:-1]) for m in used}
    print(f"   {'rank':>5}{'pos':>6}{'median ddG':>12}{'used?':>8}")
    for i, p in enumerate(order, 1):
        print(f"   {i:>5}{p:>6}{np.median([x['ddG'] for x in bypos[p]]):>+12.2f}"
              f"{('YES' if p in prod else '-'):>8}")
    json.dump({"auc": auc, "recall": res}, open(os.path.join(OUT, "ddg_per_residue.json"), "w"),
              indent=1)
    print(f"\nwritten to {OUT}/ddg_per_residue.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
