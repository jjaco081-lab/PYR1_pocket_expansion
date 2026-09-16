#!/usr/bin/env python
r"""
207_ensemble_marginals.py -- CALIBRATION. Can per-substitution MARGINAL effects,
averaged over combinatorial backgrounds, shrink a library while retaining hits?

JANNIS'S IDEA: compute many side-chain combinations and average, rather than
scoring substitutions one at a time.

WHY IT IS WELL-POSED HERE. Every scoring method in this project evaluated
substitutions in the WILD-TYPE background, but a library is combinatorial and a
substitution's value alone is not its value in context. Two measured facts say
context is where the signal is:
  * 69: Cartesian FastRelax ddG RANKS COMBINATIONS at AUC 0.251 on the coumarin
    arm, whose target ligand is 162.7 A^3 -- eugenol is 162 A^3, the same regime.
    (It cannot choose the MENU, 0/11; that is a different job.)
  * 68b: the separation WIDENS with depth, 0.350 at 4 substitutions -> 0.288 at
    5 -> 0.145 at 6. Deep backgrounds are more informative than singles.

⚠ 33's trap, and why an ensemble escapes it: a SUMMED ligand-overlap objective
"exposes nothing -- relief of every top pair equalled the sum of its singles".
An additive objective cannot show synergy. This does not assume additivity; it
MEASURES the departure from it.

CALIBRATION, on data that already exists (results/cavity_match/coumarin_*.json):
  P1 reproduce the published separation: REAL vs LIBRARY on ddG, depth-matched,
     must land near AUC 0.25. If not, the machinery is not reading the same
     numbers and nothing below counts.
  P2 marginal effect per substitution, mean(ddG | s present) - mean(ddG | absent),
     must rank the substitutions REAL sensors actually use above the rest.
  P3 the retention curve: pruning the menu by marginal effect, what fraction of
     REAL sensors survive at what library reduction? The bar is 69d's ~2x at
     90 % retention.

⚠ REAL sensors are HITS FOUND (positive-unlabeled), so this is recall only.
⚠ ddG's known pathology -- +10.1 REU to desolvate a buried salt bridge, which is
   why K59R is unreachable -- biases against buried POLAR substitutions. That is
   awkward for eugenol, whose one small-ligand signal (L117->N/D) is polar, so
   the latch set should be FORCED rather than voted on.
"""
import glob, json, os, sys
from collections import defaultdict
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "ensemble_marginals")


def load(klass):
    rows, meta = [], None
    for f in sorted(glob.glob(os.path.join(ROOT, "results", "cavity_match",
                                           f"{klass}_*.json"))):
        d = json.load(open(f))
        meta = meta or dict(target_vol=d.get("target_vol"),
                            wt_cavity=d.get("wt_cavity"))
        rows.extend(d["rows"])
    return rows, meta


def auc(pos, neg):
    """P(pos < neg) -- ddG LOWER is better, so <0.5 means REAL separates."""
    if not pos or not neg:
        return float("nan")
    p, n = np.array(pos), np.array(neg)
    gt = (p[:, None] < n[None, :]).sum() + 0.5 * (p[:, None] == n[None, :]).sum()
    return float(gt / (len(p) * len(n)))


def main():
    os.makedirs(OUT, exist_ok=True)
    klass = sys.argv[1] if len(sys.argv) > 1 else "coumarin"
    rows, meta = load(klass)
    rows = [r for r in rows if r.get("ddG") is not None]
    real = [r for r in rows if r["set"] == "REAL"]
    lib = [r for r in rows if r["set"] == "LIBRARY"]
    print(f"{klass}: {len(rows)} scored variants "
          f"({len(real)} REAL, {len(lib)} LIBRARY, "
          f"{sum(1 for r in rows if r['set']=='WILD')} WILD)")
    print(f"  target ligand volume {meta['target_vol']} A^3   "
          f"WT cavity {meta['wt_cavity']:.1f} A^3")
    print(f"  (eugenol is 162 A^3 -- the same regime)\n")

    # ---- P1 ----------------------------------------------------------------
    a_raw = auc([r["ddG"] for r in real], [r["ddG"] for r in lib])
    ns = sorted({r["n_sub"] for r in real} & {r["n_sub"] for r in lib})
    per, wts = [], []
    for n in ns:
        rr = [r["ddG"] for r in real if r["n_sub"] == n]
        ll = [r["ddG"] for r in lib if r["n_sub"] == n]
        if rr and ll:
            per.append(auc(rr, ll)); wts.append(len(rr))
    a_m = float(np.average(per, weights=wts)) if per else float("nan")
    print(f"P1 REAL vs LIBRARY on ddG: raw AUC {a_raw:.3f}, "
          f"depth-matched {a_m:.3f}  (published: 0.280 raw, 0.251 matched)")
    print(f"   {'HOLDS' if a_m < 0.40 else 'FAILS -- not reading the same numbers'}\n")

    # ---- P2 marginal effects over backgrounds -------------------------------
    subs = sorted({s for r in rows for s in r["subs"]})
    marg = {}
    for s in subs:
        wi = [r["ddG"] for r in rows if s in r["subs"]]
        wo = [r["ddG"] for r in rows if s not in r["subs"]]
        if len(wi) >= 3:
            marg[s] = (float(np.mean(wi) - np.mean(wo)), len(wi))
    real_subs = {s for r in real for s in r["subs"]}
    print(f"P2 marginal effect, {len(marg)} substitutions seen in >=3 backgrounds")
    ranked = sorted(marg, key=lambda s: marg[s][0])
    ranks = [i for i, s in enumerate(ranked) if s in real_subs]
    print(f"   substitutions used by REAL sensors: {len(real_subs & set(marg))}"
          f" of {len(marg)}")
    if ranks:
        print(f"   their median rank {np.median(ranks):.0f} of {len(ranked)}"
              f"   (chance = {len(ranked)/2:.0f})")
        from scipy.stats import mannwhitneyu
        other = [i for i, s in enumerate(ranked) if s not in real_subs]
        if other:
            p = mannwhitneyu(ranks, other, alternative="less").pvalue
            print(f"   Mann-Whitney (REAL subs rank better) p = {p:.4g}"
                  f"   {'HOLDS' if p < 0.05 else 'NULL'}")
    print(f"\n   best 8 by marginal ddG:  "
          f"{', '.join(f'{s}({marg[s][0]:+.2f})' for s in ranked[:8])}")
    print(f"   worst 5:                 "
          f"{', '.join(f'{s}({marg[s][0]:+.2f})' for s in ranked[-5:])}")

    # ---- P3 retention curve -------------------------------------------------
    print(f"\nP3 prune the menu by marginal effect, keep the best K substitutions")
    print(f"   {'K':>4}{'menu kept':>11}{'REAL retained':>15}{'LIBRARY kept':>14}"
          f"{'enrichment':>12}")
    for K in (len(ranked), 30, 25, 20, 16, 12, 10, 8):
        keep = set(ranked[:K])
        rr = sum(1 for r in real if set(r["subs"]) <= keep) / max(len(real), 1)
        ll = sum(1 for r in lib if set(r["subs"]) <= keep) / max(len(lib), 1)
        e = rr / ll if ll > 0 else float("inf")
        print(f"   {K:>4}{K/len(ranked):>10.0%}{rr:>14.0%}{ll:>13.0%}{e:>11.2f}x")

    # ---- pairwise interaction (the non-additivity the ensemble can see) -----
    print(f"\nPAIRWISE INTERACTION  mean(both) - mean(a) - mean(b) + mean(neither)")
    inter = []
    for a, b in combinations(subs, 2):
        both = [r["ddG"] for r in rows if a in r["subs"] and b in r["subs"]]
        oa = [r["ddG"] for r in rows if a in r["subs"] and b not in r["subs"]]
        ob = [r["ddG"] for r in rows if b in r["subs"] and a not in r["subs"]]
        nn = [r["ddG"] for r in rows if a not in r["subs"] and b not in r["subs"]]
        if len(both) >= 3 and oa and ob and nn:
            inter.append((a, b, float(np.mean(both) - np.mean(oa)
                                      - np.mean(ob) + np.mean(nn)), len(both)))
    inter.sort(key=lambda x: x[2])
    print(f"   {len(inter)} pairs with >=3 joint observations")
    for a, b, v, n in inter[:5]:
        print(f"     SYNERGY    {a:>7} + {b:<7} {v:+7.2f} REU  (n={n})")
    for a, b, v, n in inter[-3:]:
        print(f"     ANTAGONISM {a:>7} + {b:<7} {v:+7.2f} REU  (n={n})")
    json.dump(dict(auc_raw=a_raw, auc_matched=a_m,
                   marginals={k: v[0] for k, v in marg.items()},
                   interactions=[(a, b, v, n) for a, b, v, n in inter]),
              open(os.path.join(OUT, f"{klass}.json"), "w"), indent=1)
    print(f"\nwrote {OUT}/{klass}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
