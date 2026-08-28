#!/usr/bin/env python
r"""
155_ddg_false_negatives.py -- if we use ddG to SHRINK a library, what do we lose?

Jannis's question: what is the rate at which ddG is right to eliminate something,
and are there real sensors carrying a substitution ddG would have thrown away?

That is the false-negative rate of the filter, and it is the number that decides
whether the filter is usable at all. Every previous report of this filter (§69,
§82d) quoted ENRICHMENT -- how much better than chance the survivors are. That
is the wrong side of the ledger for library design: a library you never build
cannot be rescued by the fact that what you did build was enriched. Recall is
the constraint; enrichment is the reward. This script reports recall first.

WHAT IS AVAILABLE. results/viability_relax/ holds 670 scored variants: 70 REAL
coumarin sensors (sd07), 300 LIBRARY draws from Tian's own 138,240-member
coumarin library, and 300 WILD draws from the full DSM-Hao menu at the same 11
positions. All are COMBINATIONS of 3-10 substitutions; there are no singles, so
the per-substitution question ("was this residue nominated?") cannot be answered
directly from this set. It is answered two ways instead:

  1. at the combination level -- at a threshold retaining X % of real sensors,
     how much of the library survives, i.e. how much shrinkage is bought
  2. by ATTRIBUTION -- which substitutions are over-represented among the real
     sensors the filter eliminates, compared with those it keeps. A substitution
     that is systematically enriched in the rejected sensors is one ddG is
     penalising, and that is the closest this data can come to naming a residue
     the filter would wrongly remove.

⚠ n_sub IS A CONFOUND AND IS CONTROLLED. ddG grows with the number of
substitutions, and REAL and LIBRARY do not have identical n_sub distributions.
Every comparison below is computed WITHIN n_sub strata and pooled, never across.

⚠ LIBRARY IS NOT A CLEAN NEGATIVE (§122 header): an unscreened library member may
be a perfectly good sensor. So "library retention" is shrinkage, not specificity,
and no precision number is quoted anywhere here.
"""
import glob
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "viability_relax")


def load():
    rows = []
    for f in glob.glob(os.path.join(OUT, "chunk_*.json")):
        d = json.load(open(f))
        k = [x for x in d if x != "nrep"][0]
        rows += d[k]
    return rows


def strata_threshold(real, keep_frac):
    """Per-n_sub ddG cutoff retaining `keep_frac` of REAL sensors in that stratum."""
    th = {}
    by = defaultdict(list)
    for r in real:
        by[r["n_sub"]].append(r["ddG"])
    for n, v in by.items():
        th[n] = float(np.quantile(v, keep_frac)) if v else float("inf")
    return th


def apply_th(rows, th):
    kept, lost = [], []
    for r in rows:
        t = th.get(r["n_sub"])
        if t is None:
            continue
        (kept if r["ddG"] <= t else lost).append(r)
    return kept, lost


def main():
    rows = load()
    real = [r for r in rows if r["set"] == "REAL"]
    lib = [r for r in rows if r["set"] == "LIBRARY"]
    wild = [r for r in rows if r["set"] == "WILD"]
    print(f"{len(real)} REAL sensors, {len(lib)} LIBRARY draws, {len(wild)} WILD; "
          f"n_sub {min(r['n_sub'] for r in rows)}-{max(r['n_sub'] for r in rows)}")
    ns = sorted({r["n_sub"] for r in real})
    print(f"n_sub strata present in REAL: {ns}")

    print("\n" + "=" * 82)
    print("HOW MUCH SHRINKAGE DOES A GIVEN SENSOR RETENTION BUY?")
    print("   (thresholds set WITHIN n_sub strata, so size cannot leak in)")
    print("=" * 82)
    print(f"{'sensor retention':>18}{'library kept':>14}{'wild kept':>12}"
          f"{'shrinkage':>12}{'enrichment':>12}")
    res = {}
    for kf in (1.00, 0.95, 0.90, 0.80, 0.70, 0.50):
        th = strata_threshold(real, kf)
        rk, _ = apply_th(real, th)
        lk, _ = apply_th(lib, th)
        wk, _ = apply_th(wild, th)
        lr = len(lk) / max(len(lib), 1)
        rr = len(rk) / max(len(real), 1)
        print(f"{rr:>17.0%}{lr:>14.0%}{len(wk)/max(len(wild),1):>12.0%}"
              f"{1/max(lr,1e-9):>11.1f}x{rr/max(lr,1e-9):>11.2f}x")
        res[f"{kf:.2f}"] = dict(real=rr, lib=lr, wild=len(wk) / max(len(wild), 1))

    print("\n   'shrinkage' is how many-fold smaller the library becomes.")
    print("   'enrichment' is sensor retention / library retention -- the number")
    print("   previously quoted. Read the FIRST column as the cost.")

    # ---- attribution: what is in the sensors we lose? ------------------
    print("\n" + "=" * 82)
    print("WHICH SUBSTITUTIONS ARE OVER-REPRESENTED IN THE SENSORS ddG REJECTS?")
    print("=" * 82)
    for kf in (0.90, 0.80):
        th = strata_threshold(real, kf)
        kept, lost = apply_th(real, th)
        ck = Counter(s for r in kept for s in r["subs"])
        cl = Counter(s for r in lost for s in r["subs"])
        nk = sum(ck.values()) or 1
        nl = sum(cl.values()) or 1
        print(f"\n   at {len(kept)}/{len(real)} sensors retained "
              f"({len(lost)} eliminated)")
        rows_ = []
        for s in set(ck) | set(cl):
            fk, fl = ck[s] / nk, cl[s] / nl
            if cl[s] >= 2:
                rows_.append((fl - fk, s, cl[s], ck[s]))
        rows_.sort(reverse=True)
        if not rows_:
            print("      no substitution appears in 2+ eliminated sensors")
            continue
        print(f"      {'sub':>8}{'in lost':>9}{'in kept':>9}{'excess':>9}")
        for d, s, a, b in rows_[:12]:
            print(f"      {s:>8}{a:>9}{b:>9}{d:>+9.3f}")
        never = sorted(s for s in cl if s not in ck)
        if never:
            print(f"      substitutions appearing ONLY in eliminated sensors: "
                  f"{', '.join(never)}")

    # ---- the honest ceiling: sensors lost at ANY useful threshold ------
    th = strata_threshold(real, 0.90)
    _, lost = apply_th(real, th)
    print(f"\n   Example eliminated sensors at 90 % retention:")
    for r in sorted(lost, key=lambda r: -r["ddG"])[:6]:
        print(f"      {r['id']:<14} n_sub {r['n_sub']:>2}  ddG {r['ddG']:+7.2f}  "
              f"{'+'.join(r['subs'])}")
    # ---- a mechanism I proposed and this REFUTES -----------------------
    # Hypothesis: ddG is a clash detector (§79), so it should penalise sensors
    # that GROW side chains into the pocket -- which is precisely what pocket
    # expansion needs. Tested with Zamyatnin side-chain volumes, controlling
    # n_sub. It does not hold: 69 of 70 real sensors net-grow, so there is no
    # contrast to explain the rejections, and the partial correlation is +0.19.
    VOL = dict(zip("AGVLIPFMWSTCYNQDEKRH",
                   [88.6, 60.1, 140.0, 166.7, 166.7, 112.7, 189.9, 162.9, 227.8,
                    89.0, 116.1, 108.5, 193.6, 114.1, 143.8, 111.1, 138.4,
                    168.6, 173.4, 153.2]))
    import re as _re

    def netvol(subs):
        t = 0.0
        for s in subs:
            m = _re.fullmatch(r"([A-Z])(\d+)([A-Z])", s)
            if m:
                t += VOL[m.group(3)] - VOL[m.group(1)]
        return t

    dv = np.array([netvol(r["subs"]) for r in real])
    gg = np.array([r["ddG"] for r in real])
    nn = np.array([r["n_sub"] for r in real], float)

    def _r(a, b):
        a, b = a - a.mean(), b - b.mean()
        return float((a * b).sum() / np.sqrt((a ** 2).sum() * (b ** 2).sum()))

    def _res(y, x):
        X = np.c_[np.ones(len(x)), x]
        return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]

    print("\n" + "=" * 82)
    print("MECHANISM TEST -- does ddG reject sensors that GROW the pocket?  NO")
    print("=" * 82)
    print(f"   ddG vs net side-chain volume change   r = {_r(dv, gg):+.3f}")
    print(f"   same, controlling n_sub               r = {_r(_res(dv, nn), _res(gg, nn)):+.3f}")
    print(f"   sensors that net-grow: {int((dv > 0).sum())}/{len(real)} -- there is")
    print("   no shrink group to contrast against, so this hypothesis of mine is")
    print("   not supported and the substitution attribution above stands on its")
    print("   own small counts (2-7 per substitution). Treat it as a lead only.")
    res["volume_partial_r"] = _r(_res(dv, nn), _res(gg, nn))

    json.dump(res, open(os.path.join(ROOT, "results", "library_recall",
                                     "ddg_false_negatives.json"), "w"), indent=1)
    print("\n⚠ These are COMBINATION-level rejections. Whether a specific RESIDUE")
    print("  would be eliminated needs single-substitution ddG, which job 27925009")
    print("  is computing for Beltran's 20 positions.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
