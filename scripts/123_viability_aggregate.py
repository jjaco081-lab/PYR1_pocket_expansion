#!/usr/bin/env python
r"""
123_viability_aggregate.py -- did protein-only Rosetta separate real sensors from
random members of the same library?

READ THE TWO COMPARISONS DIFFERENTLY (see 122's header)
  REAL vs WILD     a sanity check. If this fails, the score has no signal at all.
  REAL vs LIBRARY  the question that matters. Tian's menu has already removed the
                   obviously broken combinations, so this asks whether Rosetta
                   adds anything ON TOP of what the wet lab already did.

⚠ LIBRARY is a contaminated negative set -- its members are unlabelled, not known
failures, and some are surely working sensors nobody screened. That biases the
comparison TOWARD a null, so a clear separation is strong evidence and a null is
weak evidence.

AUC is reported as the probability that a random REAL scores better (lower ddG)
than a random member of the null, with a Mann-Whitney U p-value.
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "viability")


def mannwhitney(a, b):
    """U statistic, AUC and a normal-approximation two-sided p (ties corrected)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    n1, n2 = len(a), len(b)
    allv = np.concatenate([a, b])
    order = np.argsort(allv)
    ranks = np.empty(len(allv), float)
    ranks[order] = np.arange(1, len(allv) + 1)
    # average ranks over ties
    _, inv, cnt = np.unique(allv, return_inverse=True, return_counts=True)
    for i, c in enumerate(cnt):
        if c > 1:
            m = inv == i
            ranks[m] = ranks[m].mean()
    R1 = ranks[:n1].sum()
    U1 = R1 - n1 * (n1 + 1) / 2
    auc = U1 / (n1 * n2)
    mu = n1 * n2 / 2
    tie = sum(c ** 3 - c for c in cnt)
    N = n1 + n2
    sd = np.sqrt(n1 * n2 / 12 * ((N + 1) - tie / (N * (N - 1))))
    z = (U1 - mu) / sd if sd > 0 else 0.0
    from math import erfc, sqrt
    p = erfc(abs(z) / sqrt(2))
    return auc, z, p


def main():
    files = sorted(glob.glob(os.path.join(OUT, "chunk_*.json")))
    if not files:
        raise SystemExit(f"no chunks in {OUT}")
    rows = []
    for f in files:
        rows += json.load(open(f))["rows"]
    by = defaultdict(list)
    for r in rows:
        by[r["set"]].append(r)
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    say("=" * 78)
    say("PROTEIN-ONLY VIABILITY: can Rosetta rank combinations inside a fixed menu?")
    say("=" * 78)
    say(f"   {len(files)} chunks, {len(rows)} variants scored")
    say("")
    say(f"   {'set':<10}{'n':>5}{'ddG median':>12}{'mean':>10}{'sd':>9}"
        f"{'n_sub':>8}{'shell':>8}")
    for k in ("REAL", "LIBRARY", "WILD"):
        v = [r["ddG"] for r in by[k]]
        if not v:
            continue
        say(f"   {k:<10}{len(v):>5}{np.median(v):>12.1f}{np.mean(v):>10.1f}"
            f"{np.std(v):>9.1f}"
            f"{np.mean([r['n_sub'] for r in by[k]]):>8.1f}"
            f"{np.mean([r['n_shell'] for r in by[k]]):>8.1f}")
    say("")
    say("   n_sub and shell are printed because they are the obvious confound:")
    say("   more mutations open a larger repacking shell. The random sets were")
    say("   drawn to match REAL's substitution-count range, so these columns")
    say("   should be close -- if they are not, the comparison is contaminated.")

    real = [r["ddG"] for r in by["REAL"]]
    say("")
    for null in ("WILD", "LIBRARY"):
        v = [r["ddG"] for r in by[null]]
        if not v or not real:
            continue
        auc, z, p = mannwhitney(real, v)
        say(f"   REAL vs {null:<8} AUC = {auc:.3f}   z = {z:+.2f}   p = {p:.2g}")
        say(f"      (AUC 0.5 = indistinguishable; > 0.5 means REAL scores BETTER)")
    say("")
    say("   How many of each set would survive a filter set at REAL's 90th"
        " percentile?")
    if real:
        thr = float(np.percentile(real, 90))
        say(f"   threshold ddG <= {thr:.1f} REU (keeps 90 % of real sensors)")
        for k in ("REAL", "LIBRARY", "WILD"):
            v = [r["ddG"] for r in by[k]]
            if v:
                say(f"     {k:<10}{100*np.mean(np.array(v) <= thr):5.1f} % survive"
                    f"   ({int(np.sum(np.array(v) <= thr))}/{len(v)})")
        lib = np.array([r["ddG"] for r in by["LIBRARY"]])
        if len(lib):
            keep = float(np.mean(lib <= thr))
            say("")
            say(f"   ==> such a filter would shrink Tian's 138,240-member library to")
            say(f"       about {138240*keep:,.0f} while keeping 90 % of the known")
            say(f"       sensors -- a {1/keep:.1f}x reduction, IF the score is real.")
    say("")
    say("   ⚠ Everything above is recall-only. LIBRARY members are UNLABELLED, not")
    say("   known failures, so the reduction figure is what the filter would")
    say("   REMOVE, not what it would remove correctly.")
    with open(os.path.join(OUT, "viability.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    say(f"\n   written to {OUT}/viability.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
