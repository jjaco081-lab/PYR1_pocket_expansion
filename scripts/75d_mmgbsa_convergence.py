#!/usr/bin/env python
"""
75d_mmgbsa_convergence.py -- is the MM-GBSA number converged enough to rank on?

Two independent checks, both of which the 50 ps run could not have passed:

  (1) BETWEEN RUNS. The same seven variants were scored at 50 ps (job 27655xxx,
      archived in data/mmgbsa_50ps/) and again at 250 ps (job 27676964). If the
      method were converged the two should agree.

  (2) WITHIN A RUN. Split each 250 ps ensemble into its first and second half.
      A converged average does not care which half you took it from.

Why this matters more than the error bars: the block standard error in 75c is a
WITHIN-run quantity. It measures scatter about the current mean, so it stays
small even while that mean is systematically marching somewhere else. Only a
drift check can see the march.

The stake: every ddG in an arm subtracts that arm's WT. A reference that is
still moving shifts every ddG in the arm by the same amount -- which flips
SIGN-based verdicts while leaving RANK-based ones untouched. The two conclusions
therefore have to be judged separately, and this script reports both.
"""
import csv
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NEW = os.path.join(ROOT, "data", "mmgbsa")
OLD = os.path.join(ROOT, "data", "mmgbsa_50ps")
ARMS = ("mandi", "aba")

#: the 50 ps run used spelled-out names, the extended run uses bare positions
MAP = {"WT": "WT", "K59R": "59R", "K59Q": "59Q", "K59N": "59N",
       "V81I": "81I", "F108A": "108A", "F159L": "159L"}


def dat(path):
    if not os.path.exists(path):
        return None
    m = re.search(r"DELTA TOTAL\s+(-?\d+\.\d+)", open(path).read())
    return float(m.group(1)) if m else None


def old_val(arm, v):
    if v == "WT":
        return dat(os.path.join(OLD, f"{arm}_WT_50ps.dat"))
    for f in ("mmgbsa_md.dat", "mmgbsa.dat"):
        r = dat(os.path.join(OLD, f"{arm}_{v}", f))
        if r is not None:
            return r
    return None


def frames(arm, var):
    p = os.path.join(NEW, f"{arm}_{var}", "mmgbsa_md_frames.csv")
    if not os.path.exists(p):
        return []
    vals, block, ti = [], False, None
    for row in csv.reader(open(p)):
        if not row:
            block = False
            continue
        head = ",".join(row).upper()
        if "DELTA" in head and "FRAME" in head:
            block = True
            cols = [c.strip().upper() for c in row]
            ti = cols.index("DELTA TOTAL") if "DELTA TOTAL" in cols else len(cols) - 1
            continue
        if block:
            try:
                vals.append(float(row[ti]))
            except (ValueError, IndexError):
                block = False
    return vals


def main():
    out = []
    P = lambda s="": (print(s), out.append(s))

    P("=" * 78)
    P("(1) BETWEEN RUNS: 50 ps -> 250 ps")
    P("=" * 78)
    P(f"  {'variant':<8}{'arm':<7}{'50ps':>9}{'250ps':>9}{'shift':>9}")
    sh = {}
    for v in MAP:
        for arm in ARMS:
            o, n = old_val(arm, v), dat(os.path.join(NEW, f"{arm}_{MAP[v]}", "mmgbsa_md.dat"))
            if o is None or n is None:
                continue
            sh[(arm, v)] = (o, n)
            P(f"  {v:<8}{arm:<7}{o:>9.2f}{n:>9.2f}{n - o:>9.2f}")

    P("")
    P("  The pre-registered K59R selectivity test, computed both ways:")
    for tag, i in (("50 ps", 0), ("250 ps", 1)):
        dm = sh[("mandi", "K59R")][i] - sh[("mandi", "WT")][i]
        da = sh[("aba", "K59R")][i] - sh[("aba", "WT")][i]
        P(f"    {tag:<7} ddG_mandi {dm:+7.2f}  ddG_ABA {da:+7.2f}  "
          f"selectivity {dm - da:+7.2f}  "
          f"{'PASS (prefers mandipropamid)' if dm - da < 0 else 'FAIL (prefers ABA)'}")
    big = max(sh.items(), key=lambda kv: abs(kv[1][1] - kv[1][0]))
    P(f"\n  largest single shift: {big[0][0]}_{big[0][1]} "
      f"{big[1][1] - big[1][0]:+.2f} kcal/mol -- and because it is a WT, it is the "
      f"reference\n  every other ddG in that arm is measured against.")

    P("")
    P("=" * 78)
    P("(2) WITHIN A RUN: first 125 ps vs last 125 ps")
    P("=" * 78)
    vs = [l.strip() for l in open(os.path.join(NEW, "variants.txt")) if l.strip()]
    dr = []
    for v in vs:
        for arm in ARMS:
            f = frames(arm, v)
            if len(f) < 20:
                continue
            h = len(f) // 2
            a, b = sum(f[:h]) / h, sum(f[h:]) / len(f[h:])
            dr.append((abs(b - a), arm, v, a, b))
    dr.sort(reverse=True)
    P(f"  {'variant':<12}{'arm':<7}{'1st':>9}{'2nd':>9}{'drift':>9}")
    for d, arm, v, a, b in dr[:8]:
        P(f"  {v:<12}{arm:<7}{a:>9.2f}{b:>9.2f}{b - a:>9.2f}")
    med = sorted(x[0] for x in dr)[len(dr) // 2]
    n1 = sum(1 for x in dr if x[0] > 1.0)
    P(f"\n  median |drift| {med:.2f} kcal/mol; {n1}/{len(dr)} runs drift "
      f"more than 1.0 kcal/mol")
    for d, arm, v, a, b in dr:
        if v == "WT":
            P(f"  WT {arm:<6} drift {b - a:+.2f}")

    P("")
    P("=" * 78)
    P("VERDICT")
    P("=" * 78)
    P("  The 250 ps ensemble is NOT converged. The WT-ABA reference moved +7.75")
    P("  kcal/mol between 50 and 250 ps and is still moving +2.62 within the 250 ps")
    P("  window, in the same direction. It has not settled, it is still climbing.")
    P("")
    P("  Consequences, kept separate because they are not equally damaged:")
    P("    SIGN-based conclusions are void. The 50 ps 'K59 flip' was produced by")
    P("      an unconverged WT-ABA reference, not by the K59R substitution. At")
    P("      250 ps the same test gives the opposite answer.")
    P("    RANK-based conclusions survive a uniform reference error, since a")
    P("      constant added to every ddG in an arm cannot reorder them. They are")
    P("      still limited by the per-variant drift, which at a median 1.05")
    P("      kcal/mol is the same size as the spacing at the top of the ranking.")
    open(os.path.join(ROOT, "results", "mmgbsa", "convergence.txt"), "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
