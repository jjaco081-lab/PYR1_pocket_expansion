#!/usr/bin/env python
r"""
199_cascade_report.py -- read cascade2.csv, print the funnel per arm and the
pre-registered arm comparisons.

The comparisons are fixed HERE, before the data exists, so they are not chosen
after seeing the answer:

  A. 1slack  vs 4rep      -- does slack in link3/link4 restore connectivity?
                             Same motif, same lead, same conditioning; only the
                             two linker ranges that carried 6-7 of the breaks
                             differ. This is the intended fix for fragmentation.
  B. 2noanch vs 4rep      -- is the ANCHOR or merely the SHORT LEAD responsible
                             for the fragmentation? 2noanch has the short lead
                             WITHOUT the A34-40 anchor, so it separates them.
  C. 3b3rasa vs batch03   -- what does partially_buried RASA do on batch 3's
                             known-good-connectivity contig (9/10 intact)?

Primary endpoint for A and B is the INTACT-CHAIN RATE (one connected piece),
because that is the defect Stage 0 exposed. Secondary endpoints: cavity rate on
the main piece, appendage length, and full-cascade PASS rate.

⚠ Stage 0's rates came from n=10 per arm, where 4/10 vs 0/10 sits at p = 0.087.
Nothing here is interpreted as established below n where Fisher can resolve it,
and every proportion is reported with a Wilson interval rather than bare.
"""
import csv, math, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CSV = os.path.join(ROOT, "results", "rfd3", "cascade2.csv")


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def fisher(a, b, c, d):
    """two-sided Fisher exact on [[a,b],[c,d]] without scipy."""
    def lc(n, k):
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    obs = lc(r1, a) + lc(n - r1, c1 - a) - lc(n, c1)
    tot = 0.0
    for i in range(max(0, c1 - (n - r1)), min(r1, c1) + 1):
        l = lc(r1, i) + lc(n - r1, c1 - i) - lc(n, c1)
        if l <= obs + 1e-9:
            tot += math.exp(l)
    return min(1.0, tot)


def main():
    if not os.path.exists(CSV):
        print(f"no {CSV}"); return 1
    rows = defaultdict(list)
    with open(CSV) as fh:
        for r in csv.DictReader(fh):
            rows[r["arm"]].append(r)

    def rate(rs, pred):
        k = sum(1 for r in rs if pred(r))
        lo, hi = wilson(k, len(rs))
        return k, len(rs), f"{k}/{len(rs)} ({100*k/max(len(rs),1):.0f}%, " \
                           f"95% CI {100*lo:.0f}-{100*hi:.0f})"

    intact = lambda r: r["n_pieces"] == "1"                          # noqa: E731
    usable = lambda r: r["motif_split"] == "0" and (                  # noqa: E731
        r["n_orphan"] == "0" or r["orphan_termini_only"] == "1")
    hascav = lambda r: r["cavity"] not in ("", "0", "0.0") and \
        float(r["cavity"] or 0) > 0                                  # noqa: E731
    passed = lambda r: r["verdict"] == "PASS"                        # noqa: E731

    print("=" * 78)
    print("CASCADE v2 FUNNEL")
    print("=" * 78)
    order = ["batch03", "arm0a_s2026", "arm0b_s2027", "arm4rep_s3004",
             "arm1slack_s3001", "arm2noanch_s3002", "arm3b3rasa_s3003"]
    for arm in [a for a in order if a in rows] + \
               [a for a in rows if a not in order]:
        rs = rows[arm]
        print(f"\n== {arm}   n={len(rs)}")
        print(f"   one intact chain     {rate(rs, intact)[2]}")
        print(f"   usable (truncatable) {rate(rs, usable)[2]}")
        print(f"   cavity on main piece {rate(rs, hascav)[2]}")
        print(f"   FULL CASCADE PASS    {rate(rs, passed)[2]}")
        fa = defaultdict(int)
        for r in rs:
            if r["verdict"] != "PASS":
                fa[r["fail_at"] or r["verdict"]] += 1
        for k, v in sorted(fa.items(), key=lambda x: -x[1]):
            print(f"      {v:>4}  {k}")
        ps = [r for r in rs if passed(r)]
        if ps:
            ps.sort(key=lambda r: -float(r["envelope"] or 0))
            print(f"   PASSING (by envelope):")
            for r in ps[:12]:
                print(f"      {r['design']:<6} env {r['envelope']:>7} "
                      f"cav {r['cavity']:>6} rmax {r['r_max']:>5} "
                      f"Rg {r['rg']:>5} app {r['appendage']:>2} "
                      f"orph {r['n_orphan']:>3} dev {r['dev_max']}")

    print("\n" + "=" * 78)
    print("PRE-REGISTERED COMPARISONS  (primary endpoint: intact-chain rate)")
    print("=" * 78)
    for lbl, x, y in [("A slack fix:      1slack vs 4rep",
                       "arm1slack_s3001", "arm4rep_s3004"),
                      ("B anchor vs lead: 2noanch vs 4rep",
                       "arm2noanch_s3002", "arm4rep_s3004"),
                      ("C RASA on b3:     3b3rasa vs batch03",
                       "arm3b3rasa_s3003", "batch03")]:
        if x not in rows or y not in rows:
            print(f"\n{lbl}: MISSING ({x if x not in rows else y})")
            continue
        print(f"\n{lbl}")
        for nm, pred in [("intact", intact), ("cavity", hascav),
                         ("PASS", passed)]:
            ka, na, sa = rate(rows[x], pred)
            kb, nb, sb = rate(rows[y], pred)
            p = fisher(ka, na - ka, kb, nb - kb)
            star = "  <-- significant" if p < 0.05 else ""
            print(f"   {nm:<7} {sa:<28} vs {sb:<28} p={p:.4f}{star}")
    print("\n  Wilson 95% intervals. Fisher two-sided, uncorrected -- three")
    print("  pre-registered comparisons, so read p<0.017 as Holm-significant.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
