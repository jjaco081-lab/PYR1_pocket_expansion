#!/usr/bin/env python
"""
84_explicit_aggregate.py -- did explicit solvent stabilise the reference?

Three questions, in the order that matters:

 1. WITHIN-run drift. First half vs second half of each 10 ns trajectory. The
    implicit runs failed here: WT-ABA moved +2.62 kcal/mol inside its own 250 ps
    window, in the same direction it had moved between 50 and 250 ps.

 2. BETWEEN-seed spread. Three independent velocity seeds. This is the honest
    error bar, and the one thing neither the 8 byte-identical repacks nor the
    250 ps block SE could ever have provided. A within-run error bar cannot see
    a drift it is centred on; a between-run one can.

 3. Does it agree with either implicit result? The 50 ps and 250 ps numbers
    disagree by 7.75 kcal/mol on WT-ABA. Explicit solvent is the tie-breaker.
"""
import csv
import glob
import math
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "mmgbsa_explicit")


def frames(d):
    p = os.path.join(d, "mmgbsa_exp_frames.csv")
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


def sd(v):
    if len(v) < 2:
        return 0.0
    m = sum(v) / len(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


IMPLICIT = {"aba": {"50ps": -32.74, "250ps": -25.00},
            "mandi": {"50ps": -36.33, "250ps": -38.00}}

print("=" * 78)
print("EXPLICIT SOLVENT (ff19SB/OPC/0.15 M KCl), 10 ns x 3 seeds, scored igb=8")
print("=" * 78)
print(f"  {'system':<14}{'mean':>9}{'1st half':>10}{'2nd half':>10}{'drift':>9}{'n':>6}")
arms = {}
for d in sorted(glob.glob(os.path.join(D, "*/"))):
    u = os.path.basename(d.rstrip("/"))
    v = frames(d)
    if not v:
        print(f"  {u:<14} no frames")
        continue
    h = len(v) // 2
    a, b = sum(v[:h]) / h, sum(v[h:]) / (len(v) - h)
    m = sum(v) / len(v)
    arms.setdefault(u.split("_")[0], []).append(m)
    print(f"  {u:<14}{m:>9.2f}{a:>10.2f}{b:>10.2f}{b - a:>9.2f}{len(v):>6}")

print()
print("=" * 78)
print("BETWEEN-SEED SPREAD  -- the honest error bar")
print("=" * 78)
for arm, ms in sorted(arms.items()):
    m = sum(ms) / len(ms)
    s = sd(ms)
    print(f"  {arm:<8} mean {m:>8.2f}   between-seed sd {s:>5.2f}   "
          f"seeds: {', '.join(f'{x:.2f}' for x in ms)}")

print()
print("=" * 78)
print("AGREEMENT WITH THE IMPLICIT RUNS")
print("=" * 78)
print(f"  {'arm':<8}{'50 ps':>9}{'250 ps':>9}{'explicit':>11}{'closer to':>12}")
for arm, ms in sorted(arms.items()):
    m = sum(ms) / len(ms)
    i50, i250 = IMPLICIT[arm]["50ps"], IMPLICIT[arm]["250ps"]
    closer = "50 ps" if abs(m - i50) < abs(m - i250) else "250 ps"
    print(f"  {arm:<8}{i50:>9.2f}{i250:>9.2f}{m:>11.2f}{closer:>12}")

print()
print("=" * 78)
print("WHAT THIS LICENSES")
print("=" * 78)
aba_sd = sd(arms.get("aba", [0, 0]))
print(f"  The between-seed sd on the ABA reference is {aba_sd:.2f} kcal/mol.")
print("  The selectivity spacing it would have to resolve at the top of the")
print("  22-variant ranking (README 40c) is 0.01-1.8 kcal/mol.")
if aba_sd > 1.0:
    print("  -> Still larger than the differences being ranked. Explicit solvent")
    print("     fixes the POSE, but a single trajectory per variant cannot resolve")
    print("     this ranking; replicates would be required for every variant.")
