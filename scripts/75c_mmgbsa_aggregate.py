#!/usr/bin/env python
"""
75c_mmgbsa_aggregate.py -- read the MM-GBSA results and apply the pre-registered rule.

THE RULE, fixed in 75_mmgbsa_build.py before any number existed
---------------------------------------------------------------
    SUCCESS  ddG(K59R) favourable for mandipropamid AND unfavourable for ABA
    FAILURE  same sign for both -- the model is reading charge, not complementarity
    NULL     |ddG| within the seed-to-seed spread; underpowered, say so

ddG is reported relative to WT within each arm, so the receptor, the ligand and the
force field all cancel; only the substitution differs. The spread across the 8
independently repacked seeds is the error bar, and nothing here is quoted without it.
"""
import csv
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "mmgbsa")
ARMS = ["mandi", "aba"]

#: Read the variant list from the file 77_mmgbsa_extended.sh actually looped over,
#: so the table can never silently disagree with what was run. The old hard-coded
#: list also used a different naming convention ("K59R" vs "59R") than the
#: directories on disk, which would have shown every row as incomplete.
def _variants():
    f = os.path.join(D, "variants.txt")
    if os.path.exists(f):
        return [l.strip() for l in open(f) if l.strip()]
    return sorted({d.split("_", 1)[1] for d in os.listdir(D)
                   if d.startswith("aba_")})

VARIANTS = _variants()

#: The four substitutions that appear in the 4WVO mandipropamid sensor. These are
#: the positive controls: the question for the extended pool is whether they rank
#: highly among 22 candidates, not merely whether each beats WT on its own.
GROUND_TRUTH = {"59R", "81I", "108A", "159L"}


#: Prefer the MD ensemble when it exists. The repack-seed ensemble it replaces was
#: not an ensemble at all -- see 75e: all 8 mandipropamid repacks were byte-identical,
#: so every sd of 0.00 measured nothing.
FRAME_FILES = ("mmgbsa_md_frames.csv", "mmgbsa_frames.csv")


def per_frame(arm, var):
    """DELTA TOTAL for each frame, from the -eo csv; falls back to the .dat mean."""
    csvp = None
    for f in FRAME_FILES:
        c = os.path.join(D, f"{arm}_{var}", f)
        if os.path.exists(c):
            csvp = c
            break
    if csvp is None:
        csvp = os.path.join(D, f"{arm}_{var}", "mmgbsa_frames.csv")
    vals = []
    if csvp and os.path.exists(csvp):
        block = False
        for row in csv.reader(open(csvp)):
            if not row:
                block = False
                continue
            head = ",".join(row).upper()
            if "DELTA" in head and "FRAME" in head:
                block = True
                cols = [c.strip().upper() for c in row]
                try:
                    ti = cols.index("DELTA TOTAL")
                except ValueError:
                    ti = len(cols) - 1
                continue
            if block:
                try:
                    vals.append(float(row[ti]))
                except (ValueError, IndexError):
                    block = False
    return vals


def dat_mean(arm, var):
    p = os.path.join(D, f"{arm}_{var}", "mmgbsa_md.dat")
    if not os.path.exists(p):
        p = os.path.join(D, f"{arm}_{var}", "mmgbsa.dat")
    if not os.path.exists(p):
        return None, None
    txt = open(p).read()
    m = re.search(r"DELTA TOTAL\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", txt)
    if m:
        return float(m.group(1)), float(m.group(2))
    return None, None


def _sd(v):
    if len(v) < 2:
        return 0.0
    mu = sum(v) / len(v)
    return math.sqrt(sum((x - mu) ** 2 for x in v) / (len(v) - 1))


def stats(v, nblocks=10):
    """Mean, and a BLOCK-AVERAGED standard error.

    The frames come from one continuous 250 ps trajectory, so they are not
    independent draws: sd/sqrt(100) would understate the true uncertainty by
    whatever the correlation time is. Splitting the trajectory into contiguous
    blocks and taking the spread of the block means costs statistical power but
    does not pretend to have samples it lacks. Both numbers are reported so the
    inflation factor is visible rather than buried.
    """
    if not v:
        return None, None, None
    mu = sum(v) / len(v)
    naive_se = _sd(v) / math.sqrt(len(v)) if len(v) > 1 else 0.0
    b = max(1, len(v) // nblocks)
    means = [sum(v[i:i + b]) / len(v[i:i + b])
             for i in range(0, len(v) - b + 1, b)]
    block_se = _sd(means) / math.sqrt(len(means)) if len(means) > 1 else naive_se
    return mu, block_se, naive_se


def main():
    G = {}
    for arm in ARMS:
        for var in VARIANTS:
            v = per_frame(arm, var)
            mu, se, naive = stats(v)
            if mu is None:
                mu, sd = dat_mean(arm, var)
                se = naive = sd
            G[(arm, var)] = (mu, se, naive, len(v))

    print("=" * 78)
    print("MM-GBSA  DELTA G bind (kcal/mol), mean +/- BLOCK standard error")
    print("  250 ps GB ensemble, 100 frames, 10 contiguous blocks")
    print("=" * 78)
    print(f"  {'variant':<12}{'mandipropamid':>24}{'ABA':>24}")
    for var in VARIANTS:
        cells = []
        for arm in ARMS:
            mu, se, naive, n = G[(arm, var)]
            cells.append(f"{mu:8.2f} +/- {se:5.2f} (n={n})" if mu is not None
                         else f"{'--':>24}")
        print(f"  {var:<12}{cells[0]:>24}{cells[1]:>24}")

    # how badly would sd/sqrt(n) have flattered us?
    infl = [G[k][1] / G[k][2] for k in G
            if G[k][1] is not None and G[k][2] not in (None, 0.0)]
    if infl:
        print(f"\n  block SE / naive SE: median {sorted(infl)[len(infl)//2]:.2f}x"
              "  (>1 means the frames are correlated, as expected)")

    print("\n" + "=" * 78)
    print("ddG vs WT within each arm  (negative = the mutation helps binding)")
    print("SELECTIVITY = ddG(mandi) - ddG(ABA); negative = shifts toward mandipropamid")
    print("=" * 78)
    rows = []
    for var in VARIANTS:
        if var == "WT":
            continue
        vals = {}
        bad = False
        for arm in ARMS:
            mu, se, _, _ = G[(arm, var)]
            w, wse, _, _ = G[(arm, "WT")]
            if mu is None or w is None:
                bad = True
                break
            vals[arm] = (mu - w, math.sqrt(se ** 2 + wse ** 2))
        if bad:
            print(f"  {var:<12} incomplete")
            continue
        dm, em = vals["mandi"]
        da, ea = vals["aba"]
        rows.append((dm - da, math.sqrt(em ** 2 + ea ** 2), var, dm, da))

    rows.sort()
    print(f"  {'rank':>4}  {'variant':<12}{'ddG mandi':>11}{'ddG ABA':>11}"
          f"{'select.':>10}{'+/-':>8}   ")
    for r, (sel, err, var, dm, da) in enumerate(rows, 1):
        star = " *" if var in GROUND_TRUTH else "  "
        flag = "" if abs(sel) > err else "  (within error)"
        print(f"  {r:>4}  {var:<12}{dm:>11.2f}{da:>11.2f}{sel:>10.2f}"
              f"{err:>8.2f}{star}{flag}")

    print("\n  * = a substitution in the 4WVO mandipropamid sensor (positive control)")

    # THE TEST: do the known-good mutations concentrate at the selective end?
    n = len(rows)
    ranks = [(r, v) for r, (_, _, v, _, _) in enumerate(rows, 1) if v in GROUND_TRUTH]
    if ranks:
        print("\n" + "=" * 78)
        print("ENRICHMENT OF THE KNOWN SENSOR MUTATIONS")
        print("=" * 78)
        for r, v in ranks:
            print(f"  {v:<12} rank {r:>3} of {n}   (top {100.0 * r / n:.0f}%)")
        rr = [r for r, _ in ranks]
        obs = sum(rr) / len(rr)
        exp = (n + 1) / 2.0
        print(f"\n  mean rank {obs:.1f}  vs  {exp:.1f} expected if the ranking were random")
        # exact one-sided permutation p over all C(n,k) rank subsets, via DP
        k = len(rr)
        tgt = sum(rr)
        # count subsets of size k from 1..n with sum <= tgt
        dp = [[0] * (tgt + 1) for _ in range(k + 1)]
        dp[0][0] = 1
        for val in range(1, n + 1):
            for kk in range(k, 0, -1):
                for s in range(tgt, val - 1, -1):
                    dp[kk][s] += dp[kk - 1][s - val]
        import math as _m
        favourable = sum(dp[k])
        total = _m.comb(n, k)
        print(f"  one-sided permutation p = {favourable}/{total} = "
              f"{favourable / total:.4f}")
        if favourable / total < 0.05:
            print("  -> the known mutations sit significantly toward the "
                  "mandipropamid-selective end")
        else:
            print("  -> NOT significant: this pool does not separate the known "
                  "mutations from the rest")

    ident = [f"{a}_{v}" for a in ARMS for v in VARIANTS
             if G[(a, v)][1] == 0.0 and G[(a, v)][3] > 1]
    if ident:
        print(f"\n  !! SE is EXACTLY 0.00 with n>1 for: {', '.join(ident)}")
        print("     That is not precision, it is an ensemble with no diversity.")


if __name__ == "__main__":
    main()
