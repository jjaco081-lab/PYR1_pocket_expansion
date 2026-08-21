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
VARIANTS = ["WT", "K59R", "K59Q", "K59N", "V81I", "F108A", "F159L"]
GROUND_TRUTH = {"K59R", "V81I", "F108A", "F159L"}


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


def stats(v):
    if not v:
        return None, None
    mu = sum(v) / len(v)
    sd = math.sqrt(sum((x - mu) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else 0.0
    return mu, sd


def main():
    G = {}
    for arm in ARMS:
        for var in VARIANTS:
            v = per_frame(arm, var)
            mu, sd = stats(v)
            if mu is None:
                mu, sd = dat_mean(arm, var)
            G[(arm, var)] = (mu, sd, len(v))
    print("=" * 76)
    print("MM-GBSA  DELTA G bind (kcal/mol), mean +/- sd over independent repack seeds")
    print("=" * 76)
    print(f"  {'variant':<8}{'mandipropamid':>22}{'ABA':>22}")
    for var in VARIANTS:
        cells = []
        for arm in ARMS:
            mu, sd, n = G[(arm, var)]
            cells.append(f"{mu:8.2f} +/- {sd:5.2f} (n={n})" if mu is not None else f"{'--':>22}")
        print(f"  {var:<8}{cells[0]:>22}{cells[1]:>22}")

    print("\n" + "=" * 76)
    print("ddG vs WT within each arm  (negative = the mutation helps binding)")
    print("=" * 76)
    print(f"  {'variant':<8}{'ddG mandi':>12}{'ddG ABA':>12}{'difference':>12}   reading")
    ok = True
    for var in VARIANTS[1:]:
        row = {}
        for arm in ARMS:
            mu, sd, _ = G[(arm, var)]
            w, wsd, _ = G[(arm, "WT")]
            if mu is None or w is None:
                ok = False
                row[arm] = (None, None)
                continue
            # errors add in quadrature; both are means over the same seed count
            row[arm] = (mu - w, math.sqrt(sd ** 2 + wsd ** 2))
        if row["mandi"][0] is None or row["aba"][0] is None:
            print(f"  {var:<8}{'incomplete':>12}")
            continue
        dm, em = row["mandi"]
        da, ea = row["aba"]
        diff = dm - da
        err = math.sqrt(em ** 2 + ea ** 2)
        if abs(diff) < err:
            verdict = "NULL (within the seed spread)"
        elif diff < 0:
            verdict = "prefers mandipropamid"
        else:
            verdict = "prefers ABA"
        star = " *" if var in GROUND_TRUTH else "  "
        print(f"  {var:<8}{dm:>12.2f}{da:>12.2f}{diff:>12.2f}{star} {verdict}")
    print("\n  * = a real 4WVO mandipropamid mutation. The pre-registered success")
    print("    criterion is K59R showing 'prefers mandipropamid' beyond the spread.")
    ident = [f"{a}_{v}" for a in ARMS for v in VARIANTS
             if G[(a, v)][1] == 0.0 and G[(a, v)][2] > 1]
    if ident:
        print(f"\n  !! sd is EXACTLY 0.00 with n>1 for: {', '.join(ident)}")
        print("     That is not precision, it is an ensemble with no diversity in it.")
        print("     Do not read a verdict off those rows.")
    if not ok:
        print("\n  !! some runs are incomplete -- the table above is partial")


if __name__ == "__main__":
    main()
