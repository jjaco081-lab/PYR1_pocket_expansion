#!/usr/bin/env python
r"""
126_cavity_aggregate.py -- verdict on §68e's cavity-match proposal.

THE PRE-REGISTERED PREDICTIONS (fixed in 125's header before the run)
  1. ddG on COUMARIN reproduces §67 :  AUC < 0.40
  2. ddG on PFAS FAILS              :  AUC >= 0.45
  3. cavity match works on BOTH     :  AUC < 0.40 in each arm

AUC convention throughout: P(a random REAL ranks ABOVE a random null) on a
lower-is-better score, so BELOW 0.5 means REAL wins. (§67a: this label was printed
backwards for two runs before it was caught.)
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "cavity_match")


def mannwhitney(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n1, n2 = len(a), len(b)
    allv = np.concatenate([a, b])
    order = np.argsort(allv)
    ranks = np.empty(len(allv), float)
    ranks[order] = np.arange(1, len(allv) + 1)
    uq, inv, cnt = np.unique(allv, return_inverse=True, return_counts=True)
    for i, c in enumerate(cnt):
        if c > 1:
            m = inv == i
            ranks[m] = ranks[m].mean()
    U1 = ranks[:n1].sum() - n1 * (n1 + 1) / 2
    auc = U1 / (n1 * n2)
    tie = sum(c ** 3 - c for c in cnt)
    N = n1 + n2
    sd = np.sqrt(n1 * n2 / 12 * ((N + 1) - tie / (N * (N - 1))))
    z = (U1 - n1 * n2 / 2) / sd if sd > 0 else 0.0
    from math import erfc, sqrt
    return auc, z, erfc(abs(z) / sqrt(2))


def main():
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    data = {}
    for klass in ("coumarin", "pfas"):
        files = sorted(glob.glob(os.path.join(OUT, f"{klass}_*.json")))
        if not files:
            continue
        rows, meta = [], None
        for f in files:
            d = json.load(open(f))
            rows += d["rows"]
            meta = d
        by = defaultdict(list)
        for r in rows:
            by[r["set"]].append(r)
        data[klass] = (by, meta)

    if not data:
        raise SystemExit(f"no results in {OUT}")

    say("=" * 78)
    say("CAVITY MATCH vs ddG -- the §68e test")
    say("=" * 78)
    for klass, (by, meta) in data.items():
        say(f"\n   {klass.upper()}: target ligand volume {meta['target_vol']:.1f} A^3, "
            f"unmutated closed PYR1 cavity {meta['wt_cavity']:.1f} A^3")
        say(f"   {'set':<10}{'n':>5}{'ddG med':>10}{'cavity med':>12}"
            f"{'|cav-target| med':>18}{'n_sub':>7}")
        for k in ("REAL", "LIBRARY", "WILD"):
            v = by.get(k, [])
            if not v:
                continue
            say(f"   {k:<10}{len(v):>5}"
                f"{np.median([r['ddG'] for r in v]):>10.2f}"
                f"{np.median([r['cav'] for r in v]):>12.1f}"
                f"{np.median([r['cav_err'] for r in v]):>18.1f}"
                f"{np.mean([r['n_sub'] for r in v]):>7.1f}")

    say("")
    say("=" * 78)
    say("AUC  (below 0.5 = REAL wins; REAL vs LIBRARY is the one that matters)")
    say("=" * 78)
    say(f"   {'class':<10}{'score':<14}{'vs LIBRARY':>22}{'vs WILD':>22}")
    res = {}
    for klass, (by, meta) in data.items():
        for lbl, key in (("ddG", "ddG"), ("cavity match", "cav_err")):
            cells = []
            for null in ("LIBRARY", "WILD"):
                if not by.get("REAL") or not by.get(null):
                    cells.append("--")
                    continue
                auc, z, p = mannwhitney([r[key] for r in by["REAL"]],
                                        [r[key] for r in by[null]])
                res[(klass, lbl, null)] = auc
                cells.append(f"{auc:.3f} (p {p:.1g})")
            say(f"   {klass:<10}{lbl:<14}{cells[0]:>22}{cells[1]:>22}")

    say("")
    say("=" * 78)
    say("PRE-REGISTERED VERDICT")
    say("=" * 78)
    checks = [
        ("1. ddG works on coumarin (AUC < 0.40)",
         res.get(("coumarin", "ddG", "LIBRARY")), lambda x: x < 0.40),
        ("2. ddG FAILS on pfas (AUC >= 0.45)",
         res.get(("pfas", "ddG", "LIBRARY")), lambda x: x >= 0.45),
        ("3a. cavity match works on coumarin (AUC < 0.40)",
         res.get(("coumarin", "cavity match", "LIBRARY")), lambda x: x < 0.40),
        ("3b. cavity match works on pfas (AUC < 0.40)",
         res.get(("pfas", "cavity match", "LIBRARY")), lambda x: x < 0.40),
    ]
    ok = {}
    for name, val, test in checks:
        if val is None:
            say(f"   {name:<52} -- no data")
            continue
        good = test(val)
        ok[name[:2]] = good
        say(f"   {name:<52} AUC {val:.3f}   {'HOLDS' if good else 'FAILS'}")
    say("")
    if ok.get("3a") and ok.get("3b"):
        if ok.get("2"):
            say("   ==> DIAGNOSIS CONFIRMED AND CAVITY MATCH IS THE FIX. ddG works")
            say("       only where the pocket must close; cavity match works in both")
            say("       directions, including the expansion regime we care about.")
        else:
            say("   ==> cavity match works in both arms, but ddG did NOT fail on")
            say("       PFAS -- §68d's volume argument is wrong and ddG is more")
            say("       general than its mechanism suggested. Cavity match is still")
            say("       usable; the reason for preferring it is weaker.")
    elif ok.get("3a") and not ok.get("3b"):
        say("   ==> CAVITY MATCH DOES NOT TRANSFER EITHER. It works where ddG works")
        say("       and fails where ddG fails, so it is not the fix -- whatever")
        say("       defeats a protein-only score in the expansion regime defeats")
        say("       both. That points at the ligand being genuinely required.")
    elif not ok.get("3a"):
        say("   ==> CAVITY MATCH IS DEAD. It cannot even reproduce the arm where")
        say("       ddG already works, so the objective swap is not the answer.")
    say("")
    say("   ⚠ LIBRARY is unlabelled, not a set of known failures, so every AUC here")
    say("   is a lower bound on the separation ([[feedback_hits_are_not_optima]]).")

    # does the cavity actually move in the direction the ligand demands?
    say("")
    say("=" * 78)
    say("DIRECTION CHECK -- do real sensors move the cavity toward their ligand?")
    say("=" * 78)
    for klass, (by, meta) in data.items():
        wt = meta["wt_cavity"]
        tgt = meta["target_vol"]
        need = "GROW" if tgt > wt else "SHRINK"
        say(f"   {klass:<10} wild-type {wt:.1f} -> target {tgt:.1f} A^3, so the "
            f"cavity must {need}")
        for k in ("REAL", "LIBRARY", "WILD"):
            v = by.get(k, [])
            if not v:
                continue
            d = np.median([r["cav"] - r["cav_wt"] for r in v])
            say(f"     {k:<9} median cavity change vs its own paired wild-type: "
                f"{d:+7.1f} A^3")
    with open(os.path.join(OUT, "verdict.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    say(f"\n   written to {OUT}/verdict.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
