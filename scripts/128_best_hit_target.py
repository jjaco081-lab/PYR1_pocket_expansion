#!/usr/bin/env python
r"""
128_best_hit_target.py -- change the target from "any known hit" to "the BEST
known hit", using the dose-response data.

WHY (Jannis, 2026-08-26)
Every metric so far asks whether a designed library contains >= 1 known round-2
hit. That has a weakness he identified: those hits came out of Tian's library, so
if that library was built suboptimally we are asking a method to rediscover a
mediocre answer. Raising the bar to the BEST observed sensor removes most of that
objection -- whatever the library's flaws, its best output is the best thing anyone
has actually measured for that ligand.

sd07 carries a full dose-response ladder (0.025 - 100 uM, 12 levels), and the
within-ligand spread is up to 50x (4-methylumbelliferone: 0.5 vs 25 uM), so "any
hit" and "best hit" are genuinely different targets.

WHAT THIS CANNOT DO, stated plainly
It cannot find the best possible SEQUENCE. That would need binding measurements
for sequences nobody made, and no scorer in this project predicts affinity
(§40-§46). What it CAN do is enumerate every possible LIBRARY exactly and ask
which is smallest subject to containing the best MEASURED sensor -- the ground
truth is still positive-unlabeled ([[feedback_hits_are_not_optima]]), but it is now
the strongest label available rather than an arbitrary one.
"""
import itertools
import json
import os
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_xlsx import table                                     # noqa: E402
import importlib.util                                          # noqa: E402

_s = importlib.util.spec_from_file_location(
    "m116", os.path.join(HERE, "116_secondary_library.py"))
M = importlib.util.module_from_spec(_s)
_s.loader.exec_module(M)

ROOT = os.path.dirname(HERE)
SD = M.SD
OUT = os.path.join(ROOT, "results", "secondary_library")
P18 = ["K59", "V81", "V83", "L87", "A89", "S92", "E94", "F108", "I110", "L117",
       "Y120", "S122", "E141", "F159", "A160", "V163", "V164", "N167"]


def sensors():
    """Round-2 coumarin sensors with their minimum responsive concentration."""
    hdr, recs = table(f"{SD}/pnas.2519924122.sd07(1).xlsx")
    doses = sorted((h for h in hdr if h.replace(".", "", 1).isdigit()), key=float)
    cols = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    out = []
    for r in recs:
        pot = None
        for d in doses:
            if (r.get(d) or "").strip() == "+":
                pot = float(d)
                break                       # doses are sorted, first + is the min
        if pot is None:
            continue
        st = {p: p[0] for p in P18}
        off = False
        for c in cols:
            v = (r.get(c) or "").strip().upper()
            if len(v) != 1 or v not in "ACDEFGHIKLMNPQRSTVWY":
                continue
            if c in P18:
                st[c] = v
            elif v != c[0]:
                off = True
        out.append({"lig": (r.get("compound") or "").strip(), "pot": pot,
                    "st": st, "off": off,
                    "name": (r.get("name") or "").strip()})
    return out


def exact_min(byl, ligs, forcing):
    best = [float("inf"), None]
    order = sorted(ligs, key=lambda l: (len(byl[l]), l))

    def size(cols):
        n = 1
        for i, s in enumerate(cols):
            n *= max(1, len(s | ({P18[i][0]} if not forcing else set())))
        return n

    def rec(i, cols):
        if size(cols) >= best[0]:
            return
        if i == len(order):
            best[0] = size(cols)
            best[1] = [set(c) for c in cols]
            return
        for cl in byl[order[i]]:
            add = []
            for k, aa in enumerate(cl):
                if aa not in cols[k]:
                    cols[k].add(aa)
                    add.append((k, aa))
            rec(i + 1, cols)
            for k, aa in add:
                cols[k].discard(aa)

    rec(0, [set() for _ in P18])
    return best[0], best[1]


def main():
    os.makedirs(OUT, exist_ok=True)
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    S = [s for s in sensors() if not s["off"]]
    ligs = sorted({s["lig"] for s in S})
    say("=" * 78)
    say("1. THE BEST MEASURED SENSOR PER LIGAND")
    say("=" * 78)
    say(f"   {len(S)} round-2 clones inside the 18 positions, {len(ligs)} ligands")
    say(f"   {'ligand':<32}{'n':>4}{'best uM':>9}{'worst':>8}{'spread':>8}"
        f"{'#tied at best':>14}")
    best_of = {}
    for l in ligs:
        v = sorted((s["pot"] for s in S if s["lig"] == l))
        b = v[0]
        tied = [s for s in S if s["lig"] == l and s["pot"] == b]
        best_of[l] = tied
        say(f"   {l:<32}{len(v):>4}{b:>9.3f}{v[-1]:>8.1f}{v[-1]/b:>7.0f}x"
            f"{len(tied):>14}")

    say("")
    say("=" * 78)
    say("2. SMALLEST LIBRARY: 'any hit' vs 'the BEST hit'")
    say("=" * 78)
    variants = {
        "any hit": {l: sorted({tuple(s["st"][p] for p in P18)
                               for s in S if s["lig"] == l}) for l in ligs},
        "BEST hit": {l: sorted({tuple(s["st"][p] for p in P18)
                                for s in best_of[l]}) for l in ligs},
    }
    tian = M.size_of(M.tian_menu("Coumarin"))
    res = {}
    say(f"   {'target':<12}{'WT offered':>14}{'vs Tian':>10}{'forcing':>12}"
        f"{'vs Tian':>10}")
    for name, byl in variants.items():
        row = []
        for forcing in (False, True):
            sz, _cols = exact_min(byl, ligs, forcing)
            row.append(sz)
        res[name] = row
        say(f"   {name:<12}{row[0]:>14,}{tian/row[0]:>9.0f}x{row[1]:>12,}"
            f"{tian/row[1]:>9.0f}x")
    say("")
    say(f"   ==> demanding the BEST sensor rather than any sensor costs "
        f"{res['BEST hit'][0]/res['any hit'][0]:.1f}x in library size")
    say(f"       ({res['any hit'][0]:,} -> {res['BEST hit'][0]:,}), and the target "
        f"is still {tian/res['BEST hit'][0]:.0f}x under Tian's {tian:,}.")

    say("")
    say("=" * 78)
    say("3. WHAT DO THE DESIGNED LIBRARIES ACTUALLY CATCH? potency of what they hold")
    say("=" * 78)
    r1, _ = M.load(f"{SD}/pnas.2519924122.sd03(1).xlsx", "library_name",
                   "canonical_smiles")
    r2raw, _ = M.load(f"{SD}/pnas.2519924122.sd07(1).xlsx", "compound")
    sim = M.similarity(r1, ligs, M.EXTRA_SMILES)
    subw = M.subw_from(r1, sim, 4.0)
    say(f"   {'library':<26}{'size':>11}{'ligands':>9}{'best-hit':>10}"
        f"   median potency caught / best available")
    for lbl, frac in (("class/profile f=0.50", 0.50), ("class/profile f=0.35", 0.35),
                      ("class/profile f=0.25", 0.25)):
        menu = M.design_profile_depth(subw, 138240, set(P18), 11, frac)[0]
        flat = {p + a for p, v in menu.items() for a in v}
        got, gotbest, ratios = set(), set(), []
        for l in ligs:
            held = [s for s in S if s["lig"] == l
                    and all(f"{p}{s['st'][p]}" in flat
                            for p in P18 if s["st"][p] != p[0])]
            if held:
                got.add(l)
                b = min(x["pot"] for x in held)
                ratios.append(b / best_of[l][0]["pot"])
                if any(x["pot"] == best_of[l][0]["pot"] for x in held):
                    gotbest.add(l)
        import numpy as np
        med = np.median(ratios) if ratios else float("nan")
        say(f"   {lbl:<26}{M.size_of(menu):>11,}{len(got):>7}/11{len(gotbest):>8}/11"
            f"   {med:.1f}x worse than the best")
    say("")
    say("   ==> a library that captures 'a hit' does NOT generally capture the best")
    say("       one, so the two metrics are not interchangeable and the harder one")
    say("       is the one worth designing against.")

    with open(os.path.join(OUT, "best_hit.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    say(f"\n   written to {OUT}/best_hit.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
