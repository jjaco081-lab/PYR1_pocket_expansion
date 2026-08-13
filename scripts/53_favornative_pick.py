#!/usr/bin/env python
"""
53_favornative_pick.py -- choose the favor-native weight from the NULL arm alone.

The ABA arm is WT protein with its native ligand, so its correct answer is zero
mutations at all 15 designable positions. Job 27421344 (favor_native=0) kept WT at
only 4 of 15 and deleted the K59 carboxylate salt bridge in 100% of trajectories,
which is why K59R and V81I came out uninterpretable rather than negative: at a
position the null also mutates, "no ligand-conditional signal" cannot be told apart
from "the protocol cannot hold a native contact".

PRE-REGISTERED CRITERION -- fixed before any sweep output was read, and stated in
52_submit_favornative_sweep.sh at submission time:

    pick the SMALLEST weight whose ABA null satisfies BOTH
      (a) keeps WT at >= 70% of the 15 designable positions, and
      (b) restores K59, i.e. >= 50% Lys at position 59
    in BOTH alphabets.

Only the null arm enters this decision. The mandipropamid arms are not run in the
sweep, so the weight cannot be tuned toward the result the benchmark is testing for.

FREEZE CHECK: a large enough bonus trivially satisfies both conditions by freezing
the packer, which would make the whole benchmark vacuous. Distinct sequences and
mean per-position entropy are reported so that failure is visible. A weight that
passes (a) and (b) but collapses to 1 distinct sequence is reported as FROZEN and
disqualified.

Usage:  python 53_favornative_pick.py
"""
import collections
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SWEEP = os.path.join(ROOT, "results", "stage1_rosetta_fnsweep")
BASE = os.path.join(ROOT, "results", "stage1_rosetta")
SUMMARY = "stage1_rosetta_summary.json"

DESIGN = [59, 81, 83, 92, 94, 108, 110, 120, 122, 141, 159, 160, 163, 164, 167]
WT = {59: "K", 81: "V", 83: "V", 92: "S", 94: "E", 108: "F", 110: "I", 120: "Y",
      122: "S", 141: "E", 159: "F", 160: "A", 163: "V", 164: "V", 167: "N"}
ALPHABETS = ["dsm_hao", "free"]

WT_FLOOR = 0.70          # (a)
K59_FLOOR = 0.50         # (b)
K59 = 59


def load(d):
    """arm -> trajectories, for one directory of 50_stage1_rosetta.py output."""
    arms = collections.defaultdict(list)
    if not os.path.isdir(d):
        return arms
    for f in sorted(os.listdir(d)):
        if not f.endswith(".json") or f == SUMMARY:
            continue
        j = json.load(open(os.path.join(d, f)))
        if "arm" in j:
            arms[j["arm"]].extend(j["trajectories"])
    return arms


def stats(rows):
    """WT recovery, K59 retention, and how much freedom is left."""
    n = len(rows)
    if not n:
        return None
    keeps = []
    for p in DESIGN:
        col = [r["seq"][DESIGN.index(p)] for r in rows]
        keeps.append(sum(c == WT[p] for c in col) / n)
    k59 = keeps[DESIGN.index(K59)]
    seqs = [r["seq"] for r in rows]
    ent = []
    for p in DESIGN:
        c = collections.Counter(r["seq"][DESIGN.index(p)] for r in rows)
        ent.append(-sum((v / n) * math.log(v / n, 20) for v in c.values()))
    return dict(n=n,
                wt_frac=sum(keeps) / len(keeps),
                wt_pos=sum(k >= 0.5 for k in keeps),
                k59=k59,
                distinct=len(set(seqs)),
                entropy=sum(ent) / len(ent))


# favor_native=0 baseline comes from the original campaign, not the sweep
runs = {0.0: load(BASE)}
if os.path.isdir(SWEEP):
    for sub in sorted(os.listdir(SWEEP)):
        if sub.startswith("w"):
            runs[float(sub[1:])] = load(os.path.join(SWEEP, sub))

print("=" * 88)
print("ABA NULL ARM vs favor-native weight   (correct answer: 15/15 WT, 100% K59)")
print("=" * 88)
print(f"{'weight':>7} {'alphabet':>9} {'n':>4} {'WT frac':>8} {'WT pos':>7} "
      f"{'K59 Lys':>8} {'distinct':>9} {'entropy':>8}  verdict")

passing = collections.defaultdict(dict)
for w in sorted(runs):
    for a in ALPHABETS:
        s = stats(runs[w].get(f"wt_aba__{a}", []))
        if not s:
            continue
        ok_wt, ok_k59 = s["wt_frac"] >= WT_FLOOR, s["k59"] >= K59_FLOOR
        frozen = s["distinct"] <= 1
        verdict = ("FROZEN -- disqualified" if frozen and ok_wt and ok_k59
                   else "pass" if ok_wt and ok_k59
                   else "fail: " + ", ".join(
                       ([] if ok_wt else [f"WT {s['wt_frac']:.0%}<{WT_FLOOR:.0%}"])
                       + ([] if ok_k59 else [f"K59 {s['k59']:.0%}<{K59_FLOOR:.0%}"])))
        passing[w][a] = ok_wt and ok_k59 and not frozen
        print(f"{w:>7g} {a:>9} {s['n']:>4} {s['wt_frac']:>7.0%} "
              f"{s['wt_pos']:>4}/{len(DESIGN)} {s['k59']:>7.0%} "
              f"{s['distinct']:>9} {s['entropy']:>8.3f}  {verdict}")
    print()

winners = [w for w in sorted(passing)
           if len(passing[w]) == len(ALPHABETS) and all(passing[w].values())]

print("=" * 88)
if winners:
    pick = winners[0]
    print(f"PICK: favor_native = {pick:g}")
    print(f"  smallest weight passing WT>={WT_FLOOR:.0%} and K59>={K59_FLOOR:.0%} "
          f"in both alphabets, without freezing")
    print(f"\nRe-run all six arms at this weight, into a fresh outdir:")
    print(f"  sbatch --export=ALL,FN={pick:g} scripts/54_submit_stage1_rosetta_fn.sh")
    print(f"  python scripts/51_stage1_rosetta_merge.py \\")
    print(f"      --indir results/stage1_rosetta_fn{pick:g}")
    print(f"\nFLOOR CHECK before believing the re-run: the mandipropamid arm must")
    print(f"  still mutate F108 at high frequency. If the bonus has quieted the")
    print(f"  mandi arm too, the benchmark is measuring the bonus, not the ligand.")
else:
    print("NO WEIGHT PASSES.")
    print("  Do not widen the sweep to make one pass -- that is tuning the control")
    print("  until it agrees. If even 1.5 leaves the null mutating, the problem is")
    print("  not the reference energies: suspect the pose, the ligand params, or the")
    print("  8 A pack shell. Report the Rosetta arm as inconclusive instead.")
    sys.exit(1)
