#!/usr/bin/env python
"""
35_variance_decomposition.py -- CORRECTS the discrimination test used in scripts
20, 23 and 33, and re-reports every arm on the fixed statistic.

The error
---------
Those scripts asked whether

    sd(per-variant means)  /  mean(within-variant replicate sd)   >=  2

That comparison is wrong. Replicate sd is the spread of INDIVIDUAL relax runs.
What limits our ability to separate variants is the spread of their ESTIMATED
MEANS, which is sd/sqrt(n) -- smaller by 4.5x at n=20. Comparing a between-
variant spread against the raw replicate sd therefore understates power by that
same factor, and buries real effects whenever n is large.

The right test is the one-way ANOVA F ratio:

    F = s2_between / (s2_within / n)

with the true between-variant variance recovered by subtracting off the part of
the observed spread that sampling error alone would produce:

    s2_true = max(0, s2_between - s2_within / n)

s2_true is the number that matters for design: it is the actual spread of the
underlying per-variant values, in physical units, independent of how hard we
sampled.

Consequences of the fix, which change conclusions:
  * Arm 2b `dE_close` was labelled "NOT DISCRIMINATING" at ratio 0.88 by the old
    rule. Its F is ~9, and its true between-variant sd is ~1.8 REU. It works.
  * Arm 2 (v2/v3) remains weak, but the correct statement is not "the assay is
    blind" -- it is that the true between-variant spread of dG_separated is
    ~0.5 REU against a -67 REU interface, i.e. the effect is real but negligible.

Usage:  python 35_variance_decomposition.py
"""
import glob, json, os, collections
import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")


def load_grouped(pattern, key="reps"):
    per = collections.defaultdict(list)
    for f in sorted(glob.glob(pattern)):
        j = json.load(open(f))
        per[j["variant"]].extend(j[key])
    for v in per:
        seen, uniq = set(), []
        for r in sorted(per[v], key=lambda x: x["rep"]):
            if r["rep"] not in seen:
                seen.add(r["rep"])
                uniq.append(r)
        per[v] = uniq
    return per


def decompose(per, key, derived=None):
    groups = []
    for v, reps in per.items():
        if derived:
            x = [derived(r) for r in reps if all(k in r for k in derived.keys)] \
                if hasattr(derived, "keys") else [derived(r) for r in reps]
        else:
            x = [r[key] for r in reps if key in r and r[key] is not None]
        if len(x) > 1:
            groups.append(np.asarray(x, float))
    if len(groups) < 2:
        return None
    n = float(np.mean([len(g) for g in groups]))
    means = np.array([g.mean() for g in groups])
    s2_within = float(np.mean([g.var(ddof=1) for g in groups]))
    s2_between = float(means.var(ddof=1))
    se2 = s2_within / n
    F = s2_between / se2 if se2 > 0 else float("nan")
    s2_true = max(0.0, s2_between - se2)
    try:
        p = float(stats.f_oneway(*groups).pvalue)
    except Exception:
        p = float("nan")
    return dict(n=n, k=len(groups), within_sd=np.sqrt(s2_within),
                between_sd=np.sqrt(s2_between), sem=np.sqrt(se2),
                true_sd=np.sqrt(s2_true), F=F, p=p,
                old_ratio=np.sqrt(s2_between) / np.sqrt(s2_within))


ARMS = []

cav = load_grouped(os.path.join(R, "cavity_scan", "*.json"))
ARMS.append(("Arm 1  cavity_A3", cav, "cavity_A3"))
ARMS.append(("Arm 1  score_ref2015", cav, "score_ref2015"))

for tag, sub in (("v2", "ratchet_v2"), ("v3", "ratchet_v3")):
    d = os.path.join(R, sub)
    if glob.glob(os.path.join(d, "*__c*.json")):
        rat = load_grouped(os.path.join(d, "*__c*.json"))
        for k in ("dG_separated", "gate_bb_rmsd", "latch_bb_rmsd",
                  "gate_latch_min_dist", "w385_aba_min_dist"):
            if any(k in r for reps in rat.values() for r in reps):
                ARMS.append((f"Arm 2 {tag}  {k}", rat, k))

sw = load_grouped(os.path.join(R, "switch", "*__r*.json"))
if sw:
    for k in ("dE_close", "score_closed", "score_open"):
        ARMS.append((f"Arm 2b  {k}", sw, k))

print("=" * 104)
print("VARIANCE DECOMPOSITION -- corrected discrimination test")
print("=" * 104)
print(f"{'arm / metric':<34}{'n':>4}{'within sd':>11}{'SEM':>8}"
      f"{'between sd':>12}{'TRUE sd':>10}{'F':>8}{'p':>10}{'old ratio':>11}")
for label, per, key in ARMS:
    d = decompose(per, key)
    if d is None:
        continue
    p = f"{d['p']:.2e}" if np.isfinite(d["p"]) else "--"
    print(f"{label:<34}{d['n']:>4.0f}{d['within_sd']:>11.3f}{d['sem']:>8.3f}"
          f"{d['between_sd']:>12.3f}{d['true_sd']:>10.3f}{d['F']:>8.2f}"
          f"{p:>10}{d['old_ratio']:>11.2f}")

print("""
  within sd  = spread of individual relax replicates
  SEM        = within sd / sqrt(n): the noise on each variant's MEAN
  between sd = observed spread of the per-variant means
  TRUE sd    = sqrt(between^2 - SEM^2): the real spread, sampling error removed.
               THIS is the number to judge an effect by.
  F          = s2_between / (s2_within/n); F >> 1 means the variants really differ
  old ratio  = the incorrect statistic used in scripts 20/23/33, shown so the
               earlier conclusions in this README can be re-read against it

Judge an arm by TRUE sd in physical units, not by F alone: with enough replicates
a negligible effect becomes statistically significant while remaining useless.""")
