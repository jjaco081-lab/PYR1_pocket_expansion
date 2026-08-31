#!/usr/bin/env python
r"""
161_tractability_analysis.py -- the ONE >80 %-accuracy question we can actually
score, because it is the only one with real negatives.

§90 inventoried every label source in the corpus and found THREE variant-level
tested negatives. At the ligand level there are 3,172, from sd01's Hit? column,
and §84 property-matched 181 hits to 2,531 of them so that the benchmark cannot
be won by polarity or size alone (worst residual |AUC - 0.5| = 0.021 after
matching, against TPSA reaching 0.298 on a random draw).

So this is the only place a genuine accuracy number can be computed, and the
question it answers is the ligand-level form of Jannis's: GIVEN A MOLECULE, WILL
PYR1 YIELD A SENSOR FOR IT AT ALL?

METRICS, fixed before the numbers are seen. §90b's criteria:
  * balanced accuracy and MCC, never raw accuracy -- the trivial "no" classifier
    scores 93 % on the unmatched set and is quoted beside every number
  * AUC on the MATCHED set is the headline; the RANDOM set is scored too, but
    only to show how much easier the rigged version looks
  * the operating point that reaches 80 % balanced accuracy is reported if one
    exists, and its threshold and its recall are reported with it

WHAT IT CANNOT SHOW (§84c). A ligand-level result does not rescue variant
ranking: §82 settled that Boltz-2 cannot rank pocket variants for a fixed ligand
(within-ligand rho -0.028). This decides whether a campaign is worth starting,
not what to put in the library.

⚠ Both classes are folded against the SAME wild-type PYR1. A hit is a ligand some
VARIANT bound, so this asks whether WT co-folding predicts library tractability
-- a proxy one step from the label, but symmetric between the classes.
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
T = os.path.join(ROOT, "data", "tractability")
OUT = os.path.join(ROOT, "results", "tractability")
FIELDS = ["affinity_probability_binary", "affinity_pred_value"]


def load_scores():
    out = {}
    for d in sorted(glob.glob(os.path.join(T, "out", "boltz_results_lig*"))):
        lig = os.path.basename(d).replace("boltz_results_", "")
        hits = glob.glob(os.path.join(d, "predictions", "*", "affinity*.json"))
        if not hits:
            continue
        try:
            out[lig] = json.load(open(hits[0]))
        except Exception:
            continue
    return out


def auc(score, lab):
    o = np.argsort(-np.asarray(score, float))
    u = np.asarray(lab, bool)[o]
    p, n = np.where(u)[0], np.where(~u)[0]
    if not len(p) or not len(n):
        return float("nan")
    return float(np.mean([[1.0 if a < b else 0.5 if a == b else 0.0
                           for b in n] for a in p]))


def best_operating_point(score, lab):
    """Threshold maximising balanced accuracy; returns (bacc, mcc, thr, sens, spec)."""
    score = np.asarray(score, float)
    lab = np.asarray(lab, bool)
    best = None
    for thr in np.unique(score):
        pred = score >= thr
        tp = int((pred & lab).sum())
        fp = int((pred & ~lab).sum())
        fn = int((~pred & lab).sum())
        tn = int((~pred & ~lab).sum())
        sens = tp / max(tp + fn, 1)
        spec = tn / max(tn + fp, 1)
        bacc = 0.5 * (sens + spec)
        den = np.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        mcc = ((tp * tn - fp * fn) / den) if den else 0.0
        if best is None or bacc > best[0]:
            best = (bacc, mcc, float(thr), sens, spec)
    return best


def main():
    os.makedirs(OUT, exist_ok=True)
    key = json.load(open(os.path.join(T, "key_matched.json")))
    scores = load_scores()
    print(f"{len(scores)} of 362 runs have affinity output")
    if len(scores) < 300:
        print("  (incomplete -- numbers below are provisional)")

    rows = []
    for lig, meta in (key.items() if isinstance(key, dict) else
                      {r["id"]: r for r in key}.items()):
        s = scores.get(lig)
        if not s:
            continue
        rec = {"lig": lig, "hit": bool(meta.get("hit", meta.get("label")))}
        for f in FIELDS:
            rec[f] = float(s.get(f, np.nan))
        rows.append(rec)
    lab = np.array([r["hit"] for r in rows])
    print(f"scored: {lab.sum()} hits, {(~lab).sum()} matched non-hits")
    if lab.sum() < 5 or (~lab).sum() < 5:
        print("not enough of one class yet")
        return 0

    triv = max(lab.mean(), 1 - lab.mean())
    print(f"\ntrivial-majority baseline: accuracy {triv:.0%}, "
          f"balanced accuracy 50 %, MCC 0.000")
    print("\n" + "=" * 72)
    print("CAN BOLTZ-2 TELL A TRACTABLE LIGAND FROM A MATCHED INTRACTABLE ONE?")
    print("=" * 72)
    print(f"{'field':<34}{'AUC':>7}{'bal.acc':>9}{'MCC':>8}{'sens':>7}{'spec':>7}")
    res = {}
    for f in FIELDS:
        v = np.array([r[f] for r in rows])
        m = ~np.isnan(v)
        if m.sum() < 10:
            continue
        a = auc(v[m], lab[m])
        bacc, mcc, thr, sens, spec = best_operating_point(v[m], lab[m])
        print(f"{f:<34}{a:>7.3f}{bacc:>9.3f}{mcc:>8.3f}{sens:>7.2f}{spec:>7.2f}")
        res[f] = dict(auc=a, bacc=bacc, mcc=mcc, thr=thr, sens=sens, spec=spec,
                      n=int(m.sum()))
    print("\n   bal.acc is at the BEST threshold chosen on this same data, so it")
    print("   is an optimistic ceiling, not a held-out estimate. Read AUC first.")
    hit80 = [f for f, d in res.items() if d["bacc"] >= 0.80]
    print(f"\n   fields reaching 80 % balanced accuracy: "
          f"{hit80 if hit80 else 'NONE'}")
    json.dump({"rows": rows, "metrics": res},
              open(os.path.join(OUT, "tractability.json"), "w"), indent=1)
    print(f"\nwritten to {OUT}/tractability.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
