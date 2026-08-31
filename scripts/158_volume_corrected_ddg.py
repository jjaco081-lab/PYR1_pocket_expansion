#!/usr/bin/env python
r"""
158_volume_corrected_ddg.py -- a PYR1-tuned correction to ref2015 ddG that works.

Jannis's proposition: the scoring function may have to be tailored to PYR1
rather than generalised. §157 gives the first evidence that this is right AND
identifies what to correct, on 380 single substitutions at Beltran's 20
positions with a measured replicate noise floor of 0.21 REU.

THE BIAS. Grouping the 380 by side-chain volume change (Zamyatnin):

    SHRINK (< -20 A^3)   n=179   median ddG  +2.40
    similar              n= 65   median ddG  +0.04
    GROW   (> +20 A^3)   n=136   median ddG  +0.54

ref2015 penalises SHRINK substitutions specifically -- removing a side chain
leaves a void it scores as destabilising. But a void is precisely what a ligand
pocket needs, and real sensors are shrink-biased: the 40 substitutions used by
Beltran's sensors have a median volume change of -16.3 A^3 against -3.7 for the
unlabelled remainder. So the filter is penalising the class it should favour,
which is why §89b found it buys only 1.9x at zero cost.

THE CORRECTION IS LABEL-FREE, which is what makes it legitimate. A quadratic in
volume change is fitted to ddG using ONLY the 340 UNLABELLED substitutions and
subtracted. No sensor label enters the fit, so evaluating on the 40 is not
circular. Within a single volume class the raw score is already better than
pooled (AUC 0.652 shrink-only, 0.683 grow-only, against 0.606 pooled), which is
the signature of a between-class offset rather than a within-class failure --
exactly what a constant-per-class correction should fix.

⚠ WHAT THIS IS NOT. It is not a new energy function and it is not validated
out-of-sample: the correction is fitted and evaluated on the same 380 (label-free
fit, but the same structural context). It must be re-tested on the Park-475
mandipropamid library, which is a different ligand, a different backbone and a
different lab. That job is 27977098.
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "library_recall")

#: Zamyatnin side-chain volumes, A^3
VOL = dict(zip("AGVLIPFMWSTCYNQDEKRH",
               [88.6, 60.1, 140.0, 166.7, 166.7, 112.7, 189.9, 162.9, 227.8,
                89.0, 116.1, 108.5, 193.6, 114.1, 143.8, 111.1, 138.4, 168.6,
                173.4, 153.2]))


def auc(score, lab):
    o = np.argsort(score)
    u = lab[o]
    p, n = np.where(u)[0], np.where(~u)[0]
    return float(np.mean([[1.0 if a < b else 0.5 if a == b else 0.0
                           for b in n] for a in p]))


def recall(score, lab, frac):
    k = int(round(frac * len(score)))
    return lab[np.argsort(score)][:k].sum() / lab.sum()


def fit_volume_trend(dv, dg, mask):
    """Quadratic ddG ~ dV fitted on `mask` only. Returns the predictor."""
    X = np.c_[np.ones(mask.sum()), dv[mask], dv[mask] ** 2]
    b = np.linalg.lstsq(X, dg[mask], rcond=None)[0]
    return b, (lambda v: np.c_[np.ones(len(v)), v, v ** 2] @ b)


def main():
    rows = []
    for f in glob.glob(os.path.join(ROOT, "results", "beltran_ddg", "chunk_*.json")):
        rows += json.load(open(f))["rows"]
    bel = json.load(open(os.path.join(ROOT, "data", "beltran", "win_sensors.json")))
    used = {m for s in bel["sensors"] if s["muts"] for m in s["muts"]}
    dv = np.array([VOL[r["mut"]] - VOL[r["wt"]] for r in rows])
    dg = np.array([r["ddG"] for r in rows])
    lab = np.array([r["sub"] in used for r in rows])
    print(f"{len(rows)} singles, {lab.sum()} used by a real sensor")

    print(f"\n{'class':<22}{'n':>5}{'median ddG':>13}{'median dV':>12}")
    for name, m in (("SHRINK (< -20 A^3)", dv < -20),
                    ("similar", (dv >= -20) & (dv <= 20)),
                    ("GROW (> +20 A^3)", dv > 20)):
        print(f"{name:<22}{m.sum():>5}{np.median(dg[m]):>+13.2f}{np.median(dv[m]):>+12.1f}")
    print(f"\nreal sensor substitutions: median dV {np.median(dv[lab]):+.1f} A^3 "
          f"({int((dv[lab] < -20).sum())} shrink, {int((dv[lab] > 20).sum())} grow)")
    print(f"unlabelled:                median dV {np.median(dv[~lab]):+.1f} A^3")

    b, pred = fit_volume_trend(dv, dg, ~lab)
    corr = dg - pred(dv)
    print(f"\ncorrection fitted on the {int((~lab).sum())} UNLABELLED only "
          f"(no label enters the fit):")
    print(f"   ddG_hat = {b[0]:+.2f} {b[1]:+.4f}*dV {b[2]:+.6f}*dV^2")

    print("\n" + "=" * 66)
    print("DOES IT HELP?")
    print("=" * 66)
    print(f"{'score':<26}{'AUC':>7}{'top10%':>9}{'top25%':>9}{'top50%':>9}")
    out = {}
    for name, s in (("raw ref2015 ddG", dg), ("volume-corrected", corr)):
        print(f"{name:<26}{auc(s, lab):>7.3f}{recall(s, lab, .10):>9.0%}"
              f"{recall(s, lab, .25):>9.0%}{recall(s, lab, .50):>9.0%}")
        out[name] = dict(auc=auc(s, lab),
                         **{f"top{int(f*100)}": recall(s, lab, f)
                            for f in (.10, .25, .50)})
    print(f"{'chance':<26}{0.5:>7.3f}{10:>8}%{25:>8}%{50:>8}%")

    print("\n   substitutions that move most, and the ones that do not:")
    ro = np.argsort(dg).argsort()
    rc = np.argsort(corr).argsort()
    idx = {r["sub"]: i for i, r in enumerate(rows)}
    print(f"   {'sub':>8}{'raw':>10}{'corrected':>12}")
    for s in ("F159G", "Y120A", "Y120G", "A160G", "V164W", "A160V", "V83L"):
        if s in idx:
            i = idx[s]
            print(f"   {s:>8}{ro[i]+1:>7}/380{rc[i]+1:>9}/380")
    print("\n   F159G is Beltran's round-1 DSM hit for WIN 55,212 (with A160V).")
    print("   Raw ddG buries it at 258; the correction lifts it to 41. Y120G, the")
    print("   most widely used substitution in the whole set, stays at 275 -- the")
    print("   correction is partial and the head of the list barely improves.")
    json.dump(out, open(os.path.join(OUT, "volume_corrected.json"), "w"), indent=1)
    print(f"\nwritten to {OUT}/volume_corrected.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
