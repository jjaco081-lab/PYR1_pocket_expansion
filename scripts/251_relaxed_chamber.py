#!/usr/bin/env python
r"""
251_relaxed_chamber.py -- THE COLLAPSE TEST: re-measure Arm 1's RELAXED
structures with the calibrated chamber probe.

Jannis asked what test should come next to show the F108/R79/E94 triple really
opens the second lobe, and whether anything should be run to show the pocket is
not predicted to collapse. The collapse experiment already exists -- I had
forgotten it, and he was right that it was done early.

`06_rosetta_cavity_scan.py` (Arm 1 of the pilot) is exactly that experiment:
  * MutateResidue, then FastRelax with ref2015_cst, 3 repeats
  * CA coordinate constraints to the closed backbone, so it measures SIDE-CHAIN
    INFILLING rather than a global conformational change
  * side chains UNCONSTRAINED and free to fall into the cavity -- the collapse
    IS the signal
  * APO relax, deliberately: no ligand propping the pocket open
  * WT relaxed through the identical protocol as the internal reference
and its relaxed coordinates are still on disk at results/cavity_scan/pdb/.

WHAT IS NEW HERE. Arm 1's cavity came from a seeded volume at the ABA centroid,
not the calibrated 1.4 A chamber probe. That matters because the same older
criterion credited F108G with +84 A^3 that turned out to be only +1 A^3 of
USABLE chamber. So the relaxed structures are re-measured with the probe that
was calibrated against ABA (18 of 19 atoms inside PYR1's own chamber), and
`n_chamber` is reported so fragmentation cannot hide inside a total.

⚠ APO AND CHAIN A ONLY -- that is what Arm 1 relaxed. The +HAB1 numbers from 250
are rigid-truncation and are NOT comparable to these.
⚠ Rigid (250) vs relaxed (here) on the SAME probe is the collapse measurement:
the difference is how much the neighbours reclaim once they are allowed to move.
"""
import glob, json, os, re, sys
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
m229 = __import__("importlib").import_module("229_cavity_visualise")
PDBD = os.path.join(ROOT, "results", "cavity_scan", "pdb")
#: rigid, calibrated-probe values from 250 (chain A, apo) for the same variants
RIGID = {"WT": 164.4, "R79A": 161.2, "E94A": 180.4, "R79A_E94A": 251.0,
         "F108A": 241.7, "F108A_R79A": 267.2, "F108A_R79A_E94A": 311.0,
         "K59A_F108A_R79A_E94A": 385.1, "K59A_F108A": None}
WANT = ["WT", "R79A", "E94A", "R79A_E94A", "F108A", "F108A_R79A",
        "F108A_R79A_E94A", "K59A_F108A", "K59A_F108A_R79A_E94A"]


def measure(path):
    X, el = [], []
    for l in open(path):
        if not l.startswith("ATOM"):
            continue
        e = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if e in ("H", "D"):
            continue
        X.append([float(l[30:38]), float(l[38:46]), float(l[46:54])]); el.append(e)
    if len(X) < 300:
        return None
    X = np.array(X)
    cs = m229.chambers_with_voxels(X, el, [True] * len(X))
    if not cs:
        return (0.0, 0.0, 0)
    return (cs[0]["vol"], sum(c["vol"] for c in cs), len(cs))


def main():
    by = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(PDBD, "*.pdb"))):
        m = re.match(r"(.+)_r(\d+)\.pdb$", os.path.basename(f))
        if m:
            by[m.group(1)].append(f)
    print(f"{len(by)} relaxed variants on disk, "
          f"{sum(len(v) for v in by.values())} structures\n")
    have = [w for w in WANT if w in by]
    missing = [w for w in WANT if w not in by]
    res = {}
    for v in have:
        vals = [measure(f) for f in by[v]]
        vals = [x for x in vals if x]
        if not vals:
            continue
        main_ = [x[0] for x in vals]; nch = [x[2] for x in vals]
        res[v] = (float(np.mean(main_)), float(np.std(main_)), len(main_),
                  float(np.mean(nch)))
    wt = res.get("WT", (None,))[0]
    print(f"{'variant':<24}{'RELAXED main':>14}{'SD':>7}{'n':>3}{'chmb':>6}"
          f"{'vs WT':>8}{'RIGID':>8}{'reclaimed':>11}")
    for v in have:
        if v not in res:
            continue
        m_, sd, n, nc = res[v]
        rg = RIGID.get(v)
        d = m_ - wt if wt else 0.0
        rec = (f"{rg - m_:+.1f}" if rg else "-")
        print(f"{v:<24}{m_:>14.1f}{sd:>7.1f}{n:>3}{nc:>6.1f}{d:>+8.1f}"
              f"{(f'{rg:.1f}' if rg else '-'):>8}{rec:>11}")
    if missing:
        print(f"\nnot relaxed in Arm 1: {missing}")
    # the collapse verdict
    t = res.get("F108A_R79A_E94A")
    if t and wt:
        rg = RIGID["F108A_R79A_E94A"]
        kept = (t[0] - wt) / (rg - RIGID["WT"]) if rg else float("nan")
        print(f"\n{'='*64}\nDOES THE TRIPLE COLLAPSE ON RELAXATION?\n{'='*64}")
        print(f"   rigid   {rg:.1f} A^3  ({rg-RIGID['WT']:+.1f} vs WT)")
        print(f"   relaxed {t[0]:.1f} +/- {t[1]:.1f} A^3  ({t[0]-wt:+.1f} vs WT), "
              f"n={t[2]}, {t[3]:.1f} chamber(s)")
        print(f"   -> {100*kept:.0f} % of the rigid gain SURVIVES repacking")
        print(f"   (side chains were free and apo; collapse was the signal)")
    json.dump({k: dict(relaxed_main=v[0], sd=v[1], n=v[2], n_chamber=v[3],
                       rigid=RIGID.get(k)) for k, v in res.items()},
              open(os.path.join(ROOT, "results", "cavity_scan",
                                "relaxed_calibrated.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
