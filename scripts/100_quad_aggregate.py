#!/usr/bin/env python
r"""
100_quad_aggregate.py -- PYR1^MANDI quadruple TI: the 13e calibration.

    dG(leg) = sum_i w_i <dV/dlambda>_i        12-point Gauss-Legendre
    ddG_bind(L)  = dG(complex with L) - dG(apo)
    selectivity  = ddG(mandi) - ddG(ABA) = dG(mandi leg) - dG(ABA leg)

The apo leg cancels exactly in selectivity, so that number depends on four fewer
error sources than the absolute ddG and is the one to read first.

TWO MANDIPROPAMID ARMS, differing only in where the pose came from:
  mandi_xtal  4WVO's crystallographic pose and its own mutant rotamers
  mandi_dock  docked into a MODEL of the quad pocket, no complex structure used
Both share the same ABA leg, so their DIFFERENCE isolates pose provenance.
Measured at build time the docked pose is 8.33 A heavy-atom RMSD from the crystal
one (best of 9 returned poses: 4.90 A), so this is not a small perturbation --
it is what TI does when handed a pose that is in the right pocket and wrong.

THE PRE-REGISTERED TEST (README 13e, 44e)
PYR1^MANDI + mandipropamid must rank above WT + mandipropamid: the crystal arm's
selectivity must come out NEGATIVE. Unlike 43d this is a real constraint -- the
quadruple is the sequence the ground truth is actually about.

Uncertainty is reported as three separate terms (stat / conv / quad) for the
reasons in 88's header: a within-run error bar cannot see a drift it is centred
on, and pooling them is how the previous two results went confidently wrong.
"""
import csv
import importlib.util
import math
import os
import sys

import numpy as np
from scipy.interpolate import CubicSpline

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TI = os.path.join(ROOT, "data", "ti_quad")
LEGS = ["aba", "mandi_xtal", "mandi_dock", "apo"]
LAM = [0.00922, 0.04794, 0.11505, 0.20634, 0.31608, 0.43738,
       0.56262, 0.68392, 0.79366, 0.88495, 0.95206, 0.99078]
W = [0.02359, 0.05347, 0.08004, 0.10158, 0.11675, 0.12457,
     0.12457, 0.11675, 0.10158, 0.08004, 0.05347, 0.02359]

# reuse the corrected parser from 89 rather than re-deriving it: it de-duplicates
# pmemd's two printed copies of every frame and drops the AVERAGES/RMS banners
_s = importlib.util.spec_from_file_location(
    "reparse", os.path.join(ROOT, "scripts", "89_ti_reparse.py"))
_m = importlib.util.module_from_spec(_s)
_s.loader.exec_module(_m)


def window(leg, i):
    wd = os.path.join(TI, leg, f"lam{i:02d}")
    x = _m.series_dir(wd)          # spans every prod*.out segment (see 105)
    if x is None or len(x) == 0:
        return None
    h = len(x) // 2
    t = _m.tau_int(x)
    neff = len(x) / (1 + 2 * t)
    return dict(mean=x.mean(), n=len(x), tau=t, n_eff=neff,
                se=x.std(ddof=1) / np.sqrt(max(neff, 1)),
                first=x[:h].mean(), second=x[h:].mean())


def integrate(leg):
    w = [window(leg, i) for i in range(12)]
    have = [i for i, v in enumerate(w) if v]
    if len(have) < 12:
        return None, len(have)
    y = np.array([w[i]["mean"] for i in range(12)])
    se = np.array([w[i]["se"] for i in range(12)])
    y1 = np.array([w[i]["first"] for i in range(12)])
    y2 = np.array([w[i]["second"] for i in range(12)])
    Wa = np.array(W)
    dG = float(Wa @ y)
    xs = np.array([0.0] + LAM + [1.0])
    lo = y[0] - (y[1] - y[0]) * (LAM[0]) / (LAM[1] - LAM[0])
    hi = y[-1] + (y[-1] - y[-2]) * (1.0 - LAM[-1]) / (LAM[-1] - LAM[-2])
    quad = abs(float(CubicSpline(xs, np.concatenate([[lo], y, [hi]])
                                 ).integrate(0.0, 1.0)) - dG)
    return dict(dG=dG, dG2=float(Wa @ y2),
                stat=float(np.sqrt(((Wa * se) ** 2).sum())),
                conv=abs(float(Wa @ y2) - float(Wa @ y1)), quad=quad,
                tau_max=max(w[i]["tau"] for i in range(12)),
                neff_min=min(w[i]["n_eff"] for i in range(12)),
                w=w), 12


def main():
    res, status = {}, {}
    for leg in LEGS:
        r, n = integrate(leg)
        res[leg], status[leg] = r, n

    print("=" * 82)
    print("PYR1^MANDI QUADRUPLE (K59R/V81I/F108A/F159L) -- dG per leg, kcal/mol")
    print("=" * 82)
    print(f"  {'leg':<12}{'dG':>9}{'stat':>7}{'conv':>7}{'quad':>7}"
          f"{'dG(2nd half)':>14}{'tau_max':>9}{'n_eff_min':>11}")
    for leg in LEGS:
        r = res[leg]
        if r is None:
            print(f"  {leg:<12} incomplete ({status[leg]}/12 windows)")
            continue
        print(f"  {leg:<12}{r['dG']:9.2f}{r['stat']:7.2f}{r['conv']:7.2f}"
              f"{r['quad']:7.2f}{r['dG2']:14.2f}{r['tau_max']:9.1f}"
              f"{r['neff_min']:11.1f}")

    apo = res["aba"] and res["apo"]
    print()
    print("=" * 82)
    print("ddG_bind   (negative = the quadruple IMPROVES binding of that ligand)")
    print("=" * 82)
    for leg in ("aba", "mandi_xtal", "mandi_dock"):
        if not res[leg] or not res["apo"]:
            print(f"  {leg:<12} incomplete")
            continue
        dd = res[leg]["dG"] - res["apo"]["dG"]
        st = math.hypot(res[leg]["stat"], res["apo"]["stat"])
        cv = math.hypot(res[leg]["conv"], res["apo"]["conv"])
        qd = math.hypot(res[leg]["quad"], res["apo"]["quad"])
        print(f"  {leg:<12} {dd:8.2f}  +/- {math.sqrt(st*st+cv*cv+qd*qd):5.2f} total"
              f"   (stat {st:.2f} | conv {cv:.2f} | quad {qd:.2f})")

    print()
    print("=" * 82)
    print("SELECTIVITY = ddG(mandi) - ddG(ABA)   [apo cancels exactly]")
    print("negative = the quadruple shifts preference TOWARD mandipropamid")
    print("=" * 82)
    sels = {}
    for leg, lab in (("mandi_xtal", "crystal pose"), ("mandi_dock", "docked pose")):
        if not res[leg] or not res["aba"]:
            print(f"  {lab:<14} incomplete")
            continue
        sel = res[leg]["dG"] - res["aba"]["dG"]
        st = math.hypot(res[leg]["stat"], res["aba"]["stat"])
        cv = math.hypot(res[leg]["conv"], res["aba"]["conv"])
        qd = math.hypot(res[leg]["quad"], res["aba"]["quad"])
        tot = math.sqrt(st * st + cv * cv + qd * qd)
        sels[leg] = (sel, tot)
        verdict = ("prefers mandipropamid" if sel < -tot else
                   "prefers ABA" if sel > tot else "WITHIN ERROR -- no call")
        print(f"  {lab:<14} {sel:8.2f} +/- {tot:5.2f}   {verdict}")
        print(f"  {'':14} (stat {st:.2f} | conv {cv:.2f} | quad {qd:.2f})"
              f"   second-half-only {res[leg]['dG2'] - res['aba']['dG2']:7.2f}")

    if len(sels) == 2:
        gap = sels["mandi_dock"][0] - sels["mandi_xtal"][0]
        e = math.hypot(sels["mandi_dock"][1], sels["mandi_xtal"][1])
        print()
        print(f"  COST OF NOT HAVING A STRUCTURE: {gap:+.2f} +/- {e:.2f} kcal/mol")
        print(f"  (docked pose is 8.33 A heavy-atom RMSD from crystal; best of 9 = 4.90 A)")

    print()
    print("  PRE-REGISTERED TEST (13e): the CRYSTAL arm must come out NEGATIVE.")
    print("  Unlike 43d this is a real constraint -- PYR1^MANDI is the sequence")
    print("  the ground truth is about. If it fails, read conv and the per-window")
    print("  drift below BEFORE concluding anything about the method.")

    print()
    print("=" * 82)
    print("PER-WINDOW <dV/dlambda>  (watch for a sharp jump: ligand-pose hysteresis")
    print("between the WT-accommodated and quad-accommodated basins)")
    print("=" * 82)
    print(f"  {'lam':>8} " + "".join(f"{l:>18}" for l in LEGS))
    for i in range(12):
        row = f"  {LAM[i]:8.5f} "
        for leg in LEGS:
            r = res[leg]
            if r is None:
                row += f"{'--':>18}"
            else:
                w = r["w"][i]
                row += f"{w['mean']:12.2f}{w['second']-w['first']:+6.2f}"
        print(row)
    print("   (each cell: <dV/dl> then first->second half drift)")

    out = os.path.join(ROOT, "results", "ti_quad")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "windows.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["leg", "window", "lam", "mean", "se", "tau", "n_eff",
                     "first_half", "second_half", "drift"])
        for leg in LEGS:
            if res[leg] is None:
                continue
            for i in range(12):
                w = res[leg]["w"][i]
                wr.writerow([leg, i, LAM[i], f"{w['mean']:.4f}", f"{w['se']:.4f}",
                             f"{w['tau']:.2f}", f"{w['n_eff']:.1f}",
                             f"{w['first']:.4f}", f"{w['second']:.4f}",
                             f"{w['second']-w['first']:.4f}"])
    print(f"\nwrote {out}/windows.csv")


if __name__ == "__main__":
    main()
