#!/usr/bin/env python
r"""
88_ti_aggregate.py -- integrate <dV/dlambda> and report ddG_bind + selectivity.

INTEGRATION
    dG = sum_i w_i * <dV/dlambda>_i     (12-point Gauss-Legendre on [0,1])
The nodes never touch lambda = 0 or 1, so the softcore endpoint singularity is
avoided by construction rather than extrapolated away.

THE CYCLE
    ddG_bind(L) = dG(mutate in complex with L) - dG(mutate in apo)
    selectivity = ddG_bind(mandi) - ddG_bind(ABA)
                = dG(complex,mandi) - dG(complex,ABA)      <- apo cancels exactly

The apo leg cancels out of the selectivity, so that number is the more trustworthy
of the two: it depends on four fewer sources of error than the absolute ddG does.

UNCERTAINTY -- three terms, reported separately on purpose
  stat   correlation-corrected sd/sqrt(N_eff) from 89_ti_reparse.py. The first
         version of this script used sd/sqrt(n) with n=4004, which was wrong
         twice: n double-counted pmemd's two printed copies of every frame, and
         no autocorrelation correction was applied. It quoted +/-0.06 on numbers
         whose real statistical error is ~20x larger.
  conv   |dG(second halves) - dG(first halves)|. A within-run error bar cannot
         see a drift it is centred on -- that is exactly how MM-GBSA produced a
         confident wrong sign twice (README 40, 42). This term is the one that
         decides whether a window is converged, so it is never folded silently
         into 'stat'.
  quad   |spline integral - Gauss-Legendre sum|, i.e. whether 12 nodes actually
         resolve the shape of this integrand.

A result is only quotable when conv and quad are both small next to the effect.
"""
import csv
import math
import os

import numpy as np
from scipy.interpolate import CubicSpline

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(ROOT, "data", "ti", "ti_windows.csv")

LAM = [0.00922, 0.04794, 0.11505, 0.20634, 0.31608, 0.43738,
       0.56262, 0.68392, 0.79366, 0.88495, 0.95206, 0.99078]
W = [0.02359, 0.05347, 0.08004, 0.10158, 0.11675, 0.12457,
     0.12457, 0.11675, 0.10158, 0.08004, 0.05347, 0.02359]
LEGS = ["V81I_aba", "V81I_mandi", "V81I_apo",
        "K59R_aba", "K59R_mandi", "K59R_apo"]


def load():
    if not os.path.exists(CSV):
        raise SystemExit(f"{CSV} not found -- run scripts/89_ti_reparse.py first")
    d = {}
    for r in csv.DictReader(open(CSV)):
        d.setdefault(r["leg"], {})[int(r["window"])] = {
            k: float(v) for k, v in r.items() if k not in ("leg", "window")}
    return d


def integrate(w):
    """Gauss-Legendre sum plus the three error terms."""
    y = np.array([w[i]["mean"] for i in range(12)])
    se = np.array([w[i]["se"] for i in range(12)])
    y1 = np.array([w[i]["first_half"] for i in range(12)])
    y2 = np.array([w[i]["second_half"] for i in range(12)])
    Wa = np.array(W)

    dG = float(Wa @ y)
    stat = float(np.sqrt(((Wa * se) ** 2).sum()))
    conv = abs(float(Wa @ y2) - float(Wa @ y1))

    # does a 12-node rule resolve this integrand? compare with a spline that is
    # free to bend between the nodes. Endpoints are linearly extrapolated from
    # the two outermost windows -- softcore keeps dV/dl finite there.
    xs = np.array([0.0] + LAM + [1.0])
    lo = y[0] - (y[1] - y[0]) * (LAM[0] - 0.0) / (LAM[1] - LAM[0])
    hi = y[-1] + (y[-1] - y[-2]) * (1.0 - LAM[-1]) / (LAM[-1] - LAM[-2])
    ys = np.concatenate([[lo], y, [hi]])
    quad = abs(float(CubicSpline(xs, ys).integrate(0.0, 1.0)) - dG)

    dG2 = float(Wa @ y2)     # second-half-only estimate
    return dict(dG=dG, stat=stat, conv=conv, quad=quad, dG2=dG2,
                tau_max=max(w[i]["tau"] for i in range(12)),
                neff_min=min(w[i]["n_eff"] for i in range(12)))


def main():
    data = load()
    res = {leg: integrate(data[leg]) for leg in LEGS if len(data.get(leg, {})) == 12}

    print("=" * 78)
    print("THERMODYNAMIC INTEGRATION -- dG of mutation in each leg  (kcal/mol)")
    print("=" * 78)
    print(f"  {'leg':<13}{'dG':>8}{'stat':>7}{'conv':>7}{'quad':>7}"
          f"{'dG(2nd half)':>14}{'tau_max':>9}{'n_eff_min':>11}")
    for leg in LEGS:
        r = res.get(leg)
        if not r:
            print(f"  {leg:<13} incomplete")
            continue
        print(f"  {leg:<13}{r['dG']:8.2f}{r['stat']:7.2f}{r['conv']:7.2f}"
              f"{r['quad']:7.2f}{r['dG2']:14.2f}{r['tau_max']:9.1f}{r['neff_min']:11.1f}")

    print()
    print("=" * 78)
    print("ddG_bind  (negative = the mutation IMPROVES binding of that ligand)")
    print("=" * 78)
    for mut in ("V81I", "K59R"):
        apo = res.get(f"{mut}_apo")
        for lig in ("aba", "mandi"):
            cx = res.get(f"{mut}_{lig}")
            if not cx or not apo:
                print(f"  {mut} {lig:<6} incomplete")
                continue
            dd = cx["dG"] - apo["dG"]
            st = math.hypot(cx["stat"], apo["stat"])
            cv = math.hypot(cx["conv"], apo["conv"])
            qd = math.hypot(cx["quad"], apo["quad"])
            tot = math.sqrt(st ** 2 + cv ** 2 + qd ** 2)
            print(f"  {mut} {lig:<6} ddG_bind = {dd:7.2f}  +/- {tot:5.2f} total"
                  f"   (stat {st:.2f} | conv {cv:.2f} | quad {qd:.2f})")

    print()
    print("=" * 78)
    print("SELECTIVITY = ddG(mandi) - ddG(ABA)   [the apo leg cancels exactly]")
    print("negative = the mutation shifts preference TOWARD mandipropamid")
    print("=" * 78)
    for mut in ("V81I", "K59R"):
        a, m = res.get(f"{mut}_aba"), res.get(f"{mut}_mandi")
        if not a or not m:
            print(f"  {mut}: incomplete")
            continue
        sel = m["dG"] - a["dG"]
        sel2 = m["dG2"] - a["dG2"]
        st = math.hypot(a["stat"], m["stat"])
        cv = math.hypot(a["conv"], m["conv"])
        qd = math.hypot(a["quad"], m["quad"])
        tot = math.sqrt(st ** 2 + cv ** 2 + qd ** 2)
        verdict = ("prefers mandipropamid" if sel < -tot else
                   "prefers ABA" if sel > tot else "WITHIN ERROR -- no call")
        print(f"  {mut}: {sel:7.2f} +/- {tot:5.2f}   {verdict}")
        print(f"        (stat {st:.2f} | conv {cv:.2f} | quad {qd:.2f})"
              f"   second-half-only: {sel2:7.2f}")

    print()
    print("  Ground truth: both V81I and K59R are substitutions in the 4WVO")
    print("  mandipropamid sensor, so both SHOULD come out negative. MM-GBSA")
    print("  (README 40c) put them at +1.36 and +1.31 -- both the wrong sign,")
    print("  with a +/-3.85 uncertainty that could not resolve either.")


if __name__ == "__main__":
    main()
