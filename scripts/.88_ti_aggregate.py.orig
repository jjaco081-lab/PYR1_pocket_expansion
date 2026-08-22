#!/usr/bin/env python
"""
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
Both are reported, and they should be read with that asymmetry in mind.

CONVERGENCE
Every window carries its own first-half/second-half drift. This is the diagnostic
MM-GBSA never had -- README 40/42 record two separate occasions where a tight
within-run error bar sat on top of a mean that was still moving. A window whose
<dV/dl> drifts is flagged, and the integral it feeds is not quoted as converged.
"""
import glob
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TI = os.path.join(ROOT, "data", "ti")

LAM = [0.00922, 0.04794, 0.11505, 0.20634, 0.31608, 0.43738,
       0.56262, 0.68392, 0.79366, 0.88495, 0.95206, 0.99078]
W = [0.02359, 0.05347, 0.08004, 0.10158, 0.11675, 0.12457,
     0.12457, 0.11675, 0.10158, 0.08004, 0.05347, 0.02359]
LEGS = ["V81I_aba", "V81I_mandi", "V81I_apo",
        "K59R_aba", "K59R_mandi", "K59R_apo"]


def read(leg):
    out = []
    for i in range(12):
        f = os.path.join(TI, leg, f"lam{i:02d}", "dvdl.txt")
        if not os.path.exists(f):
            out.append(None)
            continue
        d = dict(l.strip().split("=") for l in open(f) if "=" in l)
        out.append({k: float(v) for k, v in d.items()})
    return out


def integrate(leg):
    w = read(leg)
    have = [i for i, x in enumerate(w) if x]
    if len(have) < 12:
        return None, w, f"{len(have)}/12 windows"
    dG = sum(W[i] * w[i]["mean"] for i in range(12))
    # propagate each window's own sd through the quadrature weights
    err = math.sqrt(sum((W[i] * w[i]["sd"] / math.sqrt(max(w[i]["n"], 1))) ** 2
                        for i in range(12)))
    return (dG, err), w, "complete"


def main():
    res = {}
    print("=" * 78)
    print("THERMODYNAMIC INTEGRATION -- dG of mutation in each leg")
    print("=" * 78)
    for leg in LEGS:
        val, w, status = integrate(leg)
        res[leg] = val
        if val is None:
            print(f"  {leg:<14} {status}")
            continue
        dG, err = val
        drift = [abs(x["second_half"] - x["first_half"]) for x in w if x]
        bad = sum(1 for d in drift if d > 1.0)
        flag = f"  <-- {bad}/12 windows drift >1 kcal/mol" if bad else ""
        print(f"  {leg:<14} dG = {dG:8.2f} +/- {err:4.2f} kcal/mol   "
              f"max window drift {max(drift):5.2f}{flag}")

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
            dd = cx[0] - apo[0]
            e = math.sqrt(cx[1] ** 2 + apo[1] ** 2)
            print(f"  {mut} {lig:<6} ddG_bind = {dd:7.2f} +/- {e:4.2f}")

    print()
    print("=" * 78)
    print("SELECTIVITY = ddG(mandi) - ddG(ABA)   [the apo leg cancels]")
    print("negative = the mutation shifts preference TOWARD mandipropamid")
    print("=" * 78)
    for mut in ("V81I", "K59R"):
        a, m = res.get(f"{mut}_aba"), res.get(f"{mut}_mandi")
        if not a or not m:
            print(f"  {mut}: incomplete")
            continue
        sel = m[0] - a[0]
        e = math.sqrt(a[1] ** 2 + m[1] ** 2)
        verdict = ("prefers mandipropamid" if sel < -e else
                   "prefers ABA" if sel > e else "within error")
        print(f"  {mut}: {sel:7.2f} +/- {e:4.2f}   {verdict}")

    print()
    print("  Ground truth: both V81I and K59R are substitutions in the 4WVO")
    print("  mandipropamid sensor, so both SHOULD come out negative. MM-GBSA")
    print("  (README 40c) put them at +1.36 and +1.31 -- both the wrong sign,")
    print("  with a +/-3.85 uncertainty that could not resolve either.")


if __name__ == "__main__":
    main()
