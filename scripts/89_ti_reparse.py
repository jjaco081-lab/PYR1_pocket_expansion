#!/usr/bin/env python
r"""
89_ti_reparse.py -- re-extract <dV/dlambda> from prod.out with an honest error bar.

WHY THIS EXISTS
Script 87 wrote each window's dvdl.txt with
    re.findall(r"DV/DL\s+=\s+(-?\d+\.\d+)", prod.out)
which is wrong in two ways that 88 then trusted:

  1. pmemd prints EVERY timestep TWICE -- one energy block per TI region, with a
     bit-identical DV/DL in each (verified: 0 of 2000 steps disagree). The regex
     collected both copies, so n was double-counted.
  2. It also swept up the trailing "A V E R A G E S" and "R M S  F L U C T U A -
     T I O N S" banners, whose DV/DL fields are a mean and a standard deviation,
     not samples. Those are the 4 extra entries in n=4004 = 2000*2 + 4.

Both defects are nearly invisible in the MEAN (36.8074 -> 36.8199 at K59R_aba
lam04) and fatal in the UNCERTAINTY, which is the number the pilot exists to
produce. 88's error was sd/sqrt(n) with n=4004; the truth is sd/sqrt(N_eff) with
N_eff = n/(1+2*tau) on 2000 correlated frames -- typically 30-80x smaller.

This is the third time in this project a tight within-run error bar has sat on
top of a mean that was still moving (README 40, 42). So this script reports, per
window: the correlation-corrected SE, a blocking SE as an independent check, the
integrated autocorrelation time, and the drift between the first and second half
of the run -- and the aggregate carries the drift as a SEPARATE error term,
because a within-run error bar cannot see a drift it is centred on.

Writes data/ti/ti_windows.csv. Read by 88_ti_aggregate.py.
"""
import csv
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TI = os.path.join(ROOT, "data", "ti")
LEGS = ["V81I_aba", "V81I_mandi", "V81I_apo",
        "K59R_aba", "K59R_mandi", "K59R_apo"]
LAM = [0.00922, 0.04794, 0.11505, 0.20634, 0.31608, 0.43738,
       0.56262, 0.68392, 0.79366, 0.88495, 0.95206, 0.99078]


def series(path):
    """Per-frame DV/DL, de-duplicated by NSTEP, averages banners excluded."""
    txt = open(path).read()
    body = re.split(r"A V E R A G E S", txt)[0]      # drop the trailing banners
    pairs = re.findall(r"NSTEP =\s*(\d+).*?DV/DL\s+=\s+(-?\d+\.\d+)", body, re.S)
    seen, order = {}, []
    for s, v in pairs:
        s = int(s)
        if s not in seen:
            seen[s] = float(v)
            order.append(s)
        elif abs(seen[s] - float(v)) > 1e-6:
            raise SystemExit(f"{path}: the two printed copies of step {s} differ "
                             f"({seen[s]} vs {v}) -- de-duplication is unsafe here")
    return np.array([seen[s] for s in sorted(order)])


def series_dir(wd):
    """Concatenate every production segment in a window directory, in order.

    105 extends production in chunks (prod.out, prod01.out, prod02.out, ...).
    De-duplication must happen WITHIN each file and never across them: Amber
    restarts NSTEP from zero in every new run, so a global de-duplication would
    silently discard all but the first segment.
    """
    import glob
    files = sorted(glob.glob(os.path.join(wd, "prod*.out")))
    if not files:
        return None
    return np.concatenate([series(f) for f in files])


def tau_int(x):
    """Integrated autocorrelation time, initial-positive-sequence estimator."""
    x = x - x.mean()
    n = len(x)
    if x.std() == 0:
        return 0.0
    c = np.correlate(x, x, "full")[n - 1:] / (np.arange(n, 0, -1) * x.var())
    t = 0.0
    for k in range(1, min(n // 2, 5000)):
        if c[k] <= 0:
            break
        t += c[k]
    return t


def block_se(x, nb=20):
    b = len(x) // nb
    if b < 2:
        return float("nan")
    m = np.array([x[i * b:(i + 1) * b].mean() for i in range(nb)])
    return m.std(ddof=1) / np.sqrt(nb)


def main():
    rows = []
    for leg in LEGS:
        for i, lam in enumerate(LAM):
            p = os.path.join(TI, leg, f"lam{i:02d}", "prod.out")
            if not os.path.exists(p):
                print(f"MISSING {leg} lam{i:02d}", file=sys.stderr)
                continue
            x = series(p)
            h = len(x) // 2
            t = tau_int(x)
            neff = len(x) / (1 + 2 * t)
            rows.append(dict(
                leg=leg, window=i, lam=lam, n=len(x),
                mean=x.mean(), sd=x.std(ddof=1),
                tau=t, n_eff=neff,
                se=x.std(ddof=1) / np.sqrt(max(neff, 1)),
                se_block=block_se(x),
                first_half=x[:h].mean(), second_half=x[h:].mean(),
                drift=x[h:].mean() - x[:h].mean(),
                mean_2nd_half=x[h:].mean(),
                se_2nd_half=x[h:].std(ddof=1) / np.sqrt(
                    max(len(x[h:]) / (1 + 2 * tau_int(x[h:])), 1)),
            ))
            print(f"  {leg:<12} lam{i:02d} n={len(x):5d} mean={x.mean():9.3f} "
                  f"tau={t:6.1f} n_eff={neff:7.1f} se={rows[-1]['se']:6.3f} "
                  f"(block {rows[-1]['se_block']:6.3f})  drift={rows[-1]['drift']:7.2f}")
    out = os.path.join(TI, "ti_windows.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}  ({len(rows)} windows)")


if __name__ == "__main__":
    main()
