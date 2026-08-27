#!/usr/bin/env python
r"""
113_us_pmf.py -- WHAM over the umbrella windows, and the calibration verdict.

Self-contained WHAM (no pymbar dependency, which is not installed here). The
iteration is the standard one:

    p(x)  = sum_i H_i(x) / sum_i N_i exp(-beta (w_i(x) - f_i))
    f_i   = -kT ln sum_x p(x) exp(-beta w_i(x))

with w_i(x) = k (x - x_i)^2 the harmonic bias, iterated to self-consistency.

WHAT IS BEING TESTED (README 59b) -- this is a CALIBRATION, not a result:

    WT + ABA  must favour CLOSED
    WT apo    must favour OPEN

If the two PMFs do not separate in that direction, the scheme is wrong and no
designed-pocket number from it is worth anything. The headline number is

    ddG = [G_open - G_closed]_holo  -  [G_open - G_closed]_apo

which must be POSITIVE (ABA stabilises closed relative to open). Reporting the
DIFFERENCE rather than either absolute PMF is deliberate: systematic error from
the choice of reaction coordinate largely cancels, which is the same reason TI
selectivity was usable when its absolute ddG was not (§43b, §46a).

THREE DIAGNOSTICS PRINTED ALONGSIDE, none of which may be skipped:
  overlap    neighbouring windows must share histogram support, or WHAM is
             interpolating across a gap it cannot see
  drift      first-half vs second-half PMF, the same convergence check §44a
             established after a tight error bar sat on a moving mean twice
  hysteresis windows seeded from the CLOSED basin vs the OPEN basin are compared
             in the region both reach. Systematic disagreement there means the
             coordinate is missing an orthogonal slow mode, which is the stated
             risk of this whole approach.
"""
import os
import sys
from collections import OrderedDict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
US = os.path.join(ROOT, "data", "umbrella")
KT = 0.001987204 * 300.0        # kcal/mol
K = 10.0                        # kcal/mol/A^2, must match 112
DISCARD = 0.25                  # fraction of each window's samples dropped


def load(tag, branch=""):
    """branch="" reads the forward run (w*/prod*.rc); branch="reverse" reads the
    reverse-seeded rerun (w*/reverse/prod*.rc). §138 keeps them side by side so the
    two can be compared without either overwriting the other."""
    out = OrderedDict()
    d = os.path.join(US, tag)
    if not os.path.isdir(d):
        return out
    import glob
    for w in sorted(os.listdir(d)):
        # production runs in chunks so each job fits short_gpu's 2 h limit, so a
        # window's samples are spread over prod01.rc, prod02.rc, ... Concatenate
        # in order; they are one continuous trajectory (irest=1 carries velocities
        # and box across, only the Langevin seed is redrawn).
        files = sorted(glob.glob(os.path.join(d, w, branch, "prod*.rc")))
        if not files:
            continue
        v = np.array([float(l.split()[1]) for f in files for l in open(f)
                      if l.strip() and not l.startswith("#")])
        if len(v) < 100:
            continue
        out[float(w[1:])] = v[int(len(v) * DISCARD):]
    return out


def wham(windows, bins, tol=1e-7, maxit=200000):
    centers = np.array(list(windows))
    data = list(windows.values())
    N = np.array([len(v) for v in data], float)
    edges = np.linspace(bins[0], bins[1], bins[2] + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    H = np.array([np.histogram(v, edges)[0] for v in data], float)
    bias = K * (mid[None, :] - centers[:, None]) ** 2          # w_i(x)
    expb = np.exp(-bias / KT)
    f = np.zeros(len(centers))
    for it in range(maxit):
        denom = (N[:, None] * expb * np.exp(f[:, None] / KT)).sum(0)
        p = H.sum(0) / np.where(denom > 0, denom, np.inf)
        s = p.sum()
        if s > 0:
            p = p / s
        fnew = -KT * np.log(np.where((expb * p[None, :]).sum(1) > 0,
                                     (expb * p[None, :]).sum(1), 1e-300))
        fnew -= fnew[0]
        if np.max(np.abs(fnew - f)) < tol:
            f = fnew
            break
        f = fnew
    with np.errstate(divide="ignore"):
        pmf = -KT * np.log(np.where(p > 0, p, np.nan))
    pmf -= np.nanmin(pmf)
    return mid, pmf, H, it


def basin(mid, pmf, lo, hi):
    m = (mid >= lo) & (mid <= hi) & np.isfinite(pmf)
    if not m.any():
        return np.nan
    return -KT * np.log(np.nansum(np.exp(-pmf[m] / KT)))


def report(tag, win):
    if not win:
        print(f"\n{tag}: no windows with data yet")
        return None
    mid, pmf, H, it = wham(win, (4.0, 21.0, 170))
    occ = (H > 0).sum(1)
    print(f"\n=== {tag}: {len(win)} windows, {int(sum(len(v) for v in win.values()))} "
          f"samples, WHAM converged in {it} iterations")
    thin = [f"{c:.1f}" for c, n in zip(win, occ) if n < 5]
    print(f"  windows with <5 occupied bins (overlap risk): "
          f"{thin if thin else 'none'}")
    gc, go = basin(mid, pmf, 4.5, 9.5), basin(mid, pmf, 13.0, 20.5)
    print(f"  G(closed, 4.5-9.5 A)  = {gc:7.2f} kcal/mol")
    print(f"  G(open, 13.0-20.5 A)  = {go:7.2f} kcal/mol")
    print(f"  G_open - G_closed     = {go-gc:+7.2f} kcal/mol   "
          f"({'CLOSED' if go > gc else 'OPEN'} favoured)")
    # convergence: half-split
    h1 = {c: v[:len(v)//2] for c, v in win.items()}
    h2 = {c: v[len(v)//2:] for c, v in win.items()}
    d1 = basin(*wham(h1, (4.0, 21.0, 170))[:2], 13.0, 20.5) - \
        basin(*wham(h1, (4.0, 21.0, 170))[:2], 4.5, 9.5)
    d2 = basin(*wham(h2, (4.0, 21.0, 170))[:2], 13.0, 20.5) - \
        basin(*wham(h2, (4.0, 21.0, 170))[:2], 4.5, 9.5)
    print(f"  half-split drift      : {d1:+.2f} -> {d2:+.2f}  (|d| = {abs(d2-d1):.2f})")
    return go - gc


def hysteresis(tag):
    """Compare the forward-seeded and reverse-seeded PMFs.

    §60 seeded from BOTH basins precisely so that "disagreement in the overlap
    region is visible"; §70's repair replaced every open seed with an outward pull
    from closed and destroyed that control; §136/§138 restored the missing
    direction by pulling INWARD from the equilibrated 20 A window. This is the
    check this file's docstring promised and never implemented.

    Reported on the windows both branches cover (10.5-19.5 A). A converged PMF is
    seed-independent; systematic divergence means the pull direction is setting the
    answer, and §81's verdict would then be about the protocol rather than the
    physics.
    """
    fwd, rev = load(tag), load(tag, "reverse")
    shared = [w for w in rev if w in fwd]
    if len(shared) < 5:
        print(f"\n{tag}: only {len(shared)} windows have BOTH branches "
              f"({len(rev)} reverse so far) -- rerun still in progress")
        return None
    print(f"\n=== HYSTERESIS {tag}: {len(shared)} windows with both seedings")
    mf, pf, _hf, _i = wham({w: fwd[w] for w in shared}, (10.0, 20.5, 105))
    mr, pr, _hr, _j = wham({w: rev[w] for w in shared}, (10.0, 20.5, 105))
    m = np.isfinite(pf) & np.isfinite(pr)
    d = pf[m] - pr[m]
    d = d - d.mean()                       # PMFs are defined up to a constant
    print(f"  mean |forward - reverse| after aligning the offset: "
          f"{np.abs(d).mean():.2f} kcal/mol   max {np.abs(d).max():.2f}")
    # the quantity that actually matters
    for nm, br in (("forward", fwd), ("reverse", rev)):
        mm, pp, _h, _k = wham({w: br[w] for w in shared}, (10.0, 20.5, 105))
        lo = basin(mm, pp, 10.5, 13.0)
        hi = basin(mm, pp, 16.0, 19.5)
        print(f"  {nm:<8} G(16-19.5) - G(10.5-13) = {hi-lo:+7.2f} kcal/mol")
    print("  READ: a converged PMF is seed-independent. If the two disagree by more")
    print("  than the half-split drift (§81a: 0.22 holo, 1.05 apo), the pull")
    print("  direction is setting the answer.")
    return float(np.abs(d).mean())


def main():
    holo = load("holo")
    apo = load("apo")
    dh = report("HOLO (WT + ABA)", holo)
    da = report("APO  (WT)", apo)
    for tag in ("holo", "apo"):
        hysteresis(tag)
    print("\n" + "=" * 72)
    print("CALIBRATION VERDICT (README 59b)")
    print("=" * 72)
    if dh is None or da is None or not np.isfinite(dh) or not np.isfinite(da):
        print("  incomplete -- both arms need windows with production data")
        return 1
    print(f"  holo  G_open - G_closed = {dh:+.2f} kcal/mol")
    print(f"  apo   G_open - G_closed = {da:+.2f} kcal/mol")
    print(f"  ddG (holo - apo)        = {dh-da:+.2f} kcal/mol")
    print()
    print("  PASS requires: holo favours CLOSED (positive), apo favours OPEN")
    print("  (negative), and ddG POSITIVE -- ABA stabilising closed vs open.")
    ok = dh > 0 and da < 0
    print(f"\n  ==> {'PASS' if ok else 'FAIL'}: "
          + ("the scheme reproduces the known answer and may be applied to a "
             "designed pocket" if ok else
             "the scheme does NOT reproduce the known answer; no designed-pocket "
             "number from it should be quoted (check overlap, drift, hysteresis "
             "above before blaming the physics)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
