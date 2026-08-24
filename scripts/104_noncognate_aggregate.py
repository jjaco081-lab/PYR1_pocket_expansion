#!/usr/bin/env python
r"""
104_noncognate_aggregate.py -- can the closed-state filter rank LIGANDS?

THE QUESTION (README 27, 43a)
S24b licensed a stability filter on the closed state, but S1 was apo AND open
while S2 was holo AND closed, so conformation was confounded with occupancy. If
closed PYR1 holds its state just as well around a ligand it was never built for,
that filter reports CONFORMATION, not ligand fit, and cannot rank designs.

Three ABA-sized ligands (19-20 heavy atoms, MW 270-276) from three chemical
classes -- furanocoumarin, nitroaromatic anilide, steroid -- each from three
INDEPENDENTLY DOCKED poses, against S2 (ABA) as reference. Same receptor file
byte-for-byte, so the only difference is the ligand.

TWO PRE-REGISTERED ASYMMETRIES (S27c), neither of which may be over-read
  1. ABSENCE OF RELEASE IS NOT EVIDENCE OF BINDING. Residence times are us-ms;
     150 ns cannot sample unbinding. Release IF SEEN is informative; staying put
     is not.
  2. A DOCKED POSE IS A HYPOTHESIS. If a ligand arm differs from ABA the cause
     could be the pose rather than the chemistry -- which is why there are three
     poses and why the per-pose spread is printed, never just the mean.

THE STATISTICS (same discipline as 59)
  - the INFERENTIAL UNIT is the RUN (n=3 per ligand), never the frame. There are
    15,000 frames per run but loop and pose motions are correlated over ns, so a
    frame-pooled test would manufacture significance for anything.
  - integrated autocorrelation time is measured and n_eff reported beside every
    frame-pooled number
  - with n=3 vs n=3 no p-value is worth printing: effect size and the
    replicate-level RANGE are what get reported
  - a further asymmetry worth stating: S2's three replicates are three SEEDS of
    ONE pose, while each non-cognate's three are three DIFFERENT poses. The
    non-cognate spread therefore includes pose variance that ABA's does not, so
    ABA's range is the more optimistic of the two.

EQUILIBRATION DISCARD is chosen from the CORE backbone RMSD only -- never from
the observables being tested -- so it cannot be tuned to flatter the answer.
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results", "noncognate")

GROUPS = {
    "ABA (S2, cognate)":      [f"S2_holo_closed_rep{i}" for i in range(3)],
    "Imperatorin (S6)":       [f"S6_imperatorin_pose{i}" for i in range(3)],
    "Flutamide (S7)":         [f"S7_flutamide_pose{i}" for i in range(3)],
    "alpha-Estradiol (S8)":   [f"S8_estradiol_pose{i}" for i in range(3)],
}
DISCARD_NS = 30.0          # see below; frames are 10 ps apart
PS_PER_FRAME = 10.0


def col(unit, name, idx=1):
    p = os.path.join(RES, unit, f"{name}.dat")
    if not os.path.exists(p):
        return None
    v = []
    for l in open(p):
        if l.startswith("#"):
            continue
        f = l.split()
        if len(f) > idx:
            try:
                v.append(float(f[idx]))
            except ValueError:
                pass
    return np.array(v) if v else None


def tau_int(x):
    x = np.asarray(x, float)
    x = x - x.mean()
    n = len(x)
    if n < 10 or x.std() == 0:
        return 0.0
    c = np.correlate(x, x, "full")[n - 1:] / (np.arange(n, 0, -1) * x.var())
    t = 0.0
    for k in range(1, min(n // 2, 4000)):
        if c[k] <= 0:
            break
        t += c[k]
    return t


def summarise(unit, discard):
    out = {}
    for key in ("core_rmsd", "lig_rmsd", "gate_rmsd", "latch_rmsd", "lb7a5_rmsd"):
        v = col(unit, key)
        if v is None:
            return None
        v = v[discard:]
        out[key] = v
    c = col(unit, "lig_contacts")
    out["contacts"] = c[discard:] if c is not None else None
    return out


def exact_p(a, b):
    """Two-sided exact permutation test on RUN-LEVEL means.

    With n=3 vs n=3 there are only C(6,3) = 20 distinct splits, so the smallest
    attainable two-sided p is 2/20 = 0.10. THE DESIGN CANNOT PRODUCE A
    CONVENTIONALLY SIGNIFICANT RESULT no matter what the data say -- which is
    worth printing next to every comparison, because "the ranges do not overlap"
    reads as strong and is exactly this p = 0.10.
    """
    from itertools import combinations
    pool = list(a) + list(b)
    n = len(a)
    obs = abs(np.mean(a) - np.mean(b))
    cnt = tot = 0
    for idx in combinations(range(len(pool)), n):
        g1 = [pool[i] for i in idx]
        g2 = [pool[i] for i in range(len(pool)) if i not in idx]
        tot += 1
        if abs(np.mean(g1) - np.mean(g2)) >= obs - 1e-12:
            cnt += 1
    return cnt / tot, 2.0 / tot


def fmt_group(vals):
    """mean of the per-run means, and the run-level range."""
    a = np.array(vals, float)
    return f"{a.mean():5.2f}  [{a.min():.2f}-{a.max():.2f}]"


def main():
    missing = [u for g in GROUPS.values() for u in g
               if not os.path.exists(os.path.join(RES, u, "lig_rmsd.dat"))]
    if missing:
        print(f"missing extraction for: {missing}")
        print("run scripts/103_noncognate_extract.sh first")
        return 1

    discard = int(DISCARD_NS * 1000 / PS_PER_FRAME)

    print("=" * 84)
    print("1. VALIDITY -- did the systems stay themselves?")
    print("=" * 84)
    print(f"  discarding the first {DISCARD_NS:.0f} ns ({discard} frames), chosen "
          "from CORE backbone RMSD only")
    print(f"  {'unit':<26}{'core RMSD':>12}{'lig RMSD':>11}{'contacts':>11}"
          f"{'tau(lig)':>10}{'n_eff':>8}")
    data = {}
    for gname, units in GROUPS.items():
        for u in units:
            s = summarise(u, discard)
            if s is None:
                print(f"  {u:<26} MISSING")
                continue
            data[u] = s
            t = tau_int(s["lig_rmsd"])
            neff = len(s["lig_rmsd"]) / (1 + 2 * t)
            cm = s["contacts"].mean() if s["contacts"] is not None else float("nan")
            print(f"  {u:<26}{s['core_rmsd'].mean():12.2f}"
                  f"{s['lig_rmsd'].mean():11.2f}{cm:11.1f}{t:10.1f}{neff:8.0f}")

    print()
    print("=" * 84)
    print("2. THE READOUT -- ligand pose stability (A, RMSD to its own start)")
    print("   run-level means, with the range across the three runs")
    print("=" * 84)
    print(f"  {'ligand':<24}{'mean':>16}{'max reached':>20}{'final 20 ns':>18}")
    lig_means = {}
    for gname, units in GROUPS.items():
        ok = [u for u in units if u in data]
        if not ok:
            continue
        m = [data[u]["lig_rmsd"].mean() for u in ok]
        mx = [data[u]["lig_rmsd"].max() for u in ok]
        fin = [data[u]["lig_rmsd"][-2000:].mean() for u in ok]
        lig_means[gname] = m
        print(f"  {gname:<24}{fmt_group(m):>16}{fmt_group(mx):>20}{fmt_group(fin):>18}")

    if "ABA (S2, cognate)" in lig_means:
        aba = np.array(lig_means["ABA (S2, cognate)"])
        print()
        print("  effect size vs ABA (difference of run-level means; n=3 vs n=3,")
        print("  so this is an effect size and a range, NOT a p-value):")
        for gname, m in lig_means.items():
            if gname.startswith("ABA"):
                continue
            m = np.array(m)
            d = m.mean() - aba.mean()
            pooled = np.sqrt((m.var(ddof=1) + aba.var(ddof=1)) / 2)
            g = d / pooled if pooled > 0 else float("nan")
            overlap = not (m.min() > aba.max() or m.max() < aba.min())
            pv, floor = exact_p(m, aba)
            print(f"    {gname:<24} delta {d:+5.2f} A   Cohen d {g:+5.2f}   "
                  f"ranges {'OVERLAP' if overlap else 'SEPARATE'}   "
                  f"exact p {pv:.2f} (floor {floor:.2f})")

    print()
    print("=" * 84)
    print("3. THE S24b FILTER -- does gate/latch stability differ BY LIGAND?")
    print("   if it does not, the filter is conformation-reporting and cannot")
    print("   rank ligands, which disqualifies it as a design filter")
    print("=" * 84)
    print(f"  {'ligand':<24}{'gate RMSD':>16}{'latch RMSD':>18}{'Lb7-a5 RMSD':>18}")
    gate_means = {}
    for gname, units in GROUPS.items():
        ok = [u for u in units if u in data]
        if not ok:
            continue
        g = [data[u]["gate_rmsd"].mean() for u in ok]
        la = [data[u]["latch_rmsd"].mean() for u in ok]
        lb = [data[u]["lb7a5_rmsd"].mean() for u in ok]
        gate_means[gname] = g
        print(f"  {gname:<24}{fmt_group(g):>16}{fmt_group(la):>18}{fmt_group(lb):>18}")

    if "ABA (S2, cognate)" in gate_means:
        aba = np.array(gate_means["ABA (S2, cognate)"])
        print()
        for gname, m in gate_means.items():
            if gname.startswith("ABA"):
                continue
            m = np.array(m)
            overlap = not (m.min() > aba.max() or m.max() < aba.min())
            pv, floor = exact_p(m, aba)
            print(f"    {gname:<24} gate delta {m.mean()-aba.mean():+5.2f} A   "
                  f"ranges {'OVERLAP' if overlap else 'SEPARATE'}   "
                  f"exact p {pv:.2f} (floor {floor:.2f})")

    print()
    print("=" * 84)
    print("4. IS THE ANSWER AN ARTEFACT OF THE DISCARD WINDOW?")
    print("=" * 84)
    print(f"  {'discard':<10}" + "".join(f"{g.split(' (')[0]:>20}" for g in GROUPS))
    for dns in (0, 15, 30, 50, 75):
        d = int(dns * 1000 / PS_PER_FRAME)
        row = f"  {dns:>3d} ns   "
        for gname, units in GROUPS.items():
            ok = [u for u in units if u in data]
            v = [col(u, "lig_rmsd")[d:].mean() for u in ok]
            row += f"{np.mean(v):20.2f}"
        print(row)
    print("  (ligand RMSD, run-level mean; the ordering should not move)")

    print()
    print("  READ THIS WITH S27c IN MIND:")
    print("   - a ligand that STAYS PUT has not been shown to bind; 150 ns cannot")
    print("     sample unbinding. Only release, if seen, is informative.")
    print("   - ABA's three runs are three SEEDS of one pose; each non-cognate's")
    print("     are three DIFFERENT poses, so ABA's range is the optimistic one.")

    os.makedirs(os.path.join(ROOT, "results", "noncognate"), exist_ok=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
