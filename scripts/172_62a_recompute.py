#!/usr/bin/env python
r"""
172_62a_recompute.py -- §62a at n=3 and n=2 instead of n=1 and n=1.

§62a is the project's only positive result from unbiased MD: the two closed
cells have equal MEANS but different DISTRIBUTION WIDTHS, giving an effective
stiffness ratio and dF = -0.33 kcal/mol in favour of the holo closed state.

An external audit (README 95a) found it rested on ONE trajectory per cell,
because 110_factorial_extract.sh read only prod.nc and ignored every preemption
continuation. 11 of 12 replicates actually finished 300 ns. With the reaction
coordinate now measured on all of them, S9_apo_closed has n=3 and
S2_holo_closed n=2 usable (rep2 was killed at ~10 ns on a Blackwell node).

WHY THIS IS THE RIGHT TEST TO REPEAT. §61b dismissed the truncated replicates as
"error bars on a quantity that carries no information" -- written before §62
established that the information IS the distribution width. A width estimated
from a single trajectory has no error bar at all, and the whole claim is a
comparison of two widths.

k_eff = kT / var is the harmonic stiffness implied by a Gaussian of that width;
dF = -kT/2 * ln(k_apo / k_holo) is the free-energy difference between two
harmonic wells of those stiffnesses. Both are reported per replicate, so the
spread between replicates is visible rather than assumed away.
"""
import glob, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RC = os.path.join(ROOT, "results", "factorial", "rc")
KT = 0.001987 * 300.0          # kcal/mol at 300 K
MIN_FRAMES = 20000             # a 10 ns fragment is not a 300 ns width


def load(tag):
    f = os.path.join(RC, tag + ".dat")
    if not os.path.exists(f):
        return None
    v = np.array([float(l.split()[1]) for l in open(f)
                  if not l.startswith("#") and len(l.split()) > 1])
    return v


def main():
    cells = {"closed / apo": "S9_apo_closed", "closed / +ABA": "S2_holo_closed"}
    print(f"{'cell':<16}{'rep':>5}{'frames':>9}{'mean':>8}{'sd':>8}"
          f"{'p99':>8}{'max':>8}{'k_eff':>9}")
    stats = {}
    for label, sysname in cells.items():
        stats[label] = []
        for r in range(3):
            v = load(f"{sysname}_rep{r}")
            if v is None:
                continue
            keff = KT / v.var(ddof=1)
            flag = "" if len(v) >= MIN_FRAMES else "   <-- TRUNCATED, excluded"
            print(f"{label if r == 0 else '':<16}{r:>5}{len(v):>9}{v.mean():>8.2f}"
                  f"{v.std(ddof=1):>8.2f}{np.percentile(v, 99):>8.2f}"
                  f"{v.max():>8.2f}{keff:>9.2f}{flag}")
            if len(v) >= MIN_FRAMES:
                stats[label].append(dict(rep=r, n=len(v), mean=float(v.mean()),
                                         sd=float(v.std(ddof=1)),
                                         keff=float(keff)))
    print()
    a, h = stats["closed / apo"], stats["closed / +ABA"]
    print(f"usable replicates: apo n={len(a)}, holo n={len(h)}   "
          f"(README 62a used n=1 and n=1)")
    print(f"\n{'cell':<16}{'sd across reps':>26}{'k_eff across reps':>26}")
    for label, S in (("closed / apo", a), ("closed / +ABA", h)):
        sd = [x["sd"] for x in S]; kf = [x["keff"] for x in S]
        print(f"{label:<16}{np.mean(sd):>12.2f} +- {np.std(sd, ddof=1) if len(sd)>1 else float('nan'):<10.2f}"
              f"{np.mean(kf):>12.2f} +- {np.std(kf, ddof=1) if len(kf)>1 else float('nan'):<10.2f}")

    print("\nEVERY apo-vs-holo REPLICATE PAIRING (the test 62a ran once):")
    dfs = []
    for x in a:
        for y in h:
            df = -KT / 2 * np.log(x["keff"] / y["keff"])
            dfs.append(df)
            print(f"   apo rep{x['rep']} (k={x['keff']:.2f}) vs "
                  f"holo rep{y['rep']} (k={y['keff']:.2f})   "
                  f"dF = {df:+.3f} kcal/mol")
    dfs = np.array(dfs)
    print(f"\n   README 62a reported dF = -0.33 kcal/mol from ONE pairing.")
    print(f"   across all {len(dfs)} pairings: mean {dfs.mean():+.3f}, "
          f"range {dfs.min():+.3f} to {dfs.max():+.3f}, sd {dfs.std(ddof=1):.3f}")
    same = (dfs < 0).all() or (dfs > 0).all()
    print(f"   sign is {'CONSISTENT' if same else 'NOT consistent'} across pairings")
    json.dump({"stats": stats, "dF": dfs.tolist()},
              open(os.path.join(ROOT, "results", "factorial", "62a_recompute.json"), "w"),
              indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
