#!/usr/bin/env python
r"""
245_orange_filter.py -- drop designs with more than N ORANGE residues.

Jannis: "I want the RFd3 structures with a large internal loop or HAB1
interacting sequence away from the core to be filtered out. Make this any with
more than 5 orange residues for now."

ORANGE is defined by `220_core_truncate.write_pdb`: B-factor 0 = kept AND core,
**1 = kept but NOT core**, 2 = truncated. So orange is exactly `kept - core` --
residues that survive truncation but the flood-fill from the motif never reached.
Those are the dangling internal loops and the stretches that wander off to touch
HAB1 away from the body.

⚠ `truncated.csv` records only `internal_run`, the LONGEST CONSECUTIVE run. That
is not what is being asked for: five separate 1-residue excursions and one
5-residue loop both give internal_run <= 5 but differ completely. This counts the
TOTAL, which is what the eye sees as orange.

⚠ Counted on the TRUNCATED structure, not the raw one -- terminal excursions are
already removed by truncation and should not be charged twice.
"""
import argparse, csv, glob, json, os, sys
import multiprocessing as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
imp = __import__("importlib").import_module
m195, m198, m220 = imp("195_arm_compare"), imp("198_cascade2"), imp("220_core_truncate")


def one(args):
    arm, design = args
    g = glob.glob(f"{ROOT}/results/rfd3/{arm}/{design}/*.cif.gz")
    if not g:
        return dict(arm=arm, design=design, error="no structure")
    try:
        # ⚠ load_design returns SIX values (at, dat, dc, by, main, motif);
        # unpacking five silently failed on all 655 designs.
        at, dat, dc, by, main, motif = m220.load_design(g[0])
        core = m220.core_of(by, main, motif)
        if not core:
            return dict(arm=arm, design=design, error="no core")
        lo, hi = min(core), max(core)
        kept = {k for k in main if lo <= k <= hi}
        orange = sorted(kept - core)
        runs, cur = [], []
        for k in range(lo, hi + 1):
            if k in orange:
                cur.append(k)
            elif cur:
                runs.append(len(cur)); cur = []
        if cur:
            runs.append(len(cur))
        return dict(arm=arm, design=design, n_kept=len(kept), n_core=len(core),
                    n_orange=len(orange), longest_run=max(runs) if runs else 0,
                    n_runs=len(runs))
    except Exception as e:                                       # noqa: BLE001
        return dict(arm=arm, design=design, error=f"{type(e).__name__}: {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-orange", type=int, default=5)
    ap.add_argument("--procs", type=int, default=48)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(f"{ROOT}/results/rfd3/truncated.csv")))
    passed = [r for r in rows if r["verdict_after"] == "PASS"]
    print(f"{len(rows)} designs, {len(passed)} PASS after truncation")
    jobs = [(r["arm"], r["design"]) for r in passed]
    with mp.Pool(a.procs) as p:
        res = p.map(one, jobs)
    ok = [r for r in res if "error" not in r]
    print(f"{len(ok)} measured, {len(res)-len(ok)} failed\n")
    cav = {(r["arm"], r["design"]): float(r["cav_after"] or 0) for r in passed}
    for r in ok:
        r["cavity"] = cav.get((r["arm"], r["design"]), 0.0)
    import statistics as st
    o = [r["n_orange"] for r in ok]
    print(f"orange residues per PASS design: median {st.median(o):.0f}  "
          f"mean {st.mean(o):.1f}  max {max(o)}")
    print(f"  (longest RUN would have said: median "
          f"{st.median([r['longest_run'] for r in ok]):.0f})\n")
    for cut in (3, 5, 8, 10):
        keep = [r for r in ok if r["n_orange"] <= cut]
        big = [r for r in keep if r["cavity"] > 164.4]
        print(f"  <= {cut:>2} orange: {len(keep):>4} of {len(ok)} kept "
              f"({100*len(keep)/len(ok):>2.0f} %)   of which cavity > PYR1: {len(big)}")
    keep = [r for r in ok if r["n_orange"] <= a.max_orange]
    keep.sort(key=lambda r: -r["cavity"])
    out = os.path.join(ROOT, "results", "rfd3", "orange_filtered.json")
    json.dump(keep, open(out, "w"), indent=1)
    print(f"\nTOP 15 SURVIVING <= {a.max_orange} ORANGE:")
    print(f"{'arm':<20}{'design':<8}{'cavity':>8}{'xPYR1':>7}{'orange':>8}"
          f"{'runs':>6}{'core':>6}")
    for r in keep[:15]:
        print(f"{r['arm']:<20}{r['design']:<8}{r['cavity']:>8.1f}"
              f"{r['cavity']/164.4:>7.2f}{r['n_orange']:>8}{r['n_runs']:>6}"
              f"{r['n_core']:>6}")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
