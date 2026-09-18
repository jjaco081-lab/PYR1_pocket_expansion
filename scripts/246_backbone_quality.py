#!/usr/bin/env python
r"""
246_backbone_quality.py -- score RFd3 backbones BEFORE sequence design and pick
the pilot set.

Jannis: "a lot of secondary structure is good and obviously no long orange
residue stretches" -- and "Is there any program that scores backbones before
sequence design?"

THE ANSWER TO THAT QUESTION, because it shapes what this script can honestly do:
there is no validated pre-sequence "designability score". The field standard is
SELF-CONSISTENCY -- design a sequence, fold it, measure RMSD back to the design
backbone -- and that needs the sequence. The cheap proxies that DO precede it,
and that this script computes, are:

  * secondary-structure content (P-SEA, already in 143) -- loop-heavy backbones
    are hard to design and fold ambiguously;
  * ORANGE residues, i.e. kept-but-not-core from 220's flood fill -- Jannis's own
    criterion, and the thing his eye caught that the metrics missed;
  * compactness (Rg against PYR1's 15.0) and core fraction.

The one genuinely predictive cheap option is the **ProteinMPNN/LigandMPNN score
itself**: running MPNN and taking its log-probability, without committing to the
designed sequence, is a fast designability proxy that precedes folding. It is
noted here rather than used, because the pilot runs MPNN properly anyway and the
score comes out of that run for free.

⚠ ALL OF THESE ARE PROXIES. None of them says the backbone is foldable. That is
what the Boltz binary/ternary predictions in the pilot are for.
"""
import argparse, csv, glob, json, os, sys
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
imp = __import__("importlib").import_module
m143, m220 = imp("143_pocket_backbone_shape"), imp("220_core_truncate")
PYR1_RG, PYR1_CAV = 15.0, 164.4


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
        kept = sorted(k for k in main if lo <= k <= hi)
        orange = sorted(set(kept) - core)
        runs, cur = [], []
        for k in range(lo, hi + 1):
            if k in orange:
                cur.append(k)
            elif cur:
                runs.append(len(cur)); cur = []
        if cur:
            runs.append(len(cur))
        # `by` holds EVERY atom per residue, not CA -- psea and Rg need CA
        cad = {a[1]: a[4] for a in dat if a[2] == 'CA'}
        ca = np.array([cad[k] for k in kept if k in cad])
        if len(ca) < 10:
            return dict(arm=arm, design=design, error='too few CA')
        ss = m143.psea(ca)
        n = len(ss)
        cen = ca.mean(0)
        rg = float(np.sqrt(((ca - cen) ** 2).sum(1).mean()))
        return dict(arm=arm, design=design, n_kept=len(kept), n_core=len(core),
                    n_orange=len(orange),
                    longest_orange=max(runs) if runs else 0,
                    frac_H=round(ss.count("H") / n, 3),
                    frac_E=round(ss.count("E") / n, 3),
                    frac_SS=round((ss.count("H") + ss.count("E")) / n, 3),
                    rg=round(rg, 2),
                    core_frac=round(len(core) / max(len(kept), 1), 3))
    except Exception as e:                                       # noqa: BLE001
        return dict(arm=arm, design=design, error=f"{type(e).__name__}: {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", nargs="+",
                    default=[f"{ROOT}/results/rfd3/truncated.csv"])
    ap.add_argument("--max-orange", type=int, default=5)
    ap.add_argument("--min-ss", type=float, default=0.45)
    ap.add_argument("--procs", type=int, default=48)
    a = ap.parse_args()
    rows = []
    for c in a.csv:
        if os.path.exists(c):
            rows += [dict(r, _src=os.path.basename(c)) for r in
                     csv.DictReader(open(c))]
    passed = [r for r in rows if r.get("verdict_after") == "PASS"]
    print(f"{len(rows)} designs from {len(a.csv)} file(s), {len(passed)} PASS")
    cav = {}
    for r in passed:
        try: cav[(r["arm"], r["design"])] = float(r["cav_after"] or 0)
        except ValueError: pass
    with mp.Pool(a.procs) as p:
        res = p.map(one, [(r["arm"], r["design"]) for r in passed])
    ok = [r for r in res if "error" not in r]
    for r in ok:
        r["cavity"] = cav.get((r["arm"], r["design"]), 0.0)
        r["stage3"] = r["arm"].startswith(("arm9aba", "arm10ext"))
    print(f"{len(ok)} scored\n")
    import statistics as st
    print(f"secondary structure (H+E) across PASS designs: "
          f"median {st.median([r['frac_SS'] for r in ok]):.2f}")
    print(f"orange residues: median {st.median([r['n_orange'] for r in ok]):.0f}\n")
    good = [r for r in ok if r["n_orange"] <= a.max_orange
            and r["frac_SS"] >= a.min_ss]
    print(f"pass BOTH (<= {a.max_orange} orange AND >= {a.min_ss} SS): "
          f"{len(good)} of {len(ok)}")
    s3 = [r for r in good if r["stage3"]]
    old = [r for r in good if not r["stage3"]]
    print(f"   Stage 3 (ABA in the pocket): {len(s3)}")
    print(f"   earlier arms (empty pocket): {len(old)}\n")
    for lbl, s, k in (("STAGE 3 -- pick 3", s3, 3),
                      ("EARLIER ARMS, largest pockets -- pick 2", old, 2)):
        s.sort(key=lambda r: -r["cavity"])
        print(f"{lbl}")
        print(f"   {'arm':<20}{'design':<8}{'cav':>7}{'xP':>6}{'SS':>6}"
              f"{'H':>6}{'E':>6}{'orng':>6}{'Rg':>6}{'core':>6}")
        for r in s[:8]:
            m = "  <--" if s.index(r) < k else ""
            print(f"   {r['arm']:<20}{r['design']:<8}{r['cavity']:>7.1f}"
                  f"{r['cavity']/PYR1_CAV:>6.2f}{r['frac_SS']:>6.2f}"
                  f"{r['frac_H']:>6.2f}{r['frac_E']:>6.2f}{r['n_orange']:>6}"
                  f"{r['rg']:>6.1f}{r['core_frac']:>6.2f}{m}")
        print()
    sel = s3[:3] + old[:2]
    json.dump(dict(selected=sel, all_scored=ok),
              open(f"{ROOT}/results/rfd3/pilot_backbones.json", "w"), indent=1)
    print(f"selected {len(sel)} RFd3 designs -> results/rfd3/pilot_backbones.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
