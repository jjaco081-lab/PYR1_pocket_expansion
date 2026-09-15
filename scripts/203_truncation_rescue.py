#!/usr/bin/env python
r"""
203_truncation_rescue.py -- does the cavity SURVIVE truncating the HAB1-wrapping
lead?

Stage 1 established (n=120 per arm, anchor absent in both, so lead length is the
only variable):

  long lead 50-70  ->  90 % one intact chain, 50 % cavity
  short lead 10-25 ->  60 % intact          , 10 % cavity
                       Fisher p = 9.1e-08 and 8.7e-12

So shortening the lead -- my Stage-0 "fix" -- was the wrong move: it destroyed
both connectivity and the pocket. The anchor, which I thought was the culprit, is
irrelevant (p = 0.14 / 0.45).

But the long lead is exactly what wraps HAB1, and it is why 40 of 120 long-lead
designs fail F4. An audit of all 40 puts the appendage in the LEAD in 40 of 40
cases -- none internal, none trailing. Jannis said this from the start about
batch2 model 7: "a really long terminal disordered section that wraps around
HAB1. This can probably just be truncated."

THE QUESTION THIS ANSWERS. If the lead is truncated, F4 is satisfied by
construction -- but the lead may be part of the pocket WALL, in which case
truncation destroys the cavity it was holding and rescues nothing.

  cavity_full  = enclosed volume of the main piece, lead included
  cavity_trunc = same, with every non-motif residue before the first motif
                 residue deleted

RESCUED means cavity_trunc > 0 AND envelope_trunc > PYR1's 1184 A^3 AND the
truncated chain is still one piece. Reported as a rate with a Wilson interval.

⚠ Measured on the LARGEST CONNECTED PIECE only, per 198's F0 -- a cavity in the
gap between two pieces is not a pocket (d009: 61 -> 1 A^3 once the main piece was
isolated).
"""
import csv, glob, json, os, sys
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m192 = import_module("192_filter_cascade")
m195 = import_module("195_arm_compare")
m198 = import_module("198_cascade2")

PYR1_ENVELOPE = 1184.0
ARM = "arm3b3rasa_s3003"


def one(design):
    cif = glob.glob(os.path.join(ROOT, "results", "rfd3", ARM, design,
                                 "*.cif.gz"))
    if not cif:
        return None
    cif = cif[0]
    at = m195.read_cif_gz(cif)
    imap = json.load(open(cif.replace(".cif.gz", ".json")))["diffused_index_map"]
    dc = m195.design_chain(at, imap)
    dat = [a for a in at if a[0] == dc]
    ca = {a[1]: a[4] for a in dat if a[2] == "CA"}
    motif = {int(v[1:]) for v in imap.values() if v[0] == dc}
    main = max(m198.pieces_of(ca), key=len)
    keep = set(main)
    first = min(motif & keep) if (motif & keep) else min(keep)

    def measure(sel):
        mat = [(a[1], None, a[2], a[3], a[4]) for a in dat if a[1] in sel]
        if len(mat) < 50:
            return 0.0, 0.0, 0
        r, wall, _ = m192.cavity_and_wall(mat)
        env = m192.envelope(mat, wall) if wall else None
        sub = {k: ca[k] for k in sel if k in ca}
        npc = len(m198.pieces_of(sub)) if sub else 0
        return ((r or {}).get("main", 0.0),
                (env or {}).get("envelope", 0.0), npc)

    cf, ef, pf = measure(keep)
    trunc = {r for r in keep if r >= first}
    ct, et, pt = measure(trunc)
    return dict(design=design, n_lead=len(keep) - len(trunc),
                cav_full=round(cf, 1), env_full=round(ef, 1),
                cav_trunc=round(ct, 1), env_trunc=round(et, 1),
                pieces_trunc=pt,
                rescued=int(ct > 0 and et > PYR1_ENVELOPE and pt == 1))


def wilson(k, n, z=1.96):
    if not n:
        return 0.0, 0.0
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def main():
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    fails = [r["design"] for r in csv.DictReader(
        open(os.path.join(ROOT, "results", "rfd3", "cascade2.csv")))
        if r["arm"] == ARM and r["fail_at"] == "F4 HAB1 appendage"]
    print(f"{len(fails)} F4 failures in {ARM}; testing lead truncation "
          f"on {procs} procs", flush=True)
    with mp.Pool(procs) as pool:
        rows = [r for r in pool.map(one, fails) if r]
    print(f"\n{'design':<7}{'lead':>6}{'cav_full':>10}{'cav_trunc':>11}"
          f"{'env_full':>10}{'env_trunc':>11}{'pieces':>8}{'':>3}")
    for r in sorted(rows, key=lambda x: -x["env_trunc"]):
        print(f"{r['design']:<7}{r['n_lead']:>6}{r['cav_full']:>10.1f}"
              f"{r['cav_trunc']:>11.1f}{r['env_full']:>10.1f}"
              f"{r['env_trunc']:>11.1f}{r['pieces_trunc']:>8}"
              f"{'  RESCUED' if r['rescued'] else ''}")
    k = sum(r["rescued"] for r in rows)
    lo, hi = wilson(k, len(rows))
    print(f"\nRESCUED {k}/{len(rows)} ({100*k/max(len(rows),1):.0f}%, "
          f"95% CI {100*lo:.0f}-{100*hi:.0f})")
    kept = sum(1 for r in rows if r["cav_trunc"] > 0)
    print(f"  cavity survives truncation at all: {kept}/{len(rows)}")
    print(f"  median envelope {np.median([r['env_full'] for r in rows]):.0f} "
          f"-> {np.median([r['env_trunc'] for r in rows]):.0f} A^3")
    json.dump(rows, open(os.path.join(ROOT, "results", "rfd3",
                                      "truncation_rescue.json"), "w"), indent=1)
    print("\nwrote results/rfd3/truncation_rescue.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
