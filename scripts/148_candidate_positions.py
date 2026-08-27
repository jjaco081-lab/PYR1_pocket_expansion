#!/usr/bin/env python
r"""
148_candidate_positions.py -- filter 146's 52 borrowed positions into a testable
set, and check that the donor consensus is not an artefact of one donor.

Three filters, applied in this order and each able to kill a position outright:

  ADJACENCY   distance from the position's C-beta to the nearest ABA heavy atom
              in the closed holo frame. This is the filter that matters most and
              it removes most of the list: the large donors' pockets are up to
              21.5 A long against PYR1's 11.6, so much of their wall sits where
              PYR1 has solid core. Opening a position 15 A from the ligand does
              not extend the pocket, it excavates a separate void.
  CONFLICT    membership of the gate, the latch, or the HAB1 interface (the four
              segments of README 9c). Those residues are the transduction
              machinery; a position that enlarges the pocket by breaking the
              switch has not helped.
  SUPPORT     number of independent donors nominating it.

LEAVE-ONE-DONOR-OUT. With only four donors, "consensus" could mean one dominant
structure plus noise. So each donor is held out in turn and the >= 2-vote
consensus of the other three is scored against it, with a hypergeometric null
over the 155 non-wall positions. This is the internal-replicate check that
caught the capsaicinoid subgroup in README 82.

⚠ The four donors are all helix-grip relatives and share ancestry, so their
agreement is not four independent observations. The LODO p-values are against a
random-position null, not against a phylogenetic one, and are correspondingly
optimistic.
"""
import json
import os
import sys

import numpy as np
from scipy.stats import hypergeom

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m143 = import_module("143_pocket_backbone_shape")                   # noqa: E402
read_pdb, THREE = m143.read_pdb, m143.THREE

OUT = os.path.join(ROOT, "results", "pocket_shape")
SEG = {"helix-60s": [60, 61, 63], "gate": [84, 85, 86, 87, 88, 89],
       "latch": [116, 117],
       "C-lobe": [148, 151, 155, 156, 158, 159, 162, 166]}
IFACE = {p: k for k, v in SEG.items() for p in v}
NEAR = 9.0            # A from ABA -- adjacency cut, fixed before looking


def main():
    bor = json.load(open(os.path.join(OUT, "borrowed_positions.json")))
    det = bor["detail"]
    shp = {r["name"]: r for r in json.load(open(os.path.join(OUT, "pocket_shape.json")))}
    wall = set(shp["PYR1"]["lining"])
    pool = set(shp["PYR1"]["ca_resnums"]) - wall
    N = len(pool)

    at = read_pdb(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb"), want_atom=False)
    aba = np.array([a[4] for a in at if a[1] == "A8S"])
    assert len(aba), "no A8S in the frame"
    by = {}
    for a in at:
        if a[1] in THREE:
            by.setdefault(a[0], {})[a[2]] = (a[4], THREE[a[1]])

    rows = []
    for p, dv in det.items():
        p = int(p)
        d = by.get(p)
        if not d:
            continue
        c = d.get("CB", d.get("CA"))
        rows.append(dict(pos=p, aa=c[1], votes=len(dv),
                         d_aba=round(float(np.linalg.norm(aba - c[0], axis=1).min()), 1),
                         conflict=IFACE.get(p, ""),
                         donors={k: v for k, v in dv.items()}))
    rows.sort(key=lambda r: (r["d_aba"], -r["votes"]))

    keep = [r for r in rows if r["d_aba"] <= NEAR and not r["conflict"]
            and r["votes"] >= 2]
    print("=" * 84)
    print(f"CANDIDATE POSITIONS -- borrowed, adjacent (<= {NEAR:.0f} A from ABA), "
          f"no switch conflict, >= 2 donors")
    print("=" * 84)
    print(f"{'PYR1':>5}{'aa':>4}{'votes':>7}{'d(ABA)':>9}   donor residues")
    for r in keep:
        print(f"{r['pos']:>5}{r['aa']:>4}{r['votes']:>7}{r['d_aba']:>9.1f}   "
              + " ".join(f"{k.split('_')[0]}:{v}" for k, v in r['donors'].items()))
    print(f"\n{len(keep)} of {len(rows)} survive. Removed: "
          f"{sum(1 for r in rows if r['d_aba'] > NEAR)} too far from the ligand, "
          f"{sum(1 for r in rows if r['conflict'])} on gate/latch/HAB1 "
          f"({sorted(r['pos'] for r in rows if r['conflict'])}), "
          f"{sum(1 for r in rows if r['votes'] < 2 and r['d_aba'] <= NEAR)} single-donor.")

    # ---------------------------------------------------------------- LODO ---
    sets = {d: {int(p) for p, v in det.items() if d in v}
            for d in sorted({d for v in det.values() for d in v})}
    print(f"\nLEAVE-ONE-DONOR-OUT over {N} non-wall positions")
    lodo = {}
    for d in sets:
        others = [e for e in sets if e != d]
        cons = {p for p in set().union(*[sets[e] for e in others])
                if sum(p in sets[e] for e in others) >= 2}
        k = len(cons & sets[d])
        exp = len(cons) * len(sets[d]) / N
        pv = float(hypergeom.sf(k - 1, N, len(cons), len(sets[d])))
        lodo[d] = dict(consensus=len(cons), recovered=k, expected=round(exp, 1), p=pv)
        print(f"   held out {d:<9} {len(cons):>2} consensus positions, "
              f"recovered {k:>2} ({k/max(len(cons),1):.0%})  "
              f"expected {exp:.1f}  p = {pv:.1e}")
    print("   ⚠ shared ancestry: these are 4 helix-grip relatives, not 4 "
          "independent samples.")

    json.dump({"ranked": rows, "candidates": keep, "lodo": lodo,
               "near_cut": NEAR},
              open(os.path.join(OUT, "candidate_positions.json"), "w"), indent=1)
    print(f"\nwritten to {OUT}/candidate_positions.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
