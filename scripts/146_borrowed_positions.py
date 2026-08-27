#!/usr/bin/env python
r"""
146_borrowed_positions.py -- which PYR1 positions line a BIG pocket elsewhere?

Scripts 144 and 145 together rule out both ways of building a bigger pocket out
of new backbone:

  insertion    every large-cavity relative reaches its volume with AT MOST 3
               inserted pocket-lining residues (2PCS: 2, 2NS9: 2, 2BK0: 3,
               6AWV: 3). 43 of 2PCS's 45 wall residues sit at positions PYR1
               also has.
  displacement the sheet-to-grip-helix separation, measured at structurally
               equivalent positions, is 14.3 A in PYR1 and 14.7 A in 2PCS, and
               across 19 relatives correlates with cavity volume at r = -0.18.
               The lobes are not further apart.

What is left is the third possibility: the same backbone positions, with
different residues facing the cavity. 2PCS lines its pocket with 45 residues to
PYR1's 24, and only 18 of those are positions PYR1 also uses as wall. The other
~25 are positions PYR1 HAS, in the same fold, at the same place, which in CoxG
face the cavity and in PYR1 do not.

This script emits that list, per donor: the PYR1 residue number, its identity,
its secondary structure, the aligned donor residue, and how many donors
independently use it. A position nominated by several unrelated big-pocket
relatives is a stronger candidate than one nominated by a single structure.

This is a hypothesis generator, not a design. It says where the room is; it says
nothing about whether opening those positions costs the gate-latch-HAB1 switch,
which README 14/24/29 say is the risk that actually matters.
"""
import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m143 = import_module("143_pocket_backbone_shape")                   # noqa: E402
read_pdb, read_cif_resnames, THREE = m143.read_pdb, m143.read_cif_resnames, m143.THREE

HC = os.path.join(ROOT, "results", "homolog_cavities")
OUT = os.path.join(ROOT, "results", "pocket_shape")
BIG = 300.0            # A^3 -- donors whose pocket is meaningfully bigger than PYR1's
MIN_MATCH = 0.90


def main():
    shp = {r["name"]: r for r in json.load(open(os.path.join(OUT, "pocket_shape.json")))}
    pyr1 = shp["PYR1"]
    p_order = pyr1["ca_resnums"]
    p_ss = {n: pyr1["ss_string"][i] for i, n in enumerate(p_order)}
    p_wall = set(pyr1["lining"])
    frame = os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")
    p_seq = {a[0]: THREE[a[1]] for a in read_pdb(frame, want_atom=False)
             if a[1] in THREE}

    hits = {}
    for line in open(os.path.join(ROOT, "results", "foldseek", "fs_hits.tsv")):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 11:
            hits[f[1]] = f

    donors = sorted((n for n, r in shp.items()
                     if n != "PYR1" and r["cavity"] >= BIG),
                    key=lambda n: -shp[n]["cavity"])
    print(f"donors with cavity >= {BIG:.0f} A^3: "
          + ", ".join(f"{n} ({shp[n]['cavity']:.0f})" for n in donors))
    print(f"PYR1's own wall: {len(p_wall)} positions\n")

    votes, detail = Counter(), defaultdict(dict)
    for name in donors:
        tgt = next((t for t in hits if t.lower().startswith(name[:4].lower())
                    and t[4] == name[5]), None)
        if tgt is None:
            continue
        f = hits[tgt]
        qs, ts, qaln, taln = int(f[5]), int(f[7]), f[9], f[10]
        names = read_cif_resnames(os.path.join(HC, "raw", tgt + ".cif"), tgt[4])
        atoms = read_pdb(os.path.join(HC, "domains", tgt + ".pdb"))
        t_order = [a[0] for a in atoms if a[2] == "CA"]
        qi, ti, tmap, okt, tt = qs - 1, ts - 1, {}, 0, 0
        for qc, tc in zip(qaln, taln):
            if qc != "-":
                qi += 1
            if tc != "-":
                ti += 1
                if ti <= len(t_order):
                    tt += 1
                    okt += (names.get(t_order[ti - 1], "X") == tc)
            if qc != "-" and tc != "-":
                tmap[ti] = qi
        if not tt or okt / tt < MIN_MATCH:
            print(f"  {name}: index check FAILED -- dropped")
            continue
        lin = set(shp[name]["lining"])
        idx = {i: rn for i, rn in enumerate(t_order, 1)}
        new = []
        for i, rn in idx.items():
            if rn not in lin or i not in tmap:
                continue
            q = tmap[i]
            if q > len(p_order):
                continue
            pr = p_order[q - 1]
            if pr in p_wall:
                continue
            new.append(pr)
            votes[pr] += 1
            detail[pr][name] = names.get(rn, "X")
        print(f"  {name}: {len(new)} PYR1 positions that line ITS pocket "
              f"but not PYR1's")

    print("\n" + "=" * 78)
    print("BORROWED POSITIONS -- PYR1 numbering, ranked by independent donor support")
    print("=" * 78)
    print(f"{'PYR1':>5} {'aa':>3} {'ss':>3} {'donors':>7}   donor residues")
    for pr, v in sorted(votes.items(), key=lambda kv: (-kv[1], kv[0])):
        d = " ".join(f"{n.split('_')[0]}:{a}" for n, a in detail[pr].items())
        print(f"{pr:>5} {p_seq.get(pr,'?'):>3} {p_ss.get(pr,'?'):>3} {v:>7}   {d}")
    hi = [p for p, v in votes.items() if v >= 3]
    print(f"\n{len(votes)} distinct positions nominated; {len(hi)} by >= 3 of "
          f"{len(donors)} donors: {sorted(hi)}")
    ss = Counter(p_ss.get(p, "?") for p in hi)
    print(f"secondary structure of the {len(hi)}: "
          + ", ".join(f"{k}={ss[k]}" for k in "HEC" if ss[k]))
    json.dump({"votes": {str(k): v for k, v in votes.items()},
               "detail": {str(k): v for k, v in detail.items()},
               "consensus": sorted(hi)},
              open(os.path.join(OUT, "borrowed_positions.json"), "w"), indent=1)
    print(f"\nwritten to {OUT}/borrowed_positions.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
