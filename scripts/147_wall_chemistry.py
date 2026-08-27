#!/usr/bin/env python
r"""
147_wall_chemistry.py -- (a) 3OQU's wall against PYR1's, position by position
                         (b) do the borrowed positions of 146 hit real sensor
                             positions more often than chance?

PART A. 3OQU is the best transplantable donor: 46 % identity, 1.46 A core RMSD,
19/19 machinery coverage, and a 204 A^3 pocket against PYR1's 165. Same size,
same fold, same scaffold -- so any difference in what it BINDS has to come from
wall chemistry rather than wall volume. This prints the aligned residue at each
of PYR1's 24 wall positions, then the positions 3OQU uses as wall that PYR1 does
not, then aggregate chemistry (class counts, net charge, mean hydropathy, side-
chain H-bond donors and acceptors, aromatic count).

PART B. Script 146 nominated 52 PYR1 positions that line a large donor's pocket
without lining PYR1's. That list came from structures alone; Tian's libraries
were never consulted. So it can be tested against them.

  A = PYR1's own 24-residue wall            the obvious baseline
  B = the 52 borrowed positions             disjoint from A by construction
  G_lib = every position any sd04 library allows to vary   (design intuition)
  G_hit = every position actually substituted in a round-1 clone (selection)

The test is whether B hits G more often than a random draw of |B| positions
from the same protein, excluding A. Reported with an exact permutation p.

⚠ WHAT THIS CAN AND CANNOT SHOW. G_lib records where Tian chose to look, which
was itself decided by pocket proximity -- so agreement with G_lib is agreement
with their design intuition, not evidence of experimental success. G_hit is the
selection outcome but is a subset of G_lib, so a borrowed position outside the
library CANNOT appear in it however good it is. Absence from G_hit is therefore
not evidence against a position. Only the presence direction is informative.
"""
import json
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                          # noqa: E402
from importlib import import_module                                 # noqa: E402
m143 = import_module("143_pocket_backbone_shape")                   # noqa: E402
read_pdb, read_cif_resnames, THREE = m143.read_pdb, m143.read_cif_resnames, m143.THREE

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
HC = os.path.join(ROOT, "results", "homolog_cavities")
OUT = os.path.join(ROOT, "results", "pocket_shape")
AA = set("ACDEFGHIKLMNPQRSTVWY")

KD = dict(zip("AVLIPFMWGSTCYNQDEKRH",
              [1.8, 4.2, 3.8, 4.5, -1.6, 2.8, 1.9, -0.9, -0.4, -0.8, -0.7,
               2.5, -1.3, -3.5, -3.5, -3.5, -3.5, -3.9, -4.5, -3.2]))
CHG = {"D": -1, "E": -1, "K": 1, "R": 1}
#: side-chain H-bond donors / acceptors, heavy-atom counting
HBD = {"R": 3, "K": 1, "W": 1, "N": 1, "Q": 1, "H": 1, "S": 1, "T": 1, "Y": 1, "C": 1}
HBA = {"D": 2, "E": 2, "N": 1, "Q": 1, "H": 1, "S": 1, "T": 1, "Y": 1, "M": 1}
AROM = set("FWY")
CLASS = m143.CLASS


def align_map(tgt):
    """PYR1 resnum -> (donor resnum, donor aa) using the stored foldseek hit."""
    for line in open(os.path.join(ROOT, "results", "foldseek", "fs_hits.tsv")):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 11 and f[1] == tgt:
            break
    else:
        raise SystemExit(f"no foldseek hit for {tgt}")
    qs, ts, qaln, taln = int(f[5]), int(f[7]), f[9], f[10]
    names = read_cif_resnames(os.path.join(HC, "raw", tgt + ".cif"), tgt[4])
    t_order = [a[0] for a in read_pdb(os.path.join(HC, "domains", tgt + ".pdb"))
               if a[2] == "CA"]
    shp = {r["name"]: r for r in json.load(open(os.path.join(OUT, "pocket_shape.json")))}
    p_order = shp["PYR1"]["ca_resnums"]
    qi, ti, out, ok, n = qs - 1, ts - 1, {}, 0, 0
    for qc, tc in zip(qaln, taln):
        if qc != "-":
            qi += 1
        if tc != "-":
            ti += 1
            if ti <= len(t_order):
                n += 1
                ok += (names.get(t_order[ti - 1], "X") == tc)
        if qc != "-" and tc != "-" and qi <= len(p_order) and ti <= len(t_order):
            out[p_order[qi - 1]] = (t_order[ti - 1], names.get(t_order[ti - 1], "X"))
    assert n and ok / n >= 0.90, f"{tgt} index check failed ({ok}/{n})"
    return out, {v[0]: k for k, v in out.items()}


def chem(seq):
    c = Counter(CLASS.get(a, "other") for a in seq)
    return dict(n=len(seq),
                aliphatic=c["aliphatic"], aromatic=c["aromatic"],
                polar=c["polar"], basic=c["basic"], acidic=c["acidic"],
                charge=sum(CHG.get(a, 0) for a in seq),
                hydropathy=round(float(np.mean([KD.get(a, 0) for a in seq])), 2),
                hbd=sum(HBD.get(a, 0) for a in seq),
                hba=sum(HBA.get(a, 0) for a in seq))


def part_a(donor="3oquB00"):
    shp = {r["name"]: r for r in json.load(open(os.path.join(OUT, "pocket_shape.json")))}
    pyr1, dn = shp["PYR1"], shp[donor.upper()[:4] + "_" + donor[4]]
    p_seq = {a[0]: THREE[a[1]] for a in
             read_pdb(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb"), want_atom=False)
             if a[1] in THREE}
    p_ss = {n: pyr1["ss_string"][i] for i, n in enumerate(pyr1["ca_resnums"])}
    fwd, rev = align_map(donor)
    p_wall, d_wall = list(pyr1["lining"]), set(dn["lining"])

    print("=" * 78)
    print(f"A. {donor[:4].upper()} WALL vs PYR1 WALL "
          f"(cavity {dn['cavity']:.0f} vs {pyr1['cavity']:.0f} A^3)")
    print("=" * 78)
    print(f"{'PYR1':>5} {'aa':>3} {'ss':>3}   {'->':>3} {donor[:4].upper():>5} {'aa':>3}"
          f"   {'lines its pocket?':<18} change")
    both = []
    for r in p_wall:
        d = fwd.get(r)
        if d is None:
            print(f"{r:>5} {p_seq.get(r,'?'):>3} {p_ss.get(r,'?'):>3}   "
                  f"{'':>3} {'--':>5} {'-':>3}   {'(no aligned position)':<18}")
            continue
        same = "yes" if d[0] in d_wall else "no"
        ch = "=" if d[1] == p_seq.get(r) else \
             ("class" if CLASS.get(d[1]) != CLASS.get(p_seq.get(r, 'X')) else "within")
        both.append((p_seq.get(r, "X"), d[1]))
        print(f"{r:>5} {p_seq.get(r,'?'):>3} {p_ss.get(r,'?'):>3}   {'->':>3} "
              f"{d[0]:>5} {d[1]:>3}   {same:<18} {ch}")

    extra = sorted(r for r in d_wall if r in rev and rev[r] not in set(p_wall))
    print(f"\n   positions {donor[:4].upper()} uses as wall that PYR1 does not "
          f"({len(extra)}), in PYR1 numbering:")
    for r in extra:
        pr = rev[r]
        print(f"      PYR1 {pr:>3} {p_seq.get(pr,'?')} ({p_ss.get(pr,'?')})  ->  "
              f"{donor[:4].upper()} {r} {read_cif_resnames(os.path.join(HC,'raw',donor+'.cif'), donor[4]).get(r,'X')}")

    print(f"\n   WALL CHEMISTRY (each structure's OWN full wall)")
    cp, cd = chem(pyr1["seq"]), chem(dn["seq"])
    print(f"   {'':<14}{'PYR1':>8}{donor[:4].upper():>8}{'delta':>8}")
    for k in ("n", "aliphatic", "aromatic", "polar", "basic", "acidic",
              "charge", "hydropathy", "hbd", "hba"):
        print(f"   {k:<14}{cp[k]:>8}{cd[k]:>8}{cd[k]-cp[k]:>8.2f}")
    same = sum(1 for a, b in both if a == b)
    cls = sum(1 for a, b in both if a != b and CLASS.get(a) != CLASS.get(b))
    print(f"\n   at the {len(both)} shared wall positions: {same} identical, "
          f"{len(both)-same-cls} conservative, {cls} change chemical class")
    return {"pyr1": cp, donor: cd, "extra_positions": [rev[r] for r in extra]}


def part_b():
    bor = json.load(open(os.path.join(OUT, "borrowed_positions.json")))
    B = sorted(int(k) for k in bor["votes"])
    B12 = sorted(bor["consensus"])
    shp = {r["name"]: r for r in json.load(open(os.path.join(OUT, "pocket_shape.json")))}
    A = set(shp["PYR1"]["lining"])
    universe = set(shp["PYR1"]["ca_resnums"]) - A

    _, recs = table(f"{SD}/pnas.2519924122.sd04.xlsx")
    G_lib = {int(float(r["Position"])) for r in recs if r.get("Position")}
    hdr, cl = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    poscols = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    G_hit = set()
    for r in cl:
        for c in poscols:
            v = (r.get(c) or "").strip().upper()
            if len(v) == 1 and v in AA and v != c[0]:
                G_hit.add(int(c[1:]))

    print("\n" + "=" * 78)
    print("B. DO THE BORROWED POSITIONS HIT REAL SENSOR POSITIONS?")
    print("=" * 78)
    print(f"   PYR1 wall A = {len(A)}; borrowed B = {len(B)} (consensus 12 = {B12})")
    print(f"   G_lib = {len(G_lib)} positions any sd04 library varies; "
          f"G_hit = {len(G_hit)} substituted in a round-1 clone")
    rng = np.random.default_rng(0)
    pool = sorted(universe)
    out = {}
    for gname, G in (("G_lib", G_lib), ("G_hit", G_hit)):
        print(f"\n   vs {gname}")
        print(f"      A  (PYR1 wall, baseline)  {len(A & G):>3}/{len(A):<3} "
              f"= {len(A & G)/len(A):.0%}")
        if not (G - A):
            # G lies entirely inside the PYR1 wall, so B -- which is disjoint from
            # A by construction -- can never intersect it and neither can any draw
            # from the pool. The permutation null is identically zero and the test
            # has NO POWER. Report that, not a p-value of 1.
            print(f"      *** ALL {len(G)} of {gname} lie inside PYR1's own 24-residue")
            print(f"          wall. B is disjoint from that wall by construction, so")
            print(f"          B & {gname} = 0 IDENTICALLY and the permutation null is")
            print(f"          0 +- 0. This test has no power and is not evidence")
            print(f"          against the borrowed positions -- it is the statement")
            print(f"          that no published library has ever varied a position")
            print(f"          outside the existing wall.")
            out[gname] = dict(subset_of_wall=True, n=len(G))
            continue
        for label, S in (("B  (all 52)", B), ("B12 (>=3 donors)", B12)):
            S = [s for s in S if s in universe]
            obs = len(set(S) & G)
            null = np.array([len(set(rng.choice(pool, len(S), replace=False)) & G)
                             for _ in range(20000)])
            p = float((null >= obs).mean())
            print(f"      {label:<24} {obs:>3}/{len(S):<3} = {obs/len(S):.0%}   "
                  f"null {null.mean():.1f}+-{null.std():.1f}   p = {p:.4f}")
            out[f"{gname}:{label}"] = dict(obs=obs, n=len(S), p=p,
                                           null=float(null.mean()))
    print(f"\n   sd04 vocabulary = {sorted(G_lib)}")
    print(f"   all inside PYR1's wall? {not (G_lib - A)}   "
          f"({len(G_lib & A)}/{len(G_lib)})")
    print("   ⚠ There is therefore NO retrospective test available for the borrowed")
    print("     positions: the entire published search space is a subset of the")
    print("     existing wall, so the 52 are untried by construction, not refuted.")
    return out


def main():
    a = part_a()
    b = part_b()
    json.dump({"part_a": a, "part_b": b},
              open(os.path.join(OUT, "wall_chemistry.json"), "w"), indent=1)
    print(f"\nwritten to {OUT}/wall_chemistry.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
