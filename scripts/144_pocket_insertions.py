#!/usr/bin/env python
r"""
144_pocket_insertions.py -- WHERE does a bigger relative put its extra wall?

Script 143 says how big the envelope is. This says what you would have to BUILD
to get it: for every candidate, each residue of its pocket wall is classified as

  ALIGNED    it sits at a position PYR1 also has -- reachable by MUTATION
  INSERT     it sits in a stretch the target has and PYR1 does not -- reachable
             only by adding backbone (loop/strand insertion or a graft)
  UNALIGNED  outside the alignment altogether (terminus or a region foldseek
             could not match) -- effectively a rebuild

and every insertion is reported with the PYR1 position it would be spliced
into, its length, PYR1's local secondary structure there, and how many of its
residues actually line the pocket. An insertion that does not touch the pocket
is not worth the stability risk; an insertion in the middle of a beta strand is
a register shift, not a loop swap, and is a different (harder) proposition than
one in a hairpin turn.

INDEX CONVENTION, ASSERTED NOT ASSUMED
--------------------------------------
Script 79 takes the foldseek target index to be the 1-based rank of the residue
in results/homolog_cavities/domains/<target>.pdb. That convention is checked
here by reconstructing the aligned target sequence from the domain's own residue
names and comparing it to `taln`; a target whose reconstruction disagrees is
dropped rather than silently mis-mapped. The same check runs on the query side
against PYR1's sequence.

Reads results/pocket_shape/pocket_shape.json, so run 143 first.
"""
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                # noqa: E402
m143 = import_module("143_pocket_backbone_shape")                  # noqa: E402
read_pdb, read_cif_resnames, THREE = m143.read_pdb, m143.read_cif_resnames, m143.THREE

HC = os.path.join(ROOT, "results", "homolog_cavities")
OUT = os.path.join(ROOT, "results", "pocket_shape")
MIN_MATCH = 0.90          # aligned-sequence agreement required to trust indices


def main():
    shp = {r["name"]: r for r in json.load(open(os.path.join(OUT, "pocket_shape.json")))}
    pyr1 = shp["PYR1"]
    p_lining = set(pyr1["lining"])
    p_ss = {n: pyr1["ss_string"][i] for i, n in enumerate(pyr1["ca_resnums"])}
    p_order = pyr1["ca_resnums"]                       # query index i -> resnum
    frame = os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")
    p_seq = {}
    for a in read_pdb(frame, want_atom=False):
        if a[1] in THREE:
            p_seq[a[0]] = THREE[a[1]]

    hits = {}
    for line in open(os.path.join(ROOT, "results", "foldseek", "fs_hits.tsv")):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 11:
            hits[f[1]] = f

    print(f"PYR1 pocket wall: {len(p_lining)} residues, "
          f"{sum(p_ss.get(r) == 'E' for r in p_lining)} on strands, "
          f"{sum(p_ss.get(r) == 'H' for r in p_lining)} on helices")

    out = []
    for name, r in shp.items():
        if name == "PYR1":
            continue
        tgt = name[:4].lower() + name[5] + "0" + ("1" if name[:4].lower() in
                                                  ("2ns9",) else "0")
        tgt = next((t for t in hits if t.lower().startswith(name[:4].lower())
                    and t[4] == name[5]), None)
        if tgt is None:
            continue
        f = hits[tgt]
        qs, ts, qaln, taln = int(f[5]), int(f[7]), f[9], f[10]

        dom = os.path.join(HC, "domains", tgt + ".pdb")
        names = read_cif_resnames(os.path.join(HC, "raw", tgt + ".cif"), tgt[4])
        t_order = [a[0] for a in read_pdb(dom) if a[2] == "CA"]
        t_seq = {n: names.get(n, "X") for n in t_order}

        # ---- verify the index convention on BOTH sides -----------------
        qi, ti, qmap, tins = qs - 1, ts - 1, {}, []
        okq = okt = tq = tt = 0
        for qc, tc in zip(qaln, taln):
            if qc != "-":
                qi += 1
                if qi <= len(p_order):
                    tq += 1
                    okq += (p_seq.get(p_order[qi - 1]) == qc)
            if tc != "-":
                ti += 1
                if ti <= len(t_order):
                    tt += 1
                    okt += (t_seq.get(t_order[ti - 1]) == tc)
        if not tq or not tt or okq / tq < MIN_MATCH or okt / tt < MIN_MATCH:
            print(f"  {name}: index check FAILED "
                  f"(query {okq/max(tq,1):.0%}, target {okt/max(tt,1):.0%}) -- dropped")
            continue

        # ---- walk again, recording the mapping and the insertions -------
        qi, ti = qs - 1, ts - 1
        aligned_t, run = {}, []
        inserts, last_q = [], None
        for qc, tc in zip(qaln, taln):
            if qc != "-":
                qi += 1
            if tc != "-":
                ti += 1
            if qc == "-" and tc != "-":
                run.append(ti)
            else:
                if run:
                    inserts.append((last_q, run))
                    run = []
                if qc != "-" and tc != "-":
                    aligned_t[ti] = qi
                if qc != "-":
                    last_q = qi
        if run:
            inserts.append((last_q, run))

        t_res = {i: t_order[i - 1] for i in range(1, len(t_order) + 1)}
        lin = set(r["lining"])
        lin_idx = {i for i, rn in t_res.items() if rn in lin}
        n_al = sum(1 for i in lin_idx if i in aligned_t)
        ins_idx = {i for _, rr in inserts for i in rr}
        n_in = sum(1 for i in lin_idx if i in ins_idx)
        n_un = len(lin_idx) - n_al - n_in

        # which PYR1 wall positions are conserved as wall in the target
        cons = sum(1 for i in lin_idx if aligned_t.get(i) and
                   p_order[aligned_t[i] - 1] in p_lining)

        rec = dict(name=name, cavity=r["cavity"], n_lining=len(lin_idx),
                   aligned=n_al, insert=n_in, unaligned=n_un, conserved_wall=cons,
                   inserts=[])
        for lq, rr in inserts:
            if lq is None or lq > len(p_order):
                continue
            nl = sum(1 for i in rr if i in lin_idx)
            if nl == 0:
                continue
            rec["inserts"].append(dict(after_pyr1=p_order[lq - 1], length=len(rr),
                                       lining=nl, pyr1_ss=p_ss.get(p_order[lq - 1], "?"),
                                       seq="".join(t_seq.get(t_res[i], "X") for i in rr)))
        out.append(rec)

    out.sort(key=lambda x: -x["cavity"])
    print("\n" + "=" * 92)
    print("HOW MUCH OF THE WALL IS REACHABLE BY MUTATION ALONE?")
    print("=" * 92)
    print(f"{'structure':<10}{'cavity':>8}{'wall':>6}{'aligned':>9}{'insert':>8}"
          f"{'unaln':>7}{'sharedWall':>12}")
    for x in out:
        print(f"{x['name']:<10}{x['cavity']:>8.0f}{x['n_lining']:>6}"
              f"{x['aligned']:>9}{x['insert']:>8}{x['unaligned']:>7}"
              f"{x['conserved_wall']:>12}")

    print("\n" + "=" * 92)
    print("POCKET-LINING INSERTIONS -- new backbone that would have to be built")
    print("=" * 92)
    for x in out:
        if not x["inserts"]:
            continue
        print(f"\n{x['name']}  (cavity {x['cavity']:.0f} A^3)")
        for d in x["inserts"]:
            print(f"   after PYR1 {d['after_pyr1']:>3} ({d['pyr1_ss']})  "
                  f"+{d['length']:>2} aa, {d['lining']} line the pocket   {d['seq']}")
    json.dump(out, open(os.path.join(OUT, "pocket_inserts.json"), "w"), indent=1)
    print(f"\nwritten to {OUT}/pocket_inserts.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
