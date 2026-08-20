#!/usr/bin/env python
"""
79_graft_geometry.py -- can PYR1's transduction machinery be grafted INTO a
large-cavity relative, rather than the cavity being grafted into PYR1?

THE INVERSION, AND WHY IT IS THE BETTER DIRECTION
--------------------------------------------------
Everything before this tried to enlarge PYR1's pocket. That is the part we have
measured to be hard: expansion does not predict switch cost (r = 0.03, README 14),
F108 gatekeeps the second lobe, and the pocket is length-limited rather than
neck-limited (README 4e). Meanwhile CoxG (2PCS) already has a 570 A^3 pocket, 3.3x
PYR1's, for free.

So: take the big pocket, and move PYR1's machinery onto it.

WHAT "THE MACHINERY" ACTUALLY IS
--------------------------------
Not a loop. The HAB1 interface is 19 residues in FOUR discontinuous segments
(README 9c), and 15 of them are also the homodimer interface, which is why the apo
dimer is a competitive off-state. Any graft has to carry all four.

WHAT THIS SCRIPT MEASURES -- geometry only, no design
-----------------------------------------------------
For each candidate, after superposing on the shared CORE (aligned positions that
are NOT part of the machinery, so the segments being judged do not drag the fit):

  coverage  how many of the 19 interface residues have a structurally equivalent
            position at all, versus falling in an alignment gap
  indels    insertions/deletions inside each segment's span - a 4-residue insert
            where the latch belongs makes it a rebuild, not a transplant
  local     per-segment CA RMSD. The global alnTM averages over the whole domain
            and will happily hide a bad gate.
  cavity    the candidate's own pocket volume, since that is the entire point

A candidate passes only if all four segments are present, indel-free and locally
close. This is a screen for what is GEOMETRICALLY POSSIBLE; it says nothing about
whether the graft would still cycle between open and closed, which is the risk our
own MD work (README 24, 29) says is the real one.

Run with the pyr1_docking env python.
"""
import csv
import json
import os
from collections import defaultdict

import numpy as np
from Bio.PDB import PDBParser

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HC = os.path.join(ROOT, "results", "homolog_cavities")
OUT = os.path.join(ROOT, "results", "graft_candidates")
os.makedirs(OUT, exist_ok=True)

#: The HAB1 interface, README 9c, grouped into the four contiguous stretches it
#: actually forms. 15 of these 19 are also the homodimer interface.
SEGMENTS = {
    "helix-60s":  [60, 61, 63],
    "gate 84-89": [84, 85, 86, 87, 88, 89],
    "latch":      [116, 117],
    "C-lobe":     [148, 151, 155, 156, 158, 159, 162, 166],
}
INTERFACE = sorted(p for seg in SEGMENTS.values() for p in seg)
MIN_CAVITY = 100.0        # A^3 - below PYR1's own 174 there is no point grafting
p = PDBParser(QUIET=True)


def ca_by_order(path, chain=None):
    """[(resnum, CA coord)] in file order -- the order foldseek indexes into."""
    st = p.get_structure("x", path)[0]
    out = []
    for ch in st:
        if chain and ch.id != chain:
            continue
        for r in ch:
            if r.id[0] == " " and "CA" in r:
                out.append((r.id[1], r["CA"].coord))
    return out


def superpose(P, Q):
    """Kabsch: rotation taking Q onto P; returns (R, muP, muQ, rmsd)."""
    mp, mq = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((Q - mq).T @ (P - mp))
    R = U @ np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]) @ Vt
    rms = float(np.sqrt((((Q - mq) @ R - (P - mp)) ** 2).sum() / len(P)))
    return R, mp, mq, rms


def main():
    pyr1 = ca_by_order(os.path.join(ROOT, "data", "pyr1_A.pdb"))
    nums = [n for n, _ in pyr1]
    seqidx = {n: i + 1 for i, n in enumerate(nums)}      # PYR1 resnum -> query index
    pyr1_xyz = {n: c for n, c in pyr1}
    missing = [x for x in INTERFACE if x not in seqidx]
    assert not missing, f"interface residues absent from the query structure: {missing}"

    cav = {r["target"]: float(r["cavity_A3"])
           for r in csv.DictReader(open(os.path.join(HC, "homolog_cavities.csv")))}

    rows = []
    for line in open(os.path.join(ROOT, "results", "foldseek", "fs_hits.tsv")):
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        tgt, fid, tm = f[1], float(f[2]), float(f[3])
        if tm < 0.5 or tgt.startswith("af_"):
            continue                      # models carry no ligand and no crystal pocket
        dom = os.path.join(HC, "domains", tgt + ".pdb")
        if not os.path.exists(dom):
            continue
        volume = cav.get(tgt)
        if volume is None or volume < MIN_CAVITY:
            continue

        qs, ts, qaln, taln = int(f[5]), int(f[7]), f[9], f[10]
        # walk the alignment: PYR1 query index -> target index, and record gaps
        qi, ti, qmap = qs - 1, ts - 1, {}
        for qc, tc in zip(qaln, taln):
            if qc != "-":
                qi += 1
            if tc != "-":
                ti += 1
            if qc != "-":
                qmap[qi] = ti if tc != "-" else None      # None = deletion in target

        tgt_ca = ca_by_order(dom)
        tgt_xyz = {i + 1: c for i, (_, c) in enumerate(tgt_ca)}

        # superpose on the shared core: aligned pairs OUTSIDE the machinery
        iface_qi = {seqidx[x] for x in INTERFACE}
        pairs = [(qi_, ti_) for qi_, ti_ in qmap.items()
                 if ti_ is not None and qi_ not in iface_qi and ti_ in tgt_xyz]
        if len(pairs) < 40:
            continue
        idx2num = {i: n for n, i in seqidx.items()}
        P = np.array([pyr1_xyz[idx2num[a]] for a, _ in pairs])
        Q = np.array([tgt_xyz[b] for _, b in pairs])
        R, mp, mq, core_rms = superpose(P, Q)

        seg_stats, covered = {}, 0
        for name, poss in SEGMENTS.items():
            got, dev = 0, []
            for pos in poss:
                ti_ = qmap.get(seqidx[pos])
                if ti_ is None or ti_ not in tgt_xyz:
                    continue
                got += 1
                moved = (tgt_xyz[ti_] - mq) @ R + mp
                dev.append(float(np.linalg.norm(moved - pyr1_xyz[pos])))
            covered += got
            # ⚠ Net insertion/deletion relative to PYR1, which means comparing the
            # TARGET span against the QUERY span - not against the number of selected
            # positions. The interface residues are not contiguous in PYR1 either
            # (the C-lobe picks 8 positions out of 148-166), so measuring against
            # len(poss) reported an 11-residue "insertion" for PYL homologs at
            # alnTM 0.97, which is nonsense.
            span = [qmap.get(seqidx[x]) for x in poss]
            real = [v for v in span if v is not None]
            if len(real) > 1:
                qspan = seqidx[max(poss)] - seqidx[min(poss)] + 1
                tspan = max(real) - min(real) + 1
                ins = tspan - qspan
            else:
                ins = 0
            seg_stats[name] = dict(present=got, of=len(poss), insert=ins,
                                   rmsd=float(np.mean(dev)) if dev else None)
        rows.append(dict(target=tgt, pdb=tgt[:4].upper(), alntm=tm, fident=fid,
                         cavity=volume, core_rmsd=core_rms, n_core=len(pairs),
                         covered=covered, of=len(INTERFACE), segments=seg_stats))

    rows.sort(key=lambda r: (-r["covered"], r["core_rmsd"]))
    print("=" * 100)
    print("GRAFT GEOMETRY SCREEN -- can PYR1's 19-residue machinery be placed on a bigger pocket?")
    print("=" * 100)
    print(f"  candidates: solved, alnTM >= 0.5, cavity >= {MIN_CAVITY:.0f} A^3 (PYR1 is 174)")
    print(f"  {'PDB':<6}{'alnTM':>6}{'cav':>7}{'core':>6}{'iface':>7}   "
          + "".join(f"{n:<16}" for n in SEGMENTS))
    print(f"  {'':<6}{'':>6}{'A^3':>7}{'RMSD':>6}{'/19':>7}   "
          + "".join(f"{'n  ins  RMSD':<16}" for _ in SEGMENTS))
    for r in rows:
        cells = ""
        for n in SEGMENTS:
            s = r["segments"][n]
            rm = f"{s['rmsd']:.1f}" if s["rmsd"] is not None else "  -"
            cells += f"{s['present']}/{s['of']} {s['insert']:>2}  {rm:>4}    "
        print(f"  {r['pdb']:<6}{r['alntm']:>6.2f}{r['cavity']:>7.0f}"
              f"{r['core_rmsd']:>6.2f}{r['covered']:>4}/19   {cells}")
    json.dump(rows, open(os.path.join(OUT, "graft_geometry.json"), "w"), indent=1)
    print(f"\n  wrote {OUT}/graft_geometry.json")
    return rows


if __name__ == "__main__":
    main()
