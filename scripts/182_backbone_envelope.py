#!/usr/bin/env python
r"""
182_backbone_envelope.py -- what does the BACKBONE give, independent of side chains?

Jannis: cavity volume is the wrong comparator for a scaffold decision, because it
conflates two things. If a donor's larger cavity is just smaller side chains,
that change is available far more cheaply by mutating PYR1 itself. The question
for a graft or a diffused scaffold is whether the BACKBONE encloses more space.

⚠ THIS METRIC HAS BEEN WRONG TWICE. §143 defined the envelope as the C-beta hull
of the LINING residues, which the cavity itself selected -- circular, and it made
cavity vs envelope look like r = +0.99 when a structurally-defined wall gives
+0.20. §143's poly-Gly column was worse: 0-3 A^3 for five structures, which was
the enclosed component CEASING TO EXIST once buriedness fell below cut, not a
small cavity, and the leak flag missed it because a vanished component never
touches the grid boundary. §145 fixed both and this script uses that version:

  wall        the residues aligned to PYR1's own 24 lining positions -- the SAME
              set in every structure, mapped by foldseek. Never cavity-derived.
  envelope    convex hull of those residues' C-beta atoms
  free        probe-accessible volume INSIDE that hull, side chains present
  polyGly     the same with every wall side chain deleted past CA. Bounded by the
              hull, so it cannot leak -- which is what broke §143.
  depth/width principal-axis extents of the polyGly free volume: L1 (depth) and
              L2, L3 (width). This is the shape of the space the BACKBONE makes.

READ IT THIS WAY. A donor whose polyGly volume matches PYR1's has no backbone
advantage -- its bigger cavity is side chains, and mutation in PYR1 gets there
more cheaply. A donor with a genuinely larger polyGly volume is offering
something a graft or a scaffold could actually import.

⚠ Structures are written for inspection. Per Jannis's standing rule they should
be looked at before anything is built on them.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m143 = import_module("143_pocket_backbone_shape")
m145 = import_module("145_wall_displacement")

OUT = os.path.join(ROOT, "results", "backbone_envelope")
HC = os.path.join(ROOT, "results", "homolog_cavities")
SD = os.path.join(ROOT, "data", "start_domains")


BB = {"N", "CA", "C", "O", "OXT"}


def measure(name, atoms, wall_res):
    """§145's measurement, inlined -- it is nested inside 145.main and not importable."""
    by = {}
    for a in atoms:
        by.setdefault(a[0], {})[a[2]] = a[4]
    cb = [by[r].get("CB", by[r].get("CA")) for r in wall_res if r in by]
    cb = np.array([c for c in cb if c is not None])
    if len(cb) < 4:
        return None
    strip = set(wall_res)
    keep = [(a[4], a[3]) for a in atoms if not (a[0] in strip and a[2] not in BB)]
    fw, hv = m145.free_volume(cb, [a[4] for a in atoms], [a[3] for a in atoms])
    fg, _ = m145.free_volume(cb, [c for c, _ in keep], [e for _, e in keep])
    E = m143.extents(cb)
    return dict(name=name, n_wall=len(cb), envelope=round(hv, 1),
                free=round(fw, 1), free_gly=round(fg, 1),
                fill=round(1 - fw / hv, 3) if hv else 0.0,
                E1=round(E[0], 1), E2=round(E[1], 1), E3=round(E[2], 1))


def aln_map(hits_path, target, pyr1_order, want_q):
    """PYR1 query index -> target residue number, from a foldseek line."""
    for line in open(hits_path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 11 and f[1].replace(".pdb", "") == target:
            break
    else:
        return None
    qs, ts, qaln, taln = int(f[5]), int(f[7]), f[9], f[10]
    qi, ti, out = qs - 1, ts - 1, {}
    for qc, tc in zip(qaln, taln):
        if qc != "-":
            qi += 1
        if tc != "-":
            ti += 1
        if qc != "-" and tc != "-":
            out[qi] = ti
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    shp = {r["name"]: r for r in json.load(open(os.path.join(
        ROOT, "results", "pocket_shape", "pocket_shape.json")))}
    pyr1 = shp["PYR1"]
    p_order = pyr1["ca_resnums"]
    WALL = list(pyr1["lining"])
    qidx = {n: i + 1 for i, n in enumerate(p_order)}
    wall_q = sorted(qidx[r] for r in WALL if r in qidx)

    rows = []
    p_at = [a for a in m143.read_pdb(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb"),
                                     want_atom=False) if a[1] in m143.THREE]
    rows.append(("PYR1", measure("PYR1", p_at, WALL)))

    SRC = [("2PCS", HC + "/domains/2pcsA00.pdb", ROOT + "/results/foldseek/fs_hits.tsv", "2pcsA00"),
           ("2NS9", HC + "/domains/2ns9B01.pdb", ROOT + "/results/foldseek/fs_hits.tsv", "2ns9B01"),
           ("2BK0", HC + "/domains/2bk0A00.pdb", ROOT + "/results/foldseek/fs_hits.tsv", "2bk0A00"),
           ("6AWV", HC + "/domains/6awvC00.pdb", ROOT + "/results/foldseek/fs_hits.tsv", "6awvC00"),
           ("4DSB", HC + "/domains/4dsbB00.pdb", ROOT + "/results/foldseek/fs_hits.tsv", "4dsbB00"),
           ("3OQU", HC + "/domains/3oquB00.pdb", ROOT + "/results/foldseek/fs_hits.tsv", "3oquB00"),
           ("1LN1", SD + "/pdb/1LN1.pdb", SD + "/aln.tsv", "1LN1"),
           ("1EM2", SD + "/pdb/1EM2.pdb", SD + "/aln.tsv", "1EM2"),
           ("2E3M", SD + "/pdb/2E3M.pdb", SD + "/aln.tsv", "2E3M")]
    for name, pdb, hits, tgt in SRC:
        if not os.path.exists(pdb):
            print(f"  {name}: missing {pdb}"); continue
        at = m143.read_pdb(pdb)
        order = [a[0] for a in at if a[2] == "CA"]
        amap = aln_map(hits, tgt, p_order, wall_q)
        if amap is None:
            print(f"  {name}: no alignment"); continue
        w = [order[amap[q] - 1] for q in wall_q if q in amap and amap[q] <= len(order)]
        if len(w) < 0.7 * len(wall_q):
            print(f"  {name}: only {len(w)}/{len(wall_q)} wall positions aligned -- skipped")
            continue
        r = measure(name, at, w)
        rows.append((name, r))

    print(f"\nBACKBONE ENVELOPE -- the same {len(wall_q)} PYR1 wall positions, mapped into each")
    print(f"{'structure':<10}{'wall':>5}{'envelope':>10}{'free':>8}{'fill':>7}"
          f"{'polyGly':>9}{'  depth':>8}{'width':>7}{'width':>7}{'  d/w':>7}")
    for name, r in rows:
        if r is None:
            continue
        e1, e2, e3 = r["E1"], r["E2"], r["E3"]
        print(f"{name:<10}{r['n_wall']:>5}{r['envelope']:>10.0f}{r['free']:>8.0f}"
              f"{r['fill']:>7.2f}{r['free_gly']:>9.0f}{e1:>8.1f}{e2:>7.1f}{e3:>7.1f}"
              f"{e1/max(e3,0.01):>7.2f}")
    print("\n  polyGly is the volume the BACKBONE encloses once every wall side chain")
    print("  is cut to CA. A donor at PYR1's polyGly offers no backbone advantage --")
    print("  its extra cavity is side chains, and mutating PYR1 is cheaper.")
    json.dump([{**r, "name": n} for n, r in rows if r],
              open(os.path.join(OUT, "backbone_envelope.json"), "w"), indent=1)
    print(f"\n  written to {OUT}/backbone_envelope.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
