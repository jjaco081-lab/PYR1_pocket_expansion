#!/usr/bin/env python
r"""
145_wall_displacement.py -- insertion or displacement? and two fixes to 143.

Script 144 found that the wall of every large-cavity relative is built almost
entirely from ALIGNED positions: 2PCS reaches 570 A^3 with 43 of its 45 wall
residues sitting at positions PYR1 also has, and 2 inserted. So the extra volume
is not made of extra backbone. This script tests the alternative -- that the
SAME backbone positions are further apart -- and fixes two defects in 143.

DEFECT 1, CIRCULARITY. 143 defined the envelope as the C-beta hull of the LINING
residues, and lining was defined from the cavity. A bigger cavity therefore
recruits more residues into its own hull, so the near-constant fill fraction
(0.88-0.95 across 20 structures) may be definitional rather than biological.
Fixed here by defining the wall STRUCTURALLY: the residues aligned to PYR1's own
24 lining positions, the same set in every structure, mapped by foldseek. The
hull is then cavity-independent and the comparison is honest.

DEFECT 2, THE poly-Gly COLUMN WAS MEANINGLESS. 143 reported 0-3 A^3 for PYR1,
2BK0, 3QRZ, 1Z94 and 3P51. That is not a small cavity: it is the enclosed
component ceasing to exist once the side chains come off, because buriedness
falls below the 0.88 cut, leaving the seeded search to pick some unrelated void
elsewhere. The leak flag did not catch it because a component that VANISHES
never touches the grid boundary. Fixed by abandoning enclosure entirely: free
volume is counted inside the wall hull, which bounds the region by construction
and so cannot leak.

WHAT IS MEASURED
  envelope    convex hull of the C-beta atoms of the structurally aligned wall
  free        solvent-probe-free volume inside that hull, side chains present
  free_gly    the same with every wall side chain removed past CA -- the ceiling
              mutation alone can reach, now leak-proof
  sheet_helix mean CA-CA distance between the wall residues PYR1 carries on its
              beta sheet and those it carries on the grip helix, measured at the
              structurally equivalent positions. In a helix-grip fold this is the
              single degree of freedom that sets pocket size; if the big pockets
              differ here, displacement is the mechanism and insertion is not.

Reads 143's JSON for PYR1's wall definition and secondary structure.
"""
import json
import os
import sys

import numpy as np
from scipy.spatial import cKDTree, ConvexHull, Delaunay

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_cavity import VDW                                          # noqa: E402
from importlib import import_module                                 # noqa: E402
m143 = import_module("143_pocket_backbone_shape")                   # noqa: E402
read_pdb, read_cif_resnames, THREE = m143.read_pdb, m143.read_cif_resnames, m143.THREE
extents = m143.extents

HC = os.path.join(ROOT, "results", "homolog_cavities")
OUT = os.path.join(ROOT, "results", "pocket_shape")
BACKBONE = {"N", "CA", "C", "O", "OXT"}
SPACING = 0.5
PROBE = 1.4
MIN_MATCH = 0.90


def free_volume(hull_pts, xyz, elem, spacing=SPACING, probe=PROBE):
    """Probe-accessible volume strictly inside the convex hull of hull_pts.

    No buriedness test and no connected-component step, so nothing can leak:
    the hull is the boundary. Returns (free_A3, hull_A3).
    """
    H = ConvexHull(hull_pts)
    tri = Delaunay(hull_pts[H.vertices])
    lo, hi = hull_pts.min(0), hull_pts.max(0)
    axes = [np.arange(lo[i], hi[i] + spacing, spacing) for i in range(3)]
    pts = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
    inside = tri.find_simplex(pts) >= 0
    pts = pts[inside]
    if len(pts) == 0:
        return 0.0, float(H.volume)
    rad = np.array([VDW.get(e, 1.70) for e in elem])
    d, i = cKDTree(np.asarray(xyz, float)).query(pts, k=1)
    return float((d - rad[i] > probe).sum() * spacing ** 3), float(H.volume)


def main():
    shp = {r["name"]: r for r in json.load(open(os.path.join(OUT, "pocket_shape.json")))}
    pyr1 = shp["PYR1"]
    p_order = pyr1["ca_resnums"]
    p_ss = {n: pyr1["ss_string"][i] for i, n in enumerate(p_order)}
    WALL = list(pyr1["lining"])
    qidx = {n: i + 1 for i, n in enumerate(p_order)}
    wall_q = {qidx[r] for r in WALL if r in qidx}
    sheet_q = {qidx[r] for r in WALL if p_ss.get(r) == "E" and r in qidx}
    helix_q = {qidx[r] for r in WALL if p_ss.get(r) == "H" and r in qidx}
    print(f"structural wall = PYR1's {len(WALL)} lining positions "
          f"({len(sheet_q)} on sheet, {len(helix_q)} on the grip helix); "
          f"the SAME set is used in every structure")

    frame = os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")
    p_at = [a for a in read_pdb(frame, want_atom=False) if a[1] in THREE]

    hits = {}
    for line in open(os.path.join(ROOT, "results", "foldseek", "fs_hits.tsv")):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 11:
            hits[f[1]] = f

    def measure(name, atoms, cav, wall_res, sheet_res, helix_res):
        by = {}
        for a in atoms:
            by.setdefault(a[0], {})[a[2]] = a[4]
        cb = [by[r].get("CB", by[r].get("CA")) for r in wall_res
              if r in by and ("CB" in by[r] or "CA" in by[r])]
        cb = np.array([c for c in cb if c is not None])
        if len(cb) < 4:
            return None
        strip = set(wall_res)
        keep = [(a[4], a[3]) for a in atoms
                if not (a[0] in strip and a[2] not in BACKBONE)]
        fw, hv = free_volume(cb, [a[4] for a in atoms], [a[3] for a in atoms])
        fg, _ = free_volume(cb, [c for c, _ in keep], [e for _, e in keep])
        S = np.array([by[r]["CA"] for r in sheet_res if r in by and "CA" in by[r]])
        Hh = np.array([by[r]["CA"] for r in helix_res if r in by and "CA" in by[r]])
        sh = (float(np.linalg.norm(S[:, None] - Hh[None], axis=-1).mean())
              if len(S) and len(Hh) else float("nan"))
        E = extents(cb)
        return dict(name=name, cavity=cav, n_wall=len(cb), envelope=round(hv, 1),
                    free=round(fw, 1), free_gly=round(fg, 1),
                    fill=round(1 - fw / hv, 3), gly_gain=round(fg - fw, 1),
                    sheet_helix=round(sh, 2),
                    E1=round(E[0], 1), E2=round(E[1], 1), E3=round(E[2], 1))

    rows = [measure("PYR1", p_at, pyr1["cavity"], WALL,
                    [r for r in WALL if p_ss.get(r) == "E"],
                    [r for r in WALL if p_ss.get(r) == "H"])]

    for name, r in shp.items():
        if name == "PYR1":
            continue
        tgt = next((t for t in hits if t.lower().startswith(name[:4].lower())
                    and t[4] == name[5]), None)
        if tgt is None:
            continue
        f = hits[tgt]
        qs, ts, qaln, taln = int(f[5]), int(f[7]), f[9], f[10]
        dom = os.path.join(HC, "domains", tgt + ".pdb")
        names = read_cif_resnames(os.path.join(HC, "raw", tgt + ".cif"), tgt[4])
        atoms = [(a[0], names.get(a[0], "UNK"), a[2], a[3], a[4])
                 for a in read_pdb(dom)]
        t_order = [a[0] for a in atoms if a[2] == "CA"]

        qi, ti, qmap = qs - 1, ts - 1, {}
        okt = tt = 0
        for qc, tc in zip(qaln, taln):
            if qc != "-":
                qi += 1
            if tc != "-":
                ti += 1
                if ti <= len(t_order):
                    tt += 1
                    okt += (names.get(t_order[ti - 1], "X") == tc)
            if qc != "-" and tc != "-":
                qmap[qi] = ti
        if not tt or okt / tt < MIN_MATCH:
            print(f"  {name}: target index check FAILED ({okt/max(tt,1):.0%}) -- dropped")
            continue
        pick = lambda S: [t_order[qmap[q] - 1] for q in sorted(S)
                          if q in qmap and qmap[q] <= len(t_order)]
        w = pick(wall_q)
        if len(w) < 0.7 * len(wall_q):
            print(f"  {name}: only {len(w)}/{len(wall_q)} wall positions aligned -- dropped")
            continue
        rr = measure(name, atoms, r["cavity"], w, pick(sheet_q), pick(helix_q))
        if rr:
            rows.append(rr)

    rows.sort(key=lambda x: -x["cavity"])
    print("\n" + "=" * 94)
    print("STRUCTURALLY EQUIVALENT WALL -- same 24 positions in every structure")
    print("=" * 94)
    print(f"{'structure':<10}{'cavity':>8}{'wall':>6}{'envelope':>10}{'free':>8}"
          f"{'fill':>7}{'freeGly':>9}{'gain':>7}{'sheet-helix':>13}")
    for x in rows:
        print(f"{x['name']:<10}{x['cavity']:>8.0f}{x['n_wall']:>6}"
              f"{x['envelope']:>10.0f}{x['free']:>8.0f}{x['fill']:>7.2f}"
              f"{x['free_gly']:>9.0f}{x['gly_gain']:>7.0f}{x['sheet_helix']:>13.2f}")
    v = [x for x in rows if x["name"] != "PYR1"]
    p = rows[[r["name"] for r in rows].index("PYR1")]
    print(f"\nPYR1 reference: envelope {p['envelope']:.0f} A^3, free {p['free']:.0f}, "
          f"sheet-helix {p['sheet_helix']:.2f} A")
    if len(v) > 2:
        c = np.array([x["cavity"] for x in rows])
        for k in ("envelope", "sheet_helix", "free"):
            y = np.array([x[k] for x in rows], float)
            m = ~np.isnan(y)
            r_ = float(np.corrcoef(c[m], y[m])[0, 1])
            print(f"   cavity vs {k:<12} r = {r_:+.3f}  (n = {m.sum()})")
    json.dump(rows, open(os.path.join(OUT, "wall_displacement.json"), "w"), indent=1)
    print(f"\nwritten to {OUT}/wall_displacement.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
