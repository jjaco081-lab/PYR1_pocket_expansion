#!/usr/bin/env python
r"""
152_truncation_scan.py -- replacement for 148's withdrawn ABA-distance filter.

148 ranked the 52 borrowed positions by distance to the crystallographic ABA and
kept 8. That criterion is invalid for this project: it scores a candidate by how
close it is to the volume the pocket ALREADY has, so it can only ever select
positions that do not expand anything. Withdrawn.

The replacement asks the question the criterion was standing in for, directly:
if this side chain were removed, would the MAIN CHAMBER get bigger? Chamber
volume comes from 151 -- the connected sub-volume that a 2.4 A sphere can reach,
so satellite voids joined through sub-ligand channels are not counted as gain.

  dV_main   change in main-chamber volume on truncation to Ala, and to Gly
  merged    a satellite that was isolated becomes part of the main chamber
  r_max     change in the largest inscribed sphere -- girth, not just volume

CONTROLS, so the numbers mean something
  * PYR1's own 24 wall positions are scanned alongside the 52. If the borrowed
    positions do not open more volume than the wall positions already known to
    line the pocket, they are adding nothing.
  * F108 is scanned to I, L, V, A and G, because 3OQU carries I there and the
    lobe reading in README 87 says F108 is the mouth. This is the cheap test of
    that specific claim.

⚠ A known artefact of the enclosure criterion (lib_cavity docstring): truncating
a residue that SEALS the cavity lets the component escape and the reported
volume FALLS. A large negative dV is therefore a seal flag, not a shrinking
pocket, and is printed as SEAL rather than averaged into anything.
"""
import json
import os
import sys

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m143 = import_module("143_pocket_backbone_shape")                   # noqa: E402
m151 = import_module("151_cavity_bottleneck")                       # noqa: E402
read_pdb, THREE = m143.read_pdb, m143.THREE
clear_map, enclosed, SPACING = m151.clear_map, m151.enclosed, m151.SPACING

OUT = os.path.join(ROOT, "results", "pocket_shape")
BACKBONE = {"N", "CA", "C", "O", "OXT", "CB"}     # Ala keeps CB; Gly drops it
CHAMBER_PROBE = 2.4


def chamber(xyz, elem):
    """(main-chamber volume, total enclosed volume, r_max) -- 151's definition."""
    clear, _, lo, shape = clear_map(xyz, elem)
    mask = enclosed(clear, lo, shape)
    tot = float(mask.sum() * SPACING ** 3)
    if tot == 0:
        return 0.0, 0.0, 0.0
    core = mask & (clear > CHAMBER_PROBE)
    lab, n = ndimage.label(core)
    if n == 0:
        return 0.0, tot, float(clear[mask].max())
    sizes = ndimage.sum(core, lab, range(1, n + 1))
    big = int(np.argmax(sizes)) + 1
    seed = np.argwhere(lab == big)
    seed = seed[np.argmax(clear[tuple(seed.T)])]
    # main chamber = all cavity voxels reachable from the seed at probe 2.4,
    # then dilated back through the enclosed mask to recover its lining shell
    sub = mask & (clear > CHAMBER_PROBE)
    l2, _ = ndimage.label(sub)
    comp = l2 == l2[tuple(seed)]
    grown = ndimage.binary_dilation(comp, iterations=3) & mask
    return (float(grown.sum() * SPACING ** 3), tot, float(clear[mask].max()))


def main():
    at = read_pdb(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb"), want_atom=False)
    prot = [a for a in at if a[1] in THREE]
    seq = {a[0]: THREE[a[1]] for a in prot}
    base = chamber([a[4] for a in prot], [a[3] for a in prot])
    print(f"WT PYR1: main chamber {base[0]:.0f} A^3, total enclosed {base[1]:.0f}, "
          f"r_max {base[2]:.2f} A", flush=True)

    shp = {r["name"]: r for r in json.load(open(os.path.join(OUT, "pocket_shape.json")))}
    wall = sorted(shp["PYR1"]["lining"])
    bor = sorted(int(k) for k in
                 json.load(open(os.path.join(OUT, "borrowed_positions.json")))["votes"])
    groups = [("wall", wall), ("borrowed", bor)]

    rows = []
    for gname, positions in groups:
        for p in positions:
            if p not in seq or seq[p] == "G":
                continue
            for mode in ("ala", "gly"):
                keep = BACKBONE - ({"CB"} if mode == "gly" else set())
                sub = [a for a in prot if not (a[0] == p and a[2] not in keep)]
                m, t, rm = chamber([a[4] for a in sub], [a[3] for a in sub])
                rows.append(dict(group=gname, pos=p, aa=seq[p], mode=mode,
                                 main=round(m, 1), total=round(t, 1),
                                 dmain=round(m - base[0], 1),
                                 dtot=round(t - base[1], 1),
                                 r_max=round(rm, 2)))
            print(f"  {gname} {p}{seq[p]}: dV_main "
                  f"{rows[-2]['dmain']:+.0f} (Ala) {rows[-1]['dmain']:+.0f} (Gly)",
                  flush=True)

    # F108 identity scan -- the 3OQU reading
    print("\nF108 identity scan (3OQU carries I here)", flush=True)
    SC = {"I": {"CG1", "CG2", "CD1"}, "L": {"CG", "CD1", "CD2"},
          "V": {"CG1", "CG2"}, "A": set(), "G": set()}
    f108 = [a for a in prot if a[0] == 108]
    assert f108 and seq[108] == "F", "residue 108 is not F in this frame"
    for aa, extra in SC.items():
        keep = (BACKBONE - ({"CB"} if aa == "G" else set())) | extra
        sub = [a for a in prot if not (a[0] == 108 and a[2] not in keep)]
        m, t, rm = chamber([a[4] for a in sub], [a[3] for a in sub])
        rows.append(dict(group="F108", pos=108, aa=aa, mode="scan",
                         main=round(m, 1), total=round(t, 1),
                         dmain=round(m - base[0], 1), dtot=round(t - base[1], 1),
                         r_max=round(rm, 2)))
        print(f"  F108{aa}: main {m:>6.0f} ({m-base[0]:+6.0f})  "
              f"total {t:>6.0f}  r_max {rm:.2f}", flush=True)

    json.dump({"wt": dict(main=base[0], total=base[1], r_max=base[2]), "rows": rows},
              open(os.path.join(OUT, "truncation_scan.json"), "w"), indent=1)

    print("\n" + "=" * 82)
    print("TRUNCATION SCAN -- change in MAIN-CHAMBER volume, Gly arm, ranked")
    print("=" * 82)
    g = [r for r in rows if r["mode"] == "gly"]
    g.sort(key=lambda r: -r["dmain"])
    print(f"{'pos':>5}{'aa':>4}{'group':>10}{'dV_main':>10}{'dV_total':>10}{'r_max':>8}")
    for r in g[:20]:
        note = "  SEAL" if r["dtot"] < -20 else ""
        print(f"{r['pos']:>5}{r['aa']:>4}{r['group']:>10}{r['dmain']:>+10.0f}"
              f"{r['dtot']:>+10.0f}{r['r_max']:>8.2f}{note}")
    for gn in ("wall", "borrowed"):
        v = [r["dmain"] for r in g if r["group"] == gn]
        if v:
            print(f"\n   {gn:<9} n={len(v):>3}  median dV_main {np.median(v):+.0f}  "
                  f"best {max(v):+.0f}  n>+20 A^3: {sum(1 for x in v if x > 20)}")
    print(f"\nwritten to {OUT}/truncation_scan.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
