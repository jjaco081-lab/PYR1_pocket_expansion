#!/usr/bin/env python
r"""
192_filter_cascade.py -- the structural filter cascade, with WILD-TYPE PYR1 AS A
HARD GATE.

Jannis's condition: PYR1 must survive every filter before any filter is applied
to a design. A filter that rejects wild-type PYR1 is miscalibrated, not strict.
This is the control that catches the class of error that has hit this project
four times in one session -- a circular envelope (§143), a poly-Gly column that
read a vanished component as 0-3 A^3 (§143), a missing connected-component step
(§145, caught by Jannis), and an ABA-distance criterion that defined the pocket
as the volume already occupied (§149, caught by Jannis).

FILTERS, cheapest first. Each reports a count in and out so the funnel is
measured rather than asserted.

  F1 has an enclosed cavity at all
  F2 largest-CONNECTED poly-Gly envelope >= PYR1's
  F3 motif integrity -- CA deviation and chain breaks, from RFd3's own metrics
  F4 no HAB1-wrapping leading segment
  F5 no pocket-spanning residue
  F6 compactness and fold plausibility -- Rg, sheet fraction, n SS elements
  F7 shape: depth and depth/width, REPORTED not gated (§106: the importable
     difference across natural donors is depth, and PYR1 is only d/w 1.25)

⚠ Every volume here is the LARGEST CONNECTED COMPONENT inside a structurally
defined C-beta hull. §145's free_volume summed every disconnected void, which
credits a loose scaffold with the total of several small pockets. Measured
inflation on batches 1-3 was 1.00x so it changed nothing there, but it is the
correct definition and it matters the moment a pocket-spanning residue splits a
cavity -- which is exactly filter F5's failure mode.

Usage:
  python 192_filter_cascade.py --calibrate          # PYR1 only, the hard gate
  python 192_filter_cascade.py --designs <glob>     # apply to designs
"""
import argparse, glob, json, os, sys
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree, ConvexHull, Delaunay

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_cavity import VDW                                          # noqa: E402
from importlib import import_module                                 # noqa: E402
m143 = import_module("143_pocket_backbone_shape")
m151 = import_module("151_cavity_bottleneck")
BB = {"N", "CA", "C", "O", "OXT"}
SPACING, PROBE = 0.5, 1.4

#: PYR1 reference, all measured in §106/§113/§190. The gate is that PYR1 passes.
PYR1_ENVELOPE = 1184.0
PYR1_RMAX = 3.21
PYR1_RG = 15.0
PYR1_SHEET = 0.27          # sheet fraction of the helix-grip fold
THRESH = dict(envelope=PYR1_ENVELOPE, rmax=PYR1_RMAX,
              ca_dev=1.0, breaks=0, hab1_lead=5, rg_mult=1.5,
              sheet_frac=0.15, n_ss=8, span_dist=4.0)


def free_volume_connected(hull_pts, xyz, elem):
    """(largest_connected_free_A3, n_components, hull_A3) -- see module docstring."""
    H = ConvexHull(hull_pts)
    tri = Delaunay(hull_pts[H.vertices])
    lo, hi = hull_pts.min(0), hull_pts.max(0)
    axes = [np.arange(lo[i], hi[i] + SPACING, SPACING) for i in range(3)]
    shape = tuple(len(a) for a in axes)
    pts = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
    inside = tri.find_simplex(pts) >= 0
    rad = np.array([VDW.get(e, 1.70) for e in elem])
    d, i = cKDTree(np.asarray(xyz, float)).query(pts, k=1)
    free = inside & ((d - rad[i]) > PROBE)
    if not free.any():
        return 0.0, 0, float(H.volume)
    lab, n = ndimage.label(free.reshape(shape))
    sizes = ndimage.sum(free.reshape(shape), lab, range(1, n + 1))
    return float(sizes.max() * SPACING ** 3), int(n), float(H.volume)


def envelope(atoms, wall):
    by = {}
    for a in atoms:
        by.setdefault(a[0], {})[a[2]] = a[4]
    cb = np.array([by[r].get("CB", by[r]["CA"]) for r in wall
                   if r in by and "CA" in by[r]])
    if len(cb) < 4:
        return None
    strip = set(wall)
    keep = [(a[4], a[3]) for a in atoms if not (a[0] in strip and a[2] not in BB)]
    big, n, hull = free_volume_connected(cb, [c for c, _ in keep],
                                         [e for _, e in keep])
    E = m143.extents(cb)
    return dict(envelope=big, ncomp=n, hull=hull,
                depth=E[0], w1=E[1], w2=E[2],
                dw=E[0] / max(E[2], 1e-9))


def cavity_and_wall(atoms):
    r = m151.analyse("x", [a[4] for a in atoms], [a[3] for a in atoms])
    if not r or r.get("main", 0) <= 0:
        return None, None, None
    clear, _, lo, shape = m151.clear_map([a[4] for a in atoms], [a[3] for a in atoms])
    mask = m151.enclosed(clear, lo, shape)
    if not mask.sum():
        return None, None, None
    pts = np.argwhere(mask) * m151.SPACING + lo
    by = {}
    for a in atoms:
        by.setdefault(a[0], {})[a[2]] = a[4]
    wall = [k for k, d in by.items() if "CA" in d
            and np.linalg.norm(d["CA"] - pts.mean(0)) < 11.0]
    return r, wall, pts


def spanning_residue(atoms, pts, wall):
    """A residue SPLITS the cavity if truncating its side chain to CB MERGES
    cavity components -- i.e. with the side chain present the pocket is divided.

    ⚠ The first version of this filter asked whether a residue's CB lies within
    4 A of the cavity's long axis, and it REJECTED WILD-TYPE PYR1, flagging five
    ordinary wall residues (61, 81, 92, 160, 163). That was the gate working:
    PYR1's chamber is 11.6 A long with r_max 3.21 A, so every lining residue is
    within 4 A of its own axis. The old filter measured "lines the pocket", not
    "spans it".

    The merge test measures the thing Jannis actually objected to -- a side chain
    that interrupts the pocket -- and it is defined against the connected-component
    count, so it is consistent with the envelope metric.
    """
    base = _ncomp(atoms)
    if base is None:
        return []
    out = []
    for r in wall:
        trimmed = [a for a in atoms if not (a[0] == r and a[2] not in BB)]
        n = _ncomp(trimmed)
        if n is not None and n < base:
            out.append(r)
    return out


def _ncomp(atoms):
    """number of enclosed cavity components"""
    try:
        clear, _, lo, shape = m151.clear_map([a[4] for a in atoms],
                                             [a[3] for a in atoms])
        mask = m151.enclosed(clear, lo, shape)
        if not mask.sum():
            return None
        return int(ndimage.label(mask)[1])
    except Exception:
        return None


def calibrate():
    at = [(a[0], a[1], a[2], a[3], a[4]) for a in
          m143.read_pdb(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb"),
                        want_atom=False) if a[1] in m143.THREE]
    shp = {r["name"]: r for r in json.load(open(os.path.join(
        ROOT, "results", "pocket_shape", "pocket_shape.json")))}
    wall = list(shp["PYR1"]["lining"])
    print("HARD GATE -- wild-type PYR1 through the cascade")
    r, w_auto, pts = cavity_and_wall(at)
    e = envelope(at, wall)
    ca = np.array([a[4] for a in at if a[2] == "CA"])
    rg = float(np.sqrt(((ca - ca.mean(0)) ** 2).sum(1).mean()))
    span = spanning_residue(at, pts, wall) if pts is not None else []
    checks = [
        ("F1 enclosed cavity exists", r is not None and r.get("main", 0) > 0, f"main {r.get('main',0):.0f} A^3"),
        ("F2 envelope >= PYR1", e["envelope"] >= THRESH["envelope"] * 0.999, f"{e['envelope']:.0f} vs {THRESH['envelope']:.0f}"),
        ("F2b r_max >= PYR1", r.get("r_max", 0) >= THRESH["rmax"] * 0.999, f"{r.get('r_max',0):.2f} vs {THRESH['rmax']}"),
        ("F5 no spanning residue", len(span) == 0, f"{len(span)} found {span if span else ''}"),
        ("F6 Rg within 1.5x", rg <= PYR1_RG * THRESH["rg_mult"], f"Rg {rg:.1f}"),
    ]
    ok = True
    for name, passed, detail in checks:
        print(f"   {'PASS' if passed else 'FAIL':<5}{name:<28} {detail}")
        ok &= passed
    print(f"\n   depth {e['depth']:.1f}  widths {e['w1']:.1f}/{e['w2']:.1f}  "
          f"d/w {e['dw']:.2f}   components {e['ncomp']}")
    print(f"\n   {'GATE PASSES -- filters are usable' if ok else 'GATE FAILS -- RECALIBRATE, do not apply to designs'}")
    json.dump(dict(envelope=e, rg=rg, main=r.get("main", 0),
                   r_max=r.get("r_max", 0), spanning=span, gate_ok=bool(ok)),
              open(os.path.join(ROOT, "results", "filter_cascade_pyr1.json"), "w"),
              indent=1)
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calibrate", action="store_true")
    a = ap.parse_args()
    if a.calibrate:
        return calibrate()
    print("run with --calibrate first")
    return 0


if __name__ == "__main__":
    sys.exit(main())
