#!/usr/bin/env python
r"""
190_envelope_connected.py -- envelope volume, LARGEST CONNECTED COMPONENT ONLY.

Jannis: "when doing the pocket or envelope calculations, please make sure to only
count the largest pocket and not any disconnected ones."

He is right and §145's `free_volume` has this flaw by construction -- its own
docstring says "no buriedness test and no connected-component step". That was
written as leak-proofing (the hull is the boundary, so nothing can escape), but
it means the number SUMS EVERY DISCONNECTED VOID inside the hull. A loose,
extended scaffold with several small pockets is therefore credited with their
total, and compared against PYR1's single chamber that is not a fair number.

This recomputes every structure with a 6-connected component step on the free
grid and reports only the largest, so all numbers are single-pocket volumes.
Both are printed so the size of the correction is visible rather than hidden.

⚠ This changes PYR1's reference value too, so the whole comparison is redone;
no number is carried over from §106/§111/§113.
"""
import glob, json, os, sys
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


def free_volume_connected(hull_pts, xyz, elem):
    """(total_free, largest_connected_free, n_components, hull_volume)"""
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
    tot = float(free.sum() * SPACING ** 3)
    if not free.any():
        return 0.0, 0.0, 0, float(H.volume)
    lab, n = ndimage.label(free.reshape(shape))
    sizes = ndimage.sum(free.reshape(shape), lab, range(1, n + 1))
    return tot, float(sizes.max() * SPACING ** 3), int(n), float(H.volume)


def measure(atoms, wall):
    by = {}
    for a in atoms:
        by.setdefault(a[0], {})[a[2]] = a[4]
    cb = np.array([by[r].get("CB", by[r]["CA"]) for r in wall
                   if r in by and "CA" in by[r]])
    if len(cb) < 4:
        return None
    strip = set(wall)
    keep = [(a[4], a[3]) for a in atoms if not (a[0] in strip and a[2] not in BB)]
    return free_volume_connected(cb, [c for c, _ in keep], [e for _, e in keep])


def read_cif(f):
    cols, rows, inl = {}, [], False
    for line in open(f):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".")[1]] = len(cols); inl = True; continue
        if inl:
            if line.startswith(("#", "loop_", "_")):
                if rows: break
                continue
            x = line.split()
            if len(x) >= len(cols): rows.append(x)
    return [(r[cols["label_asym_id"]], int(r[cols["label_seq_id"]]),
             r[cols["label_atom_id"]], r[cols["type_symbol"]],
             np.array([float(r[cols[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")]))
            for r in rows if r[cols["type_symbol"]] != "H"]


def wall_from_cavity(atoms):
    r = m151.analyse("x", [a[4] for a in atoms], [a[3] for a in atoms])
    if not r or r.get("main", 0) <= 0:
        return None
    clear, _, lo, shape = m151.clear_map([a[4] for a in atoms], [a[3] for a in atoms])
    mask = m151.enclosed(clear, lo, shape)
    if not mask.sum():
        return None
    cen = (np.argwhere(mask) * m151.SPACING + lo).mean(0)
    by = {}
    for a in atoms:
        by.setdefault(a[0], {})[a[2]] = a[4]
    return [k for k, d in by.items() if "CA" in d
            and np.linalg.norm(d["CA"] - cen) < 11.0]


def main():
    print(f"{'structure':<14}{'hull':>8}{'free(all)':>11}{'LARGEST':>9}"
          f"{'ncomp':>7}{'inflation':>11}")
    # PYR1, on its own 24 lining positions
    shp = {r["name"]: r for r in json.load(open(os.path.join(
        ROOT, "results", "pocket_shape", "pocket_shape.json")))}
    pat = [a for a in m143.read_pdb(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb"),
                                    want_atom=False) if a[1] in m143.THREE]
    pw = [(a[0], a[1], a[2], a[3], a[4]) for a in pat]
    r = measure([(a[0], a[1], a[2], a[3], a[4]) for a in pat], list(shp["PYR1"]["lining"]))
    tot, big, n, hull = r
    print(f"{'PYR1':<14}{hull:>8.0f}{tot:>11.0f}{big:>9.0f}{n:>7}"
          f"{tot/max(big,1):>10.2f}x")
    ref = big
    for b in ("batch01", "batch02", "batch03"):
        for f in sorted(glob.glob(os.path.join(ROOT, "results", "rfd3", b, "*.cif"))):
            at = read_cif(f)
            ch = {}
            for a in at:
                ch.setdefault(a[0], set()).add(a[1])
            dc = max(ch, key=lambda c: len(ch[c]) if len(ch[c]) < 260 else 0)
            dat = [(a[1], None, a[2], a[3], a[4]) for a in at if a[0] == dc]
            w = wall_from_cavity(dat)
            nm = b[-1] + "d" + os.path.basename(f).split("_model_")[1].replace(".cif", "")
            if not w:
                print(f"{nm:<14}{'':>8}{'':>11}{0:>9}{0:>7}      no cavity")
                continue
            r = measure(dat, w)
            if not r:
                continue
            tot, big, n, hull = r
            flag = "  <- vs PYR1 %+.0f" % (big - ref)
            print(f"{nm:<14}{hull:>8.0f}{tot:>11.0f}{big:>9.0f}{n:>7}"
                  f"{tot/max(big,1):>10.2f}x{flag}")
    print(f"\n  'free(all)' is §145's number (every disconnected void summed);")
    print(f"  'LARGEST' is the single connected pocket. 'inflation' is their ratio.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
