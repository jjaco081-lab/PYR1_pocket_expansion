"""
lib_cavity.py -- grid-based enclosed-cavity volume calculator.

Method
------
1. Build a cubic grid (default 0.5 A spacing) over a box centred on a reference
   point (normally the crystallographic ligand centroid).
2. A grid point is FREE if its distance to the nearest protein atom centre
   exceeds that atom's van der Waals radius + a solvent probe radius (1.4 A).
   This is a solvent-excluded-surface-like criterion.
3. A FREE point is BURIED if >= `bur_cut` (default 0.88) of 26 evenly spread
   ray directions hit protein within 15 A. This discriminates interior cavity
   from bulk solvent and from surface grooves.
4. FREE & BURIED points are grouped into connected components; the component
   containing (or nearest to) the reference point is reported.

Absolute volumes from this method are conservative relative to fpocket/CASTp
because of the strict buriedness cut and the single-component restriction.
They are intended for RELATIVE comparison between structures treated
identically -- do not quote the absolute number against literature values
computed with another tool.

Known behaviour: truncating a residue that SEALS the cavity from bulk solvent
lowers the reported volume (the component leaks and buriedness drops below
cut-off). A negative delta therefore flags a sealing residue, not a shrinking
cavity. This is diagnostic, not a bug.
"""
import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
BACKBONE = {"N", "CA", "C", "O", "CB"}

_DIRS = []
for _v in [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1),
           (1,-1,0),(1,0,-1),(0,1,-1),(1,1,-1),(1,-1,1),(-1,1,1)]:
    _v = np.array(_v, float); _v /= np.linalg.norm(_v)
    _DIRS += [_v, -_v]
DIRS = np.array(_DIRS)


def cavity_volume(atom_xyz, atom_elem, ref_point, spacing=0.5, probe=1.4,
                  bur_cut=0.88, box=14.0, ray_max=15.0, ray_step=0.75):
    """Return (volume_A3, n_gridpoints) of the enclosed cavity nearest ref_point.

    atom_xyz  : (N,3) float array of protein heavy-atom coordinates
    atom_elem : list/array of N element symbols
    ref_point : (3,) seed point, normally the ligand centroid
    """
    atom_xyz = np.asarray(atom_xyz, float)
    rad = np.array([VDW.get(e, 1.70) for e in atom_elem])
    ref = np.asarray(ref_point, float)

    lo, hi = ref - box, ref + box
    axes = [np.arange(lo[i], hi[i], spacing) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    shape = grid.shape[:3]
    pts = grid.reshape(-1, 3)

    tree = cKDTree(atom_xyz)
    d, i = tree.query(pts, k=1)
    clearance = d - rad[i]
    free = clearance > probe
    occ = (clearance < 0).reshape(shape)

    cand = np.where(free)[0]
    if cand.size == 0:
        return 0.0, 0
    cp = pts[cand]

    steps = np.arange(1.0, ray_max, ray_step)
    hits = np.zeros(len(cp))
    shp_arr = np.array(shape)
    for dv in DIRS:
        blocked = np.zeros(len(cp), bool)
        for s in steps:
            q = cp + dv * s
            idx = ((q - lo) / spacing).astype(int)
            ok = np.all((idx >= 0) & (idx < shp_arr), axis=1)
            b = np.zeros(len(cp), bool)
            v = idx[ok]
            b[ok] = occ[v[:, 0], v[:, 1], v[:, 2]]
            blocked |= b
        hits += blocked
    buried = hits / len(DIRS)

    mask = np.zeros(len(pts), bool)
    mask[cand[buried >= bur_cut]] = True
    m3 = mask.reshape(shape)
    lab, n = ndimage.label(m3)
    if n == 0:
        return 0.0, 0

    ci = np.clip(((ref - lo) / spacing).astype(int), 0, shp_arr - 1)
    l = lab[ci[0], ci[1], ci[2]]
    if l == 0:
        nz = np.array(np.nonzero(m3)).T
        l = lab[tuple(nz[np.argmin(np.linalg.norm(nz - ci, axis=1))])]
    npts = int((lab == l).sum())
    return npts * spacing ** 3, npts


def atoms_from_pose_pdb(path, chain=None, exclude_resnames=("HOH", "MN", "A8S"),
                        truncate_to_ala=()):
    """Parse a PDB file into (xyz, elements). truncate_to_ala is a set of
    residue numbers whose side chains beyond CB are dropped (in-silico Ala)."""
    xyz, el = [], []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resn = line[17:20].strip()
        if resn in exclude_resnames:
            continue
        ch = line[21]
        if chain and ch != chain:
            continue
        name = line[12:16].strip()
        rnum = int(line[22:26])
        if rnum in truncate_to_ala and name not in BACKBONE:
            continue
        e = line[76:78].strip() or name[0]
        if e == "H":
            continue
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        el.append(e)
    return np.array(xyz), el
