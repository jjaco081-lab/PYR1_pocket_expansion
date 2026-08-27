#!/usr/bin/env python
r"""
143_pocket_backbone_shape.py -- what does the FOLD give us, before side chains?

THE QUESTION
------------
README 78/79 left a categorical tension: every relative with a pocket bigger than
~300 A^3 sits at 3.4-4.0 A core RMSD and ~10 % identity from PYR1, while every
relative close enough to transplant (4DSB 1.03 A, 3OQU 1.46 A) has a pocket the
same size as PYR1's. Jannis asks the right follow-up: for those close relatives,
how do the pocket's AMINO ACID DISTRIBUTION and its BACKBONE SHAPE differ, since
either could be a route to a different binder even at constant volume.

Volume is one number and it hides two independent things:

  ENVELOPE   how much room the backbone encloses -- a property of the fold,
             changed only by moving the backbone (grafting, indels, sheet
             register). Measured here as the convex hull of the C-beta atoms of
             the lining residues.
  FILL       how much of that envelope the side chains occupy. Changed by
             mutation alone. Measured as 1 - cavity/envelope.

A pocket can be small because the fold is small (envelope-limited: mutation
cannot fix it) or because the fold is filled (fill-limited: mutation can). These
have opposite engineering consequences and the cavity volume alone does not
distinguish them. That is the whole point of this script.

WHAT IS MEASURED, per structure
-------------------------------
  cavity          the same unbiased largest-enclosed-component volume as
                  script 08, so numbers reproduce results/homolog_cavities.csv
  shape           extents of the cavity point cloud along its own principal
                  axes: L1 >= L2 >= L3, and elongation L1/L3. A 200 A^3 tube and
                  a 200 A^3 sphere accept completely different ligands.
  envelope        C-beta convex hull of the lining residues (side-chain free)
  fill            1 - cavity/envelope
  poly-Gly        cavity recomputed with every lining side chain deleted past
                  CA, seeded at the WT cavity centroid. This is the ceiling that
                  mutation alone can reach.
  composition     lining residues by chemical class, and how many line the wall
                  through BACKBONE atoms only (those are envelope, not fill --
                  no mutation moves them)
  secondary str.  which secondary structure the wall is built from, by P-SEA
                  CA-geometry criteria. "Add residues to the beta sheet" is only
                  a coherent idea if the sheet is actually lining the pocket.

CAVEATS, STATED BEFORE THE NUMBERS
----------------------------------
* lib_cavity's docstring warns that truncating a SEALING residue makes the
  component leak to bulk and the reported volume FALL. The poly-Gly arm is
  exactly that situation, so it carries an explicit `leaked` flag: the component
  is checked for contact with the grid boundary. A leaked cell is reported as
  LEAK, never as a volume.
* Homolog cavities are measured APO (script 08's reader keeps group_PDB == ATOM
  only, so ligands are excluded). PYR1 is treated identically -- A8S is stripped
  -- or the comparison would be against a filled pocket.
* Absolute volumes are conservative (strict buriedness cut). Compare within this
  table, not against fpocket/CASTp literature values.
* P-SEA is a CA-distance heuristic, not DSSP. It is validated here against PYR1's
  known C-terminal grip helix before any homolog is trusted.

Run: python scripts/143_pocket_backbone_shape.py   (numpy + scipy only, ~1 GB)
"""
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree, ConvexHull

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_cavity import VDW, DIRS                                  # noqa: E402

HC = os.path.join(ROOT, "results", "homolog_cavities")
OUT = os.path.join(ROOT, "results", "pocket_shape")
os.makedirs(OUT, exist_ok=True)

BACKBONE = {"N", "CA", "C", "O", "OXT"}
CLASS = {}
for _aas, _c in [("AVLIMPG", "aliphatic"), ("FWY", "aromatic"),
                 ("STNQCH", "polar"), ("KR", "basic"), ("DE", "acidic")]:
    for _a in _aas:
        CLASS[_a] = _c
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
         "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
         "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V", "MSE": "M"}

#: PYR1 wall positions with their asserted identities (README 4). Never trust a
#: residue number: the script stops if the frame does not carry these residues.
PYR1_WALL = {59: "K", 79: "R", 94: "E", 108: "F", 120: "Y", 141: "E"}
#: PYR1's C-terminal grip helix, used to validate the P-SEA assignment.
GRIP_HELIX = (160, 178)

#: the graft-geometry set (script 79), plus PYR1 itself
TARGETS = ["2pcsA00", "2ns9B01", "2bk0A00", "6awvC00", "3qrzA00", "3tfzE00",
           "4xrwA02", "3oquB00", "2flhA00", "3oh8A01", "4dsbB00", "1z94E00",
           "4n3eY00", "4xrtB02", "5nonC00", "2resA00", "3p51A00", "4xrwA01",
           "3qrzB00"]


# ---------------------------------------------------------------- readers ----
def read_pdb(path, want_atom=True):
    """-> list of (resnum, resname, atomname, elem, xyz). Heavy atoms only."""
    out = []
    for line in open(path):
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        if want_atom and rec != "ATOM":
            continue
        alt = line[16]
        if alt not in (" ", "A"):
            continue
        e = (line[76:78].strip() or line[12:16].strip()[0])
        if e == "H":
            continue
        out.append((int(line[22:26]), line[17:20].strip(), line[12:16].strip(), e,
                    np.array([float(line[30:38]), float(line[38:46]),
                              float(line[46:54])])))
    return out


def read_cif_resnames(path, chain):
    """resnum -> one-letter code, for the protein atoms of one chain."""
    cols, rows, in_loop = {}, [], False
    for line in open(path):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".")[1]] = len(cols)
            in_loop = True
            continue
        if in_loop:
            if line.startswith(("#", "loop_", "_")):
                if rows:
                    break
                continue
            f = line.split()
            if len(f) >= len(cols):
                rows.append(f)
    if "auth_seq_id" not in cols:
        return {}
    out = {}
    for f in rows:
        if f[cols["group_PDB"]] != "ATOM" or f[cols["auth_asym_id"]] != chain:
            continue
        try:
            out[int(f[cols["auth_seq_id"]])] = THREE.get(
                f[cols["label_comp_id"]].strip(), "X")
        except ValueError:
            continue
    return out


# ---------------------------------------------------------------- cavity -----
def cavity(xyz, elem, spacing=0.6, probe=1.4, bur_cut=0.88, ray_max=15.0,
           ray_step=0.75, seed=None):
    """Enclosed cavity as script 08. Returns (volume, points, leaked).

    seed=None  -> the largest enclosed component anywhere (script 08 behaviour)
    seed=(3,)  -> the component nearest that point, so a truncation arm keeps
                  reporting on the SAME pocket instead of jumping to another one
    leaked     -> the chosen component touches the grid boundary, i.e. it is no
                  longer enclosed. Its volume is meaningless; see the docstring.
    """
    xyz = np.asarray(xyz, float)
    rad = np.array([VDW.get(e, 1.70) for e in elem])
    lo = xyz.min(0) - 3
    axes = [np.arange(lo[i], xyz.max(0)[i] + 3, spacing) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    shape = grid.shape[:3]
    pts = grid.reshape(-1, 3)
    d, i = cKDTree(xyz).query(pts, k=1)
    clear = d - rad[i]
    occ = (clear < 0).reshape(shape)
    cand = np.where(clear > probe)[0]
    if cand.size == 0:
        return 0.0, np.zeros((0, 3)), False
    cp = pts[cand]
    hits = np.zeros(len(cp))
    shp = np.array(shape)
    for dv in DIRS:
        blocked = np.zeros(len(cp), bool)
        for s in np.arange(1.0, ray_max, ray_step):
            idx = ((cp + dv * s - lo) / spacing).astype(int)
            ok = np.all((idx >= 0) & (idx < shp), axis=1)
            v = idx[ok]
            b = np.zeros(len(cp), bool)
            b[ok] = occ[v[:, 0], v[:, 1], v[:, 2]]
            blocked |= b
        hits += blocked
    mask = np.zeros(len(pts), bool)
    mask[cand[(hits / len(DIRS)) >= bur_cut]] = True
    lab, n = ndimage.label(mask.reshape(shape))
    if n == 0:
        return 0.0, np.zeros((0, 3)), False
    labf = lab.reshape(-1)
    if seed is None:
        sizes = ndimage.sum(mask.reshape(shape), lab, range(1, n + 1))
        pick = int(np.argmax(sizes)) + 1
    else:
        sel = np.where(labf > 0)[0]
        pick = int(labf[sel[np.argmin(((pts[sel] - seed) ** 2).sum(1))]])
    keep = labf == pick
    P = pts[keep]
    # leak test: does this component reach the outer shell of the box?
    idx = np.argwhere(lab == pick)
    leaked = bool((idx.min(0) <= 0).any() or (idx.max(0) >= shp - 1).any())
    return float(keep.sum() * spacing ** 3), P, leaked


def extents(P):
    """principal-axis extents L1>=L2>=L3 of a point cloud."""
    if len(P) < 4:
        return (0.0, 0.0, 0.0)
    Q = P - P.mean(0)
    V = np.linalg.svd(Q, full_matrices=False)[2]
    pr = Q @ V.T
    return tuple(sorted((pr.max(0) - pr.min(0)).tolist(), reverse=True))


# ------------------------------------------------------------------ P-SEA ----
def psea(ca):
    """CA-only secondary structure (Labesse et al. 1997 distance criteria).

    ca : (N,3) in chain order. -> string of H / E / C, length N.
    Helix  d(i,i+2)=5.5+-0.5  d(i,i+3)=5.3+-0.5  d(i,i+4)=6.4+-0.6
    Strand d(i,i+2)=6.7+-0.6  d(i,i+3)=9.9+-0.9  d(i,i+4)=12.4+-1.1
    """
    n = len(ca)
    s = ["C"] * n
    def dist(i, k):
        return float(np.linalg.norm(ca[i + k] - ca[i])) if i + k < n else 1e9
    for i in range(n - 4):
        d2, d3, d4 = dist(i, 2), dist(i, 3), dist(i, 4)
        if abs(d2 - 5.5) <= 0.5 and abs(d3 - 5.3) <= 0.5 and abs(d4 - 6.4) <= 0.6:
            for j in range(i, i + 5):
                s[j] = "H"
    for i in range(n - 4):
        d2, d3, d4 = dist(i, 2), dist(i, 3), dist(i, 4)
        if (abs(d2 - 6.7) <= 0.6 and abs(d3 - 9.9) <= 0.9
                and abs(d4 - 12.4) <= 1.1 and s[i] != "H"):
            for j in range(i, i + 5):
                if s[j] != "H":
                    s[j] = "E"
    return "".join(s)


# ------------------------------------------------------------- per-structure -
def analyse(name, atoms, resname, verbose=False):
    """atoms: [(resnum, resname3, atomname, elem, xyz)] protein heavy atoms."""
    xyz = np.array([a[4] for a in atoms])
    elem = [a[3] for a in atoms]
    vol, P, leak = cavity(xyz, elem)
    if len(P) == 0:
        return None
    cen = P.mean(0)

    # lining residues: any heavy atom within 4.0 A of a cavity grid point
    tree = cKDTree(P)
    near = tree.query_ball_point(xyz, 4.0)
    lining, bb_only = set(), {}
    for a, nb in zip(atoms, near):
        if not nb:
            continue
        lining.add(a[0])
        bb_only.setdefault(a[0], True)
        if a[2] not in BACKBONE:
            bb_only[a[0]] = False
    lining = sorted(lining)

    # C-beta envelope of the lining wall (CA where there is no CB)
    cb = {}
    for a in atoms:
        if a[0] in bb_only and a[2] in ("CA", "CB"):
            if a[2] == "CB" or a[0] not in cb:
                cb[a[0]] = a[4]
    E = np.array([cb[r] for r in lining if r in cb])
    env = float(ConvexHull(E).volume) if len(E) >= 4 else float("nan")

    # poly-Gly ceiling: delete every lining side chain past CA, same pocket
    keep = [(a, e) for a, e in zip(atoms, elem)
            if not (a[0] in set(lining) and a[2] not in BACKBONE)]
    gvol, _, gleak = cavity(np.array([a[4] for a, _ in keep]),
                            [e for _, e in keep], seed=cen)

    # composition
    comp = Counter()
    for r in lining:
        aa = resname.get(r, "X")
        comp[CLASS.get(aa, "other")] += 1
    seq = "".join(resname.get(r, "X") for r in lining)

    # secondary structure of the wall
    ca_order = [(a[0], a[4]) for a in atoms if a[2] == "CA"]
    ss = psea(np.array([c for _, c in ca_order]))
    ssn = {rn: ss[i] for i, (rn, _) in enumerate(ca_order)}
    sscount = Counter(ssn.get(r, "C") for r in lining)

    L = extents(P)
    LE = extents(E) if len(E) >= 4 else (0, 0, 0)
    return dict(name=name, cavity=round(vol, 1), leaked=leak,
                n_lining=len(lining),
                n_bb_only=sum(1 for r in lining if bb_only[r]),
                envelope=round(env, 1),
                fill=round(1 - vol / env, 3) if env == env and env else None,
                polygly=round(gvol, 1), polygly_leaked=gleak,
                L1=round(L[0], 1), L2=round(L[1], 1), L3=round(L[2], 1),
                elong=round(L[0] / L[2], 2) if L[2] else None,
                E1=round(LE[0], 1), E2=round(LE[1], 1), E3=round(LE[2], 1),
                comp={k: comp[k] for k in
                      ("aliphatic", "aromatic", "polar", "basic", "acidic")},
                ss={k: sscount.get(k, 0) for k in "HEC"},
                lining=lining, seq=seq, ss_string=ss,
                ca_resnums=[rn for rn, _ in ca_order])


def main():
    # ---- PYR1 reference, treated exactly like a homolog (ligand stripped) ----
    frame = os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")
    allat = read_pdb(frame, want_atom=False)
    assert any(a[1] == "A8S" for a in allat), "wt_aba.pdb has lost its ligand"
    prot = [a for a in allat if a[1] in THREE]
    rn1 = {}
    for a in prot:
        rn1[a[0]] = THREE[a[1]]
    bad = {p: (rn1.get(p), aa) for p, aa in PYR1_WALL.items() if rn1.get(p) != aa}
    assert not bad, f"PYR1 frame residue identity check FAILED: {bad}"
    print(f"PYR1 frame {os.path.basename(frame)}: {len(rn1)} residues, "
          f"wall identities {''.join(PYR1_WALL[p] for p in sorted(PYR1_WALL))} "
          f"at {sorted(PYR1_WALL)} CONFIRMED; A8S stripped for the cavity")

    rows = []
    r = analyse("PYR1", prot, rn1)
    # validate P-SEA against the known grip helix before trusting any homolog
    ss, nums = r["ss_string"], r["ca_resnums"]
    idx = [i for i, n in enumerate(nums) if GRIP_HELIX[0] <= n <= GRIP_HELIX[1]]
    fh = sum(ss[i] == "H" for i in idx) / max(len(idx), 1)
    print(f"P-SEA validation: PYR1 grip helix {GRIP_HELIX[0]}-{GRIP_HELIX[1]} "
          f"called {fh:.0%} helix")
    assert fh >= 0.70, "P-SEA failed on a known helix; do not trust the SS column"
    rows.append(r)

    for t in TARGETS:
        dom = os.path.join(HC, "domains", t + ".pdb")
        cif = os.path.join(HC, "raw", t + ".cif")
        if not (os.path.exists(dom) and os.path.exists(cif)):
            print(f"  {t}: missing input, skipped")
            continue
        atoms = read_pdb(dom)                      # UNK resnames, right coords
        names = read_cif_resnames(cif, t[4])       # 2pcsA00 -> chain A
        cov = sum(1 for a in atoms if a[2] == "CA" and a[0] in names)
        nca = sum(1 for a in atoms if a[2] == "CA")
        if cov < 0.95 * nca:
            print(f"  {t}: resname map covers {cov}/{nca} CA, skipped")
            continue
        atoms = [(a[0], names.get(a[0], "UNK"), a[2], a[3], a[4]) for a in atoms]
        rr = analyse(t.upper()[:4] + "_" + t[4], atoms, names)
        if rr:
            rows.append(rr)
        print(f"  {t}: cavity {rr['cavity']:.0f} A^3, {rr['n_lining']} lining",
              flush=True)

    json.dump(rows, open(os.path.join(OUT, "pocket_shape.json"), "w"), indent=1)

    # ------------------------------------------------------------- report ----
    rows.sort(key=lambda x: -x["cavity"])
    print("\n" + "=" * 96)
    print("ENVELOPE vs FILL -- is the pocket small because of the fold, or the side chains?")
    print("=" * 96)
    print(f"{'structure':<10}{'cavity':>8}{'envelope':>10}{'fill':>7}"
          f"{'polyGly':>10}{'lin':>5}{'bbOnly':>8}   {'L1':>5}{'L2':>5}{'L3':>5}{'elong':>7}")
    for x in rows:
        pg = "LEAK" if x["polygly_leaked"] else f"{x['polygly']:.0f}"
        print(f"{x['name']:<10}{x['cavity']:>8.0f}{x['envelope']:>10.0f}"
              f"{x['fill']:>7.2f}{pg:>10}{x['n_lining']:>5}{x['n_bb_only']:>8}"
              f"   {x['L1']:>5.1f}{x['L2']:>5.1f}{x['L3']:>5.1f}{x['elong']:>7.2f}")

    print("\n" + "=" * 96)
    print("WALL COMPOSITION -- lining residues by chemical class, and by secondary structure")
    print("=" * 96)
    print(f"{'structure':<10}{'alip':>6}{'arom':>6}{'polar':>7}{'basic':>7}{'acid':>6}"
          f"    {'helix':>6}{'sheet':>6}{'loop':>6}   lining sequence")
    for x in rows:
        c, s = x["comp"], x["ss"]
        print(f"{x['name']:<10}{c['aliphatic']:>6}{c['aromatic']:>6}{c['polar']:>7}"
              f"{c['basic']:>7}{c['acidic']:>6}    {s['H']:>6}{s['E']:>6}{s['C']:>6}"
              f"   {x['seq']}")
    print(f"\nwritten to {OUT}/pocket_shape.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
