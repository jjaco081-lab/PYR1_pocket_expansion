#!/usr/bin/env python
r"""
164_pose_geometry.py -- do tractable ligands sit differently in the WT pocket?

Jannis's idea. §82 showed Boltz-2's AFFINITY head cannot rank pocket variants,
and §84's tractability benchmark asks the ligand-level question through that same
head. This asks it through GEOMETRY instead: co-folded with wild-type PYR1, does
a ligand that PYR1 can be evolved to bind end up positioned differently from a
matched ligand that it cannot?

That is worth asking precisely because it does not go through the affinity head.
The position is a property of the structure module, and it is free -- 348 models
already exist.

FEATURES, fixed before any result is seen. Each is computed after superposing the
model's core backbone onto the 3QN1-derived PYR1 frame (166 residues, gate 85-89
and latch 116-117 excluded so ligand-driven gate motion cannot drive the fit):

  d_site      ligand centroid to the crystallographic ABA centroid. Small = the
              ligand is where ABA binds.
  d_mouth     ligand centroid to the gate/latch mouth (centroid of 85-89 +
              116-117). Jannis's hypothesis is about THIS one: tractable ligands
              may sit closer to the opening.
  buried      fraction of ligand heavy atoms with >= 16 protein heavy atoms
              within 8 A -- how enclosed it is, independent of where.
  frac_in     fraction of ligand atoms within 8 A of the ABA centroid.

⚠ FOLD IS ASSERTED FIRST. §92c was retracted because 30 structures with median
pLDDT 49.4 and a 14 A core RMSD were measured as if they were folded. Every model
here is checked for median pLDDT >= 80 and core RMSD <= 3 A, and the count that
fails is reported rather than silently dropped.

⚠ Both classes are folded against the SAME wild-type PYR1 (§84b). A hit is a
ligand some VARIANT bound, so a geometric difference in the WT pocket is a proxy
one step from the label -- but it is symmetric between the classes.
"""
import glob, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
T = os.path.join(ROOT, "data", "tractability")
BB = ("N", "CA", "C", "O")
GATE, LATCH = [85, 86, 87, 88, 89], [116, 117]


def kabsch(P, Q):
    mp, mq = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((Q - mq).T @ (P - mp))
    R = U @ np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))]) @ Vt
    return R, mp, mq


def frame():
    bb, het, allat = {}, [], []
    for l in open(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")):
        if l.startswith("ATOM"):
            xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
            allat.append(xyz)
            if l[12:16].strip() in BB:
                bb.setdefault(int(l[22:26]), {})[l[12:16].strip()] = xyz
        elif l.startswith("HETATM") and l[17:20].strip() == "A8S":
            het.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return bb, np.array(het), np.array(allat)


def read(path):
    cols, rows, inl = {}, [], False
    for line in open(path):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".")[1]] = len(cols); inl = True; continue
        if inl:
            if line.startswith(("#", "loop_", "_")):
                if rows: break
                continue
            f = line.split()
            if len(f) >= len(cols): rows.append(f)
    bb, lig, plddt, prot = {}, [], [], []
    for f in rows:
        if f[cols["type_symbol"]] == "H": continue
        xyz = [float(f[cols[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")]
        if f[cols["group_PDB"]] == "ATOM":
            prot.append(xyz)
            an = f[cols["label_atom_id"]]
            if an in BB:
                bb.setdefault(int(f[cols["label_seq_id"]]), {})[an] = xyz
            if an == "CA":
                plddt.append(float(f[cols["B_iso_or_equiv"]]))
        else:
            lig.append(xyz)
    return bb, np.array(lig), np.array(prot), (float(np.median(plddt)) if plddt else 0.0)


def paired(a, b, resl):
    P, Q = [], []
    for r in sorted(resl):
        if r in a and r in b:
            for x in BB:
                if x in a[r] and x in b[r]:
                    P.append(a[r][x]); Q.append(b[r][x])
    return np.array(P), np.array(Q)


def auc(v, lab):
    v = np.asarray(v, float); lab = np.asarray(lab, bool)
    p, n = v[lab], v[~lab]
    if not len(p) or not len(n): return float("nan")
    return float(np.mean([[1.0 if a < b else .5 if a == b else 0. for b in n] for a in p]))


def main():
    fbb, aba, fall = frame()
    core = [r for r in fbb if 6 <= r <= 180 and r not in GATE + LATCH]
    mouth = np.array([fbb[r]["CA"] for r in GATE + LATCH if r in fbb]).mean(0)
    key = json.load(open(os.path.join(T, "key_matched.json")))

    rows, bad = [], {"plddt": 0, "core": 0, "nolig": 0}
    for lig, meta in key.items():
        g = sorted(glob.glob(os.path.join(T, "out", f"boltz_results_{lig}",
                                          "predictions", "*", "*.cif")))
        if not g: continue
        bb, L, prot, pl = read(g[0])
        if len(L) == 0: bad["nolig"] += 1; continue
        if pl < 80: bad["plddt"] += 1; continue
        P, Q = paired(fbb, bb, core)
        if len(P) < 400: bad["core"] += 1; continue
        R, mp, mq = kabsch(P, Q)
        crms = float(np.sqrt((((Q - mq) @ R + mp - P) ** 2).sum(1).mean()))
        if crms > 3.0: bad["core"] += 1; continue
        Lf = (L - mq) @ R + mp
        Pf = (prot - mq) @ R + mp
        c = Lf.mean(0)
        d = np.linalg.norm(Lf[:, None] - Pf[None], axis=-1)
        rows.append(dict(lig=lig, name=meta["name"], hit=bool(meta["label"]),
                         heavy=meta["heavy"], plddt=pl, core_rmsd=crms,
                         d_site=float(np.linalg.norm(c - aba.mean(0))),
                         d_mouth=float(np.linalg.norm(c - mouth)),
                         buried=float(((d < 8).sum(1) >= 16).mean()),
                         frac_in=float((np.linalg.norm(Lf - aba.mean(0), axis=1) < 8).mean())))
    lab = [r["hit"] for r in rows]
    print(f"{len(rows)} models pass the fold assertions "
          f"({sum(lab)} hits, {len(lab)-sum(lab)} non-hits)")
    print(f"  rejected: pLDDT<80 {bad['plddt']}, core RMSD>3 A {bad['core']}, "
          f"no ligand {bad['nolig']}")
    if len(set(lab)) < 2:
        print("  only one class present"); return 0
    print(f"\n{'feature':<12}{'hits (median)':>15}{'non-hits':>12}{'AUC':>8}"
          f"{'  direction'}")
    for f, lo_is_hit, desc in (("d_site", True, "closer to the ABA site"),
                               ("d_mouth", True, "closer to the mouth"),
                               ("buried", False, "more enclosed"),
                               ("frac_in", False, "more atoms in the site"),
                               ("heavy", True, "smaller")):
        v = np.array([r[f] for r in rows], float)
        a = auc(v if lo_is_hit else -v, lab)
        h, n = np.median(v[np.array(lab)]), np.median(v[~np.array(lab)])
        print(f"{f:<12}{h:>15.2f}{n:>12.2f}{a:>8.3f}   {desc}")
    print("\n  AUC > 0.5 means the stated direction is associated with being a HIT.")
    json.dump(rows, open(os.path.join(ROOT, "results", "tractability",
                                      "pose_geometry.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
