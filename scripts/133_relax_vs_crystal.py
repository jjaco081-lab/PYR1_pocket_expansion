#!/usr/bin/env python
r"""
133_relax_vs_crystal.py -- does FastRelax reproduce the EXPERIMENTAL pocket when
handed the right sequence and the right ligand?

THE DIAGNOSTIC (Jannis, 2026-08-26). Everything ligand-present we have run --
FastDesign (§23j), coupled moves (§36), MM-GBSA (§37-42), TI (§43-46) -- asked the
protocol to RANK or DESIGN. None asked whether the modelling reproduces the
crystal. That separates two failure modes:

  A  start from the CRYSTAL rotamers, relax -- do they STAY?
       drift -> the SCORE prefers something else          (a scoring problem)
       stay  -> the score is content with the real answer
  B  start from REPACKED rotamers, relax -- do they REACH the crystal?
       no    -> a SAMPLING problem (§33b: crystal Arg59 sits at chi3 ~100 deg,
                which the backbone-independent rotamer library never proposes)
       yes   -> modelling is adequate and our failures live in RANKING

FRAME. 4WVO cannot be loaded directly: its 3UZ carries CCD atom names (CAC, CAB,
...) while the §49 params use C1/C12/N1/..., and hand-renaming is the silent
mismatch §97 exists to prevent. §49 already married the CRYSTAL coordinates to the
params names, so `stage1/params/3UZ_0001.pdb` IS the crystal ligand in the stage-1
frame -- verified here by comparing its first atom against wt_mandi.pdb's CAC.

⚠ PDB LINES ARE SPLICED FROM REAL ONES, never formatted by hand (§101). A
hand-written writer put resName in columns 17-19 and chainID in 21 -- everything
after the atom name shifted one left -- and Rosetta died with "too many tries in
fill_missing_atoms" on a file that looked entirely normal in an editor.

⚠ This measures STRUCTURAL accuracy only. A protocol can rebuild a pocket
faithfully and still rank mutations badly; those are different claims.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_rosetta as LR                                        # noqa: E402

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "relax_vs_crystal")
FRAME = os.path.join(ROOT, "data", "stage1", "wt_mandi.pdb")
LIGPDB = os.path.join(ROOT, "data", "stage1", "params", "3UZ_0001.pdb")
LIGPARAMS = os.path.join(ROOT, "data", "stage1", "params", "3UZ.params")
CIF = os.path.join(ROOT, "data", "4WVO.cif")
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
BB = {"N", "CA", "C", "O", "OXT"}
MUT = {59: "ARG", 81: "ILE", 108: "ALA", 159: "LEU"}
POCKET = [59, 81, 83, 87, 89, 92, 94, 108, 110, 117, 120, 122, 141, 159, 160,
          163, 164, 167]


def read_cif(path):
    prot, seq = {}, {}
    cols, inloop = {}, False
    for line in open(path):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".", 1)[1]] = len(cols)
            inloop = True
            continue
        if inloop and line.startswith("#"):
            inloop = False
            continue
        if not inloop or not cols or not line.startswith("ATOM"):
            continue
        p = line.split()
        if p[cols["auth_asym_id"]] != "A" or p[cols["label_alt_id"]] not in (".", "A"):
            continue
        atom = p[cols["label_atom_id"]].strip('"')
        if atom.startswith("H"):
            continue
        num = int(p[cols["auth_seq_id"]])
        prot[(num, atom)] = np.array([float(p[cols["Cartn_x"]]),
                                      float(p[cols["Cartn_y"]]),
                                      float(p[cols["Cartn_z"]])])
        seq[num] = p[cols["label_comp_id"]]
    return prot, seq


def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, _S, Wt = np.linalg.svd(Pc.T @ Qc)
    D = np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))])
    return V @ D @ Wt, P.mean(0), Q.mean(0)


def main():
    os.makedirs(OUT, exist_ok=True)
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    frame = [l for l in open(FRAME) if l.startswith("ATOM")]
    ligl = [l for l in open(LIGPDB) if l.startswith("HETATM")]
    say("=" * 78)
    say("0. FRAME AND IDENTITY")
    say("=" * 78)
    fx = np.array([float(frame[0][30:38]), float(frame[0][38:46]),
                   float(frame[0][46:54])])
    lx = np.array([float(ligl[0][30:38]), float(ligl[0][38:46]),
                   float(ligl[0][46:54])])
    wm = [l for l in open(FRAME) if l.startswith("HETATM")]
    wx = np.array([float(wm[0][30:38]), float(wm[0][38:46]), float(wm[0][46:54])])
    if np.linalg.norm(lx - wx) > 0.01:
        raise SystemExit(f"params ligand is not the crystal pose: {lx} vs {wx}")
    say(f"   stage-1 frame: {sum(1 for l in frame if l[12:16].strip()=='CA')} residues")
    say(f"   params ligand matches wt_mandi's 3UZ to "
        f"{np.linalg.norm(lx-wx):.4f} A -- it IS the crystal pose")

    xt, xseq = read_cif(CIF)
    bad = [(n, w, xseq.get(n)) for n, w in MUT.items() if xseq.get(n) != w]
    if bad:
        raise SystemExit(f"4WVO mutation check failed: {bad}")
    # superpose 4WVO backbone onto the frame
    fbb = {}
    for l in frame:
        nm = l[12:16].strip()
        if nm in ("N", "CA", "C", "O"):
            fbb[(int(l[22:26]), nm)] = np.array([float(l[30:38]), float(l[38:46]),
                                                 float(l[46:54])])
    com = sorted(set(fbb) & set(xt))
    if len(com) < 400:
        raise SystemExit(f"only {len(com)} shared backbone atoms")
    P = np.array([xt[k] for k in com])
    Q = np.array([fbb[k] for k in com])
    R, mP, mQ = kabsch(P, Q)
    rms = float(np.sqrt((((P - mP) @ R + mQ - Q) ** 2).sum(1).mean()))
    say(f"   4WVO -> frame: {len(com)} backbone atoms, RMSD {rms:.2f} A")
    say(f"   mutations confirmed in 4WVO: "
        + ", ".join(f"{n}{AA3[MUT[n]]}" for n in sorted(MUT)))

    def xtal_sc(num):
        return {a: (v - mP) @ R + mQ for (n, a), v in xt.items()
                if n == num and a not in BB}

    # ---- build arm A: frame backbone + 4WVO side chains at the 4 positions ----
    # every line is SPLICED from a real ATOM line; only coordinates and resName
    # are substituted, and the column layout of the source line is preserved
    def build(with_xtal_sc):
        out = []
        n = 0
        for l in frame:
            num = int(l[22:26])
            nm = l[12:16].strip()
            if num in MUT:
                if nm not in BB:
                    continue                     # drop the wild-type side chain
                n += 1
                out.append(l[:17] + f"{MUT[num]:>3s}" + l[20:])
            else:
                n += 1
                out.append(l)
        if with_xtal_sc:
            for num in sorted(MUT):
                sc = xtal_sc(num)
                # reuse a real side-chain line from this residue as the template
                tmpl = next(l for l in frame
                            if int(l[22:26]) == num and l[12:16].strip() not in BB)
                for a, v in sorted(sc.items()):
                    nm4 = f" {a:<3s}" if len(a) < 4 else a
                    out.append(tmpl[:12] + nm4 + tmpl[16:17] + f"{MUT[num]:>3s}"
                               + tmpl[20:30]
                               + f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}" + tmpl[54:])
        out.sort(key=lambda l: (int(l[22:26]), l[12:16]))
        return out

    import pyrosetta
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ex1 -ex2aro "
                   f"-extra_res_fa {LIGPARAMS}")
    sfxn = pyrosetta.create_score_function("ref2015_cart")

    pa_file = os.path.join(OUT, "quad_xtal.pdb")
    with open(pa_file, "w") as fh:
        fh.writelines(build(True))
        fh.write("TER\n")
        fh.writelines(ligl)
        fh.write("END\n")
    pose = pyrosetta.pose_from_pdb(pa_file)
    info = pose.pdb_info()
    idx = {info.number(i): i for i in range(1, pose.total_residue() + 1)}
    ligres = pose.total_residue()
    if pose.residue(ligres).name3().strip() != "LIG":
        raise SystemExit(f"last residue {pose.residue(ligres).name3()}, not LIG")
    say(f"   pose built: {pose.total_residue()} residues incl. ligand")
    for num, want in MUT.items():
        got = pose.residue(idx[num]).name3()
        if got != want:
            raise SystemExit(f"position {num} is {got}, expected {want}")
    say(f"   all four mutations present in the pose")

    sel = ResidueIndexSelector(",".join(str(idx[p]) for p in POCKET if p in idx))
    nb = NeighborhoodResidueSelector(sel, 6.0, True)
    allowed = sorted({i + 1 for i, b in enumerate(nb.apply(pose)) if b} | {ligres})
    say(f"   relaxable shell: {len(allowed)} residues, ligand free")

    def sc(pose_, num):
        r = pose_.residue(idx[num])
        d = {}
        for j in range(1, r.natoms() + 1):
            nm = r.atom_name(j).strip()
            if nm in BB or r.atom_type(j).element().strip() == "H":
                continue
            v = r.xyz(j)
            d[nm] = np.array([v.x, v.y, v.z])
        return d

    ref = {p: sc(pose, p) for p in POCKET if p in idx}

    def dev(pose_, p):
        a, b = ref[p], sc(pose_, p)
        k = set(a) & set(b)
        return float(np.sqrt(np.mean([np.sum((a[x] - b[x]) ** 2) for x in k]))) if k else None

    def relax(p0):
        mm = MoveMap()
        mm.set_bb(False)
        mm.set_chi(False)
        for i in allowed:
            mm.set_bb(i, True)
            mm.set_chi(i, True)
        tf, _ = LR.restrict_packing(p0, allowed)
        fr = FastRelax(sfxn, 1)
        fr.cartesian(True)
        fr.min_type("lbfgs_armijo_nonmonotone")
        fr.set_movemap(mm)
        fr.set_task_factory(tf)
        fr.apply(p0)
        return p0

    say("")
    say("=" * 78)
    say("RESULT: side-chain heavy-atom RMSD to the crystal conformation")
    say("=" * 78)
    pa = relax(pose.clone())
    pb = pose.clone()
    tf, _ = LR.restrict_packing(pb, allowed)
    pk = PackRotamersMover(sfxn)
    pk.task_factory(tf)
    pk.apply(pb)
    b0 = {p: dev(pb, p) for p in ref}
    pb = relax(pb)
    ra = {p: dev(pa, p) for p in ref}
    rb = {p: dev(pb, p) for p in ref}
    say(f"   {'pos':<6}{'res':<5}{'A: from xtal':>14}{'B: repacked':>13}"
        f"{'B relaxed':>12}")
    for p in sorted(ref):
        tag = "  <-- PYR1^MANDI" if p in MUT else ""
        say(f"   {p:<6}{pose.residue(idx[p]).name3():<5}{ra[p]:>14.2f}"
            f"{b0[p]:>13.2f}{rb[p]:>12.2f}{tag}")
    mm_ = [p for p in MUT if p in ref]
    say("")
    say(f"   mean over all {len(ref)} pocket positions:  "
        f"A {np.mean([ra[p] for p in ref]):.2f}   B {np.mean([rb[p] for p in ref]):.2f}")
    say(f"   mean over the four mutated positions:  "
        f"A {np.mean([ra[p] for p in mm_]):.2f}   B {np.mean([rb[p] for p in mm_]):.2f}")
    with open(os.path.join(OUT, "relax_vs_crystal.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    with open(os.path.join(OUT, "rmsd.json"), "w") as fh:
        json.dump({"A_from_crystal": ra, "B_repacked": b0, "B_relaxed": rb}, fh,
                  indent=1)
    say(f"\n   written to {OUT}/relax_vs_crystal.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
