#!/usr/bin/env python
r"""
90_ti_validate.py -- does each TI leg actually represent what it claims?

Written after 88 returned a confident answer with the wrong sign. The lesson from
README 40/42 is that the failure is usually upstream of the estimator, so this
checks the SETUP before the number is allowed to mean anything:

  1. residue identity  -- sequential vs native numbering, and NOT merely "is it a
     VAL" : native 81 and native 83 are BOTH valine (85's assertion passes on
     either), so the mapping is confirmed against ligand contacts instead.
  2. ligand identity   -- copies, atom count and total charge per leg.
  3. starting sterics  -- minimum protein-ligand heavy-atom distance.
  4. final pose        -- ligand RMSD AFTER superposing on CA. Measuring it in the
     raw frame reports whole-box drift, not ligand motion (the same mistake this
     project already made once on gate/latch RMSD).

Writes results/ti/validation.txt.
"""
import os
import sys

import numpy as np
import parmed as pmd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
LEGS = ["V81I_aba", "V81I_mandi", "V81I_apo", "K59R_aba", "K59R_mandi", "K59R_apo"]
EXPECT = {"aba": (38, -1.0), "mandi": (51, 0.0), "apo": (0, 0.0)}


def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    D = np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))])
    return V @ D @ Wt, P.mean(0), Q.mean(0)


def main():
    out = []
    w = out.append

    w("=" * 78)
    w("1. RESIDUE IDENTITY  (native 81 and native 83 are both VAL -- the")
    w("   identity assertion in 85 cannot tell them apart, so anchor on contacts)")
    w("=" * 78)
    src = "results/stage1_rosetta/_input_wt_aba_A8S_anion.pdb"
    prot = [l for l in open(src) if l.startswith("ATOM")]
    het = [l for l in open(src) if l.startswith("HETATM")]
    order, seen, name = [], set(), {}
    for l in prot:
        r = int(l[22:26])
        name[r] = l[17:20].strip()
        if r not in seen:
            seen.add(r)
            order.append(r)
    lig = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in het])

    def mind(nat):
        c = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                      for l in prot if int(l[22:26]) == nat])
        return np.linalg.norm(c[:, None, :] - lig[None, :, :], axis=2).min()

    for seq, nat in ((79, 81), (81, 83), (59, 59)):
        got = order[seq - 1]
        ok = "OK" if got == nat else "MISMATCH"
        w(f"  sequential {seq:3d} -> native {got:3d} {name[got]:4s}  "
          f"(expected native {nat})  {ok}   min dist to ABA {mind(got):.2f} A")
    w("  -> the TI used sequential 79 = native 81, which is 4.85 A from ABA")
    w("     (second shell). Native 83 at 2.99 A is the one in direct contact.")

    w("")
    w("=" * 78)
    w("2. LIGAND IDENTITY PER LEG")
    w("=" * 78)
    for leg in LEGS:
        t = pmd.load_file(f"data/ti/{leg}/ti.prmtop")
        L = [r for r in t.residues if r.name == "LIG"]
        n = sum(len(r.atoms) for r in L)
        q = sum(a.charge for r in L for a in r.atoms)
        en, eq = EXPECT[leg.split("_")[1]]
        ok = "OK" if n == en and abs(q - eq) < 0.02 else "MISMATCH"
        ions = {i: sum(1 for r in t.residues if r.name == i) for i in ("K+", "Cl-")}
        w(f"  {leg:<12} LIG atoms {n:3d} (expect {en:3d})  q {q:+.3f} "
          f"(expect {eq:+.1f})  {ok}   {ions}")
    w("  NOTE: Cl- = 0 everywhere. These systems were NEUTRALISED only; they do")
    w("  NOT carry the project-canonical 0.15 M KCl. Largely cancels in the")
    w("  selectivity difference, but it is a deviation and is recorded as one.")

    w("")
    w("=" * 78)
    w("3. STARTING STERICS  and  4. FINAL POSE (ligand RMSD after CA superposition)")
    w("=" * 78)
    for leg in LEGS:
        if leg.endswith("apo"):
            continue
        top = f"data/ti/{leg}/ti.prmtop"
        t = pmd.load_file(top)
        ligh = [a.idx for a in t.atoms if a.residue.name == "LIG" and a.atomic_number > 1]
        ca = [a.idx for a in t.atoms if a.name == "CA"
              and a.residue.name not in ("WAT", "LIG", "K+", "Cl-")]
        poc = [a.idx for a in t.atoms
               if a.residue.name not in ("WAT", "HOH", "K+", "Cl-", "LIG")
               and a.atomic_number > 1]
        x0 = np.array(pmd.load_file(top, f"data/ti/{leg}/ti.inpcrd").coordinates)
        D0 = np.linalg.norm(x0[poc][:, None, :] - x0[ligh][None, :, :], axis=2)
        w(f"  {leg}:  START min prot-lig heavy = {D0.min():.2f} A   "
          f"contacts <2.0 A = {int((D0 < 2.0).sum())}   <2.5 A = {int((D0 < 2.5).sum())}")
        for i in range(0, 12, 2):
            f = f"data/ti/{leg}/lam{i:02d}/prod.rst7"
            if not os.path.exists(f):
                continue
            x = np.array(pmd.load_file(top, f).coordinates)
            R, mP, mQ = kabsch(x[ca], x0[ca])
            xf = (x - mP) @ R + mQ
            rl = np.sqrt(((xf[ligh] - x0[ligh]) ** 2).sum(1).mean())
            D = np.linalg.norm(xf[poc][:, None, :] - xf[ligh][None, :, :], axis=2)
            w(f"      lam{i:02d}  CA RMSD {np.sqrt(((xf[ca]-x0[ca])**2).sum(1).mean()):4.2f}"
              f"   ligand RMSD {rl:5.2f}   contacts <4 A {int((D < 4.0).sum()):4d}")

    txt = "\n".join(out)
    print(txt)
    os.makedirs("results/ti", exist_ok=True)
    open("results/ti/validation.txt", "w").write(txt + "\n")


if __name__ == "__main__":
    main()
