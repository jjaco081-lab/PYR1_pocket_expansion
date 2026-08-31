#!/usr/bin/env python
r"""
163_pose_stability.py -- is a predicted pose stable enough to score against?

§91 scored 475 variants against a CRYSTAL mandipropamid pose and got AUC 0.868.
The three other compounds in Park's screen have no crystal structure, so the
whole question of whether this method reaches novel ligands (Goal 2/3) reduces
to one measurement: does Boltz-2 put the ligand in the SAME PLACE across seeds?

Ten seeds per compound. For each pair of seeds the protein is superposed on CA
atoms and the ligand heavy-atom RMSD is taken WITHOUT re-fitting the ligand, so
the number reports where the ligand sits in the pocket frame rather than whether
its internal geometry matches.

INTERPRETATION, fixed before the numbers:
  < 1.0 A   one pose; the rigid-ligand protocol of 91 can be applied directly
  1-2.5 A   same site, uncertain orientation; scoring must be repeated per pose
            and only conclusions stable across poses may be quoted
  > 2.5 A   no consensus pose. The method needs a crystal structure and that is
            the finding -- §46 already showed a wrong pose does not fail loudly,
            it returns a plausible number.

⚠ Seed-stability is NOT correctness. Ten seeds agreeing on the wrong pose would
look identical to ten agreeing on the right one. This measures precision only;
the mandipropamid arm is the only place we have accuracy.
"""
import glob, itertools, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POSES = os.path.join(ROOT, "data", "poses", "out")


def read_cif(path):
    """-> (CA coords by residue index, ligand heavy-atom coords, ligand elements)"""
    ca, lig, el = [], [], []
    cols, rows, inloop = {}, [], False
    for line in open(path):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".")[1]] = len(cols); inloop = True; continue
        if inloop:
            if line.startswith(("#", "loop_", "_")):
                if rows: break
                continue
            f = line.split()
            if len(f) >= len(cols): rows.append(f)
    for f in rows:
        e = f[cols["type_symbol"]]
        if e == "H": continue
        xyz = [float(f[cols[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")]
        chain = f[cols["label_asym_id"]]
        if f[cols["group_PDB"]] == "ATOM" and f[cols["label_atom_id"]] == "CA":
            ca.append(xyz)
        elif f[cols["group_PDB"]] == "HETATM" or chain == "B":
            lig.append(xyz); el.append(e)
    return np.array(ca), np.array(lig), el


def kabsch(P, Q):
    """rotation+translation taking Q onto P"""
    mp, mq = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((Q - mq).T @ (P - mp))
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return R, mp, mq


def main():
    comps = sorted({os.path.basename(d).rsplit("_s", 1)[0]
                    for d in glob.glob(os.path.join(POSES, "*_s*"))})
    print(f"{'compound':<20}{'seeds':>6}{'lig atoms':>10}"
          f"{'median pairwise RMSD':>22}{'max':>8}{'  verdict'}")
    for c in comps:
        got = []
        for s in range(10):
            g = glob.glob(os.path.join(POSES, f"{c}_s{s}", "boltz_results_*",
                                       "predictions", "*", "*.cif"))
            if g:
                got.append(read_cif(sorted(g)[0]))
        if len(got) < 2:
            print(f"{c:<20}{len(got):>6}   too few"); continue
        nel = {tuple(x[2]) for x in got}
        if len(nel) != 1:
            print(f"{c:<20} ligand atom order differs between seeds -- skipped")
            continue
        ref_ca, ref_lig, _ = got[0]
        rms = []
        for i, j in itertools.combinations(range(len(got)), 2):
            cai, ligi, _ = got[i]
            caj, ligj, _ = got[j]
            n = min(len(cai), len(caj))
            R, mp, mq = kabsch(cai[:n], caj[:n])
            lj = (ligj - mq) @ R + mp
            rms.append(float(np.sqrt(((ligi - lj) ** 2).sum(1).mean())))
        med, mx = float(np.median(rms)), float(np.max(rms))
        verdict = ("ONE POSE" if med < 1.0 else
                   "same site, uncertain orientation" if med < 2.5 else
                   "NO CONSENSUS -- needs a crystal structure")
        print(f"{c:<20}{len(got):>6}{len(ref_lig):>10}{med:>22.2f}{mx:>8.2f}  {verdict}")
    print("\n⚠ seed-stability is precision, not correctness: ten seeds agreeing")
    print("  on the wrong pose look exactly like ten agreeing on the right one.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
