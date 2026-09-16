#!/usr/bin/env python
r"""
209_anchor_place.py -- Whitehead-style PHARMACOPHORE-ANCHORED ligand placement.

THE METHOD, from the nitazene paper's Methods (read verbatim): they do NOT dock.
They fix orientation with an SVD superposition in which "the nitro group at R2
[is] superimposed on the target molecule abscisic acid ketone, with the nitro
group acting as the required hydrogen bond acceptor", then run a "fast-grid
search of PYR1 subjected to a poly-glycine shave" for backbone clashes, repack
with the PyRosetta Packer, and keep complexes below -300 REU.

Anchoring an acceptor atom AND its attached heavy atom leaves exactly ONE
rotational degree of freedom, so the whole pose space is a conformer list times a
1-D spin -- enumerable exhaustively, no sampling and no docking.

ABA'S ANCHOR ATOM, IDENTIFIED NOT ASSUMED. 3QN1 chain A A8S has four oxygens and
this script re-derives which is which from bond lengths and connectivity:
  O11/O12  carboxylate  -- salt-bridges K59 NZ at 2.85/3.02 A
  O7       ring hydroxyl
  O10      ring KETONE (its carbon's other two neighbours are both carbon)
O10 is Whitehead's anchor. ⚠ It sits 11.8 A from K59 NZ, i.e. at the DEEP end of
the pocket, opposite the carboxylate.

⚠ TWO ANCHORS ARE RUN, NOT ONE. The depth analysis found small-ligand sensors
mutate MOUTH positions (Spearman(depth, small/large ratio) = +0.498, p = 0.036;
L87 2.45x, L117 1.91x) while large ones mutate the DEEP end (Y120 0.33x). So
anchoring a small ligand's acceptor onto the deep ketone may be exactly wrong.
Both the ketone (O10, theirs) and the carboxylate (O12, mouth-ward) are used and
the results are compared.

⚠ AND AN UNANCHORED ARM. Our own record warns against hard H-bond constraints:
mandipropamid ignores the W385 latch water entirely (4WVO/8EY0), so a protocol
that requires a specific acceptor contact would have rejected a real sensor. The
unanchored arm enumerates the same conformers over a full rotation grid with no
acceptor constraint; if the anchored pose is also the best unanchored pose, the
anchor is free information, and if it is not, the anchor is an assumption.

⚠ POSE SCORING IS DELIBERATELY CRUDE HERE. This step only has to produce a small
set of non-clashing candidate poses; ranking them is 210's job (Rosetta) and
211's (MM-GBSA). 46b measured docking at 8.33 A from the crystal pose, so nothing
in this project justifies trusting a placement score.

Usage:
  python 209_anchor_place.py --ligand anthrone --arm ketone
  python 209_anchor_place.py --all --dry-run
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "anchor_place")
BB = {"N", "CA", "C", "O", "OXT"}
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "F": 1.47, "CL": 1.75,
       "BR": 1.85, "I": 1.98, "P": 1.80}
CLASH = 0.6            # A of allowed vdW overlap
NSPIN = 36             # 10 degree steps about the anchor axis

#: the control panel. anchor_smarts picks the acceptor atom; eugenol is HELD OUT
#: of calibration and is only placed once the method is frozen.
PANEL = {
    "anthrone":       dict(smiles="O=C1c2ccccc2Cc2ccccc21", uM=1.0,
                           anchor="[OX1]=[CX3]", hit=["F108M", "F159M"]),
    "3-phenylphenol": dict(smiles="Oc1cccc(-c2ccccc2)c1", uM=10.0,
                           anchor="[OX2H]", hit=["A89M"]),
    "chloroxylenol":  dict(smiles="Cc1cc(Cl)c(O)c(C)c1", uM=10.0,
                           anchor="[OX2H]", hit=["V83I", "L87M", "N167Q"]),
    "eugenol":        dict(smiles="COc1cc(CC=C)ccc1O", uM=100.0,
                           anchor="[OX2H]", hit=["V164L", "N167V"], holdout=True),
}
#: the 18 sd03 positions, asserted on every parse
POS18 = {59: "LYS", 81: "VAL", 83: "VAL", 87: "LEU", 89: "ALA", 92: "SER",
         94: "GLU", 108: "PHE", 110: "ILE", 117: "LEU", 120: "TYR", 122: "SER",
         141: "GLU", 159: "PHE", 160: "ALA", 163: "VAL", 164: "VAL", 167: "ASN"}


def read_3qn1():
    cols, rows, inl = {}, [], False
    for line in open(os.path.join(ROOT, "data", "3QN1.cif")):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".")[1]] = len(cols); inl = True; continue
        if inl:
            if line.startswith(("#", "loop_", "_")):
                if rows:
                    break
                continue
            x = line.split()
            if len(x) >= len(cols):
                rows.append(x)
    g = lambda r, k: r[cols[k]]                                  # noqa: E731
    xyz = lambda r: np.array([float(g(r, "Cartn_x")), float(g(r, "Cartn_y")),
                              float(g(r, "Cartn_z"))])           # noqa: E731
    prot, lig, seq = [], {}, {}
    for r in rows:
        if g(r, "type_symbol") == "H":
            continue
        if g(r, "label_comp_id") == "A8S" and g(r, "auth_asym_id") == "A":
            lig[g(r, "label_atom_id")] = (g(r, "type_symbol").upper(), xyz(r))
            continue
        if g(r, "group_PDB") != "ATOM" or g(r, "auth_asym_id") != "A":
            continue
        try:
            n = int(g(r, "auth_seq_id"))
        except ValueError:
            continue
        seq[n] = g(r, "label_comp_id")
        prot.append((n, g(r, "label_atom_id"), g(r, "type_symbol").upper(), xyz(r)))
    bad = [f"{n}:{seq.get(n)}!={a}" for n, a in POS18.items() if seq.get(n) != a]
    assert not bad, f"3QN1 identity assertion FAILED: {bad}"
    return prot, lig


def aba_anchors(lig):
    """Re-derive which ABA oxygen is which from geometry. Returns dict of
    (acceptor_xyz, attached_carbon_xyz) keyed by chemical role."""
    O = {k: v for k, v in lig.items() if v[0] == "O"}
    C = {k: v for k, v in lig.items() if v[0] == "C"}
    out = {}
    for on, (_, op) in O.items():
        d = sorted(((np.linalg.norm(op - cp), cn, cp) for cn, (_, cp) in C.items()))
        dist, cn, cp = d[0]
        nO = sum(1 for xn, (xe, xp) in lig.items()
                 if xe == "O" and xn != on and np.linalg.norm(cp - xp) < 1.8)
        if nO >= 1:
            out.setdefault("carboxylate", []).append((on, op, cp))
        elif dist < 1.40:
            out.setdefault("ketone", []).append((on, op, cp))
        else:
            out.setdefault("hydroxyl", []).append((on, op, cp))
    assert "ketone" in out, f"no ketone found in A8S: {list(O)}"
    assert "carboxylate" in out, "no carboxylate found in A8S"
    return out


def rot_to(a, b):
    """rotation matrix taking unit vector a onto unit vector b."""
    v = np.cross(a, b); c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * (1 / (1 + c))


def spin(axis, ang):
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * (K @ K)


def conformers(smiles, anchor_smarts):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    assert m is not None, f"bad SMILES {smiles}"
    patt = Chem.MolFromSmarts(anchor_smarts)
    hits = m.GetSubstructMatches(patt)
    assert hits, f"anchor {anchor_smarts} not found in {smiles}"
    acc = hits[0][0]
    att = [n.GetIdx() for n in m.GetAtomWithIdx(acc).GetNeighbors()
           if n.GetSymbol() != "H"]
    assert att, "acceptor has no heavy neighbour to define the axis"
    mh = Chem.AddHs(m)
    ids = list(AllChem.EmbedMultipleConfs(mh, numConfs=300, randomSeed=42,
                                          pruneRmsThresh=0.5))
    res = AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=800)
    mh = Chem.RemoveHs(mh)
    E = np.array([e for _, e in res]); E = E - E.min()
    p = np.exp(-E / 0.5924); p = p / p.sum()
    keep = [i for i, pi in zip(ids, p) if pi > 0.02]      # populated only
    el = [a.GetSymbol().upper() for a in mh.GetAtoms()]
    return [(el, np.array(mh.GetConformer(i).GetPositions())) for i in keep], \
        acc, att[0], m.GetNumHeavyAtoms()


def place(conf_xyz, acc, att, target_O, target_C, nspin=NSPIN, anchored=True):
    """yield candidate poses. anchored: match acceptor atom AND its axis."""
    out = []
    v_t = target_C - target_O; v_t /= np.linalg.norm(v_t)
    P0 = conf_xyz - conf_xyz[acc]
    if anchored:
        v_l = P0[att] / np.linalg.norm(P0[att])
        R0 = rot_to(v_l, v_t)
        base = P0 @ R0.T
        for k in range(nspin):
            out.append(base @ spin(v_t, 2 * np.pi * k / nspin).T + target_O)
    else:
        rng = np.random.default_rng(0)
        for _ in range(nspin * 8):
            A = rng.normal(size=(3, 3))
            Q, _ = np.linalg.qr(A)
            if np.linalg.det(Q) < 0:
                Q[:, 0] *= -1
            out.append(P0 @ Q.T + target_O)
    return out


def clashes(P, el, prot, polygly=False):
    sel = [(e, x) for n, a, e, x in prot if (not polygly) or a in BB or a == "CB"]
    X = np.array([x for _, x in sel])
    R = np.array([VDW.get(e, 1.7) for e, _ in sel])
    lr = np.array([VDW.get(e, 1.7) for e in el])
    d = np.linalg.norm(X[None, :, :] - P[:, None, :], axis=2)
    ov = (lr[:, None] + R[None, :]) - d
    return int((ov > CLASH).sum()), float(ov.max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ligand", default=None)
    ap.add_argument("--arm", default="ketone",
                    choices=("ketone", "carboxylate", "unanchored"))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    prot, lig = read_3qn1()
    print(f"3QN1 chain A: {len(prot)} heavy protein atoms, "
          f"A8S {len(lig)} atoms; 18/18 sd03 identities asserted")
    anch = aba_anchors(lig)
    for k, v in anch.items():
        print(f"   {k:<12} {[x[0] for x in v]}")
    ligs = list(PANEL) if a.all else [a.ligand]
    arms = ("ketone", "carboxylate", "unanchored") if a.all else (a.arm,)
    summary = []
    for name in ligs:
        spec = PANEL[name]
        confs, acc, att, nheavy = conformers(spec["smiles"], spec["anchor"])
        print(f"\n{name}: {len(confs)} populated conformers, {nheavy} heavy atoms, "
              f"anchor atom idx {acc} (axis to {att})"
              + ("   [HELD OUT]" if spec.get("holdout") else ""))
        for arm in arms:
            if arm == "unanchored":
                tO, tC = anch["ketone"][0][1], anch["ketone"][0][2]
            else:
                tO, tC = anch[arm][0][1], anch[arm][0][2]
            kept = []
            for ci, (el, cx) in enumerate(confs):
                for pi, P in enumerate(place(cx, acc, att, tO, tC,
                                             anchored=(arm != "unanchored"))):
                    npg, _ = clashes(P, el, prot, polygly=True)
                    if npg > 0:
                        continue
                    nfa, worst = clashes(P, el, prot, polygly=False)
                    kept.append(dict(conf=ci, pose=pi, fa_clashes=nfa,
                                     worst_overlap=round(worst, 2),
                                     xyz=P.tolist(), el=el))
            kept.sort(key=lambda r: (r["fa_clashes"], r["worst_overlap"]))
            ntried = len(confs) * (NSPIN if arm != "unanchored" else NSPIN * 8)
            print(f"   {arm:<12} {len(kept):>4} / {ntried:<5} poses clear the "
                  f"poly-Gly shave; best full-atom clashes "
                  f"{kept[0]['fa_clashes'] if kept else '-'}")
            summary.append(dict(ligand=name, arm=arm, uM=spec["uM"],
                                n_conf=len(confs), n_tried=ntried,
                                n_polygly_ok=len(kept),
                                best_fa_clashes=kept[0]["fa_clashes"] if kept else None))
            if kept and not a.dry_run:
                json.dump(kept[:20], open(os.path.join(
                    OUT, f"{name}_{arm}.json"), "w"))
    json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    print(f"\nwrote {OUT}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
