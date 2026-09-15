#!/usr/bin/env python
r"""
204_latch_advance.py -- how far can the LATCH BACKBONE advance toward a ligand,
and does Leu117's bulk set the floor?

THE MECHANISM, measured in 3QN1 chain A (identity asserted 18/18, auth numbering):
  L117 side chain -> ABA   3.80 A     and it points TOWARD the ligand
  L117 backbone   -> ABA   5.72 A
So the Leu side chain bridges the gap between the latch backbone and the ligand.
ABA (249 A^3 vdW) fills that space. Eugenol (162 A^3) does not, so for a small
ligand the latch must come further in to staple -- and Leu117 is the thing in the
way. That predicts SHRINK at 117 for small ligands, which is what the sd03
responsive clones do (dVol -45.6 vs +14.6, p = 0.0085; L117->N x4, D x3, G x2, A,
versus H/W/M in the weak clones).

⚠ WHY NOT LIGAND-ONLY MD. The panel is rigid: median 2 rotatable bonds and 6
distinct conformers at 0.5 A pruning over 53 small ligands, and ANTHRONE -- the
most responsive of all at 1 uM -- has ZERO rotatable bonds and exactly ONE
conformer. RDKit ETKDG enumerates that space exhaustively in seconds, so MD adds
no conformational information, and with no protein in the box it cannot measure
the latch advance that is actually unknown.

WHAT THIS COMPUTES. For each ligand, every ETKDG conformer is placed in the
closed pocket by maximising overlap with the volume ABA occupies (deterministic:
principal-axis alignment plus a fixed rotational grid, best non-clashing pose
kept). Then, holding the pose fixed, residues 115-117 are translated as a rigid
body along the ligand-centroid direction until any heavy atom clashes, under two
117 identities:

  advance_WT   = permitted translation with Leu117 present
  advance_GLY  = permitted translation with 117 truncated to CB (a Gly/Ala proxy)
  gain         = advance_GLY - advance_WT

PREDICTIONS, FIXED BEFORE THE RUN:
  1. SANITY. gain is ~0 for ABA itself -- ABA fills the space, so removing
     Leu117 buys no advance. If ABA shows a large gain the placement is wrong.
  2. THE TEST. gain is LARGER for small ligands than for large ones.
  3. THE TEST. among small ligands, gain correlates with responsiveness
     (lower min_conc), Spearman p < 0.05.

⚠ The placement is the weak link and is stated, not hidden: it is a shape-overlap
placement against ABA's volume, NOT a docked pose, and it assumes small ligands
sit where ABA sits. Prediction 1 is the control that the placement is not simply
manufacturing gaps.
"""
import json, os, re, sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                          # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "latch_advance")
AA = set("ACDEFGHIKLMNPQRSTVWY")
BB = {"N", "CA", "C", "O", "OXT"}
LATCH = [115, 116, 117]
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "F": 1.47, "CL": 1.75,
       "BR": 1.85, "I": 1.98, "P": 1.80, "NA": 2.27, "H": 1.09}
CLASH = 0.5             # A of allowed vdW overlap before calling it a clash
STEP = 0.1              # A per translation step
MAXADV = 4.0
POS18 = {59: "LYS", 81: "VAL", 83: "VAL", 87: "LEU", 89: "ALA", 92: "SER",
         94: "GLU", 108: "PHE", 110: "ILE", 117: "LEU", 120: "TYR", 122: "SER",
         141: "GLU", 159: "PHE", 160: "ALA", 163: "VAL", 164: "VAL", 167: "ASN"}


def load_3qn1():
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
    res, lig = defaultdict(list), []
    for r in rows:
        if g(r, "type_symbol") == "H":
            continue
        if g(r, "label_comp_id") == "A8S":
            lig.append((g(r, "type_symbol").upper(), xyz(r))); continue
        if g(r, "group_PDB") != "ATOM" or g(r, "auth_asym_id") != "A":
            continue
        try:
            n = int(g(r, "auth_seq_id"))
        except ValueError:
            continue
        res[n].append((g(r, "label_atom_id"), g(r, "label_comp_id"),
                       g(r, "type_symbol").upper(), xyz(r)))
    # ASSERT identity before using any number
    bad = [f"{n}:{res[n][0][1]}!={aa}" for n, aa in POS18.items()
           if n in res and res[n][0][1] != aa]
    assert not bad, f"3QN1 identity assertion failed: {bad}"
    return res, lig


def conformers(smiles, cap=60):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    mh = Chem.AddHs(m)
    ids = AllChem.EmbedMultipleConfs(mh, numConfs=300, randomSeed=42,
                                     pruneRmsThresh=0.5)
    if not len(ids):
        return None
    AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=300)
    mh = Chem.RemoveHs(mh)
    out = []
    el = [a.GetSymbol().upper() for a in mh.GetAtoms()]
    for i in list(ids)[:cap]:
        out.append((el, np.array(mh.GetConformer(i).GetPositions())))
    return out


def principal(P):
    C = P - P.mean(0)
    _, _, V = np.linalg.svd(C, full_matrices=False)
    return C, V


def place(conf_el, conf_xyz, abaP, prot_xyz, prot_r):
    """Best non-clashing pose maximising overlap with ABA's volume.
    Deterministic: principal-axis alignment x a fixed 4x4x4 sign/permutation grid."""
    cA, VA = principal(abaP)
    cen = abaP.mean(0)
    cL, VL = principal(conf_xyz)
    lr = np.array([VDW.get(e, 1.7) for e in conf_el])
    best = None
    for perm in ((0, 1, 2), (1, 0, 2), (0, 2, 1)):
        for s0 in (1, -1):
            for s1 in (1, -1):
                S = np.diag([s0, s1, s0 * s1])
                R = VA.T @ S @ VL[list(perm), :]
                P = cL @ R.T + cen
                d = np.linalg.norm(prot_xyz[None, :, :] - P[:, None, :], axis=2)
                over = (lr[:, None] + prot_r[None, :]) - d
                if over.max() > CLASH:
                    continue
                # overlap with ABA's volume = mean proximity to ABA atoms
                sc = -np.linalg.norm(abaP[None, :, :] - P[:, None, :],
                                     axis=2).min(1).mean()
                if best is None or sc > best[0]:
                    best = (sc, P, lr)
    return best


def advance(latch_atoms, lig_xyz, lig_r, other_xyz, other_r, direction):
    """max translation of the latch along `direction` before a clash."""
    L = np.array([a[3] for a in latch_atoms])
    lr = np.array([VDW.get(a[2], 1.7) for a in latch_atoms])
    tgt_xyz = np.vstack([lig_xyz, other_xyz])
    tgt_r = np.concatenate([lig_r, other_r])
    adv = 0.0
    while adv < MAXADV:
        Q = L + direction * (adv + STEP)
        d = np.linalg.norm(tgt_xyz[None, :, :] - Q[:, None, :], axis=2)
        if ((lr[:, None] + tgt_r[None, :]) - d).max() > CLASH:
            break
        adv += STEP
    return round(adv, 2)


def main():
    os.makedirs(OUT, exist_ok=True)
    res, lig = load_3qn1()
    print("3QN1 identity assertion passed (18/18 sd03 positions)")
    abaP = np.array([x for _, x in lig])
    abaR = np.array([VDW.get(e, 1.7) for e, _ in lig])

    # protein atoms, excluding the latch (it moves) -- used for clash during placement
    prot = [(a[2], a[3]) for n, v in res.items() for a in v if n not in LATCH]
    prot_xyz = np.array([x for _, x in prot])
    prot_r = np.array([VDW.get(e, 1.7) for e, _ in prot])
    latch_full = [a for n in LATCH for a in res[n]]
    latch_trunc = [a for n in LATCH for a in res[n]
                   if n != 117 or a[0] in BB | {"CB"}]
    other = [(a[2], a[3]) for n, v in res.items() for a in v if n not in LATCH]
    other_xyz = np.array([x for _, x in other])
    other_r = np.array([VDW.get(e, 1.7) for e, _ in other])
    # closure direction: latch centroid -> ligand centroid
    lc = np.array([a[3] for a in latch_full]).mean(0)
    direction = abaP.mean(0) - lc
    direction /= np.linalg.norm(direction)
    print(f"closure direction (latch centroid -> ABA centroid): "
          f"{direction.round(2)}")

    # ---- ligand panel -------------------------------------------------------
    hdr, recs = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    pc = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    best = {}
    for r in recs:
        lg = (r.get("library_name") or "").strip()
        sm = (r.get("canonical_smiles") or "").strip()
        mc = (r.get("min_conc") or "").strip()
        s = [c + (r.get(c) or "").strip().upper() for c in pc
             if len((r.get(c) or "").strip()) == 1
             and (r.get(c) or "").strip().upper() in AA
             and (r.get(c) or "").strip().upper() != c[0]]
        if not (lg and sm and s and mc):
            continue
        try:
            v = float(mc)
        except ValueError:
            continue
        m = Chem.MolFromSmiles(sm)
        if not m:
            continue
        if lg not in best or v < best[lg][0]:
            best[lg] = (v, sm, m.GetNumHeavyAtoms(), s)
    # ABA itself is prediction 1's control
    panel = [("ABA (control)", 0.0, "CC1=CC(=O)CC(C)(C)C1/C=C/C(C)=C/C(=O)O", 18, [])]
    panel += [(lg, v, sm, ha, s) for lg, (v, sm, ha, s) in best.items()]
    print(f"\npanel: {len(panel)-1} ligands with min_conc + hit sequence, "
          f"plus ABA as the control\n")

    rows = []
    print(f"{'ligand':<30}{'uM':>6}{'HA':>4}{'nconf':>6}"
          f"{'advWT':>7}{'advGLY':>8}{'gain':>7}")
    for lg, v, sm, ha, s in panel:
        cs = conformers(sm)
        if not cs:
            continue
        pl = None
        for el, cx in cs:
            b = place(el, cx, abaP, prot_xyz, prot_r)
            if b and (pl is None or b[0] > pl[0]):
                pl = b
        if pl is None:
            rows.append(dict(lig=lg, mc=v, ha=ha, nconf=len(cs),
                             advWT=None, advGLY=None, gain=None,
                             note="no non-clashing placement"))
            print(f"{lg[:29]:<30}{v:>6.0f}{ha:>4}{len(cs):>6}"
                  f"{'-':>7}{'-':>8}{'-':>7}  no placement")
            continue
        _, P, lr = pl
        aw = advance(latch_full, P, lr, other_xyz, other_r, direction)
        ag = advance(latch_trunc, P, lr, other_xyz, other_r, direction)
        rows.append(dict(lig=lg, mc=v, ha=ha, nconf=len(cs),
                         advWT=aw, advGLY=ag, gain=round(ag - aw, 2),
                         has117=any(x.startswith("L117") for x in s)))
        print(f"{lg[:29]:<30}{v:>6.0f}{ha:>4}{len(cs):>6}"
              f"{aw:>7.2f}{ag:>8.2f}{ag-aw:>+7.2f}")

    json.dump(rows, open(os.path.join(OUT, "latch_advance.json"), "w"), indent=1)
    ok = [r for r in rows if r.get("gain") is not None and r["lig"] != "ABA (control)"]
    aba = [r for r in rows if r["lig"] == "ABA (control)"]
    from scipy.stats import spearmanr, mannwhitneyu
    print("\n" + "=" * 62)
    if aba:
        print(f"P1 control: ABA gain {aba[0]['gain']:+.2f} A   "
              f"{'HOLDS (ABA fills the space)' if abs(aba[0]['gain']) < 0.6 else 'FAILS -- placement is manufacturing a gap'}")
    sm_ = [r for r in ok if r["ha"] <= 16]
    lg_ = [r for r in ok if r["ha"] >= 17]
    if sm_ and lg_:
        p = mannwhitneyu([r["gain"] for r in sm_], [r["gain"] for r in lg_]).pvalue
        print(f"P2 gain small {np.mean([r['gain'] for r in sm_]):+.2f} A (n={len(sm_)}) "
              f"vs large {np.mean([r['gain'] for r in lg_]):+.2f} (n={len(lg_)})  "
              f"MWU p={p:.4g}")
    if len(sm_) > 5:
        r_, p_ = spearmanr([np.log10(r["mc"]) for r in sm_], [r["gain"] for r in sm_])
        print(f"P3 small ligands: Spearman(log10 min_conc, gain) = {r_:+.3f}, "
              f"p={p_:.4g}   {'HOLDS' if p_ < 0.05 and r_ < 0 else 'null'}")
    print(f"\nwrote {OUT}/latch_advance.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
