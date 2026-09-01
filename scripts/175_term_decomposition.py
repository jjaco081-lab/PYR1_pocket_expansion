#!/usr/bin/env python
r"""
175_term_decomposition.py -- is the missing signal in a different ENERGY TERM?

§99c's gap: we can tell which mutation makes ROOM for a ligand too big for the
pocket (mandipropamid, AUC 0.868) and we cannot tell which mutation makes a
pocket bind a ligand that ALREADY FITS (fludioxonil 0.550, benzothiadiazole
0.528, benoxacor 0.469). §79 says the working score is essentially unrelieved
`fa_rep` -- a clash detector. So the obvious question is whether the information
for the fitting ligands lives in a term the total is drowning out.

ref2015 is a weighted sum, so this needs no new physics: rescore the same 475
variants and record dG_bind PER TERM instead of only the total.

  fa_rep                    steric repulsion -- the one that works, as control
  fa_atr                    dispersion attraction: shape complementarity
  fa_sol, lk_ball_wtd       DESOLVATION -- the cost of burying polar groups
  fa_elec                   electrostatics
  hbond_sc, hbond_bb_sc     H-BOND SATISFACTION, side-chain and backbone-to-side-chain
  n_hbond_lig               explicit count of hydrogen bonds made to the ligand

Each is scored the same way as the total: dG_bind = E(complex) - E(protein) -
E(ligand) for that term alone, paired against a wild-type control repacked in
the identical shell, 3 replicates, minimum taken. Protocol otherwise identical
to 156/166 so the per-term AUCs are comparable with the totals already recorded.

PRE-REGISTERED, so it cannot be read as post-hoc:
  * fa_rep should reproduce the known pattern -- high AUC on mandipropamid, chance
    on the three that fit. If it does not, the decomposition is wrong.
  * a term that carries the fitting-ligand signal must beat chance on
    fludioxonil / benzothiadiazole / benoxacor, where the total does not.
  * 8 terms x 4 compounds = 32 tests. A single nominal p < 0.05 means nothing;
    Holm-corrected significance is what counts, and is reported.
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import lib_rosetta as LR                                            # noqa: E402

AP = os.path.join(ROOT, "data", "agro_params")
OUT = os.path.join(ROOT, "results", "term_decomp")
FRAMES = {"mandipropamid_crystal": (os.path.join(ROOT, "data", "stage1", "wt_mandi.pdb"),
                                    os.path.join(ROOT, "data", "stage1", "params", "3UZ.params"),
                                    os.path.join(ROOT, "data", "stage1", "params", "3UZ_0001.pdb"))}
PARK25 = {55: "P", 61: "F", 62: "I", 81: "V", 83: "V", 87: "L", 88: "P", 89: "A",
          92: "S", 94: "E", 141: "E", 108: "F", 110: "I", 115: "H", 116: "R",
          117: "L", 120: "Y", 122: "S", 158: "M", 159: "F", 160: "A", 162: "T",
          163: "V", 164: "V", 167: "N"}
THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
         "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
         "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
         "W": "TRP", "Y": "TYR"}
AA20 = "ACDEFGHIKLMNPQRSTVWY"
TERMS = ["fa_rep", "fa_atr", "fa_sol", "lk_ball_wtd", "fa_elec",
         "hbond_sc", "hbond_bb_sc"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compound", required=True)
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--nchunks", type=int, default=1)
    ap.add_argument("--nrep", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    if a.compound in FRAMES:
        frame, prm, ligpdb = FRAMES[a.compound]
    else:
        meta = json.load(open(os.path.join(AP, "summary.json")))[a.compound]
        code = meta["code"]; prm = meta["params"]
        ligpdb = os.path.join(AP, f"{code}_0001.pdb")
        frame = os.path.join(AP, f"frame_{code}.pdb")

    import pyrosetta
    from pyrosetta.rosetta.core.scoring import ScoreType
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ex1 -ex2aro "
                   f"-extra_res_fa {prm}")
    sfxn = pyrosetta.create_score_function("ref2015_cart")
    base = pyrosetta.pose_from_pdb(frame)
    info = base.pdb_info()
    idx = {info.number(i): i for i in range(1, base.total_residue() + 1)}
    lig = base.total_residue()
    bad = {n: base.residue(idx[n]).name1() for n, w in PARK25.items()
           if base.residue(idx[n]).name1() != w}
    assert not bad, f"Park identities disagree: {bad}"
    MutateResidue(idx[59], "ARG").apply(base)

    st = {t: getattr(ScoreType, t) for t in TERMS
          if hasattr(ScoreType, t)}
    print(f"{a.compound}: terms {sorted(st)}", flush=True)
    ligpose = pyrosetta.pose_from_pdb(ligpdb)
    sfxn(ligpose)
    e_lig = {t: ligpose.energies().total_energies()[s] for t, s in st.items()}

    def terms_of(p):
        sfxn(p)
        E = {t: p.energies().total_energies()[s] for t, s in st.items()}
        q = p.clone(); q.delete_residue_slow(q.total_residue()); sfxn(q)
        Q = {t: q.energies().total_energies()[s] for t, s in st.items()}
        out = {t: E[t] - Q[t] - e_lig[t] for t in st}
        # explicit H-bonds involving the ligand
        hs = pyrosetta.rosetta.core.scoring.hbonds.HBondSet()
        p.update_residue_neighbors()
        pyrosetta.rosetta.core.scoring.hbonds.fill_hbond_set(p, False, hs)
        out["n_hbond_lig"] = -float(sum(
            1 for i in range(1, hs.nhbonds() + 1)
            if hs.hbond(i).don_res() == lig or hs.hbond(i).acc_res() == lig))
        return out

    def repack(p0, allowed, nrep):
        best = None
        for _ in range(nrep):
            p = p0.clone()
            tf, _ = LR.restrict_packing(p, allowed)
            pk = PackRotamersMover(sfxn); pk.task_factory(tf); pk.apply(p)
            v = terms_of(p)
            if best is None or v["fa_rep"] < best["fa_rep"]:
                best = v
        return best

    jobs = [(n, mu) for n in sorted(PARK25) for mu in AA20 if mu != PARK25[n]]
    jobs = [j for i, j in enumerate(jobs) if i % a.nchunks == a.chunk]
    shell, bgc, rows = {}, {}, []
    for num, mu in jobs:
        if num not in shell:
            sel = ResidueIndexSelector(str(idx[num]))
            nb = NeighborhoodResidueSelector(sel, 6.0, True)
            shell[num] = sorted({i + 1 for i, b in enumerate(nb.apply(base)) if b}
                                | {lig})
            bgc[num] = repack(base.clone(), shell[num], a.nrep)
        p = base.clone()
        MutateResidue(idx[num], THREE[mu]).apply(p)
        v = repack(p, shell[num], a.nrep)
        rows.append({"sub": f"{PARK25[num]}{num}{mu}", "pos": num,
                     **{f"d_{t}": v[t] - bgc[num][t] for t in v}})
        json.dump(rows, open(os.path.join(
            OUT, f"{a.compound}_{a.chunk:03d}.json"), "w"))
    print(f"wrote {len(rows)}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
