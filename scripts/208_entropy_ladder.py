#!/usr/bin/env python
r"""
208_entropy_ladder.py -- does LIGAND CONFORMATIONAL ENTROPY predict sensor
sensitivity, independent of size and permeability?

THE HYPOTHESIS. Eugenol (min_conc 100 uM, the weakest band) has 3 rotatable
bonds and 5 populated conformers spread over 1.0 kcal/mol, with ALL the freedom
in one allyl tail; the more sensitive small-ligand sensors -- anthrone (1 uM) and
chloroxylenol (10 uM) -- have ZERO rotatable bonds and exactly ONE conformer. A
floppy tail that gains no contacts costs entropy on binding without paying it
back. Jannis: the 4-ligand ladder is not enough data, so this runs all 181.

ENTROPY MEASURE. Not a rotatable-bond count: the Boltzmann conformational
entropy over the MMFF-optimised ETKDG ensemble,
    S = -R * sum(p_i ln p_i),  p_i from relative MMFF energies at 298 K
which counts only conformers that are actually POPULATED. Eugenol has 6
conformers but one is 8.3 kcal/mol up and contributes nothing.

⚠ TWO CONFOUNDS, both measured in this project and both controlled here:
  * SIZE -- bigger ligands have more rotatable bonds, and hit rate is already
    known to be two-sided in heavy-atom count (4.0 % at <=12, 10.0 % at 13-20,
    2.0 % at 28-50).
  * PERMEABILITY -- min_conc is an in vivo yeast readout and tracks cLogP at
    Spearman -0.374, p = 2.1e-7 (README 120a).
Both are removed by rank-partial correlation, and the small-ligand stratum is
reported separately.

PREDICTIONS, FIXED BEFORE THE RUN:
  1. Spearman(conformational entropy, log10 min_conc) > 0 -- more entropy, worse
     sensor -- at p < 0.05.
  2. It SURVIVES partialling out heavy atoms and cLogP. If it does not, the
     ladder was a size or permeability effect wearing a different label.
  3. It holds within the small-ligand stratum (<=16 HA) where eugenol sits.

⚠ min_conc is a screening readout over HITS FOUND (positive-unlabeled), so this
describes what was recovered, not what is achievable.
"""
import os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                          # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
AA = set("ACDEFGHIKLMNPQRSTVWY")
RT = 0.5924            # kcal/mol at 298 K
R_GAS = 0.0019872      # kcal/mol/K


def rank(x):
    from scipy.stats import rankdata
    return rankdata(x)


def partial(x, y, covs):
    """Spearman partial correlation: rank-transform, regress out, correlate."""
    from scipy.stats import pearsonr
    X, Y = rank(x), rank(y)
    C = np.column_stack([rank(c) for c in covs] + [np.ones(len(X))])
    bx = np.linalg.lstsq(C, X, rcond=None)[0]
    by = np.linalg.lstsq(C, Y, rcond=None)[0]
    return pearsonr(X - C @ bx, Y - C @ by)


def main():
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Crippen, rdMolDescriptors as D
    from scipy.stats import spearmanr
    RDLogger.DisableLog("rdApp.*")
    hdr, recs = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    pc = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    best = {}
    for r in recs:
        lg = (r.get("library_name") or "").strip()
        sm = (r.get("canonical_smiles") or "").strip()
        mc = (r.get("min_conc") or "").strip()
        s = [c for c in pc if len((r.get(c) or "").strip()) == 1
             and (r.get(c) or "").strip().upper() in AA
             and (r.get(c) or "").strip().upper() != c[0]]
        if not (lg and sm and s and mc):
            continue
        try:
            v = float(mc)
        except ValueError:
            continue
        if lg not in best or v < best[lg][0]:
            best[lg] = (v, sm)
    print(f"{len(best)} ligands with a best min_conc", flush=True)

    rows = []
    for i, (lg, (v, sm)) in enumerate(sorted(best.items()), 1):
        m = Chem.MolFromSmiles(sm)
        if m is None:
            continue
        mh = Chem.AddHs(m)
        ids = list(AllChem.EmbedMultipleConfs(mh, numConfs=300, randomSeed=42,
                                              pruneRmsThresh=0.5))
        if not ids:
            continue
        try:
            res = AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=800)
        except Exception:                                    # noqa: BLE001
            continue
        E = np.array([e for ok, e in res if ok == 0])
        if len(E) == 0:
            continue
        E = E - E.min()
        p = np.exp(-E / RT); p = p / p.sum()
        S = float(-R_GAS * np.sum(p * np.log(p + 1e-300)))
        rows.append(dict(lig=lg, mc=v, S=S, nconf=len(E),
                         neff=float(np.exp(-np.sum(p * np.log(p + 1e-300)))),
                         rotb=D.CalcNumRotatableBonds(m),
                         ha=m.GetNumHeavyAtoms(), logp=Crippen.MolLogP(m)))
        if i % 40 == 0:
            print(f"  {i}/{len(best)}", flush=True)

    print(f"\n{len(rows)} ligands with a conformer ensemble\n")
    mc = np.log10([r["mc"] for r in rows])
    S = np.array([r["S"] for r in rows])
    ne = np.array([r["neff"] for r in rows])
    rb = np.array([r["rotb"] for r in rows], float)
    ha = np.array([r["ha"] for r in rows], float)
    lp = np.array([r["logp"] for r in rows])
    print("P1  raw correlations with log10(min_conc)  (+ = more of it, WORSE sensor)")
    for nm, x in (("conformational entropy S", S), ("effective n conformers", ne),
                  ("rotatable bonds", rb), ("heavy atoms", ha), ("cLogP", lp)):
        r_, p_ = spearmanr(x, mc)
        print(f"   {nm:<26} rho {r_:+.3f}  p={p_:.3g}"
              + ("  <-- " if p_ < 0.05 else ""))
    print("\nP2  partial, removing heavy atoms AND cLogP")
    for nm, x in (("conformational entropy S", S), ("effective n conformers", ne),
                  ("rotatable bonds", rb)):
        r_, p_ = partial(x, mc, [ha, lp])
        print(f"   {nm:<26} rho {r_:+.3f}  p={p_:.3g}"
              + ("  SURVIVES" if p_ < 0.05 else "  dies"))
    print("\nP3  small-ligand stratum only (<=16 heavy atoms)")
    k = ha <= 16
    print(f"   n = {int(k.sum())}")
    for nm, x in (("conformational entropy S", S), ("rotatable bonds", rb)):
        r_, p_ = spearmanr(x[k], mc[k])
        rp, pp = partial(x[k], mc[k], [ha[k], lp[k]])
        print(f"   {nm:<26} raw rho {r_:+.3f} p={p_:.3g} | partial rho {rp:+.3f} p={pp:.3g}")
    print("\n  the 4-ligand ladder, for reference:")
    for nm in ("Anthrone", "Chloroxylenol", "3-Phenylphenol", "Eugenol"):
        e = [r for r in rows if r["lig"].lower() == nm.lower()]
        if e:
            r0 = e[0]
            print(f"    {nm:<16}{r0['mc']:>6.0f} uM   S {r0['S']:.5f}  "
                  f"n_eff {r0['neff']:.2f}  rotB {r0['rotb']}  HA {r0['ha']}")
    import json
    os.makedirs(os.path.join(ROOT, "results", "entropy_ladder"), exist_ok=True)
    json.dump(rows, open(os.path.join(ROOT, "results", "entropy_ladder",
                                      "rows.json"), "w"), indent=1)
    print(f"\nwrote results/entropy_ladder/rows.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
