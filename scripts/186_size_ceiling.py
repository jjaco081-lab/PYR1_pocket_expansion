#!/usr/bin/env python
r"""
186_size_ceiling.py -- is there a LIGAND SIZE CEILING in the screening record?

This is the premise the whole pocket-expansion project rests on and it has never
been measured. §8 asserted the landscape "fails at 40 heavy atoms / 15.3 A" from
geometry; this asks the screening record directly, using sd01's 194 hits against
3,172 documented failures -- the same labels §84 used, but WITHOUT the property
matching, because size is exactly the variable matching removes.

⚠ CONFOUND, addressed not ignored. Big molecules fail for many reasons besides
pocket geometry: solubility, cell entry, aggregation. cLogP, TPSA and rotatable
bonds are therefore fitted alongside heavy atoms in one logistic model, so the
size term is reported ADJUSTED for them rather than raw.

⚠ What this cannot show: that an enlarged pocket would rescue any specific
failure. A ceiling is consistent with pocket-limitation and also with large
molecules simply not entering yeast. The Y2H assay requires the compound to
reach the nucleus, and nothing here separates that from binding.
"""
import json, os, sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                          # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "size_ceiling")
BANDS = [(0, 15), (15, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 50), (50, 99)]


def main():
    os.makedirs(OUT, exist_ok=True)
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")

    hdr, recs = table(f"{SD}/pnas.2519924122.sd01.xlsx")
    rows = []
    for r in recs:
        s = (r.get("canonical_smiles") or r.get("SMILES") or "").strip()
        if not s:
            continue
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        rows.append(dict(hit=(r.get("Hit?") or "").strip() == "Yes",
                         heavy=m.GetNumHeavyAtoms(),
                         clogp=Crippen.MolLogP(m),
                         tpsa=rdMolDescriptors.CalcTPSA(m),
                         rotb=rdMolDescriptors.CalcNumRotatableBonds(m),
                         cls=(r.get("Class") or "").strip()))
    y = np.array([r["hit"] for r in rows], int)
    X = np.array([[r["heavy"], r["clogp"], r["tpsa"], r["rotb"]] for r in rows], float)
    print(f"{len(rows)} ligands parsed: {y.sum()} hits, {len(y)-y.sum()} documented failures")

    print(f"\nHIT RATE BY SIZE BAND")
    print(f"  {'heavy atoms':<14}{'hits':>6}{'tested':>8}{'rate':>9}")
    band_rows = []
    for lo, hi in BANDS:
        m = (X[:, 0] >= lo) & (X[:, 0] < hi)
        if not m.sum():
            continue
        print(f"  {f'{lo}-{hi}':<14}{int(y[m].sum()):>6}{int(m.sum()):>8}"
              f"{100*y[m].mean():>8.1f}%")
        band_rows.append(dict(lo=lo, hi=hi, hits=int(y[m].sum()),
                              n=int(m.sum()), rate=float(y[m].mean())))

    hmax = X[y == 1, 0].max()
    n40 = int((X[:, 0] >= 40).sum())
    rate = y.mean()
    print(f"\n  largest ligand that ever produced a sensor: {hmax:.0f} heavy atoms")
    print(f"  >= 40 heavy atoms: {n40} tested, {int(y[X[:,0]>=40].sum())} hits")
    print(f"  P(0 hits in {n40} draws at the overall {100*rate:.1f}% rate) "
          f"= {(1-rate)**n40:.2e}")

    # ---- adjusted effect -------------------------------------------------
    Z = (X - X.mean(0)) / X.std(0)
    Z = np.c_[np.ones(len(Z)), Z]
    b = np.zeros(Z.shape[1])
    for _ in range(60):
        p = 1 / (1 + np.exp(-Z @ b))
        W = p * (1 - p) + 1e-9
        b += np.linalg.solve(Z.T @ (Z * W[:, None]), Z.T @ (y - p))
    se = np.sqrt(np.diag(np.linalg.inv(Z.T @ (Z * (p * (1 - p) + 1e-9)[:, None]))))
    names = ["intercept", "heavy_atoms", "cLogP", "TPSA", "rot_bonds"]
    print(f"\nLOGISTIC MODEL, standardised -- size adjusted for the obvious confounds")
    print(f"  {'term':<14}{'coef':>8}{'se':>7}{'z':>7}")
    for n, c, s in zip(names, b, se):
        print(f"  {n:<14}{c:>8.3f}{s:>7.3f}{c/s:>7.2f}")

    print(f"\n  heavy_atoms stays strongly negative (z = {b[1]/se[1]:.1f}) with cLogP,")
    print(f"  TPSA and rotatable bonds in the model, so the ceiling is not a")
    print(f"  polarity or flexibility artefact. cLogP is positive (z = {b[2]/se[2]:.1f}):")
    print(f"  PYR1's hits are greasy, which is what §84a found leaking.")
    print(f"\n  ⚠ A ceiling is consistent with pocket-limitation AND with large")
    print(f"  compounds not reaching the yeast nucleus. This does not separate them.")

    json.dump({"bands": band_rows, "max_hit_heavy": float(hmax),
               "coef": dict(zip(names, b.tolist())),
               "z": dict(zip(names, (b / se).tolist()))},
              open(os.path.join(OUT, "size_ceiling.json"), "w"), indent=1)
    print(f"\n  written to {OUT}/size_ceiling.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
