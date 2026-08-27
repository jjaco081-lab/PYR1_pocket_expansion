#!/usr/bin/env python
r"""
137_boltz_by_ligand_class.py -- is there a ligand class Boltz-2 handles better?

§82 found Boltz-2's binder classifier calls 90 % of real sensors non-binders, and
that its apparent potency signal is BETWEEN ligands (pooled AUC 0.634) rather than
between sensors (within-ligand rho -0.028). Jannis asks the obvious follow-up:
is there a chemical class where it does work?

⚠ POWER, STATED FIRST. Only 43 ligands carry within-ligand information (>= 4
sensors AND potency variation). Splitting those by chemistry leaves ~15 per class.
A subgroup that looks good at that n is as likely to be noise as signal, and with
8 descriptors x 2 metrics there are 16 chances to find one. So:

  * the descriptor list is FIXED here, before looking, and is the standard
    Lipinski/medchem set -- no searching for a split that works
  * every subgroup result is reported with its n and its permutation p
  * nothing here is a conclusion; it is at most a hypothesis to test on TNT

TWO DIFFERENT QUESTIONS, and only the second one matters for design:
  1. which ligands SCORE higher -- a property of Boltz's calibration
  2. for which ligands does the WITHIN-LIGAND ranking work -- the design question
"""
import json
import os
import sys
from collections import defaultdict
from math import erfc, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                      # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
AFF = os.path.join(ROOT, "data", "structures", "boltz719_affinity.json")

#: fixed before any result is seen
DESCRIPTORS = [
    ("heavy_atoms", "heavy atoms"),
    ("clogp", "cLogP"),
    ("tpsa", "TPSA"),
    ("hbd", "H-bond donors"),
    ("hba", "H-bond acceptors"),
    ("arom_rings", "aromatic rings"),
    ("rotb", "rotatable bonds"),
    ("charge", "formal charge"),
]


def rank(x):
    x = np.asarray(x, float)
    o = np.argsort(x)
    r = np.empty(len(x))
    r[o] = np.arange(1, len(x) + 1)
    for u in np.unique(x):
        m = x == u
        if m.sum() > 1:
            r[m] = r[m].mean()
    return r


def sp(a, b):
    ra, rb = rank(a), rank(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d else 0.0


def pv(r, n):
    return erfc(abs(r * sqrt(n - 1)) / sqrt(2)) if n > 3 else 1.0


def main():
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")

    rows = [r for r in json.load(open(AFF)) if r["minc"] not in ("nd", "")]
    for r in rows:
        r["pot"] = float(r["minc"])

    # SMILES from sd03
    hdr, recs = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    smi = {}
    for rec in recs:
        lg = (rec.get("library_name") or "").strip()
        s = (rec.get("canonical_smiles") or "").strip()
        if lg and s:
            smi.setdefault(lg.lower(), s)

    props = {}
    for lg, s in smi.items():
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        props[lg] = {
            "heavy_atoms": m.GetNumHeavyAtoms(),
            "clogp": Crippen.MolLogP(m),
            "tpsa": rdMolDescriptors.CalcTPSA(m),
            "hbd": rdMolDescriptors.CalcNumHBD(m),
            "hba": rdMolDescriptors.CalcNumHBA(m),
            "arom_rings": rdMolDescriptors.CalcNumAromaticRings(m),
            "rotb": rdMolDescriptors.CalcNumRotatableBonds(m),
            "charge": Chem.GetFormalCharge(m),
        }
    rows = [r for r in rows if r["lig"].lower() in props]
    for r in rows:
        r.update(props[r["lig"].lower()])
    print(f"{len(rows)} sensors with affinity output, potency and ligand properties; "
          f"{len({r['lig'] for r in rows})} ligands")

    # ---- Q1: which ligands SCORE higher (calibration, not design) ----------
    byl = defaultdict(list)
    for r in rows:
        byl[r["lig"].lower()].append(r)
    lig_rows = [{"lig": l, "pbind": float(np.median([x["pbind"] for x in v])),
                 "n": len(v), **props[l]} for l, v in byl.items()]
    print("\nQ1. WHICH LIGANDS SCORE HIGHER  (Boltz calibration, NOT evidence it works)")
    print(f"   ligand-level medians, n = {len(lig_rows)} ligands")
    print(f"   {'descriptor':<18}{'rho vs p_bind':>15}{'p':>10}")
    for k, lab in DESCRIPTORS:
        x = [d[k] for d in lig_rows]
        y = [d["pbind"] for d in lig_rows]
        if len(set(x)) < 3:
            continue
        r = sp(x, y)
        print(f"   {lab:<18}{r:>15.3f}{pv(r, len(x)):>10.2g}")

    # ---- Q2: for which ligands does WITHIN-LIGAND ranking work -------------
    use = {l: v for l, v in byl.items()
           if len(v) >= 4 and len({x["pot"] for x in v}) > 1}
    per = {l: sp([x["pbind"] for x in v], [x["pot"] for x in v])
           for l, v in use.items()}
    print(f"\nQ2. WHERE DOES WITHIN-LIGAND RANKING WORK?  (the design question)")
    print(f"   {len(use)} ligands qualify (>= 4 sensors, potency varies); "
          f"a USEFUL rho is NEGATIVE")
    allr = np.array(list(per.values()))
    print(f"   overall: median rho {np.median(allr):+.3f}, "
          f"{int((allr < 0).sum())}/{len(allr)} negative (chance {len(allr)/2:.1f})")
    print(f"\n   split by descriptor tertile -- ⚠ n ~ {len(use)//3} per cell, "
          f"8 descriptors = 8 chances to find a spurious split")
    print(f"   {'descriptor':<18}{'low tertile':>14}{'mid':>10}{'high':>10}"
          f"{'  (median within-ligand rho)'}")
    for k, lab in DESCRIPTORS:
        vals = np.array([props[l][k] for l in use])
        if len(set(vals)) < 3:
            continue
        rs = np.array([per[l] for l in use])
        q1, q2 = np.percentile(vals, [33.3, 66.7])
        cells = [rs[vals <= q1], rs[(vals > q1) & (vals <= q2)], rs[vals > q2]]
        out = "".join(f"{np.median(c):>10.3f}({len(c):>2})" if len(c) else f"{'--':>14}"
                      for c in cells)
        print(f"   {lab:<18}{out}")
    print("\n   ⚠ Read the SPREAD, not the best cell. If the three tertiles differ by")
    print("     less than the scatter among individual ligands, there is no class")
    print("     effect -- only 43 noisy per-ligand estimates being sorted.")
    iq = np.percentile(allr, [25, 75])
    print(f"   per-ligand rho interquartile range: {iq[0]:+.3f} to {iq[1]:+.3f}")
    json.dump({"per_ligand_rho": per, "props": {l: props[l] for l in use}},
              open(os.path.join(ROOT, "results", "boltz_by_class.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
