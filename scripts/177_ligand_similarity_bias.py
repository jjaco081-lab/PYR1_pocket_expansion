#!/usr/bin/env python
r"""
177_ligand_similarity_bias.py -- an MSA-style ligand prior, as a BIASER.

Jannis's proposal: by analogy with an MSA guiding structure prediction, aggregate
everything known about every screened ligand into a similarity-weighted prior
over substitutions, and use it not as a predictor but to bias the choice among
sequences that other filters already allow.

WHY THIS IS NOT §52 AGAIN. §52 tested chemical-similarity LOOKUP and killed it:
leave-one-ligand-out over 89 ligands beat the ligand-blind frequency menu by
+0.010 recall@20 (23 win / 13 loss / 53 TIE), copying the single nearest ligand
was far worse (0.20 vs 0.35), and only 0.4 % of ligand pairs reach Tanimoto 0.5
(median 0.110). But an MSA does not work by finding one close homolog -- it
aggregates many weak, distant relationships. §52 searched for close neighbours
and correctly found none. GRADED AGGREGATION OVER EVERY LIGAND, WITH NO
THRESHOLD, is what is untested, and it is what Jannis is describing.

There is also precedent for it paying: §68 found ECFP4 class-weighting moves
position recovery to 11/11 in the top 12 and the best residue to rank 1 at 7 of
11 positions (blind: 4/11) -- the one place in this project where knowing the
ligand has helped. It was fragile, which a biaser tolerates better than a
predictor.

DESIGN, fixed before any result.
  prior      P(sub | q) = sum_j sim(q,j)^alpha * freq(sub | ligand j), over ALL
             training ligands, no threshold
  similarity ECFP4 Tanimoto, and a descriptor-space similarity built from the
             same z-scored medchem set §84 matched on -- Tanimoto alone is
             degenerate here (median 0.11) and cannot spread the ligands
  split      LIGAND-level, 50/50, repeated 50x. It must be ligands, not clones:
             §57 established the effective sample size is ~208 ligands against
             ~1,150 clones, because clones sharing a ligand are repeats.
  bar        the LIGAND-BLIND frequency prior, which §52 measured at recall@20 =
             0.35 and recall@40 = 0.52. Beating chance is meaningless; this is
             the number to beat, and it is advantaged.
  ⚠ alpha sensitivity is reported, not tuned away: §68's coverage objective went
    from 9/11 to 2/11 on an exponent change from 1 to 3.

⚠ CEILING, stated first. §87b showed the global prior saturates at the 67 %
VOCABULARY ceiling: a third of Beltran's sensors use positions Tian never varied.
Re-weighting can only reorder within the vocabulary, never add to it, so no
version of this can pass that ceiling.
"""
import json, os, re, sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                          # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "ligand_prior")
AA = set("ACDEFGHIKLMNPQRSTVWY")
NSPLIT = 50


def load_clones():
    """ligand -> list of substitution-sets, from sd03 round-1 clones"""
    hdr, recs = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    poscols = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    by, smi = defaultdict(list), {}
    for r in recs:
        lg = (r.get("library_name") or "").strip()
        if not lg:
            continue
        subs = [c + (r.get(c) or "").strip().upper() for c in poscols
                if len((r.get(c) or "").strip()) == 1
                and (r.get(c) or "").strip().upper() in AA
                and (r.get(c) or "").strip().upper() != c[0]]
        if subs:
            by[lg].append(subs)
        s = (r.get("canonical_smiles") or "").strip()
        if s:
            smi.setdefault(lg, s)
    return {k: v for k, v in by.items() if k in smi}, smi


def recall_at(ranked, truth, n):
    return len(set(ranked[:n]) & truth) / max(len(truth), 1)


def main():
    os.makedirs(OUT, exist_ok=True)
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")

    by, smi = load_clones()
    ligs = sorted(by)
    mols = {l: Chem.MolFromSmiles(smi[l]) for l in ligs}
    ligs = [l for l in ligs if mols[l] is not None]
    print(f"{len(ligs)} ligands with clones and parseable SMILES; "
          f"{sum(len(by[l]) for l in ligs)} clones")

    fp = {l: AllChem.GetMorganFingerprintAsBitVect(mols[l], 2, 2048) for l in ligs}
    D = np.array([[m.GetNumHeavyAtoms(), Crippen.MolLogP(m),
                   rdMolDescriptors.CalcTPSA(m), rdMolDescriptors.CalcNumHBD(m),
                   rdMolDescriptors.CalcNumHBA(m),
                   rdMolDescriptors.CalcNumRotatableBonds(m),
                   rdMolDescriptors.CalcNumAromaticRings(m),
                   float(Chem.GetFormalCharge(m))]
                  for m in (mols[l] for l in ligs)], float)
    D = (D - D.mean(0)) / (D.std(0) + 1e-9)
    n = len(ligs)
    T = np.eye(n)
    for i in range(n):
        sims = DataStructs.BulkTanimotoSimilarity(fp[ligs[i]], [fp[l] for l in ligs])
        T[i] = sims
    dist = np.linalg.norm(D[:, None] - D[None], axis=-1)
    Dsim = 1.0 / (1.0 + dist / np.median(dist[dist > 0]))
    print(f"  ECFP4 Tanimoto: median off-diagonal {np.median(T[~np.eye(n,dtype=bool)]):.3f}; "
          f"descriptor sim median {np.median(Dsim[~np.eye(n,dtype=bool)]):.3f}")

    subs_of = {l: Counter(s for c in by[l] for s in c) for l in ligs}
    truth = {l: set(subs_of[l]) for l in ligs}
    idx = {l: i for i, l in enumerate(ligs)}
    rng = np.random.default_rng(0)

    SIMS = {"ECFP4": T, "descriptor": Dsim, "ECFP4 x descriptor": T * Dsim}
    ALPHAS = [1, 2, 4]
    res = defaultdict(list)
    for rep in range(NSPLIT):
        perm = rng.permutation(n)
        tr, te = set(perm[: n // 2]), set(perm[n // 2:])
        blind = Counter()
        for j in tr:
            blind.update(subs_of[ligs[j]])
        rb = [s for s, _ in blind.most_common()]
        for l in (ligs[j] for j in te):
            if not truth[l]:
                continue
            for nn in (20, 40):
                res[("ligand-blind freq", None, nn)].append(recall_at(rb, truth[l], nn))
            for sname, S in SIMS.items():
                for a in ALPHAS:
                    w = defaultdict(float)
                    for j in tr:
                        sim = S[idx[l], j] ** a
                        if sim <= 0:
                            continue
                        for s, c in subs_of[ligs[j]].items():
                            w[s] += sim * c
                    rw = [s for s, _ in sorted(w.items(), key=lambda kv: -kv[1])]
                    for nn in (20, 40):
                        res[(sname, a, nn)].append(recall_at(rw, truth[l], nn))

    print(f"\nheld-out recall over {NSPLIT} random 50/50 LIGAND splits "
          f"(~{n//2} held-out ligands each)")
    print(f"{'prior':<22}{'alpha':>6}{'recall@20':>22}{'recall@40':>22}")
    base = {}
    for nn in (20, 40):
        v = np.array(res[("ligand-blind freq", None, nn)])
        base[nn] = v.mean()
    v20 = np.array(res[("ligand-blind freq", None, 20)])
    v40 = np.array(res[("ligand-blind freq", None, 40)])
    print(f"{'ligand-blind freq':<22}{'-':>6}"
          f"{v20.mean():>14.3f} +-{v20.std()/np.sqrt(len(v20)):<6.3f}"
          f"{v40.mean():>14.3f} +-{v40.std()/np.sqrt(len(v40)):<6.3f}   <- the bar")
    for sname in SIMS:
        for a in ALPHAS:
            x20 = np.array(res[(sname, a, 20)]); x40 = np.array(res[(sname, a, 40)])
            d20 = x20.mean() - base[20]; d40 = x40.mean() - base[40]
            print(f"{sname:<22}{a:>6}"
                  f"{x20.mean():>14.3f} ({d20:+.3f})"
                  f"{x40.mean():>14.3f} ({d40:+.3f})")
    print("\n  (delta) is against the ligand-blind bar on the SAME splits.")
    print(f"  §52 measured the blind prior at recall@20 = 0.35, @40 = 0.52.")
    json.dump({str(k): float(np.mean(v)) for k, v in res.items()},
              open(os.path.join(OUT, "similarity_bias.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
