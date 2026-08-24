#!/usr/bin/env python
r"""
107_lookup_baseline.py -- the trivial baseline this project never established.

THE QUESTION (§49c)
Before building anything from structure, find out what a lookup table achieves:
for a new ligand, take the chemically nearest ligand that has ALREADY been
screened and reuse its substitution menu. If that recalls as much as a
physics-based method, then every structural result here has been competing with
table lookup, and that has to be known before more effort is spent.

THE NULL THAT ACTUALLY MATTERS
The obvious trap is to compare similarity-lookup against nothing and declare
success. Common substitutions are common: a menu of the globally most frequent
substitutions will recall a lot of any ligand's chemistry without using the
ligand at all. So the comparison here is against a FREQUENCY baseline that is
blind to the target ligand. Similarity only earns its keep if it beats that.

    freq       rank substitutions by how often they appear across OTHER ligands
    sim        rank by sum over other ligands of Tanimoto(target, other) x count
    1nn        take the single most similar ligand's observed substitutions
    oracle     the target's own substitutions (ceiling; not achievable)

EVALUATION
Leave-one-LIGAND-out over sd03, which is the only set carrying SMILES. For each
held-out ligand with enough clones to have a ground truth, rank candidate
substitutions with each scorer, take the top K, and measure recall of the
substitutions actually observed in that ligand's functional clones.

CONFOUND, HANDLED: a clone can only contain substitutions its library allowed, so
a target screened in `tsm` cannot display a `dsm`-only residue. Evaluation is
therefore run per mut_lib as well as pooled, and both are reported. Pooling alone
would let a neighbour from a richer library appear to "predict" residues the
target could never have shown.

HELD OUT: sd08 (TNT) and sd09 (PFAS) are not touched here.
Needs RDKit: conda_envs/esmfold2/bin/python.
"""
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_xlsx import table                                   # noqa: E402
from rdkit import Chem, DataStructs, RDLogger                # noqa: E402
from rdkit.Chem import rdFingerprintGenerator                # noqa: E402

RDLogger.DisableLog("rdApp.*")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
AA = set("ACDEFGHIKLMNPQRSTVWY")
MIN_CLONES = 3            # a ligand needs this many clones to have a ground truth


def load():
    hdr, recs = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    pc = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    per = defaultdict(lambda: {"subs": Counter(), "clones": 0, "smi": None, "lib": Counter()})
    for r in recs:
        lig = (r.get("library_name") or "?").strip()
        smi = (r.get("canonical_smiles") or "").strip()
        lib = (r.get("mut_lib") or "?").strip()
        d = per[lig]
        d["clones"] += 1
        d["lib"][lib] += 1
        if smi and not d["smi"]:
            d["smi"] = smi
        for c in pc:
            v = (r.get(c) or "").strip().upper()
            if len(v) == 1 and v in AA and v != c[0]:
                d["subs"][f"{c[0]}{c[1:]}{v}"] += 1
    return per


def fingerprints(per):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = {}
    for lig, d in per.items():
        if not d["smi"]:
            continue
        m = Chem.MolFromSmiles(d["smi"])
        if m is None:
            continue
        fp[lig] = gen.GetFingerprint(m)
    return fp


def evaluate(per, fp, ligs, Ks=(5, 10, 20, 40)):
    """Leave-one-ligand-out recall@K for each scorer."""
    res = {s: {k: [] for k in Ks} for s in ("freq", "sim", "1nn", "oracle")}
    nn_sim = []
    for tgt in ligs:
        truth = set(per[tgt]["subs"])
        if not truth:
            continue
        others = [l for l in ligs if l != tgt and l in fp]
        if not others:
            continue
        sims = {l: DataStructs.TanimotoSimilarity(fp[tgt], fp[l]) for l in others}
        nn = max(sims, key=sims.get)
        nn_sim.append(sims[nn])

        freq, sim = Counter(), Counter()
        for l in others:
            for s, n in per[l]["subs"].items():
                freq[s] += n
                sim[s] += sims[l] * n
        ranked = {
            "freq": [s for s, _ in freq.most_common()],
            "sim": [s for s, _ in sim.most_common()],
            "1nn": [s for s, _ in per[nn]["subs"].most_common()],
            "oracle": [s for s, _ in per[tgt]["subs"].most_common()],
        }
        for name, order in ranked.items():
            for K in Ks:
                got = set(order[:K]) & truth
                res[name][K].append(len(got) / len(truth))
    return res, nn_sim


def show(title, res, n, nn_sim, Ks=(5, 10, 20, 40)):
    print(f"\n{'='*76}\n{title}  ({n} held-out ligands)\n{'='*76}")
    if nn_sim:
        import statistics as st
        print(f"  nearest-neighbour Tanimoto: median {st.median(nn_sim):.2f}, "
              f"min {min(nn_sim):.2f}, max {max(nn_sim):.2f}")
    print(f"  {'scorer':<10}" + "".join(f"{'recall@'+str(k):>12}" for k in Ks))
    for name in ("freq", "sim", "1nn", "oracle"):
        row = f"  {name:<10}"
        for K in Ks:
            v = res[name][K]
            row += f"{(sum(v)/len(v) if v else 0):>12.3f}"
        print(row)
    d = [a - b for a, b in zip(res["sim"][20], res["freq"][20])]
    if d:
        win = sum(1 for x in d if x > 0)
        print(f"\n  sim - freq at K=20: mean {sum(d)/len(d):+.3f}; "
              f"similarity wins on {win}/{len(d)} ligands "
              f"({100*win/len(d):.0f}%)")


def main():
    per = load()
    fp = fingerprints(per)
    ligs_all = [l for l, d in per.items()
                if d["clones"] >= MIN_CLONES and l in fp and d["subs"]]
    print(f"sd03: {len(per)} ligands, {len(fp)} with usable SMILES, "
          f"{len(ligs_all)} with >= {MIN_CLONES} clones and >=1 substitution")

    res, nn = evaluate(per, fp, ligs_all)
    show("POOLED across mut_lib (see header: this is the OPTIMISTIC reading)",
         res, len(ligs_all), nn)

    # per-library, which is the honest comparison
    for lib in ("tsm", "dsm"):
        sub = [l for l in ligs_all if per[l]["lib"].most_common(1)[0][0] == lib]
        if len(sub) < 8:
            continue
        r2, nn2 = evaluate(per, fp, sub)
        show(f"mut_lib = {lib} only (target and neighbours share allowed residues)",
             r2, len(sub), nn2)

    print("\n  READ: 'sim' only earns its keep where it beats 'freq', which does not")
    print("  look at the ligand at all. 'oracle' is the ceiling, not a method.")


if __name__ == "__main__":
    main()
