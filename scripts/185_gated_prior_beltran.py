#!/usr/bin/env python
r"""
185_gated_prior_beltran.py -- the gated similarity prior on an UNTOUCHED set.

§103b left one thing unconfirmed. The whole-molecule similarity prior beats the
ligand-blind bar by +0.05 recall@20, but §103a's control showed the entire gain
comes from ligands WITH a close analog (+0.191 at Tanimoto >= 0.5) and it is
NEGATIVE where none exists (-0.036 at < 0.2). Gating on max Tanimoto >= 0.25
keeps +0.056 across 69 % of ligands.

⚠ That threshold was chosen AFTER seeing the stratification, so it is fitted to
the sd03 data and the +0.056 is optimistic. This runs the FROZEN rule -- prior if
max Tanimoto >= 0.25, ligand-blind otherwise -- on Beltran-45, which was never
used to build or tune anything (§87: 20 design positions, only 10 overlapping
Tian's 18).

NOTHING IS FITTED HERE. The training set is all 194 sd03 ligands, the threshold
is 0.25 as recorded in §103b, the weighting exponent is 2 as recorded, and the
test set is Beltran's 14 cannabinoid ligands. One shot.

BARS, both pre-registered:
  ligand-blind frequency over the same 194 training ligands -- the §52/§103 bar
  the UNGATED prior -- so the gate has to earn its place, not just beat blind

⚠ §87b's ceiling still applies: Tian's whole vocabulary captures only 30/45 of
Beltran's sensors, because 15 use positions Tian never varied. Re-weighting
reorders within the vocabulary and cannot exceed that.
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
GATE = 0.25          # frozen in §103b
ALPHA = 2            # frozen in §103b

#: Beltran's 14 cannabinoid ligands. SMILES from PubChem canonical entries; the
#: names are exactly as they appear in data/beltran/win_sensors.json so the join
#: is asserted rather than assumed.
BEL_SMILES = {
    "WIN 55,212": "CN1CC(COc2cccc3c2c(C(=O)c2ccc4ccccc4c2)cn3CCN2CCOCC2)O1",
    "JWH-015": "CCCn1cc(C(=O)c2ccc3ccccc3c2)c2ccccc21",
    "JWH-072": "CCCn1cc(C(=O)c2ccccc2)c2ccccc21",
    "JWH-016": "CCCCn1cc(C(=O)c2ccc3ccccc3c2)c2ccccc21",
    "JWH-167": "CCCn1cc(C(=O)Cc2ccccc2)c2ccccc21",
    "JWH-018": "CCCCCn1cc(C(=O)c2ccc3ccccc3c2)c2ccccc21",
    "JWH-007": "CCCCCn1c(C)c(C(=O)c2ccc3ccccc3c2)c2ccccc21",
    "JWH-030": "CCCCCn1cc(C(=O)c2cccc3ccccc23)c2ccccc21",
    "JWH-193": "CCCCCn1cc(C(=O)C2(O)CCCCC2)c2ccccc21",
    "CBDA": "CCCCCc1cc(O)c(C(=O)O)c(O)c1C1CC(=C)CCC1C(=C)C",
    "AB-PINACA": "CCCCCn1nc(C(=O)NC(C(C)C)C(N)=O)c2ccccc21",
    "4F-MDMB-BUTINACA": "COC(=O)C(NC(=O)c1nn(CCCCF)c2ccccc12)C(C)(C)C",
    "(±)-CP 47,497": "CCCCCCC(C)(C)c1ccc(O)c(C2CCCC(O)C2)c1",
    "∆9-THC": "CCCCCc1cc(O)c2c(c1)OC(C)(C)C1CCC(C)=CC21",
}


def main():
    os.makedirs(OUT, exist_ok=True)
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")

    # ---- training: all 194 sd03 ligands -------------------------------
    hdr, recs = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    pc = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    clone, smi = defaultdict(list), {}
    for r in recs:
        lg = (r.get("library_name") or "").strip()
        if not lg:
            continue
        s = [c + (r.get(c) or "").strip().upper() for c in pc
             if len((r.get(c) or "").strip()) == 1
             and (r.get(c) or "").strip().upper() in AA
             and (r.get(c) or "").strip().upper() != c[0]]
        if s:
            clone[lg].append(s)
        sm = (r.get("canonical_smiles") or "").strip()
        if sm:
            smi.setdefault(lg, sm)
    train = [l for l in sorted(clone) if l in smi and Chem.MolFromSmiles(smi[l])]
    sub = {l: Counter(s for c in clone[l] for s in c) for l in train}
    fp = {l: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi[l]), 2, 2048)
          for l in train}
    blind = Counter()
    for l in train:
        blind.update(sub[l])
    rb = [s for s, _ in blind.most_common()]
    print(f"training: {len(train)} sd03 ligands, {len(blind)} distinct substitutions")

    # ---- test: Beltran-45, asserted join ------------------------------
    bel = json.load(open(os.path.join(ROOT, "data", "beltran", "win_sensors.json")))
    sensors = [s for s in bel["sensors"] if s["muts"]]
    bel_ligs = sorted({s["ligand"] for s in sensors})
    missing = [l for l in bel_ligs if l not in BEL_SMILES]
    assert not missing, f"no SMILES for Beltran ligands: {missing}"
    truth = defaultdict(set)
    for s in sensors:
        truth[s["ligand"]] |= set(s["muts"])
    print(f"test: {len(bel_ligs)} Beltran ligands, "
          f"{sum(len(v) for v in truth.values())} substitution-ligand pairs")

    def recall(rank, t, n):
        return len(set(rank[:n]) & t) / max(len(t), 1)

    rows = []
    for lg in bel_ligs:
        m = Chem.MolFromSmiles(BEL_SMILES[lg])
        assert m is not None, f"unparseable SMILES for {lg}"
        q = AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)
        sims = np.array(DataStructs.BulkTanimotoSimilarity(q, [fp[l] for l in train]))
        mx = float(sims.max())
        w = defaultdict(float)
        for j, l in enumerate(train):
            s = sims[j] ** ALPHA
            if s > 0:
                for x, c in sub[l].items():
                    w[x] += s * c
        rw = [x for x, _ in sorted(w.items(), key=lambda kv: -kv[1])]
        gated = rw if mx >= GATE else rb
        rows.append(dict(lig=lg, maxT=mx, n_truth=len(truth[lg]),
                         blind20=recall(rb, truth[lg], 20),
                         ungated20=recall(rw, truth[lg], 20),
                         gated20=recall(gated, truth[lg], 20),
                         blind40=recall(rb, truth[lg], 40),
                         ungated40=recall(rw, truth[lg], 40),
                         gated40=recall(gated, truth[lg], 40),
                         used_prior=mx >= GATE))
    print(f"\n{'ligand':<20}{'maxT':>6}{'n':>4}{'blind@20':>10}{'ungated':>9}"
          f"{'gated':>8}   prior used?")
    for r in sorted(rows, key=lambda x: -x["maxT"]):
        print(f"{r['lig']:<20}{r['maxT']:>6.3f}{r['n_truth']:>4}"
              f"{r['blind20']:>10.2f}{r['ungated20']:>9.2f}{r['gated20']:>8.2f}"
              f"   {'yes' if r['used_prior'] else 'no (blind)'}")
    print(f"\n{'':<20}{'':>6}{'':>4}{'recall@20':>10}{'':>9}{'':>8}")
    for nm, k in (("ligand-blind (bar)", "blind"), ("ungated prior", "ungated"),
                  ("GATED prior", "gated")):
        a = np.mean([r[k + "20"] for r in rows]); b = np.mean([r[k + "40"] for r in rows])
        d20 = a - np.mean([r["blind20"] for r in rows])
        print(f"  {nm:<22}recall@20 {a:.3f} ({d20:+.3f})   recall@40 {b:.3f}")
    nu = sum(r["used_prior"] for r in rows)
    print(f"\n  gate fired for {nu}/{len(rows)} Beltran ligands at max Tanimoto >= {GATE}")
    print(f"  §103b predicted +0.056 recall@20 from the gate on sd03; the number")
    print(f"  above is the held-out version, with nothing fitted.")
    json.dump(rows, open(os.path.join(OUT, "gated_beltran.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
