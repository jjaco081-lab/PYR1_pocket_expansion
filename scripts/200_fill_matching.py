#!/usr/bin/env python
r"""
200_fill_matching.py -- OPTION A. Does the volume DEFICIT a ligand leaves in the
pocket predict which substitutions its sensors carry?

THE IDEA. ABA, the cognate ligand, has a vdW volume of 248.8 A^3 (RDKit MMFF,
this script). Eugenol has 162 A^3 -- it under-fills ABA's own footprint by ~87
A^3. Its single observed hit, V164L+N167V, GROWS by +52.6 A^3, filling about half
that deficit. So the hypothesis is that the deficit sets the DIRECTION and rough
magnitude of the substitutions a ligand needs, and that ranking substitutions by
deficit-match beats the ligand-blind frequency prior.

⚠ WHY THIS IS NOT THE DEAD CAVITY-MATCH TEST (README 69). Cavity match was
scored REAL vs LIBRARY *inside* Tian's focused menus (AUC 0.540 coumarin, 0.461
PFAS) and its epitaph was "it re-measures the direction the menu has already
fixed" -- every member of a pocket-opening menu opens the pocket. Here there IS
no menu: the eugenol problem is to CHOOSE one from a single nonspecific hit, and
69c's own conclusion was that the menu supplies the DIRECTION. This asks whether
the deficit can supply that direction when no menu exists. Different question,
opposite side of the same finding.

⚠ WHY THIS IS NOT THE RETRACTED SIZE-CONDITIONED PRIOR (README 115d). That one
SUBSET the training clones to size-matched ligands, and 400 random subsets of
equal size reproduced the gain (P = 0.21). Nothing is subset here: the prior is
trained on every clone and only RE-ORDERED, so no clone is discarded and the
equal-size-subset artefact cannot arise. The analogous control for a re-ordering
is to permute the conditioning variable, which is prediction 3.

PREDICTIONS, FIXED BEFORE THE RUN:
  1. POSITIVE CONTROL. Deficit D(L) = V_ABA - V_L correlates positively with the
     mean dVolume of L's clones, p < 0.05. This is 69c's Spearman(size, dVol) =
     -0.30 restated in the deficit frame; if it fails, the pipeline is broken.
  2. THE TEST. Sign-matched re-ranking beats the blind prior at recall@20 in
     leave-one-LIGAND-out on sd03, by a paired sign test at p < 0.05.
  3. THE CONTROL. Permuting D(L) across ligands (n=200) does NOT reproduce the
     gain: observed delta must exceed the 95th percentile of the null.

⚠ Recall only, never precision -- clones are HITS FOUND, positive-unlabeled.
"""
import json, os, re, sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                          # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "fill_matching")
AA = set("ACDEFGHIKLMNPQRSTVWY")
#: side-chain volumes, A^3 (Zamyatnin), same table as 117_hit_signature.py
VOL = {"G": 60.1, "A": 88.6, "S": 89.0, "C": 108.5, "D": 111.1, "P": 112.7,
       "N": 114.1, "T": 116.1, "E": 138.4, "V": 140.0, "Q": 143.8, "H": 153.2,
       "M": 162.9, "I": 166.7, "L": 166.7, "K": 168.6, "R": 173.4, "F": 189.9,
       "Y": 193.6, "W": 227.8}
V_ABA = 248.8           # computed in-session, RDKit MMFF, ABA canonical SMILES
SMALL_HA = 15           # the eugenol cohort, fixed in the session before this run


def dvol(sub):
    """sub is like 'V163W' -> volume change of the side chain."""
    return VOL[sub[-1]] - VOL[sub[0]]


def load_clones():
    """[(ligand_name, smiles, [subs])] from sd03, one row per clone."""
    hdr, recs = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    pc = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    # sd03's ligand column is `library_name` (verified against the header:
    # library_name, parent_name, chem_library, cat, mut_lib, ... canonical_smiles)
    assert "library_name" in hdr, f"no library_name in {hdr[:6]}"
    name_col = "library_name"
    out = []
    for r in recs:
        sm = (r.get("canonical_smiles") or "").strip()
        s = [c + (r.get(c) or "").strip().upper() for c in pc
             if len((r.get(c) or "").strip()) == 1
             and (r.get(c) or "").strip().upper() in AA
             and (r.get(c) or "").strip().upper() != c[0]]
        if sm and s:
            out.append(((r.get(name_col) or sm).strip(), sm, s))
    return out


def volumes(smiles_set, cache):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    from rdkit.Chem.AllChem import ComputeMolVolume
    RDLogger.DisableLog("rdApp.*")
    if os.path.exists(cache):
        v = json.load(open(cache))
    else:
        v = {}
    new = 0
    for sm in smiles_set:
        if sm in v:
            continue
        m = Chem.MolFromSmiles(sm)
        if m is None:
            v[sm] = None; continue
        mh = Chem.AddHs(m)
        if AllChem.EmbedMolecule(mh, randomSeed=0xf00d) != 0:
            v[sm] = None; continue
        try:
            AllChem.MMFFOptimizeMolecule(mh, maxIters=400)
            v[sm] = float(ComputeMolVolume(mh))
        except Exception:                                   # noqa: BLE001
            v[sm] = None
        new += 1
        if new % 25 == 0:
            print(f"    volumes: {new} computed", flush=True)
            json.dump(v, open(cache, "w"))
    json.dump(v, open(cache, "w"))
    return v


def blind_prior(clones):
    c = Counter(s for _, _, subs in clones for s in subs)
    return [x for x, _ in c.most_common()]


def fill_prior(clones, deficit):
    """Blind frequency order, stably partitioned so that substitutions whose
    volume change MATCHES the sign of the deficit come first. No free
    parameters: the only input is sign(deficit)."""
    rank = blind_prior(clones)
    want = 1.0 if deficit > 0 else -1.0
    match = [s for s in rank if np.sign(dvol(s)) == want]
    rest = [s for s in rank if np.sign(dvol(s)) != want]
    return match + rest


def tiebreak_prior(clones, deficit):
    """Blind frequency order, with deficit-match used ONLY to break TIES.

    The hard partition in fill_prior() is net harmful (-0.065 recall@20,
    p = 0.0003) because it demotes high-frequency substitutions, and recall@20
    lives at the top of the blind ranking. This version can never reorder two
    substitutions of different frequency, so it cannot destroy the blind signal;
    it only decides the order inside frequency ties, which are large -- most
    substitutions are seen once or twice. Still parameter-free.
    """
    c = Counter(s for _, _, subs in clones for s in subs)
    want = 1.0 if deficit > 0 else -1.0
    return [s for s, _ in sorted(
        c.items(), key=lambda kv: (-kv[1], 0 if np.sign(dvol(kv[0])) == want else 1))]


def recall(rank, truth, n):
    return len(set(rank[:n]) & truth) / max(len(truth), 1)


def main():
    os.makedirs(OUT, exist_ok=True)
    clones = load_clones()
    print(f"sd03: {len(clones)} clones with >=1 substitution")
    vols = volumes({sm for _, sm, _ in clones},
                   os.path.join(OUT, "ligand_volumes.json"))
    clones = [(n, sm, s) for n, sm, s in clones if vols.get(sm)]
    print(f"  {len(clones)} clones with a computed ligand volume")

    by_lig = defaultdict(list)
    for n, sm, s in clones:
        by_lig[sm].append((n, sm, s))
    print(f"  {len(by_lig)} distinct ligands")

    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    ha = {sm: Chem.MolFromSmiles(sm).GetNumHeavyAtoms() for sm in by_lig}

    # ---- prediction 1: positive control --------------------------------------
    from scipy.stats import spearmanr, wilcoxon
    D, MD = [], []
    for sm, cl in by_lig.items():
        D.append(V_ABA - vols[sm])
        MD.append(float(np.mean([dvol(s) for _, _, subs in cl for s in subs])))
    r, p = spearmanr(D, MD)
    print(f"\nP1 positive control: Spearman(deficit, mean dVolume) = "
          f"{r:+.3f}, p = {p:.2g}   "
          f"{'HOLDS' if (r > 0 and p < 0.05) else 'FAILS'}")

    # ---- prediction 2: leave-one-ligand-out ----------------------------------
    rows = []
    for sm, cl in by_lig.items():
        train = [c for c in clones if c[1] != sm]
        truth = {s for _, _, subs in cl for s in subs}
        d = V_ABA - vols[sm]
        b = blind_prior(train)
        f = fill_prior(train, d)
        t = tiebreak_prior(train, d)
        rows.append(dict(name=cl[0][0], ha=ha[sm], vol=round(vols[sm], 1),
                         deficit=round(d, 1), n=len(truth),
                         b20=recall(b, truth, 20), f20=recall(f, truth, 20),
                         b40=recall(b, truth, 40), f40=recall(f, truth, 40),
                         t20=recall(t, truth, 20), t40=recall(t, truth, 40)))
    json.dump(rows, open(os.path.join(OUT, "loo_rows.json"), "w"), indent=1)

    def report(sel, lab):
        if len(sel) < 3:
            print(f"  {lab}: n={len(sel)}, too few"); return None
        for k in ("20", "40"):
          for arm in ("f", "t"):
            b = np.array([x["b" + k] for x in sel])
            f = np.array([x[arm + k] for x in sel])
            d = f - b
            w = sum(1 for x in d if x > 0); l = sum(1 for x in d if x < 0)
            try:
                p = wilcoxon(f, b).pvalue if np.any(d) else 1.0
            except Exception:                               # noqa: BLE001
                p = 1.0
            nm = "partition" if arm == "f" else "TIEBREAK "
            print(f"  {lab:<20}@{k} {nm}: blind {b.mean():.3f} -> "
                  f"{f.mean():.3f}  delta {d.mean():+.3f}  "
                  f"win/loss {w}/{l}  p={p:.3g}")
        return None

    print(f"\nP2 leave-one-ligand-out ({len(rows)} ligands)")
    report(rows, "ALL ligands")
    report([x for x in rows if x["ha"] <= SMALL_HA], f"small (<={SMALL_HA} HA)")
    report([x for x in rows if x["ha"] > SMALL_HA], f"large (>{SMALL_HA} HA)")

    # ---- prediction 3: permute the conditioning variable ---------------------
    rng = np.random.default_rng(0)
    obs = float(np.mean([x["f20"] - x["b20"] for x in rows]))
    defs = np.array([x["deficit"] for x in rows])
    keys = list(by_lig)
    null = []
    for it in range(200):
        perm = rng.permutation(len(keys))
        tot = 0.0
        for i, sm in enumerate(keys):
            cl = by_lig[sm]
            truth = {s for _, _, subs in cl for s in subs}
            train = [c for c in clones if c[1] != sm]
            b = blind_prior(train)
            f = fill_prior(train, defs[perm[i]])
            tot += recall(f, truth, 20) - recall(b, truth, 20)
        null.append(tot / len(keys))
        if (it + 1) % 50 == 0:
            print(f"    permutation {it+1}/200", flush=True)
    null = np.array(null)
    q95 = float(np.percentile(null, 95))
    pperm = float((null >= obs).mean())
    print(f"\nP3 permutation control (n=200): observed delta@20 {obs:+.4f}, "
          f"null mean {null.mean():+.4f}, null 95th {q95:+.4f}, "
          f"P(null >= obs) = {pperm:.3f}")
    print(f"   {'PASSES -- not reproducible by permuting the deficit' if pperm < 0.05 else 'FAILS -- permuted deficits do as well'}")

    # ---- eugenol, named ------------------------------------------------------
    eu = [x for x in rows if "eugenol" in x["name"].lower()]
    if eu:
        x = eu[0]
        print(f"\nEUGENOL: {x['ha']} heavy atoms, vol {x['vol']} A^3, "
              f"deficit {x['deficit']:+.1f} A^3, {x['n']} observed subs, "
              f"blind@20 {x['b20']:.2f} -> fill@20 {x['f20']:.2f}")
    json.dump(dict(p1_r=r, p1_p=p, obs_delta20=obs, null_mean=float(null.mean()),
                   null_q95=q95, p_perm=pperm, n_ligands=len(rows)),
              open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    print(f"\nwrote {OUT}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
