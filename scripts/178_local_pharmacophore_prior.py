#!/usr/bin/env python
r"""
178_local_pharmacophore_prior.py -- similarity computed PER POCKET POSITION.

Jannis's refinement of §103. Whole-molecule Tanimoto failed for a specific
reason: it is GLOBAL. A novel ligand can have no global analog while every one of
its functional groups has many local ones -- which is exactly the §103a row where
the prior went NEGATIVE (max T < 0.2, delta -0.036). The fix is to stop asking
"which ligand is this like?" and start asking, per position, "what does this
ligand present HERE, and which ligands presented the same thing?"

That is also what an MSA actually does: position-specific scoring, not
whole-sequence matching.

DATA, and it already exists: 181 of the 194 sd03 ligands with clones have a
Boltz co-folded pose from §84, covering 652 of 691 clones, all superposable into
one PYR1 frame and all past §93's fold assertions.

METHOD
  typing      each ligand heavy atom is typed by ELEMENT -- polar (N, O), sulfur,
              halogen (F, Cl, Br, I), carbon. Element typing is used rather than
              an RDKit pharmacophore because it needs no atom-order mapping
              between the SMILES and the CIF, which is the step most likely to
              fail silently.
  feature     for each library position p, the vector of type counts within
              SHELL A of p's C-beta -- what the ligand PRESENTS at p
  local sim   cosine between the query's and a training ligand's feature vector
              AT p. A ligand can be locally similar at 108 and globally unlike.
  prior       P(sub at p | q) = sum_j localsim_p(q,j)^alpha * freq(sub at p | j)

⚠ LEAVE-ONE-LIGAND-OUT IS ASSERTED, not assumed. Jannis's condition: a held-out
ligand must never contribute to its own prior. Splits are on LIGANDS, and the
code asserts the held-out ligand's index is absent from the training set before
每 prior is built.

BAR, fixed before results. It must beat BOTH:
  * ligand-blind frequency          recall@20 = 0.345  (§52, §103)
  * whole-molecule ECFP4 prior      recall@20 = 0.398  (§103)
Beating blind alone is not a result; §103 already does that.

THE DECISIVE CONTROL is §103a's: delta stratified by max whole-molecule Tanimoto
to the training set. If the gain again lives only in the T >= 0.5 stratum, this
is lookup in new clothing. The claim holds only if it gains where T < 0.2, the
stratum where whole-molecule similarity actively HURT.

⚠ Unchanged by any of this: §87b's 67 % vocabulary ceiling. Re-weighting reorders
within the vocabulary and cannot add a position no library ever varied.
"""
import glob, json, os, re, sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                          # noqa: E402
from importlib import import_module                                 # noqa: E402
G = import_module("164_pose_geometry")

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
T = os.path.join(ROOT, "data", "tractability")
OUT = os.path.join(ROOT, "results", "ligand_prior")
AA = set("ACDEFGHIKLMNPQRSTVWY")
SHELL = float(os.environ.get("SHELL","7.0"))
TYPES = ["polar", "halogen", "sulfur", "carbon"]
NSPLIT = 50


def etype(e):
    if e in ("N", "O"): return "polar"
    if e in ("F", "CL", "BR", "I", "Cl", "Br"): return "halogen"
    if e == "S": return "sulfur"
    return "carbon"


def main():
    os.makedirs(OUT, exist_ok=True)
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")

    # ---- clones and SMILES ------------------------------------------------
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
    POSITIONS = sorted({int(s[1:-1]) for v in clone.values() for c in v for s in c})

    # ---- poses ------------------------------------------------------------
    key = json.load(open(os.path.join(T, "key_matched.json")))
    byname = {v["name"].strip().lower(): k for k, v in key.items()}
    fbb, aba, _ = G.frame()
    core = [r for r in fbb if 6 <= r <= 180 and r not in G.GATE + G.LATCH]
    cb = {}
    for l in open(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")):
        if l.startswith("ATOM") and l[12:16].strip() in ("CB", "CA"):
            n = int(l[22:26])
            if l[12:16].strip() == "CB" or n not in cb:
                cb[n] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])

    ligs, feat = [], {}
    for lg in sorted(clone):
        k = byname.get(lg.strip().lower())
        if not k or lg not in smi:
            continue
        g = sorted(glob.glob(os.path.join(T, "out", f"boltz_results_{k}",
                                          "predictions", "*", "*.cif")))
        if not g:
            continue
        bb, L, _, pl = G.read(g[0])
        if len(L) == 0 or pl < 80:
            continue
        el = _elements(g[0])
        if len(el) != len(L):
            continue
        P, Q = G.paired(fbb, bb, core)
        if len(P) < 400:
            continue
        R, mp, mq = G.kabsch(P, Q)
        Lf = (L - mq) @ R + mp
        v = {}
        for p in POSITIONS:
            if p not in cb:
                continue
            d = np.linalg.norm(Lf - cb[p], axis=1)
            c = Counter(etype(el[i]) for i in np.where(d < SHELL)[0])
            v[p] = np.array([c.get(t, 0) for t in TYPES], float)
        feat[lg] = v
        ligs.append(lg)
    print(f"{len(ligs)} ligands with pose + clones; {len(POSITIONS)} library positions; "
          f"shell {SHELL} A")

    fp = {l: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi[l]), 2, 2048)
          for l in ligs}
    n = len(ligs)
    Tan = np.array([DataStructs.BulkTanimotoSimilarity(fp[a], [fp[b] for b in ligs])
                    for a in ligs])
    # substitutions split by position
    subp = {l: defaultdict(Counter) for l in ligs}
    allsub = {l: Counter() for l in ligs}
    for l in ligs:
        for c in clone[l]:
            for s in c:
                subp[l][int(s[1:-1])][s] += 1
                allsub[l][s] += 1
    truth = {l: set(allsub[l]) for l in ligs}
    idx = {l: i for i, l in enumerate(ligs)}

    def cos(a, b):
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        return float(a @ b / (na * nb)) if na > 0 and nb > 0 else 0.0

    def rec(rk, t, k):
        return len(set(rk[:k]) & t) / max(len(t), 1)

    rng = np.random.default_rng(0)
    out = defaultdict(list); strat = []
    for rep in range(NSPLIT):
        perm = rng.permutation(n)
        tr, te = list(perm[: n // 2]), list(perm[n // 2:])
        trs = set(tr)
        blind = Counter()
        for j in tr:
            blind.update(allsub[ligs[j]])
        rb = [s for s, _ in blind.most_common()]
        for j in te:
            l = ligs[j]
            assert j not in trs, "held-out ligand leaked into training"
            if not truth[l]:
                continue
            # whole-molecule prior (103)
            w = defaultdict(float)
            for k in tr:
                s2 = Tan[j, k] ** 2
                if s2 > 0:
                    for s, c in allsub[ligs[k]].items():
                        w[s] += s2 * c
            rw = [s for s, _ in sorted(w.items(), key=lambda kv: -kv[1])]
            # LOCAL prior -- similarity recomputed per position
            wl = defaultdict(float)
            for p in POSITIONS:
                if p not in feat[l]:
                    continue
                for k in tr:
                    if p not in feat[ligs[k]]:
                        continue
                    sim = cos(feat[l][p], feat[ligs[k]][p])
                    if sim <= 0:
                        continue
                    for s, c in subp[ligs[k]][p].items():
                        wl[s] += (sim ** 2) * c
            rl = [s for s, _ in sorted(wl.items(), key=lambda kv: -kv[1])]
            for nn in (20, 40):
                out[("blind", nn)].append(rec(rb, truth[l], nn))
                out[("whole-molecule", nn)].append(rec(rw, truth[l], nn))
                out[("LOCAL pharmacophore", nn)].append(rec(rl, truth[l], nn))
            strat.append((max(Tan[j, k] for k in tr),
                          rec(rl, truth[l], 20) - rec(rb, truth[l], 20),
                          rec(rw, truth[l], 20) - rec(rb, truth[l], 20)))

    print(f"\nheld-out recall over {NSPLIT} 50/50 LIGAND splits "
          f"(leave-one-ligand-out asserted)")
    print(f"{'prior':<24}{'recall@20':>14}{'recall@40':>14}")
    for nm in ("blind", "whole-molecule", "LOCAL pharmacophore"):
        a = np.array(out[(nm, 20)]); b = np.array(out[(nm, 40)])
        print(f"{nm:<24}{a.mean():>14.3f}{b.mean():>14.3f}")

    S = np.array(strat)
    print(f"\nTHE DECISIVE CONTROL -- delta vs blind, by max whole-molecule Tanimoto")
    print(f"{'stratum':<26}{'n':>6}{'LOCAL':>10}{'whole-mol':>12}")
    for lo, hi, lab in ((0, .2, "< 0.2  (no analog)"), (.2, .3, "0.2-0.3"),
                        (.3, .5, "0.3-0.5"), (.5, 1.01, ">= 0.5 (close analog)")):
        m = (S[:, 0] >= lo) & (S[:, 0] < hi)
        if m.sum():
            print(f"  {lab:<24}{int(m.sum()):>6}{S[m,1].mean():>10.4f}{S[m,2].mean():>12.4f}")
    print("\n  The claim holds only if LOCAL gains in the < 0.2 stratum, where")
    print("  whole-molecule similarity was NEGATIVE (README 103a).")
    json.dump({str(k): float(np.mean(v)) for k, v in out.items()},
              open(os.path.join(OUT, "local_prior.json"), "w"), indent=1)
    return 0


def _elements(path):
    cols, rows, inl = {}, [], False
    for line in open(path):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".")[1]] = len(cols); inl = True; continue
        if inl:
            if line.startswith(("#", "loop_", "_")):
                if rows: break
                continue
            f = line.split()
            if len(f) >= len(cols): rows.append(f)
    return [f[cols["type_symbol"]] for f in rows
            if f[cols["group_PDB"]] != "ATOM" and f[cols["type_symbol"]] != "H"]


if __name__ == "__main__":
    sys.exit(main())
