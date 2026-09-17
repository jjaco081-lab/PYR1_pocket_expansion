#!/usr/bin/env python
r"""
230_guaiacol_positions.py -- WHICH POCKET POSITIONS DO EUGENOL-LIKE LIGANDS USE?

Jannis: "can we maybe compare among structures with a phenol and an ether in a
similar geometry to eugenol and see what their residues near the mouth of the
pocket are if we assume it binds there?"

Eugenol is 4-allyl-2-METHOXYPHENOL: an aromatic ring carrying a hydroxyl and an
ether on ADJACENT carbons. That ortho arrangement is the guaiacol motif, and it
is what makes eugenol's two candidate anchors (phenol OH, ether O) sit 2.8 A
apart on a rigid ring rather than at opposite ends of a flexible molecule.

THE IDEA. Eugenol has exactly one known sensor (V164L+N167V, 100 uM) and it is a
GENERIC small-ligand solution shared with five unrelated ligands, so on its own it
says nothing about where eugenol sits. But the Tian screen contains other ligands
carrying the same motif, each with its own evolved sensors. If those sensors
converge on a particular set of positions, that is independent evidence about
where this chemistry is read -- evidence that does not depend on docking eugenol,
which 209 showed is pose-underdetermined (117 of 180 placements clash-free).

FOUR NESTED CLASSES, so the comparison is graded rather than binary:
  guaiacol   phenol + ether on ADJACENT ring carbons  (eugenol's own motif)
  phenol_ether  phenol and an ether, not necessarily ortho
  phenol_only   a phenol, no ether
  ether_only    an ether, no phenol
  neither       the rest of the screen -- the background

⚠ THE CONFOUND THAT KILLED SIX PREVIOUS ATTEMPTS. Subsetting a library by
chemistry and finding "enrichment" is not evidence: a smaller sample looks
different by chance, and §115d / §206 both died on exactly this. So every class
is compared against a RANDOM SUBSET OF EQUAL SIZE drawn from the same sensor
pool, not against the pooled remainder. A class only counts if it beats its own
size-matched null.

⚠ RECALL ONLY. Sensors are hits FOUND, so a position not mutated may simply not
have been sampled. Positions that ARE mutated are informative; absences are not.

⚠ DEPTH. Mouth-vs-deep is computed here from the structure, not taken from the
numbers quoted in 209's header, so it can be checked. Depth = distance from the
residue's closest-to-ABA atom to the pocket MOUTH, defined as the centroid of the
gate and latch loops (the end ABA's carboxylate points toward and the end HAB1
seals). Larger = deeper.
"""
import json, os, sys, random
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
lx = import_module("lib_xlsx")
BENCH = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark/"
SD03 = BENCH + "pnas.2519924122.sd03(1).xlsx"
OUT = os.path.join(ROOT, "results", "guaiacol")
POS = ["K59", "V81", "V83", "L87", "A89", "S92", "E94", "F108", "I110", "L117",
       "Y120", "S122", "E141", "F159", "A160", "V163", "V164", "N167"]
#: gate and latch, from project convention; identities asserted below
GATE = list(range(85, 95))
LATCH = list(range(112, 122))
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
         "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
         "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V"}
SMARTS = {
    "guaiacol":     "[OX2H]c1ccccc1[OX2][CX4]",
    "phenol":       "[OX2H]c",
    "ether":        "[OX2]([#6])[#6]",
    "phenol_arom_ether": "[OX2H]c1ccccc1",       # used with ether separately
}


def depths():
    """{resnum: (identity, depth_A)} for the 18 lining positions."""
    atoms, lig = defaultdict(list), []
    for l in open(os.path.join(ROOT, "data", "3QN1_complex_auth_aba.pdb")):
        if l.startswith("ATOM") and l[21] == "A":
            atoms[int(l[22:26])].append(
                (l[17:20].strip(),
                 np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
        elif l.startswith("HETATM") and l[17:20] == "A8S":
            lig.append(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    lig = np.array(lig)
    mouth = np.mean([x for n in GATE + LATCH for _, x in atoms.get(n, [])], axis=0)
    out = {}
    for p in POS:
        n = int(p[1:])
        assert n in atoms, f"residue {n} absent"
        ident = THREE.get(atoms[n][0][0], "?")
        assert ident == p[0], f"residue {n} IS {ident}, name says {p[0]}"
        near = min(atoms[n], key=lambda a: np.linalg.norm(lig - a[1], axis=1).min())
        out[p] = (ident, float(np.linalg.norm(near[1] - mouth)))
    return out


def classify(smiles, cache={}):
    if smiles in cache:
        return cache[smiles]
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        cache[smiles] = None; return None
    has = {k: m.HasSubstructMatch(Chem.MolFromSmarts(v))
           for k, v in SMARTS.items()}
    if has["guaiacol"]:
        c = "guaiacol"
    elif has["phenol"] and has["ether"]:
        c = "phenol_ether"
    elif has["phenol"]:
        c = "phenol_only"
    elif has["ether"]:
        c = "ether_only"
    else:
        c = "neither"
    cache[smiles] = c
    return c


def main():
    os.makedirs(OUT, exist_ok=True)
    D = depths()
    print("LINING POSITIONS BY DEPTH (identities asserted; mouth = gate+latch centroid)")
    for p, (i, d) in sorted(D.items(), key=lambda x: x[1][1]):
        print(f"   {p:<6} {i}  depth {d:5.2f} A")
    med = float(np.median([d for _, d in D.values()]))
    MOUTH = {p for p, (_, d) in D.items() if d < med}
    print(f"\n   median depth {med:.2f} A -> MOUTH half: "
          f"{' '.join(sorted(MOUTH, key=lambda p: D[p][1]))}\n")

    rows = lx.table(SD03)
    hdr = rows[0]
    recs = rows[1] if len(rows) > 1 else []
    print(f"sd03: {len(recs)} sensor clones")
    byclass = defaultdict(list)
    ligclass = {}
    for r in recs:
        smi = r.get("canonical_smiles", "")
        if not smi:
            continue
        c = classify(smi)
        if c is None:
            continue
        ligclass[r["library_name"]] = c
        muts = [p for p in POS if r.get(p, "")]
        if muts:
            byclass[c].append((r["library_name"], muts))
    nlig = defaultdict(set)
    for lg, c in ligclass.items():
        nlig[c].add(lg)
    print(f"{'class':<14}{'ligands':>9}{'sensors':>9}")
    for c in ("guaiacol", "phenol_ether", "phenol_only", "ether_only", "neither"):
        print(f"{c:<14}{len(nlig[c]):>9}{len(byclass[c]):>9}")
    print(f"\nguaiacol ligands: {sorted(nlig['guaiacol'])}")
    print(f"phenol_ether    : {sorted(nlig['phenol_ether'])[:14]}")

    allsens = [s for c in byclass for s in byclass[c]]
    rng = random.Random(0)

    def profile(sens):
        f = defaultdict(int)
        for _, muts in sens:
            for m in muts:
                f[m] += 1
        n = max(len(sens), 1)
        return {p: f[p] / n for p in POS}

    print(f"\n{'position':<9}{'depth':>7}", end="")
    for c in ("guaiacol", "phenol_ether", "phenol_only", "neither"):
        print(f"{c[:9]:>11}", end="")
    print("   size-matched null p (guaiacol)")
    results = {}
    for c in ("guaiacol", "phenol_ether", "phenol_only", "neither"):
        results[c] = profile(byclass[c])
    gp = results["guaiacol"]
    ng = len(byclass["guaiacol"])
    nullp = {}
    if ng >= 3:
        draws = [profile(rng.sample(allsens, ng)) for _ in range(5000)]
        for p in POS:
            obs = gp[p]
            nullp[p] = (sum(1 for d in draws if d[p] >= obs) + 1) / 5001
    for p in sorted(POS, key=lambda p: D[p][1]):
        print(f"{p:<9}{D[p][1]:>7.2f}", end="")
        for c in ("guaiacol", "phenol_ether", "phenol_only", "neither"):
            print(f"{results[c][p]*100:>10.0f}%", end="")
        star = ""
        if p in nullp:
            star = f"   p={nullp[p]:.4f}" + ("  <--" if nullp[p] < 0.05 else "")
        print(star)

    if ng >= 3:
        mo = sum(gp[p] for p in MOUTH) / max(sum(gp[p] for p in POS), 1e-9)
        nulls = [sum(d[p] for p in MOUTH) / max(sum(d[p] for p in POS), 1e-9)
                 for d in draws]
        pm = (sum(1 for x in nulls if x >= mo) + 1) / (len(nulls) + 1)
        print(f"\nfraction of guaiacol mutations in the MOUTH half: {mo:.2f}"
              f"   size-matched null {np.mean(nulls):.2f}   p = {pm:.4f}")
    else:
        print(f"\n⚠ only {ng} guaiacol sensors -- too few for a null. "
              f"Report the ligand list, not a statistic.")
    json.dump(dict(depths={p: D[p][1] for p in POS},
                   n_sensors={c: len(byclass[c]) for c in byclass},
                   ligands={c: sorted(nlig[c]) for c in nlig},
                   profiles=results, null_p=nullp),
              open(os.path.join(OUT, "guaiacol_positions.json"), "w"), indent=1)
    print(f"\nwrote {OUT}/guaiacol_positions.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
