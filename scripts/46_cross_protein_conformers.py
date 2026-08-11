#!/usr/bin/env python
"""
46_cross_protein_conformers.py -- does a ligand adopt the SAME conformation in
unrelated receptors, or is its bound conformation receptor-specific?

WHY THIS IS THE DECISIVE HALF
-----------------------------
§45 measured conformational spread between copies inside ONE crystal, which
isolates conformational variability from receptor identity. That is the wrong
control for the design question. What matters for building libraries around known
bound conformations is whether a conformation TRANSFERS:

  * same conformation across unrelated proteins -> it is intrinsically preferred,
    a genuine low-energy bioactive conformer, and a sound anchor for a PYR1 design
  * receptor-specific conformation             -> a conformer seen in receptor X
    is no evidence it is reachable or favourable in PYR1, and anchoring on it
    imports another protein's induced fit

So this compares copies drawn from DIFFERENT PDB entries with DIFFERENT UniProt
accessions, and contrasts that with the within-crystal spread for the same ligand.

SCOPE
-----
Restricted to ligands with >=40 heavy atoms -- the size range this project cares
about, and where §45 found the conformational tail. Up to MAX_ENTRIES structures
per ligand.

METRIC
------
For each ligand: pool every heavy-atom copy, tag it with its UniProt accession,
and compute Kabsch RMSD (matched by PDB atom name) for
  * within_max   -- worst pair sharing an accession
  * across_med / across_max -- pairs spanning different accessions
`across` much larger than `within` means induced fit dominates. Comparable values
mean the conformation is intrinsic and transferable.

LIMITS
------
  * Atom-name matching ignores symmetry, so all RMSDs are upper bounds.
  * No resolution or occupancy filter.
  * A ligand seen in only one protein family cannot be tested and is reported as
    n_proteins=1 rather than silently dropped.

Usage:  python 46_cross_protein_conformers.py [--min-heavy 40] [--max-entries 4]
"""
import argparse, csv, itertools, json, os, time, urllib.request
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "survey_cache")
SURVEY = os.path.join(ROOT, "results", "44_floppy_survey.csv")
OUT = os.path.join(ROOT, "results", "46_cross_protein_conformers.csv")
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"

ap = argparse.ArgumentParser()
ap.add_argument("--min-heavy", type=int, default=40)
ap.add_argument("--max-entries", type=int, default=4)
args = ap.parse_args()
os.makedirs(CACHE, exist_ok=True)

_uni = {}


def post(url, payload, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                         headers={"Content-Type": "application/json"})
            return json.load(urllib.request.urlopen(req, timeout=30))
        except Exception:
            time.sleep(1 + i)
    return {}


def get(url, tries=2):
    for i in range(tries):
        try:
            return json.load(urllib.request.urlopen(url, timeout=30))
        except Exception:
            time.sleep(1 + i)
    return {}


def entries_for(comp, rows):
    q = {"query": {"type": "terminal", "service": "text_chem",
                   "parameters": {"attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                                  "operator": "exact_match", "value": comp}},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": rows}}}
    return [r["identifier"] for r in post(SEARCH, q).get("result_set", [])]


def uniprot(entry):
    """Primary UniProt accession of entity 1; falls back to the entry id."""
    if entry in _uni:
        return _uni[entry]
    d = get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry}/1")
    acc = None
    ids = (d.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids")
    if ids:
        acc = ids[0]
    _uni[entry] = acc or f"ENTRY:{entry}"
    return _uni[entry]


def fetch(pdb):
    p = os.path.join(CACHE, f"{pdb}.cif")
    if not os.path.exists(p):
        try:
            with urllib.request.urlopen(
                    f"https://files.rcsb.org/download/{pdb}.cif", timeout=60) as r:
                open(p, "wb").write(r.read())
        except Exception:
            return None
    return p if os.path.getsize(p) > 5000 else None


def copies(cif, comp):
    out = {}
    for line in open(cif, errors="ignore"):
        if not line.startswith("HETATM"):
            continue
        f = line.split()
        if len(f) < 15 or f[5] != comp or f[2] in ("H", "D"):
            continue
        try:
            xyz = (float(f[10]), float(f[11]), float(f[12]))
        except ValueError:
            continue
        out.setdefault((f[6], f[8], f[4]), {})[f[3]] = xyz
    return {k: v for k, v in out.items() if len(v) >= 8}


def kabsch(P, Q):
    P = P - P.mean(0); Q = Q - Q.mean(0)
    V, S, W = np.linalg.svd(P.T @ Q)
    if np.linalg.det(V) * np.linalg.det(W) < 0:
        V[:, -1] = -V[:, -1]
    return float(np.sqrt((((P @ (V @ W)) - Q) ** 2).sum(1).mean()))


base = {r["comp"]: r for r in csv.DictReader(open(SURVEY))
        if int(r["n_heavy"]) >= args.min_heavy}
print(f"{len(base)} ligands with >= {args.min_heavy} heavy atoms\n", flush=True)

recs = []
for i, (comp, r) in enumerate(sorted(base.items())):
    ents = entries_for(comp, args.max_entries) or [r["entry"]]
    pool = []                       # (accession, entry, atoms)
    for e in ents[:args.max_entries]:
        c = fetch(e)
        if not c:
            continue
        acc = uniprot(e)
        for k, atoms in copies(c, comp).items():
            pool.append((acc, e, atoms))
    if len(pool) < 2:
        continue
    accs = {a for a, _, _ in pool}
    within, across = [], []
    for (a1, e1, x), (a2, e2, y) in itertools.combinations(pool, 2):
        common = sorted(set(x) & set(y))
        if len(common) < 8:
            continue
        v = kabsch(np.array([x[n] for n in common], float),
                   np.array([y[n] for n in common], float))
        (within if a1 == a2 else across).append(v)
    recs.append(dict(comp=comp, n_heavy=int(r["n_heavy"]),
                     shape="", buried=float(r["buried"]),
                     n_entries=len({e for _, e, _ in pool}), n_proteins=len(accs),
                     n_copies=len(pool),
                     within_max=round(max(within), 3) if within else None,
                     across_med=round(float(np.median(across)), 3) if across else None,
                     across_max=round(max(across), 3) if across else None))
    if (i + 1) % 25 == 0:
        print(f"  {i+1}/{len(base)} processed", flush=True)

# attach shape from 45 if available
try:
    sh = {x["comp"]: x["shape"] for x in
          csv.DictReader(open(os.path.join(ROOT, "results", "45_ligand_conformations.csv")))}
    for x in recs:
        x["shape"] = sh.get(x["comp"], "")
except OSError:
    pass

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(recs[0]))
    w.writeheader(); w.writerows(recs)

# ------------------------------- report -------------------------------
import statistics as st
multi = [x for x in recs if x["n_proteins"] > 1 and x["across_med"] is not None]
print(f"\n{len(recs)} ligands analysed; {len(multi)} seen in >1 distinct protein\n")
if multi:
    print("CONFORMATION TRANSFER ACROSS UNRELATED PROTEINS")
    print(f"{'heavy':>10} {'n':>4} {'within_max':>11} {'across_med':>11} {'across_max':>11} {'across<1A':>10}")
    for lo, hi in [(40, 55), (55, 999)]:
        g = [x for x in multi if lo <= x["n_heavy"] < hi]
        if len(g) < 3:
            continue
        wi = [x["within_max"] for x in g if x["within_max"] is not None]
        print(f"{lo:>4}-{hi:<5} {len(g):>4} "
              f"{(st.median(wi) if wi else float('nan')):>11.2f} "
              f"{st.median(x['across_med'] for x in g):>11.2f} "
              f"{st.median(x['across_max'] for x in g):>11.2f} "
              f"{sum(x['across_med']<1.0 for x in g)/len(g):>10.0%}")
    for s in ("rod", "disc"):
        g = [x for x in multi if x["shape"] == s]
        if len(g) >= 3:
            print(f"  {s:>5}: n={len(g):>3} across_med median {st.median(x['across_med'] for x in g):.2f}, "
                  f"{sum(x['across_med']<1.0 for x in g)/len(g):.0%} under 1 A")
    pairs = [(x["within_max"], x["across_med"]) for x in multi if x["within_max"] is not None]
    if len(pairs) > 5:
        a = np.array([p for p, _ in pairs]); b = np.array([q for _, q in pairs])
        print(f"\n  median within_max {np.median(a):.2f} A  vs  median across_med {np.median(b):.2f} A")
        print(f"  across exceeds within for {np.mean(b > a):.0%} of ligands")
    print("\nINTERPRETATION")
    print("  across ~ within  -> conformation is intrinsic and TRANSFERS; anchor on it")
    print("  across >> within -> induced fit dominates; a foreign conformer is not evidence")
print(f"\nwrote {OUT}")
