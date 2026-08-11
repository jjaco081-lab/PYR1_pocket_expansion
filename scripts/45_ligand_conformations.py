#!/usr/bin/env python
"""
45_ligand_conformations.py -- what conformation are bound ligands actually in,
and is it one conformation or many?

WHY THIS EXISTS
---------------
§23d established that receptors do engulf large ligands. That says nothing about
WHICH conformation gets engulfed. The proposal this script tests is to build
libraries around **experimentally observed bound conformations** rather than
around generated conformer ensembles.

That proposal is only sound if bound conformations are reproducible. If a ligand
adopts essentially one conformation wherever it is bound, that conformation is a
privileged design anchor -- it is demonstrably achievable, already "paid for"
entropically, and collapses the conformer-enumeration problem that makes floppy
rods intractable (the step Leonard et al. needed replica-exchange MD for). If
instead the same ligand is found in many different conformations, then no single
observed conformer is privileged and designing around one is arbitrary.

WHAT IS MEASURED
----------------
For every heavy-atom copy of each ligand in its cached entry:

  * SHAPE, from the inertia tensor -- NPR1 = I1/I3, NPR2 = I2/I3 with I1<=I2<=I3.
    rod NPR1~0 NPR2~1; disc ~0.5/~0.5; sphere ~1/~1. Classified on the standard
    triangle so "what conformation" gets a concrete answer.
  * EXTENSION -- max interatomic distance, and radius of gyration normalised by
    that of a compact sphere of the same atom count, so elongation is comparable
    across sizes.
  * CONFORMATIONAL SPREAD (the decisive number) -- for ligands present in more
    than one copy, pairwise heavy-atom RMSD after Kabsch superposition, matched
    by PDB atom name. This is WITHIN one crystal, so it isolates conformational
    variability from any difference in receptor.

LIMITS
------
  * Atom-name matching ignores molecular symmetry, so a symmetric group flipped
    180 degrees inflates RMSD. Reported spreads are therefore upper bounds.
  * Copies within one entry may be crystallographically related rather than
    independent, which biases spread DOWN. The two biases oppose.
  * No resolution or occupancy filter; ligand geometry is weakly restrained at
    poor resolution.
  * Within-entry only. Across-entry comparison needs the other entries downloaded.

Usage:  python 45_ligand_conformations.py
"""
import csv, itertools, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "survey_cache")
SURVEY = os.path.join(ROOT, "results", "44_floppy_survey.csv")
OUT = os.path.join(ROOT, "results", "45_ligand_conformations.csv")


def copies(cif, comp):
    """{(chain, resseq, altloc): {atom_name: xyz}} for every copy of comp."""
    out = {}
    for line in open(cif, errors="ignore"):
        if not line.startswith("HETATM"):
            continue
        f = line.split()
        if len(f) < 15 or f[5] != comp:
            continue
        el, name = f[2], f[3]
        if el == "H" or el == "D":
            continue
        try:
            xyz = (float(f[10]), float(f[11]), float(f[12]))
        except ValueError:
            continue
        alt = f[4] if f[4] not in (".", "?") else ""
        key = (f[6], f[8], alt)          # auth chain, auth seq, altloc
        out.setdefault(key, {})[name] = xyz
    return {k: v for k, v in out.items() if len(v) >= 8}


def shape(xyz):
    """NPR1, NPR2, max interatomic distance, normalised radius of gyration."""
    c = xyz - xyz.mean(0)
    # inertia tensor with unit masses
    I = np.zeros((3, 3))
    for r in c:
        I += np.dot(r, r) * np.eye(3) - np.outer(r, r)
    ev = np.sort(np.linalg.eigvalsh(I))
    if ev[2] <= 0:
        return None
    npr1, npr2 = ev[0] / ev[2], ev[1] / ev[2]
    d = np.sqrt(((c[:, None, :] - c[None, :, :]) ** 2).sum(-1))
    rg = np.sqrt((c ** 2).sum(1).mean())
    # Rg of a uniform sphere holding n atoms at ~20.1 A^3 each
    rg_sphere = np.sqrt(3 / 5) * (3 * len(c) * 20.1 / (4 * np.pi)) ** (1 / 3)
    return npr1, npr2, float(d.max()), float(rg / rg_sphere)


def classify(npr1, npr2):
    """Corner of the NPR triangle this conformation sits nearest."""
    d = {"rod": np.hypot(npr1 - 0.0, npr2 - 1.0),
         "disc": np.hypot(npr1 - 0.5, npr2 - 0.5),
         "sphere": np.hypot(npr1 - 1.0, npr2 - 1.0)}
    return min(d, key=d.get)


def kabsch_rmsd(P, Q):
    P = P - P.mean(0); Q = Q - Q.mean(0)
    V, S, W = np.linalg.svd(P.T @ Q)
    if np.linalg.det(V) * np.linalg.det(W) < 0:
        V[:, -1] = -V[:, -1]
    P = P @ (V @ W)
    return float(np.sqrt(((P - Q) ** 2).sum(1).mean()))


rows = list(csv.DictReader(open(SURVEY)))
recs = []
for r in rows:
    cif = os.path.join(CACHE, f"{r['entry']}.cif")
    if not os.path.exists(cif):
        continue
    cps = copies(cif, r["comp"])
    if not cps:
        continue
    shapes, names = [], list(cps)
    for k in names:
        a = cps[k]
        s = shape(np.array([a[n] for n in sorted(a)], float))
        if s:
            shapes.append(s)
    if not shapes:
        continue
    npr1, npr2, maxd, rgn = np.mean(shapes, axis=0)

    # pairwise RMSD between copies, matched by atom name
    spread = np.nan
    if len(names) > 1:
        vals = []
        for x, y in itertools.combinations(names, 2):
            common = sorted(set(cps[x]) & set(cps[y]))
            if len(common) < 8:
                continue
            vals.append(kabsch_rmsd(np.array([cps[x][n] for n in common], float),
                                    np.array([cps[y][n] for n in common], float)))
        if vals:
            spread = float(np.max(vals))       # worst pair = full range seen
    recs.append(dict(comp=r["comp"], entry=r["entry"], n_heavy=int(r["n_heavy"]),
                     n_rot=r["n_rot"], buried=float(r["buried"]),
                     n_copies=len(names), npr1=round(npr1, 3), npr2=round(npr2, 3),
                     shape=classify(npr1, npr2), max_dim=round(maxd, 2),
                     rg_norm=round(rgn, 3),
                     spread_rmsd=None if np.isnan(spread) else round(spread, 3)))

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(recs[0]))
    w.writeheader()
    w.writerows(recs)

# ------------------------------- report -------------------------------
import collections, statistics as st
print(f"{len(recs)} ligands with usable coordinates\n")

print("SHAPE of the bound conformation, by ligand size")
print(f"{'heavy atoms':>13} {'n':>4}   {'rod':>5} {'disc':>5} {'sphere':>6}   {'rg/rg_sph':>9} {'max dim':>8}")
for lo, hi in [(0, 20), (20, 30), (30, 40), (40, 55), (55, 999)]:
    g = [r for r in recs if lo <= r["n_heavy"] < hi]
    if len(g) < 3:
        continue
    c = collections.Counter(r["shape"] for r in g)
    n = len(g)
    print(f"{lo:>6}-{hi:<6} {n:>4}   {c['rod']/n:>5.0%} {c['disc']/n:>5.0%} {c['sphere']/n:>6.0%}   "
          f"{st.median(r['rg_norm'] for r in g):>9.2f} {st.median(r['max_dim'] for r in g):>7.1f} A")

multi = [r for r in recs if r["spread_rmsd"] is not None]
print(f"\nCONFORMATIONAL SPREAD within one crystal ({len(multi)} ligands, >1 copy)")
print(f"{'heavy atoms':>13} {'n':>4} {'median':>8} {'p90':>7} {'max':>7}   {'<0.5A':>7} {'>1.5A':>7}")
for lo, hi in [(0, 20), (20, 30), (30, 40), (40, 55), (55, 999)]:
    g = [r["spread_rmsd"] for r in multi if lo <= r["n_heavy"] < hi]
    if len(g) < 3:
        continue
    g = sorted(g)
    print(f"{lo:>6}-{hi:<6} {len(g):>4} {st.median(g):>8.2f} "
          f"{g[int(.9*(len(g)-1))]:>7.2f} {max(g):>7.2f}   "
          f"{sum(v<0.5 for v in g)/len(g):>7.0%} {sum(v>1.5 for v in g)/len(g):>7.0%}")

if len(multi) > 8:
    x = np.array([r["n_heavy"] for r in multi], float)
    y = np.array([r["spread_rmsd"] for r in multi], float)
    b = np.array([r["buried"] for r in multi], float)
    print(f"\ncorr(n_heavy, spread) = {np.corrcoef(x, y)[0,1]:+.3f}")
    print(f"corr(buried,  spread) = {np.corrcoef(b, y)[0,1]:+.3f}")
    nr = [(r["n_rot"], r["spread_rmsd"]) for r in multi
          if r["n_rot"] not in (None, "", "None")]
    if len(nr) > 8:
        a = np.array([float(v) for v, _ in nr]); c = np.array([v for _, v in nr])
        print(f"corr(n_rot,   spread) = {np.corrcoef(a, c)[0,1]:+.3f}   (n={len(nr)})")

print(f"\nwrote {OUT}")
