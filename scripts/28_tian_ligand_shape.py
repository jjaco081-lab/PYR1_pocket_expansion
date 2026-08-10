#!/usr/bin/env python
"""
28_tian_ligand_shape.py -- LONGER tunnel or BIGGER chamber? Ask the screening data.

The question (user, 2026-08-07)
-------------------------------
Script 27 established that Tian's hits stop at 15.3 A longest dimension while
compounds up to 29.1 A were screened -- i.e. the existing landscape has a LENGTH
ceiling. But that only matters for design if the large compounds that FAILED were
elongated. If the large failures were mostly globular, then a longer tunnel is
the wrong build and we should be widening/deepening instead.

So: among Tian's compounds, does increasing size arrive as LENGTH (rods) or as
BULK (spheres/discs), and which shape class did the libraries struggle with?

Method
------
3D conformer per compound (ETKDG + MMFF), then the standard PMI shape descriptors:

    NPR1 = I1/I3, NPR2 = I2/I3   (I1 <= I2 <= I3, normalised principal moments)

    rod    -> NPR1 small, NPR2 large   (canonical corner ~ (0, 1))
    disc   -> ~ (0.5, 0.5)
    sphere -> ~ (1, 1)

Classified as: ROD if NPR1 < 0.35 and NPR2 > 0.70; SPHERE if NPR1 > 0.55;
otherwise DISC/intermediate. Also reports length per heavy atom, a simple
intuition check that needs no conventions.

The decisive number is HIT RATE within each (size bin x shape class) cell. A
library that can bind bulky-but-large things and not long thin things will show
it there, and that is what decides tunnel vs chamber.

Run with the dockenv python (openpyxl + rdkit).
"""
import os, warnings, collections
warnings.filterwarnings("ignore")
import numpy as np
import openpyxl
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Descriptors3D
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
XLSX = os.path.join(SD, "pnas.2519924122.sd01.xlsx")
OUT = os.path.join(ROOT, "results", "28_tian_ligand_shape.csv")

LADDER = {
    "ABA": "CC1(C)C(=CC(=O)CC1(O)/C=C/C(C)=C/C(O)=O)C",
    "crocetin": "CC(=CC=CC=C(C)C=CC=C(C)C=CC(=O)O)C=CC(=O)O",
    "b-apo-8'-carotenal": "CC1=C(C(CCC1)(C)C)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=O",
    "b-carotene (C40)": "CC1=C(C(CCC1)(C)C)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C2=C(CCCC2(C)C)C",
}


def shape(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    nh, mw = m.GetNumHeavyAtoms(), Descriptors.MolWt(m)
    mh = Chem.AddHs(m)
    if AllChem.EmbedMolecule(mh, randomSeed=0xC0FFEE, useRandomCoords=True) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=400)
    except Exception:
        pass
    m3 = Chem.RemoveHs(mh)
    X = m3.GetConformer().GetPositions()
    length = float(np.linalg.norm(X[:, None] - X[None], axis=-1).max())
    try:
        npr1, npr2 = Descriptors3D.NPR1(m3), Descriptors3D.NPR2(m3)
        asph = Descriptors3D.Asphericity(m3)
    except Exception:
        return None
    return nh, mw, length, float(npr1), float(npr2), float(asph)


def klass(npr1, npr2):
    if npr1 < 0.35 and npr2 > 0.70:
        return "ROD"
    if npr1 > 0.55:
        return "SPHERE"
    return "DISC"


print("=" * 84)
print("REFERENCE LADDER")
print("=" * 84)
print(f"{'ligand':<22}{'heavy':>6}{'len A':>8}{'NPR1':>7}{'NPR2':>7}{'class':>8}{'len/heavy':>11}")
for n, s in LADDER.items():
    r = shape(s)
    if r:
        print(f"{n:<22}{r[0]:>6}{r[2]:>8.1f}{r[3]:>7.2f}{r[4]:>7.2f}"
              f"{klass(r[3], r[4]):>8}{r[2]/r[0]:>11.2f}")

wb = openpyxl.load_workbook(XLSX, read_only=True)
ws = wb.worksheets[0]
it = ws.iter_rows(values_only=True)
hdr = [str(h) for h in next(it)]
i_nm, i_hit = hdr.index("Library_name"), hdr.index("Hit?")
i_smi = hdr.index("canonical_smiles")
rows = [r for r in it]
wb.close()

recs = []
for r in rows:
    smi = r[i_smi]
    if not smi or not isinstance(smi, str):
        continue
    v = shape(smi)
    if v is None:
        continue
    nh, mw, ln, n1, n2, asph = v
    recs.append(dict(name=str(r[i_nm]), hit=str(r[i_hit]).strip().lower() == "yes",
                     heavy=nh, mw=mw, length=ln, npr1=n1, npr2=n2,
                     asph=asph, cls=klass(n1, n2)))

with open(OUT, "w") as fh:
    fh.write("name,hit,heavy_atoms,mw,length_A,npr1,npr2,asphericity,shape_class\n")
    for d in recs:
        fh.write(f"\"{d['name']}\",{int(d['hit'])},{d['heavy']},{d['mw']:.1f},"
                 f"{d['length']:.2f},{d['npr1']:.3f},{d['npr2']:.3f},"
                 f"{d['asph']:.3f},{d['cls']}\n")

print(f"\nparsed {len(recs)} compounds, {sum(d['hit'] for d in recs)} hits")

print("\n" + "=" * 84)
print("DOES SIZE ARRIVE AS LENGTH OR AS BULK?")
print("=" * 84)
BINS = [(0, 20), (20, 25), (25, 30), (30, 40), (40, 200)]
print(f"{'heavy atoms':<14}{'n':>6}{'med len':>9}{'len/heavy':>11}   shape mix (ROD/DISC/SPHERE)")
for lo, hi in BINS:
    sel = [d for d in recs if lo <= d["heavy"] < hi]
    if not sel:
        continue
    c = collections.Counter(d["cls"] for d in sel)
    n = len(sel)
    print(f"{f'{lo}-{hi}':<14}{n:>6}{np.median([d['length'] for d in sel]):>9.1f}"
          f"{np.median([d['length']/d['heavy'] for d in sel]):>11.2f}   "
          f"{100*c['ROD']//n:>3}% / {100*c['DISC']//n:>3}% / {100*c['SPHERE']//n:>3}%")

print("\n" + "=" * 84)
print("HIT RATE BY SIZE x SHAPE  -- the number that decides tunnel vs chamber")
print("=" * 84)
print(f"{'heavy atoms':<14}" + "".join(f"{c:>18}" for c in ("ROD", "DISC", "SPHERE")))
for lo, hi in BINS:
    line = f"{f'{lo}-{hi}':<14}"
    for c in ("ROD", "DISC", "SPHERE"):
        sel = [d for d in recs if lo <= d["heavy"] < hi and d["cls"] == c]
        if not sel:
            line += f"{'--':>18}"
        else:
            h = sum(d["hit"] for d in sel)
            line += f"{f'{h}/{len(sel)} ({100*h/len(sel):.1f}%)':>18}"
    print(line)

print("\n" + "=" * 84)
print("LENGTH, split by outcome")
print("=" * 84)
for lo, hi in BINS:
    sel = [d for d in recs if lo <= d["heavy"] < hi]
    hits = [d["length"] for d in sel if d["hit"]]
    miss = [d["length"] for d in sel if not d["hit"]]
    if hits and miss:
        print(f"  {lo}-{hi} heavy: HITS median {np.median(hits):5.1f} A (max {max(hits):5.1f}) | "
              f"MISSES median {np.median(miss):5.1f} A (max {max(miss):5.1f})")

big = [d for d in recs if d["heavy"] >= 30]
if big:
    print("\nLARGE COMPOUNDS (>=30 heavy atoms) that FAILED, longest first:")
    for d in sorted([d for d in big if not d["hit"]], key=lambda x: -x["length"])[:12]:
        print(f"   len {d['length']:5.1f} A  {d['heavy']:>3} heavy  {d['cls']:<6}  {d['name']}")

print(f"\nwrote {OUT}")
