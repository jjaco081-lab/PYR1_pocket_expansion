#!/usr/bin/env python
"""
27_tian_ligand_size.py -- did the published PYR1 engineering work ever test
ligands LARGER than the pocket we are trying to build?

Why this matters
----------------
The project's win condition (set by the user, 2026-08-07) is not just "open the
pocket". It is: produce a receptor that could not have come out of the existing
libraries. Five of our six wall positions (K59, E94, F108, Y120, E141) are in
Mosquna 2011's 39-residue site-saturation set AND in Tian's randomised set; only
R79 is in neither. So mutating those five is re-treading known ground -- UNLESS
the ligands they were screened against were all small, in which case "same
position, genuinely larger ligand" is new territory.

This measures the SIZE of Tian's screening library and of its hits.

MW IS THE WRONG METRIC and this script does not rely on it. Tian's largest hit by
mass is Tetrac at 747.8 Da, but four iodines contribute ~508 Da of that; it is a
small molecule wearing heavy atoms. What determines whether something fits a
pocket is heavy-atom count and longest dimension, so both are computed from the
deposited SMILES with RDKit (ETKDG embed + MMFF).

Outputs results/27_tian_ligand_size.csv and a summary table.

Run with the dockenv python (has openpyxl AND rdkit).
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import openpyxl
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
RDLogger.DisableLog("rdApp.*")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
XLSX = os.path.join(SD, "pnas.2519924122.sd01.xlsx")
OUT = os.path.join(ROOT, "results", "27_tian_ligand_size.csv")

# our ligand ladder, for comparison on the identical metric
LADDER = {
    "ABA (the native ligand)": "CC1(C)C(=CC(=O)CC1(O)/C=C/C(C)=C/C(O)=O)C",
    "retinoic acid": "CC1=C(C(CCC1)(C)C)/C=C/C(C)=C/C=C/C(C)=C/C(O)=O",
    "crocetin": "CC(=CC=CC=C(C)C=CC=C(C)C=CC(=O)O)C=CC(=O)O",
    "beta-apo-8'-carotenal": "CC1=C(C(CCC1)(C)C)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=O",
    "beta-carotene (C40)": "CC1=C(C(CCC1)(C)C)/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C2=C(CCCC2(C)C)C",
}


def geom(smiles):
    """heavy atoms, MW, and longest interatomic distance of a 3D conformer."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    n_heavy = m.GetNumHeavyAtoms()
    mw = Descriptors.MolWt(m)
    mh = Chem.AddHs(m)
    if AllChem.EmbedMolecule(mh, randomSeed=0xC0FFEE, useRandomCoords=True) != 0:
        return n_heavy, mw, None
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=400)
    except Exception:
        pass
    X = Chem.RemoveHs(mh).GetConformer().GetPositions()
    d = float(np.linalg.norm(X[:, None] - X[None], axis=-1).max())
    return n_heavy, mw, d


print("=" * 78)
print("OUR LIGAND LADDER on this metric")
print("=" * 78)
print(f"{'ligand':<28}{'heavy':>7}{'MW':>9}{'length A':>10}")
ref = {}
for name, smi in LADDER.items():
    g = geom(smi)
    ref[name] = g
    print(f"{name:<28}{g[0]:>7}{g[1]:>9.1f}"
          + (f"{g[2]:>10.1f}" if g[2] else f"{'--':>10}"))

ABA_HEAVY = ref["ABA (the native ligand)"][0]

print("\n" + "=" * 78)
print("TIAN SCREENING LIBRARY (sd01)")
print("=" * 78)
wb = openpyxl.load_workbook(XLSX, read_only=True)
ws = wb.worksheets[0]
it = ws.iter_rows(values_only=True)
hdr = [str(h) for h in next(it)]
i_nm, i_hit = hdr.index("Library_name"), hdr.index("Hit?")
i_smi, i_mw = hdr.index("canonical_smiles"), hdr.index("mw_parent")
rows = [r for r in it]
wb.close()

recs = []
for r in rows:
    smi = r[i_smi]
    if not smi or not isinstance(smi, str):
        continue
    g = geom(smi)
    if g is None:
        continue
    recs.append(dict(name=str(r[i_nm]), hit=str(r[i_hit]).strip().lower() == "yes",
                     heavy=g[0], mw=g[1], length=g[2]))

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as fh:
    fh.write("name,hit,heavy_atoms,mw,length_A\n")
    for d in recs:
        fh.write(f"\"{d['name']}\",{int(d['hit'])},{d['heavy']},{d['mw']:.1f},"
                 f"{d['length'] if d['length'] else ''}\n")

allr = np.array([d["heavy"] for d in recs])
hitr = np.array([d["heavy"] for d in recs if d["hit"]])
print(f"parsed {len(recs)} compounds, {len(hitr)} hits")
print(f"\nHEAVY ATOMS      ABA = {ABA_HEAVY}")
print(f"  screened: median {np.median(allr):.0f}, 90th {np.percentile(allr,90):.0f}, "
      f"99th {np.percentile(allr,99):.0f}, max {allr.max():.0f}")
print(f"  HITS    : median {np.median(hitr):.0f}, 90th {np.percentile(hitr,90):.0f}, "
      f"max {hitr.max():.0f}")

print(f"\n{'threshold':<34}{'screened':>10}{'HITS':>7}")
for name in ("crocetin", "beta-apo-8'-carotenal", "beta-carotene (C40)"):
    n = ref[name][0]
    print(f"heavy atoms >= {n:<3} ({name:<22}){int((allr>=n).sum()):>10}"
          f"{int((hitr>=n).sum()):>7}")

L = [d for d in recs if d["length"]]
Lh = [d for d in L if d["hit"]]
if Lh:
    la = np.array([d["length"] for d in L]); lh = np.array([d["length"] for d in Lh])
    print(f"\nLONGEST DIMENSION (A)   ABA = {ref['ABA (the native ligand)'][2]:.1f}")
    print(f"  screened: median {np.median(la):.1f}, 90th {np.percentile(la,90):.1f}, max {la.max():.1f}")
    print(f"  HITS    : median {np.median(lh):.1f}, 90th {np.percentile(lh,90):.1f}, max {lh.max():.1f}")
    bc = ref["beta-carotene (C40)"][2]
    print(f"  hits longer than beta-carotene ({bc:.1f} A): {int((lh>bc).sum())}")

print("\nLARGEST HITS BY HEAVY-ATOM COUNT:")
for d in sorted([d for d in recs if d["hit"]], key=lambda x: -x["heavy"])[:12]:
    ln = f"{d['length']:.1f}" if d["length"] else "--"
    print(f"   {d['heavy']:>3} heavy  MW {d['mw']:7.1f}  len {ln:>5}  {d['name']}")

print(f"\nwrote {OUT}")
