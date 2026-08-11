#!/usr/bin/env python
"""
44_floppy_ligand_survey.py -- do protein receptors actually engulf large flexible
ligands, or do they grip a head group and let the tail hang out?

WHY THIS EXISTS
---------------
This project's premise is that PYR1's pocket can be enlarged to admit ligands far
bigger than ABA (20 heavy atoms). The measured landscape wall sits near 40 heavy
atoms / 15.3 A, and the 229-compound probe set is 217/229 rods at 15-22 A. The
standing objection is thermodynamic rather than steric: ordering a floppy rod
inside a channel costs conformational entropy that no scoring function in this
pipeline models. That objection is not directly computable here.

It is, however, answerable empirically as a database question. If receptors across
the PDB systematically bury only PART of a large flexible ligand -- the head -- and
leave the tail solvent-exposed, then full engulfment is not how nature solves this,
and the design objective should change from "enclose the ligand" to "grip a head
group and tolerate an exit vector." That is a materially different target, and it
would reframe the pocket-expansion goal.

WHAT IS MEASURED
----------------
For each sampled ligand, in one representative complex:
  * n_heavy, n_rotatable        (RDKit, from the RCSB canonical SMILES)
  * buried fraction   = 1 - SASA(ligand in complex) / SASA(ligand alone)
  * AXIAL BURIAL PROFILE -- the decisive measurement. Ligand atoms are projected
    onto their first principal axis and split into three equal-length bins. If
    burial is uniform the three agree; if receptors grip a head and expose a tail,
    the terminal bins diverge. Bins are ordered so `head` is the more buried end,
    so the reported gap is head-minus-tail by construction and its SIGN carries no
    information -- only its magnitude does.

CONTROLS AND LIMITS
-------------------
  * Only protein chains are kept. Other ligands, cofactors and waters are removed,
    so burial reported here is burial BY PROTEIN, not by everything nearby. A
    ligand stacked against a cofactor will read as less buried than it feels.
  * One entry per ligand, no resolution filter, and crystal-packing neighbours are
    not modelled -- a tail that looks exposed may contact a symmetry mate.
  * Selection is biased toward ligands that crystallised at all, which under-counts
    exactly the floppy chemistry in question. This biases the survey TOWARD finding
    engulfment, so a negative result (tails exposed) is the conservative one.

Usage:  python 44_floppy_ligand_survey.py [--n-per-bin 25] [--out results/44_floppy_survey.csv]
"""
import argparse, csv, io, json, os, sys, time, urllib.request

import numpy as np
from Bio.PDB import MMCIFParser, Structure, Model, Chain
from Bio.PDB.SASA import ShrakeRupley

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "data", "survey_cache")
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
CHEMCOMP = "https://data.rcsb.org/rest/v1/core/chemcomp/"

# formula-weight bins spanning ABA-sized to well past the probe set
MW_BINS = [(250, 320), (320, 400), (400, 500), (500, 620), (620, 800), (800, 1100)]

ap = argparse.ArgumentParser()
ap.add_argument("--n-per-bin", type=int, default=25)
ap.add_argument("--out", default=os.path.join(ROOT, "results", "44_floppy_survey.csv"))
args = ap.parse_args()
os.makedirs(CACHE, exist_ok=True)
os.makedirs(os.path.dirname(args.out), exist_ok=True)


def post(url, payload, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(
                url, data=json.dumps(payload).encode(),
                headers={"Content-Type": "application/json"})
            return json.load(urllib.request.urlopen(req, timeout=30))
        except Exception:
            time.sleep(2 * (i + 1))
    return {}


def get(url, tries=3):
    for i in range(tries):
        try:
            return json.load(urllib.request.urlopen(url, timeout=30))
        except Exception:
            time.sleep(2 * (i + 1))
    return {}


def comps_in_bin(lo, hi, rows):
    q = {"query": {"type": "terminal", "service": "text_chem",
                   "parameters": {"attribute": "chem_comp.formula_weight",
                                  "operator": "range",
                                  "value": {"from": lo, "to": hi}}},
         "return_type": "mol_definition",
         "request_options": {"paginate": {"start": 0, "rows": rows}}}
    return [r["identifier"] for r in post(SEARCH, q).get("result_set", [])]


def entries_for(comp):
    q = {"query": {"type": "terminal", "service": "text_chem",
                   "parameters": {"attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                                  "operator": "exact_match", "value": comp}},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": 3}}}
    return [r["identifier"] for r in post(SEARCH, q).get("result_set", [])]


def smiles_of(comp):
    d = get(CHEMCOMP + comp)
    for k in ("rcsb_chem_comp_descriptor", "pdbx_chem_comp_descriptor"):
        v = d.get(k)
        if isinstance(v, dict) and v.get("smiles"):
            return v["smiles"]
        if isinstance(v, list):
            for e in v:
                if e.get("type", "").upper() == "SMILES_CANONICAL":
                    return e.get("descriptor")
    return None


def fetch_cif(pdb):
    p = os.path.join(CACHE, f"{pdb}.cif")
    if not os.path.exists(p):
        try:
            with urllib.request.urlopen(
                    f"https://files.rcsb.org/download/{pdb}.cif", timeout=60) as r:
                open(p, "wb").write(r.read())
        except Exception:
            return None
    return p if os.path.getsize(p) > 5000 else None


AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL MSE".split())
sr = ShrakeRupley()


def analyse(cif, comp):
    """buried fraction overall and per axial third, for the first copy of `comp`."""
    try:
        model = MMCIFParser(QUIET=True).get_structure("x", cif)[0]
    except Exception:
        return None
    lig = prot = None
    ligres, protres = None, []
    for ch in model:
        for r in ch:
            if r.get_resname() == comp and ligres is None:
                ligres = r
            elif r.get_resname() in AA:
                protres.append(r)
    if ligres is None or len(protres) < 30:
        return None
    if len(protres) > 2500:
        return None      # Shrake-Rupley cost; skip ribosome-scale assemblies
    if len([a for a in ligres if a.element != "H"]) < 12:
        return None

    def build(residues):
        """Fresh Structure holding heavy-atom copies; ligand is always residue 1.

        ShrakeRupley only accepts real Biopython entities, and the residues here
        come from different chains with colliding ids, so they are copied into one
        chain and renumbered.
        """
        s = Structure.Structure("s"); m = Model.Model(0); c = Chain.Chain("A")
        s.add(m); m.add(c)
        for i, r in enumerate(residues):
            rc = r.copy()
            rc.id = (rc.id[0], i + 1, " ")
            for a in [x for x in rc if x.element == "H"]:
                rc.detach_child(a.get_id())
            if len(rc):
                c.add(rc)
        return s, c

    s_free, c_free = build([ligres])
    sr.compute(s_free, level="A")
    free_res = list(c_free)[0]
    names = [a.get_id() for a in free_res]
    free = np.array([a.sasa for a in free_res])

    s_cx, c_cx = build([ligres] + protres)
    sr.compute(s_cx, level="A")
    cx_res = list(c_cx)[0]
    bound = np.array([cx_res[n].sasa for n in names])
    latoms = [free_res[n] for n in names]

    tot_free = free.sum()
    if tot_free <= 0:
        return None
    buried = 1.0 - bound.sum() / tot_free

    # axial profile
    xyz = np.array([a.coord for a in latoms], dtype=float)
    c = xyz - xyz.mean(0)
    axis = np.linalg.svd(c, full_matrices=False)[2][0]
    t = c @ axis
    length = float(t.max() - t.min())
    edges = np.linspace(t.min(), t.max() + 1e-6, 4)
    thirds = []
    for i in range(3):
        m = (t >= edges[i]) & (t < edges[i + 1])
        if m.sum() == 0 or free[m].sum() <= 0:
            thirds.append(np.nan)
        else:
            thirds.append(1.0 - bound[m].sum() / free[m].sum())
    ends = [thirds[0], thirds[2]]
    if not np.isnan(ends).any() and ends[0] < ends[1]:
        thirds = thirds[::-1]          # orient so index 0 is the buried "head"

    # Permutation null. Orienting the bins by which end is more buried makes the
    # head-tail gap positive by construction, so its magnitude is biased upward
    # even for a ligand buried perfectly uniformly. Recomputing the same statistic
    # on atoms shuffled between bins -- destroying axial structure but keeping bin
    # sizes and the max-orientation -- measures that bias directly. Only the gap
    # in EXCESS of null_gap is evidence of real head/tail asymmetry.
    rng = np.random.default_rng(0)
    idx = np.arange(len(t))
    sizes = [int(((t >= edges[i]) & (t < edges[i + 1])).sum()) for i in range(3)]
    nulls = []
    for _ in range(40):
        rng.shuffle(idx)
        g, o = [], 0
        for n in sizes:
            sel = idx[o:o + n]; o += n
            g.append(1.0 - bound[sel].sum() / free[sel].sum()
                     if n and free[sel].sum() > 0 else np.nan)
        if not np.isnan([g[0], g[2]]).any():
            nulls.append(abs(g[0] - g[2]))
    return dict(n_heavy=len(latoms), buried=buried, length=length,
                head=thirds[0], mid=thirds[1], tail=thirds[2],
                null_gap=float(np.mean(nulls)) if nulls else np.nan)


def rot_bonds(smi):
    try:
        sys.path.insert(0, "")
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None, None
        return rdMolDescriptors.CalcNumRotatableBonds(m), Descriptors.MolWt(m)
    except Exception:
        return None, None


rows = []
seen = set()
for lo, hi in MW_BINS:
    comps = comps_in_bin(lo, hi, args.n_per_bin * 6)
    kept = 0
    for comp in comps:
        if kept >= args.n_per_bin or comp in seen:
            continue
        seen.add(comp)
        ents = entries_for(comp)
        if not ents:
            continue
        smi = smiles_of(comp)
        nrot, mw = rot_bonds(smi) if smi else (None, None)
        cif = fetch_cif(ents[0])
        if not cif:
            continue
        a = analyse(cif, comp)
        if not a:
            continue
        a.update(comp=comp, entry=ents[0], mw_bin=f"{lo}-{hi}", n_rot=nrot, mw=mw)
        rows.append(a)
        kept += 1
    print(f"bin {lo}-{hi} Da: {kept} ligands analysed", flush=True)

if not rows:
    sys.exit("no ligands analysed")

cols = ["comp", "entry", "mw_bin", "mw", "n_heavy", "n_rot", "length",
        "buried", "head", "mid", "tail", "null_gap"]
with open(args.out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow(r)

# ---- summary ----
import statistics as st
print(f"\n{len(rows)} ligands\n")
print(f"{'heavy atoms':>14} {'n':>4} {'buried':>8} {'head':>8} {'tail':>8} {'gap':>8} {'null':>7} {'excess':>8}")
for lo, hi in [(0, 20), (20, 30), (30, 40), (40, 55), (55, 999)]:
    g = [r for r in rows if lo <= r["n_heavy"] < hi
         and not any(np.isnan([r["head"], r["tail"]]))]
    if len(g) < 3:
        continue
    print(f"{lo:>7}-{hi:<6} {len(g):>4} "
          f"{st.median(r['buried'] for r in g):>8.2f} "
          f"{st.median(r['head'] for r in g):>8.2f} "
          f"{st.median(r['tail'] for r in g):>8.2f} "
          f"{st.median(r['head']-r['tail'] for r in g):>8.2f} "
          f"{st.median(r['null_gap'] for r in g):>7.2f} "
          f"{st.median(r['head']-r['tail']-r['null_gap'] for r in g):>8.2f}")

ok = [r for r in rows if not any(np.isnan([r["head"], r["tail"]]))]
if len(ok) > 8:
    x = np.array([r["n_heavy"] for r in ok], float)
    for lab in ("buried", "head", "tail"):
        y = np.array([r[lab] for r in ok], float)
        print(f"\ncorr(n_heavy, {lab}) = {np.corrcoef(x, y)[0,1]:+.3f}")
    gap = np.array([r["head"] - r["tail"] - r["null_gap"] for r in ok], float)
    print(f"corr(n_heavy, EXCESS head-tail gap over null) = {np.corrcoef(x, gap)[0,1]:+.3f}")
    print("\nINTERPRETATION")
    print("  buried flat vs size, small gap    -> receptors DO engulf; entropy objection weakens")
    print("  buried falls, gap grows with size -> head-gripping; change the design objective")
print(f"\nwrote {args.out}")
