#!/usr/bin/env python
r"""
181_start_domains.py -- the >1000 A^3 START lipid-transfer domains, which §7
named as the attractive graft donors and which were NEVER SURVEYED.

§104's check: STARD1, MLN64/STARD3, CERT/STARD11, PCTP and relatives are all
absent from the 261-domain foldseek scan behind §56. That scan topped out at
634 A^3, and only 66 of its 261 entries are experimental structures -- its two
largest are AlphaFold models, which §79 excluded for carrying no ligand and no
crystal pocket. So the graft direction was closed on a candidate set that never
included the biggest members of the fold.

These are fetched by PDB ID rather than by foldseek, because the point is to
test named structures whose cavities are documented in the literature, not to
re-run a search that has already failed to retrieve them. Foldseek is then run
PYR1-vs-these to get identity, alnTM and the alignment needed for the §79
machinery-coverage question.

Measured identically to §86c so the numbers are comparable:
  total       largest enclosed component (script 08 criterion)
  main        the sub-volume a 2.4 A sphere can reach -- what a ligand can use
  r_max       largest inscribed sphere
PYR1 is 164 total / 119 main / 3.21; 2PCS, the best donor so far, is 570/570/3.80.

⚠ A big cavity is necessary and not sufficient. §56a's verdict on the whole graft
direction was "rebuild, not transplant", and §7's objection to any rebuild is
that the gate and latch are a two-state switch. Nothing here tests the switch.
"""
import os, subprocess, sys, urllib.request
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m151 = import_module("151_cavity_bottleneck")

OUT = os.path.join(ROOT, "data", "start_domains")
FOLDSEEK = "/bigdata/cutlerlab/jjaco081/tools/foldseek/foldseek/bin/foldseek"

#: named START/StARkin domains with documented lipid cavities, plus their chain
TARGETS = {
    "3P0L": ("A", "STARD1 / StAR, cholesterol transfer"),
    "1EM2": ("A", "MLN64 / STARD3 START domain, cholesterol"),
    "2E3M": ("A", "CERT / STARD11 START domain, ceramide"),
    "2E3P": ("A", "CERT / STARD11 with ceramide bound"),
    "1LN1": ("A", "PCTP / STARD2, phosphatidylcholine"),
    "1JSS": ("A", "STARD4-like / GLTP fold check"),
    "2R55": ("A", "STARD5, bile acid"),
    "3QZT": ("A", "STARD13-family START domain"),
}


def fetch(pdb):
    p = os.path.join(OUT, f"{pdb}.cif")
    if os.path.exists(p) and os.path.getsize(p) > 2000:
        return p
    try:
        req = urllib.request.Request(f"https://files.rcsb.org/download/{pdb}.cif",
                                     headers={"User-Agent": "PYR1-expansion/1.0"})
        open(p, "wb").write(urllib.request.urlopen(req, timeout=60).read())
        return p
    except Exception as e:
        print(f"  {pdb}: fetch failed ({e})")
        return None


def chain_atoms(path, chain):
    """protein heavy atoms of one chain -> (xyz, elements, n_res)"""
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
    xyz, el, res = [], [], set()
    for f in rows:
        if f[cols["group_PDB"]] != "ATOM" or f[cols["auth_asym_id"]] != chain:
            continue
        if f[cols["type_symbol"]] == "H":
            continue
        alt = f[cols["label_alt_id"]] if "label_alt_id" in cols else "."
        if alt not in (".", "?", "A"):
            continue
        xyz.append([float(f[cols[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")])
        el.append(f[cols["type_symbol"]])
        res.add(f[cols["auth_seq_id"]])
    return np.array(xyz), el, len(res)


def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"{'PDB':<6}{'chain':>6}{'res':>5}{'total':>8}{'main':>8}{'r_max':>7}   description")
    rows = []
    for pdb, (ch, desc) in TARGETS.items():
        p = fetch(pdb)
        if not p:
            continue
        xyz, el, nres = chain_atoms(p, ch)
        if len(xyz) < 500:
            print(f"{pdb:<6}{ch:>6}{nres:>5}   too few atoms ({len(xyz)}) -- chain id may be wrong")
            continue
        r = m151.analyse(pdb, xyz, el)
        if r is None:
            print(f"{pdb:<6}{ch:>6}{nres:>5}   no enclosed cavity found")
            continue
        print(f"{pdb:<6}{ch:>6}{nres:>5}{r['total']:>8.0f}{r.get('main',0):>8.0f}"
              f"{r.get('r_max',0):>7.2f}   {desc}")
        rows.append((pdb, ch, nres, r, desc))
    print(f"\n  reference: PYR1 164 total / 119 main / 3.21 ; 2PCS (best donor so far) "
          f"570 / 570 / 3.80")
    print(f"\n  structures written to {OUT}/ -- inspect before anything is built on them")
    return 0


if __name__ == "__main__":
    sys.exit(main())
