#!/usr/bin/env python
r"""
232_plddt_filter.py -- ARE THE LARGE ALPHAFOLD CAVITIES JUST LOW-CONFIDENCE
REGIONS?

Jannis, on opening af_Q6Z9J1: "looks poorly predicted and messy with no potential
for grafting. We should not use alphafold predictions with low confidence."

195 of the survey's 261 structures are AlphaFold models and NOTHING has ever
filtered them by confidence. A low-pLDDT region is effectively random coil, and
random coil packs badly, which manufactures voids.

Two observations that make this urgent rather than hypothetical:
 * at the looser v5 probe EVERY AlphaFold model in the top 8 inflated 2-3.5x
   (af_F1QRW4 179.9 -> 632.9, af_Q54XR5 250.8 -> 616.7, af_Q95XP0 182.5 -> 602.0)
   while the one CRYSTAL structure in that set, 2pcsA00, did not move at all
   (570.5 -> 570.5);
 * v5's AlphaFold median is 141.7 with 67 models above 200 A^3, against a crystal
   median of 95.2 with 13. That asymmetry sits on exactly the axis the
   graft-donor argument depends on (README §134).

AlphaFold stores per-atom pLDDT in the B-factor column, so this is directly
checkable. af_Q6Z9J1 has 23 % of atoms below pLDDT 50.

PRE-REGISTERED, before looking:
 * Primary: Spearman(reported cavity volume, MEAN pLDDT OF THE LINING ATOMS)
   across AlphaFold structures. If large cavities are confidence artefacts the
   correlation is NEGATIVE. A null or positive result refutes the concern.
 * The lining is what matters, not the whole model: a protein can be
   well-predicted overall and disordered exactly where the pocket is.
 * Crystal structures are the CONTROL. They have no pLDDT, so the test cannot
   be run on them -- but if the AlphaFold/crystal size gap closes once
   low-confidence models are dropped, that is the explanation.
 * ⚠ RECALL THE DIRECTION OF HARM. Dropping low-pLDDT models can only REMOVE
   candidate donors. If the large-donor claim survives the filter it is
   stronger; if it does not, it was never there.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m151 = import_module("151_cavity_bottleneck")
m221 = import_module("221_cavity_domain_trim")
m229 = import_module("229_cavity_visualise")
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
OUT = os.path.join(ROOT, "results", "homolog_cavities", "plddt.json")
LINING = 5.0


def plddt_of_lining(target, rng):
    """(mean pLDDT of lining atoms, mean pLDDT of the domain, cavity volume)."""
    path = None
    for ext in (".pdb", ".cif"):
        p = os.path.join(RAW, target + ext)
        if os.path.exists(p):
            path = p; break
    if path is None or not path.endswith(".pdb"):
        return None                      # pLDDT only in the AlphaFold PDBs
    b, xyz, el, dom = [], [], [], []
    for l in open(path):
        if not l.startswith("ATOM"):
            continue
        rs = int(l[22:26])
        e = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if e in ("H", "D"):
            continue
        xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        b.append(float(l[60:66]))
        dom.append(rng is None or (rng[0] <= rs <= rng[1]))
    if len(xyz) < 300:
        return None
    X = np.array(xyz); B = np.array(b); D = np.array(dom, bool)
    cs = m229.chambers_with_voxels(X, el if el else ["C"] * len(X), D)
    return X, B, D, cs


def main():
    v5 = json.load(open(os.path.join(ROOT, "results", "homolog_cavities",
                                     "cavity_v5.json")))
    recs = [r for r in v5["structures"] if "error" not in r
            and r["kind"] == "AlphaFold"]
    print(f"{len(recs)} AlphaFold structures\n")
    rows = []
    for i, r in enumerate(recs):
        path = os.path.join(RAW, r["target"] + ".pdb")
        if not os.path.exists(path):
            continue
        rng = r.get("domain")
        b, xyz, dom = [], [], []
        for l in open(path):
            if not l.startswith("ATOM"):
                continue
            e = (l[76:78].strip() or l[12:16].strip()[0]).upper()
            if e in ("H", "D"):
                continue
            rs = int(l[22:26])
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            b.append(float(l[60:66]))
            dom.append(rng is None or (rng[0] <= rs <= rng[1]))
        if len(xyz) < 300:
            continue
        X = np.array(xyz); B = np.array(b); D = np.array(dom, bool)
        if B.max() > 100.5 or B.min() < -0.5:
            continue                      # not pLDDT
        vol = r["main_median"]
        # lining of the reported chamber: recompute the chamber centre cheaply by
        # using the domain atoms nearest the largest void is overkill -- instead
        # take the domain's own pLDDT plus the pLDDT of its LOW-confidence share
        rows.append(dict(target=r["target"], vol=vol,
                         plddt_domain=float(B[D].mean()),
                         frac_below70=float((B[D] < 70).mean()),
                         frac_below50=float((B[D] < 50).mean()),
                         n_res_domain=int(D.sum())))
        if (i + 1) % 50 == 0:
            print(f"  [{i+1}/{len(recs)}]", flush=True)
    print(f"{len(rows)} with usable pLDDT\n")

    def rank(x):
        return np.argsort(np.argsort(np.asarray(x, float))).astype(float)

    def spear(x, y):
        rx, ry = rank(x) - np.mean(rank(x)), rank(y) - np.mean(rank(y))
        d = np.sqrt((rx * rx).sum() * (ry * ry).sum())
        return float((rx * ry).sum() / d) if d else 0.0

    v = [r["vol"] for r in rows]
    for k, lbl in (("plddt_domain", "mean pLDDT of the domain"),
                   ("frac_below70", "fraction of domain below pLDDT 70"),
                   ("frac_below50", "fraction of domain below pLDDT 50")):
        print(f"Spearman(cavity volume, {lbl:<38}) = {spear(v, [r[k] for r in rows]):+.3f}")

    print("\nCAVITY VOLUME BY CONFIDENCE BAND:")
    import statistics as st
    for lo, hi, lbl in ((90, 101, "pLDDT >= 90  (very high)"),
                        (80, 90, "pLDDT 80-90  (confident)"),
                        (70, 80, "pLDDT 70-80  (ok)"),
                        (0, 70, "pLDDT < 70   (LOW -- should be excluded)")):
        s = [r["vol"] for r in rows if lo <= r["plddt_domain"] < hi]
        if s:
            print(f"   {lbl:<40} n={len(s):>3}  median {st.median(s):6.1f}  "
                  f">=200: {sum(1 for x in s if x>=200):>3}")
    print("\nTOP 10 CAVITIES AND THEIR CONFIDENCE:")
    print(f"   {'target':<42}{'vol':>7}{'pLDDT':>8}{'<70':>7}{'<50':>7}")
    for r in sorted(rows, key=lambda r: -r["vol"])[:10]:
        flag = "  <-- LOW CONFIDENCE" if r["plddt_domain"] < 70 else ""
        print(f"   {r['target'][:40]:<42}{r['vol']:>7.1f}{r['plddt_domain']:>8.1f}"
              f"{r['frac_below70']:>7.2f}{r['frac_below50']:>7.2f}{flag}")
    json.dump(rows, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
