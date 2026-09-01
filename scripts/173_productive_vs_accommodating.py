#!/usr/bin/env python
r"""
173_productive_vs_accommodating.py -- what separates a pocket that BINDS from one
that merely FITS?

§97b is the sharpest result the project has: azoxystrobin and lufenuron are
mandipropamid-sized, clash just as hard (+1841 and +915 strain against +2066),
and the score nominates clash-relieving substitutions for them just as
confidently -- lufenuron has 16 variants below -100 REU. **Zero of 475 responded
to either compound.** So relieving the clash is necessary and nowhere near
sufficient.

The assay is not a binding assay. It scores LIGAND-DEPENDENT RECRUITMENT OF HAB1,
which requires the gate (85-89) and latch (115-117) to close over the ligand and
present the HAB1 surface. A ligand can occupy the pocket, relieve every clash,
and still never drive that closure -- and nothing in a ddG of the bound complex
would notice.

THE HYPOTHESIS, stated before the numbers: the productive compound contacts the
GATE AND LATCH; the accommodating ones sit in the chamber without touching them.

Measured per compound on its own pose in the common PYR1 frame:
  gate/latch contacts   ligand heavy atoms within 4.5 A of residues 85-89, 115-117
  depth                 ligand centroid to the gate/latch mouth
  buried fraction       ligand atoms with >= 16 protein heavy atoms within 8 A
  W385-proxy            distance to the latch water position (README W385)

⚠ Four compounds carry positives and two do not, so n is small and this is a
descriptive comparison, not a test. It can only suggest which quantity is worth
building a filter on.
"""
import glob, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from importlib import import_module
G = import_module("164_pose_geometry")

AP = os.path.join(ROOT, "data", "agro_params")
GATE = [85, 86, 87, 88, 89]
LATCH = [115, 116, 117]
POS = {"mandipropamid": 20, "fludioxonil": 28, "benzothiadiazole": 22,
       "benoxacor": 10, "azoxystrobin": 0, "lufenuron": 0}
STRAIN = {"mandipropamid": 2066, "fludioxonil": -6, "benzothiadiazole": 6,
          "benoxacor": 36, "azoxystrobin": 1841, "lufenuron": 915}


def frame_atoms():
    """protein heavy atoms by residue, from the common PYR1 frame"""
    by = {}
    for l in open(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")):
        if l.startswith("ATOM") and l[76:78].strip() != "H":
            by.setdefault(int(l[22:26]), []).append(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return {k: np.array(v) for k, v in by.items()}


def main():
    prot = frame_atoms()
    allp = np.vstack(list(prot.values()))
    gate = np.vstack([prot[r] for r in GATE if r in prot])
    latch = np.vstack([prot[r] for r in LATCH if r in prot])
    mouth = np.vstack([prot[r] for r in GATE + LATCH if r in prot]).mean(0)
    meta = json.load(open(os.path.join(AP, "summary.json")))

    def lig_coords(code):
        p = os.path.join(AP, f"{code}_0001.pdb")
        return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                         for l in open(p) if l.startswith(("ATOM", "HETATM"))
                         and l[76:78].strip() != "H"])

    rows = []
    # mandipropamid crystal pose, from the stage1 frame
    xt = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                   for l in open(os.path.join(ROOT, "data", "stage1", "wt_mandi.pdb"))
                   if l.startswith("HETATM") and l[76:78].strip() != "H"])
    cands = [("mandipropamid (crystal)", xt, "mandipropamid")]
    for name, m in meta.items():
        cands.append((f"{name} (predicted)", lig_coords(m["code"]), name))

    print(f"{'compound / pose':<30}{'pos':>5}{'strain':>8}"
          f"{'gate':>6}{'latch':>7}{'g+l':>6}{'depth':>8}{'buried':>8}")
    for label, L, key in cands:
        dg = np.linalg.norm(L[:, None] - gate[None], axis=-1).min(1)
        dl = np.linalg.norm(L[:, None] - latch[None], axis=-1).min(1)
        ng, nl = int((dg < 4.5).sum()), int((dl < 4.5).sum())
        d = np.linalg.norm(L[:, None] - allp[None], axis=-1)
        bur = float(((d < 8).sum(1) >= 16).mean())
        depth = float(np.linalg.norm(L.mean(0) - mouth))
        rows.append(dict(label=label, key=key, pos=POS.get(key, 0),
                         gate=ng, latch=nl, depth=depth, buried=bur))
        print(f"{label:<30}{POS.get(key,0):>5}{STRAIN.get(key,0):>8}"
              f"{ng:>6}{nl:>7}{ng+nl:>6}{depth:>8.2f}{bur:>8.2f}")

    print("\n  gate = ligand heavy atoms within 4.5 A of residues 85-89")
    print("  latch = the same for 115-117; these are the residues whose closure")
    print("  the assay actually reads (ligand-dependent HAB1 recruitment).")
    prod = [r for r in rows if r["pos"] > 0 and "crystal" not in r["label"]]
    dead = [r for r in rows if r["pos"] == 0]
    if prod and dead:
        print(f"\n  compounds WITH responders  (n={len(prod)}): "
              f"gate+latch contacts median {np.median([r['gate']+r['latch'] for r in prod]):.0f}")
        print(f"  compounds WITHOUT responders (n={len(dead)}): "
              f"gate+latch contacts median {np.median([r['gate']+r['latch'] for r in dead]):.0f}")
    json.dump(rows, open(os.path.join(ROOT, "results", "agro_scan",
                                      "productive_vs_accommodating.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
