#!/usr/bin/env python
r"""
120_start_ligand_census.py -- how many START/SRPBCC-fold structures have a
genuine POCKET ligand?

WHY THE NAIVE COUNT IS USELESS
Counting structures with any HETATM returns 478 of 478 (100 %), because the field
is dominated by things that are not pocket ligands: NAG/BMA/MAN/FUC glycans on the
surface, TYS/CSO/TPO modified residues, ACE caps, and cryoprotectants. A ligand
census has to be geometric, not compositional.

WHAT COUNTS HERE
A HETATM component is a pocket ligand if it
  * is not water, an ion, a glycan, a modified residue or a common cryo additive
  * has >= 8 heavy atoms (below that it is almost always a buffer fragment)
  * is BURIED: >= 70 % of its heavy atoms have a protein heavy atom within 4.5 A,
    and it makes >= 20 such contacts in total

Burial is what separates a ligand sitting in the helix-grip cavity from a sugar
stuck to the outside, and it is the property that matters for the question being
asked -- whether the fold has characterised cavity chemistry to learn from.

⚠ This counts STRUCTURES and COMPONENTS, not distinct chemistry. Several entries
of the same protein with the same ligand inflate a structure count, so distinct
(component, fold-family) pairs are reported alongside.
"""
import glob
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#: ⚠ THE RIGHT DIRECTORY MATTERS AND THE OBVIOUS ONE IS WRONG.
#: data/survey_cache/ holds 478 CIFs, but it was populated by scripts 44-46 for the
#: FLOPPY-LIGAND survey -- structures chosen because they contain a particular
#: flexible ligand, not because they share PYR1's fold. Running the census there
#: returns 97 % with a buried ligand and puts bacteriochlorophyll, chlorophyll,
#: carotenoid and lauryl maltoside at the top, i.e. it is measuring
#: light-harvesting complexes. The actual START/SRPBCC foldseek hit set is
#: results/homolog_cavities/raw/.
CACHE = os.path.join(ROOT, "results", "homolog_cavities", "raw")
OUT = os.path.join(ROOT, "results", "start_census")

AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
      "MSE", "SEC", "PYL"}
#: things that are never the pocket ligand of interest
GLYCAN = {"NAG", "BMA", "MAN", "FUC", "GAL", "GLC", "XYS", "SIA", "NDG", "BGC",
          "A2G", "GLA", "RAM", "XYP", "FUL", "NGA"}
MODRES = {"TYS", "CSO", "TPO", "SEP", "PTR", "MLY", "KCX", "CME", "OCS", "CSD",
          "LLP", "PCA", "FME", "ACE", "NH2", "SAH", "NLE", "ABA", "AIB"}
JUNK = {"HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "PGE",
        "MPD", "ACT", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "FE", "FE2",
        "NI", "CD", "CO", "CU", "CU1", "IOD", "BR", "NO3", "TRS", "IMD", "DMS",
        "FMT", "EPE", "MES", "CIT", "FLC", "TAR", "ACY", "BME", "NH4", "UNX",
        "DTT", "SCN", "AZI", "BCT", "SR", "CS", "RB", "HG", "AU", "PT", "LI",
        "BA", "OXY", "PEO", "P6G", "1PE", "2PE", "12P", "15P", "XPE", "MRD",
        "BU3", "IPA", "ETA", "NHE", "SIN", "MLI", "MLA", "UNL", "UNK", "SO3",
        "MOH", "EOH", "ACN", "BCN", "DIO", "TLA", "MAE", "CAC", "BEZ", "PYR"}
MIN_HEAVY = 8
BURIED_FRAC = 0.70
MIN_CONTACTS = 20


def parse_cif(path):
    """Minimal mmCIF atom_site reader -- column positions are read from the loop
    header rather than assumed, because the cached files do not all use the same
    column order."""
    prot, het = [], defaultdict(list)
    with open(path) as fh:
        cols, inloop = {}, False
        for line in fh:
            if line.startswith("_atom_site."):
                cols[line.strip().split(".", 1)[1]] = len(cols)
                inloop = True
                continue
            if inloop and (line.startswith("#") or line.startswith("loop_")):
                if cols:
                    inloop = False
                continue
            if not inloop or not cols:
                continue
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            p = line.split()
            if len(p) < len(cols):
                continue
            try:
                comp = p[cols["label_comp_id"]]
                el = p[cols["type_symbol"]].upper()
                xyz = (float(p[cols["Cartn_x"]]), float(p[cols["Cartn_y"]]),
                       float(p[cols["Cartn_z"]]))
            except (KeyError, ValueError, IndexError):
                continue
            if el == "H" or el == "D":
                continue
            if line.startswith("ATOM") or comp in AA:
                prot.append(xyz)
            else:
                key = (comp, p[cols.get("auth_asym_id", cols["label_asym_id"])],
                       p[cols.get("auth_seq_id", cols["label_seq_id"])])
                het[key].append(xyz)
    return np.array(prot, float) if prot else np.zeros((0, 3)), het


def main():
    os.makedirs(OUT, exist_ok=True)
    files = sorted(glob.glob(os.path.join(CACHE, "*.cif")))
    if not files:
        raise SystemExit(f"no CIFs in {CACHE}")
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    af = [f for f in files if os.path.basename(f).startswith("af_")]
    xr = [f for f in files if not os.path.basename(f).startswith("af_")]
    say(f"START/SRPBCC foldseek hit set: {len(files)} structures "
        f"({len(xr)} experimental, {len(af)} AlphaFold models)")
    say("AlphaFold models carry no ligands by construction, so only the")
    say("experimental subset can answer the question.")
    files = xr
    say(f"criteria: >= {MIN_HEAVY} heavy atoms, >= {BURIED_FRAC:.0%} of atoms "
        f"contacting protein within 4.5 A, >= {MIN_CONTACTS} contacts")
    say("")
    withlig, comps, per = [], Counter(), {}
    rejected = Counter()
    for f in files:
        pdb = os.path.basename(f)[:-4]
        try:
            prot, het = parse_cif(f)
        except Exception as e:                    # noqa: BLE001
            rejected["parse error"] += 1
            continue
        if len(prot) < 200:
            rejected["too few protein atoms"] += 1
            continue
        keep = []
        for (comp, ch, seq), xyz in het.items():
            if comp in JUNK or comp in GLYCAN or comp in MODRES or comp in AA:
                continue
            L = np.array(xyz, float)
            if len(L) < MIN_HEAVY:
                continue
            D = np.linalg.norm(L[:, None, :] - prot[None, :, :], axis=2)
            near = (D < 4.5)
            frac = float((near.any(1)).mean())
            ncon = int(near.sum())
            if frac >= BURIED_FRAC and ncon >= MIN_CONTACTS:
                keep.append((comp, len(L), frac, ncon))
        if keep:
            best = max(keep, key=lambda x: x[1])
            withlig.append((pdb, best))
            per[pdb] = keep
            # count STRUCTURES per component, not copies -- a homodimer with the
            # same ligand in both chains is one piece of evidence, not two
            for c in {k[0] for k in keep}:
                comps[c] += 1
    say(f"EXPERIMENTAL structures with >= 1 buried pocket ligand: {len(withlig)} of {len(files)}"
        f"  ({100*len(withlig)/len(files):.0f}%)")
    say(f"distinct ligand chemical components: {len(comps)}")
    say(f"skipped: {dict(rejected) if rejected else 'none'}")
    sizes = [b[1] for _p, b in withlig]
    say("")
    say(f"largest ligand per structure, heavy-atom count: median "
        f"{int(np.median(sizes))}, range {min(sizes)}-{max(sizes)}")
    say(f"  PYR1's ABA is 19 heavy atoms; mandipropamid 29")
    for lo, hi, lbl in ((0, 20, "<= 20  (ABA-sized)"),
                        (21, 28, "21-28"),
                        (29, 50, "29-50  (the §58c expansion band)"),
                        (51, 10**6, "> 50   (lipids, cofactors)")):
        n = sum(1 for s in sizes if lo <= s <= hi)
        say(f"    {lbl:<32} {n:>4}  ({100*n/len(sizes):4.1f}%)")
    say("")
    say("every structure with a buried pocket ligand:")
    say(f"   {'target':<12}{'ligand':<8}{'heavy':>7}{'buried':>8}{'contacts':>10}")
    for pdb, (comp, n, fr, nc) in sorted(withlig):
        say(f"   {pdb:<12}{comp:<8}{n:>7}{fr:>8.2f}{nc:>10}")
    say("")
    say("components, counted once per structure:")
    for c, n in comps.most_common(20):
        say(f"   {c:<6} {n}")
    say("")
    say("⚠ Structure count is not chemistry count -- the same protein solved many")
    say("  times with the same ligand inflates it. Distinct components is the")
    say("  conservative number.")
    with open(os.path.join(OUT, "census.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    with open(os.path.join(OUT, "census.json"), "w") as fh:
        json.dump({"n_structures": len(files), "with_ligand": len(withlig),
                   "components": dict(comps),
                   "per_structure": {k: v for k, v in per.items()}}, fh, indent=1)
    say(f"\nwritten to {OUT}/census.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
