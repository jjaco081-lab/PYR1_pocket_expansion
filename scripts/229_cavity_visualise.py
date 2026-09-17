#!/usr/bin/env python
r"""
229_cavity_visualise.py -- write the MEASURED cavity out as atoms so it can be
checked by eye in PyMOL instead of taken on trust.

Jannis: "For any claims that the pocket calculation has been corrected or that
pockets don't exist within a cavity that are easy for me to check on pymol please
have me do so."

Every number in README §134 comes from a voxel grid that nobody has ever looked
at. This dumps that grid. For each structure it writes:

  <name>_source.pdb    THE ORIGINAL FILE, copied verbatim -- real residue and
                       atom names, so PyMOL can draw a cartoon and so the LIGAND
                       is visible. ⚠ The first version rewrote the measured atoms
                       as UNK with element-only atom names; PyMOL could not build
                       a cartoon from it and the render showed cavity spheres
                       floating in empty space, which is useless for an eye check.
                       The measured-atom SELECTION is applied in the .pml instead.
  <name>_cavity.pdb    one HETATM per cavity voxel, as element He so PyMOL draws
                       spheres. CHAIN = chamber rank (A = largest, B = 2nd ...),
                       B-factor = that chamber's volume in A^3, occupancy = the
                       fraction of its lining that belongs to the CATH domain
  <name>.pml           loads both, colours chamber A red and the rest grey,
                       shows the domain in cartoon and non-domain in thin lines

⚠ The point of the occupancy column is the claim that needs checking hardest:
"this protein has a big chamber, but the SRPBCC domain does not line it". Chain A
with occupancy near 0 means the code found a large void somewhere else in the
protein. If it visibly sits in the helix-grip pocket, v4 is wrong.

WHAT TO LOOK AT, and what each case is supposed to show:
  3QN1_A      119.0 A^3 in 3 chambers -- the largest should sit on ABA.
  3QN1_AB     182.3 A^3 in ONE chamber -- HAB1 has no pocket residues but SEALS
              the openings, merging the satellites. This is the claim that PYR1's
              functional pocket is 182 and not 119.
  2pcsA00     570.5 A^3, the CoxG donor. Is that genuinely ONE chamber?
  af_A0A0R0IEX8_94_244   v2 said 318.6, v4 says 0.0. There IS a 277.6 A^3
              chamber; the claim is that it is NOT in the SRPBCC domain.
  af_Q6Z9J1_229_401      v2 216.6 -> v4 9.3 at domain fraction 0.76.
  2lf2A00     an NMR ensemble: ~20 models stacked into one chain. Before the fix
              this read 0.0 because the overlaid copies fill the pocket.
"""
import gzip, json, os, sys
import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m151 = import_module("151_cavity_bottleneck")
m221 = import_module("221_cavity_domain_trim")
m224 = import_module("224_cavity_domain_lined")
SPACING = m151.SPACING
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
OUT = os.path.join(ROOT, "results", "cavity_check")
CH = "ABCDEFGH"

TARGETS = [
    ("3QN1_A",  None, {"A"}, "PYR1 apo monomer: 119.0 A^3, largest chamber should sit on ABA"),
    ("3QN1_AB", None, None,  "PYR1+HAB1: 182.3 A^3 as ONE chamber -- HAB1 seals it"),
    ("2pcsA00", None, None,  "CoxG donor 570.5 A^3 -- is it really one chamber?"),
    ("af_A0A0R0IEX8_94_244_3.30.530.20", None, None,
     "v2 318.6 -> v4 0.0; a 277.6 A^3 chamber exists but should NOT be in the domain"),
    ("af_Q6Z9J1_229_401_3.30.530.20", None, None,
     "v2 216.6 -> v4 9.3, domain fraction 0.76"),
    ("2lf2A00", None, None, "NMR ensemble: read 0.0 before the first-model fix"),
]


def find(name):
    if name.startswith("3QN1"):
        return os.path.join(ROOT, "data", "3QN1_complex_auth_aba.pdb")
    for ext in (".cif", ".pdb"):
        p = os.path.join(RAW, name + ext)
        if os.path.exists(p):
            return p
    return None


def chambers_with_voxels(xyz, elem, dom):
    clear, _, lo, shape = m151.clear_map(xyz, elem)
    mask = m151.enclosed(clear, lo, shape)
    if mask.sum() == 0:
        return []
    core = mask & (clear > 2.4)
    lab, n = ndimage.label(core)
    if n == 0:
        return []
    sizes = ndimage.sum(core, lab, range(1, n + 1))
    order = np.argsort(sizes)[::-1][:8]
    seeds = []
    for k in order:
        idx = np.argwhere(lab == k + 1)
        seeds.append(idx[np.argmax(clear[tuple(idx.T)])])
    cvox = np.argwhere(mask)
    own = np.linalg.norm(cvox[:, None] - np.array(seeds)[None], axis=-1).argmin(1)
    X = np.asarray(xyz, float); D = np.asarray(dom, bool)
    out = []
    for j in range(len(seeds)):
        sel = cvox[own == j]
        if len(sel) == 0:
            continue
        pts = lo + sel * SPACING
        step = max(1, len(pts) // 4000)
        d = np.linalg.norm(X[None, :, :] - pts[::step][:, None, :], axis=2)
        near = (d < 5.0).any(0)
        nl = int(near.sum())
        frac = float(D[near].sum() / nl) if nl else 0.0
        out.append(dict(vol=len(sel) * SPACING ** 3, frac=frac, pts=pts))
    out.sort(key=lambda r: -r["vol"])
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    manifest = []
    for name, _u, keep, why in TARGETS:
        path = find(name)
        if not path:
            print(f"!! {name}: not found"); continue
        rng = m221.domain_range(name)
        ch = m221.read_trimmed(path, None)
        at = [a for c, v in ch.items() if (keep is None or c in keep) for a in v]
        if len(at) < 300:
            print(f"!! {name}: only {len(at)} atoms"); continue
        xyz = [x for _, x, _ in at]; el = [e for e, _, _ in at]
        dom = ([True] * len(at) if rng is None
               else [rng[0] <= r <= rng[1] for _, _, r in at])
        cs = chambers_with_voxels(xyz, el, dom)
        base = os.path.join(OUT, name.replace(".", "_"))

        # the ORIGINAL structure, verbatim, so cartoon and ligand render
        import shutil
        src = base + "_source" + os.path.splitext(path)[1]
        shutil.copyfile(path, src)
        # cavity voxels
        nkept = 0
        with open(base + "_cavity.pdb", "w") as fh:
            k = 0
            for ci, c in enumerate(cs[:len(CH)]):
                for pt in c["pts"]:
                    k += 1; nkept += 1
                    fh.write(f"HETATM{k%99999:5d} HE   CAV {CH[ci]}{(ci+1)%9999:4d}    "
                             f"{pt[0]:8.3f}{pt[1]:8.3f}{pt[2]:8.3f}"
                             f"{c['frac']:6.2f}{min(c['vol'],999.99):6.2f}"
                             f"          HE\n")
            fh.write("END\n")
        dsel = (f"resi {rng[0]}-{rng[1]}" if rng else "all")
        with open(base + ".pml", "w") as fh:
            fh.write(f"# {name}: {why}\n")
            fh.write(f"load {os.path.basename(src)}, prot\n")
            fh.write(f"load {os.path.basename(base)}_cavity.pdb, cav\n")
            fh.write("hide everything\n")
            fh.write(f"select domain, prot and polymer and ({dsel})\n")
            fh.write("select other, prot and polymer and not domain\n")
            fh.write("select ligand, prot and not polymer and not resn HOH\n")
            fh.write("show cartoon, domain\ncolor grey80, domain\n")
            fh.write("show lines, other\ncolor palecyan, other\n")
            fh.write("show sticks, ligand\ncolor yellow, ligand\n")
            fh.write("set stick_radius, 0.28, ligand\n")
            fh.write("util.cnc ligand\n")
            fh.write("show spheres, cav\nset sphere_scale, 0.20, cav\n")
            # translucent, or the blob hides the very ligand it should overlay
            fh.write("set sphere_transparency, 0.55, cav\n")
            fh.write("color firebrick, cav and chain A\n")
            for c_ in CH[1:]:
                fh.write(f"color grey50, cav and chain {c_}\n")
            fh.write("set cartoon_transparency, 0.55\nbg_color white\n")
            fh.write("deselect\n")
            fh.write("orient cav and chain A\nzoom cav and chain A, 7\n")
            fh.write(f'print "CHAMBERS (chain = rank, red = largest):"\n')
            for ci, c in enumerate(cs[:len(CH)]):
                fh.write(f'print "  {CH[ci]}  {c["vol"]:7.1f} A^3   '
                         f'domain-lined fraction {c["frac"]:.2f}"\n')
            fh.write('print "grey cartoon = CATH domain; pale lines = rest of '
                     'the file; YELLOW STICKS = ligand; red spheres = reported '
                     'chamber"\n')
            fh.write('print "the check: does the RED blob sit on the yellow '
                     'ligand, inside the grey cartoon?"\n')
        print(f"{name}")
        print(f"   {why}")
        for ci, c in enumerate(cs[:6]):
            tag = "  <-- REPORTED as the pocket" if ci == 0 and c["frac"] >= 2/3 else \
                  ("  <-- largest, but NOT domain-lined" if ci == 0 else "")
            print(f"   chamber {CH[ci]}  {c['vol']:7.1f} A^3   "
                  f"domain fraction {c['frac']:.2f}{tag}")
        print(f"   -> {base}.pml   ({nkept} cavity points)\n")
        manifest.append(dict(target=name, why=why, pml=base + ".pml",
                             chambers=[dict(vol=round(c["vol"], 1),
                                            frac=round(c["frac"], 3)) for c in cs[:8]]))
    json.dump(manifest, open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
    print(f"wrote {OUT}/  -- open any .pml with:  pymol <file>.pml")
    return 0


if __name__ == "__main__":
    sys.exit(main())
