#!/usr/bin/env python
r"""
149_pocket_selections.py -- audit the lining calls against each structure's OWN
crystallographic ligand, and emit PyMOL selections for inspection.

WHY. Script 143 reported 32 lining residues for 3OQU against 24 for PYR1, from a
cavity detected geometrically with all ligands stripped. Jannis, looking at the
two pockets in PyMOL, is not convinced 3OQU really has that many more. That is a
falsifiable objection and this script tests it the only honest way: both
structures carry ABA (A8S) in the crystal, so every residue called "lining" can
be scored by its distance to the ligand that is actually there.

If the extra 3OQU residues sit far from A8S, the detected component has run into
a neighbouring void and the count is inflated. If they cluster at 4-6 A like the
rest, the count stands.

The two numbering frames differ and both are stated explicitly:
  PYR1   data/stage1/wt_aba.pdb, chain A of 3QN1 -- so the residue numbers below
         are 3QN1 auth numbering and can be pasted against a fresh 3QN1
  3OQU   chain B, auth numbering straight from the deposited mmCIF

Emits copy-pasteable PyMOL for both.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m143 = import_module("143_pocket_backbone_shape")                   # noqa: E402
read_pdb, read_cif_resnames, THREE = m143.read_pdb, m143.read_cif_resnames, m143.THREE

HC = os.path.join(ROOT, "results", "homolog_cavities")
OUT = os.path.join(ROOT, "results", "pocket_shape")
NEAR = 4.5      # A -- the conventional contact shell, fixed before looking


def cif_atoms(path, chain, group="HETATM", comp=None):
    cols, rows, in_loop = {}, [], False
    for line in open(path):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".")[1]] = len(cols)
            in_loop = True
            continue
        if in_loop:
            if line.startswith(("#", "loop_", "_")):
                if rows:
                    break
                continue
            f = line.split()
            if len(f) >= len(cols):
                rows.append(f)
    out = []
    for f in rows:
        if f[cols["group_PDB"]] != group or f[cols["auth_asym_id"]] != chain:
            continue
        if comp and f[cols["label_comp_id"]] != comp:
            continue
        if f[cols["type_symbol"]] == "H":
            continue
        out.append((f[cols["label_comp_id"]], f[cols["auth_seq_id"]],
                    np.array([float(f[cols["Cartn_x"]]), float(f[cols["Cartn_y"]]),
                              float(f[cols["Cartn_z"]])])))
    return out


def report(label, resnums, coords_by_res, lig, ss):
    rows = []
    for r in resnums:
        d = coords_by_res.get(r)
        if d is None:
            continue
        allat = np.array(list(d["atoms"].values()))
        dm = float(np.linalg.norm(lig[:, None] - allat[None], axis=-1).min())
        rows.append((r, d["aa"], ss.get(r, "?"), dm))
    rows.sort(key=lambda x: x[3])
    near = [x for x in rows if x[3] <= NEAR]
    print(f"\n{label}: {len(rows)} residues called lining, "
          f"{len(near)} within {NEAR} A of the crystal ligand "
          f"({len(near)/max(len(rows),1):.0%})")
    print(f"   {'res':>5}{'aa':>4}{'ss':>4}{'d(lig)':>9}")
    for r, aa, s, dm in rows:
        flag = "" if dm <= NEAR else ("  <-- outside the contact shell"
                                      if dm > 6.0 else "  (second shell)")
        print(f"   {r:>5}{aa:>4}{s:>4}{dm:>9.2f}{flag}")
    return rows


def main():
    shp = {r["name"]: r for r in json.load(open(os.path.join(OUT, "pocket_shape.json")))}

    # ---------------- PYR1, 3QN1 chain A numbering ----------------------
    at = read_pdb(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb"), want_atom=False)
    lig1 = np.array([a[4] for a in at if a[1] == "A8S"])
    assert len(lig1), "no A8S in the PYR1 frame"
    p_by = {}
    for a in at:
        if a[1] in THREE:
            p_by.setdefault(a[0], {"aa": THREE[a[1]], "atoms": {}})["atoms"][a[2]] = a[4]
    p = shp["PYR1"]
    p_ss = {n: p["ss_string"][i] for i, n in enumerate(p["ca_resnums"])}
    print("=" * 70)
    print("PYR1  --  data/stage1/wt_aba.pdb = 3QN1 chain A, auth numbering")
    print(f"          ligand A8S, {len(lig1)} heavy atoms")
    print("=" * 70)
    r1 = report("PYR1 pocket", p["lining"], p_by, lig1, p_ss)

    # ---------------- 3OQU chain B -------------------------------------
    tgt = "3oquB00"
    cif = os.path.join(HC, "raw", tgt + ".cif")
    names = read_cif_resnames(cif, "B")
    dom = read_pdb(os.path.join(HC, "domains", tgt + ".pdb"))
    d_by = {}
    for a in dom:
        d_by.setdefault(a[0], {"aa": names.get(a[0], "X"), "atoms": {}})["atoms"][a[2]] = a[4]
    het = cif_atoms(cif, "B", comp="A8S")
    if not het:
        het = [h for h in cif_atoms(cif, "B") if h[0] not in ("HOH", "MN", "SO4", "GOL")]
    lig2 = np.array([h[2] for h in het])
    print("\n" + "=" * 70)
    print(f"3OQU  --  chain B, auth numbering from the deposited mmCIF")
    print(f"          ligand {het[0][0] if het else '?'} {len(lig2)} heavy atoms")
    print("=" * 70)
    d = shp["3OQU_B"]
    d_ss = {n: d["ss_string"][i] for i, n in enumerate(d["ca_resnums"])}
    r2 = report("3OQU pocket", d["lining"], d_by, lig2, d_ss)

    # ---------------- the verdict on the extra residues -----------------
    n1 = sum(1 for x in r1 if x[3] <= NEAR)
    n2 = sum(1 for x in r2 if x[3] <= NEAR)
    print("\n" + "=" * 70)
    print("VERDICT on the extra lining count")
    print("=" * 70)
    print(f"   geometric cavity call:  PYR1 {len(r1)}  vs  3OQU {len(r2)}   "
          f"(+{len(r2)-len(r1)})")
    print(f"   within {NEAR} A of ABA:  PYR1 {n1}  vs  3OQU {n2}   (+{n2-n1})")
    print(f"   so of 3OQU's {len(r2)-len(r1)} extra residues, "
          f"{n2-n1} are genuine ligand-contact positions and "
          f"{(len(r2)-len(r1))-(n2-n1)} sit further out.")

    # ---------------- PyMOL -------------------------------------------
    def sel(rs):
        return "+".join(str(r) for r in sorted(rs))
    lines = [
        "# ---- PYR1 / 3QN1 (chain A) ----",
        "fetch 3qn1, async=0",
        "hide everything, 3qn1",
        "show cartoon, 3qn1 and chain A",
        "select pyr1_lig, 3qn1 and chain A and resn A8S",
        f"select pyr1_wall, 3qn1 and chain A and resi {sel(x[0] for x in r1)}",
        f"select pyr1_wall_contact, 3qn1 and chain A and resi "
        f"{sel(x[0] for x in r1 if x[3] <= NEAR)}",
        "show sticks, pyr1_wall or pyr1_lig",
        "color grey80, pyr1_wall",
        "color salmon, pyr1_wall_contact",
        "color yellow, pyr1_lig",
        "",
        "# ---- 3OQU (chain B) ----",
        "fetch 3oqu, async=0",
        "hide everything, 3oqu",
        "show cartoon, 3oqu and chain B",
        "select oqu_lig, 3oqu and chain B and resn A8S",
        f"select oqu_wall, 3oqu and chain B and resi {sel(x[0] for x in r2)}",
        f"select oqu_wall_contact, 3oqu and chain B and resi "
        f"{sel(x[0] for x in r2 if x[3] <= NEAR)}",
        "show sticks, oqu_wall or oqu_lig",
        "color grey80, oqu_wall",
        "color skyblue, oqu_wall_contact",
        "color yellow, oqu_lig",
        "",
        "# ---- superpose and compare ----",
        "align 3oqu and chain B, 3qn1 and chain A",
        "# cavity surfaces:",
        "set surface_cavity_mode, 1",
        "set surface_cavity_radius, 5",
        "show surface, 3qn1 and chain A",
        "show surface, 3oqu and chain B",
    ]
    txt = "\n".join(lines)
    with open(os.path.join(OUT, "pocket_selections.pml"), "w") as fh:
        fh.write(txt + "\n")
    print("\n" + "=" * 70)
    print("PYMOL  (also written to results/pocket_shape/pocket_selections.pml)")
    print("=" * 70)
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
