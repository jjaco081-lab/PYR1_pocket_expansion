#!/usr/bin/env python
r"""
97_pose_map.py -- map a docked pose back onto the parameterised molecule by
GRAPH MATCHING rather than by atom order.

    97_pose_map.py <sybyl_template.mol2> <docked.mol2> <gaff_template.mol2> <out.mol2>

WHY NOT JUST COPY COORDINATES IN ORDER
Measured on this pipeline, the atom order survives openbabel's --gen3d intact but
does NOT survive smina:

    template : C  C  C  Cl C  C  C  C  O  ...
    docked   : C  O  N  H  C  C  C  C  Cl ...

smina also discards every nonpolar hydrogen (51 atoms in, 30 out), so neither the
count nor the order can be relied on. An order-based transfer would put a
chlorine's coordinates on a carbon and nothing downstream would complain -- the
system would build, run, and return a confident number for the wrong molecule.

So the two molecules are matched as GRAPHS. Chirality cannot be part of the match
(smina removed the stereocentre's hydrogen, so RDKit sees no centre in the docked
file at all); topology alone is unambiguous here, and the stereocentre is verified
separately by 98_check_stereo.py once the hydrogens are back.

Only the 29 heavy atoms are transferred. Every hydrogen -- polar included -- is
left at its template position and relaxed afterwards by 91_quad_build.sh with the
heavy atoms restrained.

Needs RDKit: run under conda_envs/esmfold2/bin/python.
"""
import math
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def heavy_view(mol):
    """RemoveHs keeps the surviving atoms in their original relative order, so
    position j of the result corresponds to the j-th heavy atom of the input."""
    idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    return Chem.RemoveHs(mol), idx


def main():
    sybyl_f, docked_f, gaff_f, out_f = sys.argv[1:5]

    tmpl = Chem.MolFromMol2File(sybyl_f, removeHs=False, sanitize=True)
    dock = Chem.MolFromMol2File(docked_f, removeHs=False, sanitize=True)
    if tmpl is None or dock is None:
        raise SystemExit("RDKit could not parse one of the mol2 files")
    Chem.AssignStereochemistryFrom3D(tmpl)
    Chem.AssignStereochemistryFrom3D(dock)

    t_h, t_map = heavy_view(tmpl)
    d_h, _ = heavy_view(dock)
    if t_h.GetNumAtoms() != d_h.GetNumAtoms():
        raise SystemExit(f"heavy-atom count differs: template {t_h.GetNumAtoms()} "
                         f"vs docked {d_h.GetNumAtoms()}")

    # Chirality CANNOT be enforced in the match: smina stripped the hydrogen on
    # the stereocentre, so RDKit perceives no chiral centre in the docked file at
    # all (template [(7,'S')], docked []). That is a loss of annotation, not an
    # inversion -- Vina only rotates rotatable bonds and never inverts a centre.
    # Topology alone is still an unambiguous map here: the four substituents at
    # C7 (amide, propargyl ether, chlorophenyl, H) are all topologically
    # distinct, so no graph automorphism can swap them. The only symmetry in the
    # molecule is the chlorophenyl ring flip, which is geometrically harmless.
    # The stereocentre is then VERIFIED on the rebuilt molecule after the
    # hydrogens are relaxed back on (91_quad_build.sh, step 3c).
    match = d_h.GetSubstructMatch(t_h, useChirality=False)
    if not match:
        raise SystemExit("no graph match -- the docked molecule is not the template")
    if len(match) != t_h.GetNumAtoms():
        raise SystemExit(f"partial match ({len(match)}/{t_h.GetNumAtoms()})")

    dconf = d_h.GetConformer()
    tconf = tmpl.GetConformer()
    new = {}
    for j, orig in enumerate(t_map):
        p = dconf.GetAtomPosition(match[j])
        new[orig] = (p.x, p.y, p.z)

    # rewrite the GAFF2 template, which carries the types and AM1-BCC charges
    lines = open(gaff_f).read().split("\n")
    fl, i, n_set = False, 0, 0
    for k, l in enumerate(lines):
        if l.startswith("@<TRIPOS>ATOM"):
            fl = True
            continue
        if l.startswith("@<TRIPOS>BOND"):
            fl = False
        if fl and l.strip():
            p = l.split()
            if i in new:
                x, y, z = new[i]
                n_set += 1
            else:
                x, y, z = float(p[2]), float(p[3]), float(p[4])
            lines[k] = (f"{int(p[0]):>7} {p[1]:<8} {x:9.4f} {y:9.4f} {z:9.4f} "
                        f"{p[5]:<8} {p[6]:>4} {p[7]:<8} {float(p[8]):9.6f}")
            i += 1
    if i != tmpl.GetNumAtoms():
        raise SystemExit(f"GAFF template has {i} atoms, sybyl template has "
                         f"{tmpl.GetNumAtoms()} -- these are not the same molecule")
    if n_set != len(t_map):
        raise SystemExit(f"set {n_set} of {len(t_map)} heavy atoms")
    open(out_f, "w").write("\n".join(lines))

    rms = math.sqrt(sum(
        (new[o][0] - tconf.GetAtomPosition(o).x) ** 2 +
        (new[o][1] - tconf.GetAtomPosition(o).y) ** 2 +
        (new[o][2] - tconf.GetAtomPosition(o).z) ** 2
        for o in t_map) / len(t_map))
    ca = [sum(new[o][j] for o in t_map) / len(t_map) for j in range(3)]
    cb = [sum(getattr(tconf.GetAtomPosition(o), ax) for o in t_map) / len(t_map)
          for ax in "xyz"]
    com = math.sqrt(sum((u - v) ** 2 for u, v in zip(ca, cb)))
    print(f"    graph match OK: {len(match)} heavy atoms (topology; stereo checked later)")
    print(f"    {tmpl.GetNumAtoms()-len(t_map)} hydrogens left at template "
          f"positions -- relaxed next")
    print(f"    DOCKED vs CRYSTAL: heavy-atom RMSD {rms:.2f} A, "
          f"centroid offset {com:.2f} A")


def all_pose_rmsd(sybyl_f, docked_f):
    """RMSD to the crystal pose for every returned pose.

    The top-scored pose is what a production run would take, but knowing whether
    a BETTER pose was sampled and merely mis-ranked is a different diagnosis from
    never sampling it at all -- the first is a scoring-function problem, the
    second a search problem.
    """
    tmpl = Chem.MolFromMol2File(sybyl_f, removeHs=False, sanitize=True)
    t_h, t_map = heavy_view(tmpl)
    tconf = tmpl.GetConformer()
    ref = [(tconf.GetAtomPosition(o).x, tconf.GetAtomPosition(o).y,
            tconf.GetAtomPosition(o).z) for o in t_map]
    out = []
    blocks = open(docked_f).read().split("@<TRIPOS>MOLECULE")[1:]
    import tempfile, os
    for n, b in enumerate(blocks, 1):
        fh = tempfile.NamedTemporaryFile("w", suffix=".mol2", delete=False)
        fh.write("@<TRIPOS>MOLECULE" + b)
        fh.close()
        d = Chem.MolFromMol2File(fh.name, removeHs=False, sanitize=True)
        os.unlink(fh.name)
        if d is None:
            out.append((n, None))
            continue
        d_h, _ = heavy_view(d)
        m = d_h.GetSubstructMatch(t_h, useChirality=False)
        if not m:
            out.append((n, None))
            continue
        c = d_h.GetConformer()
        r = math.sqrt(sum((c.GetAtomPosition(m[j]).x - ref[j][0]) ** 2 +
                          (c.GetAtomPosition(m[j]).y - ref[j][1]) ** 2 +
                          (c.GetAtomPosition(m[j]).z - ref[j][2]) ** 2
                          for j in range(len(ref))) / len(ref))
        out.append((n, r))
    return out


if __name__ == "__main__":
    main()
    rs = all_pose_rmsd(sys.argv[1], sys.argv[2])
    good = [(n, r) for n, r in rs if r is not None]
    print("    RMSD to crystal pose, all returned poses: " +
          ", ".join(f"#{n} {r:.2f}" for n, r in good))
    if good:
        bn, br = min(good, key=lambda x: x[1])
        print(f"    best sampled pose is #{bn} at {br:.2f} A "
              f"({'search' if br > 2.5 else 'scoring'} is the limiting step)")
