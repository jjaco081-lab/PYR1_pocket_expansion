#!/usr/bin/env python
r"""96_mol2_setcoords.py -- write Amber restart coordinates into a mol2 template.

    96_mol2_setcoords.py <template.mol2> <coords.rst7> <out.mol2>

tleap's loadmol2 preserves atom order, so the restart written after minimising
that unit is in the same order as the mol2 it came from. Asserted on count
anyway -- a silent off-by-one here would shift every charge one atom along.

Read through ParmEd rather than by hand: Amber writes NetCDF restarts by default
(ntxo=2), so a text parser hits 'utf-8 codec can't decode byte 0xb8' on a file
that is perfectly valid. ParmEd handles both encodings.
"""
import sys

import parmed as pmd


def main():
    tpl_f, rst_f, out_f = sys.argv[1], sys.argv[2], sys.argv[3]
    xyz = pmd.load_file(rst_f).coordinates
    if xyz is None:
        raise SystemExit(f"{rst_f}: no coordinates")
    xyz = xyz.reshape(-1, 3)
    n = len(xyz)

    lines, fl, i = open(tpl_f).read().split("\n"), False, 0
    for k, l in enumerate(lines):
        if l.startswith("@<TRIPOS>ATOM"):
            fl = True
            continue
        if l.startswith("@<TRIPOS>BOND"):
            fl = False
        if fl and l.strip():
            p = l.split()
            if i >= n:
                raise SystemExit(f"mol2 has more atoms than the restart ({n})")
            x, y, z = xyz[i]
            lines[k] = (f"{int(p[0]):>7} {p[1]:<8} {x:9.4f} {y:9.4f} {z:9.4f} "
                        f"{p[5]:<8} {p[6]:>4} {p[7]:<8} {float(p[8]):9.6f}")
            i += 1
    if i != n:
        raise SystemExit(f"mol2 has {i} atoms, restart has {n} -- refusing to write")
    open(out_f, "w").write("\n".join(lines))
    print(f"    hydrogens relaxed; wrote {out_f.split('/')[-1]} ({i} atoms)")


if __name__ == "__main__":
    main()
