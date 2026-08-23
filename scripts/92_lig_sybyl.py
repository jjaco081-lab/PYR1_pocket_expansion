#!/usr/bin/env python
r"""
92_lig_sybyl.py -- rewrite a gaff2-typed mol2 with Sybyl atom types.

WHY THIS EXISTS
antechamber writes mol2 files carrying GAFF2 atom types (ca, c3, os, ha, ...).
Those are not Sybyl types, and neither openbabel nor smina can read them:
openbabel warns "Cannot perform atom type translation" and smina silently types
the whole molecule as Du (dummy) and Ca (calcium). Docking then runs on a
molecule that is not the ligand -- it returns a plausible-looking score for a
structure made of dummies. That is a silent failure of exactly the kind
README 33 already records three of.

The molecule itself is fine: the mol2 BOND block carries real bond orders
(1/2/3/ar), so only the type column needs translating. Atom order, names,
coordinates and AM1-BCC charges are preserved untouched, which is what lets the
docked coordinates be transferred back onto the parameterised molecule later.
"""
import sys

# complete for the types present in A8S and 3UZ; anything else fails loudly
GAFF2SYBYL = {
    "ca": "C.ar", "c": "C.2", "c1": "C.1", "c2": "C.2", "c3": "C.3",
    "cl": "Cl", "br": "Br", "f": "F", "i": "I",
    "n": "N.am", "ns": "N.am", "nh": "N.pl3", "n3": "N.3", "n4": "N.4",
    "na": "N.ar", "nb": "N.ar", "nc": "N.2", "nd": "N.2",
    "o": "O.2", "os": "O.3", "oh": "O.3", "o2": "O.co2",
    "s": "S.3", "ss": "S.3", "sy": "S.o2", "s6": "S.o2",
    "p5": "P.3",
    "h1": "H", "h2": "H", "h3": "H", "h4": "H", "h5": "H",
    "ha": "H", "hc": "H", "hn": "H", "ho": "H", "hs": "H", "hp": "H",
}


def convert(src, dst):
    lines = open(src).read().split("\n")
    out, fl, n, seen = [], False, 0, set()
    for l in lines:
        if l.startswith("@<TRIPOS>ATOM"):
            fl = True; out.append(l); continue
        if l.startswith("@<TRIPOS>BOND"):
            fl = False
        if fl and l.strip():
            p = l.split()
            g = p[5]
            seen.add(g)
            if g not in GAFF2SYBYL:
                raise SystemExit(f"{src}: gaff2 type '{g}' has no Sybyl mapping "
                                 "-- add it rather than letting it through")
            # fixed-width so downstream readers stay happy
            out.append(f"{int(p[0]):>7} {p[1]:<8} {float(p[2]):9.4f} {float(p[3]):9.4f} "
                       f"{float(p[4]):9.4f} {GAFF2SYBYL[g]:<8} {p[6]:>4} {p[7]:<8} "
                       f"{float(p[8]):9.6f}")
            n += 1
            continue
        out.append(l)
    open(dst, "w").write("\n".join(out))
    print(f"    {src.split('/')[-1]} -> {dst.split('/')[-1]}: {n} atoms retyped "
          f"({len(seen)} distinct gaff2 types)")
    return n


if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2])
