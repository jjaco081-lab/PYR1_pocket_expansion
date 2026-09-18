#!/usr/bin/env python
r"""
250_triple_chamber.py -- does the F108/R79/E94 triple open a USABLE chamber, or
just enclosed space?

⚠ WHY THIS HAS TO BE RUN BEFORE ANY MD. §3's numbers come from the RIGID-BACKBONE
truncation scan and are TOTAL ENCLOSED VOLUME:

    R79A/E94A            251.1   (+84.8)
    F108A/R79A/E94A      308.5  (+142.1)
    K59A/F108A/R79A/E94A 375.5  (+209.1)

The same method gave F108G +84 A^3 of enclosed volume and only **+1 A^3 of
usable chamber** once a ligand-sized probe was required, and 0 of 50 borrowed
positions opened any chamber volume at all. So +142 may be satellite space a
ligand cannot occupy. Re-measured here with the CALIBRATED chamber probe (1.4 A,
the radius at which 18 of 19 ABA atoms sit inside PYR1's own chamber).

⚠ SIDE CHAINS ARE TRUNCATED IN PLACE, NOT REPACKED -- deliberately, to reproduce
§3's construction so the two numbers are comparable. That makes this an UPPER
BOUND: it is the volume before any relaxation closes it, which is exactly why the
collapse question needs MD and cannot be answered here.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
imp = __import__("importlib").import_module
m221, m229 = imp("221_cavity_domain_trim"), imp("229_cavity_visualise")
PDB = os.path.join(ROOT, "data", "3QN1_complex_auth_aba.pdb")
#: atoms kept when a residue is truncated to alanine
ALA = {"N", "CA", "C", "O", "CB", "OXT"}
THREE = {"ARG": "R", "GLU": "E", "PHE": "F", "LYS": "K", "VAL": "V", "ALA": "A"}
WANT = {59: "LYS", 79: "ARG", 94: "GLU", 108: "PHE"}


def load(keep_chains, trunc):
    """heavy atoms of 3QN1; residues in `trunc` cut back to alanine."""
    out = []
    for l in open(PDB):
        if not l.startswith("ATOM"):
            continue
        if l[21] not in keep_chains:
            continue
        e = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if e in ("H", "D"):
            continue
        r = int(l[22:26])
        nm = l[12:16].strip()
        if l[21] == "A" and r in trunc and nm not in ALA:
            continue
        out.append((e, np.array([float(l[30:38]), float(l[38:46]),
                                 float(l[46:54])]), r))
    return out


def measure(at):
    X = np.array([x for _, x, _ in at])
    el = [e for e, _, _ in at]
    cs = m229.chambers_with_voxels(X, el, [True] * len(at))
    if not cs:
        return 0.0, 0.0, 0
    tot = sum(c["vol"] for c in cs)
    return cs[0]["vol"], tot, len(cs)


def main():
    # identity assertion FIRST -- never trust the numbers
    seq = {}
    for l in open(PDB):
        if l.startswith("ATOM") and l[21] == "A" and l[12:16].strip() == "CA":
            seq[int(l[22:26])] = l[17:20].strip()
    for r, want in WANT.items():
        assert seq.get(r) == want, f"residue {r} is {seq.get(r)}, expected {want}"
    print("identities asserted: K59, R79, E94, F108 all confirmed\n")

    #: ⚠ R79A ALONE WAS NEVER MEASURED. §3's scan covered first-shell residues
    #: only and R79 is second shell, so the "R79/E94 are synergistic not
    #: additive" claim rests on E94A alone (+18.4) and the pair (+84.8) with the
    #: other single missing. If R79A alone gives ~+66 the pair is additive.
    #: Both singles are included here so the interaction term is computable.
    VARIANTS = [("WT", set()),
                ("R79A", {79}),
                ("E94A", {94}),
                ("R79A/E94A", {79, 94}),
                ("F108A", {108}),
                ("F108A/R79A", {79, 108}),
                ("F108A/R79A/E94A", {79, 94, 108}),
                ("K59A/F108A/R79A/E94A", {59, 79, 94, 108})]
    #: §3's numbers, for the side-by-side
    OLD = {"WT": 166.4, "R79A": None, "E94A": 184.8, "F108A": 239.9,
           "R79A/E94A": 251.1, "F108A/R79A": 265.1,
           "F108A/R79A/E94A": 308.5, "K59A/F108A/R79A/E94A": 375.5}
    for lbl, keep in (("PYR1 chain A alone", {"A"}),
                      ("PYR1 + HAB1 (closed complex)", {"A", "B"})):
        print(f"=== {lbl} ===")
        print(f"{'variant':<24}{'MAIN chamber':>14}{'vs WT':>8}"
              f"{'all enclosed':>14}{'n':>4}{'§3 total':>10}")
        base = None
        for name, trunc in VARIANTS:
            at = load(keep, trunc)
            main, tot, n = measure(at)
            if base is None:
                base = main
            o = OLD.get(name)
            print(f"{name:<24}{main:>14.1f}{main-base:>+8.1f}{tot:>14.1f}{n:>4}"
                  f"{(f'{o:.1f}' if o else '-'):>10}")
        print()
    # ---- the interaction term the original scan could not compute ----
    print("=" * 66)
    print("IS R79/E94 ACTUALLY SYNERGISTIC? (chain A alone, calibrated probe)")
    print("=" * 66)
    at = load({"A"}, set()); wt, _, _ = measure(at)
    d = {}
    for nm, tr in (("R79A", {79}), ("E94A", {94}), ("R79A/E94A", {79, 94})):
        m_, _, _ = measure(load({"A"}, tr)); d[nm] = m_ - wt
    add = d["R79A"] + d["E94A"]
    obs = d["R79A/E94A"]
    print(f"   R79A alone      {d['R79A']:+7.1f} A^3   <-- never measured before")
    print(f"   E94A alone      {d['E94A']:+7.1f} A^3")
    print(f"   sum (additive)  {add:+7.1f} A^3")
    print(f"   pair observed   {obs:+7.1f} A^3")
    if abs(add) > 1e-6:
        print(f"   observed / additive = {obs/add:.2f}"
              + ("   SUPER-ADDITIVE" if obs > add * 1.3 else
                 ("   SUB-ADDITIVE" if obs < add * 0.77 else "   ~ADDITIVE")))
    print()
    print("⚠ MAIN chamber is the usable, ligand-accessible volume (probe 1.4 A).")
    print("  '§3 total' is the old rigid-scan TOTAL enclosed volume. If the main")
    print("  chamber does not grow with the triple, the +142 was satellite space")
    print("  and there is no pocket to test for collapse.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
