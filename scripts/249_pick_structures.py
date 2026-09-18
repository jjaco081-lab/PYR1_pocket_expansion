#!/usr/bin/env python
r"""
249_pick_structures.py -- write named RFd3 designs out for inspection, coloured
by core / orange / truncated, with the ligand carried through if present.

Jannis: "I want any with consecutive oranges above 5 to be removed, not total.
Either way those two sound interesting, please give me the structures."

⚠ THE FILTER CHANGES THE ANSWER. Filtering on TOTAL orange keeps 42 of 655
designs and 2 above PYR1; filtering on the LONGEST CONSECUTIVE RUN -- which is
what was asked for, and what "no long orange stretches" means -- keeps 444 and
38. `arm4rep_s3004/d169` at 294.8 A^3 (1.79x PYR1) has 11 orange residues in
total but its longest run is 5, so it survives the correct filter and failed the
wrong one.

COLOURING (from 220.write_pdb): B-factor 0 = core, 1 = kept but NOT core
(orange), 2 = truncated. A `.pml` is written per design.
⚠ HAB1 is written as its FOUR separate chains; flattening them made their residue
numbers collide and PyMOL drew disconnected lines.
"""
import argparse, glob, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
m220 = __import__("importlib").import_module("220_core_truncate")
OUT = os.path.join(ROOT, "results", "rfd3", "picked")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--picks", nargs="+", required=True,
                    help="arm/design, e.g. arm6helix_s6001/d309")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    sc = {}
    p = os.path.join(ROOT, "results", "rfd3", "pilot_backbones.json")
    if os.path.exists(p):
        for r in json.load(open(p))["all_scored"]:
            sc[(r["arm"], r["design"])] = r
    print(f"{'design':<30}{'cav':>8}{'xPYR1':>7}{'SS':>6}{'tot':>5}{'run':>5}"
          f"{'kept':>6}{'cutN':>6}{'cutC':>6}")
    rows = []
    for pick in a.picks:
        arm, d = pick.split("/")
        g = glob.glob(os.path.join(ROOT, "results", "rfd3", arm, d, "*.cif.gz"))
        if not g:
            print(f"!! {pick}: no structure"); continue
        at, dat, dc, by, main, motif = m220.load_design(g[0])
        core = m220.core_of(by, main, motif)
        kept, cN, cC, run = m220.truncate(by, main, motif)
        tag = f"{arm}_{d}"
        pdb = os.path.join(OUT, f"{tag}.pdb")
        m220.write_pdb(at, dc, kept, core, pdb, f"{tag}  core-growth truncation")
        s = sc.get((arm, d), {})
        print(f"{pick:<30}{s.get('cavity',0):>8.1f}"
              f"{s.get('cavity',0)/164.4:>7.2f}{s.get('frac_SS',0):>6.2f}"
              f"{s.get('n_orange','-'):>5}{s.get('longest_orange','-'):>5}"
              f"{len(kept):>6}{cN:>6}{cC:>6}")
        pml = os.path.join(OUT, f"{tag}.pml")
        with open(pml, "w") as fh:
            fh.write(f"# {pick}\n")
            fh.write(f"# cavity {s.get('cavity',0):.1f} A^3 "
                     f"({s.get('cavity',0)/164.4:.2f}x PYR1), "
                     f"SS {s.get('frac_SS',0):.2f}, "
                     f"orange {s.get('n_orange','?')} total / "
                     f"{s.get('longest_orange','?')} longest run\n")
            fh.write(f"load {os.path.basename(pdb)}, {tag}\n")
            fh.write("hide everything\n")
            fh.write(f"show cartoon, {tag} and chain A\n")
            fh.write(f"spectrum b, grey70 orange red, {tag} and chain A, 0, 2\n")
            fh.write(f"show cartoon, {tag} and not chain A\n")
            fh.write(f"color palecyan, {tag} and not chain A\n")
            fh.write(f"set cartoon_transparency, 0.55, {tag} and not chain A\n")
            fh.write(f"show sticks, {tag} and not polymer\n")
            fh.write(f"color yellow, {tag} and not polymer\n")
            fh.write("bg_color white\n")
            fh.write(f"orient {tag} and chain A\n")
            fh.write('print "GREY = core   ORANGE = kept but not core   '
                     'RED = truncated"\n')
            fh.write('print "pale cyan = HAB1 context (four chains)   '
                     'yellow = ligand if present"\n')
        rows.append(dict(pick=pick, pdb=pdb, pml=pml, kept=len(kept),
                         cut_N=cN, cut_C=cC, longest_run=run, **
                         {k: s.get(k) for k in
                          ("cavity", "frac_SS", "n_orange", "longest_orange")}))
    json.dump(rows, open(os.path.join(OUT, "picked.json"), "w"), indent=1)
    print(f"\nwrote {OUT}/  ({len(rows)} designs)")
    print("  PyMOL:  pymol <design>.pml")
    return 0


if __name__ == "__main__":
    sys.exit(main())
