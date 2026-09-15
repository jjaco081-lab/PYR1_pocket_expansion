#!/usr/bin/env python
r"""
196_break_contact_audit.py -- WHERE are the chain breaks and WHERE does the
design touch HAB1, per segment.

Jannis inspected arm0b d000/d002/d009 and reported the opposite ranking to my
filters: d000 looks clean to him, while d002 and d009 -- the two my cascade
preferred -- have "weird appendages still reaching up to HAB1".

That is a specific, checkable claim that my F4 filter CANNOT see. F4 counts
design residues within 4.5 A of HAB1 *in the leading 25 residues only*, because
the defect it was written for was the unanchored N-terminal range in batch 1-3.
A tail off the C-terminus, or a long inter-motif linker, contacts HAB1 without
ever entering F4's window. So "leading HAB1 = 0" is not "does not touch HAB1".

This reports, per design:
  * every chain break (consecutive residue numbers, CA-CA > 4.5 A) with the
    segment it falls in and the gap distance
  * HAB1 contacts broken out by segment: leading, each inter-motif linker,
    trailing -- so an appendage is attributed to the segment that carries it
  * the longest run of consecutive non-motif residues contacting HAB1, which is
    the "appendage" in the sense Jannis means it

Usage:
  python 196_break_contact_audit.py results/rfd3/arm0b_s2027/d000 [more dirs]
"""
import glob, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m195 = import_module("195_arm_compare")

CONTACT = 4.5


def segments(motif_out):
    """[(label, lo, hi)] over the design chain: lead, linkers, trail."""
    runs = []
    for r in sorted(motif_out):
        if runs and r == runs[-1][1] + 1:
            runs[-1][1] = r
        else:
            runs.append([r, r])
    return runs


def audit(d):
    cif = (glob.glob(os.path.join(d, "*.cif.gz")) or
           glob.glob(os.path.join(d, "*.cif")))[0]
    js = cif.replace(".cif.gz", ".json").replace(".cif", ".json")
    at = m195.read_cif_gz(cif)
    imap = json.load(open(js)).get("diffused_index_map", {})
    dc = m195.design_chain(at, imap)
    dat = [a for a in at if a[0] == dc]
    hab1 = np.array([a[4] for a in at if a[0] != dc])

    ca, by = {}, {}
    for a in dat:
        by.setdefault(a[1], []).append(a[4])
        if a[2] == "CA":
            ca[a[1]] = a[4]
    motif_out = {int(v[1:]) for v in imap.values() if v[0] == dc}
    runs = segments(motif_out)
    res = sorted(ca)

    # --- label every residue by segment -------------------------------------
    lab = {}
    for r in res:
        if r in motif_out:
            lab[r] = "motif"
        elif r < runs[0][0]:
            lab[r] = "lead"
        elif r > runs[-1][1]:
            lab[r] = "trail"
        else:
            k = sum(1 for a, b in runs if b < r)
            lab[r] = f"link{k}"

    # --- contacts -----------------------------------------------------------
    touch = set()
    for r in res:
        xyz = np.array(by[r])
        if len(hab1) and np.linalg.norm(
                hab1[None, :, :] - xyz[:, None, :], axis=2).min() < CONTACT:
            touch.add(r)

    per = {}
    for r in touch:
        per[lab[r]] = per.get(lab[r], 0) + 1
    tot = {}
    for r in res:
        tot[lab[r]] = tot.get(lab[r], 0) + 1

    # longest consecutive NON-MOTIF run touching HAB1 = the appendage
    best, cur = [], []
    for r in res:
        if r in touch and lab[r] != "motif":
            cur.append(r)
            if len(cur) > len(best):
                best = list(cur)
        else:
            cur = []

    # --- breaks -------------------------------------------------------------
    brk = []
    for i in range(len(res) - 1):
        a, b = res[i], res[i + 1]
        if b == a + 1:
            gap = float(np.linalg.norm(ca[b] - ca[a]))
            if gap > 4.5:
                brk.append((a, b, gap, lab[a], lab[b]))

    print(f"\n== {os.path.basename(d)}   design chain {dc}, {len(res)} residues")
    print(f"   motif runs (output numbering): "
          f"{', '.join(f'{a}-{b}' for a, b in runs)}")
    print(f"   HAB1 contacts by segment (contacting / total residues):")
    for k in ["lead"] + [f"link{i}" for i in range(1, len(runs))] + ["trail", "motif"]:
        if k in tot:
            n = per.get(k, 0)
            mark = "   <-- APPENDAGE" if k != "motif" and n >= 5 else ""
            print(f"      {k:<8} {n:>3} / {tot[k]:<4}{mark}")
    print(f"   TOTAL non-motif residues touching HAB1: "
          f"{sum(v for k, v in per.items() if k != 'motif')}"
          f"   (F4 only counted 'lead' within the first 25)")
    if best:
        print(f"   longest contiguous non-motif run touching HAB1: "
              f"{len(best)} residues, {best[0]}-{best[-1]} ({lab[best[0]]})")
    print(f"   chain breaks: {len(brk)}")
    for a, b, gap, la, lb in brk:
        print(f"      {a}->{b}  CA-CA {gap:5.1f} A   {la} | {lb}")
    return dict(name=os.path.basename(d), per=per, tot=tot,
                nonmotif_touch=sum(v for k, v in per.items() if k != "motif"),
                appendage=len(best), breaks=brk)


def main():
    ds = sys.argv[1:]
    if not ds:
        print(__doc__); return 1
    out = [audit(os.path.join(ROOT, d) if not os.path.isabs(d) else d) for d in ds]
    print("\n  SUMMARY")
    print(f"  {'design':<8}{'nonmotif->HAB1':>16}{'longest run':>13}{'breaks':>8}")
    for o in out:
        print(f"  {o['name']:<8}{o['nonmotif_touch']:>16}{o['appendage']:>13}"
              f"{len(o['breaks']):>8}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
