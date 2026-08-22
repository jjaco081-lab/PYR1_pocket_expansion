#!/usr/bin/env python
"""
86_ti_merge.py -- collapse the two protein copies into one TI-ready topology.

tleap builds a box holding BOTH endpoint proteins superimposed. tiMerge matches
their atoms BY COORDINATE and keeps one copy of everything the two states share,
leaving only the perturbed sidechains duplicated. It then prints the timask/scmask
strings pmemd needs, renumbered for the merged topology -- which is why the masks
must be read from its output rather than written by hand.

The softcore region is the WHOLE sidechain from CB outward, not just the atoms
that appear or disappear. Measured on the real topologies, 11 of the 15 atoms VAL
and ILE nominally share change partial charge, and their CG1 hydrogens do not sit
at matching coordinates (methyl in VAL, methylene in ILE). A common-core mask
would be quietly wrong; this one is merely slower.
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TI = os.path.join(ROOT, "data", "ti")
PARMED = "/bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/parmed"

NPROT = 179
SC = {
    "V81I": (79, "CB,HB,CG1,HG11,HG12,HG13,CG2,HG21,HG22,HG23",
                 "CB,HB,CG2,HG21,HG22,HG23,CG1,HG12,HG13,CD1,HD11,HD12,HD13"),
    "K59R": (59, "CB,HB2,HB3,CG,HG2,HG3,CD,HD2,HD3,CE,HE2,HE3,NZ,HZ1,HZ2,HZ3",
                 "CB,HB2,HB3,CG,HG2,HG3,CD,HD2,HD3,NE,HE,CZ,NH1,HH11,HH12,NH2,HH21,HH22"),
}


def audit(d, mut):
    """Measure how far the two copies' shared atoms really are apart.

    tiMerge matches by coordinate at a 0.01 A default and refuses if anything
    misses. Rather than raise the tolerance until it passes -- which would hide a
    genuine mismatch -- measure the discrepancy first and set the tolerance from
    it, asserting that what is being tolerated is a rebuilt hydrogen and not a
    displaced heavy atom. Returns (tol, report) or raises.
    """
    import numpy as np
    import parmed as pmd
    seq, sc1, _ = SC[mut]
    sc = set(sc1.split(","))
    p = pmd.load_file(os.path.join(d, "dual.prmtop"), os.path.join(d, "dual.inpcrd"))
    xyz = np.array(p.coordinates)
    r1, r2 = p.residues[:NPROT], p.residues[NPROT:2 * NPROT]
    bad = []
    for i, (a, b) in enumerate(zip(r1, r2)):
        A = {at.name: xyz[at.idx] for at in a.atoms}
        B = {at.name: xyz[at.idx] for at in b.atoms}
        for nm, va in A.items():
            if i + 1 == seq and nm in sc:
                continue
            if nm not in B:
                raise SystemExit(f"{d}: {nm} in resid {i+1} has no counterpart")
            dd = float(np.linalg.norm(va - B[nm]))
            if dd > 0.01:
                bad.append((dd, i + 1, nm))
    bad.sort(reverse=True)
    heavy = [x for x in bad if not x[2].startswith("H")]
    if heavy:
        raise SystemExit(f"{d}: HEAVY atoms displaced between copies: {heavy[:3]} "
                         "-- the copies are not the same structure, do not merge")
    mx = bad[0][0] if bad else 0.0
    if mx > 0.10:
        raise SystemExit(f"{d}: max mismatch {mx:.3f} A is too large to tolerate")
    tol = max(0.02, round(mx * 2, 3))
    return tol, f"{len(bad)} H rebuilt, max {mx:.3f} A -> tol {tol}"


def merge(mut, leg):
    d = os.path.join(TI, f"{mut}_{leg}")
    if not os.path.exists(os.path.join(d, "dual.prmtop")):
        return f"{mut}/{leg}: no dual.prmtop"
    if os.path.exists(os.path.join(d, "ti.prmtop")):
        return f"{mut}/{leg}: already merged"
    seq, sc1, sc2 = SC[mut]
    tol, rep = audit(d, mut)
    m1, m2 = f":1-{NPROT}", f":{NPROT+1}-{2*NPROT}"
    s1, s2 = f":{seq}@{sc1}", f":{seq+NPROT}@{sc2}"
    script = (f"tiMerge {m1} {m2} {s1} {s2}\n"
              f"outparm ti.prmtop ti.inpcrd\nquit\n").replace(
                  f"{s2}\n", f"{s2} tol {tol}\n")
    open(os.path.join(d, "merge.parmed"), "w").write(script)
    r = subprocess.run([PARMED, "-p", "dual.prmtop", "-c", "dual.inpcrd",
                        "-i", "merge.parmed"],
                       cwd=d, capture_output=True, text=True)
    log = r.stdout + r.stderr
    open(os.path.join(d, "merge.log"), "w").write(log)
    if not os.path.exists(os.path.join(d, "ti.prmtop")):
        err = [l for l in log.splitlines() if re.search(r"error|Error|ERROR", l)]
        return f"{mut}/{leg}: MERGE FAILED -- {err[:2] if err else log.splitlines()[-3:]}"
    # tiMerge prints the renumbered masks; they are the only correct ones
    masks = {}
    for k in ("timask1", "timask2", "scmask1", "scmask2"):
        m = re.search(rf"{k}\s*=\s*'([^']*)'", log) or re.search(rf"{k}\s*[:=]\s*(\S+)", log)
        if m:
            masks[k] = m.group(1).strip().rstrip(",")
    if len(masks) == 4:
        with open(os.path.join(d, "masks.txt"), "w") as fh:
            for k in ("timask1", "timask2", "scmask1", "scmask2"):
                fh.write(f"{k}={masks[k]}\n")
        return f"{mut}/{leg}: merged OK ({rep})"
    return f"{mut}/{leg}: merged but MASKS NOT PARSED -- read merge.log by hand"


if __name__ == "__main__":
    for mut in ("V81I", "K59R"):
        for leg in ("aba", "mandi", "apo"):
            print(" ", merge(mut, leg))
