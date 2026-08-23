#!/usr/bin/env python
r"""
93_quad_merge.py -- tiMerge for the FOUR-mutation PYR1^MANDI system.

Differs from 86 (single mutations) in two ways that matter:

  1. The softcore region spans FOUR residues, so the mask is written
     generically as ":59,79,106,157&!@N,H,CA,HA,C,O" -- "these residues, minus
     the backbone" -- rather than by enumerating sidechain atom names for every
     residue pair. Enumerating 8 different sidechains by hand is exactly the
     kind of transcription that produced the masks bug 87 was written to avoid.
  2. The audit skips softcore atoms for all four residues, not one.

Everything else follows 86: measure the coordinate mismatch between the two
copies FIRST and set the tolerance from the measurement, asserting that what is
tolerated is a rebuilt hydrogen and not a displaced heavy atom.
"""
import os
import re
import subprocess

import numpy as np
import parmed as pmd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Q = os.path.join(ROOT, "data", "ti_quad")
PARMED = "/bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/parmed"

NPROT = 179
SEQ = [59, 79, 106, 157]          # sequential indices, asserted in 91
# Only the HEAVY backbone stays common. The amide H and HA are put INSIDE the
# softcore region on the mutated residues: tleap rebuilds them against CB, and
# the crystal arm's CB sits 0.43-1.08 A from the WT CB, so as common atoms they
# missed their partners by up to 0.37 A and tiMerge refused to match them. As
# perturbed atoms they need no match, and being bonded to common N/CA is no
# different from CB, which was always in the region.
BACKBONE = {"N", "CA", "C", "O", "OXT"}
LEGS = ["aba", "mandi_xtal", "mandi_dock", "apo"]


def audit(d):
    p = pmd.load_file(os.path.join(d, "dual.prmtop"), os.path.join(d, "dual.inpcrd"))
    xyz = np.array(p.coordinates)
    r1, r2 = p.residues[:NPROT], p.residues[NPROT:2 * NPROT]
    if len(r1) != NPROT or len(r2) != NPROT:
        raise SystemExit(f"{d}: expected two {NPROT}-residue copies, got "
                         f"{len(r1)}/{len(r2)}")
    bad = []
    for i, (a, b) in enumerate(zip(r1, r2), start=1):
        A = {at.name: xyz[at.idx] for at in a.atoms}
        B = {at.name: xyz[at.idx] for at in b.atoms}
        mutated = i in SEQ
        if mutated and a.name == b.name:
            raise SystemExit(f"{d}: residue {i} is {a.name} in BOTH copies "
                             "-- the mutation was not applied")
        for nm, va in A.items():
            # On a mutated residue only the SIDECHAIN is softcore. Its BACKBONE
            # is nonsoftcore and tiMerge still has to coordinate-match it, so
            # skipping the whole residue here (as the first version did) hid the
            # one atom that actually fails: tleap rebuilds HA from the residue
            # template, and VAL's and ILE's differ by 0.024 A.
            if mutated and nm not in BACKBONE:
                continue
            if nm not in B:
                raise SystemExit(f"{d}: {nm} in resid {i} has no counterpart")
            dd = float(np.linalg.norm(va - B[nm]))
            if dd > 0.01:
                bad.append((dd, i, nm))
    bad.sort(reverse=True)
    heavy = [x for x in bad if not x[2].startswith("H")]
    if heavy:
        raise SystemExit(f"{d}: HEAVY atoms displaced between copies: {heavy[:3]}")
    mx = bad[0][0] if bad else 0.0
    if mx > 0.10:
        raise SystemExit(f"{d}: max mismatch {mx:.3f} A too large to tolerate")
    muts = [(i, r1[i - 1].name, r2[i - 1].name) for i in SEQ]
    return max(0.02, round(mx * 2, 3)), f"{len(bad)} H rebuilt, max {mx:.3f} A", muts


def merge(leg):
    d = os.path.join(Q, leg)
    if not os.path.exists(os.path.join(d, "dual.prmtop")):
        return f"{leg}: no dual.prmtop"
    if os.path.exists(os.path.join(d, "ti.prmtop")):
        return f"{leg}: already merged"
    tol, rep, muts = audit(d)
    bb = ",".join(sorted(BACKBONE))
    s1 = ":" + ",".join(str(s) for s in SEQ) + f"&!@{bb}"
    s2 = ":" + ",".join(str(s + NPROT) for s in SEQ) + f"&!@{bb}"
    script = (f"tiMerge :1-{NPROT} :{NPROT+1}-{2*NPROT} {s1} {s2} tol {tol}\n"
              "outparm ti.prmtop ti.inpcrd\nquit\n")
    open(os.path.join(d, "merge.parmed"), "w").write(script)
    # `module load amber/22` exports a PYTHONPATH that shadows this conda ParmEd
    # 4.3 with amber's own 3.4.1. Inheriting it silently runs a different tiMerge
    # than the one tested here, so strip it for the child.
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    r = subprocess.run([PARMED, "-p", "dual.prmtop", "-c", "dual.inpcrd",
                        "-i", "merge.parmed"], cwd=d, capture_output=True,
                       text=True, env=env)
    log = r.stdout + r.stderr
    open(os.path.join(d, "merge.log"), "w").write(log)
    if not os.path.exists(os.path.join(d, "ti.prmtop")):
        err = [l for l in log.splitlines() if re.search(r"[Ee]rror|ERROR", l)]
        return f"{leg}: MERGE FAILED -- {err[:3] if err else log.splitlines()[-4:]}"
    masks = {}
    for k in ("timask1", "timask2", "scmask1", "scmask2"):
        m = re.search(rf"{k}\s*=\s*'([^']*)'", log) or re.search(rf"{k}\s*[:=]\s*(\S+)", log)
        if m:
            masks[k] = m.group(1).strip().rstrip(",")
    if len(masks) != 4:
        return f"{leg}: merged but MASKS NOT PARSED -- read merge.log by hand"
    with open(os.path.join(d, "masks.txt"), "w") as fh:
        for k in ("timask1", "timask2", "scmask1", "scmask2"):
            fh.write(f"{k}={masks[k]}\n")
    mt = " ".join(f"{a}{i}{b}" for i, a, b in muts)
    return f"{leg}: merged OK ({rep}; {mt})"


if __name__ == "__main__":
    fail = 0
    for leg in LEGS:
        msg = merge(leg)
        print("   ", msg)
        if "FAILED" in msg or "NOT PARSED" in msg:
            fail = 1
    raise SystemExit(fail)
