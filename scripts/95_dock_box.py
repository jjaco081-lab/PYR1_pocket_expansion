#!/usr/bin/env python
r"""95_dock_box.py -- docking box centre from pocket-lining PROTEIN atoms only.

Placing the box on the ligand would defeat the purpose of the docked arm, which
exists to measure what happens when no complex structure is available. Sequential
numbering (the pocket set below is sequential, matching the 179-residue frame).
"""
import sys
import numpy as np

POCKET = {59, 61, 62, 79, 81, 83, 85, 87, 88, 89, 90, 92, 106, 108, 111, 113,
          115, 118, 119, 139, 141, 155, 157, 159, 161, 162}
BB = {"N", "CA", "C", "O"}

X = [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
     for l in open(sys.argv[1])
     if l.startswith("ATOM") and int(l[22:26]) in POCKET
     and l[12:16].strip() not in BB]
if len(X) < 50:
    raise SystemExit(f"only {len(X)} pocket atoms found -- numbering may be wrong")
c = np.array(X).mean(0)
print(f"{c[0]:.3f} {c[1]:.3f} {c[2]:.3f}")
