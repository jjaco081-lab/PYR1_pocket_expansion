#!/usr/bin/env python
"""
lib_sterics.py -- one steric-clash rule, shared, that does not score hydrogen bonds
as clashes.

THE BUG THIS EXISTS TO KILL
---------------------------
A hard-sphere van der Waals test treats every close approach as a clash. Hydrogen
bonds are close approaches BY DEFINITION: a donor-acceptor pair sits at 2.6-3.2 A
while the Bondi radii sum to 1.55 + 1.52 = 3.07, so a perfectly ordinary N-H...O
bond registers as ~0.4 A of "overlap".

That produced three separate wrong answers in this project:

  * README 31e -- R79, V83 and H115 failed the wild-type positive control: the
    method said the crystal residue could not be placed at its own position.
  * README 33b -- the crystal Arg59 from 4WVO, which makes NE-O2 at 2.64 A, came
    back with a 1.02 A "overlap" and was rejected as infeasible.
  * and it biases every menu against polar residues, silently, in the direction
    that makes a steric method look more ligand-specific than it is.

Confirmed independently from the literature (README 34b): R79 hydrogen bonds the
main-chain carbonyl of F52 at >90% occupancy in both apo and holo PYR1 and never
touches the ligand. Its "clash" was a real hydrogen bond all along.

THE RULE
--------
A pair is exempt from the hard-sphere test when one atom can donate and the other
can accept. Exempt pairs are clashes only below `HBOND_MIN`, which is shorter than
any real hydrogen bond. Everything else keeps the ordinary overlap rule.

Acceptor-acceptor and donor-donor pairs are NOT exempt -- two carbonyl oxygens at
2.6 A really are clashing. That is why donor and acceptor are tracked separately
instead of collapsing to "is polar".
"""
from __future__ import annotations

import numpy as np

#: Bondi radii. Hydrogen is present so that including hydrogens cannot silently
#: give them a carbon-sized 1.70 A (README 33c).
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "H": 1.20,
       "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}

#: Shorter than any real hydrogen bond. Donor-acceptor pairs closer than this are
#: genuinely clashing; between this and touching they are bonding, not colliding.
#: 2.5 A is below the 2.6 A low end of the observed N/O-N/O hydrogen-bond range and
#: below the 2.64 A of the 4WVO Arg59 contact this rule exists to admit.
HBOND_MIN = 2.5

DONOR, ACCEPTOR = 1, 2


def atom_arrays(residue, side_chain_only=False, include_h=False, backbone=("N", "CA", "C", "O")):
    """(xyz, radii, flags) for one Rosetta residue; flags are DONOR|ACCEPTOR bits."""
    xyz, rad, flag = [], [], []
    for a in range(1, residue.natoms() + 1):
        t = residue.atom_type(a)
        e = t.element().strip().upper()
        nm = residue.atom_name(a).strip()
        if not include_h and e == "H":
            continue
        if side_chain_only and nm in backbone:
            continue
        v = residue.xyz(a)
        xyz.append([v.x, v.y, v.z])
        rad.append(VDW.get(e, 1.7))
        flag.append((DONOR if t.is_donor() else 0) | (ACCEPTOR if t.is_acceptor() else 0))
    return (np.array(xyz).reshape(-1, 3), np.array(rad), np.array(flag, dtype=int))


def _pair_matrices(xyzA, radA, flagA, xyzB, radB, flagB):
    d = np.linalg.norm(xyzA[:, None, :] - xyzB[None, :, :], axis=2)
    ov = (radA[:, None] + radB[None, :]) - d
    # exempt where one side can donate and the other can accept
    da = ((flagA[:, None] & DONOR) > 0) & ((flagB[None, :] & ACCEPTOR) > 0)
    ad = ((flagA[:, None] & ACCEPTOR) > 0) & ((flagB[None, :] & DONOR) > 0)
    return d, ov, (da | ad)


def max_overlap(xyzA, radA, flagA, xyzB, radB, flagB, hbond_min=HBOND_MIN):
    """Largest steric overlap, with hydrogen-bondable pairs exempted.

    An exempt pair contributes only if it is closer than `hbond_min`, and then it
    contributes its true overlap -- so a genuinely collapsed contact is still caught.
    """
    if len(xyzA) == 0 or len(xyzB) == 0:
        return 0.0
    d, ov, exempt = _pair_matrices(xyzA, radA, flagA, xyzB, radB, flagB)
    ov = np.where(exempt & (d >= hbond_min), -np.inf, ov)
    m = ov.max()
    return float(m) if np.isfinite(m) else 0.0


def sum_overlap(xyzA, radA, flagA, xyzB, radB, flagB, hbond_min=HBOND_MIN):
    """Summed positive overlap, same exemption. Additive over atoms by construction."""
    if len(xyzA) == 0 or len(xyzB) == 0:
        return 0.0
    d, ov, exempt = _pair_matrices(xyzA, radA, flagA, xyzB, radB, flagB)
    ov = np.where(exempt & (d >= hbond_min), 0.0, ov)
    return float(np.clip(ov, 0, None).sum())


def good_contacts(xyzA, radA, flagA, xyzB, radB, flagB, far=1.0, near=0.25):
    """Pairs touching without interpenetrating (README 32b).

    Hydrogen bonds count as contacts rather than being excluded for being too close,
    which is the same correction applied from the other side.
    """
    if len(xyzA) == 0 or len(xyzB) == 0:
        return 0
    d, ov, exempt = _pair_matrices(xyzA, radA, flagA, xyzB, radB, flagB)
    ok = (ov >= -far) & (ov <= near)
    ok |= exempt & (d >= HBOND_MIN) & (ov >= -far)
    return int(ok.sum())


def hbonds(xyzA, flagA, xyzB, flagB, lo=HBOND_MIN, hi=3.5):
    """Count donor-acceptor pairs in hydrogen-bonding range (heavy-atom distance)."""
    if len(xyzA) == 0 or len(xyzB) == 0:
        return 0
    d = np.linalg.norm(xyzA[:, None, :] - xyzB[None, :, :], axis=2)
    da = ((flagA[:, None] & DONOR) > 0) & ((flagB[None, :] & ACCEPTOR) > 0)
    ad = ((flagA[:, None] & ACCEPTOR) > 0) & ((flagB[None, :] & DONOR) > 0)
    return int(((da | ad) & (d >= lo) & (d <= hi)).sum())
