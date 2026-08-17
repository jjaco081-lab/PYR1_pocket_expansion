#!/usr/bin/env python
"""
lib_rosetta.py -- the one thing about FastRelax that this project must never get wrong
again, plus the assertion that proves it.

THE TRAP
--------
`FastRelax.set_movemap(mm)` controls MINIMISATION degrees of freedom ONLY. It does not
touch the packer. FastRelax repacks side chains through a TaskFactory, and if you do
not give it one, **it repacks every residue in the pose** no matter what the MoveMap
says.

That single fact produced both structural defects in the 191-aa rebuild, and both got
past a backbone check reporting 0.000 A displacement:

  * protomer B of the dimer was relaxed alone. R116 and M158 -- nowhere near the
    rebuilt regions in the MoveMap -- were repacked anyway, into the volume protomer A
    occupies. A87-B116 came out at 0.36 A where the crystal has 3.09 A.
  * the closed monomer was relaxed with no ligand in the pose (tleap adds ABA later).
    K59 was explicitly placed in a "frozen" set and repacked regardless, collapsing
    into the empty cavity: its closest approach to ABA's crystal position went from
    2.85 A to 1.68 A.

Backbone RMSD saw neither, because it watches N, CA and C.

THE FIX
-------
`restrict_packing` builds a TaskFactory that (a) forbids design outright and (b)
prevents repacking on everything outside an explicit allowed set. `assert_frozen` then
checks ALL HEAVY ATOMS of the residues that were supposed to stay put, so a regression
here cannot be silent again.

Use both. The TaskFactory is what makes the freeze real; the assertion is what proves
it stayed real.
"""
from __future__ import annotations

import numpy as np


def restrict_packing(pose, allowed):
    """TaskFactory that repacks ONLY `allowed` (1-based pose indices) and designs nothing.

    Everything else keeps its input rotamer exactly. Pass the result to
    `FastRelax.set_task_factory()`; passing only a MoveMap is not enough.
    """
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import (
        RestrictToRepacking, OperateOnResidueSubset, PreventRepackingRLT)
    from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector

    allowed = sorted(set(allowed))
    if not allowed:
        raise ValueError("restrict_packing called with an empty allowed set")
    sel = ResidueIndexSelector(",".join(map(str, allowed)))
    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())          # repack, never design
    # flip_subset=True -> the operation applies to residues NOT in `sel`
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), sel, True))

    # prove it before handing it over: build the task and count what may move
    task = tf.create_task_and_apply_taskoperations(pose)
    packable = [i for i in range(1, pose.total_residue() + 1)
                if task.residue_task(i).being_packed()]
    designable = [i for i in range(1, pose.total_residue() + 1)
                  if task.residue_task(i).being_designed()]
    if designable:
        raise AssertionError(f"task would DESIGN residues {designable[:10]}")
    if packable != allowed:
        raise AssertionError(
            f"packer task does not match the allowed set: "
            f"{len(packable)} packable vs {len(allowed)} requested; "
            f"unexpected {sorted(set(packable) - set(allowed))[:10]}, "
            f"missing {sorted(set(allowed) - set(packable))[:10]}")
    return tf, packable


#: Tolerance for "did not move", in A, chosen from the measured distribution rather
#: than picked to make a test pass. Cartesian minimisation relaxes bonded neighbours of
#: movable residues slightly, so residues just outside the movable set pick up a little
#: geometric leakage. Measured over both PYR1 builds, every residue outside the
#: rebuilt/repacked set is either EXACTLY 0.000 A (344 of 350) or below 0.07 A, while
#: the smallest legitimately rebuilt residue moves 2.19 A. The threshold sits in a gap
#: spanning a factor of ~32, so it is not delicately placed -- and it is still an order
#: of magnitude tighter than the defects it exists to catch (K59 moved 1.17 A, the
#: interface rotamers 2.5-3.0 A).
FROZEN_TOL = 0.1


def assert_frozen(before_resmap, after_resmap, residues, label, tol=FROZEN_TOL):
    """Assert that every HEAVY ATOM of `residues` is unmoved between two structures.

    Deliberately all-atom, not backbone: the defects this guards against were
    side-chain-only and a backbone check reported 0.000 A while they happened.
    Hydrogens are ignored because the input crystal structures have none.

    Returns (n_checked, max_drift) so the margin is reported and a regression that
    creeps toward the threshold is visible before it crosses it.
    """
    moved, worst = [], 0.0
    for num in sorted(residues):
        if num not in before_resmap or num not in after_resmap:
            continue
        a = {x.get_id(): x.coord for x in before_resmap[num] if x.element != "H"}
        b = {x.get_id(): x.coord for x in after_resmap[num] if x.element != "H"}
        shared = set(a) & set(b)
        if not shared:
            continue
        d = max(float(np.linalg.norm(a[k] - b[k])) for k in shared)
        worst = max(worst, d)
        if d > tol:
            moved.append((num, before_resmap[num].get_resname(), round(d, 3)))
    if moved:
        raise AssertionError(
            f"{label}: {len(moved)} residue(s) that were supposed to be frozen moved "
            f"(heavy atoms, tol {tol} A): {moved[:10]}")
    return len(residues), worst
