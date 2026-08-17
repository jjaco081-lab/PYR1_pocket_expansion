#!/usr/bin/env python
"""
lib_structqc.py -- the quality checks every rebuilt structure must pass, in one place
so they cannot be skipped by whoever writes the next builder.

`lib_resnumber.verify_build` answers "is this the right protein, numbered right, with
no broken bonds". That is necessary and not sufficient: a build can satisfy all of it
and still have silently wrecked the crystal core, buried a clash, or inflated the
solvent box to the point that the MD costs twice as much for nothing. Those are the
three failures this module looks for.

None of these are hypothetical. The relax step is given a MoveMap that is supposed to
freeze everything except the rebuilt regions; `core_preservation` is what proves it
actually did, rather than trusting the flag.
"""
from __future__ import annotations

import numpy as np

from lib_resnumber import load_model, chain_residues, select_altloc

BB = ("N", "CA", "C")
CLASH_CUTOFF = 2.2      # A between heavy atoms of non-adjacent residues
BOX_BUFFER = 12.0       # A, matching lib_solvate's BUFFER


def _resmap(path, chain_id):
    model = load_model(path)
    chain = model[chain_id]
    select_altloc(chain)
    return chain_residues(chain)


def core_preservation(before_pdb, after_pdb, chain_id, rebuilt, tol=0.01):
    """Backbone displacement of every residue that was supposed to be held fixed.

    `rebuilt` is the set of residue numbers the builder deliberately allowed to move.
    Anything outside it that moves is a bug in the MoveMap, not a modelling choice.

    Returns (max_shift, n_moved, offenders). Residues absent from `before` (i.e. newly
    built) are skipped rather than counted as movement.
    """
    a = _resmap(before_pdb, chain_id)
    b = _resmap(after_pdb, chain_id)
    offenders = []
    worst = 0.0
    for i in sorted(a):
        if i in rebuilt or i not in b:
            continue
        d = max(float(np.linalg.norm(a[i][x].coord - b[i][x].coord))
                for x in BB if x in a[i] and x in b[i])
        worst = max(worst, d)
        if d > tol:
            offenders.append((i, a[i].get_resname(), round(d, 3)))
    return worst, len(offenders), offenders


def clashes(path, chain_ids, cutoff=CLASH_CUTOFF):
    """Heavy-atom contacts below `cutoff` between residues that are not sequence
    neighbours. Adjacent residues are excluded because their backbone atoms are
    legitimately close.

    Returns a list of (chain_a, res_a, chain_b, res_b, distance).
    """
    model = load_model(path)
    coords, tags = [], []
    for cid in chain_ids:
        chain = model[cid]
        select_altloc(chain)
        for num, res in chain_residues(chain).items():
            for atom in res:
                if atom.element == "H":
                    continue
                coords.append(atom.coord)
                tags.append((cid, num))
    crd = np.asarray(coords)
    d = np.linalg.norm(crd[:, None, :] - crd[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    out = {}
    for i, j in np.argwhere(d < cutoff):
        if i >= j:
            continue
        ca, na = tags[i]
        cb, nb = tags[j]
        if ca == cb and abs(na - nb) <= 1:
            continue
        key = (ca, na, cb, nb)
        out[key] = min(out.get(key, np.inf), round(float(d[i, j]), 2))
    return sorted((k[0], k[1], k[2], k[3], v) for k, v in out.items())


def solvent_cost(path, chain_ids, exclude=None, buffer=BOX_BUFFER):
    """How much periodic box the structure demands, with and without `exclude`.

    Used to price a disordered segment before committing to it: an extended 10-residue
    tail can double the water count, and every extra water is paid for on every one of
    the 300 ns. Returns (volume_with, volume_without, fractional_increase). Volumes are
    the bounding box plus `buffer` on each side, in nm^3 -- a proxy for the solvateoct
    cell, not the cell itself, and only meaningful as a ratio.
    """
    exclude = set(exclude or ())
    model = load_model(path)
    keep, drop = [], []
    for cid in chain_ids:
        chain = model[cid]
        select_altloc(chain)
        for num, res in chain_residues(chain).items():
            for atom in res:
                (drop if (cid, num) in exclude or num in exclude else keep).append(atom.coord)
    if not drop:
        return None, None, 0.0
    full = np.asarray(keep + drop)
    part = np.asarray(keep)
    v_full = float(np.prod(full.max(0) - full.min(0) + 2 * buffer) / 1000.0)
    v_part = float(np.prod(part.max(0) - part.min(0) + 2 * buffer) / 1000.0)
    return v_full, v_part, v_full / v_part - 1.0


def bond_geometry(path, chain_id):
    """Every C-N peptide distance, so junction quality is visible rather than assumed.

    Returns (min, max, mean) over all consecutive pairs. Ideal is ~1.33 A; a value
    above 1.5 is a break and `verify_build` will already have refused it, but a build
    that squeaks in at 1.48 is worth seeing.
    """
    rm = _resmap(path, chain_id)
    nums = sorted(rm)
    d = []
    for a, b in zip(nums, nums[1:]):
        if "C" in rm[a] and "N" in rm[b]:
            d.append(float(np.linalg.norm(rm[a]["C"].coord - rm[b]["N"].coord)))
    d = np.asarray(d)
    return float(d.min()), float(d.max()), float(d.mean())


def report(before_pdb, after_pdb, chain_ids, rebuilt, exclude_for_cost=None,
           log=print):
    """Run every check and print a verdict. Returns a dict for build_report.json."""
    primary = chain_ids[0]
    worst, n_moved, offenders = core_preservation(before_pdb, after_pdb, primary, rebuilt)
    log(f"  core preservation: {n_moved} held-fixed residue(s) moved >0.01 A, "
        f"max {worst:.3f} A" + (f"  {offenders[:5]}" if offenders else ""))
    if n_moved:
        raise AssertionError(
            f"{n_moved} residues moved that the MoveMap was supposed to freeze: "
            f"{offenders[:10]} -- the crystal core is not preserved")

    cl = clashes(after_pdb, chain_ids)
    log(f"  steric clashes <{CLASH_CUTOFF} A between non-adjacent residues: {len(cl)}"
        + (f"  {cl[:5]}" if cl else ""))
    if cl:
        raise AssertionError(f"{len(cl)} steric clash(es) introduced: {cl[:10]}")

    lo, hi, mean = bond_geometry(after_pdb, primary)
    log(f"  peptide C-N bonds: min {lo:.3f}, max {hi:.3f}, mean {mean:.3f} A")

    v_full, v_part, frac = solvent_cost(after_pdb, chain_ids, exclude_for_cost)
    if v_full is not None:
        log(f"  solvent box: {v_part:.0f} -> {v_full:.0f} nm^3 ({frac:+.0%}) "
            "from the rebuilt disordered segment")

    return {"core_max_shift_A": round(worst, 4), "core_residues_moved": n_moved,
            "clashes": len(cl), "cn_min": round(lo, 3), "cn_max": round(hi, 3),
            "cn_mean": round(mean, 3),
            "box_nm3_with": v_full, "box_nm3_without": v_part,
            "box_increase_frac": round(frac, 4)}
