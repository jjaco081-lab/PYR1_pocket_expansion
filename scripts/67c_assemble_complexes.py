#!/usr/bin/env python
"""
67c_assemble_complexes.py -- assemble the rebuilt 191-residue protomers into the
multi-chain proteins the MD systems actually need, and prove the assembly is sound.

Run with the tier1_analysis env (only for consistency; this script needs no PyRosetta):
    .../conda_envs/pyr1_docking/bin/python scripts/67c_assemble_complexes.py

WHY A SEPARATE STEP
-------------------
Script 67 builds one protomer at a time, each in its own crystal frame. That is
deliberate -- the two protomers of 3K3K are crystallographically independent and
forcing them to be identical would invent a symmetry the data does not have. But it
means the rebuilt regions of each protomer were sampled WITHOUT the other one
present, so nothing yet guarantees that protomer A's new C-terminal tail does not run
straight through protomer B.

That is exactly the failure this script exists to catch. It is cheap to check and
expensive to miss: tleap will happily build a topology with two chains interpenetrating
and the first minimisation step will either explode or quietly push the fold apart.

WHY THE CRYSTAL FRAME STILL HOLDS
---------------------------------
Script 67 freezes every residue it did not rebuild -- verified at 0.000 A backbone
displacement over 175 of 176 held-fixed residues. So the protomers are still in their
deposited relative orientation, the dimer interface is the crystallographic one, and
ABA's crystal coordinates are still valid against the rebuilt pocket. The ligand
therefore does NOT need re-docking, and this script asserts that rather than assuming
it: `check_ligand_contacts` recomputes the pocket contacts and compares them with the
same contacts measured on the original CIF.

WHAT THE CLASH CHECK CAUGHT, AND WHAT IT TURNED OUT TO MEAN
-----------------------------------------------------------
The first version of this script found three inter-chain clashes, the worst being
A87 LEU against B116 ARG at **0.36 A** -- heavy atoms essentially superimposed, where
the raw crystal has 3.09 A. So the clash was ours, not the PDB's.

The cause was NOT that the protomers were built independently, which is what this
script was written to catch. It was that `FastRelax.set_movemap()` restricts
MINIMISATION only: FastRelax repacks through a TaskFactory and, given none, repacks
every residue in the pose. R116 and M158 were never in protomer B's MoveMap and were
repacked anyway -- into the volume protomer A occupies. See `lib_rosetta.py`.

With the packer properly restricted, the assembled dimer now comes out of script 67
with **zero** inter-chain clashes before any refinement here. Building the protomers
separately was never the problem.

The in-context refinement is kept regardless. It is cheap, it is the physically right
place to relax a region that sits in an interface, and it means this script does not
silently depend on script 67 having got the packer right. The clash check remains a
gate, not a warning.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
from Bio.PDB import PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_structqc  # noqa: E402
from lib_resnumber import (  # noqa: E402
    PYR1_LANDMARKS, load_model, select_altloc, chain_residues, verify_build,
    write_residue_map,
)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
STRUCT = os.path.join(ROOT, "data", "structures_191")
CONTACT_CUTOFF = 4.5     # A, heavy atom, matching how contacts are counted elsewhere


def log(msg):
    print(msg, flush=True)


def combine(parts, out_path):
    """Write several (chain_id, source_pdb, source_chain) into one multi-chain PDB.

    Coordinates are copied verbatim -- no superposition, no re-centring. The whole
    point is to preserve the crystal frame the protomers were built in.
    """
    st = Structure("complex")
    mo = Model(0)
    st.add(mo)
    for new_chain, src_pdb, src_chain in parts:
        model = load_model(src_pdb)
        chain = model[src_chain]
        select_altloc(chain)
        ch = Chain(new_chain)
        mo.add(ch)
        for num, res in chain_residues(chain).items():
            r = res.copy()
            r.id = (" ", num, " ")
            r.detach_parent()
            ch.add(r)
    io = PDBIO()
    io.set_structure(st)
    io.save(out_path)
    return out_path


def refine_in_context(pdb_path, rebuilt_by_chain, chain_order, seed, repeats=2,
                      protected_by_chain=None):
    """Relax the rebuilt regions of an assembled complex with every chain present.

    `rebuilt_by_chain` maps chain id -> set of residue numbers script 67 was allowed to
    move in that protomer. Only those, plus their repack shell, are movable here; every
    other residue stays exactly where the crystal put it, and that is asserted
    afterwards by lib_structqc.core_preservation.

    Pose numbering is derived from the chains rather than assumed: each protomer is a
    verified 1..191, but hard-coding 191 here would be exactly the kind of assumption
    that has cost this project time before.
    """
    import pyrosetta
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from lib_rosetta import restrict_packing

    protected_by_chain = protected_by_chain or {}
    pose = pyrosetta.pose_from_pdb(pdb_path)
    info = pose.pdb_info()

    # native (chain, resnum) -> pose index, read back from the pose itself
    movable = []
    for i in range(1, pose.total_residue() + 1):
        cid, num = info.chain(i), info.number(i)
        if num in rebuilt_by_chain.get(cid, ()):
            movable.append(i)
    expected = sum(len(v) for v in rebuilt_by_chain.values())
    assert len(movable) == expected, \
        (f"found {len(movable)} rebuilt residues in the pose but the build report "
         f"lists {expected}; chain/numbering mapping is wrong")

    # ligand-lining residues, resolved into pose indices the same way
    protected = [i for i in range(1, pose.total_residue() + 1)
                 if info.number(i) in protected_by_chain.get(info.chain(i), ())]

    sel = ResidueIndexSelector(",".join(map(str, movable)))
    nbr = NeighborhoodResidueSelector(sel, 6.0, True)
    mv, pk = sel.apply(pose), nbr.apply(pose)
    mm = MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    for i in range(1, pose.total_residue() + 1):
        if i in protected:
            continue
        if mv[i]:
            mm.set_bb(i, True)
            mm.set_chi(i, True)
        elif pk[i]:
            mm.set_chi(i, True)
    n_pk = sum(1 for i in range(1, pose.total_residue() + 1)
               if pk[i] and i not in protected)
    log(f"  refining {len(movable)} rebuilt residues and repacking {n_pk} neighbours, "
        "this time with every chain present")

    # ⚠ Same trap as script 67: set_movemap governs minimisation only. Without a
    # TaskFactory FastRelax repacks the entire complex, which is exactly how protomer
    # B's R116 ended up inside protomer A in the first place. See lib_rosetta.
    allowed_pack = [i for i in range(1, pose.total_residue() + 1)
                    if (mv[i] or pk[i]) and i not in protected]
    tf, packable = restrict_packing(pose, allowed_pack)
    log(f"  packer restricted to {len(packable)} of {pose.total_residue()} residues")

    sfxn = pyrosetta.create_score_function("ref2015_cart")
    relax = FastRelax(sfxn, repeats)
    relax.cartesian(True)
    relax.set_movemap(mm)
    relax.set_task_factory(tf)
    relax.min_type("lbfgs_armijo_nonmonotone")
    relax.apply(pose)
    pose.dump_pdb(pdb_path)
    log(f"  complex score after in-context refinement: {sfxn(pose):.1f} REU")
    return {"n_movable": len(movable), "n_repacked": n_pk,
            "n_packable": len(packable), "n_protected": len(protected),
            "score_REU": float(sfxn(pose))}


def ligand_contacts(model, lig_resname, protein_chain, cutoff=CONTACT_CUTOFF):
    """{residue number: min heavy-atom distance} for every residue contacting the ligand."""
    lig = [r for ch in model for r in ch
           if r.get_resname().upper() == lig_resname.upper()]
    if not lig:
        return None
    lc = np.array([a.coord for r in lig for a in r if a.element != "H"])
    chain = model[protein_chain]
    select_altloc(chain)
    out = {}
    for num, res in chain_residues(chain).items():
        pc = np.array([a.coord for a in res if a.element != "H"])
        d = float(np.min(np.linalg.norm(pc[:, None, :] - lc[None, :, :], axis=-1)))
        if d <= cutoff:
            out[num] = round(d, 2)
    return out


def check_ligand_contacts(built_pdb, built_chain, ref_cif, ref_chain, lig, tol=0.05):
    """Assert the rebuilt receptor still presents the SAME pocket to the crystal ligand.

    This is the check that licenses reusing A8S.mol2's crystal pose instead of
    re-docking. If script 67 had moved any pocket residue, the contact set or the
    distances would shift and this would say so.
    """
    ref = ligand_contacts(load_model(ref_cif), lig, ref_chain)
    if ref is None:
        raise AssertionError(f"{ref_cif}: no {lig} found to compare against")
    # rebuild the same measurement against the built receptor, using the reference
    # ligand coordinates (the built file holds protein only)
    merged = load_model(ref_cif)
    built = load_model(built_pdb)[built_chain]
    select_altloc(built)
    built_map = chain_residues(built)
    lig_res = [r for ch in merged for r in ch if r.get_resname().upper() == lig.upper()]
    lc = np.array([a.coord for r in lig_res for a in r if a.element != "H"])

    now = {}
    for num, res in built_map.items():
        pc = np.array([a.coord for a in res if a.element != "H"])
        d = float(np.min(np.linalg.norm(pc[:, None, :] - lc[None, :, :], axis=-1)))
        if d <= CONTACT_CUTOFF:
            now[num] = round(d, 2)

    gained = sorted(set(now) - set(ref))
    lost = sorted(set(ref) - set(now))
    moved = {n: (ref[n], now[n]) for n in set(ref) & set(now)
             if abs(ref[n] - now[n]) > tol}
    log(f"  ligand {lig}: {len(ref)} contacting residues in the crystal, "
        f"{len(now)} against the rebuilt receptor")
    if gained or lost or moved:
        raise AssertionError(
            f"pocket changed during the rebuild -- gained {gained}, lost {lost}, "
            f"shifted {moved}. The crystal ligand pose is no longer valid here; "
            "either re-dock or find what moved.")
    log(f"  pocket is unchanged to within {tol} A -- the crystal {lig} pose stays valid")
    return {"n_contacts": len(now), "residues": sorted(now)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--structdir", default=STRUCT)
    ap.add_argument("--outdir", default=STRUCT)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Same seeding as script 67: the in-context refinement is FastRelax, i.e. Monte
    # Carlo, so without -constant_seed the assembled dimer would differ run to run.
    with open(os.path.join(args.structdir, "build_report.json")) as fh:
        seed = json.load(fh)["seed"]
    import pyrosetta
    pyrosetta.init(f"-mute all -ex1 -ex2aro -detect_disulf false "
                   f"-run:constant_seed -run:jran {seed}")

    report = {"assembled": datetime.now(timezone.utc).isoformat(), "seed": seed,
              "builder": "67c_assemble_complexes.py", "complexes": {}}

    # ---- the open dimer (S3) ----------------------------------------------------
    log("\n== assembling pyr1_open_dimer_191 (S3) ==")
    parts = [("A", os.path.join(args.structdir, "pyr1_open_191.pdb"), "A"),
             ("B", os.path.join(args.structdir, "pyr1_open_191_protomerB.pdb"), "A")]
    dimer = combine(parts, os.path.join(args.outdir, "pyr1_open_dimer_191.pdb"))
    for cid in ("A", "B"):
        v = verify_build(dimer, cid, label=f"dimer chain {cid}")
        log(f"  chain {cid}: {v['n_residues']} residues, md5 {v['sequence_md5']}, "
            "0 breaks")

    before = lib_structqc.clashes(dimer, ["A", "B"])
    log(f"  inter-chain clashes before refinement: "
        f"{len([c for c in before if c[0] != c[2]])}")

    # which residues script 67 was allowed to move in each protomer
    with open(os.path.join(args.structdir, "build_report.json")) as fh:
        br = json.load(fh)
    rebuilt = {"A": set(br["systems"]["pyr1_open_191"]["relaxed_residues"]),
               "B": set(br["systems"]["pyr1_open_191_protomerB"]["relaxed_residues"])}
    # 3K3K chain B carries ABA (§28g), so its pocket must stay frozen here too
    protected = {"A": set(br["systems"]["pyr1_open_191"]["frozen_ligand_shell"]),
                 "B": set(br["systems"]["pyr1_open_191_protomerB"]["frozen_ligand_shell"])}
    log(f"  ligand-lining residues held frozen: chain A {len(protected['A'])}, "
        f"chain B {len(protected['B'])}")
    ref = refine_in_context(dimer, rebuilt, ["A", "B"], seed=br["seed"],
                            protected_by_chain=protected)

    for cid in ("A", "B"):
        verify_build(dimer, cid, label=f"dimer chain {cid} (post-refinement)")
    after = lib_structqc.clashes(dimer, ["A", "B"])
    inter = [c for c in after if c[0] != c[2]]
    log(f"  inter-chain clashes after refinement: {len(inter)}"
        + (f"  {inter[:8]}" if inter else ""))
    if inter:
        raise AssertionError(
            f"the two protomers still interpenetrate at {inter[:10]} after refining in "
            "context. Do not build MD on this -- resample the tails or widen the "
            "movable set.")
    rm = os.path.join(args.outdir, "pyr1_open_dimer_191_residue_map.json")
    write_residue_map(dimer, rm, ["A", "B"])
    report["complexes"]["pyr1_open_dimer_191"] = {
        "pdb": dimer, "chains": ["A", "B"],
        "clashes_before_refinement": len([c for c in before if c[0] != c[2]]),
        "clashes_after_refinement": len(inter), "refinement": ref,
        "residue_map": rm}
    log(f"  wrote {os.path.relpath(dimer, ROOT)}")

    # ---- the closed monomer, checked against the crystal ABA (S2, S6-S9) ---------
    log("\n== checking pyr1_closed_191 against the crystal ABA (S2, S6-S9) ==")
    closed = os.path.join(args.structdir, "pyr1_closed_191.pdb")
    lig = check_ligand_contacts(closed, "A", os.path.join(ROOT, "data", "3QN1.cif"),
                                "A", "A8S")
    report["complexes"]["pyr1_closed_191"] = {"pdb": closed, "aba_contacts": lig}

    path = os.path.join(args.outdir, "assembly_report.json")
    with open(path, "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    log(f"\nwrote {os.path.relpath(path, ROOT)}")


if __name__ == "__main__":
    main()
