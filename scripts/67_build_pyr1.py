#!/usr/bin/env python
"""
67_build_pyr1.py -- build the COMPLETE 191-residue PYR1 monomer, open and closed, from
the crystals, with no chain breaks and no silent choices.

Run with the tier1_analysis env (the only one carrying PyRosetta):
    /bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python scripts/67_build_pyr1.py

WHY
---
Every MD system built before this script inherited whatever the crystal happened to
model. That gave the systems artificial buried termini (tleap makes each fragment its
own chain, so a gapped build carries extra charged NH3+/COO- pairs inside the fold),
a C-terminus at 181 where the real protein runs to 191, and -- worse -- it gave the
open and closed forms DIFFERENT defects, so the open/closed contrast was partly a
contrast between two different construct definitions. See README 28.

WHAT THE CRYSTALS ACTUALLY CONTAIN (audited 2026-08-17, not assumed)
--------------------------------------------------------------------
    3K3K chain A (open)    residues 1-183, ZERO internal gaps, ZERO broken bonds
    3K3K chain B           residues 2-184, ZERO internal gaps
    3QN1 chain A (closed)  residues 1-181, ONE gap (69-70), C(68)-N(71) = 4.64 A
    3QN1 chain B (HAB1)    residues 185-505, 26 missing over 3 gaps  -> script 67b

An earlier note in this project claimed 3QN1 also broke at C(1)-N(3) and that 3K3K
shared its gaps. Both were wrong: that 3.02 A break belongs to our own derived
data/md/S1_apo_open/protein.pdb, which script 30 produced by intersecting the two
crystals. 3QN1 models residue 2 -- as ALA, which is the P2A substitution, a mutation
rather than a gap.

THE FOUR BUILD DECISIONS, EACH SETTLED BY MEASUREMENT
------------------------------------------------------
1. ALTERNATE CONFORMATIONS. 3K3K carries 0.5-occupancy alternates at >=8 positions
   including R116, a LATCH residue; 3QN1 carries none. Taking each parser's default
   would mean the open and closed systems differ by an unrecorded coin-flip at exactly
   the positions this project measures. Altloc A is selected explicitly and counted.

2. THE 69-70 LOOP IS GRAFTED, NOT INVENTED. 3K3K models 69/70 fully. Superposing the
   flanking backbone (66-68, 71-73) of 3K3K onto 3QN1 gives an RMSD of 0.85 A -- and
   it is still 0.85 A using only residues 68 and 71, so the graft is determined by the
   flanks rather than by the modeller. The loop also sits 19.1 A from the gate and
   20.3 A from ABA, well outside the switching machinery. Real coordinates from the
   same protein beat a de novo loop.

3. BOTH TAILS ARE REBUILT FROM 181. 3K3K models two extra C-terminal residues, but
   their B-factors are 142 and 151 against a core mean of 49 -- 2.9x and 3.1x, i.e.
   present but barely ordered. Keeping them would give the open form a 2-residue head
   start built from real-but-meaningless density while the closed form got 10 modelled
   residues, reintroducing exactly the open/closed protocol asymmetry this script
   exists to remove. So 3K3K is truncated at 181 and BOTH forms get an identically
   generated 182-191.

4. P2 IS RESTORED. 3QN1's A2 is not wild-type. phi(A2) = -67.4 deg, inside proline's
   -63 +/- 15 window, so this is a ring closure and not a backbone rebuild.

THE TAIL IS A SAMPLE, NOT A PREDICTION
--------------------------------------
182-191 reads SGDGSGSQVT -- three glycines and three serines, a natural GS-linker.
It is disordered in every structure ever solved of this protein and no starting
conformation for it is meaningful. This script therefore SAMPLES `--n-tail` seeded
conformers and keeps the best-scoring one, rather than dictating an extended chain
that would also inflate the solvent box. The seed is recorded in build_report.json.
Downstream analysis must exclude 182-191 from core-fit and RMSF masks.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
from Bio.PDB import PDBIO, Superimposer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_structqc  # noqa: E402
from lib_rosetta import assert_frozen  # noqa: E402
from lib_resnumber import (  # noqa: E402
    PYR1_SEQ, PYR1_LANDMARKS, THREE2ONE, load_model, select_altloc,
    chain_residues, chain_sequence, peptide_breaks, assert_identity,
    verify_build, write_residue_map,
)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BB = ("N", "CA", "C", "O")

# Flanks used for the 69-70 graft. Two independent sets: the wide one is used, the
# narrow one is a control -- if they disagree the graft is not flank-determined.
GRAFT_FLANK = [66, 67, 68, 71, 72, 73]
GRAFT_FLANK_NARROW = [68, 71]
GRAFT_RMSD_MAX = 1.5      # A; measured 0.85, so this trips only on a real problem

TAIL_START = 182
TAIL_END = 191
KEEP_TO = 181             # both crystals truncated here before the tail is built


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------------------------
# stage A -- crystals to a gap-free 1..181 backbone, in Biopython
# ---------------------------------------------------------------------------------
def graft_loop(donor_chain, acceptor_resmap, residues, flank, label):
    """Copy `residues` from `donor_chain` into `acceptor_resmap` after superposing flanks.

    Returns (rmsd_wide, rmsd_narrow). Raises if the graft is not flank-determined.
    """
    donor = chain_residues(donor_chain)
    for i in residues + flank:
        if i not in donor:
            raise AssertionError(f"{label}: donor lacks residue {i}; cannot graft")
    for i in flank:
        if i not in acceptor_resmap:
            raise AssertionError(f"{label}: acceptor lacks flank residue {i}")
        d, a = donor[i].get_resname(), acceptor_resmap[i].get_resname()
        if THREE2ONE.get(d.upper()) != THREE2ONE.get(a.upper()):
            raise AssertionError(
                f"{label}: flank residue {i} is {a} in acceptor but {d} in donor -- "
                "the two structures are not in register")

    def fit(fl):
        sup = Superimposer()
        sup.set_atoms([acceptor_resmap[i][a] for i in fl for a in BB],
                      [donor[i][a] for i in fl for a in BB])
        return sup

    sup_wide, sup_narrow = fit(flank), fit(GRAFT_FLANK_NARROW)
    log(f"  graft {label}: flank RMSD {sup_wide.rms:.2f} A (wide, n={len(flank)}), "
        f"{sup_narrow.rms:.2f} A (narrow, n={len(GRAFT_FLANK_NARROW)})")
    if sup_wide.rms > GRAFT_RMSD_MAX:
        raise AssertionError(f"{label}: graft flank RMSD {sup_wide.rms:.2f} A exceeds "
                             f"{GRAFT_RMSD_MAX} A -- the loop is NOT flank-determined, "
                             "so it must be modelled de novo, not grafted")
    if abs(sup_wide.rms - sup_narrow.rms) > 0.5:
        raise AssertionError(f"{label}: wide and narrow flank fits disagree "
                             f"({sup_wide.rms:.2f} vs {sup_narrow.rms:.2f} A)")

    # apply the WIDE transform to the whole donor chain, then lift out the loop
    sup_wide.apply(list(donor_chain.get_atoms()))
    for i in residues:
        acceptor_resmap[i] = donor[i]
    return float(sup_wide.rms), float(sup_narrow.rms)


def write_resmap(resmap, out_path, chain_id="A"):
    """Write an ordered {num: Residue} map as a PDB, KEEPING native numbering.

    Deliberately not renumbered to 1..N. lib_structqc.core_preservation compares the
    stage A and stage B structures residue number by residue number, so renumbering
    here would silently misalign that comparison for any chain that does not start at
    residue 1 -- e.g. 3K3K chain B, which starts at 2. The comparison would then
    "pass" while comparing each residue against its neighbour.
    """
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain

    st = Structure("built")
    mo = Model(0)
    ch = Chain(chain_id)
    st.add(mo)
    mo.add(ch)
    for num in sorted(resmap):
        res = resmap[num].copy()
        res.id = (" ", num, " ")
        res.detach_parent()
        ch.add(res)
    io = PDBIO()
    io.set_structure(st)
    io.save(out_path)
    return out_path


def ligand_shell(model, protein_chain, cutoff=4.5):
    """Residues of `protein_chain` within `cutoff` of any non-water heteroatom.

    These must be FROZEN during relax, and the reason is the bug that motivated this
    function. Script 67 builds the protein alone -- the ligand is added later by tleap
    from A8S.mol2. So when FastRelax repacks the shell around a rebuilt region, any
    pocket residue caught in that shell is being repacked into an EMPTY cavity. K59
    sits within 6 A of the grafted 68-72 loop, and it duly collapsed inward: its
    closest approach to ABA's crystal position went from 2.85 A (a salt bridge) to
    1.68 A (a hard clash), while the backbone check stayed at 0.000 A because it only
    watches N/CA/C.

    Their crystal positions were determined WITH the ligand present, which makes those
    coordinates better than anything a ligand-free repack can produce. Freezing them is
    not a compromise; it is the correct answer.
    """
    het = [r for ch in model for r in ch
           if r.id[0].strip() and r.get_resname() not in ("HOH", "WAT")]
    if not het:
        return set()
    lc = np.array([a.coord for r in het for a in r if a.element != "H"])
    chain = model[protein_chain]
    select_altloc(chain)
    out = set()
    for num, res in chain_residues(chain).items():
        pc = np.array([a.coord for a in res if a.element != "H"])
        if np.min(np.linalg.norm(pc[:, None, :] - lc[None, :, :], axis=-1)) <= cutoff:
            out.add(num)
    return out


def stage_a(spec, outdir):
    """Crystal -> gap-free, altloc-resolved, 1..181 PDB. Returns a provenance dict."""
    log(f"\n== stage A: {spec['name']} from {spec['cif']} chain {spec['chain']} ==")
    model = load_model(os.path.join(ROOT, spec["cif"]))
    protected = ligand_shell(model, spec["chain"])
    if protected:
        log(f"  {len(protected)} residue(s) line a ligand and will be FROZEN, not "
            f"repacked: {sorted(protected)}")
    chain = model[spec["chain"]]
    n_alt = select_altloc(chain, keep="A")
    log(f"  altloc A selected at {n_alt} atom position(s)")

    resmap = chain_residues(chain)
    nums = sorted(resmap)
    log(f"  modelled {len(nums)} residues {nums[0]}..{nums[-1]}")

    # identity check BEFORE anything is modified, excluding residues we know differ
    pre = {k: v for k, v in PYR1_LANDMARKS.items() if k not in spec.get("known_diff", {})}
    assert_identity(resmap, pre, f"{spec['name']} (pre-build)")
    for num, want in spec.get("known_diff", {}).items():
        got = resmap[num].get_resname().upper() if num in resmap else "MISSING"
        if got != want:
            raise AssertionError(f"{spec['name']}: expected the known variant {want} at "
                                 f"residue {num}, found {got}")
        log(f"  confirmed known variant: residue {num} is {got} (will be corrected)")
    log(f"  identity assertions pass ({len(pre)} landmarks)")

    graft_info = None
    if spec.get("graft_from"):
        donor_model = load_model(os.path.join(ROOT, spec["graft_from"]["cif"]))
        donor_chain = donor_model[spec["graft_from"]["chain"]]
        select_altloc(donor_chain, keep="A")
        w, n = graft_loop(donor_chain, resmap, spec["graft_from"]["residues"],
                          GRAFT_FLANK, spec["name"])
        graft_info = {"donor": spec["graft_from"]["cif"],
                      "chain": spec["graft_from"]["chain"],
                      "residues": spec["graft_from"]["residues"],
                      "flank": GRAFT_FLANK, "rmsd_wide": w, "rmsd_narrow": n}

    # truncate to 1..181 so both forms enter stage B identically
    dropped = [n for n in sorted(resmap) if n > KEEP_TO]
    for n in dropped:
        del resmap[n]
    if dropped:
        log(f"  truncated at {KEEP_TO}, dropping {dropped} "
            "(poorly ordered; both forms rebuild 182-191 identically)")

    # An N-terminal stretch may also be unmodelled: 3K3K chain B starts at residue 2,
    # so the dimer's second protomer is short its initiator methionine. That is built
    # in stage B by the same sampling used for the C-terminal tail -- INTERNAL gaps,
    # by contrast, must never be handled this way, which is why the check below
    # distinguishes the two rather than reporting one count of "missing residues".
    nums = sorted(resmap)
    nterm_missing = list(range(1, nums[0]))
    internal_missing = sorted(set(range(nums[0], KEEP_TO + 1)) - set(nums))
    if internal_missing:
        raise AssertionError(
            f"{spec['name']}: residues {internal_missing} are INTERNAL gaps still "
            "unfilled after grafting. An internal gap cannot be built by terminal "
            "extension -- supply a graft donor or model the loop explicitly.")
    if nterm_missing:
        log(f"  N-terminal residues {nterm_missing} unmodelled; will be built in stage B")

    breaks = peptide_breaks(resmap, tol=2.0)   # graft closes to ~1 A; stage B fixes it
    log(f"  residues {nums[0]}..{KEEP_TO} complete; {len(breaks)} bond(s) >2.0 A "
        f"pre-relax: {[(a, b, d) for a, _, b, _, d in breaks]}")

    seq, _ = chain_sequence(resmap)
    want = PYR1_SEQ[nums[0] - 1:KEEP_TO]
    diffs = [(i + nums[0], e, g) for i, (e, g) in enumerate(zip(want, seq)) if e != g]
    log(f"  sequence {nums[0]}..{KEEP_TO} differs from canonical at: {diffs} "
        "(expected: the P2A to be corrected in stage B, if present)")

    out = os.path.join(outdir, f"{spec['name']}_stageA.pdb")
    write_resmap(resmap, out)
    return {"stage_a_pdb": out, "n_altloc_positions": n_alt, "graft": graft_info,
            "dropped_residues": dropped, "pre_relax_long_bonds": breaks,
            "seq_diffs_vs_canonical": diffs, "nterm_missing": nterm_missing,
            "first_modelled": nums[0], "protected": sorted(protected)}


# ---------------------------------------------------------------------------------
# stage B -- PyRosetta: restore P2, build 182-191, relax only what was rebuilt
# ---------------------------------------------------------------------------------
def _residue(rts, aa1):
    """Ideal-geometry residue for a one-letter code."""
    from pyrosetta.rosetta.core.conformation import ResidueFactory
    from lib_resnumber import ONE2THREE
    return ResidueFactory.create_residue(rts.name_map(ONE2THREE[aa1]))


def stage_b(stage_a_pdb, spec, outdir, n_tail, seed, relax_repeats, nterm_missing=(),
            protected=()):
    import pyrosetta
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from lib_rosetta import restrict_packing

    log(f"\n== stage B: {spec['name']} (seed {seed}, {n_tail} tail conformers) ==")
    pose0 = pyrosetta.pose_from_pdb(stage_a_pdb)
    rts = pose0.residue_type_set_for_pose()
    rng = np.random.default_rng(seed)

    # --- N-terminal extension FIRST, so pose numbering == native numbering ---------
    # Everything downstream (the P2 mutation, the graft's rebuilt set, the tail) is
    # written in NATIVE residue numbers. Rosetta always numbers a pose 1..N, so a
    # chain that starts at native residue 2 would silently shift every one of those
    # references by one. Prepending before anything else removes the offset entirely
    # rather than threading it through the rest of the function.
    for num in sorted(nterm_missing, reverse=True):
        pose0.prepend_polymer_residue_before_seqpos(
            _residue(rts, PYR1_SEQ[num - 1]), 1, True)
        pose0.set_phi(1, float(rng.uniform(-150, -60)))
        pose0.set_psi(1, float(rng.uniform(-60, 160)))
        pose0.set_omega(1, 180.0)
        log(f"  prepended residue {num} ({PYR1_SEQ[num - 1]}) at the N-terminus")

    assert pose0.total_residue() == KEEP_TO, \
        f"stage A + N-term extension gave {pose0.total_residue()} residues, expected {KEEP_TO}"
    for probe in (1, 2, 59, 85, 108, 115):
        assert pose0.residue(probe).name1() in (PYR1_SEQ[probe - 1],
                                                THREE2ONE.get(spec.get("known_diff", {})
                                                              .get(probe, "___"), "_")), \
            (f"pose position {probe} is {pose0.residue(probe).name1()}, expected "
             f"{PYR1_SEQ[probe - 1]} -- pose numbering is NOT native numbering")

    # --- restore wild-type P2 where the crystal carries the A2 substitution --------
    for num, want in spec.get("known_diff", {}).items():
        native = PYR1_SEQ[num - 1]
        have = pose0.residue(num).name1()
        assert have == THREE2ONE[want], \
            f"residue {num} is {have}, expected the crystal's {THREE2ONE[want]}"
        MutateResidue(num, {"P": "PRO", "A": "ALA"}[native]).apply(pose0)
        log(f"  residue {num}: {have} -> {pose0.residue(num).name1()} (wild-type restored)")

    rebuilt = set(spec.get("known_diff", {})) | set(nterm_missing)
    if nterm_missing:
        # the residue following the built N-terminus also moves, its bond is new
        rebuilt.add(max(nterm_missing) + 1)
    if spec.get("graft_from"):
        # the graft itself plus one flanking residue either side, whose bonds moved
        g = spec["graft_from"]["residues"]
        rebuilt |= set(range(min(g) - 1, max(g) + 2))

    # --- build the tail: sample seeded conformers, keep the best ------------------
    sfxn = pyrosetta.create_score_function("ref2015_cart")
    tail_seq = PYR1_SEQ[TAIL_START - 1:TAIL_END]
    log(f"  building tail {TAIL_START}-{TAIL_END} = {tail_seq}")

    best, best_score, best_idx = None, None, None
    for k in range(n_tail):
        pose = pose0.clone()
        for aa1 in tail_seq:
            pose.append_polymer_residue_after_seqpos(
                _residue(rts, aa1), pose.total_residue(), True)
            i = pose.total_residue()
            # coil-region torsions rather than a fully extended chain: an extended
            # 10-mer projects ~35 A into solvent and would inflate the periodic box
            # (and therefore the water count and the cost of every ns) for a segment
            # that is disordered anyway.
            pose.set_phi(i, float(rng.uniform(-150, -60)))
            pose.set_psi(i, float(rng.uniform(-60, 160)))
            pose.set_omega(i, 180.0)
        sc = sfxn(pose)
        if best_score is None or sc < best_score:
            best, best_score, best_idx = pose, sc, k
    log(f"  {n_tail} conformers scored; kept #{best_idx} at {best_score:.1f} REU")
    pose = best
    rebuilt |= set(range(TAIL_START - 1, TAIL_END + 1))

    # --- relax ONLY what was rebuilt, plus repack its neighbours ------------------
    # Ligand-lining residues are excluded from BOTH sets. The pose has no ligand in it
    # -- tleap adds ABA later from A8S.mol2 -- so repacking a pocket residue here means
    # optimising it into an empty cavity. See ligand_shell() for what that cost.
    protected = set(protected)
    overlap = protected & rebuilt
    if overlap:
        raise AssertionError(
            f"residues {sorted(overlap)} are both rebuilt and ligand-lining. A rebuilt "
            "pocket residue cannot be relaxed without the ligand present; add ligand "
            "context to this step rather than silently freezing or moving it.")
    sel = ResidueIndexSelector(",".join(str(i) for i in sorted(rebuilt)))
    nbr = NeighborhoodResidueSelector(sel, 6.0, True)
    move = sel.apply(pose)
    pack = nbr.apply(pose)
    mm = MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    frozen_pocket = []
    for i in range(1, pose.total_residue() + 1):
        if i in protected:
            if pack[i]:
                frozen_pocket.append(i)
            continue
        if move[i]:
            mm.set_bb(i, True)
            mm.set_chi(i, True)
        elif pack[i]:
            mm.set_chi(i, True)
    n_bb = sum(1 for i in range(1, pose.total_residue() + 1)
               if move[i] and i not in protected)
    n_chi = sum(1 for i in range(1, pose.total_residue() + 1)
                if pack[i] and i not in protected)
    log(f"  relaxing {n_bb} residues (backbone+sidechain), repacking {n_chi} neighbours; "
        "the rest of the protein is held fixed so the crystal core is preserved")
    if frozen_pocket:
        log(f"  {len(frozen_pocket)} ligand-lining residue(s) fell in the repack shell "
            f"and were frozen instead: {frozen_pocket}")

    # ⚠ THE MOVEMAP IS NOT ENOUGH. set_movemap restricts MINIMISATION only; FastRelax
    # repacks through a TaskFactory and, given none, repacks the WHOLE pose regardless.
    # That is how K59 collapsed into the empty ABA cavity while sitting in the frozen
    # set, and how the dimer's R116 crossed into the other protomer. See lib_rosetta.
    allowed_pack = [i for i in range(1, pose.total_residue() + 1)
                    if (move[i] or pack[i]) and i not in protected]
    tf, packable = restrict_packing(pose, allowed_pack)
    log(f"  packer restricted to {len(packable)} residues (design disabled); "
        f"{pose.total_residue() - len(packable)} keep their input rotamers")

    relax = FastRelax(sfxn, relax_repeats)
    relax.cartesian(True)          # cartesian: the graft junction has non-ideal bond
    relax.set_movemap(mm)          # lengths, which torsion-space relax cannot repair
    relax.set_task_factory(tf)
    relax.min_type("lbfgs_armijo_nonmonotone")
    relax.apply(pose)
    log(f"  post-relax score {sfxn(pose):.1f} REU")

    out = os.path.join(outdir, f"{spec['name']}.pdb")
    pose.dump_pdb(out)
    return {"stage_b_pdb": out, "tail_seed": seed, "tail_conformers": n_tail,
            "tail_chosen": best_idx, "tail_score_REU": float(best_score),
            "relaxed_residues": sorted(rebuilt), "n_repacked": n_chi,
            "packable_residues": sorted(packable),
            "frozen_ligand_shell": sorted(protected),
            "frozen_in_repack_shell": frozen_pocket,
            "final_score_REU": float(sfxn(pose)), "relax_repeats": relax_repeats}


# ---------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=os.path.join(ROOT, "data", "structures_191"))
    ap.add_argument("--n-tail", type=int, default=12,
                    help="seeded conformers sampled for the disordered 182-191 tail")
    ap.add_argument("--seed", type=int, default=20260817,
                    help="recorded in build_report.json; the tail is a sample, so the "
                         "seed is part of the provenance")
    ap.add_argument("--relax-repeats", type=int, default=2)
    ap.add_argument("--only", default=None, help="build just this system")
    args = ap.parse_args()

    specs = [
        {"name": "pyr1_open_191", "cif": "data/3K3K.cif", "chain": "A",
         "known_diff": {}, "graft_from": None},
        {"name": "pyr1_closed_191", "cif": "data/3QN1.cif", "chain": "A",
         "known_diff": {2: "ALA"},          # the P2A substitution, corrected in stage B
         "graft_from": {"cif": "data/3K3K.cif", "chain": "A", "residues": [69, 70]}},
        # S3's dimer needs BOTH protomers. 3K3K chain B is modelled 2..184, so it is
        # gap-free like chain A but short its initiator methionine, which stage B
        # prepends. Built as its own system rather than by copying chain A: the two
        # protomers are crystallographically independent and need not be identical.
        {"name": "pyr1_open_191_protomerB", "cif": "data/3K3K.cif", "chain": "B",
         "known_diff": {}, "graft_from": None},
    ]
    if args.only:
        specs = [s for s in specs if s["name"] == args.only]
        if not specs:
            sys.exit(f"no spec named {args.only}")

    os.makedirs(args.outdir, exist_ok=True)
    import pyrosetta
    # -constant_seed/-jran: FastRelax is a Monte Carlo protocol. Without pinning
    # Rosetta's own RNG the same input gave -195.4 and -202.4 REU on two identical
    # runs, so the structure that reached MD would not have been reproducible and the
    # numpy seed below would have been recording only half the randomness.
    pyrosetta.init(f"-mute all -ex1 -ex2aro -detect_disulf false "
                   f"-run:constant_seed -run:jran {args.seed}")

    report = {"built": datetime.now(timezone.utc).isoformat(),
              "builder": "67_build_pyr1.py", "seed": args.seed,
              "canonical_sequence_md5": "5ff74f5d34", "systems": {}}

    for spec in specs:
        a = stage_a(spec, args.outdir)
        b = stage_b(a["stage_a_pdb"], spec, args.outdir, args.n_tail,
                    args.seed, args.relax_repeats,
                    nterm_missing=a["nterm_missing"], protected=a["protected"])
        log(f"\n== stage C: verifying {spec['name']} ==")
        v = verify_build(b["stage_b_pdb"], "A", label=spec["name"])
        log(f"  PASS: {v['n_residues']} residues {v['first']}..{v['last']}, "
            f"sequence md5 {v['sequence_md5']}, 0 peptide-bond breaks, "
            f"{len(v['landmarks_checked'])} landmarks verified")
        # identity and connectivity are necessary but not sufficient: prove the
        # MoveMap actually froze the crystal core, that nothing clashes, and that the
        # rebuilt tail did not inflate the solvent box we then pay for every ns.
        qc = lib_structqc.report(
            before_pdb=a["stage_a_pdb"], after_pdb=b["stage_b_pdb"], chain_ids=["A"],
            rebuilt=set(b["relaxed_residues"]),
            exclude_for_cost=set(range(TAIL_START, TAIL_END + 1)), log=log)
        # All-atom, not backbone. The two defects this catches were side-chain-only
        # and the backbone check reported 0.000 A straight through them.
        before = chain_residues(load_model(a["stage_a_pdb"])["A"])
        after = chain_residues(load_model(b["stage_b_pdb"])["A"])
        # Untouched = neither rebuilt NOR allowed to repack. Repacked neighbours are
        # supposed to move their side chains, so including them here would make the
        # assertion fire on correct behaviour.
        frozen = set(a["protected"])
        untouched = (set(before) - set(b["relaxed_residues"])
                     - set(b["packable_residues"]))
        n_f, d_f = assert_frozen(before, after, frozen, f"{spec['name']} ligand shell")
        n_u, d_u = assert_frozen(before, after, untouched,
                                 f"{spec['name']} untouched core")
        log(f"  all-atom freeze verified: {n_f} ligand-lining (max drift {d_f:.3f} A) "
            f"and {n_u} untouched (max {d_u:.3f} A); "
            f"{len(b['packable_residues'])} residues repacked")
        qc["frozen_verified_ligand_shell"] = n_f
        qc["frozen_verified_core"] = n_u
        qc["max_drift_ligand_shell_A"] = round(d_f, 4)
        qc["max_drift_core_A"] = round(d_u, 4)
        rmap = os.path.join(args.outdir, f"{spec['name']}_residue_map.json")
        write_residue_map(b["stage_b_pdb"], rmap, ["A"])
        log(f"  wrote {os.path.relpath(rmap, ROOT)}")
        report["systems"][spec["name"]] = {**a, **b, "verification": v, "qc": qc,
                                           "residue_map": rmap}

    path = os.path.join(args.outdir, "build_report.json")
    with open(path, "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    log(f"\nwrote {os.path.relpath(path, ROOT)}")


if __name__ == "__main__":
    main()
