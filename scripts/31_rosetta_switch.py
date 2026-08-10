#!/usr/bin/env python
"""
31_rosetta_switch.py -- ARM 2b: does the variant still PREFER the closed state?

The gap this fills
------------------
Arms 1, 2 and 2-v3 all start from the closed complex and ask whether it stays
put. None of them tests the SWITCH:

    apo/open  --ligand-->  closed  --> HAB1 binds  --> readout

Enlarging a pocket generally destabilises the CLOSED state relative to the open
one, because there is less to pack against. So the failure mode we were least
equipped to see is precisely the one that pocket expansion most plausibly causes.
Worse, 5-FOA counterselection in the Y2H removes constitutive (always-on)
receptors cheaply but does NOTHING about receptors that never close -- those are
silent, indistinguishable from "no ligand in this pool binds". A conformationally
dead library shrinks invisibly.

What is computed
----------------
Both conformers are APO, so this measures the protein's INTRINSIC conformational
preference:

    closed : data/pyr1_closed_A.pdb   (3QN1 chain A)
    open   : data/pyr1_open_A.pdb     (3K3K chain A, the ligand-free protomer)

Both were built by 30_prepare_open_closed.py with an IDENTICAL 178-residue set --
Rosetta's score is extensive, so one extra residue in one conformer would swamp
the effect being measured. The open pose is superposed into the closed frame
(rigid body; changes no energy) so both share a coordinate frame.

Per variant and replicate:
    E_closed, E_open  : ref2015 after FastRelax of each conformer
    dE_close          : E_closed - E_open   (same sequence, so composition cancels)

and after merging, versus WT:
    ddE_switch = dE_close(variant) - dE_close(WT)

    ddE >> 0  closing destabilised    -> DEAD receptor risk (silent in the screen)
    ddE ~  0  switch preserved        -> candidate
    ddE << 0  closed state favoured   -> CONSTITUTIVE risk (removable by 5-FOA,
                                         but it wastes library capacity)

Two-sided by construction, which matches how the assay actually fails.

CONSTRAINTS -- deliberately the OPPOSITE choice from Arm 2 v3. Here every CA is
constrained to its own starting conformer (-relax:constrain_relax_to_start_coords).
v3 released the gate/latch because it was asking whether they move. This script
is asking for the ENERGY OF EACH CONFORMER AS GIVEN, so each must be held in its
own state; letting them interconvert would destroy the comparison. `gate_bb_rmsd`
and `latch_bb_rmsd` are reported per state as a check that this held.

CAVEATS, which matter for how far this can be pushed:
  * This is a fixed-backbone energy difference, NOT a free energy. It ignores
    conformational entropy and says nothing about the barrier between states.
  * Both states are apo. The real cycle is apo-open -> holo-closed; modelling
    that needs a ligand per variant, which does not exist for novel ligands.
  * Rosetta comparisons across different backbones are less reliable than across
    sequences on one backbone. The WT subtraction removes the systematic part;
    what remains should be read as rank order, not as kcal/mol.

Usage:  python 31_rosetta_switch.py --index <0..40> [--rep 0] [--nstruct 1]
        (tier1_analysis conda env; must run under SLURM)
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from variants import PANEL   # noqa: E402

AA3 = {"A": "ALA", "G": "GLY", "S": "SER", "V": "VAL", "L": "LEU",
       "Y": "TYR", "D": "ASP", "W": "TRP", "F": "PHE", "E": "GLU",
       "K": "LYS", "R": "ARG"}
GATE, LATCH = list(range(85, 90)), list(range(115, 118))
BB = ("N", "CA", "C", "O")

ap = argparse.ArgumentParser()
ap.add_argument("--index", type=int, required=True)
ap.add_argument("--rep", type=int, default=0)
ap.add_argument("--nstruct", type=int, default=1)
ap.add_argument("--outdir", default=os.path.join(ROOT, "results", "switch"))
args = ap.parse_args()

name, muts = PANEL[args.index]
os.makedirs(os.path.join(args.outdir, "pdb"), exist_ok=True)

import pyrosetta                                                    # noqa: E402
from pyrosetta.rosetta.protocols.relax import FastRelax             # noqa: E402
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue  # noqa: E402
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory     # noqa: E402

pyrosetta.init(
    "-ex1 -ex2 -use_input_sc -no_optH false -flip_HNQ "
    "-relax:constrain_relax_to_start_coords "
    "-relax:coord_constrain_sidechains false "
    "-relax:ramp_constraints false "
    "-ignore_zero_occupancy false -mute all"
)

sfxn_cst = ScoreFunctionFactory.create_score_function("ref2015_cst")
sfxn = ScoreFunctionFactory.create_score_function("ref2015")

STATES = {"closed": os.path.join(ROOT, "data", "pyr1_closed_A.pdb"),
          "open": os.path.join(ROOT, "data", "pyr1_open_A.pdb")}
start = {k: pyrosetta.pose_from_pdb(v) for k, v in STATES.items()}

n_res = {k: p.size() for k, p in start.items()}
assert n_res["closed"] == n_res["open"], (
    f"pose sizes differ ({n_res}) -- energies would not be comparable; "
    "re-run 30_prepare_open_closed.py")
seq = {k: p.sequence() for k, p in start.items()}
assert seq["closed"] == seq["open"], "sequences differ between conformers"
print(f"both conformers: {n_res['closed']} residues, identical sequence",
      flush=True)


def sel_bb(pose, resnums):
    pi = pose.pdb_info()
    out = []
    for rn in resnums:
        i = pi.pdb2pose("A", rn)
        if i <= 0:
            continue
        r = pose.residue(i)
        for a in range(1, r.natoms() + 1):
            if r.atom_name(a).strip() in BB and not r.atom_is_hydrogen(a):
                out.append((i, a))
    return out


def coords(pose, sel):
    return np.array([list(pose.residue(i).xyz(a)) for i, a in sel])


def rmsd(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


sels = {k: dict(gate=sel_bb(p, GATE), latch=sel_bb(p, LATCH))
        for k, p in start.items()}
ref0 = {k: {n: coords(start[k], s) for n, s in d.items()}
        for k, d in sels.items()}
# cross-reference: how far is each state's gate from the OTHER state's gate
cross = {n: rmsd(ref0["open"][n], ref0["closed"][n]) for n in ("gate", "latch")}
print(f"starting stroke  gate {cross['gate']:.2f} A, latch {cross['latch']:.2f} A",
      flush=True)

rows = []
for k in range(args.nstruct):
    rep = args.rep + k
    rec = dict(variant=name, rep=rep)
    for state in ("closed", "open"):
        pose = start[state].clone()
        pi = pose.pdb_info()
        for resnum, aa in muts:
            idx = pi.pdb2pose("A", resnum)
            assert idx > 0, f"PDB residue A{resnum} absent from the {state} pose"
            MutateResidue(idx, AA3[aa]).apply(pose)

        FastRelax(sfxn_cst, 3).apply(pose)

        tag = f"{name}_{state}_r{rep}"
        out_pdb = os.path.join(args.outdir, "pdb", f"{tag}.pdb")
        pose.dump_pdb(out_pdb)

        rec[f"score_{state}"] = float(sfxn(pose))
        rec[f"gate_bb_rmsd_{state}"] = rmsd(ref0[state]["gate"],
                                            coords(pose, sels[state]["gate"]))
        rec[f"latch_bb_rmsd_{state}"] = rmsd(ref0[state]["latch"],
                                             coords(pose, sels[state]["latch"]))
        # did the relaxed pose drift toward the OTHER conformer?
        other = "open" if state == "closed" else "closed"
        rec[f"gate_to_{other}_{state}"] = rmsd(ref0[other]["gate"],
                                               coords(pose, sels[state]["gate"]))
        rec[f"pdb_{state}"] = out_pdb

    rec["dE_close"] = rec["score_closed"] - rec["score_open"]
    rows.append(rec)
    print(f"{name}_r{rep}: E_closed={rec['score_closed']:.2f} "
          f"E_open={rec['score_open']:.2f} dE_close={rec['dE_close']:.2f} "
          f"gate_c={rec['gate_bb_rmsd_closed']:.3f} "
          f"gate_o={rec['gate_bb_rmsd_open']:.3f}", flush=True)

out = os.path.join(args.outdir, f"{name}__r{args.rep}.json")
json.dump(dict(variant=name, mutations=muts, reps=rows), open(out, "w"), indent=2)
print("wrote", out)
