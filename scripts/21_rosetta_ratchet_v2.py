#!/usr/bin/env python
"""
21_rosetta_ratchet_v2.py -- ARM 2, corrected. Supersedes 07_rosetta_ratchet.py.

Why this script exists
----------------------
The first Arm-2 run (script 07, nstruct=3) produced a table in which no variant
differed from WT on any readout metric. The aggregator's variance decomposition
(20_aggregate_pilot.py, README section 5b) showed that this was not a result:
between-variant spread was at or BELOW within-variant replicate noise for every
metric (dG_separated ratio 0.73, gate_rmsd 0.65, latch_rmsd 0.92). The assay was
blind. Two causes, both fixed here.

FIX 1 -- sampling. nstruct=3 against a replicate sd of ~2.8 REU on dG_separated
cannot resolve differences smaller than the noise. Raised to 20 (chunked across
array tasks; see 22_submit_ratchet_v2.sh -- one replicate is ~50 min on the
474-residue ternary complex, so 20 serial replicates per variant would be ~17 h).

FIX 2 -- the constraints were measuring themselves. Script 07 ran FastRelax with
the global flag -relax:constrain_relax_to_start_coords, which restrains EVERY
CA to the input complex, including the gate and the latch. Gate and latch CA
positions were therefore pinned by construction, which is why gate_rmsd came out
at 0.18 +- 0.004 A for all 41 variants: that number measured the restraint, not
the protein. Asking "does the gate move?" while holding the gate still is not a
test.

  Here the global flag is OFF. Coordinate constraints are instead applied
  explicitly, CA-only, to every protein residue EXCEPT the gate (85-89) and the
  latch (115-117) of chain A. The rest of the scaffold and all of HAB1 stay
  pinned, so variants remain comparable and the pose cannot drift wholesale --
  but the gate and latch backbone is genuinely free to move, respond to the
  mutation, and be measured.

  This is the controlled version of the experiment: everything held fixed except
  the thing being asked about.

New metrics that the constraints do not bound
---------------------------------------------
  gate_bb_rmsd / latch_bb_rmsd   backbone (N, CA, C, O) RMSD to the input. Now
                                 free to vary; these are the primary readout.
  gate_latch_min_dist            closest approach of any gate heavy atom to any
                                 latch heavy atom. The gate closing ONTO the
                                 latch is the physical event the ratchet is; a
                                 variant that opens this has lost the switch.
  gate_latch_com_dist            gate centroid to latch centroid, a coarser,
                                 less noise-prone version of the same thing.
  lig_gate_min_dist              ABA to gate -- does the ligand still contact the
                                 loop it is supposed to hold shut?

Retained for continuity with the script-07 run: dG_separated, dSASA, sc_value,
gate_rmsd, latch_rmsd (all-heavy-atom), w385_aba_min_dist, total_score.

Chunking
--------
--rep-offset shifts the replicate numbering so several array tasks can fill
different replicates of the same variant without colliding. Each task writes
results/ratchet_v2/<variant>__c<chunk>.json; 23_merge_ratchet_v2.py concatenates
them into <variant>.json.

Usage:  python 21_rosetta_ratchet_v2.py --index <0..40> [--nstruct 4]
                                        [--rep-offset 0] [--chunk 0]
        (tier1_analysis conda env; must run under SLURM -- FastRelax exceeds the
         1 GB login-node cgroup cap)
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

GATE = list(range(85, 90))     # S85-A89, beta3-beta4 gate loop
LATCH = list(range(115, 118))  # H115-L117, beta5-beta6 latch
BB = ("N", "CA", "C", "O")

ap = argparse.ArgumentParser()
ap.add_argument("--index", type=int, required=True)
ap.add_argument("--nstruct", type=int, default=4)
ap.add_argument("--rep-offset", type=int, default=0)
ap.add_argument("--chunk", type=int, default=0)
ap.add_argument("--outdir", default=os.path.join(ROOT, "results", "ratchet_v2"))
args = ap.parse_args()

name, muts = PANEL[args.index]
os.makedirs(os.path.join(args.outdir, "pdb"), exist_ok=True)

import pyrosetta                                                        # noqa: E402
from pyrosetta.rosetta.protocols.relax import FastRelax                 # noqa: E402
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue      # noqa: E402
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover  # noqa: E402
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory         # noqa: E402
from pyrosetta.rosetta.protocols.constraint_generator import (           # noqa: E402
    CoordinateConstraintGenerator, AddConstraints)
from pyrosetta.rosetta.core.select.residue_selector import (             # noqa: E402
    ResidueIndexSelector, NotResidueSelector)

PARAMS = os.path.join(ROOT, "data", "aba_params", "A8S.params")
# NOTE: -relax:constrain_relax_to_start_coords is deliberately ABSENT (see FIX 2).
pyrosetta.init(
    f"-extra_res_fa {PARAMS} "
    "-ex1 -ex2 -use_input_sc -no_optH false -flip_HNQ "
    "-ignore_zero_occupancy false -mute all"
)

sfxn_cst = ScoreFunctionFactory.create_score_function("ref2015_cst")  # coord_cst = 1.0
sfxn = ScoreFunctionFactory.create_score_function("ref2015")

start = pyrosetta.pose_from_pdb(os.path.join(ROOT, "data", "complex_AB_ABA.pdb"))
pi = start.pdb_info()


def pose_idx(chain, resnums):
    out = []
    for rn in resnums:
        i = pi.pdb2pose(chain, rn)
        if i > 0:
            out.append(i)
    return out


gate_idx = pose_idx("A", GATE)
latch_idx = pose_idx("A", LATCH)
assert len(gate_idx) == len(GATE) and len(latch_idx) == len(LATCH), \
    "gate/latch residues missing from the complex"

# ---- selective coordinate constraints: everything EXCEPT gate and latch ----
free_sel = ResidueIndexSelector(",".join(str(i) for i in gate_idx + latch_idx))
constrain_sel = NotResidueSelector(free_sel)

cg = CoordinateConstraintGenerator()
cg.set_residue_selector(constrain_sel)
cg.set_ca_only(True)
cg.set_sidechain(False)
cg.set_sd(0.5)          # matches the stiffness of the old global flag
cg.set_bounded(False)
add_cst = AddConstraints()
add_cst.add_generator(cg)

n_free = len(gate_idx) + len(latch_idx)
print(f"coordinate constraints: CA-only on all residues EXCEPT {n_free} "
      f"gate/latch residues (pose idx {gate_idx + latch_idx}) -- these are FREE",
      flush=True)


def sel_atoms(pose, idxs, names=None):
    """(residue, atom) pairs: all heavy atoms, or just the named ones."""
    out = []
    for i in idxs:
        r = pose.residue(i)
        for a in range(1, r.natoms() + 1):
            if r.atom_is_hydrogen(a):
                continue
            if names is None or r.atom_name(a).strip() in names:
                out.append((i, a))
    return out


def coords(pose, sel):
    return np.array([list(pose.residue(i).xyz(a)) for i, a in sel])


gate_hv = sel_atoms(start, gate_idx)
latch_hv = sel_atoms(start, latch_idx)
gate_bb = sel_atoms(start, gate_idx, BB)
latch_bb = sel_atoms(start, latch_idx, BB)
gate0, latch0 = coords(start, gate_hv), coords(start, latch_hv)
gate_bb0, latch_bb0 = coords(start, gate_bb), coords(start, latch_bb)

w385 = pi.pdb2pose("B", 385)
lig_idx = [i for i in range(1, start.size() + 1)
           if start.residue(i).name3() == "A8S"][0]


def heavy_xyz(pose, i):
    r = pose.residue(i)
    return np.array([list(r.xyz(a)) for a in range(1, r.natoms() + 1)
                     if not r.atom_is_hydrogen(a)])


def min_dist(A, B):
    return float(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=-1).min())


def rmsd(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


rows = []
for k in range(args.nstruct):
    rep = args.rep_offset + k
    pose = start.clone()
    for resnum, aa in muts:
        idx = pi.pdb2pose("A", resnum)
        assert idx > 0, f"PDB residue A{resnum} not found"
        MutateResidue(idx, AA3[aa]).apply(pose)

    add_cst.apply(pose)          # re-applied per replicate; mutation may change idx
    FastRelax(sfxn_cst, 3).apply(pose)

    tag = f"{name}_r{rep}"
    pdb_out = os.path.join(args.outdir, "pdb", f"{tag}.pdb")
    pose.dump_pdb(pdb_out)

    ia = InterfaceAnalyzerMover("A_B")
    ia.set_pack_separated(True)
    ia.set_scorefunction(sfxn)
    ia.apply(pose)
    d = ia.get_all_data()

    g_hv, l_hv = coords(pose, gate_hv), coords(pose, latch_hv)
    gate_xyz = np.vstack([heavy_xyz(pose, i) for i in gate_idx])
    latch_xyz = np.vstack([heavy_xyz(pose, i) for i in latch_idx])
    lig_xyz = heavy_xyz(pose, lig_idx)

    rows.append(dict(
        variant=name, rep=rep,
        total_score=float(sfxn(pose)),
        dG_separated=float(d.dG[1] + d.dG[2]) if hasattr(d, "dG") else float(ia.get_interface_dG()),
        dSASA=float(ia.get_interface_delta_sasa()),
        sc_value=float(d.sc_value),
        # retained (bounded in v1, free here)
        gate_rmsd=rmsd(gate0, g_hv),
        latch_rmsd=rmsd(latch0, l_hv),
        # new, unbounded -- the primary readout
        gate_bb_rmsd=rmsd(gate_bb0, coords(pose, gate_bb)),
        latch_bb_rmsd=rmsd(latch_bb0, coords(pose, latch_bb)),
        gate_latch_min_dist=min_dist(gate_xyz, latch_xyz),
        gate_latch_com_dist=float(np.linalg.norm(
            gate_xyz.mean(0) - latch_xyz.mean(0))),
        lig_gate_min_dist=min_dist(lig_xyz, gate_xyz),
        w385_aba_min_dist=min_dist(heavy_xyz(pose, w385), lig_xyz),
        pdb=pdb_out))
    r = rows[-1]
    print(f"{tag}: dG_sep={r['dG_separated']:.2f} sc={r['sc_value']:.3f} "
          f"gate_bb={r['gate_bb_rmsd']:.3f} latch_bb={r['latch_bb_rmsd']:.3f} "
          f"g-l_min={r['gate_latch_min_dist']:.2f} "
          f"W385-ABA={r['w385_aba_min_dist']:.2f}", flush=True)

out = os.path.join(args.outdir, f"{name}__c{args.chunk}.json")
json.dump(dict(variant=name, mutations=muts, chunk=args.chunk, reps=rows),
          open(out, "w"), indent=2)
print("wrote", out)
