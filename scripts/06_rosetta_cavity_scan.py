#!/usr/bin/env python
"""
06_rosetta_cavity_scan.py -- ARM 1 of the pilot: does the cavity stay open?

Question
--------
The rigid-backbone truncation scan (02) is an UPPER BOUND: it deletes a side
chain and measures the hole, but never lets the neighbours move. The real
question is whether the enlarged cavity survives repacking, or whether
second-shell side chains simply collapse into it. That is a repacking question,
which is why this uses Rosetta FastRelax rather than a structure predictor.
ESMFold/AF2 are near-blind to point substitutions and would return a WT-like
backbone with high pLDDT for every variant in this panel.

Protocol (per variant, nstruct replicates)
------------------------------------------
  input  : data/pyr1_A.pdb  -- PYR1 chain A, apo, CLOSED conformation from 3QN1
  mutate : MutateResidue at each panel position
  relax  : FastRelax, ref2015_cst, 3 repeats, with CA coordinate constraints to
           the starting coordinates (-relax:constrain_relax_to_start_coords).
           Constraints hold the closed backbone so we measure side-chain
           infilling, not a global conformational change. Side chains are
           unconstrained and free to collapse into the cavity -- that collapse
           is the signal we want.
  score  : re-score with plain ref2015 (no constraint term) for dREU
  cavity : lib_cavity.cavity_volume seeded at the crystallographic ABA centroid

APO relax is deliberate: the ligand is removed so the cavity is free to close.
WT is relaxed through the identical protocol and is the reference for both
dREU and dVolume, so the comparison is internally controlled.

CAVEAT: dREU is the difference of independently relaxed total scores. It is a
coarse stability proxy, NOT a benchmarked ddG (no cartesian_ddg, no reference
correction). Use it to flag grossly destabilising variants, not to rank
near-neutral ones.

Usage:  python 06_rosetta_cavity_scan.py --index <0..40> [--nstruct 3]
        (run under the tier1_analysis conda env, which has PyRosetta)
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from variants import PANEL                      # noqa: E402
from lib_cavity import cavity_volume, atoms_from_pose_pdb  # noqa: E402

AA3 = {"A": "ALA", "G": "GLY", "S": "SER", "V": "VAL", "L": "LEU",
       "Y": "TYR", "D": "ASP", "W": "TRP", "F": "PHE", "E": "GLU",
       "K": "LYS", "R": "ARG", "N": "ASN", "Q": "GLN", "I": "ILE",
       "M": "MET", "T": "THR", "H": "HIS", "C": "CYS", "P": "PRO"}

ap = argparse.ArgumentParser()
ap.add_argument("--index", type=int, required=True)
ap.add_argument("--nstruct", type=int, default=3)
ap.add_argument("--outdir", default=os.path.join(ROOT, "results", "cavity_scan"))
args = ap.parse_args()

name, muts = PANEL[args.index]
os.makedirs(args.outdir, exist_ok=True)
os.makedirs(os.path.join(args.outdir, "pdb"), exist_ok=True)

import pyrosetta                                 # noqa: E402
from pyrosetta.rosetta.protocols.relax import FastRelax          # noqa: E402
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue  # noqa: E402
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory  # noqa: E402

pyrosetta.init(
    "-ex1 -ex2 -use_input_sc -no_optH false -flip_HNQ "
    "-relax:constrain_relax_to_start_coords "
    "-relax:coord_constrain_sidechains false "
    "-relax:ramp_constraints false "
    "-ignore_zero_occupancy false -mute all"
)

# crystallographic ABA centroid -- the cavity seed, in the input frame
aba = atoms_from_pose_pdb(os.path.join(ROOT, "data", "aba_xtal.pdb"),
                          exclude_resnames=())
REF = aba[0].mean(0)

sfxn_cst = ScoreFunctionFactory.create_score_function("ref2015_cst")
sfxn = ScoreFunctionFactory.create_score_function("ref2015")

start = pyrosetta.pose_from_pdb(os.path.join(ROOT, "data", "pyr1_A.pdb"))
pi = start.pdb_info()

rows = []
for rep in range(args.nstruct):
    pose = start.clone()
    for resnum, aa in muts:
        idx = pi.pdb2pose("A", resnum)
        assert idx > 0, f"PDB residue A{resnum} not found in pose"
        MutateResidue(idx, AA3[aa]).apply(pose)

    fr = FastRelax(sfxn_cst, 3)
    fr.apply(pose)

    tag = f"{name}_r{rep}"
    pdb_out = os.path.join(args.outdir, "pdb", f"{tag}.pdb")
    pose.dump_pdb(pdb_out)

    xyz, el = atoms_from_pose_pdb(pdb_out, chain="A")
    vol, npts = cavity_volume(xyz, el, REF)

    # CA drift vs input, to confirm the closed backbone was held.
    # FastRelax with coordinate constraints appends a virtual root residue
    # (VRT) that has no CA, so restrict to protein residues.
    def ca_coords(p):
        return np.array([list(p.residue(i).xyz("CA"))
                         for i in range(1, p.size() + 1)
                         if p.residue(i).is_protein()])
    ca_in, ca_out = ca_coords(start), ca_coords(pose)
    ca_rmsd = float(np.sqrt(((ca_in - ca_out) ** 2).sum(1).mean()))

    rows.append(dict(variant=name, rep=rep, score_ref2015=float(sfxn(pose)),
                     cavity_A3=float(vol), cavity_gridpts=npts,
                     ca_rmsd_to_input=ca_rmsd, pdb=pdb_out))
    print(f"{tag}: score={rows[-1]['score_ref2015']:.2f} "
          f"cavity={vol:.1f} A^3  CA-RMSD={ca_rmsd:.3f}", flush=True)

out = os.path.join(args.outdir, f"{name}.json")
json.dump(dict(variant=name, mutations=muts, reps=rows), open(out, "w"), indent=2)
print("wrote", out)
