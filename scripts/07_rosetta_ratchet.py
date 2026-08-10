#!/usr/bin/env python
"""
07_rosetta_ratchet.py -- ARM 2 of the pilot: does the switch survive?

Question
--------
Arm 1 (script 06) asks only whether the cavity opens. But R79, E94, K59 and
E141 are INVARIANT across close PYR/PYL homologs while variable across the
wider superfamily -- the signature of family-specific functional constraint.
The real risk is therefore not a collapsed pocket but a broken ratchet:
a PYR1 that binds ligand yet no longer presents the closed-state epitope HAB1
reads, or one that presents it constitutively (ligand-independent, a false
positive in the yeast readout).

This arm relaxes the full ternary complex and measures whether the PYR1-HAB1
interface is preserved.

Protocol (per variant, nstruct replicates)
------------------------------------------
  input  : data/complex_AB_ABA.pdb  -- PYR1 chain A + HAB1 chain B + ABA chain X
           ABA is the deprotonated carboxylate in its crystal pose (see 00).
  mutate : MutateResidue at panel positions on chain A only
  relax  : FastRelax, ref2015_cst, 3 repeats, CA constrained to start coords
  score  : InterfaceAnalyzerMover across the A_B jump ->
             dG_separated   : interface binding energy (the ratchet readout)
             dSASA          : buried interface area
             sc_value       : shape complementarity
           plus gate/latch heavy-atom RMSD to the input, and the
           HAB1 W385 -> ABA distance (the wedge that locks the closed state).

Interpretation
--------------
  cavity opens (06) AND dG_separated ~ WT  -> candidate: bigger pocket, intact switch
  cavity opens AND dG_separated much worse -> ratchet broken, pocket gain is useless
  cavity opens AND dG_separated much better-> possible CONSTITUTIVE binder;
                                              check the apo arm before trusting it

CAVEAT: dG_separated is a Rosetta interface energy on a relaxed model, not an
affinity. Treat it as a rank-order signal for triage, not a predicted Kd.

Usage:  python 07_rosetta_ratchet.py --index <0..40> [--nstruct 3]
        (tier1_analysis conda env)
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

# 3QN1 chain A numbering
GATE = list(range(85, 90))     # S85-A89, beta3-beta4 gate loop
LATCH = list(range(115, 118))  # H115-L117, beta5-beta6 latch

ap = argparse.ArgumentParser()
ap.add_argument("--index", type=int, required=True)
ap.add_argument("--nstruct", type=int, default=3)
ap.add_argument("--outdir", default=os.path.join(ROOT, "results", "ratchet"))
args = ap.parse_args()

name, muts = PANEL[args.index]
os.makedirs(os.path.join(args.outdir, "pdb"), exist_ok=True)

import pyrosetta                                  # noqa: E402
from pyrosetta.rosetta.protocols.relax import FastRelax           # noqa: E402
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue   # noqa: E402
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover  # noqa: E402
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory   # noqa: E402

PARAMS = os.path.join(ROOT, "data", "aba_params", "A8S.params")
pyrosetta.init(
    f"-extra_res_fa {PARAMS} "
    "-ex1 -ex2 -use_input_sc -no_optH false -flip_HNQ "
    "-relax:constrain_relax_to_start_coords "
    "-relax:coord_constrain_sidechains false "
    "-relax:ramp_constraints false "
    "-ignore_zero_occupancy false -mute all"
)

sfxn_cst = ScoreFunctionFactory.create_score_function("ref2015_cst")
sfxn = ScoreFunctionFactory.create_score_function("ref2015")

start = pyrosetta.pose_from_pdb(os.path.join(ROOT, "data", "complex_AB_ABA.pdb"))
pi = start.pdb_info()


def sel_heavy(pose, chain, resnums):
    out = []
    for rn in resnums:
        i = pi.pdb2pose(chain, rn)
        if i <= 0:
            continue
        r = pose.residue(i)
        for a in range(1, r.natoms() + 1):
            if not r.atom_is_hydrogen(a):
                out.append((i, a))
    return out


def coords(pose, sel):
    return np.array([list(pose.residue(i).xyz(a)) for i, a in sel])


gate_sel = sel_heavy(start, "A", GATE)
latch_sel = sel_heavy(start, "A", LATCH)
gate0, latch0 = coords(start, gate_sel), coords(start, latch_sel)

# HAB1 W385 -- the only HAB1 residue within 5 A of ABA (4.67 A in 3QN1)
w385 = pi.pdb2pose("B", 385)
lig_idx = [i for i in range(1, start.size() + 1)
           if start.residue(i).name3() == "A8S"][0]


def w385_to_aba(pose):
    r, l = pose.residue(w385), pose.residue(lig_idx)
    return float(min(
        np.linalg.norm(np.array(list(r.xyz(a))) - np.array(list(l.xyz(b))))
        for a in range(1, r.natoms() + 1) if not r.atom_is_hydrogen(a)
        for b in range(1, l.natoms() + 1) if not l.atom_is_hydrogen(b)))


rows = []
for rep in range(args.nstruct):
    pose = start.clone()
    for resnum, aa in muts:
        idx = pi.pdb2pose("A", resnum)
        assert idx > 0, f"PDB residue A{resnum} not found"
        MutateResidue(idx, AA3[aa]).apply(pose)

    FastRelax(sfxn_cst, 3).apply(pose)

    tag = f"{name}_r{rep}"
    pdb_out = os.path.join(args.outdir, "pdb", f"{tag}.pdb")
    pose.dump_pdb(pdb_out)

    ia = InterfaceAnalyzerMover("A_B")
    ia.set_pack_separated(True)
    ia.set_scorefunction(sfxn)
    ia.apply(pose)
    d = ia.get_all_data()

    rows.append(dict(
        variant=name, rep=rep,
        total_score=float(sfxn(pose)),
        dG_separated=float(d.dG[1] + d.dG[2]) if hasattr(d, "dG") else float(ia.get_interface_dG()),
        dSASA=float(ia.get_interface_delta_sasa()),
        sc_value=float(d.sc_value),
        gate_rmsd=float(np.sqrt(((gate0 - coords(pose, gate_sel)) ** 2).sum(1).mean())),
        latch_rmsd=float(np.sqrt(((latch0 - coords(pose, latch_sel)) ** 2).sum(1).mean())),
        w385_aba_min_dist=w385_to_aba(pose),
        pdb=pdb_out))
    r = rows[-1]
    print(f"{tag}: dG_sep={r['dG_separated']:.2f} dSASA={r['dSASA']:.0f} "
          f"sc={r['sc_value']:.3f} gate={r['gate_rmsd']:.3f} "
          f"latch={r['latch_rmsd']:.3f} W385-ABA={r['w385_aba_min_dist']:.2f}",
          flush=True)

out = os.path.join(args.outdir, f"{name}.json")
json.dump(dict(variant=name, mutations=muts, reps=rows), open(out, "w"), indent=2)
print("wrote", out)
