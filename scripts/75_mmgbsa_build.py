#!/usr/bin/env python
"""
75_mmgbsa_build.py -- build the variant structures for MM-GBSA rescoring.

WHY MM-GBSA, AND WHAT IT IS BEING ASKED
---------------------------------------
README 36 established, with two independent samplers, that K59R is unreachable for
a REASON IN THE ENERGY FUNCTION: ref2015 pays +10.1 REU of Lazaridis-Karplus
desolvation to bury K59's ammonium, and coupled moves -- which adds backbone and
ligand flexibility -- still retains K59 in 0% of trajectories even with the COGNATE
ligand. More sampling cannot fix that.

So change the solvation model. MM-GBSA replaces the pairwise LK approximation with
a generalised-Born treatment, which is the standard alternative for exactly this
class of problem.

⚠ Where MM-GBSA is itself weak is charge CHANGES -- but the decisive comparison here,
**K59 vs R59, is charge-preserving** (+1 either way), so the largest error terms
cancel between the two states. That is a favourable case, not a lucky one, and it is
the reason this test is worth running rather than a general endorsement of MM-GBSA.

PRE-REGISTERED READING, fixed before any number is computed
-----------------------------------------------------------
The ligand-swap rule from README 23 still governs: only the DIFFERENCE between arms
is evidence.

    SUCCESS  ddG(K59R) favourable for mandipropamid AND unfavourable for ABA
    FAILURE  same sign for both ligands -- the model is reading charge, not
             complementarity, and MM-GBSA inherits the problem in a new coat
    NULL     |ddG| within the seed-to-seed spread; underpowered, report as such

The spread across independent repack seeds IS the error bar. A single-structure
MM-GBSA number for a charged residue is not interpretable on its own, so no result
here is reported without it.

HOW THE STRUCTURES ARE MADE
---------------------------
Backbone fixed, ligand fixed, only the mutated residue and its shell repacked --
with `lib_rosetta.restrict_packing`, because a MoveMap alone does not restrain the
packer (README 28d). Rosetta is used to BUILD, never to score: every energy reported
downstream comes from Amber.

Run with the tier1_analysis env python.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "data", "mmgbsa")
os.makedirs(OUT, exist_ok=True)

import lib_rosetta as lr            # noqa: E402

ARMS = {
    "mandi": (os.path.join(ROOT, "results", "stage1_rosetta", "_input_wt_mandi_3UZ.pdb"),
              os.path.join(ROOT, "data", "stage1", "params", "3UZ.params")),
    "aba":   (os.path.join(ROOT, "results", "stage1_rosetta", "_input_wt_aba_A8S_anion.pdb"),
              os.path.join(ROOT, "data", "stage1", "params", "A8S_anion.params")),
}

#: WT plus the ground truth, plus the two neutral substitutions Beltran's WIN
#: sensors actually use at 59 (K59Q in 4 sensors, K59N in 4). Those are the
#: substitutions ref2015 SHOULD find easy -- no buried charge -- so they are the
#: internal check that the comparison is working at all.
VARIANTS = {
    "WT":    {},
    "K59R":  {59: "ARG"},
    "K59Q":  {59: "GLN"},
    "K59N":  {59: "ASN"},
    "V81I":  {81: "ILE"},
    "F108A": {108: "ALA"},
    "F159L": {159: "LEU"},
}
N_SEED = 8
SHELL = 6.0


def log(m):
    print(m, flush=True)


def build(arm, variant, seeds=N_SEED):
    import pyrosetta
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

    pdb, params = ARMS[arm]
    muts = VARIANTS[variant]
    d = os.path.join(OUT, f"{arm}_{variant}")
    os.makedirs(d, exist_ok=True)

    for seed in range(seeds):
        pyrosetta.init(f"-mute all -extra_res_fa {params} -detect_disulf false "
                       f"-load_PDB_components false -ex1 -ex2 "
                       f"-run:constant_seed -run:jran {20260819 + seed}")
        pose = pyrosetta.pose_from_pdb(pdb)
        pi = pose.pdb_info()
        idx = {pi.number(i): i for i in range(1, pose.total_residue() + 1)
               if pose.residue(i).is_protein()}
        # identity before mutation, never position alone
        for num, want in ((59, "LYS"), (81, "VAL"), (108, "PHE"), (159, "PHE")):
            got = pose.residue(idx[num]).name3().strip()
            assert got == want, f"{arm}: residue {num} is {got}, expected {want}"
        for num, aa in muts.items():
            MutateResidue(idx[num], aa).apply(pose)

        touched = sorted(idx[n] for n in muts) or [idx[59]]
        sel = ResidueIndexSelector(",".join(map(str, touched)))
        shell = NeighborhoodResidueSelector(sel, SHELL, True).apply(pose)
        allow = [i for i in range(1, pose.total_residue() + 1)
                 if shell[i] and pose.residue(i).is_protein()]
        tf, packable = lr.restrict_packing(pose, allow)
        sf = pyrosetta.create_score_function("ref2015")
        PackRotamersMover(sf, tf.create_task_and_apply_taskoperations(pose)).apply(pose)

        out = os.path.join(d, f"seed{seed}.pdb")
        pose.dump_pdb(out)
        # protein only -- the ligand comes from its own mol2, already in this frame
        prot = os.path.join(d, f"seed{seed}_protein.pdb")
        with open(prot, "w") as fh:
            for line in open(out):
                if line.startswith(("ATOM", "TER")) and line[17:20].strip() != "LIG":
                    fh.write(line)
            fh.write("END\n")
        os.remove(out)
        if seed == 0:
            log(f"  {arm}/{variant}: {len(packable)} residues repacked, {seeds} seeds")
    return d


def main():
    if len(sys.argv) > 2:
        build(sys.argv[1], sys.argv[2])
        return
    made = []
    for arm in ARMS:
        for v in VARIANTS:
            build(arm, v)
            made.append(f"{arm}_{v}")
    json.dump(dict(arms=list(ARMS), variants=VARIANTS, n_seed=N_SEED, built=made),
              open(os.path.join(OUT, "build_manifest.json"), "w"), indent=1)
    log(f"\nbuilt {len(made)} variant/arm combinations x {N_SEED} seeds")


if __name__ == "__main__":
    main()
