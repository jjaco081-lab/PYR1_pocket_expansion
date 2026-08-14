#!/usr/bin/env python
"""
66_mutant_receptors.py -- build the POSITIVE controls: each non-cognate ligand
docked into the evolved PYR1 mutant that actually senses it.

WHY
---
§27's negative control asks whether closed PYR1 holds its state around a ligand it
was never built for. On its own that is only half a test: if the dynamics look the
same as with ABA, we cannot tell whether the readout is blind to ligands or whether
these particular complexes happen to be fine. The positive control supplies the
other half -- the SAME ligand in the receptor Tian isolated for it, which
demonstrably works in yeast.

    ligand            sensor clone   min_conc   mutations
    Imperatorin       54_6  (tsm)      10 uM    K59A V163W V164M
    Flutamide         21_6  (tsm)      10 uM    V83I A160M V164M
    Alpha-Estradiol   124_10 (tsm)      1 uM    F159V V163W V164G

Chosen as the LOWEST min_conc clone for each ligand, i.e. the most sensitive sensor
Tian characterised. All three carry exactly three mutations, so the comparison is
not confounded by how much the receptor was changed. Provenance: Dataset S3
(`pnas.2519924122.sd03(1).xlsx`), read in script 62's export.

The contrast that matters is then WT+ligand vs MUTANT+ligand, with ligand held
fixed -- and ABA+WT (S2) and apo-closed (S9) as the anchors.

⚠️ RESIDUE NUMBERING -- CHECKED, NOT ASSUMED
--------------------------------------------
Tian quotes mutations in NATIVE PYR1 numbering. `data/md/S2_holo_closed/protein.pdb`
also carries native numbering (with gaps: residue 2 is absent), which is NOT the
same as the sequential 1..178 that tleap will later assign. So mutations are applied
in protein.pdb space, and every position is asserted to hold the expected WT residue
BEFORE anything is mutated. A silent off-by-one here would build the wrong protein
and nothing downstream would complain.

Run with the tier1_analysis env python (PyRosetta). Submit via 66b, not interactively
-- PyRosetta peaks near 900 MB and a plain srun is capped at 1 GB.
"""
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
NC = os.path.join(DATA, "noncognate")
OUT = os.path.join(NC, "mutants")
os.makedirs(OUT, exist_ok=True)

WT_PDB = os.path.join(DATA, "md", "S2_holo_closed", "protein.pdb")
ABA_REF = os.path.join(DATA, "md", "S2_holo_closed", "aba.pdb")

# native PYR1 numbering, straight from the sensor clones
SENSORS = {
    "Imperatorin":     dict(clone="54_6",   min_conc=10.0,
                            muts=[(59, "K", "A"), (163, "V", "W"), (164, "V", "M")]),
    "Flutamide":       dict(clone="21_6",   min_conc=10.0,
                            muts=[(83, "V", "I"), (160, "A", "M"), (164, "V", "M")]),
    "Alpha-Estradiol": dict(clone="124_10", min_conc=1.0,
                            muts=[(159, "F", "V"), (163, "V", "W"), (164, "V", "G")]),
}

THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
         "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
         "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
         "W": "TRP", "Y": "TYR"}

import pyrosetta  # noqa: E402
from pyrosetta import pose_from_pdb  # noqa: E402
from pyrosetta.rosetta.core.pack.task import TaskFactory  # noqa: E402
from pyrosetta.rosetta.core.pack.task.operation import (  # noqa: E402
    RestrictAbsentCanonicalAASRLT, OperateOnResidueSubset,
    PreventRepackingRLT, RestrictToRepackingRLT)
from pyrosetta.rosetta.core.select.residue_selector import (  # noqa: E402
    ResidueIndexSelector, NeighborhoodResidueSelector, NotResidueSelector)
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory  # noqa: E402
from pyrosetta.rosetta.protocols.minimization_packing import (  # noqa: E402
    PackRotamersMover)

pyrosetta.init("-mute all -ex1 -ex2", silent=True)

pose = pose_from_pdb(WT_PDB)
pdb_info = pose.pdb_info()
print(f"loaded {WT_PDB}: {pose.total_residue()} residues")

# native resnum -> Rosetta pose index, built from the PDB numbering the file carries
nat2pose = {}
for i in range(1, pose.total_residue() + 1):
    nat2pose[pdb_info.number(i)] = i

print("\nresidue-identity check (native numbering, before any mutation):")
ok = True
for lig, spec in SENSORS.items():
    for nat, wt1, new1 in spec["muts"]:
        if nat not in nat2pose:
            print(f"  !! native {nat} absent from {WT_PDB}")
            ok = False
            continue
        idx = nat2pose[nat]
        got = pose.residue(idx).name3()
        want = THREE[wt1]
        flag = "ok" if got == want else "MISMATCH"
        if got != want:
            ok = False
        print(f"  {lig:<16} native {nat:>4} (pose {idx:>4}): {got} expected {want}   {flag}")
assert ok, ("a mutation position does not hold the residue Tian's table says it "
            "does -- the numbering is wrong; do NOT build anything from this")
print("  -> every position verified")

sfxn = ScoreFunctionFactory.create_score_function("ref2015")
prot_heavy = None


def mutate(base, muts, label):
    """apply mutations, repack an 8 A shell, return the new pose"""
    p = base.clone()
    idx = [nat2pose[n] for n, _, _ in muts]
    sel = ResidueIndexSelector(",".join(str(i) for i in idx))

    tf = TaskFactory()
    # design only at the target positions, to the intended residue
    for (nat, _, new1), i in zip(muts, idx):
        one = ResidueIndexSelector(str(i))
        restrict = RestrictAbsentCanonicalAASRLT()
        restrict.aas_to_keep(new1)
        tf.push_back(OperateOnResidueSubset(restrict, one))
    shell = NeighborhoodResidueSelector(sel, 8.0, False)
    tf.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), shell))
    outside = NotResidueSelector(NeighborhoodResidueSelector(sel, 8.0, True))
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), outside))

    packer = PackRotamersMover(sfxn, tf.create_task_and_apply_taskoperations(p))
    packer.apply(p)

    for nat, _, new1 in muts:
        got = p.residue(nat2pose[nat]).name3()
        assert got == THREE[new1], (
            f"{label}: native {nat} is {got} after packing, wanted {THREE[new1]}")
    return p


records = {}
for lig, spec in SENSORS.items():
    label = "_".join(f"{w}{n}{m}" for n, w, m in spec["muts"])
    print(f"\n=== {lig}: {label}  (clone {spec['clone']}, {spec['min_conc']} uM) ===")
    mp = mutate(pose, spec["muts"], lig)
    out_pdb = os.path.join(OUT, f"{lig}_receptor.pdb")
    mp.dump_pdb(out_pdb)

    # the mutated file must still carry NATIVE numbering for the next stage
    got_nums = set()
    for l in open(out_pdb):
        if l.startswith("ATOM"):
            got_nums.add(int(l[22:26]))
    for nat, _, _ in spec["muts"]:
        assert nat in got_nums, f"{lig}: native {nat} missing from {out_pdb}"
    print(f"  wrote {out_pdb}")
    records[lig] = dict(clone=spec["clone"], min_conc=spec["min_conc"],
                        mutations=label, pdb=out_pdb)

    records[lig]["sdf"] = os.path.join(NC, f"{lig}.sdf")

json.dump(records, open(os.path.join(OUT, "sensors.json"), "w"), indent=1)
print(f"\nwrote {OUT}/sensors.json")
print("NEXT: 63_dock_noncognate.py --mutants   (docking lives there, in the rdkit env)")
