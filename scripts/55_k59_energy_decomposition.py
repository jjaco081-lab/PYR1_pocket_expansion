#!/usr/bin/env python
"""
55_k59_energy_decomposition.py -- why does ref2015 delete the ABA-K59 salt bridge?

THE OBSERVATION THIS EXPLAINS
Job 27421344's ABA null arm is WT protein with its native ligand, so its correct
answer is zero mutations. It replaced K59 with Asn in 82% of trajectories (Ile in
the rest) at EVERY favor-native weight tested, including 1.5, where 13 of the 15
designable positions were otherwise held at WT (job 27438220).

WHAT HAS ALREADY BEEN RULED OUT
  - Pose. K59 NZ sits 2.85 A from ABA carboxylate O4 and 3.02 A from O3: a textbook
    bidentate salt bridge, measured directly from the Rosetta input.
  - Ionisation. A8S_anion.params partial charges sum to -0.970 over 38 atoms, i.e.
    the intended -1 to within 2-dp rounding. (An earlier reading of "exactly 0.000"
    was an off-by-one in the reader -- field 4 of an ATOM line is the MM type `X`,
    not the charge. 49_stage1_ligand_params.py now validates this properly, and the
    regenerated params are byte-identical to the ones the campaign used.)
  - A sequence prior being too weak. favor-native up to 1.5 REU/residue does not
    recover K59, so the deficit is not "ref2015 mildly prefers something else".

WHAT THIS SCRIPT MEASURES
Score the WT pose and each single substitution at 59 with ref2015, repacking only
the shell around 59 so the comparison is like-for-like, and decompose the per-residue
energy into the terms that decide a buried salt bridge:
    fa_elec  -- Coulomb, distance-damped
    fa_sol   -- Lazaridis-Karplus desolvation, the term that punishes burying charge
    fa_rep / fa_atr / hbond_sc
If Asn beats Lys mainly on fa_sol, this is the documented ref2015 buried-salt-bridge
pathology and the Rosetta arm simply cannot speak to K59R -- a limitation to report,
not a bug to fix. If Lys wins on the total, the 100% deletion is a sampling artifact
and the protocol is at fault instead.

NOTE ON WHAT WOULD *NOT* BE LEGITIMATE
Upweighting fa_elec until K59 survives would be tuning the scorefunction until the
control agrees, which is the same error as widening the favor-native sweep. The
output here is a diagnosis, not a knob.

Usage (needs ~2 GB, so submit it -- an interactive srun defaults to 1 GB):
  sbatch scripts/55b_submit_k59.sh
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARAMS = os.path.join(ROOT, "data", "stage1", "params", "A8S_anion.params")
INP = os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")
LIGPDB = os.path.join(ROOT, "data", "stage1", "params", "A8S_anion_0001.pdb")
OUT = os.path.join(ROOT, "results", "stage1_rosetta", "k59_decomposition.json")

POS = 59
# what the null actually chose, plus the ground-truth answer and a charge control
CANDIDATES = ["K", "N", "I", "R", "Q", "M"]
SHELL = 8.0

import pyrosetta                                                    # noqa: E402
pyrosetta.init(f"-extra_res_fa {PARAMS} "
               "-ex1 -ex2 -use_input_sc -no_optH false -flip_HNQ "
               "-ignore_zero_occupancy false -ignore_unrecognized_res false "
               "-mute all")

import numpy as np                                                  # noqa: E402
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory, ScoreType  # noqa: E402
from pyrosetta.rosetta.core.pack.task import TaskFactory, operation  # noqa: E402
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover  # noqa: E402

sfxn = ScoreFunctionFactory.create_score_function("ref2015")

# Build the same input the campaign used: protein from wt_aba.pdb, ligand spliced
# from the molfile_to_params output so atom names match the params.
lig_block = [l for l in open(LIGPDB) if l.startswith(("ATOM  ", "HETATM"))]
tmp = os.path.join(os.path.dirname(OUT), "_k59_input.pdb")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(tmp, "w") as fh:
    for l in open(INP):
        if l.startswith(("ATOM  ", "HETATM")):
            fh.write(l)
    for l in lig_block:
        fh.write(l)
    fh.write("END\n")

pose0 = pyrosetta.pose_from_pdb(tmp)
pi = pose0.pdb_info()
idx = pi.pdb2pose("A", POS)
assert idx > 0, f"residue A{POS} not in pose"
lig = [i for i in range(1, pose0.size() + 1) if pose0.residue(i).name3() == "LIG"]
assert len(lig) == 1, f"expected 1 LIG, got {len(lig)}"
LIG = lig[0]


def heavy(pose, i):
    r = pose.residue(i)
    return np.array([[*r.xyz(a)] for a in range(1, r.natoms() + 1)
                     if not r.atom_is_hydrogen(a)])


def shell_of(pose, i, cutoff):
    """All-heavy-atom distance, not neighbour-atom -- a 19-atom ligand breaks NBR."""
    X = heavy(pose, i)
    out = []
    for j in range(1, pose.size() + 1):
        if j == i:
            continue
        Y = heavy(pose, j)
        if len(Y) and np.linalg.norm(X[:, None] - Y[None], axis=-1).min() <= cutoff:
            out.append(j)
    return out


NEIGH = shell_of(pose0, idx, SHELL)
print(f"position A{POS} -> pose {idx}; ligand {LIG}; "
      f"repack shell {len(NEIGH)} residues within {SHELL:.0f} A", flush=True)

# the salt-bridge geometry, restated from the pose Rosetta actually scores
nz = pose0.residue(idx).xyz("NZ")
for at in ("O3", "O4"):
    d = (np.array([*nz]) - np.array([*pose0.residue(LIG).xyz(at)]))
    print(f"  K{POS} NZ -- LIG {at}  {np.linalg.norm(d):.2f} A")

TERMS = {"fa_elec": ScoreType.fa_elec, "fa_sol": ScoreType.fa_sol,
         "fa_rep": ScoreType.fa_rep, "fa_atr": ScoreType.fa_atr,
         "hbond_sc": ScoreType.hbond_sc}

Sel = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector
FROZEN = [i for i in range(1, pose0.size() + 1) if i != idx and i not in NEIGH]

rows = {}
for aa in CANDIDATES:
    pose = pose0.clone()
    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())

    # freeze everything outside the shell
    tf.push_back(operation.OperateOnResidueSubset(
        operation.PreventRepackingRLT(),
        Sel(",".join(str(i) for i in FROZEN))))
    # repack, but do not design, the shell
    tf.push_back(operation.OperateOnResidueSubset(
        operation.RestrictToRepackingRLT(),
        Sel(",".join(str(i) for i in NEIGH))))
    # position 59 is forced to exactly this one amino acid
    rlt = operation.RestrictAbsentCanonicalAASRLT()
    rlt.aas_to_keep(aa)
    tf.push_back(operation.OperateOnResidueSubset(rlt, Sel(str(idx))))

    pr = PackRotamersMover(sfxn)
    pr.task_factory(tf)
    pr.apply(pose)

    assert pose.residue(idx).name1() == aa, (
        f"packer produced {pose.residue(idx).name1()} at {idx}, expected {aa}")

    total = float(sfxn(pose))
    e = pose.energies().residue_total_energies(idx)
    row = {"total_pose": round(total, 3),
           "res59_total": round(float(pose.energies().residue_total_energy(idx)), 3)}
    for nm, t in TERMS.items():
        row[nm] = round(float(e[t]), 3)
    rows[aa] = row
    print(f"  {aa}: pose {total:10.3f}   res59 {row['res59_total']:8.3f}", flush=True)

ref = rows["K"]["total_pose"]
print("\n" + "=" * 76)
print(f"ref2015, position {POS}, relative to WT Lys (negative = better than Lys)")
print("=" * 76)
print(f"{'aa':>3} {'d_pose':>9} {'res59':>9} {'fa_elec':>9} {'fa_sol':>9} "
      f"{'fa_rep':>8} {'fa_atr':>8} {'hbond_sc':>9}")
for aa, r in rows.items():
    print(f"{aa:>3} {r['total_pose']-ref:>9.3f} {r['res59_total']:>9.3f} "
          f"{r.get('fa_elec',0):>9.3f} {r.get('fa_sol',0):>9.3f} "
          f"{r.get('fa_rep',0):>8.3f} {r.get('fa_atr',0):>8.3f} "
          f"{r.get('hbond_sc',0):>9.3f}")

json.dump(rows, open(OUT, "w"), indent=1)
print(f"\nwrote {OUT}")
