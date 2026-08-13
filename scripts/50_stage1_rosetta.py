#!/usr/bin/env python
"""
50_stage1_rosetta.py -- stage 1, second opinion: can Rosetta FastDesign recover
PYR1^MANDI from WT plus a perfectly placed ligand?

WHY A SECOND METHOD IS REQUIRED, NOT OPTIONAL
---------------------------------------------
Stage 1's pre-registered gate (section 23f) says a single method's failure is not
a statement about the architecture. A learned designer and a physics-based one
fail for unrelated reasons -- LigandMPNN can only propose what its training
distribution supports, while Rosetta can only find what its rotamer set and
energy function reward -- so agreement between them is informative in a way that
either alone is not. Rosetta is also the baseline Leonard et al. used, which
makes this arm comparable to the published work rather than only to itself.

GROUND TRUTH (4WVO vs 3QN1, 174 residues resolved in both): K59R V81I F108A F159L

WHAT DIFFERS FROM THE LigandMPNN ARM
------------------------------------
  * Poly-Gly is a REAL arm here. LigandMPNN masks designable side chains
    regardless, so truncating them changed nothing it could see and
    `polygly_mandi` was a duplicate of `wt_mandi` (section 23g). Rosetta packs
    explicit rotamers, so removing the WT side chains genuinely removes the
    clash signal and forces the pocket to be rebuilt rather than trimmed.
  * Rosetta returns structures and energies, not a probability distribution.
    The analogue of LigandMPNN's `sampling_probs` is the FREQUENCY of the true
    residue across independent trajectories, which is why N is 50 rather than 1.
  * An interface energy (bound minus separated) is recorded per trajectory. It
    is the quantity that would actually rank library members, and unlike a raw
    total score it is comparable across sequences of different composition.

SCORING -- unchanged from section 23f, deliberately
---------------------------------------------------
The primary statistic is still the ligand-swap delta:

    delta = freq_mandi(true residue) - freq_aba(true residue)

The ABA arm runs an identical protocol with only the ligand changed, so it
inherits whatever this method's intrinsic mutational bias is, and no assumption
about independence of trajectories is needed. A high frequency that the ABA arm
matches is a ligand-independent preference, not recovery.

Any result must also beat the trivial baseline from 47_stage1_inputs.py: ranking
WT side chains by steric overlap with the ligand already identifies 3 of the 4
true POSITIONS. Position finding is nearly free; identity selection is the task.

THE PREMISE, AND ITS LIMIT
--------------------------
Rosetta optimises a fixed backbone against a fixed ligand pose, so it can only
relieve the clash it is handed. That is the same perfect-pose assumption
LigandMPNN was given, which is what keeps the comparison fair -- but neither
method is being asked stage 2's harder question, where the pose is unknown. A
success here is necessary for the cascade to continue, not sufficient for it.

Ligand internal geometry and the ligand jump are both held fixed: the crystal
conformer IS the hypothesis under test (scripts 45/46), so letting Rosetta
minimise it away would answer a different question.

Usage (one arm-block per invocation; see 50b_submit_stage1_rosetta.sh):
  python 50_stage1_rosetta.py --unit <0..29> [--n-per-unit 10]
  (tier1_analysis conda env, which is the only one with PyRosetta)
"""
import argparse, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
IN = os.path.join(ROOT, "data", "stage1")
PARAMS = os.path.join(IN, "params")
OUT = os.path.join(ROOT, "results", "stage1_rosetta")

TRUE = {59: ("K", "R"), 81: ("V", "I"), 108: ("F", "A"), 159: ("F", "L")}
DESIGN = [59, 81, 83, 92, 94, 108, 110, 120, 122, 141, 159, 160, 163, 164, 167]
DSM_HAO = {
    59: "ADEFGHILMNQRSTVWY", 81: "ILMRTY", 83: "AFGILMSTVWY",
    92: "ADEFGHIKLMNQRSTVWY", 94: "ADFGHIKLMNQRSTVWY",
    108: "ADEGHIKLMNQRSTVWY", 110: "AFGILMSTVWY", 120: "ADEFGHIKLMNQRSTVW",
    122: "ADEFGHIKLMNQRTVWY", 141: "ADEFGHIKLMNQRSTVW",
    159: "ADEGHIKLMNQRSTWY", 160: "DEFGHIKLMNQRSTVWY",
    163: "ADEFGHIKLMNQRSTWY", 164: "ADEFGHIKLMNQRSTWY",
    167: "ADEFGHIKLMQRSTVWY",
}
ALL20 = "ACDEFGHIKLMNPQRSTVWY"

# arm -> (source pdb, ligand CCD code, params stem)
ARMS = [
    ("wt_mandi", "3UZ", "3UZ"),
    ("polygly_mandi", "3UZ", "3UZ"),
    ("wt_aba", "A8S", "A8S_anion"),
]
ALPHABETS = ["dsm_hao", "free"]
UNITS_PER_ARM = 5                       # 5 blocks x 10 trajectories = 50 per arm

ap = argparse.ArgumentParser()
ap.add_argument("--unit", type=int, required=True)
ap.add_argument("--n-per-unit", type=int, default=10)
ap.add_argument("--outdir", default=OUT)
ap.add_argument("--dump-pdb", action="store_true",
                help="write the designed structure for every trajectory")
ap.add_argument("--pack-shell", type=float, default=8.0,
                help="repack residues within this many A of the ligand; "
                     "everything further out is frozen (default 8)")
ap.add_argument("--rounds", type=int, default=1,
                help="FastRelax ramp rounds (default 1)")
ap.add_argument("--favor-native", type=float, default=0.0,
                help="per-residue res_type_constraint bonus for the starting "
                     "identity (0 = off, the behaviour of job 27421344). See the "
                     "WHY block below before changing this.")
args = ap.parse_args()

n_arms = len(ARMS) * len(ALPHABETS)
assert 0 <= args.unit < n_arms * UNITS_PER_ARM, f"--unit must be 0..{n_arms*UNITS_PER_ARM-1}"
arm_idx, block = divmod(args.unit, UNITS_PER_ARM)
(struct, comp, stem), alphabet = (ARMS[arm_idx // len(ALPHABETS)],
                                  ALPHABETS[arm_idx % len(ALPHABETS)])
tag = f"{struct}__{alphabet}"
os.makedirs(args.outdir, exist_ok=True)
if args.dump_pdb:
    os.makedirs(os.path.join(args.outdir, "pdb"), exist_ok=True)

params_file = os.path.join(PARAMS, f"{stem}.params")
assert os.path.exists(params_file), (
    f"missing {params_file} -- run 49_stage1_ligand_params.py first")

# ---------------------------------------------------------------- input prep
def heavy_xyz(lines):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in lines if l[76:78].strip() != "H"])


def rosetta_input(struct, comp, stem):
    """Protein from the stage-1 PDB + ligand from molfile_to_params' own output.

    The ligand block is taken from `<stem>_0001.pdb` rather than from the
    stage-1 HETATM records, because molfile_to_params RENAMES atoms when the
    input molfile has duplicate or absent names (it does, since MDL molfiles
    carry no atom names). Rosetta matches PDB atoms to the residue type by NAME,
    so feeding it the crystal names against a params file written with C1/C2/...
    leaves every atom unmatched and dies with "too many tries in
    fill_missing_atoms". That file also already carries the hydrogens the params
    expects, on chain X as residue LIG 1.

    Its coordinates are the crystal coordinates -- 49_stage1_ligand_params.py
    verifies that one-to-one to 0.01 A -- but it is re-verified HERE against the
    stage-1 file, because the two scripts can drift independently and a ligand
    that is present but in the wrong place would look like a clean negative
    result rather than an error. Section 23h is what that costs.
    """
    src = os.path.join(IN, f"{struct}.pdb")
    prot = [l for l in open(src) if l.startswith("ATOM")]
    crystal = [l for l in open(src) if l.startswith("HETATM")]
    assert crystal, f"{src} has no HETATM records"
    for l in crystal:
        assert l[16] == " " and l[17:20].strip() == comp, (
            f"{src}: malformed ligand columns (altLoc={l[16]!r}, "
            f"resName={l[17:20]!r}) -- re-run 47_stage1_inputs.py")

    lig_pdb = os.path.join(PARAMS, f"{stem}_0001.pdb")
    assert os.path.exists(lig_pdb), (
        f"missing {lig_pdb} -- run 49_stage1_ligand_params.py first")
    lig = [l for l in open(lig_pdb) if l.startswith(("ATOM", "HETATM"))]
    assert lig, f"{lig_pdb} has no atom records"

    a, b = heavy_xyz(crystal), heavy_xyz(lig)
    assert a.shape == b.shape, (
        f"{stem}: params ligand has {len(b)} heavy atoms, crystal has {len(a)}")
    d = np.linalg.norm(b[:, None, :] - a[None, :, :], axis=-1)
    j = d.argmin(1)
    worst = float(d[np.arange(len(b)), j].max())
    assert len(set(j.tolist())) == len(b) and worst <= 0.01, (
        f"{stem}: ligand pose differs from the crystal by {worst:.4f} A -- "
        f"the params conformer is not the bound conformer")

    # per-unit filename: every array task builds its own copy, and 30 tasks
    # writing one shared path raced badly enough to leave a 16 MB file (the
    # asserts above caught no bad read, but only by luck of timing)
    dst = os.path.join(args.outdir, f"_input_u{args.unit:02d}_{struct}_{stem}.pdb")
    with open(dst, "w") as fh:
        fh.writelines(prot)
        fh.writelines(lig)
        fh.write("END\n")
    return dst, len(a), worst


inp, n_lig_atoms, pose_dev = rosetta_input(struct, comp, stem)

import pyrosetta                                                    # noqa: E402
from pyrosetta.rosetta.core.pack.task import TaskFactory            # noqa: E402
from pyrosetta.rosetta.core.pack.task import operation              # noqa: E402
from pyrosetta.rosetta.core.select.residue_selector import (        # noqa: E402
    ResidueIndexSelector, NotResidueSelector, OrResidueSelector,
    NeighborhoodResidueSelector)
from pyrosetta.rosetta.core.scoring import (                        # noqa: E402
    ScoreFunctionFactory, ScoreType)
from pyrosetta.rosetta.protocols.relax import FastRelax             # noqa: E402
# NOT protocols.simple_moves -- in PyRosetta 2026.06 FavorNativeResidue lives in
# protocols.protein_interface_design. Verified by import before submitting.
from pyrosetta.rosetta.protocols.protein_interface_design import (  # noqa: E402
    FavorNativeResidue)
from pyrosetta.rosetta.numeric import xyzVector_double_t            # noqa: E402

pyrosetta.init(
    f"-extra_res_fa {params_file} "
    "-ex1 -ex2 -use_input_sc -no_optH false -flip_HNQ "
    "-relax:constrain_relax_to_start_coords "
    "-relax:coord_constrain_sidechains false "
    "-relax:ramp_constraints false "
    "-ignore_zero_occupancy false -ignore_unrecognized_res false -mute all")

sfxn_cst = ScoreFunctionFactory.create_score_function("ref2015_cst")
sfxn = ScoreFunctionFactory.create_score_function("ref2015")

# WHY --favor-native EXISTS (added 2026-08-13, after job 27421344)
#
# Job 27421344 ran with no favor-native term, i.e. bare ref2015 reference energies.
# The ABA arm is WT protein with its native ligand, so its correct answer is ZERO
# mutations -- it is the WT-recovery control. It kept WT at only 4 of 15 positions
# and deleted the K59 carboxylate salt bridge in 100% of trajectories (to Asn 82%,
# Ile 18%), despite ABA being modelled as the anion, the protonation state most
# favourable to that salt bridge. ref2015's reference energies are fit for soluble
# monomer design and do not hold a native complex.
#
# The consequence for the benchmark: at a position the null also mutates, "no
# ligand-conditional signal" cannot be told apart from "the protocol cannot hold a
# native contact". That is what made K59R and V81I uninterpretable rather than
# negative. F108A survived only because the null happens to be quiet there.
#
# res_type_constraint scores the FavorNativeResidue bonus. It is set on sfxn_cst
# (used for design) but deliberately NOT on sfxn (used for interface energy), so
# reported energies stay on the same scale as job 27421344 and remain comparable.
if args.favor_native > 0:
    sfxn_cst.set_weight(ScoreType.res_type_constraint, 1.0)
assert sfxn.get_weight(ScoreType.res_type_constraint) == 0.0, (
    "interface energy must not include the favor-native bonus")

start = pyrosetta.pose_from_pdb(inp)
pi = start.pdb_info()

# the ligand must have survived the load, with all its atoms
lig_idx = [i for i in range(1, start.size() + 1)
           if start.residue(i).name3() == "LIG"]
assert len(lig_idx) == 1, f"expected exactly 1 LIG residue, found {len(lig_idx)}"
LIG = lig_idx[0]
n_lig_heavy = sum(1 for a in range(1, start.residue(LIG).natoms() + 1)
                  if not start.residue(LIG).atom_is_hydrogen(a))
assert n_lig_heavy == n_lig_atoms, (
    f"ligand loaded with {n_lig_heavy} heavy atoms, input had {n_lig_atoms}")

# Applied to `start`, so every cloned trajectory inherits the same constraints.
# NOTE for the poly-Gly arm: "native" here means the identity in the INPUT pose,
# which for polygly_mandi is Gly at all 15 designable positions -- so a non-zero
# bonus there rewards staying Gly, which is not what the arm is for. The sweep
# therefore calibrates on wt_aba only, and the poly-Gly arm is re-run at whatever
# weight the null selects purely so all six arms stay on one protocol.
if args.favor_native > 0:
    FavorNativeResidue(start, args.favor_native)

pose_idx = {}
for p in DESIGN:
    i = pi.pdb2pose("A", p)
    assert i > 0, f"PDB residue A{p} absent from the pose"
    pose_idx[p] = i
print(f"{tag} block {block}: {start.size()} residues, ligand {LIG} "
      f"({n_lig_heavy} heavy atoms, pose dev {pose_dev:.4f} A), "
      f"{len(pose_idx)} designable", flush=True)


def res_heavy(pose, i):
    r = pose.residue(i)
    return np.array([list(r.xyz(a)) for a in range(1, r.natoms() + 1)
                     if not r.atom_is_hydrogen(a)])


def pocket_shell(pose, lig, cutoff):
    """Residues with any heavy atom within `cutoff` A of any ligand heavy atom.

    Computed explicitly rather than with NeighborhoodResidueSelector, which
    measures between residue NEIGHBOUR atoms (CB, or a single nominated atom for
    a ligand). molfile_to_params nominates one NBR atom for the whole ligand, so
    an 8 A neighbour-atom sphere around a 12 A-long molecule selected just THREE
    residues here -- the pocket was effectively frozen solid and every
    trajectory converged to the same answer. An all-heavy-atom criterion is
    what "lines the pocket" actually means, and it is worth the extra loop to
    have it be the thing the name says.
    """
    L = res_heavy(pose, lig)
    out = []
    for i in range(1, pose.size() + 1):
        if i == lig:
            continue
        X = res_heavy(pose, i)
        if len(X) and np.linalg.norm(X[:, None, :] - L[None, :, :],
                                     axis=-1).min() <= cutoff:
            out.append(i)
    return out


def task_factory():
    """Design at DESIGN; repack the pocket shell; freeze everything else.

    WT identity is always added to the allowed set. Without it a position whose
    WT residue is absent from the library alphabet could not stay WT, which
    would force a mutation the experiment never forced and inflate apparent
    recovery at every such position.

    WHY THE PACKER IS RESTRICTED TO A SHELL. The first version let every residue
    in the protein repack. One trajectory did not finish in 28 minutes, because
    -ex1 -ex2 over 174 repackable positions plus 15 designable ones, through
    3 FastRelax rounds, is an enormous rotamer problem. Restricting the packer to
    residues within --pack-shell A of the ligand is also the better experiment,
    not merely the cheaper one: the question is whether the POCKET can be
    rebuilt around mandipropamid, and letting surface side chains 30 A away flip
    between trajectories injects variance into the sequence frequencies that has
    nothing to do with the ligand. Both arms use the identical shell definition,
    so the ligand-swap null is unaffected either way.

    The ligand itself is frozen -- see movemap().
    """
    lig_sel = ResidueIndexSelector(str(LIG))
    design_sel = ResidueIndexSelector(",".join(str(pose_idx[p]) for p in DESIGN))
    idx = sorted(set(pocket_shell(start, LIG, args.pack_shell))
                 | set(pose_idx.values()))
    mobile = ResidueIndexSelector(",".join(str(i) for i in idx))

    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.IncludeCurrent())
    # outside the shell, and the ligand: no repacking at all
    tf.push_back(operation.OperateOnResidueSubset(
        operation.PreventRepackingRLT(), NotResidueSelector(mobile)))
    tf.push_back(operation.OperateOnResidueSubset(
        operation.PreventRepackingRLT(), lig_sel))
    # inside the shell but not designable: repack only
    tf.push_back(operation.OperateOnResidueSubset(
        operation.RestrictToRepackingRLT(), NotResidueSelector(design_sel)))
    for p in DESIGN:
        wt = start.residue(pose_idx[p]).name1()
        aas = (set(DSM_HAO.get(p, ALL20)) if alphabet == "dsm_hao"
               else set(ALL20)) | {wt}
        rlt = operation.RestrictAbsentCanonicalAASRLT()
        rlt.aas_to_keep("".join(sorted(aas)))
        tf.push_back(operation.OperateOnResidueSubset(
            rlt, ResidueIndexSelector(str(pose_idx[p]))))
    return tf, [i for i in idx if i != LIG]


def movemap():
    """Backbone and the ligand jump both fixed; only mobile side chains minimise.

    The premise of stage 1 is a perfectly placed ligand on the WT backbone. If
    the ligand were allowed to translate or the backbone to move, a failure to
    recover the mutations could always be blamed on the pose drifting, and a
    success could come from Rosetta relocating the ligand to somewhere WT
    already accommodates.

    Chi is enabled only where the packer is allowed to act, so minimisation
    cannot quietly relax the 100+ residues the task factory froze -- that would
    put back most of the cost the shell restriction removes, and would let the
    total score drift for reasons unrelated to the pocket.
    """
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    mm.set_jump(False)
    for i in MOBILE:
        mm.set_chi(i, True)
    return mm


def interface_energy(pose):
    """ref2015 bound minus separated: the ligand's contribution to the score.

    Composition cancels because both terms use the same sequence, which a raw
    total score does not give -- ref2015 is extensive and its reference energies
    differ per amino acid, so total scores across different designed sequences
    are not directly comparable.
    """
    bound = float(sfxn(pose))
    apart = pose.clone()
    v = xyzVector_double_t(500.0, 0.0, 0.0)
    for a in range(1, apart.residue(LIG).natoms() + 1):
        apart.set_xyz(pyrosetta.rosetta.core.id.AtomID(a, LIG),
                      apart.residue(LIG).xyz(a) + v)
    return bound, bound - float(sfxn(apart))


tf, MOBILE = task_factory()
n_shell = len(MOBILE)
mm = movemap()
print(f"  packer: {len(DESIGN)} designable + {n_shell - len(DESIGN)} repackable "
      f"within {args.pack_shell:.0f} A of the ligand; "
      f"{start.size() - n_shell - 1} residues frozen; "
      f"FastRelax rounds={args.rounds}; "
      f"favor_native={args.favor_native:g}", flush=True)

rows = []
t_all = time.time()
for k in range(args.n_per_unit):
    traj = block * args.n_per_unit + k
    seed = 1000 + 17 * traj
    # reseeds every generator, not just the default one; verified reproducible
    pyrosetta.rosetta.basic.random.init_random_generators(seed, "mt19937")
    pose = start.clone()

    t0 = time.time()
    fd = FastRelax(sfxn_cst, args.rounds)
    fd.set_task_factory(tf)
    fd.set_movemap(mm)
    fd.apply(pose)
    dt = time.time() - t0

    total, dG = interface_energy(pose)
    seq = {p: pose.residue(pose_idx[p]).name1() for p in DESIGN}
    muts = {p: seq[p] for p in DESIGN
            if seq[p] != start.residue(pose_idx[p]).name1()}
    rec = dict(arm=tag, struct=struct, alphabet=alphabet, traj=traj, seed=seed,
               favor_native=args.favor_native,
               total=round(total, 3), interface_dG=round(dG, 3),
               seconds=round(dt, 1),
               seq="".join(seq[p] for p in DESIGN),
               mutations={str(p): a for p, a in muts.items()},
               recovered=[f"{w}{p}{t}" for p, (w, t) in TRUE.items()
                          if seq.get(p) == t])
    rows.append(rec)
    if args.dump_pdb:
        pose.dump_pdb(os.path.join(args.outdir, "pdb", f"{tag}_t{traj}.pdb"))
    print(f"  traj {traj:>3}  {dt/60:5.1f} min  dG={dG:8.2f}  "
          f"n_mut={len(muts):>2}  "
          f"recovered={','.join(rec['recovered']) or '-'}", flush=True)
    # write incrementally: a walltime kill should not discard finished work
    json.dump(dict(arm=tag, struct=struct, alphabet=alphabet, block=block,
                   favor_native=args.favor_native,
                   ligand=comp, params=params_file, n_lig_heavy=n_lig_heavy,
                   pack_shell=args.pack_shell, rounds=args.rounds,
                   n_shell=n_shell, design_positions=DESIGN, trajectories=rows),
              open(os.path.join(args.outdir, f"{tag}__b{block}.json"), "w"),
              indent=2)

print(f"wrote {tag}__b{block}.json  "
      f"({len(rows)} trajectories, {(time.time()-t_all)/60:.1f} min total)")
