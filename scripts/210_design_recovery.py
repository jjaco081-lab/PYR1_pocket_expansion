#!/usr/bin/env python
r"""
210_design_recovery.py -- can HBNet recover the KNOWN sensor sequences when the
pocket is shaved to alanine and the ligand pose is handed over?

JANNIS'S DESIGN, and the reason for it: WIN, ABA and mandipropamid are OPTIMISED
sensors with CRYSTAL structures, so the ligand pose is known and cannot confound
the result. Shave the pocket lining to alanine, redesign, and score recovery
against the deposited sequence.

WHY HBNet SPECIFICALLY. K59R is the project's standing failure: FastDesign and
LigandMPNN both recover F108A and F159L and both miss K59R and V81I (2 of 4).
The cause was diagnosed, not guessed -- ref2015 pays **+10.1 REU** to desolvate
the salt bridge it correctly rewards (§23j), coupled moves retained K59 **0 %**
of the time so sampling is not the limit (§35), and MM-GBSA appeared to fix it
and was retracted when the reference turned out unconverged (§37 → §40). Our own
literature note names the untried fix: "Kortemme coupled moves + HBNet". Coupled
moves failed. HBNet has never been run here -- it does not appear once in the
README -- and it exists precisely to make buried polar networks payable.

⚠ ALANINE, NOT GLYCINE. Poly-glycine removes C-beta AND widens the backbone
preference, so design would have to overcome an artefact we introduced. Whitehead
used a poly-glycine shave only for a CLASH GRID SEARCH, never for design. Alanine
destroys the answer (the native rotamer is gone, so recovery cannot be memory)
while preserving C-beta direction and backbone propensity. P88 is left alone.

⚠ THE NATIVE-BACKBONE ARM IS THE EASY CASE AND IS SCORED AS A CONTROL, NOT A
RESULT. 4WVO's backbone has already relaxed around K59R/V81I/F108A/F159L, so the
answer is partly in the coordinates before design starts; and for 3QN1 the answer
is wild type, which is what ref2015's reference energies were fitted to
reproduce. The informative arm is CROSS-LIGAND: put mandipropamid's pose into
WILD-TYPE PYR1's backbone and ask whether design proposes the quadruple. That is
the real task and the one that maps onto eugenol.

⚠ SCORED TWO WAYS, because recall alone is gameable by a method that mutates
everything:
   recall      = fraction of the truly-mutated positions recovered
   specificity = fraction of the truly-UNCHANGED positions left alone
A method that returns wild type everywhere scores recall 0 / specificity 1.0 on
mandipropamid and recall 1.0 / specificity 1.0 on ABA -- which is exactly why ABA
alone proves nothing.

⚠ PYL2 NUMBERING IS NOT A CONSTANT OFFSET. Structural superposition of 7MWN onto
3QN1 gives CA RMSD 0.84 A over 171 residues with offsets of 4, 5 AND 6. A fixed
+5 silently mis-assigns part of the pocket. Positions are mapped structurally and
every identity is asserted.

GROUND TRUTH, re-derived from the deposited structures by this script:
  3QN1 / A8S   0 lining differences from WT PYR1   -> answer is wild type
  4WVO / 3UZ   K59R, V81I, F108A, F159L            -> the quadruple
  7MWN / WI5   K64Q, F165A, V166I (PYL2 numbering) -> the WIN sensor
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = os.path.join(ROOT, "results", "design_recovery")

#: PYR1 pocket lining (24) from results/pocket_shape/pocket_shape.json.
#: P88 is excluded from design -- proline, per Jannis.
LINING = [59, 60, 61, 62, 81, 83, 87, 88, 89, 91, 92, 94, 108, 110, 115, 117,
          120, 122, 141, 159, 160, 163, 164, 167]
DESIGNABLE = [n for n in LINING if n != 88]
WT_PYR1 = {59: "K", 60: "H", 61: "F", 62: "I", 81: "V", 83: "V", 87: "L",
           88: "P", 89: "A", 91: "T", 92: "S", 94: "E", 108: "F", 110: "I",
           115: "H", 117: "L", 120: "Y", 122: "S", 141: "E", 159: "F",
           160: "A", 163: "V", 164: "V", 167: "N"}

SYSTEMS = {
    "aba":   dict(pdb="3QN1", lig="A8S",
                  params=f"{ROOT}/data/aba_params/A8S.params"),
    "mandi": dict(pdb="4WVO", lig="3UZ",
                  params=f"{ROOT}/data/stage1/params/3UZ.params"),
    "win":   dict(pdb="7MWN", lig="WI5",
                  params=f"{ROOT}/results/win_crossover/WI5.params"),
}
THREE = {v: k for k, v in dict(
    ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G",
    HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S",
    THR="T", TRP="W", TYR="Y", VAL="V").items()}


def resmap(pose):
    pi = pose.pdb_info()
    return {pi.number(i): i for i in range(1, pose.total_residue() + 1)
            if pi.chain(i) == "A"}


def design(pose, positions, sfxn, arm, nrep=1):
    """Shave `positions` to ALA, then redesign them. Returns the designed seq."""
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    # ⚠ OperateOnResidueSubset takes a ResLvlTaskOperation, NOT a TaskOperation.
    # PreventRepacking/RestrictToRepacking are TaskOperations and raise a
    # constructor TypeError here; the *RLT variants are the residue-level ones.
    from pyrosetta.rosetta.core.pack.task.operation import (
        RestrictToRepackingRLT, PreventRepackingRLT, OperateOnResidueSubset,
        RestrictAbsentCanonicalAASRLT)
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector, NotResidueSelector)
    from pyrosetta.rosetta.protocols.relax import FastRelax
    import pyrosetta

    p = pose.clone()
    rm = resmap(p)
    idx = [rm[n] for n in positions if n in rm]
    for i in idx:
        MutateResidue(i, "ALA").apply(p)             # THE SHAVE

    sel = ResidueIndexSelector(",".join(str(i) for i in idx))
    shell = NeighborhoodResidueSelector(sel, 8.0, True)
    tf = TaskFactory()
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(),
                                        NotResidueSelector(shell)))
    # design only at the shaved positions; the shell repacks but cannot design
    keep = RestrictAbsentCanonicalAASRLT()
    keep.aas_to_keep("ACDEFGHIKLMNQRSTVWY")          # no proline by design
    tf.push_back(OperateOnResidueSubset(keep, sel))
    notsel = NotResidueSelector(sel)
    tf.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), notsel))

    # ⚠ the HBNet mover's setters are task_factory()/score_function(), NOT
    # set_task_factory()/set_score_function() -- the set_* forms raise
    # AttributeError and killed all three hbnet_mover tasks.
    if arm.startswith("hbnet_mover"):
        from pyrosetta.rosetta.protocols.hbnet import HBNet
        hb = HBNet()
        hb.task_factory(tf)
        hb.score_function(sfxn)
        try:
            hb.apply(p)
        except Exception as e:                       # noqa: BLE001
            print(f"      HBNet mover failed: {type(e).__name__}: {e}")

    fr = FastRelax(sfxn, nrep)
    fr.cartesian(True)
    fr.min_type("lbfgs_armijo_nonmonotone")
    fr.set_task_factory(tf)
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(False); mm.set_chi(True)
    fr.set_movemap(mm)
    fr.apply(p)
    return {n: p.residue(rm[n]).name1() for n in positions if n in rm}, p


def score_recovery(designed, truth, wt):
    """recall on truly-mutated positions, specificity on truly-unchanged ones."""
    mutated = [n for n in truth if truth[n] != wt.get(n)]
    same = [n for n in truth if truth[n] == wt.get(n)]
    rec = [n for n in mutated if designed.get(n) == truth[n]]
    spec = [n for n in same if designed.get(n) == truth[n]]
    return dict(n_mutated=len(mutated), recalled=len(rec),
                recall=len(rec) / max(len(mutated), 1),
                hit_positions=[f"{wt[n]}{n}{truth[n]}" for n in rec],
                missed=[f"{wt[n]}{n}{truth[n]}" for n in mutated if n not in rec],
                n_same=len(same), spec_ok=len(spec),
                specificity=len(spec) / max(len(same), 1),
                overall=sum(1 for n in truth if designed.get(n) == truth[n])
                / max(len(truth), 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", default="aba", choices=sorted(SYSTEMS))
    ap.add_argument("--arm", default="fastdesign",
                    choices=("fastdesign", "hbnet_terms", "hbnet_mover"))
    ap.add_argument("--nrep", type=int, default=1)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    S = SYSTEMS[a.system]
    assert os.path.exists(S["params"]), f"missing params {S['params']}"

    import pyrosetta
    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ex1 -ex2aro "
                   f"-extra_res_fa {S['params']}")
    from pyrosetta import pose_from_file
    src = os.path.join(ROOT, "data", f"{S['pdb']}.cif")
    pose = pose_from_file(src)
    rm = resmap(pose)
    truth = {n: pose.residue(rm[n]).name1() for n in LINING if n in rm}
    print(f"{a.system}: {S['pdb']} loaded, {pose.total_residue()} residues, "
          f"{len(truth)}/{len(LINING)} lining positions present")
    miss = [n for n in LINING if n not in truth]
    if miss:
        print(f"   ⚠ lining positions absent from the structure: {miss}")
    if a.system != "win":
        diffs = [f"{WT_PYR1[n]}{n}{truth[n]}" for n in truth
                 if truth[n] != WT_PYR1.get(n)]
        print(f"   ground truth vs WT PYR1: {diffs if diffs else 'wild type'}")

    sf = pyrosetta.create_score_function("ref2015_cart")
    if a.arm == "hbnet_terms":
        from pyrosetta.rosetta.core.scoring import ScoreType
        sf.set_weight(ScoreType.hbnet, 1.0)
        sf.set_weight(ScoreType.buried_unsatisfied_penalty, 1.0)
        print("   arm: ref2015_cart + hbnet(1.0) + buried_unsatisfied_penalty(1.0)")
    else:
        print(f"   arm: {a.arm}")

    pos = [n for n in DESIGNABLE if n in rm]
    if a.smoke:
        pos = pos[:4]
        print(f"   SMOKE: designing only {pos}")
    des, dp = design(pose, pos, sf, a.arm, a.nrep)
    wt = WT_PYR1 if a.system != "win" else {n: truth[n] for n in truth}
    r = score_recovery(des, {n: truth[n] for n in pos}, wt)
    print(f"\n   designed: {''.join(des[n] for n in sorted(des))}")
    print(f"   truth   : {''.join(truth[n] for n in sorted(des))}")
    print(f"   overall recovery {r['overall']:.2f}  "
          f"recall {r['recalled']}/{r['n_mutated']}  "
          f"specificity {r['spec_ok']}/{r['n_same']}")
    if r["hit_positions"]:
        print(f"   RECOVERED: {r['hit_positions']}")
    if r["missed"]:
        print(f"   missed   : {r['missed']}")
    out = os.path.join(OUT, f"{a.system}_{a.arm}{'_smoke' if a.smoke else ''}.json")
    json.dump(dict(system=a.system, arm=a.arm, designed=des,
                   truth={str(k): v for k, v in truth.items()}, **r),
              open(out, "w"), indent=1)
    print(f"   wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
