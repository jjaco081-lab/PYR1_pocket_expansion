#!/usr/bin/env python
"""
74_coupled_moves.py -- stage 1 again, with the method the literature says fixes it.

WHY
---
README 23j: stage 1 does not clear. LigandMPNN and Rosetta FastDesign both recover
only the severe-clash positions, and both miss K59R and V81I -- the two where the
ligand barely touches the residue.

README 34 turned up the direct precedent. Ollikainen, de Jong & Kortemme (PLOS
Comput Biol 2015) built a benchmark of exactly this task -- predicting
specificity-altering mutations -- showed that FIXED-BACKBONE design performs poorly
on it, and introduced **coupled moves**, which samples sequence, side chains,
backbone AND the ligand's pose together. It gave a **5.75x** increase in correct
predictions and won **16 of 17** benchmark cases. Our stage 1 was fixed-backbone.

So this runs the same arms, the same designable positions and the same scoring rule
as §23j, changing only the sampler.

PRE-REGISTERED READING, fixed before the run
--------------------------------------------
The two diagnoses this project has made point in different directions, and the
outcome distinguishes them:

  * §23j says the K59 failure is SCORING -- ref2015 pays +10.1 REU to bury the
    ammonium, and coupled moves uses the same score function. On that account
    K59R stays missed however well the backbone is sampled.
  * §33b says the miss at the electrostatic stage was SAMPLING -- the crystal Arg59
    is non-rotameric and was never proposed. On that account backbone and ligand
    flexibility should reach it.

**K59R recovered => sampling was the limit. K59R still missed => the scoring
diagnosis stands and no sampler will fix it.** Either is a result. What is NOT
allowed is to re-tune afterwards until it agrees.

SCORING RULE, unchanged from §23
--------------------------------
Only `f_mandi(true) - f_aba(true)` is evidence. The null is a LIGAND SWAP, not a
random draw, so it inherits the method's own biases and needs no independence claim.
Frequencies, never recall-at-N: uniform chance at N=50 is 3.84 of 4.

⚠ ONE LIGAND PER PROCESS. Both params files declare `NAME LIG` and PyRosetta's
residue type set is global, so loading both in one process silently gives the second
ligand the first one's chemistry (README 33e).

Run with the tier1_analysis env python.
"""
import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "data", "coupled_moves")
os.makedirs(OUT, exist_ok=True)

#: The stage-1 designable set, unchanged: Tian's 18 pocket positions minus the gate
#: (87, 89) and the latch (117), which carry the transduction mechanism. No
#: ground-truth mutation sits at the excluded three.
DESIGN = [59, 81, 83, 92, 94, 108, 110, 120, 122, 141, 159, 160, 163, 164, 167]
GROUND_TRUTH = {59: "R", 81: "I", 108: "A", 159: "L"}     # 4WVO vs 3QN1

ARMS = {
    "mandi": (os.path.join(ROOT, "results", "stage1_rosetta", "_input_wt_mandi_3UZ.pdb"),
              os.path.join(ROOT, "data", "stage1", "params", "3UZ.params")),
    "aba":   (os.path.join(ROOT, "results", "stage1_rosetta", "_input_wt_aba_A8S_anion.pdb"),
              os.path.join(ROOT, "data", "stage1", "params", "A8S_anion.params")),
}

NTRIALS = int(os.environ.get("CM_NTRIALS", 1000))        # Rosetta default; the per-trajectory Monte Carlo length
LIGAND_WEIGHT = 1.0   # Rosetta default, recorded rather than tuned
BACKBONE_MOVER = "backrub"
ONE = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
       "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
       "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
       "TRP": "W", "TYR": "Y"}


def log(m):
    print(m, flush=True)


def build_task(pose, design_idx):
    """Design at DESIGN, repack their neighbours, freeze the rest. Never design the ligand."""
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task.operation import (
        RestrictToRepacking, OperateOnResidueSubset, PreventRepackingRLT,
        RestrictToRepackingRLT, IncludeCurrent)
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector, NotResidueSelector,
        OrResidueSelector)
    des = ResidueIndexSelector(",".join(map(str, sorted(design_idx))))
    shell = NeighborhoodResidueSelector(des, 6.0, False)
    tf = TaskFactory()
    tf.push_back(IncludeCurrent())
    # everything outside design+shell is frozen; the shell repacks but cannot design
    tf.push_back(OperateOnResidueSubset(
        PreventRepackingRLT(), OrResidueSelector(des, shell), True))
    tf.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), shell))
    return tf


def run_arm(arm, rep):
    import pyrosetta
    pdb, params = ARMS[arm]
    seed = 20260819 + rep
    pyrosetta.init(f"-mute all -extra_res_fa {params} -detect_disulf false "
                   f"-load_PDB_components false -ex1 -ex2 "
                   f"-run:constant_seed -run:jran {seed}")
    from pyrosetta.rosetta.protocols.coupled_moves import CoupledMovesProtocol

    pose = pyrosetta.pose_from_pdb(pdb)
    pi = pose.pdb_info()
    idx = {pi.number(i): i for i in range(1, pose.total_residue() + 1)
           if pose.residue(i).is_protein()}
    lig = [i for i in range(1, pose.total_residue() + 1) if not pose.residue(i).is_protein()]
    assert len(lig) == 1, f"{arm}: expected one ligand, found {len(lig)}"
    # identity, not position -- the standing rule in this project
    for num, aa in ((59, "LYS"), (81, "VAL"), (108, "PHE"), (159, "PHE")):
        got = pose.residue(idx[num]).name3().strip()
        assert got == aa, f"{arm}: residue {num} is {got}, expected {aa}"

    design_idx = [idx[n] for n in DESIGN]
    tf = build_task(pose, design_idx)
    sf = pyrosetta.create_score_function("ref2015")

    cm = CoupledMovesProtocol()
    cm.set_main_task_factory(tf)
    cm.set_score_fxn(sf)
    cm.set_ntrials(NTRIALS)
    cm.set_ligand_mode(True)          # the ligand moves with the sequence -- the point
    cm.set_number_ligands(1)
    cm.set_ligand_weight(LIGAND_WEIGHT)
    cm.set_backbone_mover(BACKBONE_MOVER)

    work = pose.clone()
    cm.apply(work)
    seq = {n: ONE.get(work.residue(idx[n]).name3().strip(), "X") for n in DESIGN}
    uniq = {}
    try:
        for s, v in cm.get_unique_sequences().items():
            uniq[str(s)] = float(v)
    except Exception:
        pass
    return dict(arm=arm, rep=rep, seed=seed, seq=seq, n_unique=len(uniq))


def main():
    if len(sys.argv) > 2:
        arm, rep = sys.argv[1], int(sys.argv[2])
        r = run_arm(arm, rep)
        path = os.path.join(OUT, f"{arm}_rep{rep:03d}.json")
        json.dump(r, open(path, "w"))
        log(f"CM_DONE arm={arm} rep={rep} unique={r['n_unique']} -> {path}")
        return
    # aggregate
    freq = {a: defaultdict(Counter) for a in ARMS}
    n = {a: 0 for a in ARMS}
    for f in sorted(os.listdir(OUT)):
        if not f.endswith(".json") or f == "summary.json":
            continue
        d = json.load(open(os.path.join(OUT, f)))
        n[d["arm"]] += 1
        for k, v in d["seq"].items():
            freq[d["arm"]][int(k)][v] += 1
    log(f"trajectories: " + ", ".join(f"{a}={n[a]}" for a in ARMS))
    if not all(n.values()):
        log("  incomplete -- run the array first"); return
    log("\n  ground truth, by the §23 ligand-swap rule (only the delta is evidence):")
    log(f"  {'mut':<7}{'f(mandi)':>10}{'f(aba)':>9}{'delta':>8}   verdict")
    res = {}
    for pos, aa in GROUND_TRUTH.items():
        fm = freq["mandi"][pos][aa] / n["mandi"]
        fa = freq["aba"][pos][aa] / n["aba"]
        v = "recovered" if fm - fa > 0.2 else ("inverted" if fm - fa < -0.2 else "missed")
        res[f"{pos}{aa}"] = dict(mandi=fm, aba=fa, delta=fm - fa, verdict=v)
        log(f"  {pos}{aa:<6}{fm:>10.2f}{fa:>9.2f}{fm-fa:>8.2f}   {v}")
    log("\n  WT retention in the ABA arm (the control stage 1 failed at 29%):")
    wt = {59: "K", 81: "V", 83: "V", 92: "S", 94: "E", 108: "F", 110: "I", 120: "Y",
          122: "S", 141: "E", 159: "F", 160: "A", 163: "V", 164: "V", 167: "N"}
    keep = sum(freq["aba"][p][wt[p]] for p in DESIGN) / (n["aba"] * len(DESIGN))
    log(f"    {keep:.0%} of designable positions keep WT with the cognate ligand")
    log(f"    K59 retained in {freq['aba'][59]['K']/n['aba']:.0%} of trajectories "
        f"(stage 1: 0%, and no favor-native weight fixed it)")
    json.dump(dict(n=n, ground_truth=res, wt_retention=keep,
                   k59_retained=freq["aba"][59]["K"] / n["aba"],
                   freq={a: {str(p): dict(c) for p, c in freq[a].items()} for a in ARMS}),
              open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    log(f"\nwrote {OUT}/summary.json")


if __name__ == "__main__":
    main()
