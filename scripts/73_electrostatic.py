#!/usr/bin/env python
"""
73_electrostatic.py -- electrostatic and hydrogen-bond complementarity at
NON-CLASHING pocket positions: the stage §32d showed is missing.

THE GAP THIS FILLS
------------------
§32's pairwise steric enumeration recovered F108A (rank 5), F159L (75) and V81I
(150) of 13,457 pairs, and missed K59R at rank 3907. That miss is a boundary of
CATEGORY, not of implementation: K59 does not clash with the ligand (2.86 A away),
its contribution is a salt bridge, and a steric method cannot see electrostatics.
The same residue, for the same reason, defeated LigandMPNN and Rosetta FastDesign
in §23j.

THE DECISIVE TEST, AND WHY IT IS A GOOD ONE
-------------------------------------------
The two ligands differ in exactly the way that matters:

    ABA            net -0.97   carboxylate O3/O4 at -0.74 each, plus O1/O2
    mandipropamid  net +0.10   amide N1/O2, two ethers, NO formal charge

ABA needs a counter-ion at 59 and WT already supplies it (K59). Mandipropamid needs
a hydrogen-bond donor, and the ground truth is K59R. So one method, run on two
ligands, must output **K at 59 for ABA and R at 59 for mandipropamid**. A single
flip, with the WT-retention control and the ligand-discrimination control fused into
one number. Stage 1's null arm failed exactly here -- it kept WT K59 in 0% of
trajectories even with the cognate ligand.

WHAT IS DELIBERATELY NOT MODELLED
---------------------------------
**Desolvation.** ref2015 does model it and gets the balance wrong here by a
quantified +10.1 REU (§23j), which is why K59R was unreachable. This script does
not attempt the balance at all: it scores COMPLEMENTARITY -- is there a partner for
each polar group, and do the charges match -- and is therefore a screen, not an
affinity predictor. That has a predictable failure mode of its own, in the opposite
direction: with no desolvation penalty, everything wants to be charged. So the
Arg-everywhere check below is not decoration, it is the control that decides whether
the score means anything.

CONTROLS
--------
  1. WT RETENTION -- with the cognate ligand, is WT still preferred where it should be?
  2. LIGAND DISCRIMINATION -- do the two ligands give different answers at all?
  3. LIGAND-BLIND BASELINE -- residues scored by their own polarity with no ligand.
     If the ranking matches that, the method is saying "polar residues are polar".
  4. ARG-EVERYWHERE -- what fraction of positions pick R? A score dominated by
     charge magnitude rather than complementarity gives itself away here.

Run with the tier1_analysis env python (PyRosetta).
"""
import json
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
STAGE1 = os.path.join(ROOT, "data", "stage1")
OUT = os.path.join(ROOT, "data", "electrostatic")
os.makedirs(OUT, exist_ok=True)

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "F": 1.47, "CL": 1.75,
       "BR": 1.85, "I": 1.98, "P": 1.80}
AA20 = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
        "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
BB = ("N", "CA", "C", "O")
TOL = 0.75            # steric feasibility, as in §32
HB_MIN, HB_MAX = 2.5, 3.5      # heavy-atom donor-acceptor separation
HB_ANGLE = 110.0               # minimum D-H...A angle, degrees
POCKET_CUTOFF = 6.0
N_ROT = 8
#: ⚠ USE THE ROSETTA-WRITTEN INPUTS, not data/stage1/wt_*.pdb.
#:
#: Rosetta matches a params file to a PDB residue by name, and then maps atoms by
#: NAME. The two things that went wrong here, in order:
#:   1. both stage-1 params call themselves `LIG` while data/stage1/wt_aba.pdb and
#:      wt_mandi.pdb contain `A8S` and `3UZ`, so neither ever matched. Rosetta
#:      silently built the ligand from its own PDB-components dictionary instead --
#:      `pdb_A8S`, NEUTRAL, charge sum 0.000 where the params sum to -0.970 -- and
#:      reported no error whatsoever.
#:   2. renaming the params did not fix it, because the ATOM NAMES also differ:
#:      those PDBs carry RCSB convention (`CAC CAB CAA OAI ...`) while the params
#:      carry generation order (`C12 N1 C11 ...`). They were never a matched pair.
#:
#: The matched pair is stage 1's own Rosetta-written input, where the ligand is
#: `LIG` with params-consistent atom names -- 38 atoms for ABA, 51 for 3UZ.
#:
#: §32 is UNAFFECTED: it used only coordinates and element-derived radii, never a
#: charge, and the coordinates are identical either way.
COMPLEXES = {
    "ABA":   (os.path.join(ROOT, "results", "stage1_rosetta",
                           "_input_wt_aba_A8S_anion.pdb"),
              os.path.join(STAGE1, "params", "A8S_anion.params"), "LIG", -0.970),
    "mandi": (os.path.join(ROOT, "results", "stage1_rosetta",
                           "_input_wt_mandi_3UZ.pdb"),
              os.path.join(STAGE1, "params", "3UZ.params"), "LIG", 0.100),
}


def log(m):
    print(m, flush=True)


def residue_atoms(r, side_chain_only=False, heavy_only=True):
    """[(name, element, xyz, charge, is_donor, is_acceptor, attached_H_xyz)]

    ⚠ `heavy_only` must be FALSE for anything electrostatic. Dropping hydrogens
    does not merely lose a little charge, it loses the POLAR hydrogens that carry
    most of the hydrogen-bond electrostatics, and it makes the ligand net charge
    nonsense: ABA summed to -3.20 instead of -0.97 and mandipropamid to -2.59
    instead of +0.10. Sterics use heavy atoms; charges use all of them.
    """
    out = []
    for a in range(1, r.natoms() + 1):
        t = r.atom_type(a)
        e = t.element().strip().upper()
        nm = r.atom_name(a).strip()
        if heavy_only and e == "H":
            continue
        if side_chain_only and nm in BB:
            continue
        v = r.xyz(a)
        hs = []
        for h in range(1, r.natoms() + 1) if e != "H" else []:
            if r.atom_type(h).element().strip().upper() != "H":
                continue
            hv = r.xyz(h)
            if (np.linalg.norm(np.array([hv.x, hv.y, hv.z]) -
                               np.array([v.x, v.y, v.z]))) < 1.3:
                hs.append(np.array([hv.x, hv.y, hv.z]))
        out.append(dict(name=nm, elem=e, xyz=np.array([v.x, v.y, v.z]),
                        q=float(r.atomic_charge(a)),
                        don=bool(t.is_donor()), acc=bool(t.is_acceptor()), hs=hs))
    return out


def hbonds(res_atoms, lig_atoms):
    """geometrically valid donor-acceptor pairs, with a D-H...A angle check"""
    n = 0
    detail = []
    for A in res_atoms:
        for B in lig_atoms:
            for don, acc in ((A, B), (B, A)):
                if not (don["don"] and acc["acc"]):
                    continue
                d = float(np.linalg.norm(don["xyz"] - acc["xyz"]))
                if not (HB_MIN <= d <= HB_MAX):
                    continue
                ok = not don["hs"]          # no explicit H -> accept on distance
                for h in don["hs"]:
                    v1, v2 = don["xyz"] - h, acc["xyz"] - h
                    c = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9)
                    if np.degrees(np.arccos(np.clip(c, -1, 1))) >= HB_ANGLE:
                        ok = True
                if ok:
                    n += 1
                    detail.append((A["name"], B["name"], round(d, 2)))
                    break
    return n, detail


def coulomb(res_atoms, lig_atoms):
    """
    Screened Coulomb complementarity, distance-dependent dielectric eps = 4r,
    so E ~ q1 q2 / (4 r^2). Negative is favourable.

    Reported in arbitrary units and NEVER added to a steric term: this is a
    complementarity ranking, not an energy, and mixing the two would recreate
    exactly the balance ref2015 gets wrong here.
    """
    if not res_atoms or not lig_atoms:
        return 0.0
    P = np.array([a["xyz"] for a in res_atoms]); QP = np.array([a["q"] for a in res_atoms])
    L = np.array([a["xyz"] for a in lig_atoms]); QL = np.array([a["q"] for a in lig_atoms])
    d = np.linalg.norm(P[:, None, :] - L[None, :, :], axis=2)
    d = np.clip(d, 1.5, None)
    return float((QP[:, None] * QL[None, :] / (4.0 * d ** 2)).sum() * 332.0)


def max_overlap(x1, r1, x2, r2):
    if len(x1) == 0 or len(x2) == 0:
        return 0.0
    d = np.linalg.norm(x1[:, None, :] - x2[None, :, :], axis=2)
    return float(((r1[:, None] + r2[None, :]) - d).max())


def analyse(tag, pdb, params, code, want_q, background=None):
    """`background` = {resnum: AA3} applied before scoring, to test whether a
    substitution only becomes reachable once its partner mutation exists."""
    import pyrosetta
    from pyrosetta.rosetta.core.chemical import ChemicalManager
    from pyrosetta.rosetta.core.pack.rotamer_set import bb_independent_rotamers
    # -load_PDB_components false: without it Rosetta prefers its OWN PDB-components
    # entry over the supplied params whenever the ligand code exists there, which is
    # how 3UZ kept loading as a neutral `pdb_3UZ` even after the NAME was corrected.
    pyrosetta.init(f"-mute all -extra_res_fa {params} -detect_disulf false "
                   f"-load_PDB_components false "
                   f"-run:constant_seed -run:jran 20260818")
    pose = pyrosetta.pose_from_pdb(pdb)
    pi = pose.pdb_info()
    rts = ChemicalManager.get_instance().residue_type_set("fa_standard")

    if background:
        # Apply the background mutation(s) first. This is the whole point of the
        # experiment below: K59R makes no hydrogen bond in the WT pocket, so the
        # question is whether it becomes reachable once F108A has opened the space.
        from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
        pinf = pose.pdb_info()
        idx = {pinf.number(i): i for i in range(1, pose.total_residue() + 1)}
        for num, aa in background.items():
            MutateResidue(idx[num], aa).apply(pose)
        log(f"  background applied: {background}")
    lig_res = None
    prot = {}
    for i in range(1, pose.total_residue() + 1):
        r = pose.residue(i)
        if r.is_protein():
            prot[pi.number(i)] = i
        else:
            lig_res = r
    # ASSERT the params were actually used. A name mismatch is silent, and the only
    # visible symptom is a ligand that has quietly become neutral.
    assert lig_res is not None, f"{tag}: no ligand in the pose"
    assert not lig_res.name().startswith("pdb_"), (
        f"{tag}: ligand loaded as '{lig_res.name()}' -- that prefix means Rosetta "
        f"auto-built it from its PDB-components dictionary and IGNORED {params}. "
        f"The params NAME field must equal the PDB residue name ({code}).")
    qsum = sum(lig_res.atomic_charge(a) for a in range(1, lig_res.natoms() + 1))
    assert abs(qsum - want_q) < 0.05, (
        f"{tag}: ligand net charge {qsum:+.3f}, expected {want_q:+.3f} from {params}")
    # same protein, same numbering as everywhere else in this project
    # Positions carrying a deliberate background mutation are excluded -- the check
    # is that this is the right WT protein, not that nothing was mutated on purpose.
    EXPECT = {59: "LYS", 81: "VAL", 108: "PHE", 115: "HIS", 120: "TYR",
              159: "PHE", 160: "ALA"}
    for num, want in EXPECT.items():
        if background and num in background:
            continue
        got = pose.residue(prot[num]).name3().strip()
        assert got == want, f"{tag}: residue {num} is {got}, expected {want}"
    lig = residue_atoms(lig_res)                              # heavy, for sterics/H-bonds
    lig_all = residue_atoms(lig_res, heavy_only=False)        # all, for charges
    LX = np.array([a["xyz"] for a in lig])
    LR = np.array([VDW.get(a["elem"], 1.7) for a in lig])

    # environment for steric feasibility
    env_x, env_r, owner = [], [], []
    for n, i in prot.items():
        for a in residue_atoms(pose.residue(i)):
            env_x.append(a["xyz"]); env_r.append(VDW.get(a["elem"], 1.7)); owner.append(n)
    EX, ER, OWN = np.array(env_x), np.array(env_r), np.array(owner)

    # pocket + which positions actually CLASH with the ligand in WT
    pocket, clashing = [], []
    for n, i in prot.items():
        A = residue_atoms(pose.residue(i))
        X = np.array([a["xyz"] for a in A]); R = np.array([VDW.get(a["elem"], 1.7) for a in A])
        d = np.linalg.norm(X[:, None, :] - LX[None, :, :], axis=2)
        if d.min() < POCKET_CUTOFF:
            pocket.append(n)
            if ((R[:, None] + LR[None, :]) - d).max() > 0.4:
                clashing.append(n)
    pocket.sort()
    log(f"\n=== {tag} === ligand net charge {sum(a['q'] for a in lig_all):+.2f}, "
        f"{sum(1 for a in lig if a['don'] or a['acc'])} polar atoms")
    log(f"  {len(pocket)} pocket positions; CLASHING in WT: {sorted(clashing)}")

    res = defaultdict(dict)
    for n in pocket:
        i = prot[n]
        target = pose.residue(i)
        keep = ~np.isin(OWN, [n, ])
        for aa in AA20:
            try:
                rots = bb_independent_rotamers(rts.name_map(aa))
            except Exception:
                continue
            best = None
            for k in range(1, min(len(rots), 60) + 1):
                r = rots[k].clone()
                r.orient_onto_residue(target)
                A = residue_atoms(r, side_chain_only=True)
                if not A:
                    A = []
                X = np.array([a["xyz"] for a in A]).reshape(-1, 3)
                R = np.array([VDW.get(a["elem"], 1.7) for a in A])
                if len(X) and max_overlap(X, R, EX[keep], ER[keep]) > TOL:
                    continue
                if len(X) and max_overlap(X, R, LX, LR) > TOL:
                    continue
                nhb, det = hbonds(A, lig)
                cl = coulomb(residue_atoms(r, side_chain_only=True, heavy_only=False),
                             lig_all)
                # rank rotamers by H-bonds first, then electrostatics
                key = (-nhb, cl)
                if best is None or key < best[0]:
                    best = (key, dict(hb=nhb, coul=cl, detail=det))
            if best:
                res[n][aa] = best[1]
    return dict(pose=pose, prot=prot, pocket=pocket, clashing=clashing, res=res,
                wt={n: pose.residue(prot[n]).name3().strip() for n in pocket},
                lig_charge=float(sum(a["q"] for a in lig_all)))


def run_one(tag, background=None):
    c = analyse(tag, *COMPLEXES[tag], background=background)
    out = dict(pocket=c["pocket"], clashing=c["clashing"], wt=c["wt"],
               lig_charge=c["lig_charge"],
               res={str(n): {aa: dict(hb=v["hb"], coul=v["coul"], detail=v["detail"])
                             for aa, v in d.items()}
                    for n, d in c["res"].items()})
    suffix = "" if not background else "_bg" + "".join(
        f"{k}{v}" for k, v in sorted(background.items()))
    path = os.path.join(OUT, f"electrostatic_{tag}{suffix}.json")
    json.dump(out, open(path, "w"), indent=0)
    log(f"  wrote {path}")


def main():
    """
    ⚠ ONE LIGAND PER PROCESS, and this is not an optimisation.

    Both params files declare `NAME LIG`. PyRosetta's residue type set is global and
    a second `init` does not replace an already-registered type, so loading ABA and
    then mandipropamid in one process silently gives mandipropamid ABA's charges
    (-0.970 instead of +0.100) -- caught only because the net charge is asserted.
    Two different ligands cannot share a name inside one process, so each gets its
    own.
    """
    import subprocess
    import sys
    if len(sys.argv) > 1:
        bg = json.loads(sys.argv[2]) if len(sys.argv) > 2 else None
        run_one(sys.argv[1], {int(k): v for k, v in bg.items()} if bg else None)
        return
    for tag in COMPLEXES:
        log(f"=== subprocess for {tag} ===")
        r = subprocess.run([sys.executable, os.path.abspath(__file__), tag])
        if r.returncode != 0:
            raise SystemExit(f"{tag} failed")
    merged = {t: json.load(open(os.path.join(OUT, f"electrostatic_{t}.json")))
              for t in COMPLEXES}
    json.dump(merged, open(os.path.join(OUT, "electrostatic.json"), "w"), indent=0)
    log(f"\nwrote {OUT}/electrostatic.json")


if __name__ == "__main__":
    main()
