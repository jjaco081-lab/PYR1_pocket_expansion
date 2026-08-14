#!/usr/bin/env python
"""
63_dock_noncognate.py -- dock three ABA-sized, non-cognate ligands into the CLOSED
PYR1 pocket, for the §27 ligand-dependence control.

WHY THESE THREE
---------------
§27 asks whether the loop-dynamics signal is ligand-dependent at all, or whether
closed PYR1 holds its state around anything that fits. One ligand cannot answer
that -- a single negative could be a property of that molecule. Three from
different chemical classes, matched in size to ABA, can.

    ABA (reference)     19 heavy atoms, MW 260.3, sesquiterpenoid acid
    Imperatorin         20 heavy, MW 270.3, furanocoumarin      (cluster 26)
    Flutamide           19 heavy, MW 276.2, nitroaromatic anilide (cluster 37)
    Alpha-Estradiol     20 heavy, MW 272.4, steroid              (cluster 1)

Matched to within 4 % on molecular weight and one heavy atom, but spanning a
planar fused aromatic, a single ring carrying strong electron-withdrawing groups,
and a rigid tetracyclic steroid with ZERO rotatable bonds. Size is held constant
so that any difference in loop dynamics is about chemistry, not bulk.

All three are Tian screen hits: functional sensors exist for them, but only with
mutations -- so they are genuinely non-cognate for WT PYR1, which is the premise.
All three are neutral, unlike ABA's -1, which removes a charge confound from
parameterisation.

RECEPTOR
--------
`data/md/S2_holo_closed/protein.pdb`, i.e. the EXACT protein used for S2. Not the
script-30 file, not a fresh 3QN1 pull. The comparison this feeds is "same closed
PYR1, different ligand", so the protein must be byte-identical to the one already
simulated with ABA, or the contrast is confounded at the outset.

POSES
-----
§27c registered that the three replicates per ligand must start from THREE
DIFFERENT DOCKED POSES rather than three seeds of one, so pose sensitivity is
measured instead of assumed. This script therefore keeps the top-scoring pose plus
the two best poses that differ from it by at least MIN_POSE_RMSD.

Run with the pyr1_docking env python (rdkit + vina); smina is taken from docking_env.
"""
import csv
import os
import subprocess
import sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(DATA, "noncognate")
os.makedirs(OUT, exist_ok=True)

RECEPTOR = os.path.join(DATA, "md", "S2_holo_closed", "protein.pdb")
ABA_REF = os.path.join(DATA, "md", "S2_holo_closed", "aba.pdb")
LIGCSV = os.path.join(DATA, "coumarin_benchmark", "tian_screen_ligands.csv")
SMINA = "/bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/smina"

LIGANDS = ["Imperatorin", "Flutamide", "Alpha-Estradiol"]
N_POSES = 3
MIN_POSE_RMSD = 1.5      # A, so the three starting poses are genuinely distinct
CLASH_CUTOFF = 2.2       # A, heavy-heavy; below this a pose is rejected outright
EXHAUSTIVENESS = 32


def heavy_coords(pdb, resname=None):
    xyz = []
    for l in open(pdb):
        if l.startswith(("ATOM", "HETATM")):
            el = l[76:78].strip().upper()
            if el == "H":
                continue
            if resname and l[17:20].strip() != resname:
                continue
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(xyz)


def build_3d(smiles, name):
    """SMILES -> single low-energy 3D conformer, written as SDF"""
    m = Chem.MolFromSmiles(smiles)
    assert m is not None, f"{name}: RDKit could not parse {smiles!r}"
    m = Chem.AddHs(m)
    cids = AllChem.EmbedMultipleConfs(m, numConfs=20, randomSeed=0xC0FFEE)
    assert len(cids) > 0, f"{name}: embedding failed"
    res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
    best = int(np.argmin([e for conv, e in res]))
    path = os.path.join(OUT, f"{name}.sdf")
    w = Chem.SDWriter(path)
    w.write(m, confId=cids[best])
    w.close()
    q = Chem.GetFormalCharge(Chem.MolFromSmiles(smiles))
    return path, q, Chem.RemoveHs(m).GetNumAtoms(), Chem.RemoveHs(Chem.Mol(m))


def rebuild_pose(pose, template):
    """
    smina writes poses as HEAVY ATOMS ONLY with unreliable bond orders, which makes
    antechamber die with "Weird atomic valence ... Possible open valence". So the
    docked coordinates are transplanted back onto a chemically correct molecule:
    bond orders come from the template by substructure match, hydrogens are added
    with coordinates, then ONLY the hydrogens are minimised while every heavy atom
    is held fixed -- so the docked pose is preserved exactly and the H geometry is
    still clean enough for AM1-BCC.
    """
    fixed = AllChem.AssignBondOrdersFromTemplate(template, pose)
    molh = Chem.AddHs(fixed, addCoords=True)
    mp = AllChem.MMFFGetMoleculeProperties(molh)
    if mp is not None:
        ff = AllChem.MMFFGetMoleculeForceField(molh, mp)
        for a in molh.GetAtoms():
            if a.GetAtomicNum() > 1:
                ff.AddFixedPoint(a.GetIdx())
        ff.Minimize(maxIts=500)
    return molh


def pose_rmsd(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


# ---------------------------------------------------------------- inputs
assert os.path.exists(RECEPTOR), RECEPTOR
assert os.path.exists(SMINA), f"smina not found at {SMINA}"
smi = {r["library_name"]: r for r in csv.DictReader(open(LIGCSV))}
for L in LIGANDS:
    assert L in smi, f"{L} not in {LIGCSV}"

aba = heavy_coords(ABA_REF)
prot = heavy_coords(RECEPTOR)
print(f"receptor : {RECEPTOR}")
print(f"           {len(prot)} protein heavy atoms")
print(f"autobox  : ABA, {len(aba)} heavy atoms, centre "
      f"{np.round(aba.mean(0), 2).tolist()}")

summary = []
for name in LIGANDS:
    print("\n" + "=" * 74)
    print(name)
    print("=" * 74)
    sdf, charge, nheavy, template = build_3d(smi[name]["canonical_smiles"], name)
    print(f"  {nheavy} heavy atoms, formal charge {charge}, cluster "
          f"{smi[name]['chem_cluster']}, {smi[name]['n_clones']} characterised clones")
    assert charge == 0, (
        f"{name} has formal charge {charge}; the runner assumes neutral ligands "
        f"and tleap would need counter-ions adjusted")

    out = os.path.join(OUT, f"{name}_docked.sdf")
    cmd = [SMINA, "-r", RECEPTOR, "-l", sdf,
           "--autobox_ligand", ABA_REF, "--autobox_add", "4",
           "-o", out, "--num_modes", "20",
           "--exhaustiveness", str(EXHAUSTIVENESS), "--seed", "1"]
    print("  docking ...")
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        print(p.stdout[-2000:]); print(p.stderr[-2000:])
        sys.exit(f"smina failed for {name}")

    poses = [m for m in Chem.SDMolSupplier(out, removeHs=False) if m is not None]
    assert poses, f"{name}: smina produced no poses"

    # keep the top pose, then the best poses that differ from those already kept
    kept = []
    for m in poses:
        xyz = np.array([list(m.GetConformer().GetAtomPosition(i))
                        for i in range(m.GetNumAtoms())
                        if m.GetAtomWithIdx(i).GetAtomicNum() > 1])
        score = float(m.GetProp("minimizedAffinity"))
        # reject anything clashing with the protein
        dmin = float(np.min(np.linalg.norm(xyz[:, None, :] - prot[None, :, :], axis=-1)))
        if dmin < CLASH_CUTOFF:
            continue
        # must sit in the ABA pocket, not on the surface
        d_aba = float(np.min(np.linalg.norm(xyz[:, None, :] - aba[None, :, :], axis=-1)))
        if d_aba > 4.0:
            continue
        if all(pose_rmsd(xyz, k[1]) >= MIN_POSE_RMSD for k in kept):
            kept.append((score, xyz, m, dmin))
        if len(kept) == N_POSES:
            break

    assert len(kept) == N_POSES, (
        f"{name}: only {len(kept)} acceptable, mutually distinct poses "
        f"(needed {N_POSES}). Loosen MIN_POSE_RMSD or inspect the docking.")

    for i, (score, xyz, m, dmin) in enumerate(kept):
        full = rebuild_pose(m, template)
        nH = sum(1 for a in full.GetAtoms() if a.GetAtomicNum() == 1)
        nheavy_out = sum(1 for a in full.GetAtoms() if a.GetAtomicNum() > 1)
        assert nheavy_out == nheavy, (
            f"{name} pose{i}: {nheavy_out} heavy atoms after rebuild, expected {nheavy}")
        assert nH > 0, f"{name} pose{i}: no hydrogens after rebuild"
        # the docked heavy-atom pose must be untouched by the H minimisation
        newxyz = np.array([list(full.GetConformer().GetAtomPosition(k))
                           for k in range(full.GetNumAtoms())
                           if full.GetAtomWithIdx(k).GetAtomicNum() > 1])
        shift = float(np.abs(np.sort(newxyz, axis=0) - np.sort(xyz, axis=0)).max())
        assert shift < 1e-3, (
            f"{name} pose{i}: heavy atoms moved {shift:.4f} A during H minimisation")
        w = Chem.SDWriter(os.path.join(OUT, f"{name}_pose{i}.sdf"))
        w.write(full)
        w.close()
        rms = [round(pose_rmsd(xyz, k[1]), 2) for k in kept]
        print(f"   pose{i}: affinity {score:6.2f} kcal/mol, closest protein contact "
              f"{dmin:.2f} A, RMSD to kept poses {rms}")
        summary.append(dict(ligand=name, pose=i, affinity=score,
                            min_protein_dist=dmin, n_heavy=nheavy))

print("\n" + "=" * 74)
print("SUMMARY")
print("=" * 74)
print(f"  {'ligand':<18}{'pose':>5}{'affinity':>10}{'min contact':>13}")
for s in summary:
    print(f"  {s['ligand']:<18}{s['pose']:>5}{s['affinity']:>10.2f}"
          f"{s['min_protein_dist']:>13.2f}")
print(f"\nwrote {OUT}/<ligand>_pose{{0,1,2}}.sdf")
print("NEXT: scripts/64_build_noncognate.sh  (antechamber + tleap per ligand/pose)")
