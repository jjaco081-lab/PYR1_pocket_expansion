#!/usr/bin/env python
r"""
131_ligand_aware_score.py -- the gate test for ligand-aware scoring.

THE QUESTION (Jannis, 2026-08-26)
§67/§69 scored variants protein-only and found real coumarin sensors separate from
random members of the same library at AUC 0.284. That score is blind to the ligand
by construction. Can we do better by docking the ligand and relaxing protein AND
ligand together -- and if so, is an iterative dock/relax/mutate loop worth building?

WHY A SINGLE PASS IS THE RIGHT GATE
An iterative loop re-applies the same energy function at every step. If one
dock+relax pass carries no discrimination, iterating cannot manufacture any: the
loop would be optimising a coordinate that does not track the answer. So the cheap
decisive experiment is one pass, scored the same way §69 was, on the same variants.

PROTOCOL PER VARIANT
  1. thread the mutations, Cartesian-relax the 8 A shell protein-only  (as §124)
  2. smina-dock the cognate coumarin into that relaxed pocket
  3. rebuild as a Rosetta pose with the ligand params, relax shell + ligand jointly
  4. score  dG_bind = E(complex) - E(protein alone) - E(ligand alone), all relaxed

LIGAND ASSIGNMENT, and why the nulls get one at all
REAL sensors are scored against their OWN target. LIBRARY and WILD variants have no
ligand of their own, so they are assigned one drawn from the same distribution as
REAL's, which keeps ligand identity balanced across the three sets. Without that
the comparison would confound "is this a sensor" with "which ligand was used".

⚠ THE POSE IS THE WEAK POINT AND IS NOT HIDDEN. No coumarin sensor has a crystal
structure, so the params carry a generated conformer (§130) and the pose comes from
docking into a designed pocket -- exactly the situation §46b measured at 8.33 A
error for mandipropamid. A NULL RESULT HERE IS THEREFORE AMBIGUOUS: it would mean
either that ligand-aware scoring adds nothing, or that the docking is too poor to
tell. A POSITIVE result is unambiguous and is what the test is for.

Usage: 131_ligand_aware_score.py --chunk i --nchunks n
"""
import argparse
import json
import os
import random
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_rosetta as LR                                        # noqa: E402
import importlib.util                                           # noqa: E402
_s = importlib.util.spec_from_file_location(
    "v122", os.path.join(HERE, "122_combination_viability.py"))
V = importlib.util.module_from_spec(_s)
_s.loader.exec_module(V)
_c = importlib.util.spec_from_file_location(
    "c125", os.path.join(HERE, "125_cavity_match.py"))
C = importlib.util.module_from_spec(_c)
_c.loader.exec_module(C)

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "ligand_aware")
PARAMS = os.path.join(ROOT, "data", "coumarin_params")
SMINA = "/bigdata/cutlerlab/jjaco081/conda_envs/docking_env/bin/smina"
THREE = V.THREE if hasattr(V, "THREE") else {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
    "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
    "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
    "W": "TRP", "Y": "TYR"}
#: ligand name -> params code, from 130
import importlib.util as _iu
_p = _iu.spec_from_file_location("p130", os.path.join(HERE, "130_coumarin_params.py"))
P130 = _iu.module_from_spec(_p)
_p.loader.exec_module(P130)
BYNAME = {v[0].lower(): k for k, v in P130.COUMARINS.items()}


def dock(receptor_pdb, code, seed, workdir):
    """smina into a box centred on the pocket seed. Returns the top pose PDB."""
    lig = os.path.join(PARAMS, f"{code}_0001.pdb")
    out = os.path.join(workdir, "docked.sdf")   # SDF keeps bond orders; PDB does not
    for f in (receptor_pdb, lig):
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            return None, f"INPUT MISSING/EMPTY: {f} exists={os.path.exists(f)}"
    cmd = [SMINA, "-r", receptor_pdb, "-l", lig, "-o", out,
           "--center_x", f"{seed[0]:.3f}", "--center_y", f"{seed[1]:.3f}",
           "--center_z", f"{seed[2]:.3f}",
           "--size_x", "16", "--size_y", "16", "--size_z", "16",
           "--exhaustiveness", "8", "--num_modes", "1", "--seed", "0",
           "--cpu", "1"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        return None, (r.stdout + r.stderr)[-160:]
    return out, None


RDKIT_PY = "/bigdata/cutlerlab/jjaco081/conda_envs/dockenv/bin/python"


def place_ligand(code, docked_sdf, out_pdb):
    """Graph-match via 131b, which runs under RDKit (not installed with PyRosetta)."""
    r = subprocess.run([RDKIT_PY, os.path.join(HERE, "131b_place_ligand.py"),
                        PARAMS, code, docked_sdf, out_pdb],
                       capture_output=True, text=True,
                       env={**os.environ, "PYTHONPATH": ""})
    if r.returncode != 0 or not os.path.exists(out_pdb):
        return None
    try:
        return int(r.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        return None


def write_clean_pdb(pose, path):
    """Dump ATOM/TER records only.

    ⚠ A pose that has been SCORED writes a `#BEGIN_POSE_ENERGIES_TABLE` block at
    the end of dump_pdb(), and Open Babel -- which is how smina reads a .pdb
    receptor -- rejects the entire file with "not a valid PDB file" and then
    "could not open ... for reading". The file exists, is non-empty, and looks
    perfectly normal in a text editor; only the consuming parser disagrees. The
    unrelaxed pose worked precisely because it had never been scored.
    """
    tmpf = path + ".raw"
    pose.dump_pdb(tmpf)
    n = 0
    with open(path, "w") as out:
        for line in open(tmpf):
            if line.startswith("#BEGIN_POSE_ENERGIES_TABLE"):
                break
            if line.startswith(("ATOM", "HETATM", "TER")):
                out.write(line)
                n += 1
        out.write("END\n")
    os.remove(tmpf)
    if n < 500:
        raise RuntimeError(f"only {n} atom records written to {path}")
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--nchunks", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    import pyrosetta
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    plist = ",".join(sorted(os.path.join(PARAMS, f)
                            for f in os.listdir(PARAMS) if f.endswith(".params")))
    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ex1 -ex2aro "
                   f"-extra_res_fa {plist.replace(',', ' ')}")
    base = pyrosetta.pose_from_pdb(V.PDB)
    sfxn = pyrosetta.create_score_function("ref2015_cart")
    seed, rms, ncom = C.pocket_seed(base)

    variants, _ = V.build_sets(base)
    # assign ligands: REAL keeps its own, nulls draw from REAL's distribution
    real_ligs = [s["lig"] for s in C.real_sensors("coumarin")]
    rng = random.Random(0)
    reals = [v for v in variants if v["set"] == "REAL"]
    for v, lg in zip(reals, real_ligs):
        v["lig"] = lg
    for v in variants:
        if "lig" not in v:
            v["lig"] = rng.choice(real_ligs)
    variants = [v for i, v in enumerate(variants) if i % a.nchunks == a.chunk]
    print(f"chunk {a.chunk}/{a.nchunks}: {len(variants)} variants; pocket seed from "
          f"{ncom} atoms (RMSD {rms:.2f} A)")

    def relax_shell(pose, allowed, with_lig):
        mm = MoveMap()
        mm.set_bb(False)
        mm.set_chi(False)
        for i in allowed:
            mm.set_bb(i, True)
            mm.set_chi(i, True)
        if with_lig:
            mm.set_chi(pose.total_residue(), True)
            mm.set_jump(1, True)
        idx = list(allowed) + ([pose.total_residue()] if with_lig else [])
        tf, _ = LR.restrict_packing(pose, sorted(set(idx)))
        fr = FastRelax(sfxn, 1)
        fr.cartesian(True)
        fr.min_type("lbfgs_armijo_nonmonotone")
        fr.set_movemap(mm)
        fr.set_task_factory(tf)
        fr.apply(pose)
        return float(sfxn(pose))

    rows = []
    tmp = tempfile.mkdtemp()
    for v in variants:
        code = BYNAME.get(v["lig"].lower())
        if code is None:
            print(f"  !! no params code for {v['lig']}")
            continue
        pose = base.clone()
        idx, ok = [], True
        for num, wt, mu in v["subs"]:
            one, i = V.wt_at(pose, num)
            if one != wt:
                ok = False
                break
            MutateResidue(i, THREE[mu]).apply(pose)
            idx.append(i)
        if not ok:
            continue
        sel = ResidueIndexSelector(",".join(str(i) for i in idx))
        nb = NeighborhoodResidueSelector(sel, 8.0, True)
        allowed = [i + 1 for i, b in enumerate(nb.apply(pose)) if b]
        e_prot = relax_shell(pose, allowed, False)
        rec = os.path.join(tmp, "rec.pdb")
        write_clean_pdb(pose, rec)
        dp, err = dock(rec, code, seed, tmp)
        if dp is None:
            print(f"  !! dock failed {v['id']} {v['lig']}: {err}")
            continue
        # ⚠ smina RENAMES the ligand to UNL, strips nonpolar hydrogens and
        # REORDERS the atoms (§97: template C1 C2 C3 C4 C5 O1 ... comes back as
        # C C C C O C ...). Copying coordinates by order would put an oxygen's
        # position on a carbon and nothing downstream would complain. So the
        # docked pose is matched to the parameterised molecule as a GRAPH and only
        # then written out with the params' own atom names.
        lig_pdb = os.path.join(tmp, "lig_placed.pdb")
        nplaced = place_ligand(code, dp, lig_pdb)
        if nplaced is None:
            print(f"  !! graph match failed {v['id']} {v['lig']}")
            continue
        comp = os.path.join(tmp, "complex.pdb")
        with open(comp, "w") as fh:
            for l in open(rec):
                if l.startswith("ATOM"):
                    fh.write(l)
            fh.write("TER\n")
            for l in open(lig_pdb):
                if l.startswith("HETATM"):
                    fh.write(l)
            fh.write("END\n")
        try:
            cpose = pyrosetta.pose_from_pdb(comp)
        except Exception as e:                                  # noqa: BLE001
            print(f"  !! pose build failed {v['id']}: {type(e).__name__} {e}")
            continue
        last = cpose.residue(cpose.total_residue()).name3().strip()
        if last != code:
            print(f"  !! last residue is {last}, expected {code}")
            continue
        e_cplx = relax_shell(cpose, allowed, True)
        # ligand alone, relaxed with the same function
        lp = pyrosetta.pose_from_pdb(os.path.join(PARAMS, f"{code}_0001.pdb"))
        e_lig = float(sfxn(lp))
        dg = e_cplx - e_prot - e_lig
        rows.append({"set": v["set"], "id": v["id"], "lig": v["lig"],
                     "code": code, "n_sub": len(v["subs"]),
                     "e_prot": e_prot, "e_cplx": e_cplx, "e_lig": e_lig,
                     "dG_bind": dg,
                     "subs": [f"{w}{n}{m}" for n, w, m in v["subs"]]})
        print(f"  {v['set']:<8}{v['id']:<12}{v['lig'][:20]:<22}"
              f"dG_bind {dg:9.2f}")
    with open(os.path.join(OUT, f"chunk_{a.chunk:03d}.json"), "w") as fh:
        json.dump({"rows": rows}, fh, indent=1)
    print(f"wrote {len(rows)} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
