#!/usr/bin/env python
r"""
166_agro_scan.py -- Park's 475 scored against each agrochemical's OWN pose.

Jannis's design. §92d ranked the three predicted poses by seed spread:

    fludioxonil       0.71 A   one pose
    benzothiadiazole  2.52 A   uncertain orientation
    benoxacor         3.27 A   no consensus

All three are scored here, not just the confident one, so that POSE CONFIDENCE
IS AN INDEPENDENT VARIABLE. If AUC tracks the ordering, a confident pose is the
precondition for the method and the ordering is the evidence. If the three come
out alike, pose confidence is not what determines success -- and either way the
answer needs all three, because one compound alone cannot distinguish them.

PROTOCOL IS IDENTICAL TO 156, deliberately. Same K59R background, same per-
position 6 A repack shell, same paired same-shell background control, same three
replicates, same dG_bind split. Only the ligand and its pose change, so the
resulting AUCs are directly comparable with mandipropamid's 0.868 (§91d).

Positives per compound, from the decoded Fig 5 (§91b): fludioxonil 28,
benzothiadiazole 22, benoxacor 10, against 475 tested variants each. Fludioxonil
is the better-powered test than mandipropamid's 20, and unlike mandipropamid its
positives are not concentrated at a single dominant-clash position.

⚠ A predicted pose can be confidently wrong (§92d). Seed agreement is precision,
not accuracy, and only the mandipropamid arm has a crystal pose to anchor it.
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import lib_rosetta as LR                                            # noqa: E402

AP = os.path.join(ROOT, "data", "agro_params")
OUT = os.path.join(ROOT, "results", "agro_scan")
FRAME = os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")
PARK25 = {55: "P", 61: "F", 62: "I", 81: "V", 83: "V", 87: "L", 88: "P", 89: "A",
          92: "S", 94: "E", 141: "E", 108: "F", 110: "I", 115: "H", 116: "R",
          117: "L", 120: "Y", 122: "S", 158: "M", 159: "F", 160: "A", 162: "T",
          163: "V", 164: "V", 167: "N"}
THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
         "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
         "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
         "W": "TRP", "Y": "TYR"}
AA20 = "ACDEFGHIKLMNPQRSTVWY"


def build_frame(code):
    """protein from the PYR1 frame + the transplanted ligand, one PDB"""
    out = os.path.join(AP, f"frame_{code}.pdb")
    with open(out, "w") as fh:
        for l in open(FRAME):
            if l.startswith("ATOM"):
                fh.write(l)
        fh.write("TER\n")
        for l in open(os.path.join(AP, f"{code}_0001.pdb")):
            if l.startswith(("ATOM", "HETATM")):
                fh.write("HETATM" + l[6:])
        fh.write("END\n")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compound", required=True)
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--nchunks", type=int, default=1)
    ap.add_argument("--nrep", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    meta = json.load(open(os.path.join(AP, "summary.json")))[a.compound]
    code, prm = meta["code"], meta["params"]
    frame = build_frame(code)

    import pyrosetta
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ex1 -ex2aro "
                   f"-extra_res_fa {prm}")
    sfxn = pyrosetta.create_score_function("ref2015_cart")
    base = pyrosetta.pose_from_pdb(frame)
    info = base.pdb_info()
    idx = {info.number(i): i for i in range(1, base.total_residue() + 1)}
    lig = base.total_residue()
    assert base.residue(lig).name3() == code, \
        f"last residue is {base.residue(lig).name3()}, expected {code}"
    bad = {n: base.residue(idx[n]).name1() for n, w in PARK25.items()
           if base.residue(idx[n]).name1() != w}
    assert not bad, f"Park identities disagree: {bad}"
    assert base.residue(idx[59]).name1() == "K"
    MutateResidue(idx[59], "ARG").apply(base)
    print(f"{a.compound} ({code}): {base.residue(lig).natoms()} ligand atoms, "
          f"K59R background, pose spread {meta['spread']:.2f} A", flush=True)

    e_lig = float(sfxn(pyrosetta.pose_from_pdb(
        os.path.join(AP, f"{code}_0001.pdb"))))

    def dg(p):
        e = float(sfxn(p)); q = p.clone(); q.delete_residue_slow(q.total_residue())
        return e - float(sfxn(q)) - e_lig

    def repack(p0, allowed, nrep):
        v = []
        for _ in range(nrep):
            p = p0.clone()
            tf, _ = LR.restrict_packing(p, allowed)
            pk = PackRotamersMover(sfxn); pk.task_factory(tf); pk.apply(p)
            v.append(dg(p))
        return min(v), float(max(v) - min(v))

    jobs = [(n, mu) for n in sorted(PARK25) for mu in AA20 if mu != PARK25[n]]
    assert len(jobs) == 475
    jobs = [j for i, j in enumerate(jobs) if i % a.nchunks == a.chunk]
    rows, shell, bgc = [], {}, {}
    for num, mu in jobs:
        if num not in shell:
            sel = ResidueIndexSelector(str(idx[num]))
            nb = NeighborhoodResidueSelector(sel, 6.0, True)
            shell[num] = sorted({i + 1 for i, b in enumerate(nb.apply(base)) if b}
                                | {lig})
            bgc[num] = repack(base.clone(), shell[num], a.nrep)[0]
        p = base.clone()
        MutateResidue(idx[num], THREE[mu]).apply(p)
        s, sp = repack(p, shell[num], a.nrep)
        rows.append({"sub": f"{PARK25[num]}{num}{mu}", "pos": num, "mut": mu,
                     "dG": s, "bg": bgc[num], "ddG": s - bgc[num], "spread": sp})
        json.dump(rows, open(os.path.join(
            OUT, f"{a.compound}_{a.chunk:03d}.json"), "w"))
    print(f"wrote {len(rows)}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
