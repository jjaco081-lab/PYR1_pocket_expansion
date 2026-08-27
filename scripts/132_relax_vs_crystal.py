#!/usr/bin/env python
r"""
132_relax_vs_crystal.py -- does FastRelax reproduce the EXPERIMENTAL pocket, given
the right sequence and the right ligand?

THE DIAGNOSTIC (Jannis, 2026-08-26)
Everything ligand-present we have run so far -- FastDesign (§23j), coupled moves
(§36), MM-GBSA (§37-42), TI (§43-46) -- asked the protocol to RANK or DESIGN. None
asked the prior question: with the answer handed to it, does relaxation reproduce
the crystal? That separates two failure modes we have never separated:

    start from the CRYSTAL rotamers, relax, do they STAY?
        they drift  -> the SCORE prefers something else. A scoring problem.
        they stay   -> the score is fine on the right answer.

    start from REPACKED rotamers, relax, do they REACH the crystal?
        no  -> a SAMPLING problem (cf. §33b: crystal Arg59 sits at chi3 ~100 deg,
               which the backbone-independent rotamer library never proposes)
        yes -> both sampling and scoring are adequate for structure, and our
               failures live in the ranking task itself, not in the modelling

Run on 4WVO directly, in its own frame, so no superposition error enters. Waters
and Mg are stripped: PYR1 needs neither, and leaving ordered waters in would let
the relax exploit contacts a designed model would never have.

⚠ This measures STRUCTURAL accuracy only. A protocol can rebuild a pocket
faithfully and still rank mutations badly -- those are different claims, and this
script only speaks to the first.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_rosetta as LR                                        # noqa: E402

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "relax_vs_crystal")
CIF = os.path.join(ROOT, "data", "4WVO.cif")
#: stage-1 frame: a Rosetta-loadable PYR1 + mandipropamid pair. 4WVO itself cannot
#: be loaded directly -- its 3UZ carries CCD atom names (CAC, CAB, ...) while the
#: params built in §49 use C12/N1/..., and renaming by hand is exactly the silent
#: mismatch §97 warns about. §49 already married the CRYSTAL coordinates to the
#: params names, so 3UZ_0001.pdb IS the crystal ligand in this frame.
FRAME = os.path.join(ROOT, "data", "stage1", "wt_mandi.pdb")
LIGPDB = os.path.join(ROOT, "data", "stage1", "params", "3UZ_0001.pdb")
LIGPARAMS = os.path.join(ROOT, "data", "stage1", "params", "3UZ.params")
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
BB = {"N", "CA", "C", "O", "OXT"}
#: the four PYR1^MANDI substitutions; identity is asserted against the CIF
MUT = {59: "ARG", 81: "ILE", 108: "ALA", 159: "LEU"}
POCKET = [59, 81, 83, 87, 89, 92, 94, 108, 110, 117, 120, 122, 141, 159, 160,
          163, 164, 167]


def read_cif(path):
    """chain A protein atoms + the 3UZ ligand, from the mmCIF atom_site loop."""
    prot, lig, seq = {}, [], {}
    cols, inloop = {}, False
    for line in open(path):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".", 1)[1]] = len(cols)
            inloop = True
            continue
        if inloop and line.startswith("#"):
            inloop = False
            continue
        if not inloop or not cols:
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        p = line.split()
        if len(p) < len(cols):
            continue
        ch = p[cols["auth_asym_id"]]
        comp = p[cols["label_comp_id"]]
        atom = p[cols["label_atom_id"]].strip('"')
        alt = p[cols["label_alt_id"]]
        if alt not in (".", "A"):
            continue
        xyz = np.array([float(p[cols["Cartn_x"]]), float(p[cols["Cartn_y"]]),
                        float(p[cols["Cartn_z"]])])
        if comp == "3UZ":
            lig.append((atom, xyz))
        elif ch == "A" and comp in AA3:
            num = int(p[cols["auth_seq_id"]])
            prot[(num, atom)] = xyz
            seq[num] = comp
    return prot, lig, seq


def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    D = np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))])
    return V @ D @ Wt, P.mean(0), Q.mean(0)


def pdb_line(kind, serial, name, resn, chain, num, xyz):
    """PDB columns, written out explicitly.

    ⚠ The first version of this put resName in columns 17-19 and chainID in 21 --
    every field after the atom name shifted one left. Rosetta then failed with
    "too many tries in fill_missing_atoms" on a file that looked entirely normal.
    Same class as the apo-bug where a shifted HETATM column put `3` in altLoc and
    ProDy silently dropped the ligand. Columns, for the record:
      1-6 record, 7-11 serial, 13-16 name, 17 altLoc, 18-20 resName,
      22 chain, 23-26 resSeq, 31-38/39-46/47-54 xyz, 77-78 element
    """
    el = name.strip()[0]
    nm = f" {name:<3s}" if len(name.strip()) < 4 else name
    return (f"{kind:<6s}{serial:5d} {nm:<4s} {resn:>3s} {chain}{num:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00"
            f"          {el:>2s}\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    prot, lig, seq = read_cif(CIF)
    say("=" * 78)
    say("0. INPUT AND IDENTITY CHECK")
    say("=" * 78)
    say(f"   4WVO chain A: {len(seq)} residues, ligand 3UZ {len(lig)} atoms")
    bad = [(n, want, seq.get(n)) for n, want in MUT.items() if seq.get(n) != want]
    if bad:
        raise SystemExit(f"4WVO does not carry the expected mutations: {bad}")
    say(f"   mutations present as expected: "
        + ", ".join(f"{n}{AA3[MUT[n]]}" for n in sorted(MUT)))

    # write a clean PDB: chain A protein + ligand, no waters, no Mg
    clean = os.path.join(OUT, "4wvo_clean.pdb")
    # ⚠ 4WVO chain A runs 3-180 with residues 69-72 disordered. Every residue that
    # IS present is complete (0 missing backbone atoms, 0 truncated side chains),
    # but the gap is unmarked, so Rosetta tries to bridge 68->73 and dies with
    # "too many tries in fill_missing_atoms". A TER at the break makes it two
    # segments instead of one impossible chain.
    nums = sorted({k[0] for k in prot})
    gaps = {a for a, b in zip(nums, nums[1:]) if b != a + 1}
    say(f"   chain breaks marked with TER after residue(s): {sorted(gaps)}")
    n = 0
    prev = None
    with open(clean, "w") as fh:
        for (num, atom), xyz in sorted(prot.items()):
            if atom.startswith("H"):
                continue
            if prev is not None and num != prev and prev in gaps:
                fh.write("TER\n")
            prev = num
            n += 1
            el = atom[0]
            fh.write(f"ATOM  {n:5d} {atom:<4s}{seq[num]:>3s} A{num:4d}    "
                     f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00"
                     f"          {el:>2s}\n")
        fh.write("TER\n")
        for atom, xyz in lig:
            if atom.startswith("H"):
                continue
            n += 1
            fh.write(f"HETATM{n:5d} {atom:<4s}3UZ X   1    "
                     f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00"
                     f"          {atom[0]:>2s}\n")
        fh.write("END\n")
    say(f"   wrote {clean} ({n} heavy atoms, waters and Mg stripped)")

    import pyrosetta
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ex1 -ex2aro "
                   f"-extra_res_fa {LIGPARAMS}")
    pose = pyrosetta.pose_from_pdb(clean)
    sfxn = pyrosetta.create_score_function("ref2015_cart")
    info = pose.pdb_info()
    idx = {}
    for i in range(1, pose.total_residue() + 1):
        idx[info.number(i)] = i
    ligres = pose.total_residue()
    if pose.residue(ligres).name3().strip() != "3UZ":
        raise SystemExit(f"last residue is {pose.residue(ligres).name3()}, not 3UZ")

    sel = ResidueIndexSelector(",".join(str(idx[p]) for p in POCKET if p in idx))
    nb = NeighborhoodResidueSelector(sel, 6.0, True)
    allowed = sorted({i + 1 for i, b in enumerate(nb.apply(pose)) if b} | {ligres})
    say(f"   relaxable shell: {len(allowed)} residues (pocket + 6 A, ligand free)")

    def sidechain(pose_, num):
        i = idx[num]
        r = pose_.residue(i)
        out = {}
        for j in range(1, r.natoms() + 1):
            nm = r.atom_name(j).strip()
            if nm in BB or r.atom_type(j).element().strip() == "H":
                continue
            v = r.xyz(j)
            out[nm] = np.array([v.x, v.y, v.z])
        return out

    xtal = {p: sidechain(pose, p) for p in POCKET if p in idx}

    def rmsd_to_xtal(pose_, p):
        a, b = xtal[p], sidechain(pose_, p)
        keys = set(a) & set(b)
        if not keys:
            return None
        return float(np.sqrt(np.mean([np.sum((a[k] - b[k]) ** 2) for k in keys])))

    def run_relax(p0):
        mm = MoveMap()
        mm.set_bb(False)
        mm.set_chi(False)
        for i in allowed:
            mm.set_bb(i, True)
            mm.set_chi(i, True)
        tf, _ = LR.restrict_packing(p0, allowed)
        fr = FastRelax(sfxn, 1)
        fr.cartesian(True)
        fr.min_type("lbfgs_armijo_nonmonotone")
        fr.set_movemap(mm)
        fr.set_task_factory(tf)
        fr.apply(p0)
        return p0

    # ---- arm A: start from the crystal rotamers ---------------------------
    say("")
    say("=" * 78)
    say("A. START FROM THE CRYSTAL ROTAMERS -- do they STAY?  (tests SCORING)")
    say("=" * 78)
    pa = run_relax(pose.clone())
    ra = {p: rmsd_to_xtal(pa, p) for p in xtal}

    # ---- arm B: discard crystal side chains, repack, then relax -----------
    say("")
    say("=" * 78)
    say("B. START FROM REPACKED ROTAMERS -- do they REACH it?  (tests SAMPLING)")
    say("=" * 78)
    pb = pose.clone()
    tf, _ = LR.restrict_packing(pb, allowed)
    pk = PackRotamersMover(sfxn)
    pk.task_factory(tf)
    pk.apply(pb)
    rb0 = {p: rmsd_to_xtal(pb, p) for p in xtal}
    pb = run_relax(pb)
    rb = {p: rmsd_to_xtal(pb, p) for p in xtal}

    say("")
    say(f"   {'pos':<6}{'wt':<5}{'A: from xtal':>14}{'B: repacked':>13}"
        f"{'B after relax':>15}   ")
    for p in sorted(xtal):
        star = " <-- mutated in PYR1^MANDI" if p in MUT else ""
        say(f"   {p:<6}{AA3[seq[p]]:<5}{ra[p]:>14.2f}{rb0[p]:>13.2f}"
            f"{rb[p]:>15.2f}{star}")
    mut = [p for p in MUT if p in xtal]
    say("")
    say(f"   mean sidechain RMSD to crystal, all {len(xtal)} pocket positions:")
    say(f"     A (from crystal, relaxed) : {np.mean([ra[p] for p in xtal]):.2f} A")
    say(f"     B (repacked, before relax): {np.mean([rb0[p] for p in xtal]):.2f} A")
    say(f"     B (repacked, relaxed)     : {np.mean([rb[p] for p in xtal]):.2f} A")
    say(f"   the four PYR1^MANDI positions only:")
    say(f"     A {np.mean([ra[p] for p in mut]):.2f} A     "
        f"B {np.mean([rb[p] for p in mut]):.2f} A")
    say("")
    say("   READ:")
    say("     A small  -> the score is happy with the experimental answer")
    say("     A large  -> the score prefers something else: a SCORING problem")
    say("     B small  -> sampling finds the crystal unaided")
    say("     B large while A is small -> a SAMPLING problem, not scoring")
    with open(os.path.join(OUT, "relax_vs_crystal.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    with open(os.path.join(OUT, "rmsd.json"), "w") as fh:
        json.dump({"from_crystal": ra, "repacked": rb0, "repacked_relaxed": rb},
                  fh, indent=1)
    say(f"\n   written to {OUT}/relax_vs_crystal.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
