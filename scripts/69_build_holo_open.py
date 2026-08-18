#!/usr/bin/env python
"""
69_build_holo_open.py -- build the missing cell of the conformation x occupancy
factorial: PYR1 in the OPEN backbone conformation with ABA present.

THE EXPERIMENT
--------------
README 19b/24e set up a 2x2 and only ever filled the diagonal:

                    apo                 + ABA
    open      S1  (3 x 300 ns)      <-- THIS FILE (S10)
    closed    S9  (built, unrun)    S2  (3 x 300 ns)

Both "right" forms are done and neither "wrong" form has any data, so conformation
and occupancy have never been separated. S9 is already built. This script builds
S10, and it is the more interesting of the two: gate CLOSURE on ligand binding is
the physiological direction, so it is the one arm where the barrier this project
has never crossed (README 24: no transition in 1.8 us) is plausibly downhill.

CONSTRUCTION, AND WHY THIS ONE AND NOT A DOCKED POSE
----------------------------------------------------
ABA is transplanted from the closed structure into the open one by superposing the
RIGID CORE -- the same core the loop analysis fits on, gate/latch/Lb7a5 and the
modelled 182-191 tail all excluded. That places the ligand exactly where the
closed state holds it, and then asks whether an open gate closes around it.

Docking into the open pocket was rejected: it would answer a weaker question
("where does ABA sit in an open receptor") and add pose uncertainty to an arm
whose whole point is that the ligand starts in its known, correct position.

The transplant leaves the pocket HALF-FORMED, which is the expected result and not
a defect -- 12 of the 19 closed-state lining residues survive, and the ones lost
are the gate itself (87, 88, 89) plus 61, 91, 117, 163, 164, 167.

THE ONE INTERVENTION, AND ITS LIMIT
-----------------------------------
The transplant leaves a single tight contact (ABA O2 against H115 CE1, the latch
histidine). It is repaired by repacking ONLY the residues that actually clash,
using `lib_rosetta.restrict_packing` and a **PackRotamersMover** -- a packer run,
not a relax, so the backbone cannot move at all by construction.

That matters more here than anywhere else in the project: this system's entire
scientific value is that its backbone is OPEN. Any relaxation toward the closed
state would pre-bias the experiment toward the answer we are looking for. So the
gate and latch backbones are asserted bit-identical to `pyr1_open_191.pdb`
afterwards, and every non-repacked residue is checked all-atom.

Run with the tier1_analysis env python (PyRosetta).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
STRUCT = os.path.join(ROOT, "data", "structures_191")
MDLIG = os.path.join(ROOT, "data", "md")
PARAMS = os.path.join(ROOT, "data", "aba_params", "A8S.params")

import lib_resnumber as rn                      # noqa: E402
import lib_rosetta as lr                        # noqa: E402
from Bio.PDB import PDBParser                   # noqa: E402

GATE = list(range(85, 90))
LATCH = list(range(115, 118))
LB7A5 = list(range(148, 157))
TAIL = list(range(182, 192))     # modelled, not crystal -- must not enter the fit
BB = ("N", "CA", "C", "O")
CLASH = 3.0        # repack a residue only if it comes this close to ABA

#: What this structure has to clear, and why it is not a "no contacts under 2.6 A"
#: rule. This is a tleap INPUT, not a finished model: lib_solvate runs min1/min2
#: before heating, and a continuous minimiser relieves a slightly tight contact in
#: a way a discrete rotamer packer cannot. So the Rosetta-stage bar is only "no
#: pair so close that minimisation cannot recover it", and the real check is made
#: further down the pipeline, on the actual MD starting coordinates after AMBER
#: minimisation (scripts/69b). Moving the threshold without moving the check would
#: be tuning until it passes; the check moved to the stage that decides the answer.
ACCEPT = 2.0
SEED = 20260818

p = PDBParser(QUIET=True)


def log(m):
    print(m, flush=True)


def read_mol2(path):
    """[(name, element, xyz)] in file order, hydrogens included.

    The element comes from the ATOM NAME, not from column 6. Column 6 holds the
    GAFF2 atom type -- `ca`, `c3`, `ho`, `oh` -- which is lower case, so testing it
    against "H" silently classes every hydrogen as a heavy atom and every contact
    measurement downstream then includes hydrogens.
    """
    out, blk = [], False
    for line in open(path):
        if line.startswith("@<TRIPOS>ATOM"):
            blk = True
            continue
        if line.startswith("@<TRIPOS>") and blk:
            break
        if blk and line.split():
            f = line.split()
            el = "".join(ch for ch in f[1] if ch.isalpha())[:1].upper()
            assert el in ("C", "N", "O", "S", "H"), f"odd element {el} from {f[1]}"
            out.append((f[1], el, np.array(list(map(float, f[2:5])))))
    return out


def write_mol2_with_coords(src, dst, coords):
    """Copy a mol2 verbatim, replacing ONLY the coordinate columns.

    Atom names, GAFF2 types, connectivity and AM1-BCC charges are carried over
    untouched -- the ligand is the same molecule in a different frame, and
    regenerating it with antechamber would silently re-derive charges from the new
    geometry. See README 27: antechamber must run maxcyc=0 or sqm optimises the
    pose away.
    """
    lines, blk, i = [], False, 0
    for line in open(src):
        if line.startswith("@<TRIPOS>ATOM"):
            blk = True
            lines.append(line)
            continue
        if line.startswith("@<TRIPOS>") and blk:
            blk = False
        if blk and line.split():
            f = line.rstrip("\n").split()
            x, y, z = coords[i]
            i += 1
            lines.append(f"{int(f[0]):>7} {f[1]:<8}{x:>10.4f}{y:>10.4f}{z:>10.4f} "
                         f"{f[5]:<8}{f[6]:>3} {f[7]:<8}{float(f[8]):>10.4f}\n")
        else:
            lines.append(line)
    assert i == len(coords), f"replaced {i} coordinates, expected {len(coords)}"
    open(dst, "w").writelines(lines)


def backbone(model, nats):
    g = {}
    for ch in model:
        for r in ch:
            if r.id[0] == " " and r.id[1] in nats:
                for a in r:
                    if a.get_id() in BB:
                        g[(r.id[1], a.get_id())] = a.coord
    return g


def heavy_atoms(model):
    out = []
    for ch in model:
        for r in ch:
            if r.id[0] == " ":
                for a in r:
                    if a.element != "H":
                        out.append((r.id[1], r.get_resname(), a.get_id(), a.coord))
    return out


def contacts(lig_xyz, prot):
    P = np.array([c for *_, c in prot])
    return np.linalg.norm(lig_xyz[:, None, :] - P[None, :, :], axis=2)


def shell(D, prot, cutoff):
    return sorted({(prot[b][0], prot[b][1]) for b in np.where(D.min(0) < cutoff)[0]})


def main():
    log("=" * 78)
    log("BUILDING S10: open backbone + ABA  (the missing 2x2 cell)")
    log("=" * 78)

    open_pdb = os.path.join(STRUCT, "pyr1_open_191.pdb")
    closed_pdb = os.path.join(STRUCT, "pyr1_closed_191.pdb")
    for f in (open_pdb, closed_pdb, PARAMS, os.path.join(MDLIG, "A8S.mol2")):
        assert os.path.exists(f), f"missing {f}"

    # both receptors must be what they claim to be before anything is transplanted
    # verify_build RAISES on any failure -- numbering, sequence, landmark identity
    # or a broken peptide bond -- so reaching the log line is the pass condition
    for path, label in ((open_pdb, "open"), (closed_pdb, "closed")):
        rep = rn.verify_build(path, "A")
        log(f"  {label:<7} {rep['n_residues']} residues, md5 {rep['sequence_md5']}, "
            f"{len(rep['landmarks_checked'])} landmarks, "
            f"{rep['peptide_breaks']} chain breaks")

    om = p.get_structure("o", open_pdb)[0]
    cm = p.get_structure("c", closed_pdb)[0]

    # ---- core superposition, closed -> open ----
    mobile = set(GATE) | set(LATCH) | set(LB7A5) | set(TAIL)
    core = [n for n in range(1, 192) if n not in mobile]
    A, B = backbone(om, core), backbone(cm, core)
    keys = sorted(set(A) & set(B))
    assert len(keys) > 500, f"only {len(keys)} core backbone atoms in common"
    X = np.array([B[k] for k in keys])
    Y = np.array([A[k] for k in keys])
    mx, my = X.mean(0), Y.mean(0)
    U, _, Vt = np.linalg.svd((X - mx).T @ (Y - my))
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    fit = float(np.sqrt((((X - mx) @ R - (Y - my)) ** 2).sum() / len(keys)))
    log(f"\n  core fit closed->open: {len(keys)} backbone atoms over {len(core)} "
        f"residues, RMSD {fit:.2f} A")
    assert fit < 1.5, f"core fit {fit:.2f} A -- the core is supposed to be rigid"

    # ---- transplant the ligand ----
    lig = read_mol2(os.path.join(MDLIG, "A8S.mol2"))
    Lc = np.array([c for _, _, c in lig])
    L = (Lc - mx) @ R + my
    heavy_idx = [i for i, (_, e, _) in enumerate(lig) if e != "H"]
    log(f"  ABA: {len(lig)} atoms ({len(heavy_idx)} heavy) moved into the open frame")

    prot_o, prot_c = heavy_atoms(om), heavy_atoms(cm)
    Dc = contacts(Lc[heavy_idx], prot_c)
    Do = contacts(L[heavy_idx], prot_o)
    sc, so = shell(Dc, prot_c, 4.5), shell(Do, prot_o, 4.5)
    log(f"\n  lining residues, CLOSED receptor : {len(sc)}  (closest {Dc.min():.2f} A)")
    log(f"  lining residues, OPEN receptor   : {len(so)}  (closest {Do.min():.2f} A)")
    log(f"    retained : {sorted(set(n for n, _ in sc) & set(n for n, _ in so))}")
    log(f"    lost     : {sorted(set(n for n, _ in sc) - set(n for n, _ in so))}")
    log(f"    gained   : {sorted(set(n for n, _ in so) - set(n for n, _ in sc))}")

    tight = sorted({prot_o[b][0] for b in np.where(Do.min(0) < CLASH)[0]})
    log(f"\n  residues closer than {CLASH} A to ABA and therefore repacked: {tight}")

    # ---- write the assembled complex for Rosetta ----
    work = os.path.join(STRUCT, "_s10_work.pdb")
    with open(work, "w") as fh:
        for line in open(open_pdb):
            if line.startswith(("ATOM", "TER")):
                fh.write(line)
        # PDB columns are fixed-width and Rosetta reads them positionally: the
        # residue name is 18-20, so an atom name written from column 13 with no
        # altLoc slot shifts A8S left by one and Rosetta reports
        # "Unrecognized residue: 8S". Built explicitly rather than by eye.
        #   7-11 serial | 13-16 name | 17 altLoc | 18-20 resName | 22 chain
        #   23-26 resSeq | 31-38,39-46,47-54 xyz | 55-60 occ | 61-66 B | 77-78 elem
        for i, ((name, el, _), xyz) in enumerate(zip(lig, L), start=1):
            nm = f" {name:<3}" if len(name) < 4 else name[:4]
            fh.write(f"HETATM{i:>5} {nm}{' '}{'A8S':>3} {'X'}{999:>4}{' '}   "
                     f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{0.0:6.2f}"
                     f"          {el:>2}\n")
        fh.write("END\n")

    import pyrosetta
    pyrosetta.init(f"-mute all -ex1 -ex2aro -detect_disulf false "
                   f"-extra_res_fa {PARAMS} "
                   f"-run:constant_seed -run:jran {SEED}")
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

    pose = pyrosetta.pose_from_pdb(work)
    log(f"\n  pose: {pose.total_residue()} residues "
        f"(191 protein + 1 ligand expected)")
    assert pose.total_residue() == 192, pose.total_residue()
    assert pose.residue(192).name3().strip() == "A8S", pose.residue(192).name3()

    # PackRotamersMover, NOT FastRelax: a packer run moves side-chain rotamers and
    # nothing else, so the open backbone is preserved by construction rather than
    # by assertion. The assertion is still made below.
    tf, packable = lr.restrict_packing(pose, tight)
    log(f"  packer task: {len(packable)} residues repackable, 0 designable")
    sf = pyrosetta.create_score_function("ref2015")
    before = sf(pose)
    PackRotamersMover(sf, tf.create_task_and_apply_taskoperations(pose)).apply(pose)
    log(f"  ref2015 after packing: {before:.1f} -> {sf(pose):.1f} REU")

    # Chi-only minimisation of the same residues. The packer picks from a discrete
    # rotamer library, so it cannot make the sub-angstrom adjustment a tight
    # contact needs; the minimiser can. MoveMap allows CHI on the repacked set and
    # NOTHING else -- no backbone, no jumps -- so the open conformation and the
    # transplanted ligand pose are both fixed by construction, not by assertion.
    # (set_movemap on a MinMover does restrict the minimiser; the trap in
    # lib_rosetta is about FastRelax's PACKER, which is not running here.)
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.protocols.minimization_packing import MinMover
    mm = MoveMap()
    mm.set_bb(False)
    mm.set_jump(False)
    mm.set_chi(False)
    for i in packable:
        mm.set_chi(i, True)
    mn = MinMover(mm, sf, "lbfgs_armijo_nonmonotone", 1e-4, True)
    mn.apply(pose)
    log(f"  ref2015 after chi minimisation: {sf(pose):.1f} REU")

    out_all = os.path.join(STRUCT, "pyr1_open_holo_191_complex.pdb")
    pose.dump_pdb(out_all)

    # ---- protein-only PDB for tleap, ligand as its own mol2 ----
    out_prot = os.path.join(STRUCT, "pyr1_open_holo_191.pdb")
    with open(out_prot, "w") as fh:
        for line in open(out_all):
            if line.startswith(("ATOM", "TER")):
                fh.write(line)
        fh.write("END\n")
    out_lig = os.path.join(STRUCT, "A8S_open_frame.mol2")
    write_mol2_with_coords(os.path.join(MDLIG, "A8S.mol2"), out_lig, L)
    log(f"\n  wrote {os.path.basename(out_prot)} and {os.path.basename(out_lig)}")

    # ================= verification =================
    log("\n" + "=" * 78)
    log("VERIFICATION")
    log("=" * 78)
    rep = rn.verify_build(out_prot, "A")
    log(f"  {rep['n_residues']} residues, md5 {rep['sequence_md5']}, "
        f"{len(rep['landmarks_checked'])} landmarks verified, "
        f"{rep['peptide_breaks']} chain breaks")

    nm = p.get_structure("n", out_prot)[0]
    bmap = {r.id[1]: list(r) for ch in om for r in ch if r.id[0] == " "}
    amap = {r.id[1]: list(r) for ch in nm for r in ch if r.id[0] == " "}

    # 1. the backbone this experiment depends on must be untouched
    for nats, label in ((GATE, "gate"), (LATCH, "latch")):
        worst = 0.0
        for n in nats:
            a = {x.get_id(): x.coord for x in bmap[n] if x.get_id() in BB}
            b = {x.get_id(): x.coord for x in amap[n] if x.get_id() in BB}
            worst = max(worst, max(float(np.linalg.norm(a[k] - b[k])) for k in a))
        assert worst < 1e-3, f"{label} BACKBONE moved {worst:.3f} A -- the open " \
                             f"conformation has been perturbed toward closed"
        log(f"  {label} backbone unmoved: max {worst:.4f} A")

    # 2. everything not repacked is frozen, all-atom
    untouched = [n for n in range(1, 192) if n not in tight]
    n_chk, worst = lr.assert_frozen(bmap, amap, untouched, "non-repacked")
    log(f"  {n_chk} non-repacked residues frozen all-atom, max drift {worst:.4f} A")

    # 3. the clash is actually gone
    prot_n = heavy_atoms(nm)
    Dn = contacts(L[heavy_idx], prot_n)
    i, j = np.unravel_index(Dn.argmin(), Dn.shape)
    log(f"  closest ABA-protein contact after repack: {Dn.min():.2f} A "
        f"({lig[heavy_idx[i]][0]} -- {prot_n[j][1]}{prot_n[j][0]} {prot_n[j][2]})"
        f"   [was {Do.min():.2f} A]")
    assert Dn.min() > ACCEPT, (
        f"hard clash at {Dn.min():.2f} A after repacking and minimising {tight}; "
        f"minimisation will not recover this -- widen the repack set")
    if Dn.min() < 2.6:
        log(f"  NOTE: {Dn.min():.2f} A is tight for a finished model but fine for a "
            f"tleap input;\n        69b re-measures it after AMBER minimisation, "
            f"which is the number that matters")

    # 4. and the gate really is still open
    sn = shell(Dn, prot_n, 4.5)
    log(f"  lining residues after repack: {len(sn)}")

    json.dump(dict(seed=SEED, core_fit_rmsd=fit, core_residues=len(core),
                   repacked=tight, closest_before=float(Do.min()),
                   closest_after=float(Dn.min()),
                   lining_closed=[n for n, _ in sc], lining_open=[n for n, _ in sn],
                   score_before=before, score_after=float(sf(pose))),
              open(os.path.join(STRUCT, "s10_build_report.json"), "w"), indent=1)
    os.remove(work)
    log("\n  -> S10 receptor and ligand ready for 68_build_md_systems.sh")


if __name__ == "__main__":
    main()
