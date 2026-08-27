#!/usr/bin/env python
r"""
135_win_crossover.py -- the WIN 55,212-2 cross-over on 7MWN, an independent test.

WHY THIS ONE IS WORTH RUNNING (Jannis supplied the structure, 2026-08-26)
§78 found the mandipropamid cross-over is correct without relaxation but that
**99.8 % of the signal is F108A alone** -- one steric collision, which §32a already
detected geometrically without any energy function. K59R contributed 0.23 %. So the
open question is whether the score can do anything where the answer is NOT a large
clash.

7MWN is "An engineered PYL2-based WIN 55,212-2 synthetic cannabinoid sensor with a
stabilized HAB1 variant". Chain A is the sensor, chain B HAB1, WI5 the ligand
(WIN 55,212-2, C27H26N2O3, 32 heavy atoms). Diffing chain A against wild-type PYL2
(UniProt O80992) gives exactly THREE substitutions:

    K64Q   F165A   V166I

K64 is the position homologous to PYR1's K59 -- the residue that has defeated every
method here -- and F165A is the same phenylalanine-to-alanine clash-relief class as
F108A. So the decomposition asks directly: does the score again find only the F->A
and miss the K->Q?

⚠ PYL2 IS NOT PYR1. Numbering is offset by about +5 and the sequences differ, so
nothing about POSITIONS transfers. What transfers is the question -- can the score
separate a real sensor from its wild-type parent given the crystal pose -- on an
independent fold instance and an independent ligand.

⚠ ASYMMETRY, and how it is handled. The crystal IS the sensor, so its side chains
are experimental while the reverted wild-type's must be modelled. That biases `raw`
toward the sensor for a trivial reason. The fair comparisons are `repack` and
`relax`, which rebuild both cells identically; `raw` is reported with this caveat
attached rather than omitted.

⚠ Chain B (HAB1) is dropped: this is a pocket test, and §78's cells had no partner
either.

Ligand params are built from the CRYSTAL coordinates (§49's route: bond orders from
the CCD SMILES via AssignBondOrdersFromTemplate) rather than a generated conformer
(§130's route), because here the bound pose exists.
"""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_rosetta as LR                                        # noqa: E402

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "win_crossover")
#: two cells of a 2x2 within PYL2, SAME numbering (both verified against O80992):
#:   7MWN  WIN sensor  (K64Q/F165A/V166I per the depositors' own annotation) + WI5
#:   3KDI  WILD-TYPE PYL2 + ABA, 181 residues 7-187, zero gaps, 100 % identity
#: In each column the crystal cell is the one that SHOULD win, so `raw` is biased
#: toward the right answer in BOTH columns -- which is why repack/relax are the
#: comparisons that count. 3KDJ was rejected: it is PYL1 + ABI1, not PYL2 (27 %
#: identity), and would have silently supplied the wrong receptor.
SYSTEMS = {
    "win": {"cif": "7MWN.cif", "comp": "WI5", "crystal_is": "sensor",
            "smiles": "Cc1c(c2cccc3c2n1[C@@H](CO3)CN4CCOCC4)C(=O)c5cccc6c5cccc6"},
    "aba": {"cif": "3KDI.cif", "comp": "A8S", "crystal_is": "wildtype",
            "smiles": "CC1=CC(=O)CC([C@]1(/C=C/C(=CC(=O)O)C)O)(C)C"},
}
CIF = os.path.join(ROOT, "data", "7MWN.cif")
FASTA = os.path.join(ROOT, "data", "pyl2_wt.fasta")
M2P = ("/opt/linux/rocky/8.x/x86_64/pkgs/rosetta/2023.45/scripts/python/public/"
       "molfile_to_params.py")
RDKIT_PY = "/bigdata/cutlerlab/jjaco081/conda_envs/dockenv/bin/python"
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
ONE3 = {v: k for k, v in AA3.items()}
BB = {"N", "CA", "C", "O", "OXT"}
#: sensor -> wild-type reversion, in PYL2 numbering
REVERT = {64: ("GLN", "LYS"), 165: ("ALA", "PHE"), 166: ("ILE", "VAL")}


def read_cif(path):
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
        if p[cols["label_alt_id"]] not in (".", "A"):
            continue
        atom = p[cols["label_atom_id"]].strip('"')
        if atom.startswith("H") or p[cols["type_symbol"]].upper() == "H":
            continue
        xyz = np.array([float(p[cols["Cartn_x"]]), float(p[cols["Cartn_y"]]),
                        float(p[cols["Cartn_z"]])])
        comp = p[cols["label_comp_id"]]
        if comp == LIGCOMP:
            lig.append((atom, xyz))
        elif p[cols["auth_asym_id"]] == "A" and comp in AA3:
            num = int(p[cols["auth_seq_id"]])
            prot[(num, atom)] = xyz
            seq[num] = comp
    return prot, lig, seq


def build_params(lig):
    """Params from the CRYSTAL pose: bond orders from SMILES, coordinates from 7MWN."""
    os.makedirs(OUT, exist_ok=True)
    pdb = os.path.join(OUT, f"{LIGCOMP}_xtal.pdb")
    with open(pdb, "w") as fh:
        for i, (a, v) in enumerate(lig, start=1):
            nm = f" {a:<3s}" if len(a) < 4 else a
            fh.write(f"HETATM{i:5d} {nm}{'':1s}LIG X   1    "
                     f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}  1.00  0.00"
                     f"          {a[0]:>2s}\n")
        fh.write("END\n")
    script = os.path.join(OUT, "_mk.py")
    with open(script, "w") as fh:
        fh.write(f'''
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
raw = Chem.MolFromPDBFile({pdb!r}, removeHs=False, sanitize=False)
tmpl = Chem.MolFromSmiles({SMILES!r})
if raw is None or tmpl is None:
    sys.exit("read failed")
m = AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
m = Chem.AddHs(m, addCoords=True)
Chem.MolToMolFile(m, {os.path.join(OUT, f"{LIGCOMP}.mol")!r}, kekulize=False)
print(m.GetNumAtoms(), sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1))
''')
    r = subprocess.run([RDKIT_PY, script], capture_output=True, text=True,
                       env={**os.environ, "PYTHONPATH": ""})
    if r.returncode != 0:
        raise SystemExit(f"molfile build failed: {r.stdout}{r.stderr}")
    r2 = subprocess.run([sys.executable, M2P, "-n", LIGCOMP, "-p",
                         os.path.join(OUT, LIGCOMP), "--keep-names",
                         os.path.join(OUT, f"{LIGCOMP}.mol")],
                        capture_output=True, text=True, cwd=OUT)
    pf = os.path.join(OUT, f"{LIGCOMP}.params")
    if not os.path.exists(pf):
        raise SystemExit(f"molfile_to_params failed: {r2.stdout}{r2.stderr}")
    # verify the params PDB still sits on the crystal coordinates
    pp = os.path.join(OUT, f"{LIGCOMP}_0001.pdb")
    got = [np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
           for l in open(pp) if l.startswith(("ATOM", "HETATM"))
           and l[76:78].strip() != "H"]
    ref = [v for _a, v in lig]
    if len(got) != len(ref):
        raise SystemExit(f"params has {len(got)} heavy atoms, crystal {len(ref)}")
    d = max(float(np.linalg.norm(np.array(g) - np.array(rr)))
            for g, rr in zip(got, ref))
    return pf, pp, d, r.stdout.strip()


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", choices=("win", "aba"), default="win")
    args = ap.parse_args()
    SYS = SYSTEMS[args.system]
    global CIF, LIGCOMP, SMILES, CRYSTAL_IS, OUT
    CIF = os.path.join(ROOT, "data", SYS["cif"])
    LIGCOMP = SYS["comp"]
    SMILES = SYS["smiles"]
    CRYSTAL_IS = SYS["crystal_is"]
    OUT = os.path.join(ROOT, "results", f"crossover_{args.system}")
    os.makedirs(OUT, exist_ok=True)
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    prot, lig, seq = read_cif(CIF)
    wt = "".join(l.strip() for l in open(FASTA) if not l.startswith(">"))
    say("=" * 76)
    say("WIN 55,212-2 CROSS-OVER on 7MWN (PYL2 sensor + WI5)")
    say("=" * 76)
    say(f"   chain A: {len(seq)} residues; ligand WI5 {len(lig)} heavy atoms")
    diff = [(n, wt[n - 1], AA3[seq[n]]) for n in sorted(seq)
            if n <= len(wt) and wt[n - 1] != AA3[seq[n]]]
    say(f"   vs wild-type PYL2 (O80992): "
        + ", ".join(f"{a}{n}{b}" for n, a, b in diff))
    for n, (sensor_aa, wt_aa) in REVERT.items():
        want = sensor_aa if CRYSTAL_IS == "sensor" else wt_aa
        if seq.get(n) != want:
            raise SystemExit(f"position {n} is {seq.get(n)}, expected {want} "
                             f"(crystal is {CRYSTAL_IS})")
    say(f"   reversion map asserted: "
        + ", ".join(f"{AA3[v[0]]}{k}->{AA3[v[1]]}" for k, v in sorted(REVERT.items())))

    pf, pp, dmax, counts = build_params(lig)
    say(f"   params built from the CRYSTAL pose; max deviation "
        f"{dmax:.4f} A over {len(lig)} heavy atoms ({counts} atoms total/heavy)")

    import pyrosetta
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ex1 -ex2aro "
                   f"-extra_res_fa {pf}")
    sfxn = pyrosetta.create_score_function("ref2015_cart")

    base_pdb = os.path.join(OUT, "sensor.pdb")
    n = 0
    with open(base_pdb, "w") as fh:
        for (num, atom), v in sorted(prot.items()):
            n += 1
            nm = f" {atom:<3s}" if len(atom) < 4 else atom
            fh.write(f"ATOM  {n:5d} {nm}{'':1s}{seq[num]:>3s} A{num:4d}    "
                     f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}  1.00  0.00"
                     f"          {atom[0]:>2s}\n")
        fh.write("TER\n")
        fh.writelines(l for l in open(pp) if l.startswith(("ATOM", "HETATM")))
        fh.write("END\n")
    pose0 = pyrosetta.pose_from_pdb(base_pdb)
    info = pose0.pdb_info()
    idx = {info.number(i): i for i in range(1, pose0.total_residue() + 1)}
    ligres = pose0.total_residue()
    if pose0.residue(ligres).name3().strip() != LIGCOMP:
        raise SystemExit(f"last residue {pose0.residue(ligres).name3()}, not WI5")
    say(f"   pose: {pose0.total_residue()} residues including WI5")

    pocket = [p for p in idx if p != info.number(ligres)]
    L = pose0.residue(ligres)
    near = []
    for p in sorted(idx):
        if idx[p] == ligres:
            continue
        r = pose0.residue(idx[p])
        d = min(((r.xyz(i).x - L.xyz(j).x) ** 2 + (r.xyz(i).y - L.xyz(j).y) ** 2
                 + (r.xyz(i).z - L.xyz(j).z) ** 2) ** 0.5
                for i in range(1, r.natoms() + 1)
                for j in range(1, L.natoms() + 1))
        if d < 5.0:
            near.append(p)
    say(f"   {len(near)} residues within 5 A of WI5: {near}")
    sel = ResidueIndexSelector(",".join(str(idx[p]) for p in near))
    nb = NeighborhoodResidueSelector(sel, 6.0, True)
    allowed = sorted({i + 1 for i, b in enumerate(nb.apply(pose0)) if b} | {ligres})

    def dg(p):
        ec = float(sfxn(p))
        q = p.clone()
        q.delete_residue_slow(q.total_residue())
        lp = pyrosetta.pose_from_pdb(pp)
        return ec - float(sfxn(q)) - float(sfxn(lp))

    def cell(muts, label):
        """`muts` = positions to CHANGE away from the crystal. In the WIN system the
        crystal is the sensor, so changing means REVERTING to wild-type; in the ABA
        system the crystal is wild-type, so it means INSTALLING the sensor
        mutations. The target residue is picked accordingly."""
        p = pose0.clone()
        for num in muts:
            tgt = REVERT[num][1] if CRYSTAL_IS == "sensor" else REVERT[num][0]
            MutateResidue(idx[num], tgt).apply(p)
        raw = dg(p.clone())
        pr = p.clone()
        tf, _ = LR.restrict_packing(pr, allowed)
        pk = PackRotamersMover(sfxn)
        pk.task_factory(tf)
        pk.apply(pr)
        rep = dg(pr)
        px = p.clone()
        mm = MoveMap()
        mm.set_bb(False)
        mm.set_chi(False)
        for i in allowed:
            mm.set_bb(i, True)
            mm.set_chi(i, True)
        tf2, _ = LR.restrict_packing(px, allowed)
        fr = FastRelax(sfxn, 1)
        fr.cartesian(True)
        fr.min_type("lbfgs_armijo_nonmonotone")
        fr.set_movemap(mm)
        fr.set_task_factory(tf2)
        fr.apply(px)
        rel = dg(px)
        say(f"   {label:<22}{raw:>12.2f}{rep:>12.2f}{rel:>12.2f}")
        return raw, rep, rel

    say("")
    say(f"   {'cell':<22}{'raw':>12}{'repack':>12}{'relax':>12}")
    lab0 = "SENSOR (crystal)" if CRYSTAL_IS == "sensor" else "wild-type (crystal)"
    lab1 = "wild-type PYL2" if CRYSTAL_IS == "sensor" else "SENSOR (modelled)"
    sensor = cell(set(), lab0)
    wtcell = cell(set(REVERT), lab1)
    singles = {k: cell({k}, f"revert {AA3[REVERT[k][0]]}{k}{AA3[REVERT[k][1]]}")
               for k in sorted(REVERT)}
    say("")
    winner = "sensor" if CRYSTAL_IS == "sensor" else "wild-type"
    say(f"   SELECTIVITY: dG_bind(crystal cell) - dG_bind(other); "
        f"{winner} should WIN (negative)")
    for i, proto in enumerate(("raw", "repack", "relax")):
        d = sensor[i] - wtcell[i]
        note = "  (raw is biased toward the sensor -- see header)" if proto == "raw" else ""
        say(f"     {proto:<8}{d:+10.2f} REU   "
            f"{'CORRECT' if d < 0 else 'WRONG'}{note}")
    say("")
    say("   DECOMPOSITION: cost of reverting ONE position (positive = that mutation")
    say("   is doing work; compare against the full reversion)")
    for k in sorted(REVERT):
        for i, proto in enumerate(("raw", "repack", "relax")):
            pass
        a, b = (REVERT[k][0], REVERT[k][1]) if CRYSTAL_IS == "sensor" \
            else (REVERT[k][1], REVERT[k][0])
        say(f"     {AA3[a]}{k}{AA3[b]:<8}"
            + "".join(f"{singles[k][i]-sensor[i]:+12.2f}" for i in range(3)))
    say(f"     {'ALL THREE':<12}"
        + "".join(f"{wtcell[i]-sensor[i]:+12.2f}" for i in range(3)))
    with open(os.path.join(OUT, "win_crossover.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    with open(os.path.join(OUT, "win.json"), "w") as fh:
        json.dump({"sensor": sensor, "wt": wtcell,
                   "singles": {str(k): v for k, v in singles.items()}}, fh, indent=1)
    say(f"\n   written to {OUT}/win_crossover.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
