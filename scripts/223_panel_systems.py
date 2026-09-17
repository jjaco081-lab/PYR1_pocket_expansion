#!/usr/bin/env python
r"""
223_panel_systems.py -- turn 209's anchored poses into MM-GBSA inputs.

209 emits poses as {xyz, element} lists. MM-GBSA needs a ligand .mol WITH BOND
ORDERS at that pose, plus a receptor frame, for WT and for the known hit.

⚠ WHY THE .mol IS REBUILT FROM SMILES RATHER THAN WRITTEN BY 209. 179's header
records six different failures building a .mol at a pose: explicit hydrogens
defeat AssignBondOrdersFromTemplate, sanitize=False defeats the substructure
match, and MDL aromatic bond type 4 makes antechamber die. The robust route is
to build the molecule from the SMILES that generated the conformer in the first
place and only replace its COORDINATES.

That is only valid if the atom order matches. 209 does MolFromSmiles -> AddHs ->
Embed -> RemoveHs, and RemoveHs preserves heavy-atom order, so a fresh
MolFromSmiles(smiles) has the same heavy-atom order as 209's `el` list.
⚠ THAT IS AN ASSUMPTION, SO IT IS ASSERTED ELEMENT BY ELEMENT and the script
fails rather than writing a scrambled ligand.

Systems written: 4 ligands x {WT, known hit} = 8, each to be run at n=8 = 64.

⚠ POSE SELECTION IS A REAL DECISION HERE, NOT A DETAIL. 209 found eugenol
accepts 117 of 180 anchored placements at zero clashes, against anthrone's 6 of
36. For eugenol the "best" pose is one draw from a large near-degenerate set, so
its MM-GBSA number carries pose variance that the n=8 replicate error bar does
NOT see. This is recorded per system as `n_polygly_ok` so the eventual ddG is
read against it. It is the same class of error as 93c (seed spread understates
true error ~3x).
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m209 = import_module("209_anchor_place")
PANEL = m209.PANEL
POSES = os.path.join(ROOT, "results", "anchor_place")
OUT = os.path.join(ROOT, "data", "panel_mmgbsa")
#: 3-letter resnames, distinct from anything already parameterised
RESN = {"anthrone": "ANT", "3-phenylphenol": "PPH",
        "chloroxylenol": "CXL", "eugenol": "EUG"}


def write_mol(name, pose, path):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    smi = PANEL[name]["smiles"]
    m = Chem.MolFromSmiles(smi)
    assert m is not None, f"bad SMILES for {name}"
    el_want = [e.upper() for e in pose["el"]]
    el_have = [a.GetSymbol().upper() for a in m.GetAtoms()]
    assert el_have == el_want, (
        f"{name}: atom order differs between 209 and a fresh MolFromSmiles\n"
        f"  209  {el_want}\n  here {el_have}")
    xyz = np.array(pose["xyz"], float)
    assert xyz.shape == (m.GetNumAtoms(), 3), \
        f"{name}: {xyz.shape} coords for {m.GetNumAtoms()} atoms"
    conf = Chem.Conformer(m.GetNumAtoms())
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, [float(p[0]), float(p[1]), float(p[2])])
    m.RemoveAllConformers()
    m.AddConformer(conf, assignId=True)
    # ⚠ 209 returns HEAVY ATOMS ONLY (it does Chem.RemoveHs before yielding the
    # conformer). antechamber refuses that outright -- "This molecule has no
    # hydrogens nor halogens ... Weird atomic valence (1) for atom C1" -- because
    # it reads the missing H as an open valence. Add them back at the pose, with
    # coordinates generated from the heavy-atom geometry.
    m = Chem.AddHs(m, addCoords=True)
    nH = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 1)
    assert nH > 0, "AddHs produced no hydrogens"
    Chem.Kekulize(m, clearAromaticFlags=True)     # antechamber rejects type 4
    Chem.MolToMolFile(m, path, kekulize=True)
    # assert on the RESULT: re-read and check the coordinates survived
    back = Chem.MolFromMolFile(path, removeHs=False, sanitize=False)
    assert back is not None, f"{name}: wrote an unreadable .mol"
    bx = back.GetConformer().GetPositions()[:len(xyz)]   # heavy atoms come first
    d = float(np.abs(bx - xyz).max())
    assert d < 1e-2, f"{name}: coords moved {d:.3f} A on round-trip"
    naro = sum(1 for b in back.GetBonds() if b.GetBondTypeAsDouble() == 1.5)
    return m.GetNumAtoms(), naro, d, nH


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="ketone")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    # receptor frame: 3QN1 chain A, protein heavy atoms, as 209 itself read it
    prot, lig = m209.read_3qn1()
    print(f"3QN1 chain A: {len(prot)} heavy protein atoms (identities asserted "
          f"by 209.read_3qn1)")

    rows, missing = [], []
    for name in PANEL:
        pj = os.path.join(POSES, f"{name}_{a.arm}.json")
        if not os.path.exists(pj):
            missing.append(pj); continue
        kept = json.load(open(pj))
        if not kept:
            missing.append(pj + " (empty)"); continue
        best = kept[0]
        d = os.path.join(OUT, name.replace("-", ""))
        os.makedirs(d, exist_ok=True)
        mp = os.path.join(d, f"{RESN[name]}.mol")
        na, naro, dev, nH = write_mol(name, best, mp)
        rows.append(dict(ligand=name, resname=RESN[name], mol=mp,
                         n_atom=na, n_H=nH, aromatic_left=naro, roundtrip_dev=dev,
                         fa_clashes=best["fa_clashes"],
                         worst_overlap=best["worst_overlap"],
                         n_pose_alternatives=len(kept),
                         hit=PANEL[name]["hit"], uM=PANEL[name]["uM"],
                         holdout=bool(PANEL[name].get("holdout"))))
        print(f"{name:<16} {RESN[name]}  {na:>2} atoms  "
              f"aromatic-left {naro}  +{nH} H  round-trip {dev:.4f} A  "
              f"clashes {best['fa_clashes']}  "
              f"alternatives {len(kept)}"
              + ("   [HELD OUT]" if PANEL[name].get("holdout") else ""))

    if missing:
        print("\nMISSING POSE FILES -- run 209 WITHOUT --dry-run first:")
        for m_ in missing:
            print("   ", m_)
        return 1

    json.dump(rows, open(os.path.join(OUT, "systems.json"), "w"), indent=1)
    print(f"\nwrote {OUT}/systems.json  ({len(rows)} ligands, "
          f"x {{WT, hit}} = {2*len(rows)} systems, x n=8 = {16*len(rows)} runs)")
    print("\n⚠ POSE DEGENERACY (read every ddG against this column):")
    for r in rows:
        print(f"   {r['ligand']:<16} {r['n_pose_alternatives']:>4} "
              f"clash-free alternatives to the pose being scored")
    return 0


if __name__ == "__main__":
    sys.exit(main())
