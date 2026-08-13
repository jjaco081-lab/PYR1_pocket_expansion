#!/usr/bin/env python
"""
49_stage1_ligand_params.py -- build Rosetta .params for the two stage-1 ligands
IN THEIR CRYSTAL POSES, so FastDesign can be run as a second opinion on whether
PYR1^MANDI is recoverable (section 23i).

WHY THIS IS A SEPARATE SCRIPT
-----------------------------
Every prior Rosetta run in this project (31_rosetta_switch.py and the arms before
it) was APO -- it compared conformer energies of the bare protein and never
needed a ligand residue type. Stage 1 is the first time Rosetta has to score a
protein-ligand interface here, and Rosetta cannot read a ligand out of a PDB the
way LigandMPNN can. It needs a residue topology file generated ahead of time.

THE PIPELINE, AND WHY EACH STEP IS THERE
----------------------------------------
A PDB ligand record carries elements and coordinates but no bond orders, and
Rosetta's typing depends on bond orders. The CCD chemical definition carries
bond orders but its `ideal` coordinates are a generated conformer, not the bound
one -- and this project's whole argument (sections 23d-e, scripts 45/46) is that
the BOUND conformer is the object of interest. So the two have to be married:

  1. crystal coordinates  <- the same HETATM block 47_stage1_inputs.py wrote,
                             i.e. 4WVO's mandipropamid already superposed into
                             the 3QN1 frame. Read from the file rather than
                             recomputed, so the Rosetta and LigandMPNN arms are
                             guaranteed to see one identical pose.
  2. bond orders          <- RDKit AssignBondOrdersFromTemplate against the CCD
                             SMILES. Preserves atom order and coordinates.
  3. hydrogens            <- AddHs(addCoords=True). Rosetta requires explicit H.
  4. aromatic perception  <- MolToMolBlock(kekulize=False), so the molfile
                             carries bond type 4. molfile_to_params warns and
                             mistypes aromatic rings if handed a Kekule
                             structure, and both of mandipropamid's rings sit
                             against the residues being designed.
  5. molfile_to_params.py <- Rosetta's own generator.

VERIFICATION -- the part that is not optional
---------------------------------------------
The params PDB is checked atom-for-atom against the crystal coordinates by
nearest-neighbour matching, and the run aborts unless every heavy atom matches
within POSE_TOL and the matching is one-to-one. molfile_to_params re-centres its
output on the input centroid, which is a no-op here but would not be if the
molfile were ever regenerated from an ideal conformer -- in that case this check
is what catches it. This is deliberately paranoid: the previous stage-1 result
was destroyed by a silent input defect that every downstream step tolerated.

ABA PROTONATION -- a real choice, not a default
-----------------------------------------------
The CCD gives ABA as the neutral acid. At assay pH it is the carboxylATE, and in
3QN1 that carboxylate salt-bridges K59. K59R is one of the four ground-truth
mutations and was the only position with apparent signal in the retracted run,
so this choice acts directly on the most informative residue: model ABA neutral
and Rosetta is scoring a hydrogen-bond where the real complex has a charged
interaction. Default here is DEPROTONATED. `--aba-neutral` builds the acid
instead, and running both is the honest sensitivity check.

Mandipropamid is neutral in any case, so the mandipropamid arms are unaffected.

Usage:
  /bigdata/cutlerlab/jjaco081/conda_envs/esmfold2/bin/python 49_stage1_ligand_params.py
"""
import argparse, os, subprocess, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = os.path.join(ROOT, "data", "stage1")
OUT = os.path.join(ROOT, "data", "stage1", "params")
M2P = ("/opt/linux/rocky/8.x/x86_64/pkgs/rosetta/2023.45/scripts/python/public/"
       "molfile_to_params.py")
POSE_TOL = 0.01          # A; molfile_to_params should be exact, so this is tight
# e; per-atom charges are written to 2 dp, so the sum carries up to 0.005 of rounding
# per atom. A fixed tolerance is wrong: 3UZ's 51 atoms sum to +0.100 while being a
# perfectly correct neutral molecule. The bound below scales with atom count and is
# still ~4x tighter than the 1.0 discrepancy a genuinely dropped formal charge gives.
def charge_tol(n_atoms):
    return 0.005 * n_atoms + 0.02

# CCD isomeric SMILES. Verified against files.rcsb.org/ligands/download/<id>.cif
SMILES = {
    "3UZ": "COc1cc(CCNC(=O)[C@@H](OCC#C)c2ccc(Cl)cc2)ccc1OCC#C",
    "A8S": r"CC(/C=C/[C@@]1(O)C(=CC(=O)CC1(C)C)C)=C/C(O)=O",
}
# deprotonated carboxylate; see the header
SMILES_ANION = {"A8S": r"CC(/C=C/[C@@]1(O)C(=CC(=O)CC1(C)C)C)=C/C([O-])=O"}
SOURCE = {"3UZ": "wt_mandi.pdb", "A8S": "wt_aba.pdb"}
N_HEAVY = {"3UZ": 29, "A8S": 19}

ap = argparse.ArgumentParser()
ap.add_argument("--aba-neutral", action="store_true",
                help="model ABA as the neutral acid instead of the carboxylate")
args = ap.parse_args()
os.makedirs(OUT, exist_ok=True)

from rdkit import Chem                                        # noqa: E402
from rdkit.Chem import AllChem                                 # noqa: E402
from rdkit import RDLogger                                     # noqa: E402
RDLogger.DisableLog("rdApp.warning")


def crystal_block(comp):
    """HETATM records for comp, straight from the stage-1 input file."""
    p = os.path.join(IN, SOURCE[comp])
    het = [l for l in open(p) if l.startswith("HETATM")]
    assert het, f"{p} contains no HETATM records"
    # guard against the 23h column bug reappearing upstream
    for l in het:
        assert l[16] == " " and l[17:20].strip() == comp, (
            f"{p}: malformed HETATM columns, altLoc={l[16]!r} "
            f"resName={l[17:20]!r} -- re-run 47_stage1_inputs.py")
    assert len(het) == N_HEAVY[comp], f"{p}: {len(het)} atoms, want {N_HEAVY[comp]}"
    return "".join(het)


def build_molfile(comp, smiles, path):
    """Crystal coordinates + CCD bond orders + explicit H + aromatic flags."""
    m = Chem.MolFromPDBBlock(crystal_block(comp), removeHs=False, sanitize=False)
    assert m is not None and m.GetNumAtoms() == N_HEAVY[comp]
    t = Chem.MolFromSmiles(smiles)
    assert t is not None, f"unparseable template SMILES for {comp}"
    assert t.GetNumAtoms() == N_HEAVY[comp], (
        f"{comp}: template has {t.GetNumAtoms()} heavy atoms, "
        f"crystal has {N_HEAVY[comp]}")
    fixed = AllChem.AssignBondOrdersFromTemplate(t, m)
    charge = Chem.GetFormalCharge(fixed)
    fh = Chem.AddHs(fixed, addCoords=True)
    Chem.SetAromaticity(fh, Chem.AromaticityModel.AROMATICITY_DEFAULT)
    blk = Chem.MolToMolBlock(fh, kekulize=False)     # bond type 4 for aromatics
    open(path, "w").write(blk)
    n_arom = sum(1 for b in fh.GetBonds() if b.GetIsAromatic())
    return dict(n_heavy=fixed.GetNumAtoms(), n_total=fh.GetNumAtoms(),
                n_arom=n_arom, charge=charge)


def params_net_charge(path):
    """
    Sum the partial charges in a .params ATOM block -- what Rosetta actually scores.

    Column order is:  ATOM <name> <rosetta_type> <mm_type> <charge>
    so the charge is field 5 (index 4). Field 4 is the MM type, which is the literal
    `X` for every ligand atom molfile_to_params writes. Summing field 4 by mistake
    yields exactly 0.000 for any ligand -- a number that looks like a real finding
    ("the formal charge was dropped") and is purely an off-by-one in the reader.
    That is why this returns the count as well: 38 atoms summing to 0.000 should read
    as suspicious, not as evidence.
    """
    tot, n = 0.0, 0
    for l in open(path):
        if l.startswith("ATOM "):
            f = l.split()
            assert len(f) >= 5, f"unexpected ATOM line in {path}: {l!r}"
            tot += float(f[4])
            n += 1
    return tot, n


def run_m2p(comp, molfile, prefix, formal_charge):
    """
    molfile_to_params DOES honour the molfile's `M  CHG` record here -- verified,
    not assumed. A8S_anion.mol carries `M  CHG  1  14  -1` and the params it produces
    sum to -0.970; 3UZ (a neutral molecule) sums to +0.100. Both are the intended
    formal charge to within the 2-decimal rounding of the per-atom charge column,
    which is why CHARGE_TOL is 0.05 and not 0.001.

    --recharge is therefore NOT passed. It would "ignore existing partial charges"
    (its own help text) and redistribute them uniformly, discarding the per-atom
    values in favour of a cruder model, to fix a problem that does not exist.

    The check below still earns its place: it confirms the ionisation state Rosetta
    will score, rather than the ionisation state of the RDKit molecule we handed in.
    Those are different objects, and only the first one affects the result.
    """
    cmd = [sys.executable, M2P, "-n", "LIG", "-p", prefix, "--clobber", molfile]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=OUT)
    if r.returncode != 0:
        sys.exit(f"molfile_to_params failed for {comp}:\n{r.stdout}\n{r.stderr}")

    out_params = os.path.join(OUT, prefix + ".params")
    net, n_at = params_net_charge(out_params)
    tol = charge_tol(n_at)
    if abs(net - formal_charge) > tol:
        sys.exit(f"{comp}: {os.path.basename(out_params)} partial charges sum to "
                 f"{net:+.3f} over {n_at} atoms, expected {formal_charge:+d} "
                 f"(tolerance {tol:.3f}). Rosetta would score the wrong ionisation "
                 f"state. Do not proceed.")
    print(f"    params net charge {net:+.3f} over {n_at} atoms "
          f"(formal {formal_charge:+d}, tol {tol:.3f})  OK")
    return r.stdout


def heavy_xyz(path, rec):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(path)
                     if l.startswith(rec) and l[76:78].strip() != "H"])


def verify_pose(comp, params_pdb):
    """Every heavy atom must land on a distinct crystal atom, within POSE_TOL."""
    a = heavy_xyz(os.path.join(IN, SOURCE[comp]), "HETATM")
    b = heavy_xyz(params_pdb, "HETATM")
    assert len(a) == len(b) == N_HEAVY[comp], (
        f"{comp}: {len(b)} params atoms vs {len(a)} crystal atoms")
    d = np.linalg.norm(b[:, None, :] - a[None, :, :], axis=-1)
    j = d.argmin(1)
    worst = float(d[np.arange(len(b)), j].max())
    assert len(set(j.tolist())) == len(b), (
        f"{comp}: params atoms do not map one-to-one onto crystal atoms")
    assert worst <= POSE_TOL, (
        f"{comp}: pose drifted {worst:.4f} A > {POSE_TOL} A -- the params "
        f"conformer is NOT the crystal conformer, which is the whole premise "
        f"of stage 1")
    return worst


print(f"molfile_to_params: {M2P}")
print(f"output: {OUT}\n")
for comp in ("3UZ", "A8S"):
    smi = SMILES[comp]
    tag = comp
    if comp == "A8S" and not args.aba_neutral:
        smi, tag = SMILES_ANION[comp], "A8S_anion"
    molfile = os.path.join(OUT, f"{tag}.mol")
    info = build_molfile(comp, smi, molfile)
    run_m2p(comp, molfile, tag, info["charge"])
    params = os.path.join(OUT, f"{tag}.params")
    pdb = os.path.join(OUT, f"{tag}_0001.pdb")
    assert os.path.exists(params), f"no params written for {tag}"
    worst = verify_pose(comp, pdb)
    print(f"{comp} -> {tag}.params")
    print(f"   {info['n_heavy']} heavy + H = {info['n_total']} atoms, "
          f"{info['n_arom']} aromatic bonds, formal charge {info['charge']:+d}")
    print(f"   pose check: max heavy-atom deviation {worst:.4f} A "
          f"(one-to-one over {N_HEAVY[comp]} atoms)\n")

print("done. Rosetta needs these on the command line as:")
for f in sorted(os.listdir(OUT)):
    if f.endswith(".params"):
        print(f"   -extra_res_fa {os.path.join(OUT, f)}")
