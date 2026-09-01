#!/usr/bin/env python
r"""
165_agro_params.py -- transplant the Boltz-2 poses into the PYR1 frame and build
Rosetta params, for all THREE agrochemicals.

Jannis's design: score all three, not just the one with a clean pose, so that
pose confidence becomes an INDEPENDENT VARIABLE rather than a filter. §92d
ordered them by median pairwise ligand RMSD across 10 seeds:

    fludioxonil       0.82 A   one pose
    benzothiadiazole  2.47 A   same site, uncertain orientation
    benoxacor         3.32 A   no consensus

If ddG discrimination tracks that ordering, the pose is what determines whether
the method works and a confident pose is a usable precondition. If the three come
out alike, pose confidence is not the operative variable and something else is.
Either answer is worth having, and only running all three can tell them apart.

The representative pose is the MEDOID over the 10 seeds -- the one with the
smallest summed RMSD to the others -- not seed 0, which would be an arbitrary
draw from a distribution whose width is the whole point.

Ligand coordinates are transplanted into the 3QN1-derived PYR1 frame by
superposing the model's core backbone (gate and latch excluded), so all four
ligands including crystal mandipropamid sit in one common frame.
"""
import glob, json, os, subprocess, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from importlib import import_module
G = import_module("164_pose_geometry")

OUT = os.path.join(ROOT, "data", "agro_params")
RDKIT_PY = "/bigdata/cutlerlab/jjaco081/conda_envs/dockenv/bin/python"
M2P = ("/opt/linux/rocky/8.x/x86_64/pkgs/rosetta/2023.45/scripts/python/public/"
       "molfile_to_params.py")
SMILES = {"benzothiadiazole": "CSC(=O)c1cccc2nsnc12",
          "benoxacor": "CC1COc2ccccc2N1C(=O)C(Cl)Cl",
          "fludioxonil": "FC1(F)Oc2cccc(c2O1)-c1c[nH]cc1C#N",
          # mandipropamid with a PREDICTED pose -- the contamination-aware arm
          # (4WVO is in Boltz's training set, so success is uninformative and
          # only failure is conclusive; README 94d)
          "mandipropamid": "CC(=C)COc1ccc(cc1OC)CCNC(=O)C(OCC#C)c1ccc(Cl)cc1",
          # the two all-negative specificity compounds: mandipropamid-sized,
          # screened against all 475 with ZERO responders
          "azoxystrobin": r"CO/C=C(\C(=O)OC)c1ccccc1Oc1cc(Oc2ccccc2C#N)ncn1",
          "lufenuron": "O=C(Nc1cc(OC(F)(F)C(F)(F)F)cc(Cl)c1Cl)NC(=O)c1c(F)cccc1F"}
CODE = {"benzothiadiazole": "BZT", "benoxacor": "BNX", "fludioxonil": "FLD",
        "mandipropamid": "MDP", "azoxystrobin": "AZO", "lufenuron": "LUF"}


def medoid_pose(name, fbb, core):
    """Ligand coords of the medoid seed, already in the PYR1 frame."""
    poses, meta = [], []
    for s in range(10):
        g = sorted(glob.glob(os.path.join(ROOT, "data", "poses", "out",
                                          f"{name}_s{s}", "boltz_results_*",
                                          "predictions", "*", "*.cif")))
        if not g:
            continue
        bb, L, _, pl = G.read(g[0])
        P, Q = G.paired(fbb, bb, core)
        R, mp, mq = G.kabsch(P, Q)
        poses.append((L - mq) @ R + mp)
        meta.append((s, pl, float(np.sqrt((((Q - mq) @ R + mp - P) ** 2).sum(1).mean()))))
    n = min(len(p) for p in poses)
    D = np.array([[float(np.sqrt(((a[:n] - b[:n]) ** 2).sum(1).mean()))
                   for b in poses] for a in poses])
    k = int(np.argmin(D.sum(1)))
    return poses[k], meta[k], float(np.median(D[np.triu_indices(len(D), 1)])), len(poses)


def elements(name):
    g = sorted(glob.glob(os.path.join(ROOT, "data", "poses", "out", f"{name}_s0",
                                      "boltz_results_*", "predictions", "*", "*.cif")))
    cols, rows, inl = {}, [], False
    for line in open(g[0]):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".")[1]] = len(cols); inl = True; continue
        if inl:
            if line.startswith(("#", "loop_", "_")):
                if rows: break
                continue
            f = line.split()
            if len(f) >= len(cols): rows.append(f)
    return [f[cols["type_symbol"]] for f in rows
            if f[cols["group_PDB"]] != "ATOM" and f[cols["type_symbol"]] != "H"]


def main():
    os.makedirs(OUT, exist_ok=True)
    fbb, aba, _ = G.frame()
    core = [r for r in fbb if 6 <= r <= 180 and r not in G.GATE + G.LATCH]
    summary = {}
    for name, smi in SMILES.items():
        lig, (seed, pl, crms), spread, nseed = medoid_pose(name, fbb, core)
        el = elements(name)
        assert len(el) == len(lig), f"{name}: {len(el)} elements vs {len(lig)} atoms"
        code = CODE[name]
        pdb = os.path.join(OUT, f"{code}_pose.pdb")
        with open(pdb, "w") as fh:
            for i, (e, v) in enumerate(zip(el, lig), 1):
                an = f"{e}{i}"
                # column 17 is altLoc; omitting it puts resName in 17-19
                # instead of 18-20 and RDKit then reads ZERO atoms. This is the
                # same shift that voided a whole arm in README 23g.
                fh.write(f"HETATM{i:5d} {an:<4s}{'':1s}{code:>3s} X   1    "
                         f"{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}  1.00  0.00"
                         f"          {e:>2s}\n")
            fh.write("END\n")
        mk = os.path.join(OUT, f"_mk_{code}.py")
        open(mk, "w").write(f'''
import sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
raw = Chem.MolFromPDBFile({pdb!r}, removeHs=False, sanitize=False)
tmpl = Chem.MolFromSmiles({smi!r})
if raw is None or tmpl is None: sys.exit("read failed")
m = AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
m = Chem.AddHs(m, addCoords=True)
Chem.MolToMolFile(m, {os.path.join(OUT, code + ".mol")!r}, kekulize=False)
print(m.GetNumAtoms(), sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1))
''')
        r = subprocess.run([RDKIT_PY, mk], capture_output=True, text=True,
                           env={**os.environ, "PYTHONPATH": ""})
        if r.returncode != 0:
            print(f"  {name}: MOL BUILD FAILED\n{r.stdout}{r.stderr}"); continue
        subprocess.run([sys.executable, M2P, "-n", code, "-p",
                        os.path.join(OUT, code), "--keep-names",
                        os.path.join(OUT, f"{code}.mol")],
                       capture_output=True, text=True, cwd=OUT)
        pf = os.path.join(OUT, f"{code}.params")
        if not os.path.exists(pf):
            print(f"  {name}: molfile_to_params FAILED"); continue
        # the params PDB must still sit on the transplanted pose
        pp = os.path.join(OUT, f"{code}_0001.pdb")
        got = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                        for l in open(pp) if l.startswith(("ATOM", "HETATM"))
                        and l[76:78].strip() != "H"])
        dev = float(np.abs(np.sort(got, 0) - np.sort(lig, 0)).max())
        assert dev < 0.01, f"{name}: params moved the pose by {dev:.3f} A"
        summary[name] = dict(code=code, seed=seed, plddt=pl, core_rmsd=crms,
                             spread=spread, n_seeds=nseed, atoms=len(lig),
                             params=pf, pose_pdb=pp,
                             d_site=float(np.linalg.norm(lig.mean(0) - aba.mean(0))))
        print(f"  {name:<18} {code}  medoid seed {seed}  pLDDT {pl:.1f}  "
              f"spread {spread:.2f} A  {len(lig)} atoms  "
              f"d(ABA site) {summary[name]['d_site']:.2f} A  params OK")
    json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    print(f"\n{len(summary)}/3 built -> {OUT}/summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
