#!/usr/bin/env python
r"""
205_predicted_cavity.py -- does the predicted pocket CLOSE onto a small ligand,
or stay open with residual cavity behind it?

JANNIS'S QUESTION, and his caveat. He asked whether small ligands sit at the
mouth with cavity retained behind them, while noting the structural predictors
have been poor at ligand LOCATION and that the CAVITY measurement is what he
actually wants.

⚠ THE POSE HALF IS ALREADY ANSWERED AND IT IS FLAT (README 93a, n = 362 Boltz-2
co-folds): distance to the gate/latch mouth is 8.07 A for hits vs 8.00 A for
non-hits (AUC 0.508), every ligand lands on the ABA site (centroid displacement
1.3 A both classes) and the enclosed fraction is 1.00 for both. Co-folding is not
docking -- the structure module puts whatever it is given into the obvious pocket
and reports no doubt. So this script does NOT re-ask where the ligand sits.

WHAT IS ACTUALLY NEW. Nobody measured the CAVITY of those 1,810 structures. Per
model, with the ligand stripped and then present:

  apo_cav    enclosed volume of the predicted protein alone
  holo_cav   enclosed volume with the ligand present = RESIDUAL unfilled space
  fill       1 - holo/apo
  n_comp     components of the residual, i.e. one void or several

They already exist on disk (data/tractability/out, 362 ligands x 5 models), so
this costs no GPU.

RECEPTOR VERIFIED BY IDENTITY, not by filename: chain A is 191-residue WILD-TYPE
PYR1 and all 13 checked sd03 positions match WT, despite a `pyr1_k59r.a3m`
sitting in the same directory.

PREDICTIONS, FIXED BEFORE THE RUN:
  1. CALIBRATION. The predicted apo cavity must land near crystal 3QN1's, which
     this script measures with the identical code. If the predictions do not
     reproduce PYR1's own pocket geometry, nothing below is interpretable and
     the arm stops here.
  2. SEED STABILITY. Cavity must be more reproducible across the 5 diffusion
     samples than the POSE is (13-18 A across seeds for ESMFold2, and Boltz pose
     spread median 1.80-2.06 A in 93b). Reported as per-ligand CV; a cavity CV
     comparable to the between-ligand spread means the measurement is noise.
  3. THE QUESTION. If the pocket closes onto whatever it is given, residual
     cavity is ~flat in ligand size. If small ligands leave unfilled volume,
     residual RISES as heavy-atom count falls.

⚠ Prediction 3 is confounded by prediction 1's failure mode in one direction: a
model that always folds the same pocket regardless of ligand would give both a
flat residual AND a tight seed spread. The apo cavity's dependence on ligand size
is therefore reported too -- if apo_cav itself varies with the ligand, the
backbone is responding; if not, the pocket is a fixed prior and residual is just
arithmetic on ligand volume.
"""
import argparse, csv, glob, json, os, sys
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m151 = import_module("151_cavity_bottleneck")

OUT = os.path.join(ROOT, "results", "predicted_cavity")
PRED = os.path.join(ROOT, "data", "tractability", "out")
KEY = os.path.join(ROOT, "data", "tractability", "key_matched.json")
FIELDS = ["lig", "name", "label", "heavy", "clogp", "model",
          "apo_cav", "holo_cav", "apo_main", "holo_main", "fill",
          "apo_rmax", "holo_rmax", "n_comp"]
#: the 13 positions whose identity is asserted on every receptor parsed
POS = {59: "LYS", 81: "VAL", 83: "VAL", 87: "LEU", 89: "ALA", 108: "PHE",
       117: "LEU", 120: "TYR", 159: "PHE", 160: "ALA", 163: "VAL",
       164: "VAL", 167: "ASN"}


def read_cif(path):
    cols, rows, inl = {}, [], False
    for line in open(path):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".")[1]] = len(cols); inl = True; continue
        if inl:
            if line.startswith(("#", "loop_", "_")):
                if rows:
                    break
                continue
            x = line.split()
            if len(x) >= len(cols):
                rows.append(x)
    g = lambda r, k: r[cols[k]]                                  # noqa: E731
    prot, lig, seq = [], [], {}
    for r in rows:
        e = g(r, "type_symbol").upper()
        if e == "H":
            continue
        xyz = np.array([float(g(r, "Cartn_x")), float(g(r, "Cartn_y")),
                        float(g(r, "Cartn_z"))])
        if g(r, "auth_asym_id") == "A":
            prot.append((e, xyz))
            try:
                seq[int(g(r, "auth_seq_id"))] = g(r, "label_comp_id")
            except ValueError:
                pass
        else:
            lig.append((e, xyz))
    return prot, lig, seq


def cav(atoms):
    """(total enclosed A^3, main chamber A^3, r_max, n_components).

    ⚠ The first version returned only `main` and read residual volume as 0.0 for
    ALL 362 ligands. That was my bug. 151.analyse defines a CHAMBER as a
    component admitting a 2.4 A sphere; when none exists it returns
    note="no 2.4 A chamber" with no `main` key, so .get("main", 0.0) collapsed
    "the leftover void is too narrow for a 2.4 A probe" into "there is no
    leftover void". Eugenol actually leaves 200 enclosed voxels = 43.2 A^3,
    against the crystal ABA residual of 39.5 A^3. Both numbers are now returned:
    `total` is the volume, `main` is whether any of it is probe-accessible.
    """
    if len(atoms) < 50:
        return 0.0, 0.0, 0.0, 0
    xyz = [x for _, x in atoms]
    el = [e for e, _ in atoms]
    r = m151.analyse("x", xyz, el)
    if not r:
        return 0.0, 0.0, 0.0, 0
    try:
        from scipy import ndimage
        clear, _, lo, shape = m151.clear_map(xyz, el)
        mask = m151.enclosed(clear, lo, shape)
        n = int(ndimage.label(mask)[1]) if mask.sum() else 0
    except Exception:                                            # noqa: BLE001
        n = 0
    return (float(r.get("total", 0.0)), float(r.get("main", 0.0)),
            float(r.get("r_max", 0.0)), n)


def one(job):
    lig_id, meta, path = job
    model = os.path.basename(path).split("_model_")[1].split(".")[0]
    row = dict.fromkeys(FIELDS, "")
    row.update(lig=lig_id, name=meta["name"], label=meta["label"],
               heavy=meta["heavy"], clogp=round(meta["clogp"], 2), model=model)
    try:
        prot, lg, seq = read_cif(path)
        bad = [n for n, aa in POS.items() if seq.get(n) != aa]
        if bad:
            row["apo_cav"] = f"IDENTITY_FAIL:{bad[:3]}"
            return row
        at_, am_, ar, _ = cav(prot)
        ht_, hm_, hr, hn = cav(prot + lg)
        row.update(apo_cav=round(at_, 1), holo_cav=round(ht_, 1),
                   apo_main=round(am_, 1), holo_main=round(hm_, 1),
                   fill=round(1 - ht_ / at_, 3) if at_ > 0 else "",
                   apo_rmax=round(ar, 2), holo_rmax=round(hr, 2), n_comp=hn)
    except Exception as e:                                       # noqa: BLE001
        row["apo_cav"] = f"ERROR:{type(e).__name__}"
    return row


def calibrate():
    """P1: crystal 3QN1 through the identical code."""
    from importlib import import_module as _im
    m143 = _im("143_pocket_backbone_shape")
    at = [a for a in m143.read_pdb(os.path.join(ROOT, "data", "stage1",
                                                "wt_aba.pdb"), want_atom=False)]
    prot = [(a[3], a[4]) for a in at if a[1] in m143.THREE]
    lg = [(a[3], a[4]) for a in at if a[1] not in m143.THREE]
    at_, a, ar, _ = cav(prot)
    ht_, h, hr, n = cav(prot + lg)
    print(f"P1 CALIBRATION -- crystal PYR1+ABA (wt_aba.pdb), identical code")
    print(f"   apo   total {at_:8.1f}  chamber {a:8.1f} A^3   r_max {ar:.2f}")
    print(f"   holo  total {ht_:8.1f}  chamber {h:8.1f} A^3   r_max {hr:.2f}   comp {n}")
    print(f"   fill {1-ht_/max(at_,1e-9):.3f}   ({len(lg)} ligand heavy atoms)")
    return dict(apo=at_, holo=ht_, apo_main=a, holo_main=h,
                apo_rmax=ar, holo_rmax=hr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=64)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    ref = calibrate()
    json.dump(ref, open(os.path.join(OUT, "crystal_reference.json"), "w"), indent=1)

    key = json.load(open(KEY))
    jobs = []
    for lig_id, meta in sorted(key.items()):
        for p in sorted(glob.glob(os.path.join(
                PRED, f"boltz_results_{lig_id}", "predictions", lig_id, "*.cif"))):
            jobs.append((lig_id, meta, p))
    if a.limit:
        jobs = jobs[:a.limit]
    print(f"\n{len(jobs)} structures ({len(key)} ligands x 5 models), "
          f"{a.procs} procs", flush=True)
    outp = os.path.join(OUT, "predicted_cavity.csv")
    with open(outp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        with mp.Pool(a.procs) as pool:
            for i, row in enumerate(pool.imap_unordered(one, jobs, 4), 1):
                w.writerow(row); fh.flush()
                if i % 200 == 0:
                    print(f"  [{i}/{len(jobs)}]", flush=True)
    print(f"wrote {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
