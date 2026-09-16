#!/usr/bin/env python
r"""
216_cavity_recompute.py -- recompute every homolog cavity as the USABLE MAIN
CHAMBER, over EVERY chain, with junction cost, for experimental and AlphaFold
structures alike.

WHY EVERYTHING BEFORE THIS IS PROVISIONAL. Three errors were found by following
Jannis's question about why 3QRZ's chains disagreed 2x:

 1. `cavity_A3` in the old survey is TOTAL ENCLOSED VOLUME, which sums voids
    separated by necks a ligand cannot pass. Ara h 8 reads 319 A^3 total but
    only 126 A^3 of usable chamber -- essentially PYR1's 119. This is the same
    error Jannis caught in the RFd3 cascade ("only count the largest pocket");
    it was fixed there in 190/198 and never propagated back here.
 2. ONE ARBITRARY CHAIN was measured. 4JDL's three copies give main chambers of
    128, 16 and 13 A^3 -- an 8x spread WITHIN ONE CRYSTAL, with chains A and B
    equally complete (14 vs 13 missing residues). Whichever chain CATH picked
    became the published number.
 3. An OBSOLETE entry is in the set: 3QRZ was superseded by 4JDL on 2013-03-13
    (Jannis). Both are present, giving three different answers for PYL5.

And PYR1's own reference is the apo MONOMER: chain A alone gives 164 total /
119 main, but chain A + HAB1 gives 182 as a SINGLE chamber. HAB1 does not line
the pocket -- it SEALS it, merging the satellites. So the functional closed-state
pocket is 182 A^3 and we have been comparing homologs against 119.

WHAT THIS WRITES, per structure and per chain:
    main_chamber   usable volume (largest 2.4 A-sphere-accessible chamber)
    total          all enclosed volume, for comparison with the old column
    spread         max-min main chamber across chains = a free error bar
⚠ AlphaFold models are single-chain by construction, so they have NO spread.
That is a limitation of the prediction, not a property of the protein, and it
must not be read as "more reproducible".
⚠ CATH domain files carry UNK residue names for BOTH experimental and predicted
entries, but both retain side chains (CB/CG/CD present), so cavities are
side-chain-lined in both cases.
"""
import csv, glob, json, os, sys
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m151 = import_module("151_cavity_bottleneck")
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
OBSOLETE = {"3qrz"}          # superseded by 4JDL on 2013-03-13


#: WHITELIST of residue names counted as protein. A whitelist, not a ligand
#: blacklist, because a blacklist fails silently on the next unfamiliar HETATM
#: code -- an unrecognised residue is then excluded and shows up as an
#: anomalously LARGE cavity, which is visible, rather than a ligand being
#: counted as protein, which is not.
#:
#: ⚠ MSE AND CME ARE HETATM BUT ARE PART OF THE CHAIN. Selenomethionine appears
#: in 1,554 atoms across 23 of the 66 experimental structures (35 %), and
#: dropping it leaves a methionine-shaped hole the cavity code reads as pocket.
#: AlphaFold models contain NO HETATM at all, so this error would have inflated
#: EXPERIMENTAL cavities relative to predicted ones -- precisely the comparison
#: the figure makes.
#:
#: ⚠ A8S (ABA) is present in the homolog set. Every structure is stripped to
#: APO, because a ligand left in place turns the measurement into leftover
#: volume rather than pocket volume, and only some homologs are ligand-bound.
#: CATH domain files also carry UNK for real residues, so UNK is kept.
AA_OK = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "UNK",                                    # CATH domain extraction artifact
    "MSE", "CME", "SEP", "TPO", "PTR", "CSO", "KCX", "LLP", "MLY", "OCS",
    "CSD", "PCA", "SAC", "HIC", "M3L", "CSX", "CAS", "SNC", "MHO", "ALY",
}


def read_any(path, chain=None):
    """returns {chain: [(element, xyz)]} for .cif or .pdb, protein only, no H"""
    out = {}
    if path.endswith(".cif"):
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
        g = lambda r, k: r[cols[k]]                              # noqa: E731
        for r in rows:
            if g(r, "type_symbol") in ("H", "D"):
                continue
            if g(r, "label_comp_id") not in AA_OK:      # whitelist, not blacklist
                continue
            out.setdefault(g(r, "auth_asym_id"), []).append(
                (g(r, "type_symbol").upper(),
                 np.array([float(g(r, "Cartn_x")), float(g(r, "Cartn_y")),
                           float(g(r, "Cartn_z"))])))
    else:
        for l in open(path):
            if not l.startswith(("ATOM", "HETATM")):
                continue
            if l[17:20].strip() not in AA_OK:
                continue
            e = (l[76:78].strip() or l[12:16].strip()[0]).upper()
            if e in ("H", "D"):
                continue
            out.setdefault(l[21], []).append(
                (e, np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out


def one(path):
    t = os.path.basename(path).rsplit(".", 1)[0]
    rec = dict(target=t, pdb=t[:4].upper(),
               kind="AlphaFold" if t.startswith("af_") else "experimental",
               obsolete=t[:4].lower() in OBSOLETE)
    try:
        ch = read_any(path)
        per = {}
        for c, at in ch.items():
            if len(at) < 300:
                continue
            r = m151.analyse("x", [x for _, x in at], [e for e, _ in at]) or {}
            per[c] = dict(main=round(float(r.get("main", 0.0)), 1),
                          total=round(float(r.get("total", 0.0)), 1),
                          r_max=round(float(r.get("r_max", 0.0)), 2),
                          n_atom=len(at))
        if not per:
            rec["error"] = "no chain with >=300 heavy atoms"
            return rec
        mains = [v["main"] for v in per.values()]
        rec.update(chains=per, n_chain=len(per),
                   main_median=float(np.median(mains)),
                   main_min=min(mains), main_max=max(mains),
                   spread=round(max(mains) - min(mains), 1),
                   total_median=float(np.median([v["total"] for v in per.values()])))
    except Exception as e:                                       # noqa: BLE001
        rec["error"] = f"{type(e).__name__}: {e}"
    return rec


def main():
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 48
    files = sorted(glob.glob(f"{RAW}/*.cif") + glob.glob(f"{RAW}/*.pdb"))
    print(f"{len(files)} structures "
          f"({sum(1 for f in files if os.path.basename(f).startswith('af_'))} AlphaFold), "
          f"{procs} procs", flush=True)
    with mp.Pool(procs) as pool:
        rows = pool.map(one, files)
    ok = [r for r in rows if "error" not in r]
    print(f"{len(ok)} measured, {len(rows)-len(ok)} failed")
    multi = [r for r in ok if r["n_chain"] > 1]
    if multi:
        sp = [r["spread"] for r in multi]
        print(f"\nWITHIN-CRYSTAL SPREAD ({len(multi)} multi-chain structures):")
        print(f"  median {np.median(sp):.1f} Å³, max {max(sp):.1f} Å³")
        print(f"  structures whose chains disagree by >50 Å³: "
              f"{sum(1 for s in sp if s > 50)}")
        for r in sorted(multi, key=lambda x: -x["spread"])[:6]:
            per = ", ".join(f"{c}:{v['main']:.0f}" for c, v in r["chains"].items())
            print(f"    {r['pdb']}  spread {r['spread']:>6.1f}   chains {per}")
    obs = [r for r in ok if r["obsolete"]]
    print(f"\nOBSOLETE entries present: {[r['pdb'] for r in obs]}")
    out = os.path.join(ROOT, "results", "homolog_cavities", "cavity_v2.json")
    json.dump(rows, open(out, "w"), indent=1)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
