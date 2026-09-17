#!/usr/bin/env python
r"""
221_cavity_domain_trim.py -- ERRORS 5 AND 6 in the homolog cavity survey.

216 fixed the CHAIN axis (measure every chain, not one arbitrary one) and the
RESIDUE-NAME axis (whitelist, so MSE/CME are not deleted). It left two axes
untouched, and both inflate exactly the tail that carries the "large donor
exists" claim.

 5. THE ALPHAFOLD FILES ARE NOT TRIMMED TO THEIR CATH DOMAIN. The filename
    `af_A0A1D6QKW9_84_319_3.30.530.20` declares the SRPBCC domain to be residues
    84-319, but the file on disk holds residues 1-575. 125 of 195 AlphaFold
    entries (64%) contain more than their named domain, by a median of 2.44x.
    The cavity we reported is therefore the largest chamber ANYWHERE in a
    multi-domain protein, which need not be the SRPBCC pocket at all.
    ⚠ This is the SAME error as the CATH multi-chain bug (3qrzB00 holding chains
    A, B and C). I fixed it on the chain axis and never checked the residue axis,
    and the residue axis is the one that bites the AlphaFold arm.
    It is not cosmetic: 21 of the 25 AlphaFold structures reading >=200 A^3 are
    in this set, including every one of the top seven.

 6. NMR ENSEMBLES ARE READ AS ONE CHAIN. 2LF2 gives 28,600 atoms in "chain A"
    because all ~20 deposited models are stacked on top of each other. Overlaid
    models FILL the pocket, so these read main = 0.0 -- nine false zeros in the
    experimental arm, which is the arm being compared against AlphaFold.

PRE-REGISTERED, before the run:
  * Trimming can only REMOVE volume, never add it, for any structure whose file
    exceeded its domain. If any such structure's main chamber GROWS, the trim is
    wrong and this script must fail rather than report.
  * PYR1 (3QN1 chain A) is not an `af_` entry and has one model, so it must come
    back BIT-IDENTICAL to 216's 119.0 A^3. That is the control.
  * After trimming, no chain may exceed 1.25x its declared domain length in
    residues; assert on the RESULT, not the exit code.
"""
import glob, json, os, re, sys
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m151 = import_module("151_cavity_bottleneck")
m216 = import_module("216_cavity_recompute")
AA_OK = m216.AA_OK
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
OUT = os.path.join(ROOT, "results", "homolog_cavities", "cavity_v3.json")
OBSOLETE = {"3qrz"}

DOM = re.compile(r"^af_[A-Za-z0-9]+_(\d+)_(\d+)_")


def domain_range(target):
    """(lo, hi) CATH domain bounds from an AlphaFold target name, else None."""
    m = DOM.match(target)
    return (int(m.group(1)), int(m.group(2))) if m else None


def read_trimmed(path, rng):
    """{chain: [(element, xyz, resseq)]}; FIRST model only; within rng if given."""
    out = {}

    def keep(resseq):
        return rng is None or (rng[0] <= resseq <= rng[1])

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
        has_model = "pdbx_PDB_model_num" in cols
        first = None
        for r in rows:
            if has_model:                       # ERROR 6: NMR ensembles
                mn = g(r, "pdbx_PDB_model_num")
                if first is None:
                    first = mn
                if mn != first:
                    continue
            if g(r, "type_symbol") in ("H", "D"):
                continue
            if g(r, "label_comp_id") not in AA_OK:
                continue
            try:
                rs = int(g(r, "auth_seq_id"))
            except (ValueError, KeyError):
                continue
            if not keep(rs):
                continue
            out.setdefault(g(r, "auth_asym_id"), []).append(
                (g(r, "type_symbol").upper(),
                 np.array([float(g(r, "Cartn_x")), float(g(r, "Cartn_y")),
                           float(g(r, "Cartn_z"))]), rs))
    else:
        stop = False
        for l in open(path):
            if l.startswith("ENDMDL"):          # ERROR 6, PDB flavour
                stop = True
                continue
            if stop:
                continue
            if not l.startswith(("ATOM", "HETATM")):
                continue
            if l[17:20].strip() not in AA_OK:
                continue
            e = (l[76:78].strip() or l[12:16].strip()[0]).upper()
            if e in ("H", "D"):
                continue
            try:
                rs = int(l[22:26])
            except ValueError:
                continue
            if not keep(rs):
                continue
            out.setdefault(l[21], []).append(
                (e, np.array([float(l[30:38]), float(l[38:46]),
                              float(l[46:54])]), rs))
    return out


def one(path):
    t = os.path.basename(path).rsplit(".", 1)[0]
    rng = domain_range(t)
    rec = dict(target=t, pdb=t[:4].upper(),
               kind="AlphaFold" if t.startswith("af_") else "experimental",
               obsolete=t[:4].lower() in OBSOLETE,
               domain=list(rng) if rng else None)
    try:
        ch = read_trimmed(path, rng)
        per = {}
        for c, at in ch.items():
            if len(at) < 300:
                continue
            nres = len({r for _, _, r in at})
            r = m151.analyse("x", [x for _, x, _ in at],
                             [e for e, _, _ in at]) or {}
            per[c] = dict(main=round(float(r.get("main", 0.0)), 1),
                          total=round(float(r.get("total", 0.0)), 1),
                          r_max=round(float(r.get("r_max", 0.0)), 2),
                          n_atom=len(at), n_res=nres)
        if not per:
            rec["error"] = "no chain with >=300 heavy atoms after trim"
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
    print(f"{len(files)} structures, {procs} procs", flush=True)
    with mp.Pool(procs) as p:
        recs = p.map(one, files)
    ok = [r for r in recs if "error" not in r]
    print(f"{len(ok)} measured, {len(recs)-len(ok)} failed\n")

    old = {r["target"]: r for r in
           json.load(open(os.path.join(ROOT, "results", "homolog_cavities",
                                       "cavity_v2.json")))}

    # -- PRE-REGISTERED CHECK 1: trimming must never GROW a trimmed structure --
    grew = []
    for r in ok:
        o = old.get(r["target"])
        if not o or "main_median" not in o or not r["domain"]:
            continue
        # only meaningful where the file really did exceed its domain
        if max(c["n_res"] for c in r["chains"].values()) < \
           (r["domain"][1] - r["domain"][0] + 1) * 0.95:
            pass
        if r["main_median"] > o["main_median"] + 1.0:
            grew.append((r["target"], o["main_median"], r["main_median"]))

    # -- PRE-REGISTERED CHECK 2: no chain exceeds 1.25x its domain length --
    fat = []
    for r in ok:
        if not r["domain"]:
            continue
        dl = r["domain"][1] - r["domain"][0] + 1
        for c, v in r["chains"].items():
            if v["n_res"] > dl * 1.25:
                fat.append((r["target"], c, v["n_res"], dl))

    # -- PRE-REGISTERED CHECK 3: PYR1 control must be unchanged --
    ctrl = None
    for r in ok:
        if r["pdb"] == "3QN1":
            ctrl = r
    print("CONTROL 3QN1 (must reproduce 119.0 exactly):")
    if ctrl:
        o = old.get(ctrl["target"], {})
        print(f"  v2 {o.get('main_median')}   v3 {ctrl['main_median']}   "
              f"{'OK' if o.get('main_median') == ctrl['main_median'] else 'CHANGED'}")
    else:
        print("  3QN1 NOT IN SET")

    print(f"\nstructures whose cavity GREW after trimming: {len(grew)} (expect 0)")
    for g in sorted(grew, key=lambda x: x[1] - x[2])[:10]:
        print(f"    {g[0][:46]:<46} {g[1]:7.1f} -> {g[2]:7.1f}")
    print(f"chains still exceeding 1.25x their domain: {len(fat)} (expect 0)")
    for f in fat[:10]:
        print(f"    {f[0][:46]:<46} ch{f[1]} {f[2]} res vs domain {f[3]}")

    # -- what actually changed --
    print("\nLARGEST DROPS (v2 -> v3):")
    dl = []
    for r in ok:
        o = old.get(r["target"])
        if o and "main_median" in o:
            dl.append((o["main_median"] - r["main_median"], r["target"],
                       o["main_median"], r["main_median"], r["kind"]))
    for d, t, a, b, k in sorted(dl, reverse=True)[:15]:
        print(f"    {t[:46]:<46} {k[:4]:<4} {a:7.1f} -> {b:7.1f}   -{d:.1f}")

    import statistics as st
    for lbl in ("experimental", "AlphaFold"):
        v2 = [old[r["target"]]["main_median"] for r in ok
              if r["kind"] == lbl and "main_median" in old.get(r["target"], {})]
        v3 = [r["main_median"] for r in ok if r["kind"] == lbl]
        print(f"\n{lbl:<13} n={len(v3):>3}  median {st.median(v2):6.1f} -> "
              f"{st.median(v3):6.1f}   max {max(v2):6.1f} -> {max(v3):6.1f}   "
              f"n>=200: {sum(1 for x in v2 if x >= 200)} -> "
              f"{sum(1 for x in v3 if x >= 200)}")

    json.dump(recs, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    if grew or fat:
        print("\n*** PRE-REGISTERED CHECK FAILED -- do not use cavity_v3 ***")
        return 1
    print("OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
