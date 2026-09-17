#!/usr/bin/env python
r"""
224_cavity_domain_lined.py -- the cavity survey, measured the way it should have
been from the start: IN FULL CONTEXT, but ATTRIBUTED TO THE DOMAIN BY ITS LINING.

WHY v3 (221) WAS ALSO WRONG. 221 trimmed each AlphaFold model to its named CATH
domain, and its own pre-registered check caught the problem: 40 structures got
BIGGER after trimming, which "removing atoms can only remove volume" says is
impossible. The check was right to fire and my premise was wrong. Two mechanisms
were operating at once and they push opposite ways:

  af_Q08058_43_193   v2 whole file   main  82.1   total 266.8
                     v3 domain only  main 297.2   total 327.0

  * In v2 the neighbouring domain's atoms SUBDIVIDE the cavity network, so
    `main` (the largest single chamber) reads 82 while `total` reads 267.
  * In v3 those atoms are gone, the fragments MERGE, and `main` jumps to 297.

But deleting a domain that packs against the fold leaves a DOMAIN-SHAPED HOLE,
and if the remaining structure caps it that hole is indistinguishable from a
ligand pocket. So v3 trades one artefact for another and neither number can be
defended.

THE FIX IS NOT TO CHOOSE BETWEEN THEM. Keep the full structure, so no artificial
interface void is ever created, and then ask of each chamber: IS THIS CHAMBER
LINED BY THE CATH DOMAIN? Report the largest chamber whose lining is mostly
domain. That answers the question the survey is actually asking -- "how big is
the SRPBCC pocket in this protein" -- rather than "how big is the largest hole
anywhere in this file".

⚠ THE NMR FIX FROM 221 IS KEPT. 2LF2 stacks ~20 deposited models into one chain
(28,600 atoms); overlaid models fill the pocket and produce false zeros. First
model only.

⚠ PYR1 IS NOT IN THE SURVEY. `results/homolog_cavities/raw` has no 3QN1 entry,
so the reference "PYR1 = 119 A^3" every homolog has been compared against was
produced by a DIFFERENT script on a different input. This script measures PYR1
through the identical code path, as `3QN1_A` (chain A alone) and `3QN1_AB`
(with HAB1), so the comparison is finally like-for-like.

PRE-REGISTERED:
 * PYR1 chain A must land near 119 A^3 and PYR1+HAB1 near 182 A^3 through this
   path. Those are the two numbers §122 established; if this code cannot
   reproduce them it is measuring something else and the run is void.
 * A chamber is DOMAIN-LINED if >= 2/3 of its lining atoms are domain residues.
   Reported per structure so the threshold can be re-examined, not buried.
 * Assert on the RESULT: every structure reports how many chambers it found and
   how many were domain-lined.
"""
import glob, json, os, re, sys
import multiprocessing as mp
import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m151 = import_module("151_cavity_bottleneck")
m216 = import_module("216_cavity_recompute")
m221 = import_module("221_cavity_domain_trim")
AA_OK = m216.AA_OK
SPACING, PROBE = m151.SPACING, m151.PROBE
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
OUT = os.path.join(ROOT, "results", "homolog_cavities", "cavity_v4.json")
LINING_CUT = 5.0          # A from a cavity voxel to count an atom as lining
DOMAIN_FRAC = 2.0 / 3.0


def chambers(xyz, elem, resseq):
    """[(volume, domain_fraction, r_max)] for every 2.4 A chamber, largest first.

    `resseq` is a per-atom flag: True if the atom belongs to the CATH domain.
    """
    clear, _, lo, shape = m151.clear_map(xyz, elem)
    mask = m151.enclosed(clear, lo, shape)
    if mask.sum() == 0:
        return [], 0.0
    tot = float(mask.sum() * SPACING ** 3)
    core = mask & (clear > 2.4)
    lab, n = ndimage.label(core)
    if n == 0:
        return [], tot
    sizes = ndimage.sum(core, lab, range(1, n + 1))
    order = np.argsort(sizes)[::-1][:8]
    seeds = []
    for k in order:
        idx = np.argwhere(lab == k + 1)
        seeds.append(idx[np.argmax(clear[tuple(idx.T)])])
    cvox = np.argwhere(mask)
    S = np.array(seeds)
    own = np.linalg.norm(cvox[:, None] - S[None], axis=-1).argmin(1)
    X = np.asarray(xyz, float)
    dom = np.asarray(resseq, bool)
    out = []
    for j in range(len(seeds)):
        sel = cvox[own == j]
        if len(sel) == 0:
            continue
        vol = float(len(sel) * SPACING ** 3)
        pts = lo + sel * SPACING
        # lining atoms: any atom within LINING_CUT of any voxel of this chamber
        step = max(1, len(pts) // 4000)
        p = pts[::step]
        d = np.linalg.norm(X[None, :, :] - p[:, None, :], axis=2)
        near = (d < LINING_CUT).any(0)
        nl = int(near.sum())
        frac = float(dom[near].sum() / nl) if nl else 0.0
        rm = float(clear[tuple(sel.T)].max())
        out.append((round(vol, 1), round(frac, 3), round(rm, 2), nl))
    out.sort(key=lambda r: -r[0])
    return out, tot


def one(path):
    t = os.path.basename(path).rsplit(".", 1)[0]
    rng = m221.domain_range(t)
    rec = dict(target=t, pdb=t[:4].upper(),
               kind="AlphaFold" if t.startswith("af_") else "experimental",
               obsolete=t[:4].lower() in m221.OBSOLETE,
               domain=list(rng) if rng else None)
    try:
        ch = m221.read_trimmed(path, None)     # FULL context, first model only
        per = {}
        for c, at in ch.items():
            if len(at) < 300:
                continue
            xyz = [x for _, x, _ in at]
            el = [e for e, _, _ in at]
            if rng is None:
                dom = [True] * len(at)         # single-domain file
            else:
                dom = [rng[0] <= r <= rng[1] for _, _, r in at]
            if not any(dom):
                continue
            cs, tot = chambers(xyz, el, dom)
            dl = [c_ for c_ in cs if c_[1] >= DOMAIN_FRAC]
            per[c] = dict(
                total=round(tot, 1), n_chamber=len(cs), n_domain_lined=len(dl),
                main=dl[0][0] if dl else 0.0,
                main_frac=dl[0][1] if dl else None,
                r_max=dl[0][2] if dl else 0.0,
                largest_any=cs[0][0] if cs else 0.0,
                largest_any_frac=cs[0][1] if cs else None,
                n_atom=len(at), n_res_domain=int(sum(dom)))
        if not per:
            rec["error"] = "no usable chain"
            return rec
        mains = [v["main"] for v in per.values()]
        rec.update(chains=per, n_chain=len(per),
                   main_median=float(np.median(mains)),
                   main_min=min(mains), main_max=max(mains),
                   spread=round(max(mains) - min(mains), 1))
    except Exception as e:                                       # noqa: BLE001
        rec["error"] = f"{type(e).__name__}: {e}"
    return rec


def pyr1_controls():
    """PYR1 through the IDENTICAL code path: chain A alone, and with HAB1."""
    import gzip
    cif = os.path.join(ROOT, "data", "3QN1.cif")
    if not os.path.exists(cif):
        for c in glob.glob(os.path.join(ROOT, "data", "*3QN1*")) + \
                 glob.glob(os.path.join(ROOT, "data", "*3qn1*")):
            cif = c; break
    if not os.path.exists(cif):
        return [{"target": "3QN1", "error": "3QN1 not found under data/"}]
    out = []
    for lbl, keep in (("3QN1_A", {"A"}), ("3QN1_AB", None)):
        ch = m221.read_trimmed(cif, None)
        at = [a for c, v in ch.items() if (keep is None or c in keep) for a in v]
        if len(at) < 300:
            out.append({"target": lbl, "error": f"only {len(at)} atoms"}); continue
        cs, tot = chambers([x for _, x, _ in at], [e for e, _, _ in at],
                           [True] * len(at))
        out.append(dict(target=lbl, kind="control", total=round(tot, 1),
                        n_chamber=len(cs),
                        main=cs[0][0] if cs else 0.0,
                        r_max=cs[0][2] if cs else 0.0,
                        n_atom=len(at)))
    return out


def main():
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 48
    print("PYR1 CONTROLS through this exact code path "
          "(expect ~119 chain A, ~182 with HAB1):", flush=True)
    ctrl = pyr1_controls()
    for c in ctrl:
        if "error" in c:
            print(f"   {c['target']}: ERROR {c['error']}")
        else:
            print(f"   {c['target']:<9} main {c['main']:7.1f}  total {c['total']:7.1f}"
                  f"  chambers {c['n_chamber']}  atoms {c['n_atom']}")
    files = sorted(glob.glob(f"{RAW}/*.cif") + glob.glob(f"{RAW}/*.pdb"))
    print(f"\n{len(files)} structures, {procs} procs", flush=True)
    with mp.Pool(procs) as p:
        recs = p.map(one, files)
    ok = [r for r in recs if "error" not in r]
    print(f"{len(ok)} measured, {len(recs)-len(ok)} failed")

    # ⚠ PERSIST BEFORE REPORTING. v4's first run measured all 261 structures in
    # 66 minutes and then died in the summary printer on `sorted()` comparing two
    # dicts after a tie -- and because the dump came last, every measurement was
    # lost. Nothing expensive may sit upstream of a print statement again.
    json.dump(dict(controls=ctrl, structures=recs), open(OUT, "w"), indent=1)
    print(f"wrote {OUT}\n")

    v2 = {r["target"]: r for r in json.load(
        open(os.path.join(ROOT, "results", "homolog_cavities", "cavity_v2.json")))}
    import statistics as st
    print("\nHOW MUCH OF THE OLD SIGNAL WAS THE WRONG CHAMBER:")
    print(f"{'target':<46}{'v2':>8}{'v4':>8}{'domfrac':>9}{'largest_any':>12}")
    dl = []
    for r in ok:
        o = v2.get(r["target"], {}).get("main_median")
        if o is None:
            continue
        dl.append((o - r["main_median"], r["target"], r))
    for d, _t, r in sorted(dl, key=lambda x: (-x[0], x[1]))[:12]:
        c = list(r["chains"].values())[0]
        print(f"   {r['target'][:43]:<43}{v2[r['target']]['main_median']:>8.1f}"
              f"{r['main_median']:>8.1f}{str(c['main_frac']):>9}"
              f"{c['largest_any']:>12.1f}")
    for lbl in ("experimental", "AlphaFold"):
        a = [v2[r["target"]]["main_median"] for r in ok if r["kind"] == lbl
             and r["target"] in v2]
        b = [r["main_median"] for r in ok if r["kind"] == lbl]
        nd = sum(list(r["chains"].values())[0]["n_domain_lined"] == 0
                 for r in ok if r["kind"] == lbl)
        print(f"\n{lbl:<13} n={len(b):>3}  median {st.median(a):6.1f} -> "
              f"{st.median(b):6.1f}   max {max(a):6.1f} -> {max(b):6.1f}   "
              f"n>=200: {sum(1 for x in a if x>=200)} -> "
              f"{sum(1 for x in b if x>=200)}   "
              f"NO domain-lined chamber: {nd}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
