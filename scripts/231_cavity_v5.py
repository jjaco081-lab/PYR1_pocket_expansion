r"""
231_cavity_v5.py -- the chamber probe radius, CALIBRATED against the ligand
instead of asserted.

⚠ THIS PARTLY REVERSES §122. I switched the survey from `total` to `main` because
total sums voids a ligand cannot pass between. That was right in principle and
the threshold was wrong: the 2.4 A sphere was my own invention and nothing ever
checked it against a known answer.

Jannis checked it by eye in PyMOL (`surface_cavity_mode 2`, radius -3, cutoff -5)
and reported that PyMOL finds ONE pocket where v4 reports three, and that ABA
pokes outside the measured surface. Both are true:

    ABA heavy atoms outside the v4 main chamber:  6 of 19
    chambers A/B/C = 119.0 / 26.8 / 18.6 = 164.4 = the `total` v4 discarded
    ABA occupies ALL THREE

The necks between those pieces are narrower than 2.4 A, but ABA threads through
them, so they are not separations. Sweeping the probe against the known answer:

    probe R   n chambers   main A^3   ABA spans
       2.4        3          119.0        3
       2.0        5          122.0        5
       1.8        2          134.8        2
       1.4        1          164.4        1   <-- ABA in ONE chamber
       1.2        1          164.4        1

**1.4 A -- the water probe -- is the smallest radius at which the cognate ligand
lies in a single chamber, and it is also what PyMOL, CASTp and fpocket use.**
At 1.4 A, 18 of 19 ABA atoms sit inside the enclosed volume; the one exception is
O10, the KETONE oxygen, 1.07 A out, which is the most solvent-exposed atom in the
ligand (it coordinates the W385 gate water) and is excluded by the burial cutoff
rather than by the chamber criterion.

PRE-REGISTERED: PYR1 chain A must come back as ONE chamber of ~164 A^3, and
PYR1+HAB1 as one chamber of ~182. If chain A still fragments, the probe change
did not take effect and the run is void.

⚠ WHAT THIS DOES NOT LICENSE. Merging at 1.4 A does NOT mean `total` was right
all along: 1.4 A still separates genuinely disconnected voids, which is the §122
error. It only stops separating voids a ligand can thread. The domain-LINING
attribution from v4 is unchanged and still does the work of rejecting chambers
that belong to a different domain.
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
OUT = os.path.join(ROOT, "results", "homolog_cavities", "cavity_v5.json")
#: Chamber-separation probe, CALIBRATED (see header): the smallest radius at
#: which ABA lies in a single chamber of PYR1. 2.4 A was asserted, never checked.
CHAMBER_PROBE = 1.4
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
    core = mask & (clear > CHAMBER_PROBE)
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
          "(expect ~164 chain A as ONE chamber, ~182 with HAB1):", flush=True)
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
        open(os.path.join(ROOT, "results", "homolog_cavities", "cavity_v4.json")))}
    import statistics as st
    print("\nv4 (2.4 A probe) -> v5 (1.4 A water probe):")
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
