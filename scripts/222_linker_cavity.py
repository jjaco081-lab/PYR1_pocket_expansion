#!/usr/bin/env python
r"""
222_linker_cavity.py -- DOES LINKER LENGTH BUY POCKET VOLUME?

THE DECISION THIS GATES. The gap analysis showed internal appendages >=5 res
concentrate in gap3 (26 of 73, median 12) and gap4 (17, median 16) -- the two
gaps where our mandated minimum EXCEEDS PYR1's native linker:

    gap  motif span        native len   mandated      excess
    1    A40  -> A58          17         15-30          ok
    2    A65  -> A81          15         15-30          ok
    3    A92  -> A111         18         20-35        +2 minimum
    4    A121 -> A146         24         25-40        +1 minimum

Jannis: "I want to have significant size increases in the pocket still if
possible so I am hesitant to [shorten linkers]." That hesitation is a testable
claim, and it has never been tested. Two rival models:

  MODEL A (Jannis's concern): extra linker residues form part of the wall that
    encloses a larger chamber. Shortening them shrinks the pocket. Prediction:
    Spearman(gap length, cavity) > 0, surviving control for design size.

  MODEL B (my expectation): the pocket is defined by the MOTIF, which is fixed
    in all arms. Extra linker residues have nowhere useful to go and leave as
    surface appendages. Prediction: rho ~ 0 against cavity, and a POSITIVE
    correlation with appendage length instead.

They make opposite predictions on the same numbers, so this is decidable.

PRE-REGISTERED, before looking:
 * Primary: partial Spearman(realised gap length, cavity) controlling n_main,
   per gap, BH-FDR over the 4 gaps. Model A needs rho > 0 at q < 0.05.
 * Secondary: Spearman(gap length, appendage length) -- Model B's signature.
 * ⚠ CONTROL FOR DESIGN SIZE. A longer design has more residues everywhere and
   a bigger envelope; the raw correlation would be confounded. n_main is the
   control variable.
 * ⚠ ONLY INTACT DESIGNS. A gap length measured across a chain break is the
   distance between two pieces, not a linker. Require n_pieces == 1 and
   motif_split == 0, and report how many designs that discards.
 * Cavity is the MAIN-CHAMBER value on the largest connected piece (cascade v2),
   not total enclosed volume -- the §122 error.

If Model B wins, gaps 3 and 4 can be shortened to native length without costing
pocket volume, and the appendage problem goes away with them.
"""
import csv, glob, gzip, json, os, sys
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m195 = import_module("195_arm_compare")
MOTIF = m195.MOTIF
CSV = os.path.join(ROOT, "results", "rfd3", "cascade2.csv")
OUT = os.path.join(ROOT, "results", "rfd3", "linker_cavity.csv")

#: native PYR1 linker lengths, from the motif definition itself
NATIVE = [MOTIF[i + 1][1] - MOTIF[i][2] - 1 for i in range(len(MOTIF) - 1)]


def gaps_of(cifp):
    """realised residue count in each of the 4 inter-motif gaps, or None."""
    js = cifp.replace(".cif.gz", ".json").replace(".cif", ".json")
    if not os.path.exists(js):
        return None
    try:
        imap = json.load(open(js)).get("diffused_index_map", {})
        at = m195.read_cif_gz(cifp)
        dc = m195.design_chain(at, imap)
    except Exception:                                            # noqa: BLE001
        return None
    idx = {}
    for c, lo, hi, _ in MOTIF:
        for r in (lo, hi):
            k = f"{c}{r}"
            if k not in imap or imap[k][0] != dc:
                return None
            idx[r] = int(imap[k][1:])
    out = []
    for i in range(len(MOTIF) - 1):
        end_prev = idx[MOTIF[i][2]]
        start_next = idx[MOTIF[i + 1][1]]
        n = start_next - end_prev - 1
        if n < 0:
            return None
        out.append(n)
    return out


def rank(x):
    o = np.argsort(np.argsort(np.asarray(x, float)))
    return o.astype(float)


def spearman(x, y):
    rx, ry = rank(x), rank(y)
    rx = rx - rx.mean(); ry = ry - ry.mean()
    d = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / d) if d else 0.0


def partial_spearman(x, y, z):
    """Spearman(x,y) with z partialled out, on ranks."""
    rx, ry, rz = rank(x), rank(y), rank(z)
    def resid(a, b):
        b1 = np.c_[b - b.mean(), np.ones(len(b))]
        coef, *_ = np.linalg.lstsq(b1, a - a.mean(), rcond=None)
        return (a - a.mean()) - b1 @ coef
    ex, ey = resid(rx, rz), resid(ry, rz)
    d = np.sqrt((ex * ex).sum() * (ey * ey).sum())
    return float((ex * ey).sum() / d) if d else 0.0


def perm_p(x, y, z, obs, n=20000, seed=0):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    c = 0
    for _ in range(n):
        if abs(partial_spearman(rng.permutation(x), y, z)) >= abs(obs):
            c += 1
    return (c + 1) / (n + 1)


def bh(ps):
    m = len(ps); order = np.argsort(ps); q = np.empty(m)
    prev = 1.0
    for rank_i, i in enumerate(reversed(order)):
        k = m - rank_i
        prev = min(prev, ps[i] * m / k)
        q[i] = prev
    return q


def main():
    rows = {}
    with open(CSV) as fh:
        for r in csv.DictReader(fh):
            rows[(r["arm"], r["design"])] = r
    jobs = []
    for (arm, d) in rows:
        for pat in (f"{ROOT}/results/rfd3/{arm}/{d}/*.cif.gz",):
            g = glob.glob(pat)
            if g:
                jobs.append((arm, d, g[0]))
                break
    print(f"{len(rows)} cascade rows, {len(jobs)} with a structure", flush=True)
    with mp.Pool(48) as p:
        res = p.map(gaps_of, [j[2] for j in jobs])

    recs = []
    for (arm, d, _), g in zip(jobs, res):
        if g is None:
            continue
        r = rows[(arm, d)]
        if r["n_pieces"] != "1" or r["motif_split"] != "0":
            continue
        try:
            cav = float(r["cavity"] or 0)
            nmain = float(r["n_main"] or 0)
            app = float(r["appendage"] or 0)
        except ValueError:
            continue
        recs.append(dict(arm=arm, design=d, cavity=cav, n_main=nmain,
                         appendage=app,
                         **{f"gap{i+1}": g[i] for i in range(4)}))
    n_disc = len(jobs) - len(recs)
    print(f"INTACT designs used: {len(recs)}   discarded (broken/split/unreadable): {n_disc}")
    if len(recs) < 30:
        print("too few designs"); return 1

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0]))
        w.writeheader(); w.writerows(recs)

    cav = [r["cavity"] for r in recs]
    nm = [r["n_main"] for r in recs]
    app = [r["appendage"] for r in recs]

    print(f"\nnative linker lengths from MOTIF: {NATIVE}")
    print(f"cavity: median {np.median(cav):.1f}  mean {np.mean(cav):.1f}  "
          f"n>0 {sum(1 for c in cav if c > 0)}/{len(cav)}")

    print("\n" + "=" * 74)
    print("PRIMARY  partial Spearman(gap length, CAVITY | n_main)")
    print("=" * 74)
    print(f"{'gap':<6}{'native':>7}{'median':>8}{'range':>12}"
          f"{'rho_partial':>13}{'p_perm':>9}{'q_BH':>8}")
    obs, ps, meta = [], [], []
    for i in range(4):
        g = [r[f"gap{i+1}"] for r in recs]
        o = partial_spearman(g, cav, nm)
        p = perm_p(g, cav, nm, o, n=4000, seed=i)
        obs.append(o); ps.append(p)
        meta.append((i + 1, NATIVE[i], np.median(g), min(g), max(g)))
    q = bh(ps)
    for (gi, nat, med, lo, hi), o, p, qq in zip(meta, obs, ps, q):
        star = "  <-- MODEL A" if (o > 0 and qq < 0.05) else ""
        print(f"gap{gi:<3}{nat:>7}{med:>8.0f}{f'{lo:.0f}-{hi:.0f}':>12}"
              f"{o:>13.3f}{p:>9.4f}{qq:>8.3f}{star}")

    print("\n" + "=" * 74)
    print("SECONDARY  Spearman(gap length, APPENDAGE length)   [Model B signature]")
    print("=" * 74)
    for i in range(4):
        g = [r[f"gap{i+1}"] for r in recs]
        print(f"gap{i+1}   rho = {spearman(g, app):+.3f}")
    tot = [sum(r[f"gap{i+1}"] for i in range(4)) for r in recs]
    print(f"\nTOTAL linker vs cavity   partial rho = "
          f"{partial_spearman(tot, cav, nm):+.3f}")
    print(f"TOTAL linker vs appendage        rho = {spearman(tot, app):+.3f}")
    excess = [sum(max(0, r[f"gap{i+1}"] - NATIVE[i]) for i in range(4))
              for r in recs]
    print(f"EXCESS over native vs cavity     partial rho = "
          f"{partial_spearman(excess, cav, nm):+.3f}")
    print(f"EXCESS over native vs appendage          rho = "
          f"{spearman(excess, app):+.3f}")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
