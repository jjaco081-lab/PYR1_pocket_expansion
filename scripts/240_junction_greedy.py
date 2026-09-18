#!/usr/bin/env python
r"""
240_junction_greedy.py -- place junctions at the NEAREST acceptable position,
not the globally best one.

⚠ THE BUG THIS FIXES. `214.scan_junctions` slides each boundary through its whole
gap and takes the position of least CA deviation. That minimises junction cost
and nothing else: every residue between the required boundary and the chosen
junction is handed to PYR1, with no term rewarding donor content. Measured:

    PYR1 keep-set (flank 2)        70 of 179 residues = 39 %
    FREE for a donor              109 of 179          = 61 %
    what the splice actually gave  126 to PYR1, 19 to the donor

so ~56 free residues were spent buying junction quality, and only 8 % of the
donor's own pocket lining survived. Jannis: "The majority of PYR1 is not asked to
be conserved, right?" Correct -- 61 % is free, and the optimiser was giving it
away.

THE FIX: walk outward from each required boundary and stop at the FIRST position
whose deviation is below `--cut` (default 1.5 A, the same threshold "clean"
already uses). That keeps the junction good enough while leaving the donor
everything past it. Only if no position in the gap clears the cut does it fall
back to the gap's minimum, which is 214's behaviour.

⚠ This trades junction cost for donor content deliberately, so BOTH are reported:
a chimera that keeps more donor but has a 3 A junction is not better, it is a
different trade. The comparison against 214's placement is printed side by side.
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
imp = __import__("importlib").import_module
m214, m237 = imp("214_schema_foldseek"), imp("237_junctions_all")
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
OUT = os.path.join(ROOT, "results", "schema_chimera")


def scan_greedy(ks, K, dev, cut):
    """as 214.scan_junctions but stopping at the FIRST acceptable position."""
    idx = {k: i for i, k in enumerate(ks)}
    runs, cur = [], []
    for k in ks:
        if k in K:
            cur.append(k)
        elif cur:
            runs.append(cur); cur = []
    if cur:
        runs.append(cur)
    parent = {k: ("PYR1" if k in K else "DON") for k in ks}
    J = []
    for run in runs:
        for end, step in ((run[0], -1), (run[-1], +1)):
            i = idx[end]
            best, bestd = end, dev.get(end, 9.9)
            chosen = None
            off, fallback, fbd = 1, end, bestd
            while True:
                j = i + step * off
                if not (0 <= j < len(ks)) or ks[j] in K:
                    break
                dj = dev.get(ks[j], 9.9)
                if dj < fbd:
                    fallback, fbd = ks[j], dj      # 214's global minimum
                if chosen is None and dj < cut:
                    chosen = ks[j]                 # FIRST acceptable -- stop here
                off += 1
            if chosen is None:
                chosen, cd = fallback, fbd         # nothing clears the cut
            else:
                cd = dev[chosen]
            lo, hi = sorted((idx[end], idx[chosen]))
            for t in range(lo, hi + 1):
                parent[ks[t]] = "PYR1"
            J.append((chosen, cd))
    return parent, J


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+", required=True)
    ap.add_argument("--flank", type=int, default=2)
    ap.add_argument("--cut", type=float, default=1.5)
    a = ap.parse_args()
    pca, ph, ps = m214.read_cif(os.path.join(ROOT, "data", "3QN1.cif"), "A")
    pres = sorted(pca)
    aln = {}
    for line in open(f"{ROOT}/results/foldseek/fs_hits.tsv"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 11 and f[1] not in aln:
            aln[f[1]] = f
    K0 = m214.keep_set(a.flank)
    print(f"PYR1 {len(pres)} residues; keep-set {len(K0 & set(pres))} "
          f"({100*len(K0 & set(pres))/len(pres):.0f} %); junction cut {a.cut} A\n")
    print(f"{'target':<12}{'':>4}{'PYR1':>6}{'donor':>7}{'clean':>7}{'meanJ':>7}"
          f"{'maxJ':>7}")
    rows = []
    for t in a.targets:
        path = None
        for ext in (".cif", ".pdb"):
            p = os.path.join(RAW, t + ext)
            if os.path.exists(p): path = p; break
        if path is None or t not in aln:
            print(f"!! {t}: missing"); continue
        dchain = (t[4] if (path.endswith(".cif") and len(t) > 4) else None)
        dca, dh, dsq = m237.read_any(path, dchain)
        dres = sorted(dca)
        f = aln[t]
        m = m214.map_from_alignment(f[9], f[10], int(f[5]), int(f[7]), pres, dres)
        m = {k: v for k, v in m.items() if k in pca and v in dca}
        ks = sorted(m)
        P = np.array([pca[k] for k in ks]); Q = np.array([dca[m[k]] for k in ks])
        R, pm, qm = m214.kabsch(P, Q)
        dev = {k: float(np.linalg.norm(((pca[k] - pm) @ R + qm) - dca[m[k]]))
               for k in ks}
        K = K0 & set(ks)
        for lbl, fn in (("214 (min)", lambda: m214.scan_junctions(ks, K, dev)),
                        ("greedy", lambda: scan_greedy(ks, K, dev, a.cut))):
            parent, J = fn()
            npy = sum(1 for k in ks if parent[k] == "PYR1")
            nd = len(ks) - npy
            jd = [d for _, d in J]
            print(f"{t[:11]:<12}{lbl:>10}{npy:>6}{nd:>7}"
                  f"{sum(1 for d in jd if d < 1.5):>4}/{len(J):<2}"
                  f"{np.mean(jd):>7.2f}{max(jd):>7.2f}")
            rows.append(dict(target=t, mode=lbl, n_pyr1=npy, n_donor=nd,
                             clean=sum(1 for d in jd if d < 1.5),
                             n_junction=len(J),
                             mean_junction=round(float(np.mean(jd)), 2),
                             max_junction=round(float(max(jd)), 2)))
        print()
    json.dump(rows, open(os.path.join(OUT, "junction_greedy.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
