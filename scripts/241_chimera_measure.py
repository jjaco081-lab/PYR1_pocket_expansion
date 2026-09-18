#!/usr/bin/env python
r"""
241_chimera_measure.py -- build the chimera keeping DONOR INSERTIONS, then
measure its cavity.

⚠ THE CAP I IMPOSED WITHOUT NOTICING. Jannis: "The total number of amino acids
does not have a cap, right?" It should not, and it did. 214's junction scan works
on `ks = sorted(m)`, the ALIGNED positions only, and 239 wrote only residues in
that set -- so every donor residue with no PYR1 partner was silently dropped.
Those insertions are exactly where a bigger pocket would come from, and a real
chimera has no reason to be PYR1's length.

WHAT THIS BUILDS. Junctions are placed on the aligned frame as before (greedy:
first position under `--cut`), but each donor SEGMENT is then taken as the full
contiguous donor residue range between its two flanking junctions -- insertions
included. PYR1 contributes its kept residues unchanged. The result can be longer
or shorter than PYR1 and that is the point.

THEN IT IS MEASURED, with the same calibrated code as everything else (chamber
probe 1.4 A, validated so 18 of 19 ABA atoms sit inside PYR1's chamber), so the
number is directly comparable to PYR1's 164.4 and to the RFd3 designs.

⚠ WHAT THIS IS NOT. Still a splice map: no loop closure, no repacking, no
minimisation, no sequence design. A cavity measured on it is an UPPER bound on
what the backbone offers, not a prediction -- side chains from the two parents
have never been packed against each other. Treat a collapse as informative and a
large value as permissive. Reported alongside the composition so the two are
never read apart.
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
imp = __import__("importlib").import_module
m214, m229, m237, m240 = (imp("214_schema_foldseek"), imp("229_cavity_visualise"),
                          imp("237_junctions_all"), imp("240_junction_greedy"))
m239 = imp("239_build_chimera")
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
OUT = os.path.join(ROOT, "results", "chimera")
PYR1_REF = 164.4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+", required=True)
    ap.add_argument("--flank", type=int, default=2)
    ap.add_argument("--cut", type=float, default=1.5)
    ap.add_argument("--no-insertions", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    pca, ph, ps = m214.read_cif(os.path.join(ROOT, "data", "3QN1.cif"), "A")
    pres = sorted(pca)
    aln = {}
    for line in open(f"{ROOT}/results/foldseek/fs_hits.tsv"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 11 and f[1] not in aln:
            aln[f[1]] = f
    K0 = m214.keep_set(a.flank)
    print(f"PYR1: {len(pres)} residues, cavity {PYR1_REF} A^3, "
          f"keep-set {len(K0 & set(pres))} ({100*len(K0&set(pres))//len(pres)} %)\n")
    print(f"{'donor':<12}{'len':>5}{'PYR1':>6}{'don':>5}{'ins':>5}"
          f"{'cavity':>8}{'xPYR1':>7}{'clean':>7}")
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
        gx = lambda Y: (Y - qm) @ R.T + pm                       # noqa: E731
        dev = {k: float(np.linalg.norm(((pca[k] - pm) @ R + qm) - dca[m[k]]))
               for k in ks}
        K = K0 & set(ks)
        parent, J = m240.scan_greedy(ks, K, dev, a.cut)

        keep_p = sorted(k for k in ks if parent[k] == "PYR1")
        # donor segments: contiguous runs of DON in the aligned frame, expanded
        # to the FULL donor range between their flanking aligned positions
        segs, cur = [], []
        for k in ks:
            if parent[k] == "DON":
                cur.append(k)
            elif cur:
                segs.append(cur); cur = []
        if cur:
            segs.append(cur)
        keep_d, n_ins = [], 0
        for seg in segs:
            lo, hi = m[seg[0]], m[seg[-1]]
            if a.no_insertions:
                rng = [m[k] for k in seg]
            else:
                rng = [r for r in dres if lo <= r <= hi]   # INSERTIONS INCLUDED
                n_ins += len(rng) - len(seg)
            keep_d.extend(rng)
        keep_d = sorted(set(keep_d))

        pa = m239.atoms_of(os.path.join(ROOT, "data", "3QN1.cif"), "A", set(keep_p))
        da = m239.atoms_of(path, dchain, set(keep_d))
        tag = t + ("_noins" if a.no_insertions else "_ins")
        out = os.path.join(OUT, f"chimera_{tag}.pdb")
        with open(out, "w") as fh:
            fh.write(f"REMARK  PYR1 / {t} chimera, greedy junctions, "
                     f"insertions {'EXCLUDED' if a.no_insertions else 'INCLUDED'}\n")
            fh.write(f"REMARK  chain A = PYR1 ({len(keep_p)}), "
                     f"chain B = donor ({len(keep_d)}, {n_ins} insertions)\n")
            i = 0
            for n, nm, el, rn, x in pa:
                i += 1; m239.w(fh, i, nm, rn, "A", n, x, el)
            for n, nm, el, rn, x in da:
                i += 1; m239.w(fh, i, nm, rn, "B", n, gx(x[None, :])[0], el)
            fh.write("END\n")

        # measure with the SAME calibrated code as everything else
        X, el = [], []
        for n, nm, e, rn, x in pa:
            X.append(x); el.append(e)
        for n, nm, e, rn, x in da:
            X.append(gx(x[None, :])[0]); el.append(e)
        X = np.array(X)
        cs = m229.chambers_with_voxels(X, el, [True] * len(X))
        cav = cs[0]["vol"] if cs else 0.0
        jd = [d for _, d in J]
        n_tot = len(keep_p) + len(keep_d)
        print(f"{t[:11]:<12}{n_tot:>5}{len(keep_p):>6}{len(keep_d):>5}{n_ins:>5}"
              f"{cav:>8.1f}{cav/PYR1_REF:>7.2f}"
              f"{sum(1 for d in jd if d<1.5):>4}/{len(J):<2}")
        rows.append(dict(target=t, insertions=not a.no_insertions,
                         n_total=n_tot, n_pyr1=len(keep_p), n_donor=len(keep_d),
                         n_inserted=n_ins, cavity=round(cav, 1),
                         x_pyr1=round(cav / PYR1_REF, 2),
                         clean=sum(1 for d in jd if d < 1.5), n_junction=len(J),
                         pdb=out))
    json.dump(rows, open(os.path.join(OUT, "chimera_measured.json"), "w"), indent=1)
    print(f"\nPYR1 = {PYR1_REF} A^3 in {len(pres)} residues.")
    print("⚠ Splice map only: no loop closure, no repacking, no sequence design.")
    print("  The cavity is an UPPER BOUND on what the backbone offers.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
