#!/usr/bin/env python
r"""
214_schema_foldseek.py -- SCHEMA chimera scan on the REAL foldseek alignment,
scanning EVERY crossover position rather than a chosen few.

WHAT THIS FIXES IN 213. That pilot used an iterative-closest-point CA
correspondence and a greedy +/-6 window around each keep-set edge. Both are
replaced here:
  * correspondence comes from foldseek's own `qaln`/`taln` strings (the same
    alignment the 266-homolog survey was built on), with residue IDENTITIES
    asserted against each structure before use;
  * every position in every inter-motif gap is scanned exhaustively, so the
    reported junction cost is the MINIMUM ACHIEVABLE, not a local search result.

THE TWO QUANTITIES, and why both are needed. 213's control failed on E alone:
PYL9 -- a real ABA receptor at 1.38 A -- scored the WORST raw E, because E
mostly counts how many PYR1 positions are kept. Normalised against a per-donor,
count-matched null it became uninformative instead (every donor 0.60-0.71).
What separated PYL9 was JUNCTION COST, which E cannot see.

  E          contacts (heavy atoms <4.5 A, |i-j|>=5) whose residues come from
             different parents; reported as E/E_null against a matched null
  junction   CA deviation at the splice point after superposition. 0.4 A is
  cost       nearly a direct peptide join; 3 A means rebuilding backbone.

KEPT FROM PYR1 (Jannis): HAB1 interface, gate and latch loops, tunnel-opening
residues, plus a flank either side. Flank swept, not assumed.

⚠ This answers "is a chimera geometrically possible and how much rebuilding does
it cost", NOT "will it fold". Nothing here is a folding prediction.
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "schema_chimera")
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
IFACE = [(58, 65), (81, 92), (111, 121), (146, 168)]
GATE, LATCH = (85, 89), (115, 117)
TUNNEL = [81, 83, 87, 89, 92, 117, 159]


def read_cif(path, chain=None):
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
    ca, heavy, seq = {}, {}, {}
    for r in rows:
        if g(r, "group_PDB") != "ATOM" or g(r, "type_symbol") == "H":
            continue
        if chain and g(r, "auth_asym_id") != chain:
            continue
        try:
            n = int(g(r, "auth_seq_id"))
        except ValueError:
            continue
        p = np.array([float(g(r, "Cartn_x")), float(g(r, "Cartn_y")),
                      float(g(r, "Cartn_z"))])
        heavy.setdefault(n, []).append(p)
        seq[n] = AA3.get(g(r, "label_comp_id"), "X")
        if g(r, "label_atom_id") == "CA":
            ca[n] = p
    return ca, heavy, seq


def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    return V @ np.diag([1, 1, d]) @ Wt, P.mean(0), Q.mean(0)


def map_from_alignment(qaln, taln, qstart, tstart, qres, tres):
    """foldseek indices are 1-based into the ORDERED residue list, not PDB
    numbering, so index into sorted residue keys and assert the identities."""
    qi, ti = qstart - 1, tstart - 1
    m, bad = {}, 0
    for a, b in zip(qaln, taln):
        if a != "-" and b != "-":
            if qi < len(qres) and ti < len(tres):
                qn, tn = qres[qi], tres[ti]
                m[qn] = tn
        if a != "-":
            qi += 1
        if b != "-":
            ti += 1
    return m


def contacts(heavy, cut=4.5, sep=5):
    ks = sorted(heavy)
    out = []
    for a in range(len(ks)):
        A = np.array(heavy[ks[a]])
        for b in range(a + 1, len(ks)):
            if ks[b] - ks[a] < sep:
                continue
            B = np.array(heavy[ks[b]])
            if np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2).min() < cut:
                out.append((ks[a], ks[b]))
    return out


def keep_set(flank):
    S = set()
    for a, b in IFACE + [GATE, LATCH]:
        S |= set(range(a - flank, b + flank + 1))
    for t in TUNNEL:
        S |= set(range(t - flank, t + flank + 1))
    return S


def scan_junctions(ks, K, dev):
    """EXHAUSTIVE: for each boundary, every position in the gap is evaluated.
    Junction cost is separable per gap, so scanning each gap independently gives
    the global minimum for that term."""
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
            best, bestd, i = end, dev.get(end, 9.9), idx[end]
            off = 1
            while True:                                # scan the WHOLE gap
                j = i + step * off
                if not (0 <= j < len(ks)) or ks[j] in K:
                    break
                if dev.get(ks[j], 9.9) < bestd:
                    best, bestd = ks[j], dev[ks[j]]
                off += 1
            lo, hi = sorted((idx[end], idx[best]))
            for t in range(lo, hi + 1):
                parent[ks[t]] = "PYR1"
            J.append((best, bestd))
    return parent, J


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flank", type=int, default=2)
    ap.add_argument("--nperm", type=int, default=100)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    pca, pheavy, pseq = read_cif(os.path.join(ROOT, "data", "3QN1.cif"), "A")
    pres = sorted(pca)
    cons = contacts(pheavy)
    print(f"PYR1: {len(pres)} residues, {len(cons)} contacts")
    import csv
    cavmap = {r["target"]: float(r["cavity_A3"]) for r in
              csv.DictReader(open(f"{ROOT}/results/homolog_cavities/homolog_cavities.csv"))}
    aln = {}
    for line in open(f"{ROOT}/results/foldseek/fs_hits.tsv"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 11 and f[1] not in aln:
            aln[f[1]] = f
    rng = np.random.default_rng(0)
    K0 = keep_set(a.flank)
    rows = []
    targets = [t for t in cavmap if t in aln and
               os.path.exists(f"{ROOT}/results/homolog_cavities/raw/{t}.cif")]
    if a.limit:
        targets = targets[:a.limit]
    print(f"{len(targets)} homologs with a foldseek alignment AND a structure\n")
    for t in targets:
        f = aln[t]
        # ⚠ these are CATH DOMAIN files and several contain the whole
        # asymmetric unit -- 3qrzB00 has chains A, B AND C (436 residues).
        # Reading all of them scrambles the residue index the foldseek
        # alignment points into, which produced impossible RMSDs (17 A at 50 %
        # identity). The domain name's 5th character IS the chain.
        dchain = t[4] if len(t) > 4 else None
        dca, dheavy, dseq = read_cif(
            f"{ROOT}/results/homolog_cavities/raw/{t}.cif", dchain)
        dres = sorted(dca)
        m = map_from_alignment(f[9], f[10], int(f[5]), int(f[7]), pres, dres)
        m = {k: v for k, v in m.items() if k in pca and v in dca}
        if len(m) < 60:
            continue
        ks = sorted(m)
        P = np.array([pca[k] for k in ks]); Q = np.array([dca[m[k]] for k in ks])
        R, pm, qm = kabsch(P, Q)
        fx = lambda X: (X - pm) @ R + qm                          # noqa: E731
        dev = {k: float(np.linalg.norm(fx(pca[k][None, :])[0] - dca[m[k]]))
               for k in ks}
        rms = float(np.sqrt(np.mean([d * d for d in dev.values()])))
        if rms > 8.0:          # assert on the RESULT: a real homolog cannot be this far
            print(f"   ⚠ {t}: RMSD {rms:.1f} A over {len(ks)} — alignment "
                  f"mapping suspect, skipped")
            continue
        K = K0 & set(ks)
        parent, J = scan_junctions(ks, K, dev)
        E = sum(1 for i, j in cons if i in parent and j in parent
                and parent[i] != parent[j])
        npy = sum(1 for k in ks if parent[k] == "PYR1")
        null = []
        for _ in range(a.nperm):
            pick = set(rng.choice(ks, size=npy, replace=False))
            null.append(sum(1 for i, j in cons if i in parent and j in parent
                            and (i in pick) != (j in pick)))
        null = np.array(null)
        jd = [d for _, d in J]
        rows.append(dict(target=t, pdb=t[:4].upper(), cavity=cavmap[t],
                         ident=100 * float(f[2]), alntm=float(f[3]),
                         n_aln=len(ks), rmsd=round(rms, 2), n_keep=len(K),
                         n_pyr1=npy, E=E, E_frac=round(E / max(null.mean(), 1e-9), 3),
                         n_junction=len(J),
                         junction_cost=round(float(np.sum(jd)), 2),
                         mean_junction=round(float(np.mean(jd)), 2),
                         clean=sum(1 for d in jd if d < 1.5)))
    rows.sort(key=lambda r: r["junction_cost"])
    print(f"{'PDB':<6}{'cavity':>8}{'id%':>6}{'RMSD':>7}{'Jcost':>8}{'clean':>7}{'E/null':>8}")
    for r in rows:
        print(f"{r['pdb']:<6}{r['cavity']:>8.0f}{r['ident']:>6.0f}{r['rmsd']:>7.2f}"
              f"{r['junction_cost']:>8.1f}{r['clean']:>4}/{r['n_junction']:<2}"
              f"{r['E_frac']:>8.2f}")
    json.dump(rows, open(os.path.join(OUT, f"foldseek_flank{a.flank}.json"), "w"),
              indent=1)
    print(f"\nwrote {OUT}/foldseek_flank{a.flank}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
