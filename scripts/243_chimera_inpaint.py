#!/usr/bin/env python
r"""
243_chimera_inpaint.py -- write the chimera as ONE chain in sequential order and
emit an RFd3 contig that inpaints the seams.

⚠ WHY THE EARLIER FILES COULD NOT BE INPAINTED. 239/241 wrote PYR1 as chain A and
the donor as chain B, which is right for LOOKING at a splice but wrong as a
model: the chimera is a SINGLE polypeptide that alternates parent along the
alignment. Split across two chains it appeared to have 12 fragments, several of
them 1-2 residues, and the sequential order was lost.

Here the chain is rebuilt in alignment order -- for each aligned position take
PYR1's residue if `parent` says PYR1, else the donor's -- and renumbered 1..N.
A SEAM is then a consecutive pair whose CA-CA distance is not ~3.8 A, i.e. a
place where the two parents' backbones do not meet. Those are what inpainting has
to close.

THE CONTIG. Segments that are continuous become fixed motif; seams become
scaffolded gaps sized from the actual CA-CA distance (about 3.3 A of span per
residue, floored at 2 and given +/- slack). ⚠ Per §127, the gap MINIMUM is never
set above what the geometry needs -- over-specified linkers buy appendages, not
structure (excess vs appendage rho = +0.689).

⚠ ABA IS IN THE INPUT. Every design before §129 was scaffolded around an empty
pocket. The ligand is carried through so the inpainted loops are built around it.
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
imp = __import__("importlib").import_module
m214, m237, m240, m239 = (imp("214_schema_foldseek"), imp("237_junctions_all"),
                          imp("240_junction_greedy"), imp("239_build_chimera"))
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
OUT = os.path.join(ROOT, "results", "chimera")
CA_IDEAL, CA_TOL = 3.80, 0.45


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--flank", type=int, default=2)
    ap.add_argument("--cut", type=float, default=1.5)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    t = a.target
    path = None
    for ext in (".cif", ".pdb"):
        p = os.path.join(RAW, t + ext)
        if os.path.exists(p): path = p; break
    assert path, f"no structure for {t}"
    pca, ph, ps = m214.read_cif(os.path.join(ROOT, "data", "3QN1.cif"), "A")
    pres = sorted(pca)
    aln = {}
    for line in open(f"{ROOT}/results/foldseek/fs_hits.tsv"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 11 and f[1] not in aln:
            aln[f[1]] = f
    f = aln[t]
    dchain = (t[4] if (path.endswith(".cif") and len(t) > 4) else None)
    dca, dh, dsq = m237.read_any(path, dchain)
    m = m214.map_from_alignment(f[9], f[10], int(f[5]), int(f[7]), pres, sorted(dca))
    m = {k: v for k, v in m.items() if k in pca and v in dca}
    ks = sorted(m)
    P = np.array([pca[k] for k in ks]); Q = np.array([dca[m[k]] for k in ks])
    R, pm, qm = m214.kabsch(P, Q)
    gx = lambda Y: (Y - qm) @ R.T + pm                            # noqa: E731
    dev = {k: float(np.linalg.norm(((pca[k] - pm) @ R + qm) - dca[m[k]]))
           for k in ks}
    K = m214.keep_set(a.flank) & set(ks)
    parent, J = m240.scan_greedy(ks, K, dev, a.cut)

    # ---- ONE chain, alignment order, renumbered ----
    chain, prov = [], []
    for k in ks:
        if parent[k] == "PYR1":
            chain.append(pca[k]); prov.append(("PYR1", k))
        else:
            chain.append(gx(dca[m[k]][None, :])[0]); prov.append(("DON", m[k]))
    X = np.array(chain)
    d = np.linalg.norm(X[1:] - X[:-1], axis=1)
    seam = [i for i, dd in enumerate(d) if abs(dd - CA_IDEAL) > CA_TOL]
    npy = sum(1 for s, _ in prov if s == "PYR1")
    print(f"{t}: {len(X)} residues in ONE chain "
          f"({npy} PYR1, {len(X)-npy} donor)")
    print(f"parent switches along the chain: "
          f"{sum(1 for i in range(1,len(prov)) if prov[i][0]!=prov[i-1][0])}")
    print(f"seams (CA-CA outside {CA_IDEAL}+/-{CA_TOL} A): {len(seam)}\n")

    # ---- contiguous runs -> contig ----
    runs, start = [], 0
    for i in seam:
        runs.append((start, i)); start = i + 1
    runs.append((start, len(X) - 1))
    runs = [(s, e) for s, e in runs if e >= s]
    #: ⚠ a 1-2 residue "fixed segment" between two gaps is not a motif, it is
    #: noise -- RFd3 cannot anchor on it and it forces two needless junctions.
    #: Absorb anything shorter than MIN_SEG into the surrounding inpainted span.
    MIN_SEG = 3
    merged, i_ = [], 0
    while i_ < len(runs):
        s_, e_ = runs[i_]
        if e_ - s_ + 1 < MIN_SEG and merged:
            merged[-1] = (merged[-1][0], e_)      # swallow into the previous run
            merged[-1] = (merged[-1][0], merged[-1][1])
            runs[i_] = None
            i_ += 1
            continue
        merged.append((s_, e_)); i_ += 1
    # recompute: a swallowed segment means the gap now spans from the previous
    # kept run's end to the next kept run's start
    keep = []
    for s_, e_ in merged:
        if e_ - s_ + 1 >= MIN_SEG:
            keep.append((s_, e_))
    dropped = len(runs) - len(keep)
    if dropped:
        print(f"absorbed {dropped} segment(s) shorter than {MIN_SEG} residues "
              f"into the adjacent inpainted span\n")
    runs = keep
    parts, spec = [], []
    for j, (s, e) in enumerate(runs):
        n = e - s + 1
        src = prov[s][0]
        parts.append(f"{src}:{s+1}-{e+1} ({n})")
        #: `first`/`last` are the PARENT's own numbers; `lo`/`hi` are the
        #: positions in the RENUMBERED chain, which is what the contig needs.
        spec.append(dict(kind="fixed", n=n, src=src, lo=s + 1, hi=e + 1,
                         first=prov[s][1], last=prov[e][1]))
        if j < len(runs) - 1:
            nxt = runs[j + 1][0]
            gap_d = float(np.linalg.norm(X[nxt] - X[e]))
            n_missing = nxt - e - 1
            need = max(2, int(np.ceil(gap_d / 3.3)) - 1, n_missing)
            spec.append(dict(kind="gap", lo=need, hi=need + 4,
                             span=round(gap_d, 2), dropped=n_missing))
            parts.append(f"[inpaint {need}-{need+4}, span {gap_d:.1f} A"
                         + (f", {n_missing} absorbed" if n_missing else "") + "]")
    print("CHAIN LAYOUT")
    for p in parts:
        print("   " + p)
    nfix = sum(s["n"] for s in spec if s["kind"] == "fixed")
    ngap = sum(s["lo"] for s in spec if s["kind"] == "gap")
    print(f"\nfixed {nfix} residues, {sum(1 for s in spec if s['kind']=='gap')} "
          f"gaps needing >= {ngap} inpainted residues")
    # ---- write the motif PDB (fixed segments only) ----
    outp = os.path.join(OUT, f"inpaint_{t}.pdb")
    keep_p = {k for k in ks if parent[k] == "PYR1"}
    keep_d = {m[k] for k in ks if parent[k] == "DON"}
    pa = m239.atoms_of(os.path.join(ROOT, "data", "3QN1.cif"), "A", keep_p)
    da = m239.atoms_of(path, dchain, keep_d)
    idx = {(s, r): i + 1 for i, (s, r) in enumerate(prov)}
    with open(outp, "w") as fh:
        fh.write(f"REMARK  chimera PYR1/{t}, ONE chain, renumbered in alignment order\n")
        fh.write(f"REMARK  {npy} PYR1 + {len(X)-npy} donor; {len(seam)} seams to inpaint\n")
        i = 0
        for n, nm, el, rn, x in pa:
            if ("PYR1", n) not in idx: continue
            i += 1; m239.w(fh, i, nm, rn, "A", idx[("PYR1", n)], x, el)
        for n, nm, el, rn, x in da:
            if ("DON", n) not in idx: continue
            i += 1; m239.w(fh, i, nm, rn, "A", idx[("DON", n)],
                           gx(x[None, :])[0], el)
        # ABA, carried through -- every design before 129 had an EMPTY pocket
        nl = 0
        for l in open(os.path.join(ROOT, "data", "3QN1_complex_auth_aba.pdb")):
            if l.startswith("HETATM") and l[17:20] == "A8S":
                nl += 1; fh.write(l)
        fh.write("END\n")
    print(f"ABA atoms carried into the input: {nl}")
    # ---- the RFd3 contig, dialect 2 ----
    cg = []
    for sp in spec:
        if sp["kind"] == "fixed":
            cg.append(f"A{sp['lo']}-{sp['hi']}")
        else:
            cg.append(f"{sp['lo']}-{sp['hi']}")
    contig = ",".join(cg)
    print(f"\nCONTIG  {contig}")
    json.dump(dict(target=t, n_res=len(X), n_pyr1=npy, n_donor=len(X)-npy,
                   seams=len(seam), spec=spec, pdb=outp, contig=contig,
                   ligand="A8S"),
              open(os.path.join(OUT, f"inpaint_{t}.json"), "w"), indent=1)
    print(f"wrote {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
