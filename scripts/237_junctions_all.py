#!/usr/bin/env python
r"""
237_junctions_all.py -- splice feasibility for ALL 261 homologs, with v5 cavities.

⚠ WHY 214 ONLY COVERED 41. Its target list was
    [t for t in cavmap if t in aln and os.path.exists(f"raw/{t}.cif")]
and **every AlphaFold model in the survey is a .pdb**, so all 195 of them were
silently dropped. 220 of 261 homologs — including most of the 99 whose pocket
beats PYR1 — were never screened for splice points at all.

⚠ AND WHY THE AXIS CHANGED. §145 plotted donor cavity against SEQUENCE IDENTITY
and concluded grafting was dominated. Jannis: "I still do not know if identity is
the correct measurement... The big pockets are supposed to be different. They're
bigger and need different backbone architecture for that. Just hopefully some of
them happen to align well enough near the start of the helix that they can be
spliced together in a manner that conserves the gate, the latch, the opening, a
few extra residues, and the HAB1 binding interface."

He is right, and the existing 41-donor table already showed identity hiding the
answer:
    2bk0A00   344.5 A^3   12.2 % identity   7 of 8 clean junctions at 0.76 A
    6awvC00   319.0 A^3   10.8 % identity   7 of 8 clean junctions at 1.16 A
Both look ungraftable on identity and are spliceable on backbone geometry. The
discrimination is not vacuous either: CoxG (2pcsA00), the biggest pocket in the
set at 569.6 A^3, gets only 3 of 8 clean junctions at 2.07 A — the junction axis
correctly rejects the one you would reach for first.

WHAT A JUNCTION IS HERE: PYR1 keeps the gate, the latch, the tunnel positions and
the HAB1 interface (plus `flank` residues either side); the donor supplies
everything else. `scan_junctions` then slides each boundary through its whole gap
and takes the position of least CA deviation after superposition. "Clean" means
that deviation is under 1.5 A — i.e. the two backbones can be joined there
without a hinge, which is Jannis's "a little flexibility ... overlap them well".

⚠ CAVITIES COME FROM cavity_v5.json, not the old homolog_cavities.csv that 214
used: the calibrated probe (1.4 A, validated against ABA) and the domain-lining
attribution. 3oquB00 reads 203.7 in the old file and differs in v5.
⚠ 3QRZ IS OBSOLETE (superseded by 4JDL, 2013-03-13) and is flagged, not dropped
silently, so it cannot re-enter a shortlist unnoticed.
"""
import argparse, csv, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
imp = __import__("importlib").import_module
m214 = imp("214_schema_foldseek")
m221 = imp("221_cavity_domain_trim")
OUT = os.path.join(ROOT, "results", "schema_chimera")
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
AA3 = m214.AA3
OBSOLETE = {"3QRZ"}


def read_any(path, chain=None):
    """ca, heavy, seq -- .cif via 214's reader, .pdb here (first model only)."""
    if path.endswith(".cif"):
        return m214.read_cif(path, chain)
    ca, heavy, seq = {}, {}, {}
    for l in open(path):
        if l.startswith("ENDMDL"):
            break
        if not l.startswith("ATOM"):
            continue
        if chain and l[21] != chain:
            continue
        e = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        if e == "H":
            continue
        try:
            n = int(l[22:26])
        except ValueError:
            continue
        p = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        heavy.setdefault(n, []).append(p)
        seq[n] = AA3.get(l[17:20].strip(), "X")
        if l[12:16].strip() == "CA":
            ca[n] = p
    return ca, heavy, seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flank", type=int, default=2)
    ap.add_argument("--nperm", type=int, default=100)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    pca, pheavy, pseq = m214.read_cif(os.path.join(ROOT, "data", "3QN1.cif"), "A")
    pres = sorted(pca)
    cons = m214.contacts(pheavy)
    print(f"PYR1: {len(pres)} residues, {len(cons)} contacts")

    v5 = {r["target"]: r for r in json.load(open(os.path.join(
        ROOT, "results", "homolog_cavities", "cavity_v5.json")))["structures"]
        if "error" not in r}
    aln = {}
    for line in open(f"{ROOT}/results/foldseek/fs_hits.tsv"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 11 and f[1] not in aln:
            aln[f[1]] = f
    rng = np.random.default_rng(0)
    K0 = m214.keep_set(a.flank)

    targets = []
    for t in v5:
        if t not in aln:
            continue
        for ext in (".cif", ".pdb"):
            p = os.path.join(RAW, t + ext)
            if os.path.exists(p):
                targets.append((t, p)); break
    print(f"{len(targets)} homologs with BOTH a structure and an alignment "
          f"(214 saw only the .cif subset)\n")

    rows, skipped = [], []
    for t, path in targets:
        f = aln[t]
        # the CATH domain name's 5th character IS the chain; AlphaFold is chain A
        dchain = (t[4] if (path.endswith(".cif") and len(t) > 4) else None)
        try:
            dca, dheavy, dseq = read_any(path, dchain)
        except Exception as e:                                   # noqa: BLE001
            skipped.append((t, f"read: {e}")); continue
        dres = sorted(dca)
        if len(dres) < 60:
            skipped.append((t, f"{len(dres)} CA")); continue
        m = m214.map_from_alignment(f[9], f[10], int(f[5]), int(f[7]), pres, dres)
        m = {k: v for k, v in m.items() if k in pca and v in dca}
        if len(m) < 60:
            skipped.append((t, f"{len(m)} mapped")); continue
        ks = sorted(m)
        P = np.array([pca[k] for k in ks]); Q = np.array([dca[m[k]] for k in ks])
        R, pm, qm = m214.kabsch(P, Q)
        fx = lambda X: (X - pm) @ R + qm                          # noqa: E731
        dev = {k: float(np.linalg.norm(fx(pca[k][None, :])[0] - dca[m[k]]))
               for k in ks}
        rms = float(np.sqrt(np.mean([d * d for d in dev.values()])))
        if rms > 8.0:                      # assert on the RESULT
            skipped.append((t, f"RMSD {rms:.1f}")); continue
        K = K0 & set(ks)
        parent, J = m214.scan_junctions(ks, K, dev)
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
        rows.append(dict(
            target=t, pdb=t[:4].upper(),
            kind=v5[t]["kind"],
            cavity=v5[t]["main_median"],
            obsolete=t[:4].upper() in OBSOLETE,
            ident=100 * float(f[2]), alntm=float(f[3]),
            n_aln=len(ks), rmsd=round(rms, 2), n_keep=len(K), n_pyr1=npy,
            E=E, E_frac=round(E / max(null.mean(), 1e-9), 3),
            n_junction=len(J),
            junction_cost=round(float(np.sum(jd)), 2),
            mean_junction=round(float(np.mean(jd)), 2),
            clean=sum(1 for d in jd if d < 1.5)))
    json.dump(rows, open(os.path.join(OUT, f"junctions_all_flank{a.flank}.json"),
                         "w"), indent=1)                  # PERSIST FIRST
    print(f"wrote {OUT}/junctions_all_flank{a.flank}.json")
    print(f"{len(rows)} scored, {len(skipped)} skipped")
    for t, why in skipped[:8]:
        print(f"   skip {t[:44]:<46} {why}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
