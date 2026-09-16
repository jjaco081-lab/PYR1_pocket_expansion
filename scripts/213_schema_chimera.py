#!/usr/bin/env python
r"""
213_schema_chimera.py -- SCHEMA recombination: put PYR1's functional parts on a
LARGER donor scaffold, and score how much structure the junctions break.

WHY THIS ARM EXISTS, and why the earlier one closed the wrong question.
README 86 measured whether a donor's EXISTING backbone already sits in PYR1's
conformation at the gate/latch/interface, found 3.5-4 A core RMSD with gates
5-9 A out of register, and closed "grafting". Jannis: that is the wrong test.

  "I suppose when we graft that means taking part of PYR1's sequence and part of
   the new protein's sequence... it would not matter what the gate RMSD is as
   long as there are good places to anchor all of the infrastructure near the
   tunnel ABA enters, the loops, and the HAB1 interface."

He is right. In a CHIMERA the donor's gate is REPLACED, not superimposed, so its
conformation is irrelevant -- what matters is whether backbone positions exist
where the two chains can be JOINED. Measured in-session, they do, in every large
donor, and they bracket every segment we need:

  2PCS (CoxG, 570 A^3)   8 junctions <1.5 A   brackets gate, latch, interface
  2BK0 (Api g 1, 344)    7                    brackets all three
  6AWV (Ara h 8, 319)    6                    brackets all three
  3OQU (PYL9, 204)      98                    brackets all three  <- the control

⚠ JUNCTION COUNT IS NECESSARY, NOT SUFFICIENT. A junction can be geometrically
placeable and still gut the fold by severing the contacts that hold it together.
That is exactly what SCHEMA's E measures, and it is what this script adds.

  contact  : residue pair with any heavy atoms within 4.5 A and |i-j| >= 5
             (the |i-j| cut removes trivial local-backbone contacts)
  E        : number of contacts whose two residues come from DIFFERENT parents
  lower E  = fewer interactions broken = more likely to fold

⚠ CONTROL, and it is the point of the pilot: 3OQU is PYL9, a real ABA receptor at
1.38 A CA RMSD. A PYL9 chimera must score LOW E and should fold and bind. If the
method cannot show that, its verdict on CoxG means nothing. A random-junction
null is also computed, because a keep-set this large may force E to be similar
whatever we do -- in which case the optimisation is doing no work.

KEPT FROM PYR1 (Jannis): HAB1 interface, gate and latch loops, the tunnel-opening
residues, plus a flank either side to bolster the motifs. The flank is SWEPT
(0-4) rather than assumed, because he asked for "1-2 ... more if that's what the
best chimeras want".
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "schema_chimera")

#: the four HAB1-interface motif segments used throughout the RFd3 arm
IFACE = [(58, 65), (81, 92), (111, 121), (146, 168)]
GATE = (85, 89)
LATCH = (115, 117)
#: tunnel / mouth-opening residues: cavity-lining positions on the MOUTH side of
#: the pocket axis, measured in-session (depth axis > 0, mouth = gate/latch end)
TUNNEL = [81, 83, 87, 89, 92, 117, 159]
DONORS = {"2PCS": ("2pcsA00", "CoxG", 570),
          "2BK0": ("2bk0A00", "Api g 1", 344),
          "6AWV": ("6awvC00", "Ara h 8", 319),
          "2NS9": ("2ns9B01", "APE2225", 455),
          "3OQU": ("3oquB00", "PYL9", 204)}


def load_ca_and_heavy(path, chain=None):
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
    ca, heavy = {}, {}
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
        if g(r, "label_atom_id") == "CA":
            ca[n] = p
    return ca, heavy


def contacts(heavy, cut=4.5, seq_sep=5):
    ks = sorted(heavy)
    out = []
    for a in range(len(ks)):
        for b in range(a + 1, len(ks)):
            i, j = ks[a], ks[b]
            if j - i < seq_sep:
                continue
            A, B = np.array(heavy[i]), np.array(heavy[j])
            if np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2).min() < cut:
                out.append((i, j))
    return out


def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    return V @ np.diag([1, 1, d]) @ Wt, P.mean(0), Q.mean(0)


def correspond(pyr_ca, don_ca, iters=12):
    kp, kd = sorted(pyr_ca), sorted(don_ca)
    m = {kp[i]: kd[int(i * len(kd) / len(kp))] for i in range(len(kp))}
    for _ in range(iters):
        ks = list(m)
        P = np.array([pyr_ca[k] for k in ks]); Q = np.array([don_ca[m[k]] for k in ks])
        R, pm, qm = kabsch(P, Q)
        f = lambda X: (X - pm) @ R + qm                           # noqa: E731
        B = np.array([don_ca[k] for k in kd]); new = {}
        for k in kp:
            d = np.linalg.norm(B - f(pyr_ca[k][None, :]), axis=1)
            j = int(np.argmin(d))
            if d[j] < 6.0:
                new[k] = kd[j]
        if new == m:
            break
        m = new
    ks = sorted(m)
    P = np.array([pyr_ca[k] for k in ks]); Q = np.array([don_ca[m[k]] for k in ks])
    R, pm, qm = kabsch(P, Q)
    f = lambda X: (X - pm) @ R + qm                               # noqa: E731
    dev = {k: float(np.linalg.norm(f(pyr_ca[k][None, :])[0] - don_ca[m[k]]))
           for k in ks}
    return m, dev


def keep_set(flank):
    S = set()
    for a, b in IFACE + [GATE, LATCH]:
        S |= set(range(a - flank, b + flank + 1))
    for t in TUNNEL:
        S |= set(range(t - flank, t + flank + 1))
    return S


def schema_E(cons, parent):
    """contacts whose two residues come from different parents"""
    return sum(1 for i, j in cons
               if parent.get(i) is not None and parent.get(j) is not None
               and parent[i] != parent[j])


def place_junctions(ks, K, dev, window=6):
    """Choose WHERE each crossover sits, instead of accepting the keep-set edge.

    ⚠ The pilot took the keep-set boundary as the junction and the control
    FAILED: PYL9 -- a real ABA receptor at 1.38 A CA RMSD, the easiest chimera in
    the set -- scored the WORST E, because raw E mostly counts how many PYR1
    positions are kept rather than how hard the splice is.

    A crossover may slide outward into the donor region (never into the kept
    motifs) to reach a position where the two backbones actually superimpose.
    Cost is the CA deviation AT the splice point: 0.4 A is nearly a direct
    peptide join, 3.4 A means rebuilding backbone either side.
    """
    parent = {k: ("PYR1" if k in K else "DON") for k in ks}
    runs, cur = [], []
    for k in ks:                                   # maximal runs of kept positions
        if k in K:
            cur.append(k)
        elif cur:
            runs.append(cur); cur = []
    if cur:
        runs.append(cur)
    idx = {k: i for i, k in enumerate(ks)}
    junctions = []
    for run in runs:
        for end, step in ((run[0], -1), (run[-1], +1)):
            best, bestd = end, dev.get(end, 9.9)
            for off in range(1, window + 1):
                i = idx[end] + step * off
                if not (0 <= i < len(ks)):
                    break
                k = ks[i]
                if k in K:                          # never eat into another motif
                    break
                if dev.get(k, 9.9) < bestd:
                    best, bestd = k, dev[k]
            # extend PYR1 out to the chosen splice point
            lo, hi = sorted((idx[end], idx[best]))
            for i in range(lo, hi + 1):
                parent[ks[i]] = "PYR1"
            junctions.append((best, bestd))
    return parent, junctions


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nperm", type=int, default=200)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    pyr_ca, pyr_heavy = load_ca_and_heavy(os.path.join(ROOT, "data", "3QN1.cif"), "A")
    cons = contacts(pyr_heavy)
    print(f"PYR1 chain A: {len(pyr_ca)} residues, {len(cons)} contacts "
          f"(heavy atoms <4.5 A, |i-j|>=5)\n")
    rng = np.random.default_rng(0)
    rows = []
    for pdb, (fn, name, cav) in DONORS.items():
        path = os.path.join(ROOT, "results", "homolog_cavities", "raw", f"{fn}.cif")
        if not os.path.exists(path):
            print(f"{pdb}: missing {path}"); continue
        don_ca, _ = load_ca_and_heavy(path)
        m, dev = correspond(pyr_ca, don_ca)
        print(f"{name} ({pdb}), cavity {cav} Å³ — {len(m)} aligned, "
              f"CA RMSD {np.sqrt(np.mean([d*d for d in dev.values()])):.2f} Å")
        ks = sorted(m)
        for flank in (0, 1, 2, 3, 4):
            K = keep_set(flank) & set(m)
            parent, J = place_junctions(ks, K, dev)
            E = schema_E(cons, parent)
            n_pyr1 = sum(1 for k in ks if parent[k] == "PYR1")
            jdev = [d for _, d in J]
            clean = sum(1 for d in jdev if d < 1.5)
            # NULL matched on THIS donor AND on the realised PYR1 count, so E is
            # comparable across donors whose alignments differ in length.
            null = []
            for _ in range(a.nperm):
                pick = set(rng.choice(ks, size=n_pyr1, replace=False))
                null.append(schema_E(cons, {k: ("PYR1" if k in pick else "DON")
                                            for k in ks}))
            null = np.array(null)
            z = (E - null.mean()) / max(null.std(), 1e-9)
            rows.append(dict(pdb=pdb, name=name, cavity=cav, flank=flank,
                             n_keep=len(K), n_pyr1=n_pyr1, E=E,
                             E_frac=round(E / max(null.mean(), 1e-9), 3),
                             z=round(float(z), 2), n_junction=len(J),
                             junction_cost=round(float(np.sum(jdev)), 2),
                             mean_junction_dev=round(float(np.mean(jdev)), 2),
                             clean_junctions=clean))
            print(f"   ±{flank}: PYR1 {n_pyr1:>3}  E {E:>4}  "
                  f"E/null {E/max(null.mean(),1e-9):>5.2f}  z {z:>6.1f}   "
                  f"junctions {len(J):>2}, {clean} clean(<1.5Å), "
                  f"total cost {np.sum(jdev):>5.1f} Å")
        print()
    json.dump(rows, open(os.path.join(OUT, "schema.json"), "w"), indent=1)
    print(f"wrote {OUT}/schema.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
