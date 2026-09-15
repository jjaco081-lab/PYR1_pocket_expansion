#!/usr/bin/env python
r"""
202_water_positions.py -- OPTION B. Do crystallographic pocket waters mark the
positions that matter?

WHY POSITIONS AND NOT SUBSTITUTIONS. Beltran-45 (87) showed the ligand-blind
frequency prior transfers at 64 % recall and then saturates at Tian's 67 %
vocabulary ceiling, so the remaining headroom is POSITIONAL. 201 has just shown
the ligand-conditioning ceiling is real (4.57x within- over between-ligand
Jaccard after deduping), so there IS structure to find -- but five substitution-
level conditioning attempts have failed. This asks a positional question with a
physical observable instead.

THE PHYSICS. ABA leaves only ~40 A^3 of probe-accessible space in the closed
pocket, holding 1.5 +- 0.7 waters (README 114e). A side chain that displaces an
ordered water pays or gains a desolvation term, and fa_sol produced the largest
non-volume AUC gap in this project (+0.133 on fludioxonil, README 101). If
ordered water marks where the pocket is thermodynamically tender, water proximity
should rank the mutable positions.

NUMBERING, verified not assumed. Chain and scheme are chosen BY ASSERTION: for
each structure every chain and both auth/label schemes are scored against the 18
sd03 position identities, and the best is used. Measured here: 3QN1 chain A
auth = 18/18, 3K3K chain A auth = 18/18, 4WVO 14/18 and 8EY0 12/18 -- the
4WVO/8EY0 shortfalls are their ENGINEERED mutations, so those two are reported
but excluded from the WT reference, and any position whose identity disagrees is
dropped for that structure.

PREDICTIONS, FIXED BEFORE THE RUN:
  1. SANITY. The 18 sd03 positions sit closer to pocket waters than a random
     equal-sized set of PYR1 surface positions, p < 0.05.
  2. THE TEST. Ranking positions by pocket-water contact beats ranking them by
     blind sd03 position frequency, at position recall@5 on Beltran-45.
  3. THE CONTROL. Permuting water counts across the 18 positions (n=2000) does
     not reproduce the gain.

⚠ POWER IS LOW BY CONSTRUCTION: there are only 18-20 candidate positions, so
recall@5 moves in steps of ~1/|truth|. A null here is weak evidence, and a win
needs the permutation control to mean anything.
"""
import json, os, re, sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                          # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "water_positions")
AA = set("ACDEFGHIKLMNPQRSTVWY")
BB = {"N", "CA", "C", "O", "OXT"}
WT_REFS = ["3QN1", "3K3K"]
ALL_REFS = ["3QN1", "3K3K", "4WVO", "8EY0"]
CONTACT = 5.0           # A, side-chain atom to water oxygen
#: sd03's 18 mutable positions and their WT identities -- the assertion target
POS = {59: "LYS", 81: "VAL", 83: "VAL", 87: "LEU", 89: "ALA", 92: "SER",
       94: "GLU", 108: "PHE", 110: "ILE", 117: "LEU", 120: "TYR", 122: "SER",
       141: "GLU", 159: "PHE", 160: "ALA", 163: "VAL", 164: "VAL", 167: "ASN"}
LIGANDS = {"A8S", "ABA", "MND", "4WV"}          # ABA and engineered-ligand codes


def load_cif(pdb):
    cols, rows, inl = {}, [], False
    for line in open(os.path.join(ROOT, "data", f"{pdb}.cif")):
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
    return cols, rows


def resolve(pdb):
    """(chain, scheme, identity_hits) chosen by asserting residue IDENTITY."""
    cols, rows = load_cif(pdb)
    g = lambda r, k: r[cols[k]]                                  # noqa: E731
    best = None
    for scheme in ("auth_seq_id", "label_seq_id"):
        for ch in sorted({g(r, "auth_asym_id") for r in rows
                          if g(r, "group_PDB") == "ATOM"}):
            seen = {}
            for r in rows:
                if g(r, "group_PDB") != "ATOM" or g(r, "auth_asym_id") != ch:
                    continue
                try:
                    seen[int(g(r, scheme))] = g(r, "label_comp_id")
                except ValueError:
                    pass
            ok = sum(1 for n, aa in POS.items() if seen.get(n) == aa)
            if best is None or ok > best[2]:
                best = (ch, scheme, ok, seen)
    return best


def water_contacts(pdb):
    """{position: n_pocket_waters_within CONTACT of its side chain}, plus info."""
    cols, rows = load_cif(pdb)
    g = lambda r, k: r[cols[k]]                                  # noqa: E731
    ch, scheme, nid, seen = resolve(pdb)
    xyz = lambda r: np.array([float(g(r, "Cartn_x")), float(g(r, "Cartn_y")),
                              float(g(r, "Cartn_z"))])           # noqa: E731
    side = defaultdict(list)
    lig, wat = [], []
    for r in rows:
        if g(r, "type_symbol") == "H":
            continue
        comp = g(r, "label_comp_id")
        if comp == "HOH":
            wat.append(xyz(r)); continue
        if comp in LIGANDS:
            lig.append(xyz(r)); continue
        if g(r, "group_PDB") != "ATOM" or g(r, "auth_asym_id") != ch:
            continue
        try:
            n = int(g(r, scheme))
        except ValueError:
            continue
        if n in POS and seen.get(n) == POS[n] and g(r, "label_atom_id") not in BB:
            side[n].append(xyz(r))
    wat = np.array(wat) if wat else np.zeros((0, 3))
    lig = np.array(lig) if lig else np.zeros((0, 3))
    # POCKET waters only: within 8 A of the ligand when a ligand is present,
    # otherwise within 8 A of any of the asserted mutable side chains
    if len(wat):
        anchor = lig if len(lig) else np.array(
            [p for v in side.values() for p in v])
        if len(anchor):
            d = np.linalg.norm(wat[:, None, :] - anchor[None, :, :], axis=2).min(1)
            pocket = wat[d < 8.0]
        else:
            pocket = wat
    else:
        pocket = wat
    out = {}
    for n, pts in side.items():
        P = np.array(pts)
        if len(pocket):
            d = np.linalg.norm(pocket[:, None, :] - P[None, :, :], axis=2).min(1)
            out[n] = int((d < CONTACT).sum())
        else:
            out[n] = 0
    return out, dict(chain=ch, scheme=scheme, identity=f"{nid}/18",
                     n_water=len(wat), n_pocket=len(pocket), n_lig=len(lig),
                     n_pos=len(side))


def blind_positions():
    hdr, recs = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    pc = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    c = Counter()
    for r in recs:
        for col in pc:
            v = (r.get(col) or "").strip().upper()
            if len(v) == 1 and v in AA and v != col[0]:
                c[int(col[1:])] += 1
    return c


def main():
    os.makedirs(OUT, exist_ok=True)
    print("Resolving chain/scheme BY IDENTITY ASSERTION")
    wc, info = {}, {}
    for pdb in ALL_REFS:
        w, i = water_contacts(pdb)
        wc[pdb], info[pdb] = w, i
        print(f"   {pdb}: chain {i['chain']} {i['scheme']} identity {i['identity']}"
              f"  waters {i['n_water']} (pocket {i['n_pocket']})"
              f"  ligand atoms {i['n_lig']}  positions asserted {i['n_pos']}")

    agg = {n: float(np.mean([wc[p][n] for p in WT_REFS if n in wc[p]]))
           for n in POS if any(n in wc[p] for p in WT_REFS)}
    print(f"\nWT reference = {WT_REFS}; {len(agg)} positions with a water count")
    blind = blind_positions()
    print(f"\n{'pos':>5}{"wt":>5}{'waters(WT mean)':>17}{'sd03 clones':>13}")
    for n in sorted(agg, key=lambda x: -agg[x]):
        print(f"{n:>5}{POS[n]:>5}{agg[n]:>17.1f}{blind.get(n,0):>13}")

    # ---- P2: Beltran position recall ---------------------------------------
    bel = json.load(open(os.path.join(ROOT, "data", "beltran", "win_sensors.json")))
    truth = set()
    for s in bel["sensors"]:
        for m in s.get("muts", []):
            mm = re.match(r"^[A-Z](\d+)[A-Z]$", m)
            if mm:
                truth.add(int(mm.group(1)))
    truth &= set(agg)
    print(f"\nBeltran-45 mutated positions inside the asserted set: "
          f"{len(truth)} of {len(agg)}  {sorted(truth)}")
    rank_w = sorted(agg, key=lambda x: (-agg[x], x))
    rank_b = sorted(agg, key=lambda x: (-blind.get(x, 0), x))

    def rec(rank, n):
        return len(set(rank[:n]) & truth) / max(len(truth), 1)

    print(f"\n{'k':>4}{'water rank':>13}{'blind rank':>13}{'delta':>9}")
    deltas = {}
    for k in (3, 5, 8, 10):
        a, b = rec(rank_w, k), rec(rank_b, k)
        deltas[k] = a - b
        print(f"{k:>4}{a:>13.3f}{b:>13.3f}{a-b:>+9.3f}")

    # ---- P3: permute the water counts ---------------------------------------
    rng = np.random.default_rng(0)
    keys = list(agg)
    vals = np.array([agg[k] for k in keys])
    null5 = []
    for _ in range(2000):
        pv = rng.permutation(vals)
        pm = dict(zip(keys, pv))
        r = sorted(keys, key=lambda x: (-pm[x], x))
        null5.append(rec(r, 5) - rec(rank_b, 5))
    null5 = np.array(null5)
    p = float((null5 >= deltas[5]).mean())
    print(f"\nP3 permutation (n=2000) at k=5: observed {deltas[5]:+.3f}, "
          f"null mean {null5.mean():+.3f}, P(null >= obs) = {p:.3f}")
    print(f"   {'water ranking is doing real work' if p < 0.05 and deltas[5] > 0 else 'NULL -- water proximity adds nothing over blind frequency'}")
    json.dump(dict(water=agg, blind={str(k): v for k, v in blind.items()},
                   truth=sorted(truth), deltas={str(k): v for k, v in deltas.items()},
                   p_perm_k5=p, info=info),
              open(os.path.join(OUT, "water_positions.json"), "w"), indent=1)
    print(f"\nwrote {OUT}/water_positions.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
