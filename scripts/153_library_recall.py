#!/usr/bin/env python
r"""
153_library_recall.py -- the Beltran-45 library-enrichment benchmark.

§47e called this "the single highest-value thing not yet done": 45 real sensors
that NOTHING in this project was built on, scored by the quantity the reframed
goal actually cares about -- top-N recall of real sensors at a fixed library
size. Everything before this was retrospective recovery of a 4-mutation answer
we already knew.

THE TASK. A library is a MENU: a set of positions, each with a set of allowed
residues, wild type always retained. Its size is prod(1 + n_i), matching sd04's
arithmetic (§53). A sensor is CAPTURED when every one of its substitutions is
inside the menu. The score is the fraction of real sensors captured at a given
size -- recall, never precision, because the ground truth is positive-unlabeled:
these are hits that were FOUND, not the best sequences that exist (§51).

WHY BELTRAN IS THE RIGHT TEST SET. It uses 20 design positions, and only 7 of
them (59, 81, 83, 120, 159, 160, 164) overlap the 11 this project has ever
scored. So a method tuned on Tian's vocabulary is not merely being asked to
re-rank familiar substitutions; it is being asked to work where it has never
looked. That is the transfer failure mode §47e point 3 warns about, and it is
measurable here for free.

METHODS, fixed before any result is seen
  oracle       exact smallest menu capturing a given number of sensors -- the
               ceiling, not a method
  frequency    substitutions ranked by their frequency among Tian's sd03 round-1
               clones. This is the ADVANTAGED baseline (§50): it encodes real
               selection outcomes, so failing to beat it is weak evidence
               against a method, and beating it is the only interesting result
  tian_vocab   the ceiling imposed by restricting to Tian's vocabulary at all --
               free to compute and it bounds every Tian-trained method
  random       menus of matched shape drawn uniformly, 2000 draws

⚠ WHAT THIS CANNOT SHOW. Recall at fixed size says a library CONTAINS the known
sensors; it does not say those sensors are findable in it, nor that the library
contains anything better. And Beltran's sensors were drawn from Beltran's own
libraries, so any position they never varied is invisible here -- absence is not
evidence against a position (the same asymmetry as §86's sd04 result).
"""
import json
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_xlsx import table                                          # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "library_recall")
os.makedirs(OUT, exist_ok=True)
AA = set("ACDEFGHIKLMNPQRSTVWY")


def menu_size(menu):
    n = 1
    for v in menu.values():
        n *= (1 + len(v))
    return n


def captured(sensors, menu):
    out = 0
    for s in sensors:
        ok = True
        for m in s["muts"]:
            p, a = m[:-1], m[-1]
            if a not in menu.get(p, ()):
                ok = False
                break
        if ok:
            out += 1
    return out


def curve(sensors, ranked, cap=1e7):
    """Add substitutions in the given order; record (size, recall) as we go."""
    menu, pts = defaultdict(set), []
    for sub in ranked:
        p, a = sub[:-1], sub[-1]
        menu[p].add(a)
        sz = menu_size(menu)
        if sz > cap:
            menu[p].discard(a)
            if not menu[p]:
                del menu[p]
            continue
        pts.append((sz, captured(sensors, menu) / len(sensors)))
    return pts


def at_size(pts, target):
    best = 0.0
    for sz, r in pts:
        if sz <= target:
            best = max(best, r)
    return best


def main():
    bel = json.load(open(os.path.join(ROOT, "data", "beltran", "win_sensors.json")))
    sensors = [s for s in bel["sensors"] if s["muts"]]
    bel_pos = {p for s in sensors for p in (m[:-1] for m in s["muts"])}
    universe = sorted({m for s in sensors for m in s["muts"]})
    print(f"Beltran-45: {len(sensors)} sensors with >=1 substitution, "
          f"{len({s['ligand'] for s in sensors})} ligands, "
          f"{len(bel_pos)} positions used, {len(universe)} distinct substitutions")

    # ---- Tian sd03 frequency, the advantaged baseline ------------------
    hdr, recs = table(f"{SD}/pnas.2519924122.sd03(1).xlsx")
    poscols = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    freq = Counter()
    for r in recs:
        for c in poscols:
            v = (r.get(c) or "").strip().upper()
            if len(v) == 1 and v in AA and v != c[0]:
                freq[c + v] += 1
    tian_vocab = set(freq)
    tian_pos = {s[:-1] for s in tian_vocab}
    print(f"Tian sd03 vocabulary: {len(tian_vocab)} substitutions at "
          f"{len(tian_pos)} positions")

    ov_pos = sorted(bel_pos & tian_pos, key=lambda x: int(x[1:]))
    print(f"position overlap Beltran n Tian: {len(ov_pos)}/{len(bel_pos)}  {ov_pos}")
    ov_sub = universe and len(set(universe) & tian_vocab)
    print(f"substitution overlap: {ov_sub}/{len(universe)}")

    # ---- ceiling from Tian's vocabulary alone --------------------------
    full_tian = {p: {s[-1] for s in tian_vocab if s[:-1] == p} for p in tian_pos}
    n_tian = captured(sensors, full_tian)
    print(f"\nCEILING -- Tian's ENTIRE vocabulary as one menu "
          f"(size {menu_size(full_tian):.3g}) captures "
          f"{n_tian}/{len(sensors)} = {n_tian/len(sensors):.0%} of Beltran's sensors")
    per = defaultdict(lambda: [0, 0])
    for s in sensors:
        per[s["ligand"]][1] += 1
        if captured([s], full_tian):
            per[s["ligand"]][0] += 1
    for lg, (a, b) in sorted(per.items(), key=lambda kv: -kv[1][1]):
        print(f"      {lg:<22}{a:>3}/{b:<3}")

    # ---- the recall-vs-size curves -------------------------------------
    ranked_freq = [s for s, _ in freq.most_common() if s in universe]
    tail = [s for s in universe if s not in set(ranked_freq)]
    print(f"\n{len(ranked_freq)} of Beltran's {len(universe)} substitutions have a "
          f"Tian frequency; {len(tail)} are unranked by it (appended at random)")
    rng = np.random.default_rng(0)
    order = list(ranked_freq) + list(rng.permutation(tail))
    c_freq = curve(sensors, order)
    c_oracle = curve(sensors, sorted(universe,
                                     key=lambda s: -sum(1 for x in sensors
                                                        if s in x["muts"])))

    # ---- random null ----------------------------------------------------
    sizes = [1e3, 1e4, 1e5, 1e6]
    print("\n" + "=" * 78)
    print("RECALL OF REAL SENSORS AT A FIXED LIBRARY SIZE")
    print("=" * 78)
    print(f"{'library size':>14}{'frequency':>12}{'greedy-cover':>14}{'random':>18}")
    res = {}
    for N in sizes:
        f_ = at_size(c_freq, N)
        o_ = at_size(c_oracle, N)
        null = []
        for _ in range(2000):
            perm = list(rng.permutation(universe))
            null.append(at_size(curve(sensors, perm, cap=N), N))
        null = np.array(null)
        p = float((null >= f_).mean())
        print(f"{N:>14.0e}{f_:>12.0%}{o_:>14.0%}"
              f"{null.mean():>10.0%} +- {null.std():.0%}   p={p:.3f}")
        res[f"{N:.0e}"] = dict(freq=f_, greedy=o_, null=float(null.mean()), p=p)

    print("\n   'greedy-cover' ranks substitutions by how many BELTRAN sensors use")
    print("   them, so it sees the answer sheet. It is the achievable ceiling for a")
    print("   ranking method of this shape, not a competitor.")
    print("   The frequency column is the transfer test: Tian's selection outcomes")
    print("   applied to a library that used different positions.")
    json.dump({"ceiling_tian": n_tian / len(sensors), "sizes": res,
               "curve_freq": c_freq, "curve_greedy": c_oracle},
              open(os.path.join(OUT, "beltran_recall.json"), "w"), indent=1)
    print(f"\nwritten to {OUT}/beltran_recall.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
