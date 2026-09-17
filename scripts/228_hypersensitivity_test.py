#!/usr/bin/env python
r"""
228_hypersensitivity_test.py -- WT PYR1 IS NOT THE OPTIMUM FOR ABA, and that is
now a testable claim rather than a caveat.

SOURCE. Elzinga, Sternburg, Sabbadin, Bartsch, Park, Vaidya, Mosquna, Kaundal,
Wendeborn, Lachia, Karginov, Cutler, "Defining and Exploiting Hypersensitivity
Hotspots to Facilitate Abscisic Acid Agonist Optimization", ACS Chem. Biol.
2019, 14(3), 332-336, Figure S1 (Y2H staining; WT PYR1 first responds at 1 uM).
data/papers/Defining and Exploiting Hypersensitivity Hotspots to Facilitate/

⚠ WHY THIS MATTERS BEYOND A NEW TEST SET. Every method here that reports "WT
ranks best" has been treating that as validation. It is not: PYR1 A160V binds ABA
with Kd 4.7 uM by ITC against WT's >50 uM, a >10x improvement. WT is a
functional compromise, not an affinity optimum. Any scorer, filter or library
rule must ALLOW the WT identity, never assert it is the best one.

THE DATA (LOD by Y2H staining; WT = 1 uM; tested only at 0/0.25/0.5/1/5 uM, so
0.25 is a FLOOR, not a measurement -- F61, V81 and A160 would likely score lower
if lower concentrations had been used).

⚠ TEST DESIGN, and the trap it avoids. These are HYPERSENSITIVE MUTATIONS FOUND,
not a complete ranking: a substitution absent from the list is UNTESTED, not
worse. So this measures RECALL ONLY and every test below is one-sided -- known-
better substitutions should rank better. Nothing here can be read as precision.
See feedback_hits_are_not_optima.

PRE-REGISTERED, before looking at the scan output:
  T1. Do hypersensitive substitutions rank ABOVE WT at their own position?
      Binomial sign test against p = (rank_wt - 1) / 18 per position, which is
      the chance of beating WT by luck GIVEN where WT itself ranks. Using a flat
      0.5 would be wrong wherever WT is already extreme.
  T2. Mean rank of the 18 hypersensitive substitutions versus a permutation null
      that redraws the same number of substitutions from the same positions.
  T3. Does the scorer put WT first at these positions? If it does, that is a
      DEMONSTRATED FAILURE now, not a success.
"""
import json, glob, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCAN = os.path.join(ROOT, "results", "position_scan")
OUT = os.path.join(ROOT, "results", "hypersensitivity")

#: ⚠ TWO TIERS. SI Table S2 gives quantitative ABA EC50 for 8 of these, and 3 of
#: the 8 do NOT improve it: F61M is 0.29x (EC50 2.98 vs WT 0.864, i.e. 3.4x
#: WORSE), I110S 0.75x, I110C 0.99x. The Y2H LOD therefore OVERSTATES. EC50 is
#: the primary benchmark; Y2H is secondary and kept only for the positions EC50
#: never measured.
EC50_CONFIRMED = {160: {"V": 8.80, "C": 5.17}, 141: {"L": 3.63},
                  81: {"I": 3.57}, 61: {"L": 1.41}}
EC50_REFUTED = {61: {"M": 0.29}, 110: {"S": 0.75, "C": 0.99}}

#: position -> {substitution: LOD in uM}. WT PYR1 = 1.0 uM.
HYPER = {
    61:  {"L": 0.25, "M": 0.25},
    81:  {"I": 0.25, "Y": 0.25},
    110: {"C": 0.50, "S": 0.25},
    160: {"C": 0.25, "I": 0.25, "V": 0.25},
    141: {"C": 0.50, "I": 0.25, "L": 0.25, "M": 0.25, "N": 0.50,
          "T": 0.50, "V": 0.50, "W": 0.50, "Y": 0.25},
}
WT_AT = {61: "F", 81: "V", 110: "I", 141: "E", 160: "A"}
WT_LOD = 1.0


def main():
    if "--ec50" in sys.argv:
        for pos, subs in list(HYPER.items()):
            keep = EC50_CONFIRMED.get(pos, {})
            if keep:
                HYPER[pos] = {k: v for k, v in subs.items() if k in keep}
            else:
                del HYPER[pos]
        print("EC50-CONFIRMED SUBSET ONLY: "
              + ", ".join(f"{WT_AT[p]}{p}{s_}" for p in HYPER for s_ in HYPER[p])
              + "\n")
    files = sorted(glob.glob(os.path.join(SCAN, "*.json")))
    if not files:
        print("no position_scan output"); return 1
    by_arm = {}
    for f in files:
        tag = os.path.basename(f).rsplit("_", 1)[0]       # system_arm
        for r in json.load(open(f)):
            by_arm.setdefault(tag, {})[r["pos"]] = r
    print(f"{len(files)} scan files -> {len(by_arm)} system/arm combinations\n")

    os.makedirs(OUT, exist_ok=True)
    summary = {}
    for tag in sorted(by_arm):
        rows = by_arm[tag]
        have = [p for p in HYPER if p in rows]
        if not have:
            continue
        print("=" * 76)
        print(f"{tag}   ({len(have)}/5 hotspot positions present)")
        print("=" * 76)
        beats, total, ranks, exp_p = 0, 0, [], []
        wt_first = 0
        for p in sorted(have):
            r = rows[p]
            order = r["order"]
            wt = WT_AT[p]
            assert r["wt"] == wt, f"pos {p}: scan says WT={r['wt']}, paper says {wt}"
            rw = order.index(wt) + 1
            if rw == 1:
                wt_first += 1
            line = []
            for sub, lod in sorted(HYPER[p].items(), key=lambda x: x[1]):
                if sub not in order:
                    line.append(f"{sub}=n/a"); continue
                rk = order.index(sub) + 1
                ranks.append(rk)
                total += 1
                better = rk < rw
                beats += better
                # P(a random non-WT substitution outranks WT) = (rw-1)/18.
                # ⚠ The first version used (19-rw)/18, which is P(ranked BELOW
                # WT) -- inverted. It made a scorer that puts WT FIRST look like
                # it was failing a test it cannot possibly pass: if WT is rank 1,
                # nothing outranks it and the expectation is 0, not 18.
                exp_p.append((rw - 1) / 18.0)
                line.append(f"{sub}{rk:>2}{'<' if better else '>'}")
            print(f"  {wt}{p:<4} WT ranks {rw:>2}/19   best={order[0]}   "
                  f"hypersensitive: {' '.join(line)}")
        # T1
        from math import comb
        pe = float(np.mean(exp_p)) if exp_p else 0.5
        pval = sum(comb(total, k) * pe**k * (1-pe)**(total-k)
                   for k in range(beats, total + 1)) if total else 1.0
        note = ""
        if pe * total < 1.0:
            note = ("   [no power: the scorer ranks WT so high that almost "
                    "nothing can outrank it]")
        print(f"\n  T1  hypersensitive substitutions ranked ABOVE WT: "
              f"{beats}/{total}   expected by chance {pe*total:.1f}   "
              f"p = {pval:.4f}{note}")
        # T2 permutation on mean rank
        rng = np.random.default_rng(0)
        obs = float(np.mean(ranks)) if ranks else 0.0
        null = []
        pool = {p: [i + 1 for i in range(19)] for p in have}
        cnt = {p: sum(1 for s in HYPER[p] if s in rows[p]["order"]) for p in have}
        for _ in range(20000):
            draw = []
            for p in have:
                draw += list(rng.choice(pool[p], size=cnt[p], replace=False))
            null.append(np.mean(draw))
        null = np.array(null)
        pperm = float((null <= obs).mean())
        print(f"  T2  mean rank of hypersensitive set {obs:.2f} vs null "
              f"{null.mean():.2f}   p = {pperm:.4f}")
        print(f"  T3  positions where the scorer puts WT FIRST: {wt_first}/{len(have)}"
              + ("   <-- DEMONSTRATED FAILURE (WT is not optimal here)"
                 if wt_first else "   (scorer does not claim WT is best)"))
        summary[tag] = dict(beats=beats, total=total, p_T1=pval,
                            mean_rank=obs, p_T2=pperm, wt_first=wt_first,
                            n_pos=len(have))
        print()
    json.dump(summary, open(os.path.join(OUT, "scan_vs_hypersensitive.json"), "w"),
              indent=1)
    print(f"wrote {OUT}/scan_vs_hypersensitive.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
