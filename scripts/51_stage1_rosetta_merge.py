#!/usr/bin/env python
"""
51_stage1_rosetta_merge.py -- score the Rosetta FastDesign arm of stage 1 on the
same footing as the LigandMPNN arm (sections 23f-i).

WHAT IS COMPARED, AND WHY IT IS NOT THE SAME NUMBER
---------------------------------------------------
LigandMPNN exposes a per-position distribution (`sampling_probs`), so its
statistic is a probability read straight off the model. Rosetta exposes only
sampled structures, so the analogue is the FREQUENCY of the true residue across
independent trajectories. These are not interchangeable: a frequency estimated
from 50 trajectories has sampling error a probability does not, and Rosetta's
trajectories are correlated through a shared starting structure. Both are
therefore compared only to their OWN ligand-swap null, never to each other's
absolute value:

    delta = freq_mandi(true residue) - freq_aba(true residue)

Reported with a two-sided Fisher exact test on the 2x2 of (true / not-true)
against (mandipropamid / ABA). With 50 trajectories per arm the smallest
detectable frequency difference is coarse, so a null result here bounds the
effect rather than excluding it, and the p-value is descriptive.

WHAT ELSE THIS PRINTS, AND WHY IT MATTERS MORE THAN RECOVERY
------------------------------------------------------------
Recovery of the four known mutations is the benchmark, but the project's actual
deliverable is LIBRARY POSITIONS. So this also reports, per designable position,
the substitutions Rosetta proposed and how often -- including at positions with
no ground truth. A method that misses K59R but converges hard on some other
position is telling us something usable, and that is invisible in a recall
number. Substitutions are reported against the ABA arm too, so a proposal that
appears regardless of ligand can be discounted.

Interface energies (bound minus separated, ref2015) are summarised per arm.
These rank designs within an arm; they are NOT comparable to the LigandMPNN
arm and not a binding affinity.

Usage:  python 51_stage1_rosetta_merge.py [--indir results/stage1_rosetta]
"""
import argparse, collections, itertools, json, math, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRUE = {59: ("K", "R"), 81: ("V", "I"), 108: ("F", "A"), 159: ("F", "L")}
DESIGN = [59, 81, 83, 92, 94, 108, 110, 120, 122, 141, 159, 160, 163, 164, 167]

ap = argparse.ArgumentParser()
ap.add_argument("--indir", default=os.path.join(ROOT, "results", "stage1_rosetta"))
ap.add_argument("--top", type=int, default=3, help="substitutions to list per position")
args = ap.parse_args()

files = sorted(f for f in os.listdir(args.indir) if f.endswith(".json"))
if not files:
    sys.exit(f"no result files in {args.indir}")

arms = collections.defaultdict(list)
for f in files:
    d = json.load(open(os.path.join(args.indir, f)))
    arms[d["arm"]].extend(d["trajectories"])

print(f"{len(files)} result files, {len(arms)} arms\n")
for a, r in sorted(arms.items()):
    print(f"  {a:>24}: {len(r):>3} trajectories")


def col(rows, pos):
    """Residue at a designable PDB position across trajectories."""
    i = DESIGN.index(pos)
    return [r["seq"][i] for r in rows if len(r["seq"]) > i]


def freq(rows, pos, aa):
    c = col(rows, pos)
    return (c.count(aa) / len(c)) if c else float("nan"), c.count(aa), len(c)


def fisher(a, b, c, d):
    """Two-sided Fisher exact on [[a,b],[c,d]]; returns p."""
    def lc(n, k):
        return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))
    n = a + b + c + d
    p_obs = math.exp(lc(a + b, a) + lc(c + d, c) - lc(n, a + c))
    tot = 0.0
    for i in range(0, min(a + b, a + c) + 1):
        j, k, l = a + b - i, a + c - i, d - a + i
        if j < 0 or k < 0 or l < 0:
            continue
        p = math.exp(lc(a + b, i) + lc(c + d, k) - lc(n, a + c))
        if p <= p_obs * (1 + 1e-9):
            tot += p
    return min(1.0, tot)


print("\n" + "=" * 84)
print("PRIMARY: frequency of the TRUE residue, minus the ABA ligand-swap null")
print("=" * 84)
for alphabet in ("dsm_hao", "free"):
    m, pg, ab = (arms.get(f"wt_mandi__{alphabet}", []),
                 arms.get(f"polygly_mandi__{alphabet}", []),
                 arms.get(f"wt_aba__{alphabet}", []))
    if not (m and ab):
        print(f"\n--- alphabet: {alphabet} --- (incomplete, skipping)")
        continue
    print(f"\n--- alphabet: {alphabet} ---")
    print(f"{'mut':>8} {'f(wt_mandi)':>12} {'f(polygly)':>11} {'f(wt_aba)':>10} "
          f"{'delta_wt':>9} {'delta_pg':>9} {'p(wt)':>8}")
    for pos, (w, t) in sorted(TRUE.items()):
        fm, km, nm = freq(m, pos, t)
        fp, _, _ = freq(pg, pos, t) if pg else (float("nan"), 0, 0)
        fa, ka, na = freq(ab, pos, t)
        p = fisher(km, nm - km, ka, na - ka)
        print(f"{w}{pos}{t:<5} {fm:>12.3f} {fp:>11.3f} {fa:>10.3f} "
              f"{fm-fa:>+9.3f} {fp-fa:>+9.3f} {p:>8.3f}")
    print("  only delta is evidence; a frequency the ABA arm matches is a")
    print("  ligand-independent preference, not recovery")

print("\n" + "=" * 84)
print("RECALL-AT-N  (uniform chance at N=50 is 3.84/4 -- an upper bound, not a target)")
print("=" * 84)
print(f"{'arm':>24} {'n':>4} {'recalled':>9}  {'which'}")
for a, rows in sorted(arms.items()):
    got = [f"{w}{p}{t}" for p, (w, t) in TRUE.items()
           if any(r["seq"][DESIGN.index(p)] == t for r in rows if len(r["seq"]) > DESIGN.index(p))]
    print(f"{a:>24} {len(rows):>4} {len(got)}/4{'':>5}  {','.join(got) or '-'}")

print("\n" + "=" * 84)
print("EFFECTIVE SAMPLE SIZE -- are N trajectories worth N draws?")
print("=" * 84)
print(f"{'arm':>24} {'n':>4} {'distinct':>9} {'mean Hamming':>13} {'mean entropy':>13}")
for a, rows in sorted(arms.items()):
    seqs = [r["seq"] for r in rows]
    if not seqs:
        continue
    pairs = list(itertools.combinations(range(len(seqs)), 2))[:2000]
    ham = (sum(sum(x != y for x, y in zip(seqs[i], seqs[j])) for i, j in pairs)
           / len(pairs)) if pairs else 0.0
    ent = 0.0
    for k in range(len(DESIGN)):
        c = collections.Counter(s[k] for s in seqs if len(s) > k)
        n = sum(c.values())
        ent += -sum((v / n) * math.log2(v / n) for v in c.values()) if n else 0
    print(f"{a:>24} {len(seqs):>4} {len(set(seqs)):>9} {ham:>13.2f} "
          f"{ent/len(DESIGN):>13.2f}")

print("\n" + "=" * 84)
print("SUBSTITUTIONS PROPOSED -- the actual deliverable: candidate library positions")
print("=" * 84)
for alphabet in ("dsm_hao", "free"):
    m, ab = arms.get(f"wt_mandi__{alphabet}", []), arms.get(f"wt_aba__{alphabet}", [])
    if not (m and ab):
        continue
    print(f"\n--- alphabet: {alphabet}   (mandipropamid vs ABA, WT excluded) ---")
    print(f"{'pos':>5} {'wt':>3} {'truth':>6}  {'top substitutions (mandi)':<34} "
          f"{'same in ABA arm'}")
    for pos in DESIGN:
        cm, ca = col(m, pos), col(ab, pos)
        if not cm:
            continue
        wt = None
        for r in m:
            wt = r["seq"][DESIGN.index(pos)] if not r["mutations"].get(str(pos)) else wt
            if wt:
                break
        cnt = collections.Counter(x for x in cm if x != wt)
        cna = collections.Counter(x for x in ca if x != wt)
        top = ", ".join(f"{a}{100*n/len(cm):.0f}%" for a, n in cnt.most_common(args.top))
        same = ", ".join(f"{a}{100*cna[a]/len(ca):.0f}%" for a, _ in cnt.most_common(args.top)
                         if cna.get(a))
        t = f"->{TRUE[pos][1]}" if pos in TRUE else ""
        print(f"{pos:>5} {wt or '?':>3} {t:>6}  {top or '(none)':<34} {same or '-'}")
    print("  a substitution at the same rate in the ABA arm is ligand-independent")

print("\n" + "=" * 84)
print("INTERFACE ENERGY (ref2015 bound - separated; within-arm ranking only)")
print("=" * 84)
print(f"{'arm':>24} {'n':>4} {'median':>9} {'best':>9} {'worst':>9}")
for a, rows in sorted(arms.items()):
    v = sorted(r["interface_dG"] for r in rows if r.get("interface_dG") is not None)
    if not v:
        continue
    med = v[len(v) // 2] if len(v) % 2 else 0.5 * (v[len(v)//2 - 1] + v[len(v)//2])
    print(f"{a:>24} {len(v):>4} {med:>9.2f} {v[0]:>9.2f} {v[-1]:>9.2f}")

out = os.path.join(args.indir, "stage1_rosetta_summary.json")
json.dump({a: r for a, r in arms.items()}, open(out, "w"), indent=1)
print(f"\nwrote {out}")
