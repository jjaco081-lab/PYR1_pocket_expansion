#!/usr/bin/env python
r"""
116_secondary_library.py -- design a coumarin-class SECONDARY library from pooled
round-1 hits, and score it against the round-2 hits Tian actually found.

THE TEST, as specified in README 54c
------------------------------------
Round-2 (sd07) is the answer sheet and is NOT an input. The inputs are:

  * sd03  -- 692 round-1 clones over 194 ligands, which is what Tian had before
             building the coumarin library
  * the target ligands' SMILES, which any prospective user has by definition

and the question is whether a library of at most Tian's 138,240 members, designed
from those inputs alone, contains the round-2 sensors.

WHY THIS IS THE RIGHT PROBLEM (README 53b, 54b)
-----------------------------------------------
The libraries are substitution-DEPTH limited, not full combinatorial. Round 1 is
broad and shallow (18 positions, 2-3 substitutions per clone); round 2 is narrow
and DEEP (11 positions, 5-8 per clone). So the design task is not "make it
smaller" -- it is which positions and residues are worth combining deeply.

Measured here before anything was designed: all 23 of Tian's coumarin
substitutions are already in round-1's 144-substitution vocabulary, and the only
round-2 substitutions missing from round-1 are 8 singletons at positions outside
the 18 (C65F, F71S, R74C, S109R, T124M, R134L, M178I, D184Y) -- PCR artefacts
carried through, not designed content. So the vocabulary is CLOSED and this is
subset selection, exactly as 54b argued.

THE METRIC (README 51c)
-----------------------
Hits are positive-unlabeled: a clone is a hit that was FOUND, not the best
sequence (feedback_hits_are_not_optima). So recall only, never precision, and the
headline is LIGAND-level:

    fraction of the 11 target ligands for which the library contains at least one
    known round-2 clone -- where a clone counts only if EVERY one of its
    substitutions is in the menu.

Clone-level capture is reported alongside but is the wrong row to read: you need
A sensor, not every sensor.

CEILINGS, so no number here is read as better than it is
    11/11 ligands and 70/78 clones are reachable within the 18 round-1 positions.
    Tian's own coumarin library captures 11/11 and 69/78 at 138,240 members --
    and that is near-automatic, because sd07 IS that library's output. The
    incumbent is therefore a CONTAINMENT ceiling, not a competitor: the only way
    to win is to reach comparable ligand capture in a materially smaller library.

LIBRARY SIZE
    size = product over positions of (1 + n_substitutions_at_that_position)
    WT is retained at every position, matching sd04's arithmetic (the Coumarin
    row multiplies out to 1.382e5 only if WT is counted).
    Adding a residue at a position already carrying k costs log((k+2)/(k+1));
    opening a NEW position costs log 2. Positions are expensive, depth is cheap,
    which is why the greedy below spends on depth.
"""
import argparse
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_xlsx import table                                     # noqa: E402

ROOT = os.path.dirname(HERE)
SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "secondary_library")

AA = "ACDEFGHIKLMNPQRSTVWY"

#: SMILES for the 11 coumarin-class targets. Eight come from sd03's own
#: canonical_smiles column; the three with no round-1 clone at all
#: (4-methylumbelliferone, 5,7-dihydroxy-4-methylcoumarin, 7-methoxycoumarin) are
#: public structures written out here. That is not leakage -- sd07 contains no
#: SMILES column, so nothing about the ANSWER is being read; a prospective user
#: knows their target's structure by definition.
EXTRA_SMILES = {
    "4-methylumbelliferone": "CC1=CC(=O)Oc2cc(O)ccc12",
    "5,7-dihydroxy-4-methylcoumarin": "CC1=CC(=O)Oc2cc(O)cc(O)c12",
    "7-methoxycoumarin": "COc1ccc2ccc(=O)oc2c1",
}


def load(path, ligcol, smicol=None):
    hdr, recs = table(path)
    poscols = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    out = []
    for r in recs:
        subs = []
        for c in poscols:
            v = (r.get(c) or "").strip().upper()
            if len(v) == 1 and v in AA and v != c[0]:
                subs.append(c + v)
        out.append({"lig": (r.get(ligcol) or "?").strip(),
                    "subs": subs,
                    "smiles": (r.get(smicol) or "").strip() if smicol else "",
                    "min_conc": (r.get("min_conc") or "").strip()})
    return out, poscols


def tian_menu(name):
    """One library definition from sd04. Used ONLY as the incumbent yardstick."""
    _, recs = table(f"{SD}/pnas.2519924122.sd04.xlsx")
    m = defaultdict(set)
    for r in recs:
        if (r.get("Library") or "").strip() == name:
            m[f"{r['WT']}{int(float(r['Position']))}"] |= set(
                r["Amino Acids allowed for mutation"].strip())
    return dict(m)


def size_of(menu):
    n = 1
    for p, aas in menu.items():
        n *= (1 + len(aas))
    return n


def capture(menu, round2):
    """Ligand-level and clone-level containment. A clone counts only if all of its
    substitutions are in the menu -- an unmatched substitution means that exact
    sequence cannot be assembled."""
    flat = {p + a for p, aas in menu.items() for a in aas}
    byl = defaultdict(lambda: [0, 0])
    for c in round2:
        byl[c["lig"]][0] += 1
        if all(s in flat for s in c["subs"]):
            byl[c["lig"]][1] += 1
    ligs = sum(1 for v in byl.values() if v[1] > 0)
    clones = sum(v[1] for v in byl.values())
    return ligs, len(byl), clones, sum(v[0] for v in byl.values()), dict(byl)



def oracle(round2, positions):
    """EXACT smallest library containing at least one round-2 clone for every
    ligand, by branch and bound.

    Not a method -- it reads the answer sheet. Its only job is to bound the prize:
    the ratio Tian_size / oracle_size is the entire headroom on the size axis, so
    it decides up front whether a 10x reduction is even available.

    ⚠ A GREEDY version of this was tried first and was wrong twice over. It
    returned 31,104 against the true optimum of 9,216 -- so it understated the
    headroom by 3.4x, which would have made the whole exercise look not worth
    doing -- and because it iterated a SET of ligand names, PYTHONHASHSEED made it
    return a different answer on different runs. Both are fixed here: the search
    is exact, and every iteration order is sorted.

    Branch and bound is valid because library size is monotone non-decreasing as
    residues are added, so a partial assignment already at or above the incumbent
    can be cut. Dominated clones (a superset of another clone for the same ligand)
    are dropped first -- they can never be the cheaper choice."""
    byl = defaultdict(set)
    for c in round2:
        if c["subs"] and all(s[:-1] in positions for s in c["subs"]):
            byl[c["lig"]].add(frozenset(c["subs"]))
    opts = {}
    for l, S in byl.items():
        opts[l] = sorted((x for x in S if not any(t < x for t in S)),
                         key=lambda z: (len(z), sorted(z)))
    ligs = sorted(opts, key=lambda l: (len(opts[l]), l))     # fewest options first
    best = [float("inf"), None]

    def rec(i, menu):
        sz = size_of({p: v for p, v in menu.items() if v})
        if sz >= best[0]:
            return
        if i == len(ligs):
            best[0] = sz
            best[1] = {p: set(v) for p, v in menu.items() if v}
            return
        for subs in opts[ligs[i]]:
            added = []
            for x in sorted(subs):
                p, a = x[:-1], x[-1]
                if a not in menu.setdefault(p, set()):
                    menu[p].add(a)
                    added.append((p, a))
            rec(i + 1, menu)
            for p, a in added:
                menu[p].discard(a)

    rec(0, {})
    return best[1], best[0], len(ligs)


def oracle_k(round2, positions, k):
    """Exact smallest library covering the EASIEST k of the ligands. Same search,
    minimised over all C(n,k) subsets -- so the size axis can be read at partial
    capture, where the reductions are much larger."""
    import itertools
    byl = defaultdict(set)
    for c in round2:
        if c["subs"] and all(s[:-1] in positions for s in c["subs"]):
            byl[c["lig"]].add(frozenset(c["subs"]))
    opts = {l: sorted((x for x in S if not any(t < x for t in S)),
                      key=lambda z: (len(z), sorted(z))) for l, S in byl.items()}
    globalbest = [float("inf")]
    for combo in itertools.combinations(sorted(opts), k):
        ligs = sorted(combo, key=lambda l: (len(opts[l]), l))

        def rec(i, menu):
            sz = size_of({p: v for p, v in menu.items() if v})
            if sz >= globalbest[0]:
                return
            if i == len(ligs):
                globalbest[0] = sz
                return
            for subs in opts[ligs[i]]:
                added = []
                for x in sorted(subs):
                    p, a = x[:-1], x[-1]
                    if a not in menu.setdefault(p, set()):
                        menu[p].add(a)
                        added.append((p, a))
                rec(i + 1, menu)
                for p, a in added:
                    menu[p].discard(a)

        rec(0, {})
    return globalbest[0]


# ---------------------------------------------------------------------------
# LIGAND SIMILARITY -- the only ligand-aware input
# ---------------------------------------------------------------------------
def similarity(r1, targets, extra):
    """max Morgan/ECFP4 Tanimoto from each round-1 ligand to the target class."""
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    def fp(smi):
        m = Chem.MolFromSmiles(smi)
        return gen.GetFingerprint(m) if m else None

    smi = {}
    for c in r1:
        if c["smiles"]:
            smi.setdefault(c["lig"].lower(), c["smiles"])
    tfp = []
    for t in targets:
        s = smi.get(t.lower()) or extra.get(t.lower())
        if not s:
            raise SystemExit(f"no SMILES for target {t}")
        f = fp(s)
        if f is None:
            raise SystemExit(f"unparseable SMILES for {t}")
        tfp.append(f)
    if len(tfp) != len(targets):
        raise SystemExit("target fingerprint count mismatch")
    out = {}
    for lig, s in smi.items():
        f = fp(s)
        if f is not None:
            out[lig] = max(DataStructs.TanimotoSimilarity(f, t) for t in tfp)
    return out


def clone_weights(r1, sim, alpha, potency=False):
    """One weight per round-1 clone: class similarity sharpened by alpha, times
    potency if asked. alpha = 0 gives the ligand-blind baseline."""
    pw = {"1.0": 3.0, "10.0": 2.0, "100.0": 1.0}
    out = []
    for c in r1:
        w = 1.0 if alpha == 0 else sim.get(c["lig"].lower(), 0.0) ** alpha
        if potency:
            w *= pw.get(c["min_conc"], 1.0)
        out.append(w)
    return out


# ---------------------------------------------------------------------------
# DESIGNERS -- each returns a menu {position: set(residues)} within the budget
# ---------------------------------------------------------------------------
def _cost(k):
    """log-cost of adding one residue to a position that already carries k.
    Opening a position (k=0) costs log 2; the second residue log 3/2; the third
    log 4/3. Depth is cheaper than breadth, but only by ~1.7x -- which is why the
    first version of this script, ranking on raw weight over this cost, bought
    eleven residues at six positions and captured 1 of 11 ligands."""
    return math.log((k + 2) / (k + 1))


def design_global(subw, budget, positions):
    """Greedy on raw substitution weight per unit log-size.

    KEPT AS A NEGATIVE CONTROL. It fails, and the failure is instructive: pooled
    counts are dominated by a few positions (F159 alone carries 317 of 1747
    round-1 substitutions), so it spends the whole budget deepening them. Round-2
    clones use 7-8 DISTINCT positions each, and a clone is lost if even one of its
    positions is shut, so breadth is what containment actually needs."""
    menu = defaultdict(set)
    order = []
    while True:
        best, bv = None, 0.0
        for s, w in subw.items():
            p, a = s[:-1], s[-1]
            if p not in positions or w <= 0 or a in menu[p]:
                continue
            if size_of({**{q: v for q, v in menu.items() if v},
                        p: menu[p] | {a}}) > budget:
                continue
            v = w / _cost(len(menu[p]))
            if v > bv:
                best, bv = s, v
        if best is None:
            break
        menu[best[:-1]].add(best[-1])
        order.append(best)
    return {p: v for p, v in menu.items() if v}, order


def design_profile(subw, budget, positions):
    """Per-POSITION normalisation, which is the fix.

    A substitution scores  W(position) * share(residue | position).  The best
    residue at a closed position therefore scores about W(position), which is
    comparable across positions and far above the third or fourth residue at an
    open one. Breadth wins until the good positions are all open, then depth
    accumulates -- the shape round 2 actually has (11 positions, ~2 residues each).

    This is also, in effect, what Tian does: build a residue profile per position
    from pooled round-1 hits and take the top of each."""
    posw = Counter()
    share = defaultdict(Counter)
    for s, w in subw.items():
        p = s[:-1]
        if p in positions:
            posw[p] += w
            share[p][s[-1]] += w
    score = {}
    for p, resid in share.items():
        tot = sum(resid.values()) or 1.0
        for a, w in resid.items():
            score[p + a] = posw[p] * (w / tot)
    return design_global(score, budget, positions)


def design_cover(r1, cw, budget, positions):
    """Maximise weighted coverage of round-1 CLONES, not substitutions.

    The metric is set containment -- a clone counts only if every substitution is
    present -- so the objective is made to have the same shape. Each round-1 clone
    contributes its weight times (fraction of its substitutions present)^2; the
    square is what gives a gradient before any clone is complete and what rewards
    finishing one that is nearly there. Round-1 clones are 2-3 deep rather than
    7-8, so this is a proxy for the round-2 shape, not a model of it."""
    clones = [(set(c["subs"]), w) for c, w in zip(r1, cw)
              if c["subs"] and w > 0 and all(s[:-1] in positions for s in c["subs"])]
    menu = defaultdict(set)
    order = []
    flat = set()

    def obj(extra):
        f = flat | {extra}
        t = 0.0
        for subs, w in clones:
            t += w * (len(subs & f) / len(subs)) ** 2
        return t

    base = sum(w * (len(subs & flat) / len(subs)) ** 2 for subs, w in clones)
    while True:
        best, bv, bo = None, 0.0, None
        for s in {x for subs, _ in clones for x in subs} - flat:
            p = s[:-1]
            if size_of({**{q: v for q, v in menu.items() if v},
                        p: menu[p] | {s[-1]}}) > budget:
                continue
            o = obj(s)
            v = (o - base) / _cost(len(menu[p]))
            if v > bv:
                best, bv, bo = s, v, o
        if best is None:
            break
        menu[best[:-1]].add(best[-1])
        flat.add(best)
        order.append(best)
        base = bo
    return {p: v for p, v in menu.items() if v}, order



def design_profile_depth(subw, budget, positions, npos, frac):
    """The designer that survives: a per-position PROFILE with variable depth.

    Two knobs, both structural rather than fitted:
      npos  how many positions to open, ranked by total class weight
      frac  a residue is kept if its weight is at least `frac` of the best
            residue at that position -- so a position with one dominant answer
            contributes one residue and an ambiguous one contributes several

    Depth is allocated by CONFIDENCE, not by a global ranking. That matters
    because pooled counts are wildly uneven across positions (F159 alone carries
    317 of 1747 round-1 substitutions), and a global ranking spends the whole
    budget deepening the loud positions while leaving quiet ones shut. Round-2
    clones use 7-8 distinct positions each and are lost if even one is shut."""
    posw = Counter()
    share = defaultdict(Counter)
    for s, w in subw.items():
        p = s[:-1]
        if p in positions:
            posw[p] += w
            share[p][s[-1]] += w
    chosen = [p for p, _ in posw.most_common(npos)]
    menu = {}
    for p in chosen:
        top = share[p].most_common(1)[0][1]
        menu[p] = {a for a, w in share[p].items() if w >= frac * top}
    while size_of(menu) > budget:
        # drop the weakest residue at the position that is currently deepest
        cand = [(share[p][a], p, a) for p, v in menu.items() for a in v
                if len(v) > 1]
        if not cand:
            break
        _, p, a = min(cand)
        menu[p].discard(a)
    return {p: v for p, v in menu.items() if v}, []


def subw_from(r1x, sim, alpha, potency=False):
    cw = clone_weights(r1x, sim, alpha, potency)
    w = Counter()
    for c, x in zip(r1x, cw):
        for s in c["subs"]:
            w[s] += x
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget", type=float, default=138240)
    ap.add_argument("--alpha", type=float, default=4.0,
                    help="sharpening exponent on ligand similarity")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    r1, pos18 = load(f"{SD}/pnas.2519924122.sd03(1).xlsx", "library_name",
                     "canonical_smiles")
    r2, _ = load(f"{SD}/pnas.2519924122.sd07(1).xlsx", "compound")
    targets = sorted({c["lig"] for c in r2})
    P18 = set(pos18)
    sim = similarity(r1, targets, EXTRA_SMILES)
    tm = tian_menu("Coumarin")
    TIAN_SUBS = {p + x for p, v in tm.items() for x in v}
    nown = {t: sum(1 for c in r1 if c["lig"].lower() == t.lower()) for t in targets}

    say("=" * 80)
    say("1. INPUTS   (round-2 sd07 is the answer sheet and is never an input)")
    say("=" * 80)
    say(f"   round-1 sd03  {len(r1)} clones / {len({c['lig'] for c in r1})} ligands / "
        f"{len({s for c in r1 for s in c['subs']})} distinct substitutions / "
        f"{len(P18)} positions")
    say(f"   targets       {len(targets)} coumarin-class ligands, "
        f"{sum(1 for t in targets if nown[t])} with round-1 clones of their own "
        f"({sum(1 for t in targets if not nown[t])} had none at all)")
    say(f"   budget        {a.budget:,.0f} members = Tian's coumarin library")

    say("")
    say("=" * 80)
    say("2. THE VOCABULARY IS CLOSED -- this is subset selection, not invention")
    say("=" * 80)
    V1 = Counter(s for c in r1 for s in c["subs"])
    inv = [s for s in sorted(TIAN_SUBS) if s in V1]
    say(f"   {len(inv)}/{len(TIAN_SUBS)} of Tian's coumarin substitutions already "
        f"appear in the pooled round-1 vocabulary")
    V2 = Counter(s for c in r2 for s in c["subs"])
    out = sorted(s for s in V2 if s not in V1)
    say(f"   round-2 substitutions absent from round-1: {len(out)} of {len(V2)}, "
        f"all singletons outside the 18 positions")
    say(f"     {', '.join(out)}")
    say("   ==> nothing has to be invented. Every residue the answer uses was")
    say("       already visible in round 1; the whole task is choosing a subset.")

    say("")
    say("=" * 80)
    say("3. CEILINGS -- read every number below against these")
    say("=" * 80)
    full = defaultdict(set)
    for c in r1:
        for s in c["subs"]:
            full[s[:-1]].add(s[-1])
    lg, nl, cl, nc, _ = capture(full, r2)
    say(f"   whole round-1 vocabulary  {lg}/{nl} ligands  {cl}/{nc} clones  "
        f"size {size_of(full):.2e}")
    lg, nl, cl, nc, _ = capture(tm, r2)
    say(f"   Tian's coumarin library   {lg}/{nl} ligands  {cl}/{nc} clones  "
        f"size {size_of(tm):,}  ({len(tm)} pos, {len(TIAN_SUBS)} subs)")
    om, osz, nlig = oracle(r2, P18)
    say(f"   EXACT ORACLE              {nlig}/{nlig} ligands            "
        f"size {osz:,}  ({len(om)} pos, "
        f"{sum(len(v) for v in om.values())} subs)")
    say(f"     menu  " + " ".join(f"{p}:{''.join(sorted(v))}" for p, v in
                                  sorted(om.items(), key=lambda x: int(x[0][1:]))))
    say(f"   ==> headroom at full capture is {size_of(tm)/osz:.1f}x. A 10x "
        f"reduction IS available in principle.")
    say("   exact minimum for the easiest k ligands:")
    okc = {}
    for k in (8, 9, 10, 11):
        okc[k] = oracle_k(r2, P18, k)
        say(f"     {k:2d}/11 ligands  {okc[k]:>8,}   "
            f"({size_of(tm)/okc[k]:5.1f}x smaller than Tian)")
    say("   NOTE the shape of the optimum: 11 positions but only "
        f"{sum(len(v) for v in om.values())} substitutions, i.e. "
        f"{sum(1 for v in om.values() if len(v)==1)} of 11 positions need exactly")
    say("   ONE residue. Positions are non-negotiable (every round-2 clone spans")
    say("   7-8 of them); the entire prize is in residue identity.")

    say("")
    say("=" * 80)
    say("4. CAN THE POSITIONS BE RECOVERED?   (yes)")
    say("=" * 80)
    need = Counter(s[:-1] for c in r2 for s in c["subs"])
    say(f"   {'position':<10}{'used in r2':>11}{'rank, blind':>13}{'rank, class':>13}")
    ranks = {}
    for alpha, key in ((0, "blind"), (a.alpha, "class")):
        w = Counter()
        for s, x in subw_from(r1, sim, alpha).items():
            w[s[:-1]] += x
        ranks[key] = {p: i + 1 for i, (p, _) in enumerate(w.most_common())}
    for p, n in need.most_common():
        if p not in P18:
            continue
        star = "*" if p in tm else " "
        say(f"  {star}{p:<9}{n:>11}{ranks['blind'].get(p,'-'):>13}"
            f"{ranks['class'].get(p,'-'):>13}")
    for key in ("blind", "class"):
        top11 = {p for p, r in ranks[key].items() if r <= 11}
        top12 = {p for p, r in ranks[key].items() if r <= 12}
        say(f"   {key:<6}: {len(top11 & set(tm))}/11 true positions in its top 11, "
            f"{len(top12 & set(tm))}/11 in its top 12")
    say("   (* = one of the 11 positions Tian actually used)")

    say("")
    say("=" * 80)
    say("5. CAN THE RESIDUES BE RECOVERED?   (partly -- and this is the bottleneck)")
    say("=" * 80)
    say("   Rank of the OPTIMAL residue inside its own position's profile:")
    say(f"   {'pos':<7}{'optimal':<9}{'rank blind':>12}{'rank class':>12}")
    agg = {}
    for key, alpha in (("blind", 0), ("class", a.alpha)):
        per = defaultdict(Counter)
        for s, x in subw_from(r1, sim, alpha).items():
            per[s[:-1]][s[-1]] += x
        agg[key] = per
    t1 = {"blind": 0, "class": 0}
    t3 = {"blind": 0, "class": 0}
    for p in sorted(om, key=lambda x: int(x[1:])):
        row = {}
        for key in ("blind", "class"):
            order = [x for x, _ in agg[key][p].most_common()]
            rk = [order.index(x) + 1 for x in om[p] if x in order]
            row[key] = min(rk) if rk else 99
            t1[key] += row[key] == 1
            t3[key] += row[key] <= 3
        say(f"   {p:<7}{''.join(sorted(om[p])):<9}{row['blind']:>12}{row['class']:>12}")
    for key in ("blind", "class"):
        say(f"   {key:<6}: optimal residue ranked 1st at {t1[key]}/11 positions, "
            f"top-3 at {t3[key]}/11")
    say("   ==> ligand-awareness moves residue identity a long way "
        f"({t1['blind']}/11 -> {t1['class']}/11 exact hits), which is the")
    say("       first place in this project where knowing the ligand has paid.")

    say("")
    say("=" * 80)
    say("6. DESIGNED LIBRARIES, re-optimised at each budget")
    say("=" * 80)
    marks = [3e3, 1e4, 3e4, 1e5, 1.3824e5]
    sw_blind = subw_from(r1, sim, 0)
    sw_class = subw_from(r1, sim, a.alpha)
    designs = {}

    def run(name, fn):
        row = []
        for b in marks:
            menu = fn(b)
            lg, _, cl, _, _ = capture(menu, r2)
            flat = {p + x for p, v in menu.items() for x in v}
            row.append((lg, cl, size_of(menu), len(flat & TIAN_SUBS)))
        designs[name] = row
        say(f"   {name:<20}" + "".join(f"{f'{l}/{c}':>12}" for l, c, _s, _r in row))

    say(f"   {'design':<20}" + "".join(f"{m:>12,.0f}" for m in marks))
    say("   " + "-" * 76)
    run("blind/global", lambda b: design_global(sw_blind, b, P18)[0])
    run("blind/profile", lambda b: design_profile_depth(sw_blind, b, P18, 11, 0.35)[0])
    run("class/global", lambda b: design_global(sw_class, b, P18)[0])
    run("class/profile", lambda b: design_profile_depth(sw_class, b, P18, 11, 0.35)[0])
    run("class/cover",
        lambda b: design_cover(r1, clone_weights(r1, sim, a.alpha), b, P18)[0])
    say("   " + "-" * 76)
    say(f"   {'ORACLE (exact)':<20}"
        + "".join(f"{('11/-' if m>=okc[11] else max([k for k in okc if okc[k]<=m], default=0)):>12}"
                  for m in marks))
    say(f"   {'Tian (138,240)':<20}{'':>36}{'':>12}{'11/69':>12}")
    say("")
    say("   ⚠ THE COMPARISON WITH TIAN IS NOT SYMMETRIC, and this bounds what the")
    say("     benchmark can ever show. sd07 clones were DRAWN FROM Tian's library,")
    say("     so it contains them by construction -- 11/11 is automatic, not earned.")
    say("     A different library of equal quality would score badly here because")

    say("     the sensors IT would have found were never screened. Hits are")
    say("     positive-unlabeled (feedback_hits_are_not_optima), so this metric")
    say("     measures REDISCOVERY OF TIAN'S CHOICES, not library quality.")

    say("")
    say("=" * 80)
    say("7. THE SYMMETRIC METRIC: recovery of Tian's 23 substitutions")
    say("=" * 80)
    say("   Tian chose their 23 from the same round-1 data we are given, so this")
    say("   is method against method with no drawn-from advantage.")
    say(f"   {'design':<20}" + "".join(f"{m:>12,.0f}" for m in marks))
    for name, row in designs.items():
        say(f"   {name:<20}" + "".join(f"{f'{r}/23':>12}" for _l, _c, _s, r in row))

    say("")
    say("=" * 80)
    say("8. LEAVE-ONE-LIGAND-OUT (class/profile) -- no memorising own round-1")
    say("=" * 80)
    say(f"   {'ligand':<34}{'own r1':>8}{'captured':>10}{'Tian subs':>11}")
    won = 0
    lolo = {}
    for t in targets:
        keep = [c for c in r1 if c["lig"].lower() != t.lower()]
        menu = design_profile_depth(subw_from(keep, sim, a.alpha),
                                    a.budget, P18, 11, 0.35)[0]
        flat = {p + x for p, v in menu.items() for x in v}
        ok = any(all(s in flat for s in c["subs"])
                 for c in r2 if c["lig"] == t and c["subs"])
        won += ok
        lolo[t] = {"own_r1": nown[t], "captured": bool(ok), "size": size_of(menu)}
        say(f"   {t:<34}{nown[t]:>8}{'YES' if ok else 'no':>10}"
            f"{len(flat & TIAN_SUBS):>8}/23")
    say(f"   ==> {won}/{len(targets)} ligands captured with their own round-1 "
        f"evidence deleted")

    with open(os.path.join(OUT, "design.json"), "w") as fh:
        json.dump({"budget": a.budget, "alpha": a.alpha,
                   "oracle": {"size": osz, "by_k": okc,
                              "menu": {p: "".join(sorted(v)) for p, v in om.items()}},
                   "tian_size": size_of(tm),
                   "designs": {k: [list(x) for x in v] for k, v in designs.items()},
                   "lolo": lolo}, fh, indent=1)
    with open(os.path.join(OUT, "design.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    say(f"\n   written to {OUT}/design.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
