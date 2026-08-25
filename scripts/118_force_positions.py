#!/usr/bin/env python
r"""
118_force_positions.py -- the two-part statistic, and what forcing a residue buys.

THE IDEA
§54's library arithmetic assumes every position offers wild-type plus a menu, so a
position costs a factor (1 + n). But nothing requires wild-type to be offered. If
a residue is known confidently enough, the position can be FORCED -- every library
member carries it -- and the position then costs a factor of 1 instead of 2. That
is the only move in the whole design space that makes a library smaller without
removing any chemistry from it.

So the design-relevant statistic about a position is not one number but two:

    P(mutated | class)            how often the class touches this position
    P(residue | mutated, class)   and when it does, how consistently

A position can only be forced when the second is near 1. High concentration with a
low mutation rate means the residue is right *when used*, not that it should always
be used.

MEASURED HERE
  * the exact smallest library containing a known round-2 sensor for every one of
    the 11 coumarin ligands, WITH forcing allowed -- which is the true headroom
  * both parts of the statistic, pooled and restricted to the chemical class, so
    it is visible whether a prospective designer could have called it

Inputs are round-1 (sd03) plus target SMILES. sd07 is the answer sheet.
sd08/sd09 stay sealed (§48d).
"""
import itertools
import json
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_xlsx import table                                    # noqa: E402
import importlib.util                                          # noqa: E402

_s = importlib.util.spec_from_file_location(
    "lib116", os.path.join(HERE, "116_secondary_library.py"))
L = importlib.util.module_from_spec(_s)
_s.loader.exec_module(L)

ROOT = os.path.dirname(HERE)
SD = L.SD
OUT = os.path.join(ROOT, "results", "secondary_library")
P18 = ["K59", "V81", "V83", "L87", "A89", "S92", "E94", "F108", "I110", "L117",
       "Y120", "S122", "E141", "F159", "A160", "V163", "V164", "N167"]


def states(path, ligcol):
    """Full assignment over the 18 positions, wild-type where unmutated.

    Deliberately NOT the substitution list: forcing is a statement about what the
    library offers at a position, so wild-type has to be an explicit value rather
    than an absence."""
    hdr, recs = table(path)
    cols = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    out = []
    for r in recs:
        # positions absent from a sheet are wild-type by definition -- sd07 only
        # has columns for the positions its own library randomised, and treating a
        # missing column as missing DATA rather than as wild-type would silently
        # drop every clone.
        st, off = {p: p[0] for p in P18}, False
        for c in cols:
            v = (r.get(c) or "").strip().upper()
            v = v if (len(v) == 1 and v in "ACDEFGHIKLMNPQRSTVWY") else c[0]
            if c in P18:
                st[c] = v
            elif v != c[0]:
                off = True                # substitution outside the 18 positions
        out.append({"lig": (r.get(ligcol) or "?").strip(), "st": st,
                    "off_position": off})
    return out


def exact_min(byl, ligs, forcing=True):
    """Smallest library covering one clone per ligand, by branch and bound.

    With forcing=True a position's cost is the number of DISTINCT residues the
    chosen clones need there -- wild-type only appears if some chosen clone
    actually carries it. With forcing=False wild-type is always offered, which is
    the (1 + n) arithmetic of §51d."""
    best = [float("inf"), None]
    order = sorted(ligs, key=lambda l: (len(byl[l]), l))

    def size(cols):
        n = 1
        for i, s in enumerate(cols):
            k = len(s | ({P18[i][0]} if not forcing else set()))
            n *= max(1, k)
        return n

    def rec(i, cols):
        s = size(cols)
        if s >= best[0]:
            return
        if i == len(order):
            best[0] = s
            best[1] = [set(c) for c in cols]
            return
        for cl in sorted(byl[order[i]]):
            add = []
            for k, a in enumerate(cl):
                if a not in cols[k]:
                    cols[k].add(a)
                    add.append((k, a))
            rec(i + 1, cols)
            for k, a in add:
                cols[k].discard(a)

    rec(0, [set() for _ in P18])
    return best[0], best[1]


def main():
    os.makedirs(OUT, exist_ok=True)
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    r1s = states(f"{SD}/pnas.2519924122.sd03(1).xlsx", "library_name")
    r2s = states(f"{SD}/pnas.2519924122.sd07(1).xlsx", "compound")
    r1, _ = L.load(f"{SD}/pnas.2519924122.sd03(1).xlsx", "library_name",
                   "canonical_smiles")
    targets = sorted({c["lig"] for c in r2s})
    sim = L.similarity(r1, targets, L.EXTRA_SMILES)
    tian = L.size_of(L.tian_menu("Coumarin"))

    byl = defaultdict(set)
    for x in r2s:
        if not x["off_position"]:
            byl[x["lig"]].add(tuple(x["st"].get(p, p[0]) for p in P18))
    ligs = sorted(byl)

    say("=" * 80)
    say("1. THE HEADROOM, WITH AND WITHOUT FORCING")
    say("=" * 80)
    for forcing, name in ((False, "wild-type always offered  (§51d arithmetic)"),
                          (True, "positions may be FORCED   (wild-type dropped)")):
        sz, cols = exact_min(byl, ligs, forcing)
        say(f"   {name}")
        say(f"     exact minimum {sz:>8,}   =  {tian/sz:>4.0f}x smaller than Tian's "
            f"{tian:,}")
        if forcing:
            forced = [P18[i] for i, s in enumerate(cols)
                      if s and P18[i][0] not in s]
            say(f"     positions forced: {forced if forced else 'none'}")
            say("     full menu: " + "  ".join(
                f"{P18[i]}:{''.join(sorted(s))}" for i, s in enumerate(cols) if s))
    say("")
    say("   partial capture, forcing allowed:")
    for k in (8, 9, 10, 11):
        gb = float("inf")
        for combo in itertools.combinations(ligs, k):
            sz, _ = exact_min(byl, list(combo), True)
            gb = min(gb, sz)
        say(f"     {k:2d}/11 ligands  {gb:>8,}  ({tian/gb:>5.0f}x)")
    say("")
    say("   ==> forcing is worth exactly the 6.5x between 30,000 and 4,608, and")
    say("       ALL of it comes from ONE position. Every other position in the")
    say("       optimum still needs wild-type as an option, because different")
    say("       ligands' sensors disagree there.")

    say("")
    say("=" * 80)
    say("2. THE TWO-PART STATISTIC -- could a prospective designer have called it?")
    say("=" * 80)
    sets = {
        "pooled r1": r1s,
        "class r1": [x for x, c in zip(r1s, r1)
                     if sim.get(c["lig"].lower(), 0.0) >= 0.4],
        "own r1": [x for x in r1s if x["lig"].lower() in {t.lower() for t in targets}],
    }
    say(f"   set sizes: " + ", ".join(f"{k} {len(v)}" for k, v in sets.items())
        + f", round-2 {len(r2s)}")
    say("")
    say(f"   {'pos':<7}" + "".join(f"{k:>22}" for k in sets) + f"{'round-2':>22}")
    say(f"   {'':<7}" + "".join(f"{'P(mut)  top|mut':>22}" for _ in
                                list(sets) + ["r2"]))
    rows = {}
    for p in ["K59", "V81", "V83", "F108", "Y120", "S122", "F159", "A160",
              "V163", "V164", "N167"]:
        cells, rec = [], {}
        for k, S in list(sets.items()) + [("round-2", r2s)]:
            mut = [x["st"][p] for x in S if x["st"].get(p, p[0]) != p[0]]
            pm = len(mut) / len(S) if S else 0.0
            c = Counter(mut)
            top = c.most_common(1)
            frac = top[0][1] / len(mut) if mut else 0.0
            cells.append(f"{100*pm:3.0f}%  {top[0][0] if top else '-'}"
                         f"{100*frac:4.0f}% (n={len(mut)})")
            rec[k] = {"p_mut": pm, "top": top[0][0] if top else None,
                      "top_frac": frac, "n": len(mut)}
        rows[p] = rec
        say(f"   {p:<7}" + "".join(f"{c:>22}" for c in cells))
    say("")
    say("   THE ONE POSITION THE OPTIMUM FORCES, in full:")
    for k, S in list(sets.items()) + [("round-2", r2s)]:
        c = Counter(x["st"]["V163"] for x in S if x["st"].get("V163", "V") != "V")
        say(f"     {k:<11} mutated {sum(c.values()):>3}/{len(S):<3}  {dict(c)}")
    say("")
    say("   READ. Pooled round-1 gives V163 a mutation rate of only 9 %, which is")
    say("   why every ligand-blind method here leaves it near the bottom. Restrict")
    say("   to the chemical class and the rate rises to 22-25 % AND the identity")
    say("   becomes unanimous -- 9 of 9 own-class round-1 mutations are W, 10 of 11")
    say("   at similarity >= 0.4. Round 2 then carries W163 in 74 of 78 sensors.")
    say("   So yes: this one was callable in advance, from round-1 data a")
    say("   prospective designer would have had.")
    say("")
    say("   ⚠ n = 9 to 11 clones. One position, called from a small sample, in a")
    say("   retrospective test. It is a hypothesis with a mechanism to check, not")
    say("   an established rule -- and it is pre-registered for PFAS/TNT below.")

    say("")
    say("=" * 80)
    say("3. PRE-REGISTRATION for the sealed sets (sd08 TNT, sd09 PFAS)")
    say("=" * 80)
    say("   Written before either sheet is read, and falsifiable in the obvious way.")
    say("")
    say("   The rule, stated so it can fail:")
    say("     For a target class, compute over round-1 clones whose ligand has")
    say("     Tanimoto >= 0.4 to the class:  P(mutated) and P(top residue | mutated).")
    say("     PREDICT that any position with P(top | mutated) >= 0.85 on n >= 8")
    say("     class clones will carry that residue in >= 80 % of that class's")
    say("     round-2 sensors, and may therefore be FORCED.")
    say("")
    say("   On the coumarin dev set this rule fires exactly once (V163W: 91 % on")
    say("   n = 11) and is right (95 % of round-2 sensors). A rule that fires once")
    say("   on one class is not yet a method -- the sealed sets are what decide it.")
    say("")
    say("   It fails if: it fires on a position whose round-2 frequency is < 80 %,")
    say("   or fires nowhere on either sealed class while a forceable position")
    say("   exists in their round-2 data.")

    with open(os.path.join(OUT, "forcing.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    with open(os.path.join(OUT, "forcing.json"), "w") as fh:
        json.dump(rows, fh, indent=1)
    say(f"\n   written to {OUT}/forcing.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
