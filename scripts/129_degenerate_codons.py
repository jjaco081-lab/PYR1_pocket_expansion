#!/usr/bin/env python
r"""
129_degenerate_codons.py -- how far do IUPAC mixed bases cut the oligo count, and
what do they cost in library size?

THE TRADE-OFF, which is the whole answer
Right now each variant within an assembly block is a SEPARATE oligo in the pool, so
the pool size is  sum over blocks of  prod over positions in that block of (1+n).
A degenerate codon collapses a whole position into ONE synthesis -- but the genetic
code does not let you choose an arbitrary amino-acid set. A codon that covers the
residues you want will usually drag in residues you did not ask for, which INFLATES
the library. So the question is not "how many oligos do we save" but "how many do
we save, at what library-size cost, and can we afford the extras".

WHAT IS COMPUTED
  * all 15^3 = 3,375 IUPAC codons, each expanded to the amino-acid set it encodes
  * for every position in every library, the best single codon covering
    {wild-type} U {menu}, ranked by how few EXTRA residues it adds
  * when no single codon is tolerable, the minimum NUMBER of codons whose union
    covers the target exactly (a set-cover) -- that many oligos are then needed at
    that position, and they multiply within a block
  * the resulting pool size and library size, against the current ones

⚠ Stop codons are disqualifying: a degenerate codon that can produce TAA/TAG/TGA
puts truncations into the pool, which is exactly the autoactivation/dead-clone
problem the DSM-Hao redesign existed to reduce. Codons containing a stop are
rejected outright rather than counted with a penalty.
"""
import itertools
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_xlsx import table                                     # noqa: E402

ROOT = os.path.dirname(HERE)
SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "oligo_design")

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
BASES = "TCAG"
AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
CODON = {}
i = 0
for b1 in BASES:
    for b2 in BASES:
        for b3 in BASES:
            CODON[b1 + b2 + b3] = AAS[i]
            i += 1


def encoded(deg):
    """Amino acids encoded by a 3-letter IUPAC codon, and whether it hits a stop."""
    out = set()
    stop = False
    for a in IUPAC[deg[0]]:
        for b in IUPAC[deg[1]]:
            for c in IUPAC[deg[2]]:
                aa = CODON[a + b + c]
                if aa == "*":
                    stop = True
                else:
                    out.add(aa)
    return frozenset(out), stop


ALL = {}
for d in itertools.product(IUPAC, repeat=3):
    deg = "".join(d)
    aas, stop = encoded(deg)
    ALL[deg] = (aas, stop)
CLEAN = {d: a for d, (a, s) in ALL.items() if not s}


def best_single(target):
    """Cheapest single stop-free codon covering `target`: fewest extra residues,
    then fewest codons in the mixture (less skew in the synthesised pool)."""
    cands = [(len(a - target), len(IUPAC[d[0]]) * len(IUPAC[d[1]]) * len(IUPAC[d[2]]), d, a)
             for d, a in CLEAN.items() if target <= a]
    if not cands:
        return None
    cands.sort()
    return cands[0]


def min_cover(target):
    """Fewest stop-free codons whose union is EXACTLY target (no extras allowed).
    Greedy over subsets that are themselves subsets of target -- exact cover of a
    <=20-element set, so a small search suffices."""
    sub = [(d, a) for d, a in CLEAN.items() if a <= target]
    if not sub:
        return None
    for k in range(1, 6):
        for combo in itertools.combinations(sub, k):
            u = set()
            for _d, a in combo:
                u |= a
            if u == target:
                return [d for d, _a in combo]
    return None


def main():
    os.makedirs(OUT, exist_ok=True)
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    _, recs = table(f"{SD}/pnas.2519924122.sd04.xlsx")
    libs = defaultdict(list)
    for r in recs:
        L = (r.get("Library") or "").strip()
        if not L:
            continue
        libs[L].append((int(float(r["Position"])), r["WT"].strip(),
                        r["Amino Acids allowed for mutation"].strip(),
                        str(r.get("Assembly block") or "1")))

    say("=" * 78)
    say("DEGENERATE CODONS: oligo count vs library size")
    say("=" * 78)
    say(f"   {len(CLEAN):,} of {len(ALL):,} IUPAC codons are stop-free and usable")
    say("")
    summary = {}
    for L in ("Coumarin", "PFAS", "TNTv1", "TNTv2"):
        if L not in libs:
            continue
        rows = sorted(libs[L])
        say("=" * 78)
        say(f"{L}  ({len(rows)} positions)")
        say("=" * 78)
        say(f"   {'pos':<6}{'want':<22}{'codon':>7}{'encodes':<24}"
            f"{'extra':<14}{'n_oligo'}")
        cur_block = defaultdict(lambda: 1)
        deg_block = defaultdict(lambda: 1)
        cur_size = 1
        deg_size = 1
        detail = []
        for pos, wt, menu, blk in rows:
            target = frozenset(set(menu) | {wt})
            bs = best_single(target)
            cov = min_cover(target)
            if bs is not None:
                extra_n, _mix, deg, aas = bs
                extras = "".join(sorted(aas - target))
                n_oligo = 1
                enc = aas
            else:
                deg = "+".join(cov) if cov else "?"
                enc = target
                extras = ""
                n_oligo = len(cov) if cov else len(target)
            detail.append({"pos": pos, "wt": wt, "want": "".join(sorted(target)),
                           "codon": deg, "encodes": "".join(sorted(enc)),
                           "extra": extras, "n_oligo": n_oligo, "block": blk})
            say(f"   {wt}{pos:<5}{''.join(sorted(target)):<22}{deg:>7} "
                f"{''.join(sorted(enc)):<24}{extras if extras else '-':<14}{n_oligo}")
            cur_block[blk] *= (1 + len(menu))
            deg_block[blk] *= n_oligo
            cur_size *= (1 + len(menu))
            deg_size *= len(enc)
        cur_oligos = sum(cur_block.values())
        deg_oligos = sum(deg_block.values())
        say("")
        say(f"   oligos to order : {cur_oligos:>12,}  ->  {deg_oligos:>10,}"
            f"   ({cur_oligos/deg_oligos:.0f}x fewer)")
        say(f"   library size    : {cur_size:>12,}  ->  {deg_size:>10,}"
            f"   ({deg_size/cur_size:.1f}x larger)")
        say(f"   per block, oligos: current {dict(cur_block)} -> degenerate "
            f"{dict(deg_block)}")
        summary[L] = {"cur_oligos": cur_oligos, "deg_oligos": deg_oligos,
                      "cur_size": cur_size, "deg_size": deg_size,
                      "positions": detail}
    # ---- three strategies, as a trade-off rather than one answer -----------
    say("")
    say("=" * 78)
    say("THE TRADE-OFF: three ways to use mixed bases")
    say("=" * 78)
    say("   EXACT     multiple codons per position, union == the intended menu.")
    say("             Library is UNCHANGED; oligos = product of cover sizes per block.")
    say("   TOLERANT  one codon per position, but only where it adds <= 2 extras;")
    say("             exact cover elsewhere.")
    say("   LOOSE     one codon per position always (what the table above used).")
    say("")
    say(f"   {'library':<10}{'strategy':<10}{'oligos':>10}{'vs now':>9}"
        f"{'library':>16}{'vs now':>10}")
    strat = {}
    for L in ("Coumarin", "PFAS", "TNTv1", "TNTv2"):
        if L not in libs:
            continue
        rows = sorted(libs[L])
        cur_o = defaultdict(lambda: 1)
        cur_s = 1
        for pos, wt, menu, blk in rows:
            cur_o[blk] *= (1 + len(menu))
            cur_s *= (1 + len(menu))
        cur_oligos = sum(cur_o.values())
        say(f"   {L:<10}{'current':<10}{cur_oligos:>10,}{'1x':>9}{cur_s:>16,}{'1x':>10}")
        for name, maxextra in (("EXACT", 0), ("TOLERANT", 2), ("LOOSE", 99)):
            ob = defaultdict(lambda: 1)
            sz = 1
            for pos, wt, menu, blk in rows:
                target = frozenset(set(menu) | {wt})
                bs = best_single(target)
                use_single = bs is not None and bs[0] <= maxextra
                if use_single:
                    ob[blk] *= 1
                    sz *= len(bs[3])
                else:
                    cov = min_cover(target)
                    k = len(cov) if cov else len(target)
                    ob[blk] *= k
                    sz *= len(target)
            o = sum(ob.values())
            strat.setdefault(L, {})[name] = (o, sz)
            say(f"   {'':<10}{name:<10}{o:>10,}{cur_oligos/o:>8.0f}x{sz:>16,}"
                f"{sz/cur_s:>9.1f}x")
        say("")

    say("")
    say("=" * 78)
    say("SUMMARY")
    say("=" * 78)
    say(f"   {'library':<10}{'oligos now':>12}{'oligos deg':>12}{'saving':>9}"
        f"{'library now':>14}{'library deg':>14}{'inflation':>11}")
    for L, v in summary.items():
        say(f"   {L:<10}{v['cur_oligos']:>12,}{v['deg_oligos']:>12,}"
            f"{v['cur_oligos']/v['deg_oligos']:>8.0f}x{v['cur_size']:>14,}"
            f"{v['deg_size']:>14,}{v['deg_size']/v['cur_size']:>10.1f}x")
    with open(os.path.join(OUT, "degenerate_codons.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    with open(os.path.join(OUT, "degenerate_codons.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    say(f"\n   written to {OUT}/degenerate_codons.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
