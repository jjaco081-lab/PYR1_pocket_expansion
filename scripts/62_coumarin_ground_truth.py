#!/usr/bin/env python
"""
62_coumarin_ground_truth.py -- build the ground truth for the "could we have
skipped Tian's round-1 screen?" benchmark, and measure whether that benchmark has
any headroom before anyone builds a method for it.

THE QUESTION (user, 2026-08-14)
-------------------------------
Tian et al. mined their primary screen for the coumarin cluster, derived a
sequence profile from the hits, and built a focused 77,327-member library that
produced sensors for three ligands the original screen had missed. Could the
positions -- and ideally the substitutions -- have been predicted computationally,
so that round-1 screen would not have to be done? It does not have to be perfect;
it has to narrow the window.

WHERE THE TRUTH LIVES
---------------------
`pnas.2519924122.sd04.xlsx` (internally "Dataset S3") holds the DESIGN of every
library: per library, per position, the exact amino acids offered. This includes
the focused Coumarin, PFAS, TNTv1 and TNTv2 libraries, not just DSM-Hao/TSM.
That is the ground truth: it is what Tian actually chose to build after seeing
round 1.

  Coumarin library = 11 positions, 24 substitutions
      [59, 81, 83, 108, 120, 122, 159, 160, 163, 164, 167]

Those are the "11 of the 19 binding-pocket sites" the paper reports.

`pnas.2519924122.sd03(1).xlsx` ("Dataset S3" as printed, 692 characterised clones)
holds the round-1 SCREEN RESULTS -- the input a predictor would be allowed to use
only in the "with screen data" variant of the task.

THE HEADROOM RESULT, WHICH DECIDES THE TASK DEFINITION
------------------------------------------------------
Positions are nearly ligand-INDEPENDENT and are therefore almost useless as a
prediction target:

    Coumarin/PFAS/TNTv2 share 9 of their 11-14 positions
    pairwise Jaccard 0.60-0.80; exactly ONE position is unique to each library
    (F108 coumarin, E94 PFAS, L117 TNTv2)

So "name the nine positions every focused library uses" already scores 9/11 recall
on coumarin with zero ligand information. A position-prediction benchmark has no
headroom, which is the same defect that sank the stage-1 benchmark (README 23j),
and it is why this script measures headroom BEFORE the task is fixed.

The substitutions are the opposite -- they are almost perfectly ligand-specific:

    mean 3-way Jaccard of the substitution menus at the 9 shared positions = 0.02
    at 8 of those 9 positions the three libraries share NO allowed residue at all
    e.g. V163 -> W (coumarin) / S (PFAS) / G,M,W (TNTv2)

CONCLUSION: the pocket positions are a fixed lining set that Tian re-randomises
every time. The ligand information lives entirely in WHICH residue goes there.
The benchmark must therefore predict SUBSTITUTION MENUS, not positions.

Run with the iggypop env python (has openpyxl).
"""
import collections
import itertools
import json
import os

import openpyxl

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "coumarin_benchmark")
os.makedirs(OUT, exist_ok=True)

DESIGN = os.path.join(SD, "pnas.2519924122.sd04.xlsx")
SCREEN = os.path.join(SD, "pnas.2519924122.sd03(1).xlsx")
AA_COL = "Amino Acids allowed for mutation"
FOCUSED = ["Coumarin", "PFAS", "TNTv2"]


def sheet(path):
    wb = openpyxl.load_workbook(path, read_only=True)
    rows = list(wb.worksheets[0].iter_rows(values_only=True))
    wb.close()
    hdr = [str(c) if c is not None else "" for c in rows[0]]
    return hdr, [dict(zip(hdr, r)) for r in rows[1:] if r and any(x is not None for x in r)]


# ---------------------------------------------------------------- library designs
hdr, recs = sheet(DESIGN)
assert AA_COL in hdr, (
    f"column {AA_COL!r} missing from the design sheet; found {hdr}. Do NOT guess "
    f"the column name -- an earlier read truncated it to 'Amino Acids allowed ' "
    f"purely because the display was clipped at 20 characters.")

libs, wt = collections.defaultdict(dict), {}
for r in recs:
    try:
        pos = int(float(r["Position"]))
    except (TypeError, ValueError):
        continue
    lib = str(r["Library"]).strip()
    aas = set(str(r[AA_COL]).strip())
    libs[lib][pos] = aas
    wt[pos] = str(r["WT"]).strip()

print("=" * 76)
print("LIBRARY DESIGNS (ground truth)")
print("=" * 76)
for lib in sorted(libs):
    n_sub = sum(len(v) for v in libs[lib].values())
    print(f"  {lib:<10} {len(libs[lib]):>2} positions, {n_sub:>3} substitutions")

cou = libs["Coumarin"]
print(f"\n  Coumarin positions: {sorted(cou)}")
print(f"  -> this is the paper's \"11 of the 19 binding-pocket sites\"")
for p in sorted(cou):
    print(f"       {wt[p]}{p:<5} {''.join(sorted(cou[p]))}")

# ---------------------------------------------------------------- headroom
print()
print("=" * 76)
print("HEADROOM: are POSITIONS ligand-specific?  (if not, do not predict them)")
print("=" * 76)
for a, b in itertools.combinations(FOCUSED, 2):
    A, B = set(libs[a]), set(libs[b])
    print(f"  {a:<9} vs {b:<7} shared {len(A & B):>2}/{len(A | B):<2}  "
          f"Jaccard {len(A & B) / len(A | B):.2f}")
common = set.intersection(*[set(libs[L]) for L in FOCUSED])
print(f"  shared by all three: {len(common)} -> {sorted(common)}")
for L in FOCUSED:
    uniq = set(libs[L]) - set().union(*[set(libs[x]) for x in FOCUSED if x != L])
    print(f"  unique to {L:<9}: {sorted(uniq) or 'none'}")
hit = len(common & set(cou))
print(f"\n  BASELINE: naming just the shared-by-all-three set scores "
      f"{hit}/{len(cou)} = {hit/len(cou):.0%} recall on Coumarin,")
print(f"  using NO ligand information at all. Position prediction has no headroom.")

print()
print("=" * 76)
print("HEADROOM: are SUBSTITUTIONS ligand-specific?  (this is the real target)")
print("=" * 76)
print(f"  {'pos':<5}{'WT':<4}" + "".join(f"{L:<13}" for L in FOCUSED) + "3-way J")
jac = []
for p in sorted(common):
    sets = [libs[L][p] for L in FOCUSED]
    j = len(set.intersection(*sets)) / len(set.union(*sets))
    jac.append(j)
    print(f"  {p:<5}{wt[p]:<4}" + "".join(f"{''.join(sorted(s)):<13}" for s in sets)
          + f"{j:.2f}")
mean_j = sum(jac) / len(jac)
print(f"\n  mean 3-way Jaccard of substitution menus = {mean_j:.2f}")
print(f"  ({sum(1 for j in jac if j == 0)} of {len(jac)} shared positions share NO "
      f"allowed residue between the three libraries)")
print("  -> the ligand signal is entirely here. THIS is what a method must predict.")

# ---------------------------------------------------------------- library sizes
print()
print("=" * 76)
print("WHAT 'NARROWING THE WINDOW' IS WORTH, in library members")
print("=" * 76)


def size(menu):
    """combinatorial size if every position varies independently (WT + allowed)"""
    n = 1
    for aas in menu.values():
        n *= (len(aas) + 1)          # +1 for keeping WT
    return n


for lib in ("DSM-Hao", "TSM", "Coumarin", "PFAS", "TNTv2"):
    if lib in libs:
        print(f"  {lib:<10} {len(libs[lib]):>2} pos  ->  {size(libs[lib]):>18,d} "
              f"combinatorial members")
print("\n  Tian's Coumarin library as BUILT was 77,327 members (paper), against a")
print("  DSM-Hao design space many orders larger. That ratio is the prize: a method")
print("  that proposes the right menu turns an unscreenable space into a screenable")
print("  library WITHOUT running round 1.")

# ---------------------------------------------------------------- save
truth = {
    "note": "Tian et al. PNAS 2519924122. Library designs from sd04 (per library, "
            "per position, allowed amino acids). Positions are near-ligand-"
            "independent (3-way Jaccard 0.60-0.80); substitution menus are not "
            "(mean 3-way Jaccard 0.02). Predict MENUS, not positions.",
    "wt": wt,
    "libraries": {k: {str(p): sorted(v) for p, v in d.items()} for k, d in libs.items()},
    "shared_positions_focused": sorted(common),
    "position_jaccard": {f"{a}|{b}":
                         len(set(libs[a]) & set(libs[b])) / len(set(libs[a]) | set(libs[b]))
                         for a, b in itertools.combinations(FOCUSED, 2)},
    "mean_substitution_jaccard_3way": mean_j,
}
path = os.path.join(OUT, "tian_library_truth.json")
json.dump(truth, open(path, "w"), indent=1)
print(f"\nwrote {path}")

# ---------------------------------------------------------------- ligand export
# The screen table is the only place the ligand SMILES live, and it needs openpyxl,
# which no RDKit env here has. So export a small CSV once; downstream scripts
# (docking, MD setup) then need RDKit only. Tracked, because it is a ground-truth
# extract rather than bulk.
hdr_s, recs_s = sheet(SCREEN)
seen, out_rows = set(), []
for r in recs_s:
    n = str(r.get("library_name") or "").strip()
    if not n or n in seen:
        continue
    seen.add(n)
    out_rows.append(dict(
        library_name=n,
        parent_name=str(r.get("parent_name") or ""),
        chem_cluster=str(r.get("chem_cluster") or ""),
        cat=str(r.get("cat") or ""),
        canonical_smiles=str(r.get("canonical_smiles") or ""),
    ))
counts = collections.Counter(str(r.get("library_name") or "").strip() for r in recs_s)
for row in out_rows:
    row["n_clones"] = counts[row["library_name"]]

import csv as _csv
lp = os.path.join(OUT, "tian_screen_ligands.csv")
with open(lp, "w", newline="") as fh:
    w = _csv.DictWriter(fh, fieldnames=list(out_rows[0]))
    w.writeheader()
    w.writerows(sorted(out_rows, key=lambda d: -d["n_clones"]))
print(f"wrote {lp}  ({len(out_rows)} distinct ligands with characterised clones)")
