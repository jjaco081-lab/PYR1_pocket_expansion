#!/usr/bin/env python
"""
11_tian_library_overlap.py -- what has Tian et al. already sampled at our positions?

Source: supplementary datasets of Tian et al., PNAS (doi 10.1073/pnas.2519924122),
downloaded to /bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark/:
  sd04  Dataset S3  library DESIGN: per-position list of allowed amino acids
  sd07  Dataset S6  characterised hits (imperatorin etc.), per-position residues
  sd08  Dataset S7  characterised hits (TNT/DNT etc.)
  sd09  Dataset S8  characterised hits (PFOA/PFOS etc.)

Purpose
-------
Their libraries randomise pocket positions to REPROGRAM specificity within the
existing cavity envelope. This project instead tries to ENLARGE the cavity.
The question is which of our six wall positions (59/79/94/108/120/141) they
already varied, to which residues, and -- critically -- whether they ever
sampled the small/charge-removing substitutions that open volume, or only
side-chain swaps of comparable bulk.

Run with the iggypop env python (has openpyxl/pandas).
"""
import os, sys, collections
import openpyxl

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OURS = [59, 79, 94, 108, 120, 141]
WT = dict(zip(OURS, "KREFYE"))
# van der Waals volumes (A^3, Richards 1974) for a "does this open space?" call
VOL = {"G": 60, "A": 89, "S": 89, "C": 109, "D": 111, "P": 113, "N": 114,
       "T": 116, "E": 138, "V": 140, "Q": 144, "H": 153, "M": 163, "I": 167,
       "L": 167, "K": 169, "R": 174, "F": 190, "Y": 194, "W": 228}


def sheet_rows(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()
    return rows


# ---------- 1. library DESIGN (sd04) ----------
print("=" * 78)
print("LIBRARY DESIGN -- Tian et al. Dataset S3 (sd04): allowed substitutions")
print("=" * 78)
rows = sheet_rows(os.path.join(SD, "pnas.2519924122.sd04.xlsx"))
hdr = [str(c).strip() if c else "" for c in rows[0]]
i_lib, i_pos, i_wt, i_aa = 0, 1, 2, 3

design = collections.defaultdict(dict)   # position -> {library: allowed}
all_pos = set()
for r in rows[1:]:
    if r[i_pos] is None:
        continue
    pos = int(float(r[i_pos]))
    lib = str(r[i_lib]).strip()
    allowed = str(r[i_aa]).strip() if r[i_aa] else ""
    design[pos][lib] = (str(r[i_wt]).strip(), allowed)
    all_pos.add(pos)

libs = sorted({l for d in design.values() for l in d})
print(f"\nlibraries: {', '.join(libs)}")
print(f"positions randomised across all libraries ({len(all_pos)}): "
      f"{sorted(all_pos)}\n")

print("--- OUR SIX WALL POSITIONS ---")
for p in OURS:
    print(f"\n  PYR1 {WT[p]}{p}:")
    if p not in design:
        print("     NOT randomised in any Tian library  <-- unexplored")
        continue
    for lib, (wt, allowed) in sorted(design[p].items()):
        aas = [a for a in allowed if a.isalpha()]
        smaller = [a for a in aas if VOL[a] < VOL[WT[p]] - 20]
        print(f"     {lib:<6} WT={wt}  allowed({len(aas)}): {allowed}")
        if smaller:
            print(f"            volume-opening subs present: {''.join(smaller)} "
                  f"(WT {WT[p]}={VOL[WT[p]]} A^3)")
        else:
            print(f"            NO volume-opening substitution "
                  f"(nothing < {VOL[WT[p]]-20} A^3)")

print("\n--- positions Tian randomised that we did NOT pick ---")
print("   ", sorted(set(all_pos) - set(OURS)))
print("--- our positions Tian never randomised ---")
print("   ", sorted(set(OURS) - set(all_pos)))


# ---------- 2. observed substitutions in characterised hits ----------
print("\n" + "=" * 78)
print("OBSERVED SUBSTITUTIONS in characterised hits (sd07 / sd08 / sd09)")
print("=" * 78)
obs = collections.defaultdict(collections.Counter)
files = ["pnas.2519924122.sd07(1).xlsx", "pnas.2519924122.sd08.xlsx",
         "pnas.2519924122.sd09.xlsx"]
n_clones = 0
for fn in files:
    rows = sheet_rows(os.path.join(SD, fn))
    hdr = [str(c).strip() if c else "" for c in rows[0]]
    colmap = {}
    for j, h in enumerate(hdr):
        if h and h[0].isalpha() and h[1:].isdigit():
            colmap[int(h[1:])] = j            # e.g. "K59" -> column index
    for r in rows[1:]:
        if all(c is None for c in r):
            continue
        n_clones += 1
        for pos, j in colmap.items():
            v = r[j]
            if v and str(v).strip():
                obs[pos][str(v).strip()[0]] += 1
    print(f"  {fn}: {len(colmap)} position columns")

print(f"\n{n_clones} characterised clones total\n")
print("--- what was actually RECOVERED at our six positions ---")
for p in OURS:
    c = obs.get(p)
    if not c:
        print(f"  {WT[p]}{p}: never mutated in any characterised hit")
        continue
    tot = sum(c.values())
    items = ", ".join(f"{a}:{n}" for a, n in c.most_common())
    smaller = [a for a in c if VOL.get(a, 999) < VOL[WT[p]] - 20]
    print(f"  {WT[p]}{p}: {tot} clones -- {items}")
    print(f"        volume-opening recovered: "
          f"{''.join(sorted(smaller)) if smaller else 'NONE'}")

print("\n--- all positions seen mutated in characterised hits ---")
print("   ", sorted(obs))
