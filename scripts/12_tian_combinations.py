#!/usr/bin/env python
"""
12_tian_combinations.py -- did Tian et al. ever break the WALL, or only swap
side chains inside the existing cavity envelope?

The single-position analysis (script 11) shows they randomised 5 of our 6 wall
positions. But pocket EXPANSION requires simultaneously removing bulk/charge
from the wall cluster, not swapping one residue at a time. This script asks:

  1. For each characterised clone, how many of the wall positions
     (59 / 79 / 94 / 108 / 120 / 141) are mutated, and what is the net change
     in side-chain volume summed over those positions?
  2. Which clones show the largest net volume RELEASE at the wall?
  3. Was the E94-R79 salt bridge ever broken? (R79 is never mutated, so the
     bridge can only be broken from the E94 side.)
  4. Are the K59 + F108 positions -- the two top restrictors -- ever
     simultaneously reduced?

Run with the iggypop env python.
"""
import os, collections
import openpyxl

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
WALL = [59, 79, 94, 108, 120, 141]
WT = dict(zip(WALL, "KREFYE"))
VOL = {"G": 60, "A": 89, "S": 89, "C": 109, "D": 111, "P": 113, "N": 114,
       "T": 116, "E": 138, "V": 140, "Q": 144, "H": 153, "M": 163, "I": 167,
       "L": 167, "K": 169, "R": 174, "F": 190, "Y": 194, "W": 228}
FILES = ["pnas.2519924122.sd07(1).xlsx", "pnas.2519924122.sd08.xlsx",
         "pnas.2519924122.sd09.xlsx"]


def load(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()
    return rows


clones = []
for fn in FILES:
    rows = load(os.path.join(SD, fn))
    hdr = [str(c).strip() if c else "" for c in rows[0]]
    colmap = {int(h[1:]): j for j, h in enumerate(hdr)
              if h and h[0].isalpha() and h[1:].isdigit()}
    # locate descriptive columns
    j_mut = next((j for j, h in enumerate(hdr) if h.lower() == "mutation"), None)
    # sd07/sd08 call it "compound"; sd09 calls it "chem_name"
    j_cmp = next((j for j, h in enumerate(hdr)
                  if h.lower() in ("compound", "chem_name")), None)
    j_nm = next((j for j, h in enumerate(hdr) if h.lower() in ("name", "clone")), None)
    for r in rows[1:]:
        if all(c is None for c in r):
            continue
        muts = {}
        for pos, j in colmap.items():
            v = r[j]
            if v and str(v).strip():
                muts[pos] = str(v).strip()[0]
        clones.append(dict(
            src=fn.split(".")[2],
            name=str(r[j_nm]).strip() if j_nm is not None and r[j_nm] else "",
            compound=str(r[j_cmp]).strip() if j_cmp is not None and r[j_cmp] else "",
            mutation=str(r[j_mut]).strip() if j_mut is not None and r[j_mut] else "",
            muts=muts))

print(f"{len(clones)} characterised clones loaded\n")

for c in clones:
    wall_muts = {p: a for p, a in c["muts"].items() if p in WALL}
    c["wall_muts"] = wall_muts
    c["n_wall"] = len(wall_muts)
    c["dvol"] = sum(VOL.get(a, VOL[WT[p]]) - VOL[WT[p]] for p, a in wall_muts.items())

print("=== how many wall positions are mutated per clone ===")
h = collections.Counter(c["n_wall"] for c in clones)
for k in sorted(h):
    print(f"  {k} wall position(s) mutated: {h[k]} clones")

print("\n=== net side-chain volume change at the wall (negative = opens space) ===")
srt = sorted(clones, key=lambda c: c["dvol"])
print(f"  most volume RELEASED: {srt[0]['dvol']:+.0f} A^3")
print(f"  most volume ADDED   : {srt[-1]['dvol']:+.0f} A^3")
neg = [c for c in clones if c["dvol"] < 0]
print(f"  clones with net release: {len(neg)}/{len(clones)}")

print("\n--- top 12 clones by wall volume released ---")
print(f"{'clone':<14}{'compound':<26}{'dVol':>7}  wall substitutions")
for c in srt[:12]:
    w = " ".join(f"{WT[p]}{p}{a}" for p, a in sorted(c["wall_muts"].items()))
    print(f"{c['name'][:13]:<14}{c['compound'][:25]:<26}{c['dvol']:>+7.0f}  {w}")

print("\n=== was the E94-R79 salt bridge ever broken? ===")
print("  R79 mutated in any clone: "
      f"{sum(1 for c in clones if 79 in c['muts'])} clones")
e94 = [c for c in clones if 94 in c["muts"]]
e94_break = [c for c in e94 if c["muts"][94] not in "DEQN"]
print(f"  E94 mutated: {len(e94)} clones; of those, to a residue that CANNOT "
      f"salt-bridge R79 (not D/E/Q/N): {len(e94_break)}")
c2 = collections.Counter(c["muts"][94] for c in e94)
print(f"  E94 substitutions recovered: "
      + ", ".join(f"{a}:{n}" for a, n in c2.most_common()))

print("\n=== K59 and F108 (the two top restrictors) simultaneously reduced? ===")
both = [c for c in clones if 59 in c["muts"] and 108 in c["muts"]]
print(f"  clones with BOTH mutated: {len(both)}")
red = [c for c in both
       if VOL.get(c["muts"][59], 999) < VOL["K"] - 20
       and VOL.get(c["muts"][108], 999) < VOL["F"] - 20]
print(f"  ...with BOTH reduced in volume: {len(red)}")
for c in both[:10]:
    print(f"     {c['name'][:16]:<17}{c['compound'][:22]:<24}"
          f"K59{c['muts'][59]}  F108{c['muts'][108]}")

print("\n=== F108: which way did selection go? ===")
f = collections.Counter(c["muts"][108] for c in clones if 108 in c["muts"])
tot = sum(f.values())
bigger = sum(n for a, n in f.items() if VOL.get(a, 0) > VOL["F"])
print(f"  {tot} clones mutate F108; {bigger} ({100*bigger/max(tot,1):.0f}%) go "
      f"LARGER than Phe (190 A^3)")
print("  distribution: " + ", ".join(
    f"{a}({VOL.get(a,'?')}):{n}" for a, n in f.most_common()))
