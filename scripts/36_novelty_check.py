#!/usr/bin/env python
"""
36_novelty_check.py -- is a given variant actually NEW relative to Tian's work?

Refined win condition (user, 2026-08-08)
----------------------------------------
Novelty is NOT simply "a position the libraries never touched". A change at a
position Tian did randomise is still new if

  * the SUBSTITUTION was not among the amino acids their oligo pools allowed at
    that position, or
  * the substitution was allowed but never RECOVERED in a characterised clone, or
  * the COMBINATION of changes never co-occurred in any characterised clone

...because their libraries were screened against a chemical library, so what has
really been tested is a set of (sequence, chemistry) pairs. What must be avoided
is re-testing combinations that were already tested.

This classifies any variant on three levels, from Tian's supplementary data:

  L1 POSITION      was the position randomised at all? (sd04)
  L2 SUBSTITUTION  was this exact residue allowed by the oligo pools? (sd04)
                   and was it actually recovered in a characterised hit?
                   (sd07 / sd08 / sd09)
  L3 COMBINATION   did this exact set of substitutions co-occur in any single
                   characterised clone? and did the SET OF POSITIONS co-occur,
                   regardless of identity?

Verdict per variant:
  NOVEL-POSITION    contains a position never randomised          (strongest)
  NOVEL-SUBST       every position sampled, but >=1 substitution outside the
                    allowed pools -- unreachable by their libraries
  NOVEL-UNRECOVERED all substitutions allowed, but >=1 never seen in a hit
  NOVEL-COMBO       all substitutions individually recovered, but never together
  RETREAD           this combination was already made and screened -- avoid

Run with the dockenv python (needs openpyxl).
"""
import os, sys, collections, json, csv
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from variants import PANEL, WT   # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "36_novelty.csv")
HIT_FILES = ["pnas.2519924122.sd07(1).xlsx", "pnas.2519924122.sd08.xlsx",
             "pnas.2519924122.sd09.xlsx"]


def sheet_rows(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()
    return rows


# ---------- L1/L2: library design (sd04) ----------
allowed = collections.defaultdict(set)     # pos -> set of allowed AAs (any library)
libs_at = collections.defaultdict(set)     # pos -> libraries randomising it
for r in sheet_rows(os.path.join(SD, "pnas.2519924122.sd04.xlsx"))[1:]:
    if r[1] is None:
        continue
    pos = int(float(r[1]))
    allowed[pos] |= set(str(r[3]).strip()) if r[3] else set()
    libs_at[pos].add(str(r[0]).strip())
print(f"sd04: {len(allowed)} randomised positions across "
      f"{len(set().union(*libs_at.values()))} libraries")

# ---------- L2/L3: characterised clones (sd07/08/09) ----------
clones = []          # list of dict pos->aa
recovered = collections.defaultdict(set)   # pos -> AAs seen in any hit
for fn in HIT_FILES:
    path = os.path.join(SD, fn)
    if not os.path.exists(path):
        print(f"  !! missing {fn}")
        continue
    rows = sheet_rows(path)
    hdr = [str(c).strip() if c is not None else "" for c in rows[0]]
    poscols = {}
    for i, h in enumerate(hdr):
        h2 = h.strip()
        if len(h2) >= 2 and h2[0].isalpha() and h2[1:].isdigit():
            poscols[i] = int(h2[1:])
    for r in rows[1:]:
        combo = {}
        for i, pos in poscols.items():
            v = r[i]
            if v is None:
                continue
            aa = str(v).strip().upper()
            if len(aa) == 1 and aa.isalpha():
                combo[pos] = aa
                recovered[pos].add(aa)
        if combo:
            clones.append(combo)
print(f"characterised clones parsed: {len(clones)}")

combo_keys = {frozenset(c.items()) for c in clones}
posset_keys = {frozenset(c) for c in clones}


def classify(muts):
    """muts: list of (pos, aa)."""
    if not muts:
        return "WT", "", ""
    notes = []
    unsampled = [p for p, _ in muts if p not in allowed]
    if unsampled:
        return ("NOVEL-POSITION",
                "position(s) never randomised: " + ",".join(map(str, unsampled)), "")
    not_allowed = [f"{WT.get(p,'?')}{p}{a}" for p, a in muts if a not in allowed[p]]
    if not_allowed:
        return ("NOVEL-SUBST",
                "substitution(s) outside the oligo pools: " + ",".join(not_allowed), "")
    not_recovered = [f"{WT.get(p,'?')}{p}{a}" for p, a in muts
                     if a not in recovered.get(p, set())]
    if not_recovered:
        return ("NOVEL-UNRECOVERED",
                "allowed but never recovered in a hit: " + ",".join(not_recovered), "")
    key = frozenset((p, a) for p, a in muts)
    if key in combo_keys:
        return ("RETREAD", "this exact combination appears in a characterised clone", "")
    pk = frozenset(p for p, _ in muts)
    same_pos = "position set also co-occurred (different residues)" if pk in posset_keys \
        else "position set never co-occurred either"
    return ("NOVEL-COMBO", "all substitutions individually recovered, never together",
            same_pos)


print("\n" + "=" * 100)
print("NOVELTY OF THE 41-VARIANT PANEL vs Tian's libraries and characterised clones")
print("=" * 100)
rows_out = []
order = {"NOVEL-POSITION": 0, "NOVEL-SUBST": 1, "NOVEL-UNRECOVERED": 2,
         "NOVEL-COMBO": 3, "RETREAD": 4, "WT": 5}
for name, muts in PANEL:
    verdict, why, extra = classify(muts)
    rows_out.append(dict(variant=name,
                         mutations="+".join(f"{WT.get(p,'?')}{p}{a}" for p, a in muts),
                         verdict=verdict, reason=why, extra=extra))
rows_out.sort(key=lambda r: (order.get(r["verdict"], 9), r["variant"]))

cur = None
for r in rows_out:
    if r["verdict"] != cur:
        cur = r["verdict"]
        print(f"\n--- {cur} ---")
    print(f"  {r['variant']:<26}{r['mutations']:<28}{r['reason']}"
          + (f" [{r['extra']}]" if r["extra"] else ""))

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["variant", "mutations", "verdict",
                                       "reason", "extra"])
    w.writeheader()
    w.writerows(rows_out)
print(f"\nwrote {OUT}")

# ---------- what the libraries allow at our six positions ----------
print("\n" + "=" * 100)
print("PER-POSITION DETAIL at the six wall positions")
print("=" * 100)
print(f"{'pos':<6}{'WT':<4}{'randomised?':<13}{'allowed by pools':<26}"
      f"{'recovered in hits':<26}{'allowed but NEVER recovered'}")
for p in sorted(WT):
    a = "".join(sorted(allowed.get(p, set())))
    rec = "".join(sorted(recovered.get(p, set())))
    gap = "".join(sorted(set(a) - set(rec)))
    print(f"{p:<6}{WT[p]:<4}{('yes: '+','.join(sorted(libs_at[p]))) if p in allowed else 'NO':<13}"
          f"{a or '-':<26}{rec or '-':<26}{gap or '-'}")

print("""
READING THIS
  "allowed but NEVER recovered" is the cheapest source of genuine novelty at an
  already-sampled position: the oligo pools could make it, the screen was run,
  and it did not come back -- so that chemistry is untested in practice rather
  than merely unexplored. Combining one of those with a NOVEL-POSITION change
  (R79) gives a variant that is new at two independent levels.""")
