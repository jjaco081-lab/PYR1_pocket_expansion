#!/usr/bin/env python
"""
33_merge_switch.py -- merge Arm 2b and turn dE_close into a two-sided verdict.

Reads results/switch/<variant>__r<N>.json, merges replicates, and reports

    ddE_switch = dE_close(variant) - dE_close(WT),   dE_close = E_closed - E_open

    ddE >> 0   closing destabilised   -> DEAD receptor (silent in the Y2H: it
                                         simply does not grow, indistinguishable
                                         from "no ligand in this pool binds")
    ddE ~  0   switch preserved       -> candidate
    ddE << 0   closed state favoured  -> CONSTITUTIVE (5-FOA removes it, but it
                                         wastes library capacity)

Runs the same discrimination check that condemned Arm 2 v1: if the between-variant
spread does not exceed replicate noise, the numbers mean nothing and nstruct must
go up before any of this is used.

Also verifies the conformers were actually held apart during relax. Each pose was
CA-constrained to its own starting state, so `gate_bb_rmsd_*` should stay small
and `gate_to_<other>_*` should stay near the starting stroke (~5.2 A). If a
relaxed open pose drifts toward closed, that replicate is not measuring what the
script claims and is flagged.

Usage:  python 33_merge_switch.py [--expect-reps 12] [--drift-max 2.0]
"""
import argparse, glob, json, os, collections
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "results", "switch")
OUT = os.path.join(ROOT, "results", "33_switch_summary.csv")

ap = argparse.ArgumentParser()
ap.add_argument("--expect-reps", type=int, default=12)
ap.add_argument("--drift-max", type=float, default=2.0,
                help="max allowed gate movement toward the other conformer (A)")
args = ap.parse_args()

files = sorted(glob.glob(os.path.join(D, "*__r*.json")))
if not files:
    raise SystemExit(f"no output in {D} -- has the array started?")

per = collections.defaultdict(list)
muts = {}
for f in files:
    j = json.load(open(f))
    per[j["variant"]].extend(j["reps"])
    muts[j["variant"]] = j["mutations"]
for v in per:
    seen, uniq = set(), []
    for r in sorted(per[v], key=lambda x: x["rep"]):
        if r["rep"] not in seen:
            seen.add(r["rep"])
            uniq.append(r)
    per[v] = uniq

print(f"{len(per)} variants, {sum(len(v) for v in per.values())} replicates")
short = [(v, len(r)) for v, r in sorted(per.items()) if len(r) < args.expect_reps]
if short:
    print(f"  incomplete ({args.expect_reps} expected): "
          + ", ".join(f"{v}:{n}" for v, n in short[:12])
          + (" ..." if len(short) > 12 else ""))

# ---- conformer integrity ----
print("\n" + "=" * 76)
print("CONFORMER INTEGRITY -- did the two states stay apart during relax?")
print("=" * 76)
bad = []
for v, reps in per.items():
    for r in reps:
        drift_o = r["gate_to_closed_open"]
        drift_c = r["gate_to_open_closed"]
        if drift_o < args.drift_max or drift_c < args.drift_max:
            bad.append((v, r["rep"], drift_o, drift_c))
allr = [r for reps in per.values() for r in reps]
for k in ("gate_bb_rmsd_closed", "gate_bb_rmsd_open",
          "gate_to_closed_open", "gate_to_open_closed"):
    x = np.array([r[k] for r in allr])
    print(f"  {k:<26} mean {x.mean():6.2f}  min {x.min():6.2f}  max {x.max():6.2f}")
if bad:
    print(f"  !! {len(bad)} replicates drifted toward the other conformer "
          f"(< {args.drift_max} A) -- these are NOT measuring two distinct states:")
    for v, rep, a, b in bad[:10]:
        print(f"       {v} r{rep}: open->closed {a:.2f}, closed->open {b:.2f}")
else:
    print(f"  all replicates kept the two conformers >= {args.drift_max} A apart")

if "WT" not in per:
    raise SystemExit("WT missing -- cannot compute ddE_switch")
wt_dE = float(np.mean([r["dE_close"] for r in per["WT"]]))
wt_sd = float(np.std([r["dE_close"] for r in per["WT"]], ddof=1)) if len(per["WT"]) > 1 else 0.0
print(f"\nWT dE_close = {wt_dE:.2f} +- {wt_sd:.2f} REU "
      f"(n={len(per['WT'])}); negative means WT intrinsically prefers CLOSED")

# ---- discrimination ----
print("\n" + "=" * 76)
print("DISCRIMINATION CHECK")
print("=" * 76)
print(f"{'metric':<22}{'within sd':>11}{'between sd':>12}{'ratio':>8}   verdict")
for k in ("dE_close", "score_closed", "score_open"):
    within = [np.std([r[k] for r in reps], ddof=1) for reps in per.values() if len(reps) > 1]
    means = [np.mean([r[k] for r in reps]) for reps in per.values()]
    if not within or len(means) < 2:
        continue
    w, b = float(np.mean(within)), float(np.std(means, ddof=1))
    ratio = b / w if w else float("nan")
    print(f"{k:<22}{w:>11.3f}{b:>12.3f}{ratio:>8.2f}   "
          + ("informative" if ratio >= 2.0 else "NOT DISCRIMINATING -- raise nstruct"))

# ---- table ----
rows = []
for v, reps in per.items():
    dE = np.array([r["dE_close"] for r in reps])
    m = float(dE.mean())
    sd = float(dE.std(ddof=1)) if len(dE) > 1 else 0.0
    sem = sd / np.sqrt(len(dE)) if len(dE) else float("nan")
    dd = m - wt_dE
    if v == "WT":
        verdict = "reference"
    elif abs(dd) <= 2 * max(sem, 0.5):
        verdict = "switch preserved (within noise of WT)"
    elif dd > 0:
        verdict = "closing DESTABILISED -- dead-receptor risk"
    else:
        verdict = "closed state FAVOURED -- constitutive risk"
    rows.append(dict(variant=v, mutations="+".join(f"{p}{a}" for p, a in muts[v]),
                     n=len(dE), dE_close=m, sd=sd, sem=sem, ddE_switch=dd,
                     score_closed=float(np.mean([r["score_closed"] for r in reps])),
                     score_open=float(np.mean([r["score_open"] for r in reps])),
                     verdict=verdict))

rows.sort(key=lambda r: r["ddE_switch"])
cols = list(rows[0])
with open(OUT, "w") as fh:
    fh.write(",".join(cols) + "\n")
    for r in rows:
        fh.write(",".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c])
                          for c in cols) + "\n")

print("\n" + "=" * 76)
print("SWITCH PREFERENCE, most constitutive-leaning first")
print("=" * 76)
print(f"{'variant':<26}{'n':>3}{'dE_close':>10}{'sd':>7}{'ddE':>8}   verdict")
for r in rows:
    print(f"{r['variant']:<26}{r['n']:>3}{r['dE_close']:>10.2f}{r['sd']:>7.2f}"
          f"{r['ddE_switch']:>+8.2f}   {r['verdict']}")
print(f"\nwrote {OUT}")
print("""
CAVEATS (see 31_rosetta_switch.py): fixed-backbone energy difference, not a free
energy -- no conformational entropy, no barrier. Both states apo. Read as rank
order, not kcal/mol. Calibrate the cutoff on Tian's 419 characterised functional
clones before using this to reject anything (README section 13d).""")
