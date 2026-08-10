#!/usr/bin/env python
"""
23_merge_ratchet_v2.py -- concatenate the chunked Arm-2 v2 output and re-run the
discrimination check that condemned v1.

Each array task of 22_submit_ratchet_v2.sh writes
results/ratchet_v2/<variant>__c<chunk>.json holding 4 replicates. This merges the
5 chunks per variant into <variant>.json (20 replicates), then reports, for every
metric, the between-variant spread against the within-variant replicate noise --
the same test that showed v1 was blind (README section 5b).

If the ratio is still below ~2 after this, the problem is not sampling and not
the constraints, and the metric itself has to change.

Usage:  python 23_merge_ratchet_v2.py [--expect-reps 20]
"""
import argparse, glob, json, os, collections
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ap = argparse.ArgumentParser()
ap.add_argument("--expect-reps", type=int, default=20)
ap.add_argument("--dir", default="ratchet_v2",
                help="results subdirectory: ratchet_v2 or ratchet_v3")
args = ap.parse_args()
D = os.path.join(ROOT, "results", args.dir)

chunks = collections.defaultdict(list)
muts = {}
for f in sorted(glob.glob(os.path.join(D, "*__c*.json"))):
    j = json.load(open(f))
    chunks[j["variant"]].append(j)
    muts[j["variant"]] = j["mutations"]

if not chunks:
    raise SystemExit(f"no chunk files in {D} -- has the array finished?")

print(f"{len(chunks)} variants with chunk files")
merged = {}
short = []
for v, js in sorted(chunks.items()):
    reps = sorted((r for j in js for r in j["reps"]), key=lambda r: r["rep"])
    seen, uniq = set(), []
    for r in reps:                      # guard against a re-run duplicating a chunk
        if r["rep"] not in seen:
            seen.add(r["rep"])
            uniq.append(r)
    merged[v] = uniq
    if len(uniq) < args.expect_reps:
        short.append((v, len(uniq)))
    json.dump(dict(variant=v, mutations=muts[v], reps=uniq),
              open(os.path.join(D, f"{v}.json"), "w"), indent=2)

if short:
    print(f"\n!! {len(short)} variants have fewer than {args.expect_reps} replicates:")
    for v, n in short:
        print(f"     {v:<26} {n}")
    print("   (missing chunks -- check logs/ratchet2_*.log before trusting the stats)")
else:
    print(f"all variants have {args.expect_reps} replicates")

KEYS = ["dG_separated", "dSASA", "sc_value",
        "gate_bb_rmsd", "latch_bb_rmsd", "gate_latch_min_dist",
        "gate_latch_com_dist", "lig_gate_min_dist", "w385_aba_min_dist",
        "gate_rmsd", "latch_rmsd",
        # v3 only -- the HAB1-contacting 149-156 loop (README section 9d)
        "iloop_bb_rmsd", "iloop_hab1_min_dist"]

print("\n" + "=" * 82)
print(f"DISCRIMINATION CHECK -- {args.dir} (nstruct={args.expect_reps})")
print("=" * 82)
print(f"{'metric':<24}{'within sd':>11}{'between sd':>12}{'ratio':>8}   verdict")
for k in KEYS:
    within, means = [], []
    for v, reps in merged.items():
        x = [r[k] for r in reps if k in r]
        if len(x) > 1:
            within.append(np.std(x, ddof=1))
        if x:
            means.append(np.mean(x))
    if not within or len(means) < 2:
        continue
    w, b = float(np.mean(within)), float(np.std(means, ddof=1))
    r = b / w if w else float("nan")
    print(f"{k:<24}{w:>11.3f}{b:>12.3f}{r:>8.2f}   "
          + ("informative" if r >= 2.0 else "still not discriminating"))

print("""
ratio = (sd of per-variant means) / (mean within-variant replicate sd).
Compare against v1 (README section 5b): dG_separated 0.73, gate_rmsd 0.65,
latch_rmsd 0.92 -- all blind. gate_bb_rmsd and gate_latch_min_dist are the
metrics the v1 constraints made impossible to measure; they are the ones to
watch here.

Next: point 20_aggregate_pilot.py at results/ratchet_v2 to rebuild the joined
decision table with a readout arm that actually reads out.""")
