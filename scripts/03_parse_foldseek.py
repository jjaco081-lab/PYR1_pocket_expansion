#!/usr/bin/env python
"""
03_parse_foldseek.py -- superfamily conservation at the PYR1 wall positions.

Reads the foldseek hit table from 05_foldseek_search.sh, walks each pairwise
structural alignment, and reports which residue each homolog presents at PYR1
positions 59 / 79 / 94 / 108 / 120 / 141.

alnTM >= 0.5 is used as the helix-grip/SRPBCC core-set cut (CATH 3.30.530.20).

Key question this answers: is the PYR1 wall residue the superfamily norm, or a
PYR/PYL-specific specialisation? A position that is variable across the
superfamily but invariant across close PYR/PYL homologs is under
FAMILY-specific functional constraint, not fold constraint -- the fold will
tolerate changing it, but the switch might not.

Run with the esmfold2 env python.
"""
import collections, csv, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HITS = os.path.join(ROOT, "results", "foldseek", "fs_hits.tsv")
POSITIONS = [59, 79, 94, 108, 120, 141]
WT = dict(zip(POSITIONS, "KREFYE"))
MIN_TM = 0.5

# PYR1 PDB residue number -> 1-based index into the query sequence foldseek saw
seen, nums = set(), []
for line in open(os.path.join(ROOT, "data", "pyr1_A.pdb")):
    if line.startswith("ATOM"):
        rn = int(line[22:26])
        if rn not in seen:
            seen.add(rn); nums.append(rn)
seqidx = {rn: i + 1 for i, rn in enumerate(nums)}
want = {p: seqidx[p] for p in POSITIONS}
print(f"query length {len(nums)}; seq-index: "
      + ", ".join(f"{p}->{want[p]}" for p in POSITIONS))

rows = []
for line in open(HITS):
    f = line.rstrip("\n").split("\t")
    if len(f) < 11:
        continue
    tgt, fid, tm = f[1], float(f[2]), float(f[3])
    qs, qaln, taln = int(f[5]), f[9], f[10]
    qi = qs - 1
    m = {}
    for qc, tc in zip(qaln, taln):
        if qc != "-":
            qi += 1
            for p, si in want.items():
                if qi == si:
                    m[p] = tc
    r = dict(target=tgt, fident=fid, alntm=tm)
    r.update({p: m.get(p, "-") for p in POSITIONS})   # int keys, so not **kwargs
    rows.append(r)

core = [r for r in rows if r["alntm"] >= MIN_TM]
print(f"\n{len(rows)} hits total; {len(core)} at alnTM >= {MIN_TM}\n")

print("=== residue identity at each PYR1 wall position (alnTM >= 0.5) ===")
summary = {}
for p in POSITIONS:
    c = collections.Counter(r[p] for r in core)
    tot = sum(v for k, v in c.items() if k != "-")
    top = [(k, v, 100 * v / max(tot, 1)) for k, v in c.most_common() if k != "-"]
    summary[p] = top
    s = ", ".join(f"{k}:{v}({pc:.0f}%)" for k, v, pc in top[:8])
    print(f"  PYR1 {p:<4}(WT {WT[p]})  {s}   [n={tot}, gaps={c.get('-',0)}]")

print("\n=== close PYR/PYL homologs (alnTM > 0.89) ===")
print(f"{'target':<34}{'fid':>6}{'alnTM':>7}   " + " ".join(f"{p:>4}" for p in POSITIONS))
for r in sorted([r for r in core if r["alntm"] > 0.89], key=lambda x: -x["alntm"]):
    print(f"{r['target']:<34}{r['fident']:>6.2f}{r['alntm']:>7.3f}   "
          + " ".join(f"{r[p]:>4}" for p in POSITIONS))

out = os.path.join(ROOT, "results", "foldseek", "wall_position_identity.csv")
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["target", "fident", "alntm"] + [f"pos{p}" for p in POSITIONS])
    for r in sorted(core, key=lambda x: -x["alntm"]):
        w.writerow([r["target"], r["fident"], r["alntm"]] + [r[p] for p in POSITIONS])
print(f"\nwrote {out}")
