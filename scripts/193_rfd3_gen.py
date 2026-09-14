#!/usr/bin/env python
r"""
193_rfd3_gen.py -- ONE parameterised RFdiffusion3 generator, with a reproducible
and manuscript-reportable sampling scheme.

WHY THIS REPLACES THE PER-ARM SCRIPTS. Arms 183/187/188/189/191 were produced by
sed-copying one another. Two of them failed on defects that only existed in the
copy: `${arm: -1}` wrote `+seed=2026091a` (not an integer), and an over-escaped
line continuation broke the override list. Neither was an RFd3 problem. A single
generator with explicit parameters removes that whole class.

REPRODUCIBLE DIVERSITY, which is Jannis's requirement. Design length and scaffold
segment lengths are the levers that change the envelope, and sampling them widens
the pool. To be reportable the sampling must be exactly replayable, so:

  * ONE integer --master-seed controls everything.
  * Per-design seeds are derived deterministically:
        seed_i = master_seed * 100000 + i
    so design i of a run is reproducible on its own, and two runs with different
    master seeds cannot collide.
  * Each design's sampled parameters are drawn from
        numpy.random.default_rng(seed_i)
    which is version-stable, and are WRITTEN OUT per design to params.jsonl.
  * The manuscript line is then: "designs were generated with master seed S;
    per-design seeds S*1e5+i; segment lengths drawn uniformly from the ranges in
    Table X" -- and params.jsonl is the machine-readable record.

⚠ The sampled ranges are declared here in SEGMENTS, not passed loosely, so the
sampling space itself is versioned with the code.

Usage:
  python 193_rfd3_gen.py --arm 0a --n 10 --master-seed 2026 [--dry-run]
"""
import argparse, json, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

#: fixed motif segments, PYR1 auth numbering, verified against
#: data/3QN1_complex_auth.pdb before every run (see assert_motif)
MOTIF = [("A34-40", "HAQRIHA"), ("A58-65", "YKHFIKSC"),
         ("A81-92", "VIVISGLPANTS"), ("A111-121", "IGGEHRLTNYK"),
         ("A146-168", "DMPEGNSEDDTRMFADTVVKLNL")]
#: HAB1 blocks -- four, because chain B is missing 26 residues
HAB1 = ["B185-221", "B232-270", "B283-461", "B466-505"]
#: scaffold segment ranges, sampled per design. (lo, hi) inclusive.
#: leading range is short and anchored now: the old 50-70 unanchored segment put
#: 46 of the first 60 residues within 4.5 A of HAB1 and hosted K57/K63.
SEGMENTS = [(10, 25), (10, 20), (15, 30), (20, 35), (25, 40), (5, 15)]
#: the 12 pocket-facing motif residues that are NOT also HAB1 interface
POCKET = "A59,A62,A81,A83,A91,A92,A115,A120,A160,A163,A164,A167"
ARMS = {"0a": None, "0b": "partially_buried", "0c": "buried"}
RFD3 = "/bigdata/cutlerlab/jjaco081/conda_envs/foundry/bin/rfd3"
INPUT = os.path.join(ROOT, "data", "3QN1_complex_auth.pdb")


def assert_motif():
    T = {"HIS": "H", "PHE": "F", "ILE": "I", "LYS": "K", "SER": "S", "GLY": "G",
         "LEU": "L", "PRO": "P", "ALA": "A", "ASN": "N", "THR": "T", "ARG": "R",
         "ASP": "D", "MET": "M", "GLU": "E", "VAL": "V", "TYR": "Y", "CYS": "C",
         "GLN": "Q", "TRP": "W"}
    seq = {}
    for l in open(INPUT):
        if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == "A":
            seq[int(l[22:26])] = T.get(l[17:20].strip(), "X")
    for span, expect in MOTIF:
        a, b = (int(x) for x in span[1:].split("-"))
        got = "".join(seq.get(i, "-") for i in range(a, b + 1))
        assert got == expect, f"{span}: frame reads {got}, expected {expect}"
    return len(seq)


def contig(rng):
    """one sampled contig; SEGMENTS are the scaffold, MOTIF the fixed parts"""
    lens = [int(rng.integers(lo, hi + 1)) for lo, hi in SEGMENTS]
    parts = [str(lens[0])]
    for i, (span, _) in enumerate(MOTIF):
        parts.append(span)
        parts.append(str(lens[i + 1]))
    return ",".join(parts) + ",/0," + ",/0,".join(HAB1), lens


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=sorted(ARMS))
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--master-seed", type=int, required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    n_res = assert_motif()
    out = os.path.join(ROOT, "results", "rfd3", f"arm{a.arm}_s{a.master_seed}")
    os.makedirs(out, exist_ok=True)
    rec = open(os.path.join(out, "params.jsonl"), "a")
    print(f"motif asserted against {os.path.basename(INPUT)} ({n_res} chain-A residues)")
    print(f"arm {a.arm}  conditioning={ARMS[a.arm]}  n={a.n}  master_seed={a.master_seed}")

    for i in range(a.n):
        seed = a.master_seed * 100000 + i
        rng = np.random.default_rng(seed)
        cg, lens = contig(rng)
        d = os.path.join(out, f"d{i:03d}")
        if os.path.exists(os.path.join(d, "done")):
            continue
        os.makedirs(d, exist_ok=True)
        cmd = [RFD3, "inputs=null", f"out_dir={d}", "n_batches=1",
               "diffusion_batch_size=1", f"+seed={seed}",
               f"+specification.input={INPUT}", f"+specification.contig={cg}",
               "+specification.dialect=2"]
        if ARMS[a.arm]:
            cmd.insert(-1, f"+specification.select_{ARMS[a.arm]}={POCKET}")
        rec.write(json.dumps(dict(design=i, seed=seed, arm=a.arm,
                                  conditioning=ARMS[a.arm], contig=cg,
                                  segment_lengths=lens)) + "\n")
        rec.flush()
        if a.dry_run:
            print(f"  d{i:03d} seed={seed} lens={lens}")
            continue
        r = subprocess.run(cmd, capture_output=True, text=True)
        ok = bool(__import__("glob").glob(os.path.join(d, "**", "*.cif*"),
                                          recursive=True))
        if ok:
            open(os.path.join(d, "done"), "w").write(str(seed))
        print(f"  d{i:03d} seed={seed} lens={lens} rc={r.returncode} "
              f"{'OK' if ok else 'NO OUTPUT'}", flush=True)
        if not ok and r.returncode != 0:
            print("    " + (r.stderr or r.stdout)[-400:].replace("\n", "\n    "))
    rec.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
