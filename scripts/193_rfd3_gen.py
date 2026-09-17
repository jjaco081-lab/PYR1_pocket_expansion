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

#: --- Stage 1 arms -------------------------------------------------------
#: Stage 0 measured two defects that trade against each other:
#:
#:   * batch 3 (no A34-40 anchor, 50-70 unanchored lead) -> 9/10 ONE intact
#:     chain, 7/10 cavities, but ALL TEN wrap HAB1 (9-21 leading contacts) and
#:     Rg 22-26 A.
#:   * arms 0a/0b (anchor + 10-25 lead) -> HAB1 wrapping fixed decisively
#:     (median 17 -> 2 contacts, p = 0.011) but only 4/10 stay one piece.
#:
#: The breaks localise to lead, link3 and link4 -- and link3/link4 carry the
#: SAME (20,35)/(25,40) ranges as batch 3, which breaks there zero times. So the
#: extra fixed segment (61 motif residues instead of 54) is being paid for by
#: breaking downstream linkers. Two things are therefore untested and both are
#: on the critical path: whether SLACK in link3/link4 buys connectivity back,
#: and whether the ANCHOR or merely the SHORT LEAD causes the fragmentation.
MOTIF_ANCHORED = MOTIF                     # 5 segments, includes A34-40
MOTIF_BARE = MOTIF[1:]                     # 4 segments, batch 3's motif
SEG_BASE = SEGMENTS                                              # arms 0a/0b
SEG_SLACK = [(10, 25), (10, 20), (15, 30), (30, 45), (35, 50), (5, 15)]
SEG_SHORTLEAD = [(10, 25), (15, 30), (20, 35), (25, 40), (5, 15)]
SEG_BATCH3 = [(50, 70), (15, 30), (20, 35), (25, 40), (10, 25)]

#: Jannis: extend the scaffolded grip helix that contacts HAB1 "a little
#: further along the length of the helix", 3-4 extra residues. Measured in 3QN1
#: chain A: A146-168 is the grip helix (74 % helical by P-SEA, 8 of 23 residues
#: within 4.5 A of HAB1, closest 2.57 A) and the helix RUNS ON to residue 177
#: (coil from 178), so a +4 extension to A146-172 stays entirely helical.
#: It also shrinks the free-scaffold fraction, which is what produces the
#: dangling arms -- 65 fixed residues instead of 61.
MOTIF_HELIXEXT = MOTIF[1:-1] + [("A146-172", "DMPEGNSEDDTRMFADTVVKLNLQKLA")]

#: PERMISSIVE distribution. `specification.length` constrains the TOTAL and the
#: per-segment ranges become 0-N, so the model decides WHERE the residues go
#: instead of being told. Verified against the installed parser
#: (foundry/utils/components.py:85 get_design_pattern_with_constraints):
#:   * `length` counts EVERY residue, HAB1 context included -- fixed = 54 PYR1
#:     motif + 295 HAB1 = 349, so a ~190-residue design chain is length ~485,
#:     not 190. Passing 300 raises "No valid selections possible".
#:   * a 0 minimum IS legal, so a flank can be omitted entirely -- but only
#:     26/300 draws contain a zero segment and 1/300 drops the lead, because the
#:     allocator takes randint per segment in order and early segments win. This
#:     LOOSENS the constraint, it does not remove it.
SEG_PERMISSIVE = [(0, 90), (0, 45), (0, 50), (0, 55), (0, 40)]
#: THE UNTESTED CELL. Every anchored arm so far used the SHORT lead and every
#: long-lead arm was unanchored, because the A34-40 anchor and the lead
#: shortening were introduced as ONE change in Stage 0 -- so they could not be
#: told apart. With proper n both help INDEPENDENTLY:
#:   long vs short lead, both unanchored : cavity 50% vs 10%  (p = 8.7e-12)
#:   anchor vs none, both short lead     : cavity 25% vs 10%  (p = 0.0015)
#: and their combination has never been run. SEG_BATCH3 carries 5 entries for a
#: 4-segment motif; the anchored motif has 5 segments and needs 6, so batch 3's
#: long lead (50-70) and trailing (10-25) are kept and SEG_BASE's A34-40 ->
#: A58-65 linker (10-20) is inserted after the lead.
SEG_ANCHLONG = [(50, 70), (10, 20), (15, 30), (20, 35), (25, 40), (10, 25)]

#: ---------------------------------------------------------------------------
#: STAGE 3: ABA IN THE POCKET, AND THE POCKET RESIDUES IN THE MOTIF.
#:
#: ⚠ EVERY DESIGN BEFORE THIS WAS SCAFFOLDED AROUND AN EMPTY POCKET. The input
#: `3QN1_complex_auth.pdb` has no ligand -- ABA was stripped and nobody noticed
#: (README §124). `3QN1_complex_auth_aba.pdb` (built by 218) restores it as 19
#: heavy atoms, CCD A8S, chain L residue 1, with all five motif segments
#: re-asserted by identity.
#:
#: Jannis: "One with also scaffolding ABA, another with adding additional
#: residues to scaffold that define the ABA pocket that are near the opening and
#: are directly next to other sections that we scaffold or can be added through
#: 1-2 extra amino acids so that we do not need to once again define the sizes
#: of more intermediary regions." -- i.e. GROW the existing segments outward
#: rather than add new ones, so no new linker range has to be invented.
#:
#: He also asked to keep R79/E94 in ("these can still be changed later but for
#: now I want to test with ABA") and to RANGE the amount added. Identities were
#: asserted against the input before these were written: 79 IS R, 94 IS E,
#: 108 IS F.
#:
#: ⚠ LINKER RANGES FOLLOW §127, NOT SEG_ANCHLONG. Linker length buys APPENDAGES,
#: not pocket volume (excess vs appendage rho = +0.689; gap4 vs cavity -0.251,
#: q = 0.002), and the mandated minimum exceeded PYR1's native length in exactly
#: gaps 3 and 4. So every gap here is bracketed at (native-4, native+7) rather
#: than the old (20,35)/(25,40).
MOTIF_EXT1 = [("A34-40", "HAQRIHA"), ("A58-65", "YKHFIKSC"),
              ("A80-93", "DVIVISGLPANTST"), ("A110-122", "IIGGEHRLTNYKS"),
              ("A144-168", "VVDMPEGNSEDDTRMFADTVVKLNL")]
MOTIF_EXT2 = [("A34-40", "HAQRIHA"), ("A58-65", "YKHFIKSC"),
              ("A79-94", "RDVIVISGLPANTSTE"), ("A109-122", "SIIGGEHRLTNYKS"),
              ("A142-168", "SYVVDMPEGNSEDDTRMFADTVVKLNL")]
MOTIF_EXT3 = [("A34-40", "HAQRIHA"), ("A58-65", "YKHFIKSC"),
              ("A79-94", "RDVIVISGLPANTSTE"), ("A108-122", "FSIIGGEHRLTNYKS"),
              ("A141-168", "ESYVVDMPEGNSEDDTRMFADTVVKLNL")]


def _seg_for(motif, lead=(50, 70), tail=(10, 25)):
    """Linker ranges bracketing each gap's NATIVE length, per §127.

    Returns lead + one range per inter-segment gap + tail. The native gap is
    read from the motif's own residue numbers, so it cannot drift out of sync
    with the segments the way a hand-written table can.
    """
    out = [lead]
    for a, b in zip(motif, motif[1:]):
        hi_a = int(a[0].split("-")[1])
        lo_b = int(b[0][1:].split("-")[0])
        nat = lo_b - hi_a - 1
        out.append((max(4, nat - 4), nat + 7))
    out.append(tail)
    return out


SEG_EXT1, SEG_EXT2, SEG_EXT3 = (_seg_for(m) for m in
                                (MOTIF_EXT1, MOTIF_EXT2, MOTIF_EXT3))
#: the ABA-only arm: 8best's motif and ranges, but §127 linkers and a ligand.
SEG_ABA = _seg_for(MOTIF)

#: arms whose input carries ABA and that pass specification.ligand.
LIGAND_ARMS = {"9aba", "10ext1", "10ext2", "10ext3"}
LIGAND_CCD = "A8S"
INPUT_ABA = os.path.join(ROOT, "data", "3QN1_complex_auth_aba.pdb")

#: arm -> specification.length. None means unconstrained (every other arm).
LENGTHS = {"7free": "465-505"}

#: (motif, segment ranges, RASA conditioning). len(segs) == len(motif) + 1.
RECIPES = {
    "0a":      (MOTIF_ANCHORED, SEG_BASE,      None),
    "0b":      (MOTIF_ANCHORED, SEG_BASE,      "partially_buried"),
    "0c":      (MOTIF_ANCHORED, SEG_BASE,      "buried"),
    # Stage 1: s1 vs s4 is the slack test; s2 vs s4 isolates the anchor from
    # the short lead; s3 puts RASA on batch 3's good-connectivity contig.
    "1slack":  (MOTIF_ANCHORED, SEG_SLACK,     "partially_buried"),
    "2noanch": (MOTIF_BARE,     SEG_SHORTLEAD, "partially_buried"),
    "3b3rasa": (MOTIF_BARE,     SEG_BATCH3,    "partially_buried"),
    "4rep":    (MOTIF_ANCHORED, SEG_BASE,      "partially_buried"),
    # Stage 2: the winning long-lead recipe, and the helix-extension variant.
    "5long":   (MOTIF_BARE,      SEG_BATCH3,    "partially_buried"),
    "6helix":  (MOTIF_HELIXEXT,  SEG_BATCH3,    "partially_buried"),
    # one-variable test against 5long: same motif, same total length envelope,
    # permissive distribution instead of mandated flanks.
    "7free":   (MOTIF_BARE,      SEG_PERMISSIVE, "partially_buried"),
    # anchor for the sheet nucleus + long lead for the pocket; the tail it
    # reintroduces is handled post hoc by truncation (35/40 rescued, cavity
    # survived 40/40), which is how batch 3 reaches 30% PASS.
    "8best":   (MOTIF_ANCHORED,  SEG_ANCHLONG,   "partially_buried"),
    # Stage 3: ABA in the pocket (9aba), then the pocket residues progressively
    # folded into the motif (10ext1/2/3). ext3 includes R79, E94 and F108.
    "9aba":    (MOTIF_ANCHORED,  SEG_ABA,        "partially_buried"),
    "10ext1":  (MOTIF_EXT1,      SEG_EXT1,       "partially_buried"),
    "10ext2":  (MOTIF_EXT2,      SEG_EXT2,       "partially_buried"),
    "10ext3":  (MOTIF_EXT3,      SEG_EXT3,       "partially_buried"),
}
ARMS = {k: v[2] for k, v in RECIPES.items()}
for _k, (_m, _s, _c) in RECIPES.items():
    assert len(_s) == len(_m) + 1, f"{_k}: {len(_s)} segments for {len(_m)} motifs"
RFD3 = "/bigdata/cutlerlab/jjaco081/conda_envs/foundry/bin/rfd3"
INPUT = os.path.join(ROOT, "data", "3QN1_complex_auth.pdb")


def assert_motif(motif=None):
    T = {"HIS": "H", "PHE": "F", "ILE": "I", "LYS": "K", "SER": "S", "GLY": "G",
         "LEU": "L", "PRO": "P", "ALA": "A", "ASN": "N", "THR": "T", "ARG": "R",
         "ASP": "D", "MET": "M", "GLU": "E", "VAL": "V", "TYR": "Y", "CYS": "C",
         "GLN": "Q", "TRP": "W"}
    seq = {}
    for l in open(INPUT):
        if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == "A":
            seq[int(l[22:26])] = T.get(l[17:20].strip(), "X")
    for span, expect in (motif or MOTIF):
        a, b = (int(x) for x in span[1:].split("-"))
        got = "".join(seq.get(i, "-") for i in range(a, b + 1))
        assert got == expect, f"{span}: frame reads {got}, expected {expect}"
    return len(seq)


def contig(rng, motif=None, segs=None, ranges=False):
    """one sampled contig; segs are the scaffold, motif the fixed parts.

    ⚠ Normally this PRE-SAMPLES each segment length here, so the contig handed
    to RFd3 carries concrete numbers and the sampling is ours (reproducible from
    master_seed). That is incompatible with `specification.length`, which asks
    RFd3's own allocator to distribute a total across the segments -- concrete
    numbers leave it nothing to allocate. With ranges=True the contig carries
    the RANGES instead and RFd3 does the distribution.
    """
    motif, segs = motif or MOTIF, segs or SEGMENTS
    if ranges:
        lens = [f"{lo}-{hi}" for lo, hi in segs]
    else:
        lens = [int(rng.integers(lo, hi + 1)) for lo, hi in segs]
    parts = [str(lens[0])]
    for i, (span, _) in enumerate(motif):
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

    motif, segs, cond = RECIPES[a.arm]
    n_res = assert_motif(motif)
    out = os.path.join(ROOT, "results", "rfd3", f"arm{a.arm}_s{a.master_seed}")
    os.makedirs(out, exist_ok=True)
    rec = open(os.path.join(out, "params.jsonl"), "a")
    print(f"motif asserted against {os.path.basename(INPUT)} ({n_res} chain-A residues)")
    print(f"arm {a.arm}  conditioning={cond}  n={a.n}  "
          f"master_seed={a.master_seed}")
    print(f"  motif   {[m[0] for m in motif]}")
    print(f"  segments {segs}")

    for i in range(a.n):
        seed = a.master_seed * 100000 + i
        rng = np.random.default_rng(seed)
        cg, lens = contig(rng, motif, segs, ranges=bool(LENGTHS.get(a.arm)))
        d = os.path.join(out, f"d{i:03d}")
        if os.path.exists(os.path.join(d, "done")):
            continue
        os.makedirs(d, exist_ok=True)
        length = LENGTHS.get(a.arm)
        cmd = [RFD3, "inputs=null", f"out_dir={d}", "n_batches=1",
               "diffusion_batch_size=1", f"+seed={seed}",
               f"+specification.input={INPUT_ABA if a.arm in LIGAND_ARMS else INPUT}",
               # Hydra treats a bare comma-separated value as ambiguous
               # ("To use it as string, quote the value"). subprocess does not go
               # through a shell, so the quotes must be literal characters here.
               f"+specification.contig='{cg}'",
               "+specification.dialect=2"]
        if length:
            cmd.insert(-1, f"+specification.length='{length}'")
        if a.arm in LIGAND_ARMS:
            cmd.insert(-1, f"+specification.ligand='{LIGAND_CCD}'")
        if cond:
            cmd.insert(-1, f"+specification.select_{cond}='{POCKET}'")
        rec.write(json.dumps(dict(design=i, seed=seed, arm=a.arm,
                                  conditioning=cond, contig=cg,
                                  motif=[m[0] for m in motif],
                                  segment_ranges=[list(x) for x in segs],
                                  segment_lengths=lens,
                                  total_length=LENGTHS.get(a.arm))) + "\n")
        rec.flush()
        if a.dry_run:
            print(f"  d{i:03d} seed={seed} lens={lens}")
            print("       " + " ".join(cmd[1:]))
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
