#!/usr/bin/env python
r"""
220_core_truncate.py -- THE truncation rule: keep what packs against the MOTIF
core, cut the rest from the termini.

TWO EARLIER RULES, BOTH WRONG, AND JANNIS DIAGNOSED BOTH BY EYE.

  HAB1-contact rule: delete the shortest terminal piece containing every
  HAB1-contacting residue. FAILS when an appendage does not reach HAB1 --
  "blue gets it right on the C terminus" of arm6helix d001.

  local-packing rule: trim termini while poorly packed against the design.
  FAILS when the appendage packs against ITSELF -- on arm3b3rasa d028 "the
  protein randomly places two alpha helices that pack against each other on the
  other side of HAB1", which are densely packed and so were kept.

THE FIX is to measure packing against the MOTIF CORE rather than against the
design as a whole. Flood-fill outward from the motif: a residue joins the core
if it contacts >=K core residues within NB A. Two helices packing against each
other are never reached, because they do not touch the motif-containing domain.
HAB1 is not mentioned anywhere in the criterion.

CALIBRATION, against Jannis's own visual call. He specified arm6helix d000
should run 52-208. K=4 contacts within 6.0 A gives 51-209 -- one residue out at
each end. On the other three it reproduces his preference: d028 cuts the
self-packing helices (61 residues, where local-packing cut 5), d074 cuts BOTH
termini (42 N and 22 C, where the HAB1 rule missed the C-term entirely), and
d001 cuts 48 at the N-terminus as he preferred.

⚠ Tuned on ONE design and sanity-checked on three. That is a small validation
set and the constants are the weakest part of this.

⚠ Only TERMINAL residues are removed. An appendage on an inter-motif linker
cannot be cut without breaking the chain; those designs are flagged, not
trimmed. The real fix for those is upstream -- we mandate 20-35 and 25-40
residue linkers to span gaps of 7.9 and 9.9 A, and the excess has to go
somewhere.

⚠ HAB1 is written as its FOUR separate chains. An earlier version flattened them
into one, so their residue numbers collided (1 1 2 2 3 3 ...) and PyMOL drew
disconnected lines.
"""
import glob, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m195 = import_module("195_arm_compare")
m198 = import_module("198_cascade2")
K_CONTACT, NB_DIST, SEQ_SEP = 4, 6.0, 3
#: minimum run of consecutive non-core residues before it counts as an INTERNAL
#: APPENDAGE. Jannis: "sometimes the orange are just 1 amino acid stretches and
#: I am not convinced that these are an excessive internal appendage." Correct --
#: an isolated residue marginally below threshold is not a structural defect,
#: and flagging it made every one of five designs look afflicted.
MIN_INTERNAL_RUN = 4

#: ⚠ RE-TUNE AFTER SEQUENCE DESIGN. Jannis: "the packing might change slightly
#: once the glycines all get side chains but that is unknown." He is right and it
#: cuts one way: RFd3 emits a low-complexity placeholder sequence, so contacts
#: are counted on backbone plus a stub. Real side chains can only ADD contacts,
#: so this rule is currently biased toward cutting slightly too much, and
#: K_CONTACT should be re-checked on LigandMPNN output before it is trusted
#: downstream of Stage 3.


def load_design(cif):
    at = m195.read_cif_gz(cif)
    imap = json.load(open(cif.replace(".cif.gz", ".json")))["diffused_index_map"]
    dc = m195.design_chain(at, imap)
    dat = [a for a in at if a[0] == dc]
    ca, by = {}, {}
    for a in dat:
        by.setdefault(a[1], []).append(a[4])
        if a[2] == "CA":
            ca[a[1]] = a[4]
    motif = {int(v[1:]) for v in imap.values() if v[0] == dc}
    main = set(max(m198.pieces_of(ca), key=len))
    return at, dat, dc, by, main, motif & main


def core_of(by, main, motif, K=K_CONTACT, nb=NB_DIST):
    P = {k: np.array(by[k]) for k in main if k in by}
    core = set(motif) & set(P)
    changed = True
    while changed:
        changed = False
        for k in sorted(P):
            if k in core:
                continue
            n = 0
            for j in core:
                if abs(j - k) < SEQ_SEP:
                    continue
                if np.linalg.norm(P[j][None, :, :] - P[k][:, None, :],
                                  axis=2).min() < nb:
                    n += 1
                    if n >= K:
                        break
            if n >= K:
                core.add(k); changed = True
    return core


def truncate(by, main, motif):
    """returns (kept, n_cut_N, n_cut_C, has_internal_gap)"""
    core = core_of(by, main, motif)
    if not core:
        return set(main), 0, 0, False
    lo, hi = min(core), max(core)
    kept = {k for k in main if lo <= k <= hi}
    # longest CONSECUTIVE run of non-core residues inside the span
    best = cur = 0
    for k in range(lo, hi + 1):
        if k in kept and k not in core:
            cur += 1; best = max(best, cur)
        else:
            cur = 0
    return (kept, sum(1 for k in main if k < lo),
            sum(1 for k in main if k > hi), best)


def write_pdb(at, dc, kept, core, path, tag):
    others = sorted({a[0] for a in at if a[0] != dc})
    cm = {dc: "A"}
    cm.update({c: ch for c, ch in zip(others, "BCDE")})
    n = 0
    with open(path, "w") as fh:
        fh.write(f"REMARK  {tag}\n")
        fh.write("REMARK  chain A = design; chains B-E = the four HAB1 blocks\n")
        fh.write("REMARK  B-factor: 0 kept core, 1 kept but not core, 2 TRUNCATED\n")
        fh.write("REMARK  PyMOL: spectrum b, grey70 orange red, chain A, 0, 2\n")
        for a in at:
            if a[0] == dc:
                b = 0.0 if a[1] in core else (1.0 if a[1] in kept else 2.0)
            else:
                b = 0.0
            n += 1
            fh.write(f"ATOM  {n:>5} {a[2]:<4}{'':1}{a[5]:>3} {cm[a[0]]}"
                     f"{a[1] % 10000:>4}{'':1}   {a[4][0]:>8.3f}{a[4][1]:>8.3f}"
                     f"{a[4][2]:>8.3f}  1.00{b:>6.2f}          {a[3]:>2}\n")
        fh.write("END\n")


def main():
    out = os.path.join(ROOT, "results", "rfd3", "truncation_compare")
    os.makedirs(out, exist_ok=True)
    PICKS = [("arm3b3rasa_s3003", "d028"), ("arm3b3rasa_s3003", "d074"),
             ("arm6helix_s6001", "d000"), ("arm6helix_s6001", "d001"),
             ("arm1slack_s3001", "d172")]
    print(f"core-growth truncation  (K={K_CONTACT} contacts within {NB_DIST} Å)\n")
    print(f"{'design':<22}{'chain':>12}{'core':>12}{'cutN':>6}{'cutC':>6}"
          f"{'kept':>6}{'longest':>10}")
    for arm, d in PICKS:
        g = glob.glob(os.path.join(ROOT, "results", "rfd3", arm, d, "*.cif.gz"))
        if not g:
            continue
        at, dat, dc, by, main, motif = load_design(g[0])
        core = core_of(by, main, motif)
        kept, cN, cC, run = truncate(by, main, motif)
        tag = f"{arm.split('_')[0]}_{d}"
        write_pdb(at, dc, kept, core, os.path.join(out, f"{tag}_core.pdb"),
                  f"{tag}  core-growth truncation")
        print(f"{tag:<22}{f'{min(main)}-{max(main)}':>12}"
              f"{f'{min(core)}-{max(core)}':>12}{cN:>6}{cC:>6}{len(kept):>6}"
              f"{(str(run) if run >= MIN_INTERNAL_RUN else '-'):>10}")
    print(f"\nwrote *_core.pdb to {out}/")
    print("  PyMOL:  spectrum b, grey70 orange red, chain A, 0, 2")
    print("    grey = core    orange = kept, not core    red = truncated")
    return 0


if __name__ == "__main__":
    sys.exit(main())
