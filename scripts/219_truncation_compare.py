#!/usr/bin/env python
r"""
219_truncation_compare.py -- write designs with BOTH truncation rules painted
into the B-factor, so the cut can be judged by eye in PyMOL.

Jannis cannot read residue numbers in his viewer, so the rules are encoded as
colour instead:

    B = 0   kept by both rules
    B = 1   cut by the HAB1-CONTACT rule only
    B = 2   cut by the PACKING rule only
    B = 3   cut by BOTH

    PyMOL:  spectrum b, grey80 red blue magenta, all, 0, 3

THE TWO RULES

  HAB1-contact (what 217 currently does): from each terminus, delete the
  shortest piece containing every HAB1-contacting residue in that segment.
  ⚠ Defined by HAB1 rather than by the design, so it leaves dangling residues
  that happen not to touch HAB1, and it only works because the appendages seen
  so far reach HAB1 at all.

  PACKING (Jannis's idea, made robust): "truncate at the start of the loop
  region that leads away from the core". Rather than secondary structure --
  which he correctly noted fails when a small helix sits in the tail -- use
  how well each terminal residue is PACKED against the design itself. A residue
  in the fold has many non-local neighbours; a dangling tail residue has few,
  AND SO DOES A SMALL HELIX FLAPPING AROUND HAB1, because it is not packed
  against the core either. Trim inward from each terminus while the neighbour
  count stays below threshold; stop at the first well-packed residue.
  Independent of HAB1, removes all dangling tail, keeps tail that is real.

Neighbour count excludes |i-j| < 3 so sequence-local contacts cannot make a
dangling residue look packed. The threshold is calibrated per design from the
distribution over the motif residues, which are by definition part of the fold.
"""
import glob, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m195 = import_module("195_arm_compare")
m198 = import_module("198_cascade2")
OUT = os.path.join(ROOT, "results", "rfd3", "truncation_compare")
CONTACT, NB_CUT, SEQ_SEP = 4.5, 8.0, 3


def neighbours(by, main):
    """non-local packing count per residue, within the design only"""
    ks = sorted(main)
    P = {k: np.array(by[k]) for k in ks if k in by}
    out = {}
    for k in ks:
        if k not in P:
            continue
        n = 0
        for j in ks:
            if j == k or abs(j - k) < SEQ_SEP or j not in P:
                continue
            if np.linalg.norm(P[j][None, :, :] - P[k][:, None, :],
                              axis=2).min() < NB_CUT:
                n += 1
        out[k] = n
    return out


def rule_hab1(main, motif, touch):
    lo, hi = min(motif), max(motif)
    cut = set()
    lead_hit = [r for r in main if r < lo and r in touch]
    tail_hit = [r for r in main if r > hi and r in touch]
    if lead_hit:
        cut |= {r for r in main if r < lo and r <= max(lead_hit)}
    if tail_hit:
        cut |= {r for r in main if r > hi and r >= min(tail_hit)}
    return cut


def rule_packing(main, motif, nb, thresh):
    """trim inward from each terminus while poorly packed; stop at the core"""
    ks = sorted(main)
    cut = set()
    for seq in (ks, ks[::-1]):
        for k in seq:
            if k in motif:
                break
            if nb.get(k, 99) < thresh:
                cut.add(k)
            else:
                break
    return cut


def process(cif, tag):
    at = m195.read_cif_gz(cif)
    imap = json.load(open(cif.replace(".cif.gz", ".json")))["diffused_index_map"]
    dc = m195.design_chain(at, imap)
    dat = [a for a in at if a[0] == dc]
    hab1 = np.array([a[4] for a in at if a[0] != dc])
    ca, by = {}, {}
    for a in dat:
        by.setdefault(a[1], []).append(a[4])
        if a[2] == "CA":
            ca[a[1]] = a[4]
    motif = {int(v[1:]) for v in imap.values() if v[0] == dc}
    main = set(max(m198.pieces_of(ca), key=len))
    motif &= main
    touch = set()
    for r in main:
        X = np.array(by[r])
        if len(hab1) and np.linalg.norm(hab1[None, :, :] - X[:, None, :],
                                        axis=2).min() < CONTACT:
            touch.add(r)
    nb = neighbours(by, main)
    mot_nb = [nb[k] for k in motif if k in nb]
    thresh = max(4, int(np.percentile(mot_nb, 10))) if mot_nb else 6
    A = rule_hab1(main, motif, touch)
    B = rule_packing(main, motif, nb, thresh)
    code = {}
    for r in main:
        code[r] = (1 if r in A else 0) + (2 if r in B else 0)
    path = os.path.join(OUT, f"{tag}.pdb")
    n = 0
    with open(path, "w") as fh:
        fh.write(f"REMARK  B=0 kept  B=1 HAB1-rule only  B=2 packing-rule only  "
                 f"B=3 both\n")
        fh.write(f"REMARK  packing threshold {thresh} neighbours "
                 f"(10th pct of motif residues)\n")
        fh.write(f"REMARK  PyMOL: spectrum b, grey80 red blue magenta, all, 0, 3\n")
        for a in at:
            if a[0] == dc and a[1] not in main:
                continue
            n += 1
            ch = "A" if a[0] == dc else "H"
            b = float(code.get(a[1], 0)) if a[0] == dc else 0.0
            fh.write(f"ATOM  {n:>5} {a[2]:<4}{'':1}{a[5]:>3} {ch}{a[1] % 10000:>4}"
                     f"{'':1}   {a[4][0]:>8.3f}{a[4][1]:>8.3f}{a[4][2]:>8.3f}"
                     f"  1.00{b:>6.2f}          {a[3]:>2}\n")
        fh.write("END\n")
    return dict(tag=tag, n_main=len(main), thresh=thresh,
                hab1_cut=len(A), pack_cut=len(B), both=len(A & B),
                only_hab1=len(A - B), only_pack=len(B - A),
                internal_appendage=int(any(
                    r in touch and r not in motif and
                    min(motif) < r < max(motif) for r in main)),
                path=path)


def main():
    os.makedirs(OUT, exist_ok=True)
    PICKS = [("arm3b3rasa_s3003", "d028"), ("arm3b3rasa_s3003", "d074"),
             ("arm6helix_s6001", "d000"), ("arm6helix_s6001", "d001"),
             ("arm1slack_s3001", "d172")]
    print(f"{'design':<26}{'res':>5}{'thr':>5}{'HAB1cut':>9}{'PACKcut':>9}"
          f"{'both':>6}{'onlyH':>7}{'onlyP':>7}{'internal':>9}")
    rows = []
    for arm, d in PICKS:
        g = glob.glob(os.path.join(ROOT, "results", "rfd3", arm, d, "*.cif.gz"))
        if not g:
            print(f"  {arm}/{d}: missing"); continue
        r = process(g[0], f"{arm.split('_')[0]}_{d}")
        rows.append(r)
        print(f"{r['tag']:<26}{r['n_main']:>5}{r['thresh']:>5}{r['hab1_cut']:>9}"
              f"{r['pack_cut']:>9}{r['both']:>6}{r['only_hab1']:>7}"
              f"{r['only_pack']:>7}{'YES' if r['internal_appendage'] else '-':>9}")
    json.dump(rows, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    print(f"\nwrote {len(rows)} structures to {OUT}/")
    print("  PyMOL:  spectrum b, grey80 red blue magenta, all, 0, 3")
    print("    grey = kept    red = HAB1-rule cuts    blue = packing-rule cuts"
          "    magenta = both")
    return 0


if __name__ == "__main__":
    sys.exit(main())
