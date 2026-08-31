#!/usr/bin/env python
r"""
159_dsm_doubles.py -- Beltran's ACTUAL WIN library, enumerated and scored.

The WIN paper's methods define the search space exactly: a DOUBLE-mutant library
at 19 positions, 15 free to mutate to anything but cysteine and proline, 4
(62, 81, 87, 110) restricted to smaller subsets. sd04 records that library as
`DSM-Hao` with the per-position amino-acid menus, at 18 of those 19 positions
(62 is absent from sd04's version). Its allowed doubles number **35,863**, which
matches the ~36,140 of §53.

Round 1 of that library yielded five WIN 55,212-2 hits, and every one is a
double at 159+160:

    F159G+A160V   F159S+A160L   F159G+A160I   F159S+A160V   F159T+A160L

all at ~10 µM, before the CB-S shuffles took it to K59Q+F159A+A160I (50 nM) and
finally E4G+Y23H+D26G+K59Q+F159A+A160I (10 nM). So the prediction task is sharp
and closed: **of 35,863 allowed doubles, where do those five rank?**

FRAME: 7MWN, REVERTED TO WILD TYPE. 7MWN is a PYL2-based WIN sensor with WI5
bound — the only crystal structure of a WIN-bound receptor, and the "exact
coordinates" this test needs. Chain A carries K64Q/F165A/V166I, so it is
reverted to wild-type PYL2 first; otherwise the background already contains the
answer at two of the three positions. The WT receptor then clashes with the
crystal ligand pose, which is the point: the question is which double relieves
it.

⚠ PYL2 IS NOT PYR1. The DSM library is defined on PYR1 numbering, so positions
are mapped by global sequence alignment and the map is ASSERTED against the
three correspondences the paper states independently (PYR1 K59/F159/A160 ↔ PYL2
K64/F165/V166, from its PYL2^WIN construct). A mapping that fails those three
stops the run. The receptor is still a homolog, and that is the main caveat on
any result here.

⚠ Cross-reactivity is not modelled (Jannis): a high rank is not proof a variant
would have been selected, only that the score sees it.

SCORING is §78's licence: repack only, ligand pose fixed, dG_bind by the
three-term split, PAIRED wild-type control repacked in the IDENTICAL shell. The
WT control depends only on the shell, and shells repeat across all ~289
substitution pairs at a given position pair, so it is cached per position pair —
that is what makes 35,863 doubles affordable.
"""
import argparse
import itertools
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import lib_rosetta as LR                                          # noqa: E402
from lib_xlsx import table                                        # noqa: E402

SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "dsm_doubles")
WIN = os.path.join(ROOT, "results", "win_crossover")
SENSOR = os.path.join(WIN, "sensor.pdb")
LIGPDB = os.path.join(WIN, "WI5_0001.pdb")
LIGPRM = os.path.join(WIN, "WI5.params")
#: 7MWN chain A is the sensor; revert to wild-type PYL2 (paper's own annotation)
REVERT = {64: ("Q", "LYS"), 165: ("A", "PHE"), 166: ("I", "VAL")}
#: stated independently by the paper's PYL2^WIN construct -- the map must satisfy these
ANCHORS = {59: 64, 159: 165, 160: 166}
THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
         "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
         "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
         "W": "TRP", "Y": "TYR"}


def dsm_menu():
    _, recs = table(f"{SD}/pnas.2519924122.sd04.xlsx")
    m = {}
    for r in recs:
        if (r.get("Library") or "").strip() == "DSM-Hao":
            m[int(float(r["Position"]))] = (
                r["WT"].strip(),
                "".join(sorted(set(r["Amino Acids allowed for mutation"].strip()))))
    return m


def nw(a, b):
    """Global alignment, identity scoring. -> dict a_index -> b_index (0-based)."""
    n, m = len(a), len(b)
    g = -4
    S = np.zeros((n + 1, m + 1))
    S[:, 0] = np.arange(n + 1) * g
    S[0, :] = np.arange(m + 1) * g
    P = np.zeros((n + 1, m + 1), int)
    P[1:, 0] = 1
    P[0, 1:] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d = S[i - 1, j - 1] + (2 if a[i - 1] == b[j - 1] else -1)
            u, l = S[i - 1, j] + g, S[i, j - 1] + g
            k = int(np.argmax([d, u, l]))
            S[i, j] = [d, u, l][k]
            P[i, j] = k
    i, j, out = n, m, {}
    while i > 0 or j > 0:
        if P[i, j] == 0 and i and j:
            out[i - 1] = j - 1
            i, j = i - 1, j - 1
        elif P[i, j] == 1 and i:
            i -= 1
        else:
            j -= 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--nchunks", type=int, default=1)
    ap.add_argument("--nrep", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    import pyrosetta
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    pyrosetta.init(f"-mute all -ignore_unrecognized_res -ex1 -ex2aro "
                   f"-extra_res_fa {LIGPRM}")
    sfxn = pyrosetta.create_score_function("ref2015_cart")
    pose = pyrosetta.pose_from_pdb(SENSOR)
    info = pose.pdb_info()
    idx = {info.number(i): i for i in range(1, pose.total_residue() + 1)}
    lig = pose.total_residue()

    # ---- revert 7MWN chain A to wild-type PYL2 -------------------------
    for num, (have, want3) in REVERT.items():
        got = pose.residue(idx[num]).name1()
        assert got == have, f"7MWN position {num} is {got}, expected sensor {have}"
        MutateResidue(idx[num], want3).apply(pose)
    print(f"7MWN chain A reverted to wild-type PYL2 at {sorted(REVERT)}", flush=True)

    # ---- map PYR1 DSM positions onto this frame, then ASSERT ------------
    pyl2 = "".join(pose.residue(i).name1() for i in range(1, lig))
    pnums = [info.number(i) for i in range(1, lig)]
    pyr1_frame = os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")
    seq1, num1 = [], []
    for line in open(pyr1_frame):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            r3 = line[17:20].strip()
            if r3 in THREE.values():
                seq1.append(next(k for k, v in THREE.items() if v == r3))
                num1.append(int(line[22:26]))
    amap = nw("".join(seq1), pyl2)
    p2l = {num1[i]: pnums[j] for i, j in amap.items()}
    bad = {p: (p2l.get(p), q) for p, q in ANCHORS.items() if p2l.get(p) != q}
    assert not bad, f"PYR1->PYL2 map fails the paper's own anchors: {bad}"
    print(f"alignment maps {len(p2l)} PYR1 positions; anchors "
          f"{ANCHORS} CONFIRMED", flush=True)

    menu = dsm_menu()
    menu = {p: v for p, v in menu.items() if p in p2l}
    jobs = []
    for pa, pb in itertools.combinations(sorted(menu), 2):
        for aa in menu[pa][1]:
            for bb in menu[pb][1]:
                jobs.append((pa, aa, pb, bb))
    print(f"DSM-Hao at {len(menu)} mapped positions -> {len(jobs)} allowed doubles",
          flush=True)
    jobs = [j for i, j in enumerate(jobs) if i % a.nchunks == a.chunk]
    print(f"chunk {a.chunk}/{a.nchunks}: {len(jobs)}", flush=True)

    e_lig = float(sfxn(pyrosetta.pose_from_pdb(LIGPDB)))

    def dg(p):
        e_c = float(sfxn(p))
        q = p.clone()
        q.delete_residue_slow(q.total_residue())
        return e_c - float(sfxn(q)) - e_lig

    def repack(p0, allowed, nrep):
        v = []
        for _ in range(nrep):
            p = p0.clone()
            tf, _ = LR.restrict_packing(p, allowed)
            pk = PackRotamersMover(sfxn)
            pk.task_factory(tf)
            pk.apply(p)
            v.append(dg(p))
        return min(v), float(max(v) - min(v))

    shell_cache, wt_cache, rows = {}, {}, []
    for pa, aa, pb, bb in jobs:
        key = (pa, pb)
        if key not in shell_cache:
            sel = ResidueIndexSelector(f"{idx[p2l[pa]]},{idx[p2l[pb]]}")
            nb = NeighborhoodResidueSelector(sel, 6.0, True)
            shell_cache[key] = sorted({i + 1 for i, b in enumerate(nb.apply(pose))
                                       if b} | {lig})
        allowed = shell_cache[key]
        if key not in wt_cache:
            wt_cache[key] = repack(pose.clone(), allowed, a.nrep)[0]
        p = pose.clone()
        MutateResidue(idx[p2l[pa]], THREE[aa]).apply(p)
        MutateResidue(idx[p2l[pb]], THREE[bb]).apply(p)
        s, spread = repack(p, allowed, a.nrep)
        rows.append({"sub": f"{menu[pa][0]}{pa}{aa}+{menu[pb][0]}{pb}{bb}",
                     "pa": pa, "aa": aa, "pb": pb, "bb": bb,
                     "dG": s, "wt": wt_cache[key], "ddG": s - wt_cache[key],
                     "spread": spread, "n_shell": len(allowed)})
        if len(rows) % 200 == 0:
            json.dump(rows, open(os.path.join(OUT, f"chunk_{a.chunk:03d}.json"), "w"))
            print(f"    {len(rows)}/{len(jobs)}", flush=True)
    json.dump(rows, open(os.path.join(OUT, f"chunk_{a.chunk:03d}.json"), "w"))
    print(f"wrote {len(rows)} to chunk_{a.chunk:03d}.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
