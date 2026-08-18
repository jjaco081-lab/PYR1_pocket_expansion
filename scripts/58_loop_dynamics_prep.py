#!/usr/bin/env python
"""
58_loop_dynamics_prep.py -- establish and VERIFY the residue-numbering map for the
gate/latch loop-dynamics analysis, PER ANALYSIS UNIT, then hand the masks to 58b.

WHY THIS SCRIPT EXISTS SEPARATELY FROM THE ANALYSIS
---------------------------------------------------
tleap renumbers residues sequentially 1..N in the prmtop. The gate (85-89) and
latch (115-117) are quoted in NATIVE PYR1 numbering everywhere in this project,
and the crystal files have gaps, so `:85-89` in a cpptraj mask does NOT select the
gate. Getting this wrong is silent: cpptraj returns RMSD values for the wrong loop
and every downstream conclusion is confidently wrong.

AND THE MAP IS NOT THE SAME FOR EVERY UNIT. The offsets actually in play:

    S1_apo_open / S2_holo_closed   178 residues, chain A, gate = sequential 82-86
    S4_ternary                     179 residues, chain A, gate = sequential 83-87
    S3 chain A                     183 residues, prmtop 1-183,   gate = 85-89
    S3 chain B                     183 residues, prmtop 184-366, gate = 268-272

S1/S2 come from script 30, which took the INTERSECTION of 3K3K and 3QN1 and so
dropped residue 2 (disordered in 3K3K). S4 was built straight from 3QN1, which HAS
residue 2 -- as ALA, the P2A that README 23 records 3QN1 handing to every pipeline
in this project. S3 was built from 3K3K directly and keeps both protomers, so its
chain A is gap-free 1-183 and its chain B is 2-184 (one residue fewer at the N
terminus, one more at the C terminus -- verified, not assumed).

WHY S3 IS ANALYSED AS TWO UNITS
-------------------------------
3K3K is a MIXED dimer (README 28g): chain A is apo-open, but chain B is CLOSED and
was ABA-bound in the crystal. S3 was built protein-only, so its chain B is a
closed, ligand-shaped protomer simulated around an EMPTY pocket -- which is the
apo-closed cell that README 24e names as the missing control. The two protomers
are therefore different experiments in one box and are never pooled.

WHAT IT CHECKS
--------------
  1. per unit: sequential index -> native resnum, built by ORDER across the whole
     protein (so a second chain gets the right offset), then verified by RESNAME
     at a dozen known positions
  2. the same landmarks are re-verified against the PRMTOP's own RESIDUE_LABEL
     block -- the file cpptraj will actually read -- so a protein.pdb that
     disagrees with the topology cannot pass
  3. the gate, latch and Lbeta7-alpha5 loops carry their expected sequences
     (gate = SER GLY LEU PRO ALA, latch = HIS ARG LEU) in every unit
  4. ONE common core superposition set, the intersection over every unit and the
     references, asserted identical everywhere. S3 chain B lacks native residue 1,
     so the core is now 160 residues where the S1/S2/S4-only run used 161; every
     unit is refitted on the same set rather than letting the sets drift apart.
  5. the open/closed references really are an open/closed pair, by recovering a
     ~5 A gate/latch stroke. A MAGNITUDE check only.

OBSERVABLES (pre-registered in README 19d, fixed before the runs finished)
-------------------------------------------------------------------------
  - gate (85-89) and latch (115-117) backbone RMSD and RMSF
  - gate-latch minimum heavy-atom distance
  - Lbeta7-alpha5 (148-156) RMSF          [Dorosh 2013, README 17a]

Nothing outside that list is computed, so the answer cannot be shopped for.

Run with the pyr1_docking env python.
"""
import json
import os
import re

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
MD = os.path.join(DATA, "md")
OUT = os.path.join(DATA, "loop_dynamics")
os.makedirs(OUT, exist_ok=True)

GATE = list(range(85, 90))        # native PYR1 numbering
LATCH = list(range(115, 118))
LB7A5 = list(range(148, 157))
BB = ("N", "CA", "C", "O")

# Expected identities from UniProt O49686, the canonical 191-aa PYR1 sequence
# established in this project. Hard-coded so a shifted map cannot pass silently.
# Residue 2 is deliberately NOT listed: it is absent in S1/S2 and is ALA (not the
# native PRO) in S4, which is exactly the discrepancy that motivates this file.
EXPECT = {85: "SER", 86: "GLY", 87: "LEU", 88: "PRO", 89: "ALA",
          115: "HIS", 116: "ARG", 117: "LEU",
          59: "LYS", 108: "PHE", 79: "ARG", 94: "GLU"}
# tleap renames histidines by protonation state
ALIAS = {"HIE": "HIS", "HID": "HIS", "HIP": "HIS", "CYX": "CYS"}

# analysis unit -> (system directory, chain holding the PYR1 protomer)
#
# S3 contributes TWO units from the SAME trajectory. They are separate rows, not a
# single "S3", because the two protomers are in different conformational states
# (README 28g) and averaging them would erase the only apo-closed data we have.
UNITS = {
    "S1_apo_open":      ("S1_apo_open",    "A"),
    "S2_holo_closed":   ("S2_holo_closed", "A"),
    "S3_dimer_openA":   ("S3_apo_dimer",   "A"),
    "S3_dimer_closedB": ("S3_apo_dimer",   "B"),
    "S4_ternary":       ("S4_ternary",     "A"),
}
REFS = {"open": os.path.join(DATA, "pyr1_open_A.pdb"),
        "closed": os.path.join(DATA, "pyr1_closed_A.pdb")}

p = PDBParser(QUIET=True)


def residues(path, chain=None):
    """
    Ordered [(seq, resnum, resname)] of standard amino acids with a full backbone.

    `seq` counts across the WHOLE protein in file order, not from 1 within the
    requested chain, because that is what tleap does: in S3 the second protomer
    starts at prmtop residue 184, and restarting the count per chain would put
    every mask on the wrong protomer while looking perfectly sensible.
    """
    st = p.get_structure("x", path)[0]
    out, i = [], 0
    for ch in st:
        for r in ch:
            if r.id[0] != " " or not is_aa(r, standard=True):
                continue
            if not all(a in r for a in BB):
                continue
            i += 1
            if chain is None or ch.id == chain:
                out.append((i, r.id[1], ALIAS.get(r.get_resname(), r.get_resname())))
    return out


def prmtop_labels(path):
    """RESIDUE_LABEL from the prmtop -- the file cpptraj actually reads."""
    txt = open(path).read()
    m = re.search(r"%FLAG RESIDUE_LABEL\s*\n%FORMAT\((\d+)a(\d+)\)\s*\n(.*?)(?=%FLAG|\Z)",
                  txt, re.S)
    assert m, f"no RESIDUE_LABEL block in {path}"
    w = int(m.group(2))
    labels = []
    for line in m.group(3).rstrip("\n").split("\n"):
        labels += [line[i:i + w].strip() for i in range(0, len(line.rstrip()), w)]
    return [x for x in labels if x]


def build_map(reslist, label, common):
    """
    sequential (1-based, in prmtop order) -> native, verified by identity.

    `common` restricts the CORE (superposition) mask to residues present in EVERY
    unit and in the references. Without it S4's core mask spans 162 residues
    against the references' 161 -- S4 keeps 3QN1's residue 2 -- and cpptraj sets
    the rms up ANYWAY, silently returning ~84 A instead of refusing.
    """
    seq2nat = {i: rn for i, rn, _ in reslist}
    nat2seq = {rn: i for i, rn, _ in reslist}
    name_by_nat = {rn: nm for _, rn, nm in reslist}
    for nat, want in sorted(EXPECT.items()):
        assert nat in nat2seq, f"{label}: native residue {nat} absent"
        got = name_by_nat[nat]
        assert got == want, (
            f"{label}: native {nat} is {got}, expected {want} -- the numbering "
            f"map is WRONG, or this is not canonical PYR1")

    def span(nats, name):
        seqs = [nat2seq[n] for n in nats]
        assert seqs == list(range(seqs[0], seqs[-1] + 1)), (
            f"{label}: {name} not contiguous in sequential numbering: {seqs}")
        return f"{seqs[0]}-{seqs[-1]}"

    masks = {"gate": span(GATE, "gate"), "latch": span(LATCH, "latch"),
             "lb7a5": span(LB7A5, "Lb7a5")}

    mobile = set(GATE) | set(LATCH) | set(LB7A5)
    core_nat = [rn for _, rn, _ in reslist if rn not in mobile and rn in common]
    core_seq = sorted(nat2seq[n] for n in core_nat)
    ranges, start, prev = [], core_seq[0], core_seq[0]
    for v in core_seq[1:]:
        if v == prev + 1:
            prev = v
            continue
        ranges.append((start, prev))
        start = prev = v
    ranges.append((start, prev))
    core_mask = ",".join(f"{a}-{b}" if a != b else str(a) for a, b in ranges)
    seqs = [i for i, _, _ in reslist]
    return dict(n_residues=len(reslist), first_seq=min(seqs), last_seq=max(seqs),
                seq_to_native=seq2nat, native_to_seq=nat2seq,
                masks_sequential=masks, core_mask_sequential=core_mask,
                core_native=core_nat)


print("=" * 78)
print("PER-UNIT RESIDUE MAPS")
print("=" * 78)

ref_res = {k: residues(v) for k, v in REFS.items()}
assert [(n, m) for _, n, m in ref_res["open"]] == [(n, m) for _, n, m in ref_res["closed"]], (
    "the open and closed references disagree on numbering or identity; script 30 "
    "guarantees they match -- re-run 30_prepare_open_closed.py")
canonical = ref_res["open"]

# ---- one common core for every unit, computed BEFORE any map is built ----
# The superposition only means anything if every unit is fitted on the SAME
# residues. S3 chain B starts at native 2, so it cannot supply residue 1 and the
# 161-residue core used for the S1/S2/S4-only run is no longer available to all.
# Intersect rather than special-case: the alternative is per-unit cores that drift
# apart and make the RMSDs quietly incomparable.
present = {lab: {rn for _, rn, _ in residues(os.path.join(MD, d, "protein.pdb"), ch)}
           for lab, (d, ch) in UNITS.items()
           if os.path.exists(os.path.join(MD, d, "protein.pdb"))}
common = {rn for _, rn, _ in canonical}
for lab, s in present.items():
    common &= s
dropped = ({rn for _, rn, _ in canonical} - common)
print(f"  common residue set: {len(common)} of {len(canonical)} reference residues")
if dropped:
    print(f"  dropped (absent from at least one unit): {sorted(dropped)}")

refmap = build_map(canonical, "references", common)
print(f"  references     : {len(canonical)} res, gate = :{refmap['masks_sequential']['gate']}")

maps = {"_reference": refmap}
for lab, (d, ch) in UNITS.items():
    path = os.path.join(MD, d, "protein.pdb")
    if not os.path.exists(path):
        print(f"  {lab:<17}: protein.pdb absent -- skipped")
        continue
    rl = residues(path, ch)
    m = build_map(rl, lab, common)
    m["dir"], m["chain"] = d, ch
    maps[lab] = m
    same = "same as references" if m["masks_sequential"] == refmap["masks_sequential"] \
        else "DIFFERS from references"
    print(f"  {lab:<17}: {m['n_residues']} res (prmtop {m['first_seq']}-{m['last_seq']}), "
          f"gate = :{m['masks_sequential']['gate']}, "
          f"latch = :{m['masks_sequential']['latch']}   <- {same}")

# ---- re-verify against the TOPOLOGY, not just the pdb it was built from ----
# protein.pdb is an input to tleap; the prmtop is what cpptraj reads. They are
# supposed to agree residue-for-residue, so check it instead of assuming it.
print("\n  Landmark re-check against each prmtop's RESIDUE_LABEL:")
for lab, m in maps.items():
    if lab == "_reference":
        continue
    top = os.path.join(MD, m["dir"], "system.prmtop")
    if not os.path.exists(top):
        print(f"    {lab:<17} prmtop absent -- NOT verified")
        continue
    labels = prmtop_labels(top)
    for nat, want in sorted(EXPECT.items()):
        seq = m["native_to_seq"][nat]
        got = ALIAS.get(labels[seq - 1], labels[seq - 1])
        assert got == want, (
            f"{lab}: prmtop residue {seq} is {got}, but the map says native {nat} "
            f"= {want}. protein.pdb and system.prmtop disagree -- every mask "
            f"derived here would address the wrong residue.")
    print(f"    {lab:<17} {len(EXPECT)}/{len(EXPECT)} landmarks agree with {os.path.basename(top)}")

# the references must share numbering with whatever they are compared against
for lab, m in maps.items():
    if lab == "_reference":
        continue
    if m["masks_sequential"] != refmap["masks_sequential"]:
        print(f"\n  !! {lab} does not share the references' numbering. cpptraj masks "
              f"for\n     the TRAJECTORY and for the REFERENCE must therefore "
              f"differ -- 58b handles\n     this with an explicit refmask.")

# The superposition only means anything if the trajectory core and the reference
# core contain the SAME residues in the same order. Assert it rather than trust it.
for lab, m in maps.items():
    if lab == "_reference":
        continue
    assert m["core_native"] == refmap["core_native"], (
        f"{lab}: core residue list differs from the references "
        f"({len(m['core_native'])} vs {len(refmap['core_native'])}) -- cpptraj "
        f"would set the fit up anyway and return garbage")
print(f"\n  core superposition set: {len(refmap['core_native'])} residues, "
      f"identical in every unit (asserted)")

print("\n  Loop identities, verified per unit:")
for lab, m in maps.items():
    rl = canonical if lab == "_reference" else residues(
        os.path.join(MD, m["dir"], "protein.pdb"), m["chain"])
    nm = {rn: n for _, rn, n in rl}
    g = " ".join(nm[n] for n in GATE)
    l = " ".join(nm[n] for n in LATCH)
    print(f"    {lab:<17} gate={g}  latch={l}")

# ---- confirm open really is open ----
print()
print("=" * 78)
print("sanity: does the open reference actually differ from the closed one?")
print("=" * 78)
core_nat = refmap["core_native"]


def coords(path, nats, atoms=BB):
    st = p.get_structure("x", path)[0]
    want = {(n, a) for n in nats for a in atoms}
    got = {}
    for ch in st:
        for r in ch:
            if r.id[1] in nats and r.id[0] == " ":
                for a in r:
                    if a.get_id() in atoms:
                        got[(r.id[1], a.get_id())] = a.coord
    missing = want - set(got)
    assert not missing, f"missing atoms in {path}: {sorted(missing)[:5]}"
    return np.array([got[k] for k in sorted(want)])


def fitted_rmsd(pa, pb, nats):
    """superpose on the core, then RMSD of the loop -- no realignment on the loop"""
    ca, cb = coords(pa, core_nat), coords(pb, core_nat)
    mu_a, mu_b = ca.mean(0), cb.mean(0)
    A, B = ca - mu_a, cb - mu_b
    V, S, W = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(V @ W))
    R = V @ np.diag([1, 1, d]) @ W
    la = (coords(pa, nats) - mu_a) @ R
    lb = coords(pb, nats) - mu_b
    return float(np.sqrt(((la - lb) ** 2).sum() / len(la)))


g = fitted_rmsd(REFS["open"], REFS["closed"], GATE)
l = fitted_rmsd(REFS["open"], REFS["closed"], LATCH)
c = fitted_rmsd(REFS["open"], REFS["closed"], core_nat)
print(f"  gate  backbone open vs closed : {g:5.2f} A   (~5 A expected; script 24")
print(f"  latch backbone open vs closed : {l:5.2f} A    said 5.23/5.34 on the")
print(f"  core  backbone open vs closed : {c:5.2f} A    untrimmed pair)")
assert g > 3.0 and l > 3.0, (
    "the references differ by <3 A at the gate/latch -- not an open/closed pair, "
    "and the whole two-reference projection is meaningless")
assert c < 1.5, f"core differs by {c:.2f} A; it is supposed to be the rigid part"
print("  -> the 5 A stroke is real and the core is rigid: the projection is valid")

out = dict(
    note="sequential = 1-based residue index in the tleap prmtop, counted across "
         "the whole protein; native = PYR1/UniProt O49686 numbering. MAPS ARE PER "
         "UNIT: S4 carries 3QN1's residue 2, and S3 contributes two protomers from "
         "one trajectory (chain B starts at prmtop 184).",
    gate_native=GATE, latch_native=LATCH, lb7a5_native=LB7A5,
    common_core_native=sorted(common),
    systems=maps,
    ref_open_vs_closed={"gate_bb_rmsd": g, "latch_bb_rmsd": l, "core_bb_rmsd": c},
)
mp = os.path.join(OUT, "residue_map.json")
json.dump(out, open(mp, "w"), indent=1)
print(f"\nwrote {mp}")
print("\nNEXT: scripts/58b_loop_dynamics_run.sh")
