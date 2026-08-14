#!/usr/bin/env python
"""
58_loop_dynamics_prep.py -- establish and VERIFY the residue-numbering map for the
gate/latch loop-dynamics analysis, PER SYSTEM, then hand the masks to 58b.

WHY THIS SCRIPT EXISTS SEPARATELY FROM THE ANALYSIS
---------------------------------------------------
tleap renumbers residues sequentially 1..N in the prmtop. The gate (85-89) and
latch (115-117) are quoted in NATIVE PYR1 numbering everywhere in this project,
and the crystal files have gaps, so `:85-89` in a cpptraj mask does NOT select the
gate. Getting this wrong is silent: cpptraj returns RMSD values for the wrong loop
and every downstream conclusion is confidently wrong.

AND THE MAP IS NOT THE SAME FOR EVERY SYSTEM. That is the trap this file was
rewritten to close:

    S1_apo_open / S2_holo_closed   178 residues, gate = sequential 82-86
    S4_ternary                     179 residues, gate = sequential 83-87

S1 and S2 come from script 30, which took the INTERSECTION of 3K3K and 3QN1 and so
dropped residue 2 (disordered in 3K3K). S4 was built straight from 3QN1, which HAS
residue 2 -- as ALA, the P2A that README 23 records 3QN1 handing to every pipeline
in this project. One extra residue near the N-terminus shifts everything after it
by one, so reusing S1's masks on S4 would measure a loop one residue out of
register and never complain.

So the map is rebuilt from each system's own protein.pdb and verified by residue
IDENTITY, never assumed or shared.

WHAT IT CHECKS
--------------
  1. per system: sequential index -> native resnum, built by ORDER within the PYR1
     chain, then verified by RESNAME at a dozen known positions
  2. the gate, latch and Lbeta7-alpha5 loops carry their expected sequences
     (gate = SER GLY LEU PRO ALA, latch = HIS ARG LEU) in every system
  3. the open/closed references really are an open/closed pair, by recovering a
     ~5 A gate/latch stroke. A MAGNITUDE check only: the 5.23/5.34 A figures in
     script 30's header come from script 24 on the UNTRIMMED structures, and this
     script uses a third core again (it also excludes Lb7a5), so exact agreement
     is not expected and would be suspicious if claimed.

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

# system -> chain holding PYR1 in that system's protein.pdb
SYSTEMS = {"S1_apo_open": "A", "S2_holo_closed": "A", "S4_ternary": "A"}
REFS = {"open": os.path.join(DATA, "pyr1_open_A.pdb"),
        "closed": os.path.join(DATA, "pyr1_closed_A.pdb")}

p = PDBParser(QUIET=True)


def residues(path, chain=None):
    """ordered [(resnum, resname)] of standard amino acids with a full backbone"""
    st = p.get_structure("x", path)[0]
    out = []
    for ch in st:
        if chain is not None and ch.id != chain:
            continue
        for r in ch:
            if r.id[0] != " " or not is_aa(r, standard=True):
                continue
            if not all(a in r for a in BB):
                continue
            out.append((r.id[1], ALIAS.get(r.get_resname(), r.get_resname())))
    return out


def build_map(reslist, label, common=None):
    """
    sequential (1-based, in prmtop order) -> native, verified by identity.

    `common` restricts the CORE (superposition) mask to residues that also exist
    in the reference structures. Without it S4's core mask spans 162 residues
    against the references' 161 -- S4 keeps 3QN1's residue 2 -- and cpptraj sets
    the rms up ANYWAY, silently returning ~84 A instead of refusing. The fit never
    happens and every loop RMSD downstream is meaningless.
    """
    seq2nat = {i: rn for i, (rn, _) in enumerate(reslist, start=1)}
    nat2seq = {rn: i for i, rn in seq2nat.items()}
    name_by_nat = {rn: nm for rn, nm in reslist}
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
    core_nat = [rn for rn, _ in reslist
                if rn not in mobile and (common is None or rn in common)]
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
    return dict(n_residues=len(reslist), seq_to_native=seq2nat,
                native_to_seq=nat2seq, masks_sequential=masks,
                core_mask_sequential=core_mask, core_native=core_nat)


print("=" * 78)
print("PER-SYSTEM RESIDUE MAPS")
print("=" * 78)

ref_res = {k: residues(v) for k, v in REFS.items()}
assert ref_res["open"] == ref_res["closed"], (
    "the open and closed references disagree on numbering or identity; script 30 "
    "guarantees they match -- re-run 30_prepare_open_closed.py")
canonical = ref_res["open"]
refmap = build_map(canonical, "references")
print(f"  references     : {len(canonical)} res, gate = :{refmap['masks_sequential']['gate']}")

maps = {"_reference": refmap}
for s, chain in SYSTEMS.items():
    path = os.path.join(MD, s, "protein.pdb")
    if not os.path.exists(path):
        print(f"  {s:<15}: protein.pdb absent -- skipped")
        continue
    rl = residues(path, chain)
    m = build_map(rl, s, common={rn for rn, _ in canonical})
    maps[s] = m
    same = "same as references" if m["masks_sequential"] == refmap["masks_sequential"] \
        else "DIFFERS from references"
    print(f"  {s:<15}: {m['n_residues']} res, gate = :{m['masks_sequential']['gate']}, "
          f"latch = :{m['masks_sequential']['latch']}   <- {same}")

# the references must share numbering with whatever they are compared against
for s in ("S1_apo_open", "S2_holo_closed", "S4_ternary"):
    if s not in maps:
        continue
    if maps[s]["masks_sequential"] != refmap["masks_sequential"]:
        print(f"\n  !! {s} does not share the references' numbering. cpptraj masks "
              f"for\n     the TRAJECTORY and for the REFERENCE must therefore "
              f"differ -- 58b handles\n     this with an explicit refmask.")

# The superposition only means anything if the trajectory core and the reference
# core contain the SAME residues in the same order. Assert it rather than trust it.
for s2, m in maps.items():
    if s2 == "_reference":
        continue
    assert m["core_native"] == refmap["core_native"], (
        f"{s2}: core residue list differs from the references "
        f"({len(m['core_native'])} vs {len(refmap['core_native'])}) -- cpptraj "
        f"would set the fit up anyway and return garbage")
print(f"\n  core superposition set: {len(refmap['core_native'])} residues, "
      f"identical in every system (asserted)")

print("\n  Loop identities, verified per system:")
for s, m in maps.items():
    nm = {rn: n for rn, n in (canonical if s == "_reference"
                              else residues(os.path.join(MD, s, "protein.pdb"),
                                            SYSTEMS[s]))}
    g = " ".join(nm[n] for n in GATE)
    l = " ".join(nm[n] for n in LATCH)
    print(f"    {s:<15} gate={g}  latch={l}")

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
    note="sequential = 1-based residue index in the tleap prmtop (PYR1 chain "
         "only); native = PYR1/UniProt O49686 numbering. MAPS ARE PER SYSTEM: "
         "S4 carries 3QN1's residue 2 and is shifted by one relative to S1/S2.",
    gate_native=GATE, latch_native=LATCH, lb7a5_native=LB7A5,
    systems=maps,
    ref_open_vs_closed={"gate_bb_rmsd": g, "latch_bb_rmsd": l, "core_bb_rmsd": c},
)
mp = os.path.join(OUT, "residue_map.json")
json.dump(out, open(mp, "w"), indent=1)
print(f"\nwrote {mp}")
print("\nNEXT: scripts/58b_loop_dynamics_run.sh")
