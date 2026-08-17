#!/usr/bin/env python
"""
lib_resnumber.py -- THE residue-identity and numbering layer for every structure this
project hands to MD.

WHY THIS EXISTS
---------------
Residue NUMBERS lie in this project, repeatedly and silently:

  * 3QN1 chain A models 1-181 with 69/70 unresolved; 3K3K models 1-183 with no gaps.
  * script 30 intersected the two and dropped residue 2 as well, so S1/S2 carry a
    THIRD gap that exists in neither crystal.
  * pdb4amber/tleap renumber 1..N, so native 85 is sequential 82 in S1/S2 but 83 in
    S4 -- and cpptraj will happily fit a 162-residue mask onto a 161-residue
    reference and return a gate RMSD of 84 A instead of erroring.
  * 3QN1's residue 2 is ALA, not the wild-type PRO. That is a MUTATION wearing the
    costume of a numbering artefact, and it propagated into AF3, BindCraft and Boltz.
  * AF3 chimeras start `GAMASEL...` so their numbering is offset by 2 from `MPSEL...`.

Every one of those cost real time. The rule this module enforces is the one from
`feedback-verify-residue-identity`: never trust a residue number, assert the residue
IS what you expect. Do it at BUILD time, once, and emit a map that downstream
analysis reads instead of re-deriving.

WHAT IT GUARANTEES
------------------
`prepare_chain()` returns a chain whose residues have been altloc-resolved and
identity-checked, and `verify_build()` refuses to pass a structure that does not have
the exact expected sequence with zero broken peptide bonds.

USAGE
    from lib_resnumber import (PYR1_SEQ, prepare_chain, chain_sequence,
                               peptide_breaks, verify_build, write_residue_map)
"""
from __future__ import annotations

import hashlib
import json
import os
from collections import OrderedDict

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

# --------------------------------------------------------------------------------
# canonical sequence
# --------------------------------------------------------------------------------
# UniProt O49686 (PYR1_ARATH) residues 1-191, verbatim. This is a VENDORED copy of
# HAB1_shrinkinator/BoltzProt-1/esmfold2_repredict/pyr1_sequence.py, which remains the
# single source of truth; the md5 guard below makes any drift between the two loud
# instead of silent. It was also re-derived independently here from the
# `_entity_poly.pdbx_seq_one_letter_code_can` records of BOTH 3K3K (211 aa = 20 aa
# His-tag/thrombin + 191) and 3QN1 (193 aa = 2 aa `GA` remnant + 191). The two
# crystals' deposited sequences agree with each other and with UniProt at all 191
# positions EXCEPT residue 2, where 3QN1 has the A2 substitution.
PYR1_SEQ = (
    "MPSELTPEERSELKNSIAEFHTYQLDPGSCSSLHAQRIHAPPELVWSIVRRFDKPQTYKH"   # 1-60
    "FIKSCSVEQNFEMRVGCTRDVIVISGLPANTSTERLDILDDERRVTGFSIIGGEHRLTNY"   # 61-120
    "KSVTTVHRFEKENRIWTVVLESYVVDMPEGNSEDDTRMFADTVVKLNLQKLATVAEAMAR"   # 121-180
    "NSGDGSGSQVT"                                                     # 181-191
)
PYR1_MD5 = "5ff74f5d34"
assert len(PYR1_SEQ) == 191
assert hashlib.md5(PYR1_SEQ.encode()).hexdigest()[:10] == PYR1_MD5, \
    "PYR1_SEQ drifted from UniProt O49686"
assert PYR1_SEQ[1] == "P", "residue 2 must be wild-type PRO, not 3QN1's A2"
assert PYR1_SEQ[68:70] == "QN", "Q69/N70 gap is back"

# Residues whose identity is asserted on EVERY new structure/pipeline combination.
# Gate, latch, and the pocket/salt-bridge positions this project actually reasons about.
PYR1_LANDMARKS = {
    85: "SER", 86: "GLY", 87: "LEU", 88: "PRO", 89: "ALA",   # gate
    115: "HIS", 116: "ARG", 117: "LEU",                      # latch
    59: "LYS", 108: "PHE", 79: "ARG", 94: "GLU",             # pocket / salt bridge
    2: "PRO",                                                # the P2A tripwire
}

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
    # protonation/derivative aliases that mean the same residue for identity purposes
    "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "CYM": "C", "ASH": "D",
    "GLH": "E", "LYN": "K", "MSE": "M",
}
# Built explicitly, NOT by inverting THREE2ONE -- that dict is many-to-one (HID/HIE/HIP
# all map to H), so inverting it would silently pick whichever alias came last.
ONE2THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}

PEPTIDE_BOND_MAX = 1.5   # A; a real C-N bond is ~1.33


# --------------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------------
def load_model(path, model_index=0):
    """Parse a .cif or .pdb and return one model. Format is chosen by extension."""
    ext = os.path.splitext(path)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
    return parser.get_structure("s", path)[model_index]


def select_altloc(chain, keep="A"):
    """Resolve alternate conformations EXPLICITLY rather than taking Biopython's default.

    3K3K carries 0.5-occupancy alternates at >=8 positions including R116 -- a LATCH
    residue -- while 3QN1 carries none. Letting each file's default win means the open
    and closed systems are built from silently different choices at exactly the
    positions this project measures. So the choice is named, applied, and reported.

    Returns the number of residues that had a choice to make.
    """
    n = 0
    for res in chain:
        for atom in res:
            if atom.is_disordered():
                ids = atom.disordered_get_id_list()
                atom.disordered_select(keep if keep in ids else ids[0])
                n += 1
    return n


def chain_residues(chain, keep_hetero=False):
    """Ordered {resseq: Residue} for the standard amino acids of one chain."""
    out = OrderedDict()
    for res in chain:
        het, num, icode = res.id
        if icode.strip():
            raise ValueError(f"insertion code {icode!r} at residue {num}; not handled")
        if het.strip() and not keep_hetero:
            continue
        out[num] = res
    return out


# --------------------------------------------------------------------------------
# assertions
# --------------------------------------------------------------------------------
def assert_identity(resmap, expect, label):
    """Assert residue NUMBER -> residue NAME for every landmark. Raises on any mismatch.

    `expect` maps residue number to a three-letter name. Missing residues are
    mismatches too -- a landmark that is not modelled is exactly as dangerous as one
    that is the wrong amino acid.
    """
    bad = []
    for num, want in sorted(expect.items()):
        got = resmap[num].get_resname().upper() if num in resmap else "MISSING"
        if THREE2ONE.get(got) != THREE2ONE.get(want.upper()):
            bad.append((num, want, got))
    if bad:
        lines = "\n".join(f"    residue {n}: expected {w}, found {g}" for n, w, g in bad)
        raise AssertionError(f"{label}: residue identity check FAILED\n{lines}")
    return len(expect)


def chain_sequence(resmap):
    """(sequence string, residue numbers) in ascending numeric order.

    Note this is the sequence of what is MODELLED. Gaps do not appear as dashes --
    use `peptide_breaks` to find them.
    """
    nums = sorted(resmap)
    seq = "".join(THREE2ONE.get(resmap[n].get_resname().upper(), "X") for n in nums)
    return seq, nums


def peptide_breaks(resmap, tol=PEPTIDE_BOND_MAX):
    """Every consecutive pair whose C-N distance exceeds `tol`.

    Deliberately measures the BOND, not the numbering: a structure can have
    consecutive residue numbers and still be broken, and can have a numbering jump
    that is genuinely bonded. Both have happened here.
    """
    nums = sorted(resmap)
    out = []
    for a, b in zip(nums, nums[1:]):
        ra, rb = resmap[a], resmap[b]
        if "C" not in ra or "N" not in rb:
            out.append((a, ra.get_resname(), b, rb.get_resname(), float("nan")))
            continue
        d = float(np.linalg.norm(ra["C"].coord - rb["N"].coord))
        if d > tol:
            out.append((a, ra.get_resname(), b, rb.get_resname(), round(d, 2)))
    return out


def verify_build(path, chain_id, expect_seq=PYR1_SEQ, first_resnum=1,
                 landmarks=PYR1_LANDMARKS, label=None):
    """The gate every built structure must pass before it is allowed near tleap.

    Checks, in order:
      1. the modelled residue numbers are exactly first_resnum .. first_resnum+len-1
      2. the one-letter sequence equals `expect_seq` EXACTLY
      3. every landmark residue is the expected amino acid
      4. there are ZERO peptide-bond breaks

    Raises AssertionError with a specific message on the first failure. Returns a dict
    suitable for writing beside the structure as provenance.
    """
    label = label or f"{os.path.basename(path)}:{chain_id}"
    model = load_model(path)
    chain = model[chain_id]
    select_altloc(chain)
    resmap = chain_residues(chain)
    seq, nums = chain_sequence(resmap)

    want_nums = list(range(first_resnum, first_resnum + len(expect_seq)))
    if nums != want_nums:
        missing = sorted(set(want_nums) - set(nums))
        extra = sorted(set(nums) - set(want_nums))
        raise AssertionError(
            f"{label}: residue numbering is not {want_nums[0]}..{want_nums[-1]}\n"
            f"    modelled {len(nums)} residues {nums[0]}..{nums[-1]}\n"
            f"    missing: {missing[:20]}\n    unexpected: {extra[:20]}")

    if seq != expect_seq:
        diffs = [(i + first_resnum, e, g)
                 for i, (e, g) in enumerate(zip(expect_seq, seq)) if e != g]
        raise AssertionError(
            f"{label}: sequence does not match the canonical target\n" +
            "\n".join(f"    residue {n}: expected {e}, built {g}" for n, e, g in diffs[:20]))

    assert_identity(resmap, landmarks, label)

    breaks = peptide_breaks(resmap)
    if breaks:
        raise AssertionError(
            f"{label}: {len(breaks)} broken peptide bond(s)\n" +
            "\n".join(f"    {a}{an} -> {b}{bn}  C-N = {d} A" for a, an, b, bn, d in breaks))

    return {
        "label": label, "file": os.path.abspath(path), "chain": chain_id,
        "n_residues": len(nums), "first": nums[0], "last": nums[-1],
        "sequence_md5": hashlib.md5(seq.encode()).hexdigest()[:10],
        "landmarks_checked": sorted(landmarks),
        "peptide_breaks": 0,
    }


# --------------------------------------------------------------------------------
# provenance: facts about source structures, checked rather than remembered
# --------------------------------------------------------------------------------
PROVENANCE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "data", "structure_provenance.json")


def assert_provenance(cif_path, chain_id, path=PROVENANCE, ligand_cutoff=4.5):
    """Check a source structure against `data/structure_provenance.json` before using it.

    This exists because a system was once built as "S3_apo_dimer" from 3K3K, whose
    chain B is closed and ABA-bound. The information was on record in several places;
    the SYSTEM NAME said otherwise, and the name won. Names summarise, and a summary
    that is wrong outranks the records it was derived from.

    So the facts are asserted here, at build time, from the structure itself. What is
    checked: modelled residue range, internal gaps, broken peptide bonds, and -- the
    one that was missed -- WHICH LIGANDS THE CHAIN ACTUALLY CARRIES.

    If this raises, do not loosen it. Either the structure is not the one expected, or
    the JSON is wrong and should be corrected with a measurement.
    """
    with open(path) as fh:
        prov = json.load(fh)
    key = os.path.splitext(os.path.basename(cif_path))[0].upper()
    if key not in prov:
        raise AssertionError(
            f"{key} has no entry in {os.path.basename(path)}. Every structure this "
            "project builds from must be characterised there first -- measure it, "
            "record it, then build.")
    want = prov[key]["chains"].get(chain_id)
    if want is None:
        raise AssertionError(f"{key} chain {chain_id} has no provenance entry")

    model = load_model(cif_path)
    chain = model[chain_id]
    select_altloc(chain)
    resmap = chain_residues(chain)
    nums = sorted(resmap)
    problems = []

    lo, hi = want["modelled"]
    if [nums[0], nums[-1]] != [lo, hi]:
        problems.append(f"modelled range is {nums[0]}-{nums[-1]}, expected {lo}-{hi}")

    gaps = [i for i in range(nums[0], nums[-1] + 1) if i not in resmap]
    flat = []
    for g in want["internal_gaps"]:
        flat.extend(range(g[0], g[1] + 1) if isinstance(g, list) else [g])
    if gaps != sorted(flat):
        problems.append(f"internal gaps are {gaps}, expected {sorted(flat)}")

    n_break = len(peptide_breaks(resmap))
    if n_break != want["peptide_breaks"]:
        problems.append(f"{n_break} broken peptide bonds, expected {want['peptide_breaks']}")

    # the check that would have caught the S3 mistake: ligands belong to a CHAIN,
    # not to a file, and a dimer's two protomers need not be occupied alike
    lc = []
    for r in chain:
        if r.id[0].strip() and r.get_resname() not in ("HOH", "WAT"):
            lc.append(r.get_resname().upper())
    for other in model:
        if other.id == chain_id:
            continue
        for r in other:
            if not r.id[0].strip() or r.get_resname() in ("HOH", "WAT"):
                continue
            hcs = np.array([a.coord for a in r if a.element != "H"])
            for num in nums:
                pc = np.array([a.coord for a in resmap[num] if a.element != "H"])
                if np.min(np.linalg.norm(pc[:, None, :] - hcs[None, :, :], axis=-1)) \
                        <= ligand_cutoff:
                    lc.append(r.get_resname().upper())
                    break
    if sorted(lc) != sorted(x.upper() for x in want["ligands"]):
        problems.append(
            f"ligands in contact with this chain are {sorted(lc)}, expected "
            f"{sorted(want['ligands'])} -- an unexpected ligand means the chain is not "
            "in the occupancy state the pipeline assumes")

    for num, name in want.get("sequence_variants", {}).items():
        got = resmap[int(num)].get_resname().upper() if int(num) in resmap else "MISSING"
        if got != name.upper():
            problems.append(f"residue {num} is {got}, expected the known variant {name}")

    if problems:
        raise AssertionError(
            f"{key} chain {chain_id} does not match its provenance record:\n" +
            "\n".join(f"    {p}" for p in problems))
    return want


# --------------------------------------------------------------------------------
# the map downstream analysis reads
# --------------------------------------------------------------------------------
def write_residue_map(path, out_json, chains, landmarks=PYR1_LANDMARKS):
    """Emit native <-> sequential residue maps at BUILD time.

    tleap renumbers every chain into one 1..N sequence, so `:85-89` in a cpptraj mask
    is NOT the gate. Downstream scripts must read this file rather than re-deriving
    the offset, and must never assume the offset is shared between systems -- it is 0
    for S3, 1 for S1/S2 and 3 for others, and it changed again once the gaps were
    filled.

    `chains` is an ordered list of chain ids, in the order tleap will concatenate them.
    """
    model = load_model(path)
    seq_index = 0
    native_to_seq, seq_to_native, per_chain = {}, {}, {}
    for cid in chains:
        chain = model[cid]
        select_altloc(chain)
        resmap = chain_residues(chain)
        entries = {}
        for num in sorted(resmap):
            seq_index += 1
            key = f"{cid}:{num}"
            native_to_seq[key] = seq_index
            seq_to_native[str(seq_index)] = key
            entries[str(num)] = {"seq": seq_index,
                                 "name": resmap[num].get_resname().upper()}
        per_chain[cid] = {"n": len(resmap), "first_seq": seq_index - len(resmap) + 1,
                          "last_seq": seq_index, "residues": entries}

    payload = {
        "source": os.path.abspath(path),
        "chain_order": list(chains),
        "total_residues": seq_index,
        "native_to_seq": native_to_seq,
        "seq_to_native": seq_to_native,
        "chains": per_chain,
        "landmarks_seq": {str(n): native_to_seq.get(f"{chains[0]}:{n}")
                          for n in sorted(landmarks)},
    }
    with open(out_json, "w") as fh:
        json.dump(payload, fh, indent=2)
    return payload
