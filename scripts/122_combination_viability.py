#!/usr/bin/env python
r"""
122_combination_viability.py -- can Rosetta tell a REAL sensor from a random
member of the SAME library, using the protein alone?

THE QUESTION THIS ANSWERS (Jannis, 2026-08-25)
Every previous attempt asked structure-based tools to pick the right RESIDUE, and
every one failed (§31 admissibility, §32's K59R miss, §33 electrostatics, §23j
LigandMPNN/FastDesign, §37-§43 MM-GBSA). §63-§65 then showed round-1 data already
narrows the menu well -- positions especially. So the tools no longer have to
generate candidates. They only have to RANK combinations inside a menu the wet lab
has already chosen.

That is a different and much easier job, and crucially it can be done POSE-FREE.
The circularity that killed everything pose-based (§49b: you need the ligand
placed, but the pocket you would place it in does not exist yet, and docking into
a hypothetical pocket measured 8.33 A error, §46b) does not apply if the score
never uses the ligand. A working sensor must at minimum FOLD and PACK. Asking
"is this combination of 5-8 mutations structurally viable" needs no ligand at all.

THE TEST
Score three sets on protein-only ref2015 after repacking the mutated shell:

  REAL     the 78 known round-2 coumarin sensors (sd07)
  LIBRARY  random members of Tian's own 138,240-member coumarin library
  WILD     random combinations at the same 11 positions drawn from the FULL
           DSM-Hao menu, i.e. residues the library did NOT select

Two nulls, deliberately. WILD asks whether Rosetta has any signal at all here.
LIBRARY asks the question that matters: once the wet lab has already picked the
menu, does Rosetta add anything on top?

⚠ RECALL ONLY, AND THE NULL IS CONTAMINATED. The 78 are hits that were FOUND, not
optima ([[feedback_hits_are_not_optima]]), and a random member of Tian's library
may well be a perfectly good sensor nobody screened. So LIBRARY is not a clean
negative set and a null result there is weak evidence. A POSITIVE result --
real sensors scoring clearly better than random library members -- would be
strong, because the contamination runs against it.

Usage: 122_combination_viability.py --chunk i --nchunks n
"""
import argparse
import json
import os
import random
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_xlsx import table                                    # noqa: E402
import lib_rosetta as LR                                       # noqa: E402

ROOT = os.path.dirname(HERE)
SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "viability")
PDB = os.path.join(ROOT, "data", "structures_191", "pyr1_closed_191.pdb")
N_RANDOM = 300          # per null set
SEED = 0
AA = "ACDEFGHIKLMNPQRSTVWY"


def sd07_sensors():
    hdr, recs = table(f"{SD}/pnas.2519924122.sd07(1).xlsx")
    cols = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    out = []
    for r in recs:
        subs = []
        for c in cols:
            v = (r.get(c) or "").strip().upper()
            if len(v) == 1 and v in AA and v != c[0]:
                subs.append((int(c[1:]), c[0], v))
        if subs:
            out.append({"id": (r.get("name") or "?").strip(),
                        "lig": (r.get("compound") or "?").strip(), "subs": subs})
    return out


def sd04_menu(name):
    _, recs = table(f"{SD}/pnas.2519924122.sd04.xlsx")
    m = defaultdict(set)
    for r in recs:
        if (r.get("Library") or "").strip() == name:
            m[int(float(r["Position"]))] = (
                m[int(float(r["Position"]))] |
                set(r["Amino Acids allowed for mutation"].strip()))
    return dict(m)


def wt_at(pose, num):
    """Wild-type one-letter code at a NATIVE residue number, read from the pose."""
    info = pose.pdb_info()
    for i in range(1, pose.total_residue() + 1):
        if info.number(i) == num:
            return pose.residue(i).name1(), i
    raise SystemExit(f"native residue {num} not found in {PDB}")


def build_sets(pose):
    cou = sd04_menu("Coumarin")
    dsm = sd04_menu("DSM-Hao")
    positions = sorted(cou)
    # assert the structure agrees with the sheet about what wild-type is
    for p in positions:
        one, _ = wt_at(pose, p)
        want = None
        _, recs = table(f"{SD}/pnas.2519924122.sd04.xlsx")
        for r in recs:
            if ((r.get("Library") or "").strip() == "Coumarin"
                    and int(float(r["Position"])) == p):
                want = r["WT"].strip()
        if want and one != want:
            raise SystemExit(f"position {p}: structure has {one}, sd04 says {want}")
    rng = random.Random(SEED)
    real = [{"set": "REAL", "id": s["id"], "subs": s["subs"]}
            for s in sd07_sensors()
            if all(n in cou for n, _w, _m in s["subs"])]
    n_sub = [len(r["subs"]) for r in real]
    lo, hi = min(n_sub), max(n_sub)

    def draw(menu, tag):
        out = []
        for i in range(N_RANDOM):
            k = rng.randint(lo, hi)
            ps = rng.sample(positions, min(k, len(positions)))
            subs = []
            for p in ps:
                one, _ = wt_at(pose, p)
                choices = sorted(menu[p] - {one})
                if choices:
                    subs.append((p, one, rng.choice(choices)))
            if subs:
                out.append({"set": tag, "id": f"{tag}_{i}", "subs": subs})
        return out

    # WILD uses the full DSM-Hao menu at the SAME positions, so the only thing
    # that differs from LIBRARY is which residues the wet lab selected
    return real + draw(cou, "LIBRARY") + draw(dsm, "WILD"), n_sub


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk", type=int, default=0)
    ap.add_argument("--nchunks", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    import pyrosetta
    from pyrosetta.rosetta.core.scoring import get_score_function
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector)
    pyrosetta.init("-mute all -ignore_unrecognized_res -ex1 -ex2aro")
    base = pyrosetta.pose_from_pdb(PDB)
    sfxn = get_score_function()

    variants, n_sub = build_sets(base)
    variants = [v for i, v in enumerate(variants) if i % a.nchunks == a.chunk]
    print(f"chunk {a.chunk}/{a.nchunks}: {len(variants)} variants "
          f"(sensor sizes {min(n_sub)}-{max(n_sub)} substitutions)")

    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    THREE = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
             "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
             "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
             "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
    wt_raw = sfxn(base)
    rows = []
    for v in variants:
        pose = base.clone()
        idx = []
        ok = True
        for num, wt, mu in v["subs"]:
            one, i = wt_at(pose, num)
            if one != wt:
                print(f"  !! {v['id']}: residue {num} is {one}, sheet says {wt}")
                ok = False
                break
            MutateResidue(i, THREE[mu]).apply(pose)
            idx.append(i)
        if not ok:
            continue
        sel = ResidueIndexSelector(",".join(str(i) for i in idx))
        nb = NeighborhoodResidueSelector(sel, 8.0, True)
        allowed = [i + 1 for i, b in enumerate(nb.apply(pose)) if b]
        tf, _packable = LR.restrict_packing(pose, allowed)  # returns (tf, verified)
        pk = PackRotamersMover(sfxn)
        pk.task_factory(tf)
        pk.apply(pose)
        s_mut = sfxn(pose)

        # ⚠ PAIRED CONTROL, and it is not optional. The first version scored the
        # mutant after repacking its shell against the RAW wild-type score, and
        # every variant came out 130-165 REU "better" than wild-type -- that is
        # the input structure's own unrelieved strain being repacked away, not an
        # effect of the mutations. It also biases toward variants with MORE
        # mutations, since they open a larger shell. Repacking the identical shell
        # on the unmutated pose removes both.
        ref = base.clone()
        tf2, _ = LR.restrict_packing(ref, allowed)
        pk2 = PackRotamersMover(sfxn)
        pk2.task_factory(tf2)
        pk2.apply(ref)
        s_wt = sfxn(ref)

        rows.append({"set": v["set"], "id": v["id"], "n_sub": len(v["subs"]),
                     "score": float(s_mut), "wt_same_shell": float(s_wt),
                     "ddG": float(s_mut - s_wt), "n_shell": len(allowed),
                     "subs": [f"{w}{n}{m}" for n, w, m in v["subs"]]})
        print(f"  {v['set']:<8}{v['id']:<14}{len(v['subs'])} subs  "
              f"shell {len(allowed):3d}  mut {s_mut:9.2f}  wt {s_wt:9.2f}  "
              f"ddG {s_mut-s_wt:+8.2f}")
    with open(os.path.join(OUT, f"chunk_{a.chunk:03d}.json"), "w") as fh:
        json.dump({"wt_raw": float(wt_raw), "rows": rows}, fh, indent=1)
    print(f"wrote {OUT}/chunk_{a.chunk:03d}.json  ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
