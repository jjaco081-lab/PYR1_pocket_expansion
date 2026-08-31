#!/usr/bin/env python
r"""
140_tractability_yaml.py -- Boltz-2 inputs for the matched tractability benchmark.

One protein (wild-type PYR1, the SAME construct the June 2025 sensor run used, so
p_bind values are comparable to §82's 637) against each ligand, with the affinity
head enabled. Both classes get the identical receptor -- see §139's header for why
that symmetry is non-negotiable.

The MSA is computed ONCE and reused: every run has the same protein, so paying the
server 362 times would be the dominant cost for no benefit.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "data", "tractability")
#: identical to the construct in the June 2025 sensor predictions (§80a consensus)
PYR1 = ("MPSELTPEERSELKNSIAEFHTYQLDPGSCSSLHAQRIHAPPELVWSIVRRFDKPQTYKHFIKSCSVEQNF"
        "EMRVGCTRDVIVISGLPANTSTERLDILDDERRVTGFSIIGGEHRLTNYKSVTTVHRFEKENRIWTVVLES"
        "YVVDMPEGNSEDDTRMFADTVVKLNLQKLATVAEAMARNSGDGSGSQVT")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "matched"
    msa = sys.argv[2] if len(sys.argv) > 2 else None
    recs = json.load(open(os.path.join(OUT, f"set_{which}.json")))
    d = os.path.join(OUT, f"yaml_{which}")
    os.makedirs(d, exist_ok=True)
    n = 0
    for i, r in enumerate(recs):
        # SINGLE-quoted YAML: no escape processing. A double-quoted scalar
        # treats a backslash as an escape, so any SMILES carrying cis/trans
        # stereochemistry (C/C=C\\C) either fails to parse or is silently
        # corrupted -- \\N is YAML's NEL character, which reached RDKit as a
        # control byte and returned None. That lost 14 of 362 runs, and it lost
        # them NON-RANDOMLY: every casualty had stereochemistry (README 92a).
        smi = r["smiles"].replace("'", "''")
        block = (f"      msa: {msa}\n" if msa else "")
        with open(os.path.join(d, f"lig{i:04d}.yaml"), "w") as fh:
            fh.write(
                "version: 1\n"
                "sequences:\n"
                "  - protein:\n"
                "      id: A\n"
                f"      sequence: {PYR1}\n"
                f"{block}"
                "  - ligand:\n"
                "      id: B\n"
                f"      smiles: '{smi}'\n"
                "properties:\n"
                "  - affinity:\n"
                "      binder: B\n")
        n += 1
    key = {f"lig{i:04d}": {"name": r["name"], "label": r["label"],
                           "heavy": r["heavy"], "clogp": r["clogp"]}
           for i, r in enumerate(recs)}
    json.dump(key, open(os.path.join(OUT, f"key_{which}.json"), "w"), indent=1)
    print(f"wrote {n} YAMLs to {d}")
    print(f"  {sum(1 for r in recs if r['label']==1)} hits, "
          f"{sum(1 for r in recs if r['label']==0)} non-hits")
    print(f"  protein length {len(PYR1)} aa; msa = {msa or 'server (slow)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
