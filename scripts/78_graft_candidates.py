#!/usr/bin/env python
"""
78_graft_candidates.py -- START/SRPBCC (helix-grip) relatives of PYR1 that carry a
real ligand, ranked as grafting candidates.

WHERE THE NUMBERS COME FROM
---------------------------
`05_foldseek_search.sh` searched PYR1 (3QN1 chain A) against **CATH50** - a
50%-sequence-identity-clustered representative set of CATH domains that INCLUDES
AlphaFold models. The superfamily is CATH 3.30.530.20, helix-grip / SRPBCC / START.

That mixture matters for this question: most of the fold's diversity in CATH50 is
predicted, not solved, and a predicted model has no ligand in it at all. So "how
many relatives are there" and "how many can show us a bound ligand" differ by a lot.

WHAT COUNTS AS A LIGAND
-----------------------
HETATM alone will not do: it also catches water, cryoprotectants and buffer. A
ligand here is a HETATM residue that is NOT in the solvent/additive list below and
has at least MIN_HEAVY heavy atoms - small enough to admit fragments, large enough
to exclude stray ions.

Run with the pyr1_docking env python.
"""
import csv
import glob
import json
import os
import re
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HC = os.path.join(ROOT, "results", "homolog_cavities")
OUT = os.path.join(ROOT, "results", "graft_candidates")
os.makedirs(OUT, exist_ok=True)

PYR1_CAVITY = 174.0        # A^3, the pocket we are trying to enlarge

#: Water, ions, cryoprotectants, buffers and crystallisation additives. Present in
#: a structure because of how it was grown, not because the protein binds them.
SOLVENT = {
    "HOH", "DOD", "WAT", "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "P6G",
    "SO4", "PO4", "NO3", "ACT", "ACY", "FMT", "CIT", "FLC", "TAR", "MES", "EPE",
    "TRS", "BTB", "IMD", "DMS", "MPD", "BME", "DTT", "TCE", "CAC", "MLI", "OXL",
    "NA", "K", "MG", "CA", "ZN", "MN", "FE", "FE2", "CD", "NI", "CO", "CU", "CU1",
    "CL", "BR", "IOD", "F", "LI", "CS", "RB", "SR", "BA", "HG", "AU", "PT", "PB",
    "SM", "EU", "GD", "YB", "AG", "AL", "XE", "KR", "NH4", "AZI", "CO3", "BCT",
    "UNX", "UNL", "UNK", "PGO", "PDO", "MRD", "IPA", "EOH", "MOH", "ACE", "NH2",
    "SIN", "SCN", "THJ", "PER", "BEZ", "BU3", "BUD", "12P", "15P", "33O",
    # buffers and cryoprotectants that are large enough to clear a heavy-atom floor
    # and therefore masquerade as ligands: CXS is CAPS, HEZ is hexane-1,6-diol
    "CXS", "HEZ", "CPS", "CHD", "C8E", "LDA", "BOG", "OGA", "HED", "DTV", "PGR",
}
#: Modified amino acids. MSE is selenomethionine, incorporated for phasing and part
#: of the POLYPEPTIDE - it is not a ligand. Left unfiltered it dominated the list,
#: putting 14 structures in whose only "ligand" was their own backbone.
MODIFIED_AA = {
    "MSE", "SEC", "PYL", "PTR", "SEP", "TPO", "CSO", "CSD", "CME", "CSS", "OCS",
    "KCX", "MLY", "MLZ", "M3L", "ALY", "HYP", "LLP", "NEP", "SAC", "CGU", "TYS",
    "PCA", "FME", "DAL", "DLE", "DVA", "ABA", "AIB", "NLE", "ORN", "SAR", "MHO",
}

#: Density that was modelled but never identified. Excluded from the ligand list
#: because we cannot name the chemistry - but reported separately, because a large
#: unidentified blob in a large pocket is itself evidence the pocket is occupied.
UNIDENTIFIED = {"UNL", "UNX", "UNK", "LIG"}

MIN_HEAVY = 6              # below this it is an ion or a fragment of buffer


def het_residues(path):
    """{resname: heavy atoms in its largest copy} from an mmCIF atom_site loop.

    ⚠ These files are mmCIF, not PDB. The first version sliced fixed PDB columns
    (line[17:20]) and every ligand came back named "." - a CIF null - which then
    passed the solvent filter and made 34 structures look ligand-bearing. The loop
    header is read here so the columns are located rather than assumed.
    """
    cols, in_loop, out = {}, False, defaultdict(lambda: defaultdict(int))
    for line in open(path, errors="ignore"):
        t = line.strip()
        if t.startswith("_atom_site."):
            cols[t.split(".", 1)[1].strip()] = len(cols)
            in_loop = True
            continue
        if in_loop and (t.startswith("HETATM") or t.startswith("ATOM")):
            f = t.split()
            if len(f) < len(cols):
                continue
            if f[cols["group_PDB"]] != "HETATM":
                continue
            if f[cols["type_symbol"]].upper() == "H":
                continue
            name = f[cols.get("auth_comp_id", cols["label_comp_id"])].strip().upper()
            if name in (".", "?"):
                name = f[cols["label_comp_id"]].strip().upper()
            # ⚠ altLoc MUST be part of the copy key. Keying on (chain, resid) alone
            # merges alternate conformations of ONE molecule into a single phantom
            # with twice the atoms: 5NON's 93H is a 21-heavy-atom ligand modelled in
            # two conformations in each of three chains, and came out as "42 heavy".
            alt = f[cols["label_alt_id"]]
            copy = (f[cols.get("auth_asym_id", cols["label_asym_id"])],
                    f[cols.get("auth_seq_id", cols["label_seq_id"])],
                    "" if alt in (".", "?") else alt)
            out[name][copy] += 1
        elif in_loop and t.startswith("#"):
            in_loop = False
    return {n: max(c.values()) for n, c in out.items()}


def main():
    # ---- the four counts ----
    hits = [l.rstrip("\n").split("\t") for l in
            open(os.path.join(ROOT, "results", "foldseek", "fs_hits.tsv"))]
    tm05 = [h for h in hits if float(h[3]) >= 0.5]
    af = [h for h in tm05 if h[1].startswith("af_")]
    exp = [h for h in tm05 if not h[1].startswith("af_")]

    rows = list(csv.DictReader(open(os.path.join(HC, "homolog_cavities.csv"))))
    measured_exp = [r for r in rows if not r["target"].startswith("af_")]

    # ---- ligands, from the raw coordinate files ----
    cav = {r["target"]: r for r in rows}
    cands, with_ligand = [], 0
    for f in sorted(glob.glob(os.path.join(HC, "raw", "*"))):
        tgt = os.path.basename(f).rsplit(".", 1)[0]
        if tgt.startswith("af_"):
            continue
        res = het_residues(f)
        ligs = {n: c for n, c in res.items()
                if n not in SOLVENT and n not in MODIFIED_AA
                and n not in UNIDENTIFIED and c >= MIN_HEAVY}
        unid = {n: c for n, c in res.items()
                if n in UNIDENTIFIED and c >= MIN_HEAVY}
        if not ligs and not unid:
            continue
        if ligs:
            with_ligand += 1
        r = cav.get(tgt, {})
        cands.append(dict(
            target=tgt, pdb=tgt[:4].upper(), chain=tgt[4:5],
            alntm=float(r.get("alntm", 0)), fident=float(r.get("fident", 0)),
            cavity=float(r.get("cavity_A3", 0)),
            ligands=", ".join(f"{n} ({c} heavy)" for n, c in
                              sorted(ligs.items(), key=lambda kv: -kv[1])) or "-",
            unidentified=", ".join(f"{n} ({c} heavy)" for n, c in
                                   sorted(unid.items(), key=lambda kv: -kv[1])) or "-",
            biggest=max(list(ligs.values()) + [0]),
            wall="".join(r.get(f"aa{p}", "-") for p in (59, 79, 94, 108, 120, 141)),
        ))

    print("=" * 74)
    print("START / SRPBCC (helix-grip) RELATIVES OF PYR1 -- what we actually have")
    print("=" * 74)
    print(f"  1. structural hits returned by foldseek vs CATH50 : {len(hits)}")
    print(f"  2. hits at alnTM >= 0.5, i.e. the same fold       : {len(tm05)}")
    print(f"       of which AlphaFold MODELS (no ligand at all) : {len(af)}")
    print(f"       of which SOLVED structures                   : {len(exp)}"
          f"  ({len({h[1][:4] for h in exp})} distinct PDB entries)")
    print(f"  3. solved structures fetched and cavity-measured  : {len(measured_exp)}"
          f"   (Arm 3, alnTM >= 0.61)")
    unid_only = [c for c in cands if c["ligands"] == "-"]
    print(f"  4. of those, carrying an IDENTIFIED ligand        : {with_ligand}")
    print(f"       plus, with unidentified density only          : {len(unid_only)}")
    print()
    print(f"  Note on 2 vs 3: Arm 3 used a stricter alnTM cut and only kept domains that")
    print(f"  fetched full-atom, so it is a subset. Note on 4: 37 of 147 raw files contain")
    print(f"  HETATM records, but most are water, glycerol, sulfate, ions or SELENOMETHIONINE")
    print(f"  (a modified residue, not a ligand) - {with_ligand} survive a solvent + modified-")
    print(f"  residue filter and a {MIN_HEAVY}-heavy-atom floor.")

    cands.sort(key=lambda c: -c["cavity"])
    print()
    print("=" * 74)
    print("GRAFTING CANDIDATES -- solved, same fold, real ligand, ranked by cavity")
    print("=" * 74)
    print(f"  PYR1's own pocket is {PYR1_CAVITY:.0f} A^3 for reference.")
    print()
    print(f"  {'PDB':<6}{'domain':<10}{'alnTM':>6}{'%id':>6}{'cavity':>8}{'x PYR1':>8}  "
          f"{'wall 59/79/94/108/120/141':<26}{'ligands':<22} unidentified")
    for c in cands:
        print(f"  {c['pdb']:<6}{c['target']:<10}{c['alntm']:>6.2f}{100*c['fident']:>6.0f}"
              f"{c['cavity']:>8.1f}{c['cavity']/PYR1_CAVITY:>8.1f}  {c['wall']:<26}{c['ligands']:<22} {c['unidentified']}")

    fetched = {os.path.basename(f).rsplit(".", 1)[0]
               for f in glob.glob(os.path.join(HC, "raw", "*"))}
    missing = sorted({h[1] for h in exp} - fetched)
    print()
    print("=" * 74)
    print("GAP -- solved domains at alnTM >= 0.5 that were never fetched")
    print("=" * 74)
    print(f"  Arm 3 kept {len(measured_exp)} of the {len(exp)} solved domains, because it used a")
    print(f"  stricter alnTM cut. {len(missing)} were never examined and may carry ligands:")
    for i in range(0, len(missing), 10):
        print("    " + "  ".join(missing[i:i + 10]))
    print("  -> fetching these is the cheapest way to extend the candidate list.")

    json.dump(dict(unfetched_solved=missing, total_hits=len(hits), same_fold=len(tm05), alphafold=len(af),
                   solved=len(exp), solved_pdbs=len({h[1][:4] for h in exp}),
                   measured=len(measured_exp), with_ligand=with_ligand,
                   candidates=cands),
              open(os.path.join(OUT, "graft_candidates.json"), "w"), indent=1)
    with open(os.path.join(OUT, "graft_candidates.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cands[0].keys()))
        w.writeheader()
        w.writerows(cands)
    print(f"\nwrote {OUT}/graft_candidates.csv and .json")


if __name__ == "__main__":
    main()
