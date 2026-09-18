#!/usr/bin/env python
r"""
238_identify_hits.py -- what ARE the donor hits?

Jannis: "tell me the identities of the best hits. I am wondering if they are all
PYLs." If they are, the graft arm is only recovering PYR1's own family and has
found nothing new; if they are not, the non-PYL donors are the interesting ones
because their pockets are genuinely different architecture.

CATH domain names resolve through RCSB, AlphaFold targets through UniProt.
⚠ 3QRZ is OBSOLETE (superseded by 4JDL 2013-03-13) and is labelled as such.
"""
import json, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PYL = re.compile(r"PYR|PYL|abscisic|ABA recept|RCAR", re.I)


def curl(url):
    try:
        r = subprocess.run(["curl", "-s", "--max-time", "25", url],
                           capture_output=True, text=True, timeout=30)
        return json.loads(r.stdout) if r.stdout.strip() else None
    except Exception:                                            # noqa: BLE001
        return None


def name_of(target):
    if target.startswith("af_"):
        acc = target.split("_")[1]
        d = curl(f"https://rest.uniprot.org/uniprotkb/{acc}.json"
                 f"?fields=protein_name,organism_name")
        if not d:
            return acc, "(lookup failed)", ""
        pd = d.get("proteinDescription", {})
        nm = (pd.get("recommendedName", {}).get("fullName", {}).get("value")
              or (pd.get("submissionNames") or [{}])[0]
                 .get("fullName", {}).get("value") or "?")
        return acc, nm, d.get("organism", {}).get("scientificName", "")
    pdb = target[:4].upper()
    d = curl(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}")
    if not d:
        return pdb, "(lookup failed)", ""
    return pdb, d.get("struct", {}).get("title", "?"), ""


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        ROOT, "results", "schema_chimera", "junctions_all_flank2.json")
    rows = json.load(open(src))
    cav_cut = float(sys.argv[2]) if len(sys.argv) > 2 else 164.4
    clean_cut = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    hits = [r for r in rows if r["cavity"] > cav_cut and r["clean"] >= clean_cut]
    hits.sort(key=lambda r: -r["cavity"])
    print(f"{len(rows)} donors scored; {len(hits)} with cavity > {cav_cut} "
          f"AND >= {clean_cut}/8 clean junctions\n")
    npyl = 0
    for r in hits:
        acc, nm, org = name_of(r["target"])
        is_pyl = bool(PYL.search(nm))
        npyl += is_pyl
        tag = "  [PYL/PYR family]" if is_pyl else ""
        obs = "  ⚠ OBSOLETE" if r.get("obsolete") else ""
        print(f"{r['target'][:34]:<36} {r['cavity']:>7.1f} A^3  "
              f"{r['clean']}/{r['n_junction']} clean  meanJ {r['mean_junction']:.2f}  "
              f"id {r['ident']:.0f}%{obs}")
        print(f"    {acc}: {nm[:88]}")
        if org:
            print(f"    {org}{tag}")
        elif tag:
            print(f"   {tag}")
        time.sleep(0.15)
    print(f"\n{npyl} of {len(hits)} are PYR/PYL-family; "
          f"{len(hits)-npyl} are NOT")
    return 0


if __name__ == "__main__":
    sys.exit(main())
