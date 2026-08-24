#!/usr/bin/env python
r"""
106_sensor_benchmark.py -- what chemistry do REAL PYR1 sensors actually use?

WHY
§47 proposed reframing the goal from "rank mutations by predicted ddG" (nothing
here reaches the needed precision) to "design a maximally enriched library,
scored by recall of known sensors in the top N". §32's pairwise steric method is
the only candidate engine, and it recovers 3 of 4 of PYR1^MANDI -- missing K59R,
whose benefit is a bidentate hydrogen bond and which does not clash at all.

The objection this script answers: if a whole CATEGORY of substitution is
invisible to a steric objective, then 75 % recall is not "pretty good", it is
blind to a class of chemistry that may be REQUIRED for binding. A library built
from a steric method alone would systematically exclude it.

So: across every characterised sensor we have, how much of the observed
chemistry is steric, and how much is charge?

THE FROZEN SPLIT -- declared here, before any result is computed
    DEV      sd03 (692 clones, mixed), sd07 (78, coumarin/imperatorin),
             Beltran 45 cannabinoid sensors
    HELD OUT sd08 (96 clones, TNT / dinitrotoluene)
             sd09 (245 clones, PFAS / perfluorooctanoic acid, PFOS)

The holdout is split by LIGAND CHEMISTRY, not at random. That is deliberately
harder: it asks whether a rule fitted on cannabinoids and coumarins transfers to
a nitroaromatic and a perfluorinated acid, which is exactly the generalisation a
design method has to make. Nothing in this project has ever been tested
prospectively (§47e), and PFAS/TNT are the only untouched sets we have -- so they
are not to be looked at while a method is being chosen.
"""
import json
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_xlsx import table                                    # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"

# Zamyatnin residue volumes, A^3
VOL = {"G": 60.1, "A": 88.6, "S": 89.0, "C": 108.5, "D": 111.1, "P": 112.7,
       "N": 114.1, "T": 116.1, "E": 138.4, "V": 140.0, "Q": 143.8, "H": 153.2,
       "M": 162.9, "I": 166.7, "L": 166.7, "K": 168.6, "R": 173.4, "F": 189.9,
       "Y": 193.6, "W": 227.8}
FORMAL = {"D": -1, "E": -1, "K": +1, "R": +1}      # H treated separately
POLAR = set("STNQYCWH")

DEV = {
    "sd03 (mixed screen)": (f"{SD}/pnas.2519924122.sd03(1).xlsx", "library_name"),
    "sd07 (coumarin)":     (f"{SD}/pnas.2519924122.sd07(1).xlsx", "compound"),
}
HELD = {
    "sd08 (TNT)":  (f"{SD}/pnas.2519924122.sd08.xlsx", "compound"),
    "sd09 (PFAS)": (f"{SD}/pnas.2519924122.sd09.xlsx", "chem_name"),
}


def clone_subs(path, ligcol):
    """Per-position columns are named like 'K59' = WT residue + native number.

    A blank cell means 'wild type at this position', NOT missing data -- these
    sheets record only what changed.
    """
    hdr, recs = table(path)
    poscols = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    out = []
    for r in recs:
        subs = []
        for c in poscols:
            v = (r.get(c) or "").strip().upper()
            if len(v) == 1 and v in VOL and v != c[0]:
                subs.append(f"{c[0]}{c[1:]}{v}")
        if subs:
            out.append({"ligand": (r.get(ligcol) or "?").strip(), "subs": subs})
    return out, poscols


def classify(sub):
    wt, pos, mu = sub[0], sub[1:-1], sub[-1]
    dv = VOL[mu] - VOL[wt]
    qa, qb = FORMAL.get(wt, 0), FORMAL.get(mu, 0)
    if qa != qb:
        kind = "CHARGE"
    elif (mu in POLAR) != (wt in POLAR):
        kind = "polar"
    else:
        kind = "steric"
    return kind, dv, int(pos)


def report(name, records, out):
    subs = [s for r in records for s in r["subs"]]
    uniq = sorted(set(subs))
    out.append(f"\n{'='*78}\n{name}: {len(records)} clones, {len(subs)} substitutions, "
               f"{len(uniq)} distinct\n{'='*78}")
    cnt, byk = Counter(), defaultdict(list)
    for s in uniq:
        k, dv, p = classify(s)
        cnt[k] += 1
        byk[k].append((s, dv, p))
    for k in ("steric", "CHARGE", "polar"):
        n = cnt[k]
        out.append(f"  {k:<8} {n:3d} distinct ({100*n/len(uniq):4.1f}%)")
    # the class a steric objective cannot see: charge change with little volume change
    invisible = [(s, dv) for s, dv, p in byk["CHARGE"] if abs(dv) < 25]
    out.append(f"\n  CHARGE change with |dVolume| < 25 A^3 -- near-isosteric, so a")
    out.append(f"  volume/overlap objective has nothing to key on: "
               f"{len(invisible)} of {len(uniq)} distinct ({100*len(invisible)/len(uniq):.1f}%)")
    if invisible:
        out.append("    " + ", ".join(f"{s}({dv:+.0f})" for s, dv in sorted(invisible)[:18]))
    # weighted by how often they actually appear
    ci = sum(subs.count(s) for s, _ in invisible)
    out.append(f"  weighted by occurrence: {ci} of {len(subs)} "
               f"({100*ci/len(subs):.1f}%) of all substitutions seen")
    # position 59 specifically
    p59 = sorted({s for s in uniq if s[1:-1] == "59"})
    if p59:
        out.append(f"\n  position 59 (WT Lys, charged): {len(p59)} distinct -- "
                   + ", ".join(p59))
        n59 = sum(subs.count(s) for s in p59)
        out.append(f"    appears in {n59} substitutions "
                   f"({100*n59/len(subs):.1f}% of all)")
    return cnt, uniq


def main():
    show_held = "--reveal-holdout" in sys.argv
    out = []
    out.append("FROZEN SPLIT (declared in this file's header, before any result)")
    out.append("  DEV      : sd03, sd07, Beltran-45")
    out.append("  HELD OUT : sd08 (TNT), sd09 (PFAS)   <- not to be inspected while")
    out.append("             a method is being chosen; pass --reveal-holdout to see it")

    allsubs = {}
    for name, (path, ligcol) in DEV.items():
        recs, _ = clone_subs(path, ligcol)
        cnt, uniq = report(name, recs, out)
        allsubs[name] = uniq

    # Beltran 45
    b = json.load(open(os.path.join(ROOT, "data", "beltran", "win_sensors.json")))
    brecs = [{"ligand": s["ligand"], "subs": s["muts"]} for s in b["sensors"]]
    report("Beltran-45 (cannabinoids)", brecs, out)

    if show_held:
        out.append(f"\n\n{'#'*78}\n# HELD-OUT SETS REVEALED\n{'#'*78}")
        for name, (path, ligcol) in HELD.items():
            recs, _ = clone_subs(path, ligcol)
            report(name, recs, out)
    else:
        out.append(f"\n\n{'#'*78}")
        out.append("# HELD OUT (sd08 TNT, sd09 PFAS) NOT READ -- "
                   f"{sum(1 for _ in HELD)} sets reserved for the prospective test")
        out.append(f"{'#'*78}")

    txt = "\n".join(out)
    print(txt)
    d = os.path.join(ROOT, "results", "sensor_benchmark")
    os.makedirs(d, exist_ok=True)
    fn = "chemistry_dev.txt" if not show_held else "chemistry_with_holdout.txt"
    open(os.path.join(d, fn), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
