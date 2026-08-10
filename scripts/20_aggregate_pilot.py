#!/usr/bin/env python
"""
20_aggregate_pilot.py -- join the two pilot arms into one decision table.

Arm 1 (06_rosetta_cavity_scan.py) asks: does the mutation actually OPEN a cavity
once Rosetta is allowed to repack the neighbourhood? The rigid truncation scan
(script 02) deleted side chains without letting anything relax, so its volumes
are upper bounds. A position only survives if the cavity is still there after
FastRelax -- i.e. the surrounding side chains do not simply collapse into the
hole. The cost of that opening is the score change.

Arm 2 (07_rosetta_ratchet.py) asks: does the mutation break the READOUT? PYR1-
HAB1 is a ratchet (see README section 1): HAB1 contributes essentially nothing to
the pocket and reads the closed-state SURFACE. So a variant is only useful if the
gate/latch stay put, the interface energy survives, and W385 still wedges against
the ligand site.

A variant must pass BOTH. Volume that costs the interface is not progress.

Outputs
-------
results/20_pilot_summary.csv   one row per variant, all metrics, mean +- sd
results/20_pilot_summary.md    the same, ranked and annotated, for the README

Metric notes (the ones that are easy to misread)
------------------------------------------------
  d_score_ref2015  Arm-1 apo score MINUS the WT apo score, in REU. POSITIVE is
                   worse (Rosetta scores are energies). This is a whole-protein
                   number dominated by the removed side chain's own terms, so a
                   few REU of penalty for deleting a large buried residue is
                   expected and is NOT evidence of instability; what matters is
                   whether it is a few REU or a few tens.
  dG_separated     interface energy from InterfaceAnalyzerMover with
                   set_pack_separated(True), so the unbound state is repacked.
                   NEGATIVE is favourable. WT is about -67 REU.
  dSASA            buried interface area, A^2. Falls if the complex loosens.
  sc_value         shape complementarity, 0-1. WT about 0.754.
  gate/latch_rmsd  CA RMSD of residues 85-89 / 115-117 to the relaxed WT
                   complex. This is the direct readout-integrity test.
  w385_aba_min_dist  HAB1 Trp385 to ABA, A. Crystal 4.67; relaxed WT about 4.6.

Run with the esmfold2 env python.
"""
import os, sys, json, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from variants import PANEL, WT  # noqa: E402

CAV_DIR = os.path.join(ROOT, "results", "cavity_scan")
RAT_DIR = os.path.join(ROOT, "results", "ratchet")
OUT_CSV = os.path.join(ROOT, "results", "20_pilot_summary.csv")
OUT_MD = os.path.join(ROOT, "results", "20_pilot_summary.md")

# tolerances for the PASS/FAIL call. Set from the WT run's own replicate spread
# plus a margin -- see README section 5. Deliberately permissive on score
# (deleting a buried Phe legitimately costs REU) and strict on the readout.
GATE_MAX = 0.60      # A; WT ~0.18, so this is ~3x the WT deviation
LATCH_MAX = 1.80     # A; WT ~1.23 (the latch is intrinsically mobile here)
DG_MAX_LOSS = 12.0   # REU of interface energy we are willing to give up
W385_MAX = 6.50      # A; beyond this W385 has lost contact with the ligand site
CAV_MIN_GAIN = 25.0  # A^3 over WT to count as "opened"


def load(d):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "*.json"))):
        j = json.load(open(f))
        out[j["variant"]] = j
    return out


def fmt_muts(muts):
    """[[59,'A'],[108,'A']] -> 'K59A+F108A' using the WT identities."""
    return "+".join(f"{WT.get(int(p), '?')}{int(p)}{a}" for p, a in muts)


def agg(reps, key):
    v = [r[key] for r in reps if key in r and r[key] is not None]
    if not v:
        return float("nan"), float("nan"), 0
    return float(np.mean(v)), float(np.std(v, ddof=1)) if len(v) > 1 else 0.0, len(v)


cav, rat = load(CAV_DIR), load(RAT_DIR)
print(f"loaded {len(cav)} cavity-scan and {len(rat)} ratchet variants")
expected = {name for name, _ in PANEL} | {"WT"}   # PANEL is [(name, [(pos, aa)...])]
missing = expected - set(cav)
if missing:
    print(f"  !! missing from Arm 1: {sorted(missing)}")
missing_r = expected - set(rat)
if missing_r:
    print(f"  !! missing from Arm 2: {sorted(missing_r)}")

if "WT" not in cav or "WT" not in rat:
    sys.exit("WT reference missing -- cannot compute deltas")

wt_cav = agg(cav["WT"]["reps"], "cavity_A3")[0]
wt_score = agg(cav["WT"]["reps"], "score_ref2015")[0]
wt_dg = agg(rat["WT"]["reps"], "dG_separated")[0]
wt_dsasa = agg(rat["WT"]["reps"], "dSASA")[0]
wt_sc = agg(rat["WT"]["reps"], "sc_value")[0]
print(f"WT reference: cavity {wt_cav:.1f} A^3, score {wt_score:.2f} REU, "
      f"dG_sep {wt_dg:.2f}, dSASA {wt_dsasa:.0f}, sc {wt_sc:.3f}")

rows = []
for name in sorted(set(cav) | set(rat)):
    c = cav.get(name, {}).get("reps", [])
    r = rat.get(name, {}).get("reps", [])
    muts = (cav.get(name) or rat.get(name) or {}).get("mutations", [])

    cv, cv_sd, n_c = agg(c, "cavity_A3")
    sc_, sc_sd, _ = agg(c, "score_ref2015")
    rms, _, _ = agg(c, "ca_rmsd_to_input")
    dg, dg_sd, n_r = agg(r, "dG_separated")
    ds, _, _ = agg(r, "dSASA")
    shc, _, _ = agg(r, "sc_value")
    gt, _, _ = agg(r, "gate_rmsd")
    lt, _, _ = agg(r, "latch_rmsd")
    w3, _, _ = agg(r, "w385_aba_min_dist")

    d_cav = cv - wt_cav
    d_score = sc_ - wt_score
    d_dg = dg - wt_dg          # positive = interface got WORSE

    fails = []
    if name != "WT":
        if not np.isnan(gt) and gt > GATE_MAX:
            fails.append(f"gate {gt:.2f}")
        if not np.isnan(lt) and lt > LATCH_MAX:
            fails.append(f"latch {lt:.2f}")
        if not np.isnan(d_dg) and d_dg > DG_MAX_LOSS:
            fails.append(f"dG +{d_dg:.1f}")
        if not np.isnan(w3) and w3 > W385_MAX:
            fails.append(f"W385 {w3:.2f}")

    rows.append(dict(
        _fails=fails,
        variant=name, mutations=fmt_muts(muts), n_cav=n_c, n_rat=n_r,
        cavity_A3=cv, cavity_sd=cv_sd, d_cavity_A3=d_cav,
        score_ref2015=sc_, score_sd=sc_sd, d_score_ref2015=d_score,
        ca_rmsd_to_input=rms,
        dG_separated=dg, dG_sd=dg_sd, d_dG_separated=d_dg,
        dSASA=ds, d_dSASA=ds - wt_dsasa, sc_value=shc, d_sc_value=shc - wt_sc,
        gate_rmsd=gt, latch_rmsd=lt, w385_aba_min_dist=w3,
        verdict=""))

# ---------------- does either arm actually discriminate? ----------------
# A metric only carries information if the spread BETWEEN variants exceeds the
# spread between replicates of the SAME variant. With nstruct=3 this is a crude
# F-like ratio, but it is the difference between "no variant breaks the readout"
# and "this assay cannot tell whether a variant breaks the readout" -- which are
# opposite conclusions drawn from the same table.
def discrimination(store, key):
    within, means = [], []
    for j in store.values():
        v = [r[key] for r in j["reps"] if key in r and r[key] is not None]
        if len(v) > 1:
            within.append(np.std(v, ddof=1))
        if v:
            means.append(np.mean(v))
    w = float(np.mean(within)) if within else float("nan")
    b = float(np.std(means, ddof=1)) if len(means) > 1 else float("nan")
    return w, b, (b / w if w else float("nan"))


print("\n" + "=" * 78)
print("DISCRIMINATION CHECK -- between-variant spread vs replicate noise")
print("=" * 78)
print(f"{'metric':<26}{'arm':>6}{'within sd':>11}{'between sd':>12}{'ratio':>8}  ")
DISCRIM = {}
for store, arm, keys in ((cav, "1", ["cavity_A3", "score_ref2015"]),
                         (rat, "2", ["dG_separated", "dSASA", "sc_value",
                                     "gate_rmsd", "latch_rmsd",
                                     "w385_aba_min_dist"])):
    for k in keys:
        w, b, r = discrimination(store, k)
        DISCRIM[k] = r
        flag = "informative" if r >= 2.0 else "NOT DISCRIMINATING (noise >= signal)"
        print(f"{k:<26}{arm:>6}{w:>11.3f}{b:>12.3f}{r:>8.2f}  {flag}")
print("""
ratio = (sd of the per-variant means) / (mean within-variant replicate sd).
Below ~2 the metric cannot separate variants at this nstruct: the differences in
the table are replicate noise, so a variant showing no change is NOT evidence
that nothing changed.""")

READOUT_KEYS = ["dG_separated", "gate_rmsd", "latch_rmsd", "w385_aba_min_dist"]
ARM2_BLIND = all(DISCRIM.get(k, 0) < 2.0 for k in READOUT_KEYS)
if ARM2_BLIND:
    print("""
*** ARM 2 IS BLIND AT nstruct=3. ***
Every readout metric has between-variant spread at or below its own replicate
noise, so the ratchet arm currently discriminates NOTHING. Two causes, both real:
  (a) nstruct=3 is far too few for a Rosetta interface energy, whose replicate
      sd here is ~2.8 REU -- larger than the entire between-variant range;
  (b) FastRelax was run with -relax:constrain_relax_to_start_coords, which pins
      every CA to the input complex. Gate and latch CA RMSD are therefore
      bounded BY CONSTRUCTION (within-variant sd 0.004 A), so a near-constant
      gate RMSD of ~0.18 A across all 41 variants measures the restraint, not
      the protein.
The Arm-2 columns below must NOT be read as "the readout survives". They are
uninformative. Arm 1 is unaffected -- its cavity ratio is ~9, i.e. genuine
signal -- so the volume ranking stands on its own.""")

# verdicts are assigned only now, because what a clean Arm-2 row MEANS depends
# on whether Arm 2 can discriminate at all.
ok = "readout untested" if ARM2_BLIND else "readout intact"
for r in rows:
    if r["variant"] == "WT":
        r["verdict"] = "reference"
    elif r["_fails"] and not ARM2_BLIND:
        r["verdict"] = "FAIL readout: " + ", ".join(r["_fails"])
    elif r["_fails"]:
        r["verdict"] = "cavity opened; flagged " + ", ".join(r["_fails"]) + " (arm 2 blind)"
    elif r["d_cavity_A3"] >= CAV_MIN_GAIN:
        r["verdict"] = f"opens cavity, {ok}"
    elif r["d_cavity_A3"] > 0:
        r["verdict"] = f"little volume, {ok}"
    else:
        r["verdict"] = "no gain"
    del r["_fails"]

cols = list(rows[0])
with open(OUT_CSV, "w") as fh:
    fh.write(",".join(cols) + "\n")
    for r in rows:
        fh.write(",".join(
            (f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]))
            for c in cols) + "\n")
print(f"\nwrote {OUT_CSV}  ({len(rows)} rows)")

# ---------------- ranked report ----------------
ranked = sorted((r for r in rows if r["variant"] != "WT"),
                key=lambda r: -r["d_cavity_A3"])
wt_row = next(r for r in rows if r["variant"] == "WT")

hdr = (f"{'variant':<26}{'cav A3':>9}{'dCav':>8}{'dREU':>8}"
       f"{'dG_sep':>9}{'ddG':>7}{'sc':>7}{'gate':>7}{'latch':>7}{'W385':>7}  verdict")
lines = [hdr, "-" * len(hdr)]


def fmt(r):
    return (f"{r['variant']:<26}{r['cavity_A3']:>9.1f}{r['d_cavity_A3']:>+8.1f}"
            f"{r['d_score_ref2015']:>+8.2f}{r['dG_separated']:>9.2f}"
            f"{r['d_dG_separated']:>+7.1f}{r['sc_value']:>7.3f}"
            f"{r['gate_rmsd']:>7.2f}{r['latch_rmsd']:>7.2f}"
            f"{r['w385_aba_min_dist']:>7.2f}  {r['verdict']}")


lines.append(fmt(wt_row))
lines += [fmt(r) for r in ranked]
report = "\n".join(lines)
print("\n" + report)

n_open = sum(1 for r in ranked if r["d_cavity_A3"] >= CAV_MIN_GAIN)
summary = (f"\n{n_open} of {len(ranked)} variants open the cavity by "
           f">={CAV_MIN_GAIN:.0f} A^3 after full repacking (Arm 1, ratio "
           f"{DISCRIM['cavity_A3']:.1f} -- real signal).")
if ARM2_BLIND:
    summary += ("\nNo variant can be called safe or unsafe for the readout: Arm 2 "
                "does not discriminate\nat nstruct=3 with CA constraints. Re-run "
                "it before any variant is ordered.")
print(summary)

with open(OUT_MD, "w") as fh:
    fh.write("# Pilot summary -- Arm 1 (cavity) x Arm 2 (ratchet)\n\n")
    fh.write(f"Generated by `scripts/20_aggregate_pilot.py` from "
             f"`results/cavity_scan/*.json` and `results/ratchet/*.json`.\n\n")
    fh.write(f"WT reference (mean of {wt_row['n_cav']} relax replicates): "
             f"cavity {wt_cav:.1f} A^3, apo score {wt_score:.2f} REU, "
             f"dG_separated {wt_dg:.2f} REU, dSASA {wt_dsasa:.0f} A^2, "
             f"sc {wt_sc:.3f}.\n\n")
    fh.write(f"PASS criteria: cavity gain >= {CAV_MIN_GAIN:.0f} A^3 AND "
             f"gate RMSD <= {GATE_MAX} A AND latch RMSD <= {LATCH_MAX} A AND "
             f"interface energy loss <= {DG_MAX_LOSS} REU AND "
             f"W385-ABA <= {W385_MAX} A.\n\n")
    fh.write("dCav = A^3 over WT (positive = bigger pocket). dREU = apo score "
             "over WT (positive = worse). ddG = interface energy over WT "
             "(positive = worse interface).\n\n```\n" + report + "\n```\n")
    fh.write(summary + "\n")
print(f"wrote {OUT_MD}")
