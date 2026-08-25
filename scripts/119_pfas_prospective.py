#!/usr/bin/env python
r"""
119_pfas_prospective.py -- the §64e forcing rule, run prospectively on PFAS.

THE PRE-REGISTERED RULE (README §64e, written 2026-08-25 before sd09 was opened)

  For a target class, over round-1 clones whose ligand has Tanimoto >= 0.40 to the
  class, compute P(mutated) and P(top residue | mutated). PREDICT that any position
  with P(top | mutated) >= 0.85 on n >= 8 class clones will carry that residue in
  >= 80 % of that class's round-2 sensors, and may therefore be FORCED.

  It fails if it fires on a position whose round-2 frequency is < 80 %, or fires
  nowhere while a forceable position exists in the round-2 data.

TNT IS NOT RUN. Jannis: TNT's round-2 library was not built from hits on
chemically related ligands -- no near neighbours had hits -- but from ligands
sharing specific structural features. The rule keys on 2D whole-molecule
similarity, so its input does not exist for TNT. Running it there would test the
wrong thing; sd08 stays sealed.

⚠ THE FIRST RUN OF THIS SCRIPT USED THE WRONG ROUND-1 INPUT AND ITS VERDICT IS
VOID. It took round 1 to be sd03, the 194 hits from the 3,366-compound
Selleck/LATCA deck, which contains no PFAS at all -- so the class set came out
empty and the rule "passed" by having nothing to say. Tian in fact ran a separate
PFAS round-1 screen: an improved DSM-Hao library against a panel of 103 PFAS,
yielding sensors for 18 targets. That data is inside sd09 itself, separated by the
`mut_lib` column:

    mut_lib = DSM-Hao      89 clones, 18 ligands   <- ROUND 1, the design input
    mut_lib = PFOS_GenWT  154 clones, 25 ligands   <- ROUND 2, the answer sheet

Seven of the 25 round-2 targets had no round-1 hit at all, the same shape as
coumarin's 4 of 11 (§54a).

TWO PHASES, AND THE ORDER IS THE POINT
  Phase 1 reads sd03 and the target SMILES ONLY, and prints the predictions.
  Phase 2 then opens sd09 and scores them.
No threshold is touched between the two. Anything computed after Phase 2 is
labelled a post-hoc diagnostic and is not part of the test.
"""
import json
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_xlsx import table                                    # noqa: E402

ROOT = os.path.dirname(HERE)
SD = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
OUT = os.path.join(ROOT, "results", "secondary_library")
TANIMOTO = 0.40          # pre-registered
P_TOP = 0.85             # pre-registered
N_MIN = 8                # pre-registered
TARGET_R2 = 0.80         # pre-registered

#: The 25 PFAS screened in sd09, written from their names. Every one is checked
#: below against the fluorine count its own name states ("hexadecafluoro" = 16,
#: "dodecafluoro" = 12, "octafluoro" = 8), which catches a mistyped chain.
PFAS = {
    "Perfluorooctanesulfonic acid":
        "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "Perfluorohexanesulfonic Acid":
        "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "Perfluorooctanoic acid":
        "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "Perfluorononanoic acid":
        "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "Perfluorodecanoic acid":
        "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "Perfluorooctanesulfonamide":
        "NS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "Perfluorohexanesulfonamide":
        "NS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "Perfluorobutylsulfonamide":
        "NS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "N-Methylperfluorobutanesulfonamide":
        "CNS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "2H-Perfluoro-2-octenoic acid":
        "OC(=O)C=C(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "Perfluoro-3,6,9-trioxatridecanoic acid":
        "OC(=O)C(F)(F)OC(F)(F)C(F)(F)OC(F)(F)C(F)(F)OC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "1H,1H,9H,9H-Perfluoro-1,9-nonanediol":
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CO",
    "1H,1H,8H,8H-Dodecafluoro-1,8-octanediol":
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CO",
    "1H,1H,10H,10H-Perfluoro-1,10-decanediol":
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CO",
    "1H,1H,6H,6H-Perfluoro-1,6-hexandiol":
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)CO",
    "1H,1H,2H,3H,3H-Perfluorononane-1,2-diol":
        "OCC(O)CC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "Perfluoroheptanal hydrate (diol)":
        "OC(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "1H,1H,9H-Hexadecafluorononan-1-ol":
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)F",
    "1H,1H,8H-Perfluorooctan-1-ol":
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)F",
    "2H,2H,3H,3H-Perfluorononanoic acid":
        "OC(=O)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "2H,2H,3H,3H-Perfluorooctanoic acid":
        "OC(=O)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "2H,2H,3H,3H-Perfluoroheptanoic acid":
        "OC(=O)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "2H,2H,3H,3H-Perfluorohexanoic acid":
        "OC(=O)CCC(F)(F)C(F)(F)C(F)(F)F",
    "3,5,6-trichlorooctafluorohexanoic acid":
        "OC(=O)C(F)(F)C(F)(Cl)C(F)(F)C(F)(Cl)C(F)(F)Cl",
    "9-Chlorohexadecafluorononanoic Acid":
        "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)Cl",
}
#: fluorine counts stated by the names themselves
F_STATED = {"1H,1H,8H,8H-Dodecafluoro-1,8-octanediol": 12,
            "1H,1H,9H-Hexadecafluorononan-1-ol": 16,
            "3,5,6-trichlorooctafluorohexanoic acid": 8,
            "9-Chlorohexadecafluorononanoic Acid": 16}

POS18 = ["K59", "V81", "V83", "L87", "A89", "S92", "E94", "F108", "I110",
         "L117", "Y120", "S122", "E141", "F159", "A160", "V163", "V164", "N167"]


def clones(path, ligcol, smicol=None):
    hdr, recs = table(path)
    cols = [h for h in hdr if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]\d+", h)]
    out = []
    for r in recs:
        lig = (r.get(ligcol) or "").strip()
        if not lig or lig.startswith("Bold:") or lig.startswith("*"):
            continue                        # sd09 carries two footnote rows
        st = {p: p[0] for p in POS18}
        off = False
        for c in cols:
            v = (r.get(c) or "").strip().upper()
            if len(v) != 1 or v not in "ACDEFGHIKLMNPQRSTVWY":
                continue
            if c in POS18:
                st[c] = v
            elif v != c[0]:
                off = True
        out.append({"lig": lig, "st": st, "off": off,
                    "smiles": (r.get(smicol) or "").strip() if smicol else "",
                    "lib": (r.get("mut_lib") or "").strip()})
    return out


def two_part(cset, positions):
    """P(mutated) and P(top residue | mutated) at each position, over a clone set."""
    out = {}
    for p in positions:
        mut = [c["st"][p] for c in cset if c["st"][p] != p[0]]
        cnt = Counter(mut)
        top = cnt.most_common(1)
        out[p] = {"n_clones": len(cset), "n_mut": len(mut),
                  "p_mut": len(mut) / len(cset) if cset else 0.0,
                  "top": top[0][0] if top else None,
                  "top_frac": (top[0][1] / len(mut)) if mut else 0.0}
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    log = []

    def say(*t):
        s = " ".join(str(x) for x in t)
        print(s)
        log.append(s)

    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    say("=" * 80)
    say("PHASE 0. TARGET SET, VALIDATED BEFORE USE")
    say("=" * 80)
    tfp, bad = [], []
    for name, smi in PFAS.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            bad.append((name, "unparseable"))
            continue
        nf = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "F")
        want = F_STATED.get(name)
        if want is not None and nf != want:
            bad.append((name, f"{nf} F, name states {want}"))
            continue
        tfp.append((name, gen.GetFingerprint(mol),
                    rdMolDescriptors.CalcMolFormula(mol)))
    if bad:
        raise SystemExit(f"target SMILES failed validation: {bad}")
    say(f"   {len(tfp)} PFAS structures parsed; every name that states a fluorine")
    say(f"   count matches it ({', '.join(sorted(F_STATED))})")

    say("")
    say("=" * 80)
    say("PHASE 1. PREDICTIONS -- from ROUND-1 hits and target SMILES ONLY")
    say("=" * 80)
    r1 = clones(f"{SD}/pnas.2519924122.sd03(1).xlsx", "library_name",
                "canonical_smiles")
    # the PFAS round-1 screen, which lives in sd09 under mut_lib = DSM-Hao
    pf1 = [c for c in clones(f"{SD}/pnas.2519924122.sd09.xlsx", "chem_name")
           if c["lib"] == "DSM-Hao"]
    for c in pf1:
        c["smiles"] = PFAS.get(c["lig"], "")
        if not c["smiles"]:
            raise SystemExit(f"no SMILES for round-1 PFAS ligand {c['lig']!r}")
    say(f"   round-1 corpus: sd03 {len(r1)} clones ({len({c['lig'] for c in r1})} "
        f"ligands) + sd09/DSM-Hao {len(pf1)} PFAS clones "
        f"({len({c['lig'] for c in pf1})} ligands)")
    r1 = r1 + pf1
    sim = {}
    for c in r1:
        if not c["smiles"] or c["lig"].lower() in sim:
            continue
        m = Chem.MolFromSmiles(c["smiles"])
        if m:
            sim[c["lig"].lower()] = max(
                DataStructs.TanimotoSimilarity(gen.GetFingerprint(m), f)
                for _n, f, _fo in tfp)
    top = sorted(sim.items(), key=lambda x: -x[1])[:10]
    say(f"   round-1 ligands fingerprinted: {len(sim)}")
    say(f"   max Tanimoto to the PFAS class, best 10 round-1 ligands:")
    for l, v in top:
        say(f"     {v:.3f}  {l}")
    say(f"   median {sorted(sim.values())[len(sim)//2]:.3f}, "
        f"max {max(sim.values()):.3f}")
    cset = [c for c in r1 if sim.get(c["lig"].lower(), 0.0) >= TANIMOTO]
    say(f"   of the class set, {sum(1 for c in cset if c['lib']=='DSM-Hao')} are "
        f"PFAS round-1 hits and "
        f"{sum(1 for c in cset if c['lib']!='DSM-Hao')} are from the general deck")
    say("")
    say(f"   clones at Tanimoto >= {TANIMOTO}: {len(cset)}")
    fired = {}
    if len(cset) == 0:
        say("")
        say("   ==> THE CLASS SET IS EMPTY. No round-1 ligand comes within "
            f"{TANIMOTO} of")
        say("       any PFAS. The rule as written therefore fires NOWHERE, and its")
        say("       Phase-1 prediction is: NO POSITION IS FORCEABLE for PFAS.")
        say("       This is a legitimate outcome of the pre-registration, not a")
        say("       failure to run -- and it is falsifiable: if sd09's round-2 data")
        say("       DOES contain a forceable position, the rule has failed by being")
        say("       too conservative.")
    else:
        tp = two_part(cset, POS18)
        for p, v in tp.items():
            if v["top_frac"] >= P_TOP and v["n_mut"] >= N_MIN:
                fired[p] = v
                say(f"     FIRES {p}->{v['top']}  P(top|mut) {v['top_frac']:.2f} "
                    f"on n = {v['n_mut']}")
        if not fired:
            say("   ==> the rule fires nowhere.")

    say("")
    say("=" * 80)
    say("PHASE 2. OPENING sd09 (PFAS round 2) -- 245 clones, sealed until now")
    say("=" * 80)
    r2 = [c for c in clones(f"{SD}/pnas.2519924122.sd09.xlsx", "chem_name")
          if c["lib"] == "PFOS_GenWT"]
    say(f"   {len(r2)} round-2 clones over {len({c['lig'] for c in r2})} PFAS "
        f"ligands (mut_lib = PFOS_GenWT)")
    novel = {c["lig"] for c in r2} - {c["lig"] for c in pf1}
    say(f"   {len(novel)} of them had NO round-1 hit: "
        f"{', '.join(sorted(novel))}")
    lib = {}
    _, sd04 = table(f"{SD}/pnas.2519924122.sd04.xlsx")
    for r in sd04:
        if (r.get("Library") or "").strip() == "PFAS":
            lib[f"{r['WT']}{int(float(r['Position']))}"] = set(
                r["Amino Acids allowed for mutation"].strip())
    say(f"   Tian's PFAS library: {len(lib)} positions, "
        f"{sum(len(v) for v in lib.values())} substitutions, size "
        f"{eval('*'.join(str(1+len(v)) for v in lib.values())):,}")
    say("")
    say("   THE ANSWER -- which positions are actually forceable in round 2?")
    say("   ⚠ P(top|mut) is CONFOUNDED where the round-2 library offered only one")
    say("   substitution -- it is then 1.00 by construction, not by selection. The")
    say("   column that carries real information in both rounds is P(mutated),")
    say("   i.e. whether the winners kept wild-type when they could have.")
    say(f"   {'pos':<7}{'P(mutated)':>12}{'top':>6}{'P(top|mut)':>12}"
        f"{'carry top':>11}   forceable (>= 80 %)?")
    tp2 = two_part(r2, POS18)
    forceable = []
    for p in POS18:
        v = tp2[p]
        carry = v["p_mut"] * v["top_frac"]
        if v["n_mut"] == 0:
            continue
        ok = carry >= TARGET_R2
        if ok:
            forceable.append((p, v["top"], carry))
        say(f"   {p:<7}{v['p_mut']:>12.2f}{str(v['top']):>6}"
            f"{v['top_frac']:>12.2f}{carry:>11.2f}   {'YES' if ok else 'no'}")
    say("")
    say("=" * 80)
    say("VERDICT")
    say("=" * 80)
    say(f"   rule fired on   : {sorted(fired) if fired else 'NOTHING'}")
    say(f"   actually forceable: "
        f"{[f'{p}->{a} ({100*c:.0f}%)' for p, a, c in forceable] or 'NOTHING'}")
    say("")
    if not fired and not forceable:
        say("   ==> PASS, in the weak sense that matters least: the rule declined")
        say("       to fire and there was nothing to fire on. It was right, but the")
        say("       test had no opportunity to catch it being wrong.")
    elif not fired and forceable:
        say("   ==> FAIL (too conservative). A forceable position exists and the")
        say("       rule missed it -- exactly the second failure mode §64e names.")
    elif fired and all(p in [f[0] for f in forceable] for p in fired):
        say("   ==> PASS. Every position the rule fired on is genuinely forceable.")
    else:
        miss = [p for p in fired if p not in [f[0] for f in forceable]]
        say(f"   ==> FAIL (false positive). Fired on {miss}, which round 2 does not")
        say("       support. Forcing those would have cost sensors.")

    # ------------------------------------------------------------- phase 3
    say("")
    say("=" * 80)
    say("PHASE 3 (post-hoc diagnostic, NOT part of the pre-registered test)")
    say("=" * 80)
    say("   The rule declined because no position's round-1 residue IDENTITY was")
    say("   concentrated enough. But the other half of the two-part statistic --")
    say("   P(mutated), i.e. whether the winners keep wild-type -- transfers well.")
    say("")
    t1 = two_part([c for c in pf1], POS18)
    say(f"   {'pos':<7}{'r1 P(mut)':>10}{'r1 n':>6}{'r1 top|mut':>12}"
        f"{'r2 P(mut)':>11}{'r2 top|mut':>12}{'options':>10}")
    X, Y = [], []
    for p in POS18:
        a, b = t1[p], tp2[p]
        if a["n_mut"] == 0 and b["n_mut"] == 0:
            continue
        X.append(a["p_mut"])
        Y.append(b["p_mut"])
        say(f"   {p:<7}{a['p_mut']:>10.2f}{a['n_mut']:>6}"
            f"{str(a['top'] or '-'):>7}{a['top_frac']:>5.2f}"
            f"{b['p_mut']:>11.2f}{str(b['top'] or '-'):>7}{b['top_frac']:>5.2f}"
            f"{len(lib.get(p, '')) or '-':>10}")
    import numpy as np

    def sp(a, b):
        ra = np.argsort(np.argsort(a)).astype(float)
        rb = np.argsort(np.argsort(b)).astype(float)
        ra -= ra.mean()
        rb -= rb.mean()
        return float((ra * rb).sum() / np.sqrt((ra ** 2).sum() * (rb ** 2).sum()))

    say("")
    say(f"   Spearman( round-1 P(mutated), round-2 P(mutated) ) = {sp(X, Y):+.2f} "
        f"over n = {len(X)} positions")
    hi = [p for p in POS18 if tp2[p]["p_mut"] >= 0.80]
    say(f"   positions where round 2 nearly always drops wild-type: {hi}")
    for p in hi:
        rk = sorted([t1[q]["p_mut"] for q in POS18], reverse=True).index(
            t1[p]["p_mut"]) + 1
        say(f"     {p}: round-2 {tp2[p]['p_mut']:.2f}, round-1 was "
            f"{t1[p]['p_mut']:.2f} (rank {rk} of 18)")
    say("")
    say("   ⚠ AND THE PAYOFF IS SMALL HERE, for a structural reason worth keeping.")
    sz = 1
    dr = 1
    for p, v in lib.items():
        sz *= (1 + len(v))
        dr *= (len(v) if p in hi else 1 + len(v))
    say(f"   Dropping wild-type at all three: {sz:,} -> {dr:,} = {sz/dr:.1f}x only.")
    say("   Forcing pays in proportion to how SHALLOW the position's menu is: a")
    say("   position offering one substitution halves the library when wild-type")
    say("   goes, one offering eight saves 11 %. Coumarin's V163 offered exactly")
    say("   one residue (W) and was worth 2x; PFAS's K59/E94/Y120 offer 3/5/8.")
    say("")
    say("   So the two classes differ in KIND, not degree:")
    say("     coumarin  one position converges on one residue      -> force it, 2x")
    say("     PFAS      three positions drop wild-type but spread  -> 1.8x, no")
    say("               single residue to collapse to")

    with open(os.path.join(OUT, "pfas_prospective.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    with open(os.path.join(OUT, "pfas_prospective.json"), "w") as fh:
        json.dump({"fired": fired, "forceable": forceable,
                   "round2": tp2, "n_class_clones": len(cset)}, fh, indent=1)
    say(f"\n   written to {OUT}/pfas_prospective.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
