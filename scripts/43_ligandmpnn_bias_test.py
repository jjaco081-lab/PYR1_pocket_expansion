#!/usr/bin/env python
"""
43_ligandmpnn_bias_test.py -- is LigandMPNN's poor recovery on the Tian coumarin
set caused by a conservative amino-acid prior, or by the docked pose being wrong?

THE HYPOTHESIS UNDER TEST
-------------------------
Baseline LigandMPNN recovers 28.7% of experimental top-1 residues (33/115), while
a ligand-blind oracle built from the same wetlab data reaches 79.1%. Inspecting
the errors showed a systematic direction: the model proposes conservative,
natural-looking substitutions where the real selected sequences install bulk.

    pos 163  WT V   LigandMPNN I (11/11)   experiment W (11/11)
    pos 164  WT V   LigandMPNN V (11/11)   experiment S / F / E
    pos 120  WT Y   LigandMPNN L (10/10)   experiment M / A

If that is a prior problem, adding a logit bias toward bulky aromatics should
recover those positions. If it is a pose problem -- the pocket was docked in the
wild-type backbone, so the model is packing against a ligand that sits where it
would not sit once the pocket is redesigned -- then biasing will move residues
toward bulk indiscriminately and NOT improve agreement.

WHY THE OBVIOUS METRIC WOULD RIG THE TEST
-----------------------------------------
Biasing toward W trivially makes position 163 return W. Scoring "did 163 become W"
would therefore confirm the hypothesis no matter what is true. The primary metric
here is instead OVERALL top-1 recovery across all designable positions, so a bias
that fixes 163 by breaking five other positions shows up as no gain -- which is
the honest outcome. Per-position deltas are reported as diagnostics only.

Secondary metric: recovery restricted to the positions where the experiment is
ligand-INVARIANT versus those where it VARIES across ligands. A prior fix should
lift the invariant positions. Only a pose fix can lift the variable ones, because
only the ligand distinguishes them.

DESIGN
------
For each ligand, LigandMPNN is run on that ligand's existing selected pose
complex -- the same input as the baseline, so bias level is the only variable.
Bias is applied equally to W/F/Y at several strengths including 0.0, which
reproduces the baseline as an internal control (same seed, so any drift means
the comparison is broken).

Usage:
    python 43_ligandmpnn_bias_test.py [--n 30] [--dry-run]
"""
import argparse, csv, json, os, re, subprocess, sys, collections

BENCH = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
LMPNN = os.path.join(BENCH, "envs", "LigandMPNN")
PY = "/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python"
MODEL = os.path.join(LMPNN, "model_params", "ligandmpnn_v_32_010_25.pt")
POSITIONS = os.path.join(BENCH, "results", "design_sets", "wetlab_library_positions.txt")
TRUTH = os.path.join(BENCH, "results", "analysis", "per_ligand_recovery.csv")
OUTROOT = os.path.join(BENCH, "results", "bias_test")

BIAS_LEVELS = [0.0, 1.0, 2.0, 3.0]
BULKY = ["W", "F", "Y"]
SEED = 37

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=30, help="sequences per run")
ap.add_argument("--dry-run", action="store_true")
args = ap.parse_args()


def truth_table():
    """{ligand: {pos: experimental_top_aa}} and the ligand-variance split."""
    rows = list(csv.DictReader(open(TRUTH)))
    t = collections.defaultdict(dict)
    for r in rows:
        t[r["ligand"]][int(r["position"])] = r["real_top_aa"]
    byposn = collections.defaultdict(set)
    for r in rows:
        byposn[int(r["position"])].add(r["real_top_aa"])
    invariant = {p for p, s in byposn.items() if len(s) == 1}
    variable = {p for p, s in byposn.items() if len(s) > 1}
    return t, invariant, variable


def find_complex(ligand):
    """The pose complex the baseline used for this ligand."""
    d = os.path.join(BENCH, "results", "by_ligand", ligand, "selected_pose_complexes")
    if not os.path.isdir(d):
        d = os.path.join(BENCH, "results", "selected_pose_complexes")
    hits = []
    for root, _, files in os.walk(d):
        hits += [os.path.join(root, f) for f in files if f.endswith(".pdb")]
    return sorted(hits)[0] if hits else None


def run_lmpnn(pdb, outdir, bias):
    os.makedirs(outdir, exist_ok=True)
    resids = open(POSITIONS).read().split()
    cmd = [PY, os.path.join(LMPNN, "run.py"),
           "--model_type", "ligand_mpnn",
           "--checkpoint_ligand_mpnn", MODEL,
           "--pdb_path", pdb,
           "--out_folder", outdir,
           "--redesigned_residues", " ".join(resids),
           "--number_of_batches", str(args.n),
           "--batch_size", "1",
           "--temperature", "0.1",
           "--seed", str(SEED),
           "--ligand_mpnn_use_side_chain_context", "1"]
    if bias > 0:
        cmd += ["--bias_AA", ",".join(f"{a}:{bias}" for a in BULKY)]
    if args.dry_run:
        print("  " + " ".join(cmd)); return True
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=LMPNN)
    if r.returncode != 0:
        print(f"  !! failed: {r.stderr.strip().splitlines()[-1] if r.stderr else '?'}")
        return False
    return True


def seq_index_map(pdb, chain="A"):
    """{residue number -> index into the emitted sequence}.

    LigandMPNN writes the FULL chain sequence, not just the designed positions,
    and PYR1 pose files do not start at residue 1. Indexing the sequence by
    anything other than the actual CA order silently reads the wrong columns --
    which is exactly what an earlier version of this script did, producing a flat
    0.9% at every bias level because characters 0-18 are not designable.
    """
    order = []
    for line in open(pdb):
        if line.startswith("ATOM") and line[12:16] == " CA " and line[21] == chain:
            n = int(line[22:26])
            if not order or order[-1] != n:
                order.append(n)
    return {n: i for i, n in enumerate(order)}


def consensus(outdir, resids, pdb):
    """Modal residue per designed position across the sampled sequences."""
    fa = os.path.join(outdir, "seqs")
    seqs = []
    for root, _, files in os.walk(fa):
        for f in files:
            if not f.endswith((".fa", ".fasta")):
                continue
            keep = False
            for line in open(os.path.join(root, f)):
                if line.startswith(">"):
                    # the first record echoes the input sequence; skip it
                    keep = "id=" in line
                elif keep:
                    seqs.append(line.strip().split(":")[0])
    if not seqs:
        return {}
    imap = seq_index_map(pdb)
    out = {}
    for x in resids:
        p = int(re.sub(r"\D", "", x))
        i = imap.get(p)
        if i is None:
            continue
        col = [s[i] for s in seqs if len(s) > i]
        if col:
            out[p] = collections.Counter(col).most_common(1)[0][0]
    return out


truth, INVAR, VARY = truth_table()
ligands = sorted(truth)
resids = open(POSITIONS).read().split()
print(f"{len(ligands)} ligands, {len(resids)} designable positions, "
      f"{len(INVAR)} ligand-invariant / {len(VARY)} ligand-variable")
print(f"bias levels {BIAS_LEVELS} applied to {'/'.join(BULKY)}, n={args.n}, seed={SEED}\n")

results = []
for bias in BIAS_LEVELS:
    hit = tot = hi = ti = hv = tv = 0
    per_pos = collections.Counter()
    for lig in ligands:
        pdb = find_complex(lig)
        if not pdb:
            print(f"  [skip] no pose complex for {lig}"); continue
        od = os.path.join(OUTROOT, f"bias{bias:g}", lig.replace(",", "_").replace(" ", "_"))
        if not run_lmpnn(pdb, od, bias):
            continue
        if args.dry_run:
            continue
        pred = consensus(od, resids, pdb)
        for p, exp in truth[lig].items():
            if p not in pred:
                continue
            ok = pred[p] == exp
            tot += 1; hit += ok
            if p in INVAR: ti += 1; hi += ok
            else:          tv += 1; hv += ok
            if ok: per_pos[p] += 1
    if args.dry_run:
        continue
    results.append(dict(bias=bias, hit=hit, tot=tot,
                        inv=(hi, ti), var=(hv, tv), per_pos=dict(per_pos)))
    f = lambda a, b: f"{a}/{b} = {a/b:.1%}" if b else "n/a"
    print(f"bias {bias:>4}:  overall {f(hit,tot):<18} "
          f"invariant {f(hi,ti):<18} variable {f(hv,tv)}")

if results:
    os.makedirs(OUTROOT, exist_ok=True)
    with open(os.path.join(OUTROOT, "bias_test_summary.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    base = results[0]
    best = max(results, key=lambda r: r["hit"] / max(r["tot"], 1))
    print(f"\nbaseline (bias 0) {base['hit']/max(base['tot'],1):.1%} -> "
          f"best (bias {best['bias']:g}) {best['hit']/max(best['tot'],1):.1%}")
    print("\nINTERPRETATION")
    print("  overall recovery rises with bias      -> conservative prior was a real cause")
    print("  invariant rises, variable flat        -> prior fixed; pose still wrong")
    print("  neither rises                         -> pose is the problem, rebuild on an ensemble")
    print(f"\nwrote {OUTROOT}/bias_test_summary.json")
