#!/usr/bin/env python
"""
48_stage1_run.py -- stage 1: can LigandMPNN recover PYR1^MANDI from WT plus a
perfectly placed ligand?

Ground truth (4WVO vs 3QN1): K59R, V81I, F108A, F159L.

ARMS
----
Three structures from 47_stage1_inputs.py x two alphabets:

  wt_mandi      WT backbone, mandipropamid in its true pose. Ligand clashes with
                WT side chains; relieving that is the design problem.
  polygly_mandi same, 15 designable positions truncated to Gly. No clash signal,
                so the pocket must be rebuilt. Harder and more honest.
  wt_aba        WT + ABA. THE NULL. Correct answer is zero mutations.

  dsm_hao       identities restricted to the DSM-Hao library (plus WT)
  free          all 20, the ceiling

SCORING -- why it is not recall-at-N
------------------------------------
Sampled sequences are NOT independent; designers mode-collapse, so N samples may
carry far fewer effective draws and a null built on N independent draws is too
harsh. Two consequences:

1. The primary statistic is the per-position probability mass LigandMPNN assigns
   to the TRUE residue, read from `sampling_probs` -- the full output distribution
   rather than a sample of it, so it is N-independent.

2. The primary null is the LIGAND-SWAP control, not a uniform draw. wt_aba runs
   the identical protocol changing only the ligand, so it inherits whatever
   convergence the method has and no independence assumption is needed anywhere.

       delta = P_mandi(true residue) - P_aba(true residue)

   ONLY delta IS EVIDENCE. A high P_mandi that the ABA run matches means the model
   is expressing a ligand-independent preference, exactly the failure the
   ligand-blind oracle exposed in section 23b.

Recall-at-N is still reported, against the uniform chance level (2.11/4 at N=10,
3.84/4 at N=50) -- an upper bound on the null under independence, useful for
bounding, not for deciding.

Effective sample size (distinct sequences, mean pairwise Hamming, mean per-position
entropy) is reported per arm. If 50 sequences carry 3 effective draws that is
itself a finding.

Note the trivial baseline this must beat: ranking WT side chains by steric clash
alone already identifies 3 of the 4 true POSITIONS. Position identification is
nearly free; identity selection is the task.

Usage:  python 48_stage1_run.py [--n 50] [--temperature 0.2] [--dry-run]
"""
import argparse, collections, itertools, json, os, re, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCH = "/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark"
LMPNN = os.path.join(BENCH, "envs", "LigandMPNN")
PY = "/bigdata/cutlerlab/jjaco081/conda_envs/mutpred/bin/python"
MODEL = os.path.join(LMPNN, "model_params", "ligandmpnn_v_32_010_25.pt")
IN = os.path.join(ROOT, "data", "stage1")
OUT = os.path.join(ROOT, "results", "stage1")
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
MPNN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"

TRUE = {59: ("K", "R"), 81: ("V", "I"), 108: ("F", "A"), 159: ("F", "L")}
DESIGN = [59, 81, 83, 92, 94, 108, 110, 120, 122, 141, 159, 160, 163, 164, 167]
DSM_HAO = {
    59: "ADEFGHILMNQRSTVWY", 81: "ILMRTY", 83: "AFGILMSTVWY",
    92: "ADEFGHIKLMNQRSTVWY", 94: "ADFGHIKLMNQRSTVWY",
    108: "ADEGHIKLMNQRSTVWY", 110: "AFGILMSTVWY", 120: "ADEFGHIKLMNQRSTVW",
    122: "ADEFGHIKLMNQRTVWY", 141: "ADEFGHIKLMNQRSTVW",
    159: "ADEGHIKLMNQRSTWY", 160: "DEFGHIKLMNQRSTVWY",
    163: "ADEFGHIKLMNQRSTWY", 164: "ADEFGHIKLMNQRSTWY",
    167: "ADEFGHIKLMQRSTVWY",
}
STRUCTS = ["wt_mandi", "polygly_mandi", "wt_aba"]

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=50)
ap.add_argument("--temperature", type=float, default=0.2)
ap.add_argument("--seed", type=int, default=37)
ap.add_argument("--dry-run", action="store_true")
args = ap.parse_args()
os.makedirs(OUT, exist_ok=True)


def wt_seq(pdb):
    """{resnum: (index_in_sequence, wt_one_letter)} from CA order of chain A."""
    t31 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
           "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
           "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
           "TYR": "Y", "VAL": "V"}
    out, i = {}, 0
    for line in open(pdb):
        if line.startswith("ATOM") and line[12:16] == " CA " and line[21] == "A":
            out[int(line[22:26])] = (i, t31.get(line[17:20].strip(), "X"))
            i += 1
    return out


def omit_json(path, wtmap):
    """Residues to OMIT per position so only DSM-Hao (plus WT) remain."""
    d = {}
    for p in DESIGN:
        if p not in wtmap:
            continue
        allowed = set(DSM_HAO.get(p, AA_ORDER)) | {wtmap[p][1]}
        d[f"A{p}"] = "".join(a for a in AA_ORDER if a not in allowed)
    json.dump(d, open(path, "w"), indent=1)
    return path


def run(struct, alphabet):
    pdb = os.path.join(IN, f"{struct}.pdb")
    tag = f"{struct}__{alphabet}"
    od = os.path.join(OUT, tag)
    os.makedirs(od, exist_ok=True)
    wtmap = wt_seq(pdb)
    cmd = [PY, os.path.join(LMPNN, "run.py"),
           "--model_type", "ligand_mpnn",
           "--checkpoint_ligand_mpnn", MODEL,
           "--pdb_path", pdb, "--out_folder", od,
           "--redesigned_residues", " ".join(f"A{p}" for p in DESIGN if p in wtmap),
           "--number_of_batches", str(args.n), "--batch_size", "1",
           "--temperature", str(args.temperature), "--seed", str(args.seed),
           "--ligand_mpnn_use_side_chain_context", "1",
           "--save_stats", "1"]
    if alphabet == "dsm_hao":
        cmd += ["--omit_AA_per_residue", omit_json(os.path.join(od, "omit.json"), wtmap)]
    if args.dry_run:
        print("  " + " ".join(cmd)); return None, wtmap
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=LMPNN)
    if r.returncode != 0:
        print(f"  !! {tag} failed: {(r.stderr or '').strip().splitlines()[-1:]}")
        return None, wtmap
    return od, wtmap


def load_probs(od):
    """mean sampling_probs over samples -> [L, 21]; None if unavailable."""
    import torch
    st = os.path.join(od, "stats")
    f = [os.path.join(st, x) for x in os.listdir(st)] if os.path.isdir(st) else []
    if not f:
        return None
    d = torch.load(f[0], map_location="cpu", weights_only=False)
    p = d.get("sampling_probs")
    if p is None:
        return None
    p = p.float()
    return (p.mean(0) if p.dim() == 3 else p).numpy()


def load_seqs(od):
    seqs, keep = [], False
    sd = os.path.join(od, "seqs")
    for root, _, files in os.walk(sd):
        for f in files:
            if f.endswith((".fa", ".fasta")):
                for line in open(os.path.join(root, f)):
                    if line.startswith(">"):
                        keep = "id=" in line
                    elif keep:
                        seqs.append(line.strip().split(":")[0])
    return seqs


results = {}
for struct, alphabet in itertools.product(STRUCTS, ["dsm_hao", "free"]):
    tag = f"{struct}__{alphabet}"
    print(f"=== {tag}", flush=True)
    od, wtmap = run(struct, alphabet)
    if args.dry_run or od is None:
        continue
    probs, seqs = load_probs(od), load_seqs(od)
    results[tag] = dict(probs=probs, seqs=seqs, wtmap=wtmap)
    print(f"    {len(seqs)} sequences, probs {'ok' if probs is not None else 'MISSING'}")

if args.dry_run or not results:
    sys.exit(0)

# ---------------------------------------------------------------- report
def pmass(tag, pos, aa):
    r = results.get(tag)
    if not r or r["probs"] is None or pos not in r["wtmap"]:
        return None
    i = r["wtmap"][pos][0]
    j = MPNN_ALPHABET.find(aa)
    return float(r["probs"][i, j]) if 0 <= i < len(r["probs"]) and j >= 0 else None


print("\n" + "=" * 78)
print("PRIMARY: probability mass on the TRUE residue, minus the ABA ligand-swap null")
print("=" * 78)
for alphabet in ("dsm_hao", "free"):
    print(f"\n--- alphabet: {alphabet} ---")
    print(f"{'mut':>8} {'P(wt_mandi)':>12} {'P(polygly)':>11} {'P(wt_aba)':>10} "
          f"{'delta_wt':>9} {'delta_pg':>9}")
    for pos, (w, t) in sorted(TRUE.items()):
        a = pmass(f"wt_mandi__{alphabet}", pos, t)
        b = pmass(f"polygly_mandi__{alphabet}", pos, t)
        c = pmass(f"wt_aba__{alphabet}", pos, t)
        f = lambda v: "   n/a" if v is None else f"{v:6.3f}"
        da = "   n/a" if (a is None or c is None) else f"{a-c:+6.3f}"
        db = "   n/a" if (b is None or c is None) else f"{b-c:+6.3f}"
        print(f"{w}{pos}{t:<5} {f(a):>12} {f(b):>11} {f(c):>10} {da:>9} {db:>9}")
    print("  only delta is evidence; P alone can be a ligand-independent preference")

print("\n" + "=" * 78)
print("SECONDARY: recall-at-N  (uniform chance at N=50 is 3.84/4 -- upper bound)")
print("=" * 78)
print(f"{'arm':>24} {'recalled':>9} {'which':>26}")
for tag, r in results.items():
    got = []
    for pos, (w, t) in TRUE.items():
        if pos not in r["wtmap"]:
            continue
        i = r["wtmap"][pos][0]
        if any(len(s) > i and s[i] == t for s in r["seqs"]):
            got.append(f"{w}{pos}{t}")
    print(f"{tag:>24} {len(got)}/4{'':>5} {','.join(got) if got else '-':>26}")

print("\n" + "=" * 78)
print("EFFECTIVE SAMPLE SIZE -- are N sequences worth N draws?")
print("=" * 78)
print(f"{'arm':>24} {'n':>4} {'distinct':>9} {'mean Hamming':>13} {'mean entropy':>13}")
import math
for tag, r in results.items():
    s = r["seqs"]
    if not s:
        continue
    idx = [r["wtmap"][p][0] for p in DESIGN if p in r["wtmap"]]
    sub = ["".join(x[i] for i in idx if i < len(x)) for x in s]
    dist = len(set(sub))
    pairs = list(itertools.combinations(range(len(sub)), 2))[:2000]
    ham = (sum(sum(c1 != c2 for c1, c2 in zip(sub[i], sub[j])) for i, j in pairs)
           / len(pairs)) if pairs else 0.0
    ent = 0.0
    for k in range(len(idx)):
        c = collections.Counter(x[k] for x in sub if len(x) > k)
        n = sum(c.values())
        ent += -sum((v / n) * math.log2(v / n) for v in c.values()) if n else 0
    print(f"{tag:>24} {len(s):>4} {dist:>9} {ham:>13.2f} {ent/max(len(idx),1):>13.2f}")

json.dump({k: {"seqs": v["seqs"]} for k, v in results.items()},
          open(os.path.join(OUT, "stage1_sequences.json"), "w"))
print(f"\nwrote {OUT}/stage1_sequences.json")
