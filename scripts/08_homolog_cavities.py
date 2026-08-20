#!/usr/bin/env python
"""
08_homolog_cavities.py -- ARM 3: is there a sequence signature for a big pocket?

Question
--------
Across the helix-grip/SRPBCC superfamily, cavity volume spans roughly 400 A3
(PYR/PYL) to >1000 A3 (START domains binding cholesterol/ceramide) in the SAME
fold. If large-cavity members share an identity signature at the PYR1 wall
positions (59/79/94/108/120/141), that signature IS the design prescription --
read off evolution rather than guessed.

WHY THIS SCRIPT DOES NOT USE `foldseek convert2pdb`
---------------------------------------------------
The foldseek CATH50 database stores only C-alpha coordinates. `convert2pdb`
therefore emits CA-ONLY backbone traces (e.g. 148 CA atoms, no side chains),
and a side-chain-free structure has no measurable cavity -- the first version
of this script returned 0.0 A3 for every domain. Additionally `createsubdb`
silently dropped all 105 AlphaFold (`af_*`) entries, extracting only the 45
experimental domains.

This version fetches genuine full-atom structures instead:

  * experimental CATH domains (e.g. `1t17A00`) -> RCSB mmCIF for PDB `1t17`,
    chain `A`. Domain residues are selected by matching the foldseek CA trace
    coordinates against the full structure (exact coordinates, 0.5 A tolerance),
    so only the CATH domain is measured, not the whole chain.
  * AlphaFold entries (e.g. `af_Q7XBY6_2_207_3.30.530.20`) -> AlphaFold DB
    model `AF-Q7XBY6-F1-model_v6.pdb`, residues 2-207 as encoded in the name.
    (v6 is current; v3/v5 return 404.)

Stages
------
  --stage fetch    download structures (network; light on memory, and compute
                   nodes may lack outbound access, so run this on the login
                   node -- curl fits well under the 1 GB cap)
  --stage compute  measure cavities and test the per-position signature.
                   ⚠ The "must run under SLURM" note here was WRONG and cost a
                   needless queue wait. Measured 2026-08-20: peak RSS is **135 MB**
                   and it runs fine on an interactive node, or even under the login
                   node's 1 GB cap. ~30 s per structure, so ~2 h for 266.
                   The 1.15 GB memory warning belongs to the FOLDSEEK search in
                   script 05 and to convert2pdb in the fetch stage - not here.

Cavity method: the same lib_cavity criterion used in Arms 1 and 2, but seeded
UNBIASEDLY -- the largest enclosed component anywhere in the domain, rather
than at a ligand centroid. Volumes are therefore comparable within this script
but only qualitatively comparable to Arm 1.

Usage:
  python 08_homolog_cavities.py --stage fetch    [--min-tm 0.5] [--max-n 150]
  python 08_homolog_cavities.py --stage compute  [--min-tm 0.5] [--max-n 150]
  (esmfold2 env)
"""
import argparse, csv, gzip, json, os, re, subprocess, sys, time, urllib.request
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from lib_cavity import VDW, DIRS, atoms_from_pose_pdb  # noqa: E402

FOLDSEEK = "/bigdata/cutlerlab/jjaco081/tools/foldseek/foldseek/bin/foldseek"
CATH50 = "/bigdata/cutlerlab/jjaco081/tools/foldseek/databases/CATH50"
POSITIONS = [59, 79, 94, 108, 120, 141]
WT = dict(zip(POSITIONS, "KREFYE"))

ap = argparse.ArgumentParser()
ap.add_argument("--stage", choices=["fetch", "compute"], required=True)
ap.add_argument("--min-tm", type=float, default=0.5)
ap.add_argument("--max-n", type=int, default=150)
ap.add_argument("--hits", default=os.path.join(ROOT, "results", "foldseek", "fs_hits.tsv"))
ap.add_argument("--outdir", default=os.path.join(ROOT, "results", "homolog_cavities"))
args = ap.parse_args()

OUT = args.outdir
RAW = os.path.join(OUT, "raw")          # downloaded full-atom structures
DOM = os.path.join(OUT, "domains")      # trimmed, full-atom domain PDBs
CA = os.path.join(OUT, "ca_traces")     # foldseek CA traces (experimental only)
for d in (OUT, RAW, DOM, CA):
    os.makedirs(d, exist_ok=True)


# ---------------------------------------------------------------- hits ------
def load_hits():
    seen, nums = set(), []
    for line in open(os.path.join(ROOT, "data", "pyr1_A.pdb")):
        if line.startswith("ATOM"):
            rn = int(line[22:26])
            if rn not in seen:
                seen.add(rn); nums.append(rn)
    seqidx = {rn: i + 1 for i, rn in enumerate(nums)}
    want = {p: seqidx[p] for p in POSITIONS}

    hits = []
    for line in open(args.hits):
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        tgt, fid, tm = f[1], float(f[2]), float(f[3])
        if tm < args.min_tm:
            continue
        qs, qaln, taln = int(f[5]), f[9], f[10]
        qi = qs - 1
        m = {}
        for qc, tc in zip(qaln, taln):
            if qc != "-":
                qi += 1
                for p, si in want.items():
                    if qi == si:
                        m[p] = tc
        h = dict(target=tgt, fident=fid, alntm=tm)
        h.update({f"aa{p}": m.get(p, "-") for p in POSITIONS})
        hits.append(h)
    return sorted(hits, key=lambda h: -h["alntm"])[:args.max_n]


hits = load_hits()
print(f"{len(hits)} homologs at alnTM >= {args.min_tm}")

AF_RE = re.compile(r"^af_([A-Z0-9]+)_(\d+)_(\d+)_")


def parse_target(t):
    m = AF_RE.match(t)
    if m:
        return dict(kind="af", uniprot=m.group(1),
                    start=int(m.group(2)), end=int(m.group(3)))
    if len(t) >= 5:
        return dict(kind="pdb", pdb=t[:4].lower(), chain=t[4])
    return None


# ---------------------------------------------------------------- fetch -----
def http_get(url, dest, retries=3):
    for k in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "PYR1-pocket-expansion/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            if url.endswith(".gz"):
                data = gzip.decompress(data)
            open(dest, "wb").write(data)
            return True
        except Exception as e:
            if k == retries - 1:
                print(f"      FAIL {url}: {e}")
                return False
            time.sleep(2 * (k + 1))
    return False


if args.stage == "fetch":
    # CA traces for the experimental entries, used to pick domain residues
    exp = [h["target"] for h in hits if parse_target(h["target"])["kind"] == "pdb"]
    if exp:
        nf = os.path.join(OUT, "targets_exp.txt")
        open(nf, "w").write("\n".join(exp) + "\n")
        sub = os.path.join(OUT, "subdb_exp")
        subprocess.run([FOLDSEEK, "createsubdb", nf, CATH50, sub,
                        "--subdb-mode", "1", "--id-mode", "1"],
                       capture_output=True)
        subprocess.run([FOLDSEEK, "convert2pdb", sub, CA,
                        "--pdb-output-mode", "1", "--threads", "4"],
                       capture_output=True)
        print(f"  {len(os.listdir(CA))} CA traces extracted (for domain selection)")

    ok = 0
    for n, h in enumerate(hits, 1):
        t = h["target"]
        info = parse_target(t)
        dest = os.path.join(RAW, t + (".cif" if info["kind"] == "pdb" else ".pdb"))
        if os.path.exists(dest) and os.path.getsize(dest) > 1000:
            ok += 1
            continue
        if info["kind"] == "pdb":
            url = f"https://files.rcsb.org/download/{info['pdb'].upper()}.cif"
        else:
            url = (f"https://alphafold.ebi.ac.uk/files/"
                   f"AF-{info['uniprot']}-F1-model_v6.pdb")
        got = http_get(url, dest)
        ok += bool(got)
        if n % 10 == 0 or not got:
            print(f"  [{n}/{len(hits)}] {t}: {'ok' if got else 'FAILED'}", flush=True)
        time.sleep(0.15)          # be polite to RCSB / EBI
    print(f"\nfetched {ok}/{len(hits)} structures into {RAW}")
    sys.exit(0)


# -------------------------------------------------------------- compute -----
def read_cif_atoms(path, chain):
    """Minimal mmCIF atom_site reader -> list of (resnum, atomname, elem, xyz)."""
    cols, rows, in_loop = {}, [], False
    for line in open(path):
        if line.startswith("_atom_site."):
            cols[line.strip().split(".")[1]] = len(cols)
            in_loop = True
            continue
        if in_loop:
            if line.startswith(("#", "loop_", "_")):
                if rows:
                    break
                continue
            f = line.split()
            if len(f) < len(cols):
                continue
            rows.append(f)
    out = []
    need = ("group_PDB", "label_atom_id", "type_symbol", "auth_seq_id",
            "auth_asym_id", "Cartn_x", "Cartn_y", "Cartn_z", "label_alt_id")
    if not all(k in cols for k in need[:8]):
        return out
    for f in rows:
        if f[cols["group_PDB"]] != "ATOM":
            continue
        if f[cols["auth_asym_id"]] != chain:
            continue
        alt = f[cols["label_alt_id"]] if "label_alt_id" in cols else "."
        if alt not in (".", "?", "A"):
            continue
        e = f[cols["type_symbol"]]
        if e == "H":
            continue
        try:
            rn = int(f[cols["auth_seq_id"]])
        except ValueError:
            continue
        out.append((rn, f[cols["label_atom_id"]].strip('"'), e,
                    np.array([float(f[cols["Cartn_x"]]),
                              float(f[cols["Cartn_y"]]),
                              float(f[cols["Cartn_z"]])])))
    return out


def read_pdb_atoms(path):
    out = []
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        e = (line[76:78].strip() or line[12:16].strip()[0])
        if e == "H":
            continue
        out.append((int(line[22:26]), line[12:16].strip(), e,
                    np.array([float(line[30:38]), float(line[38:46]),
                              float(line[46:54])])))
    return out


def domain_atoms(h):
    """Return full-atom coordinates for the CATH domain only."""
    t = h["target"]
    info = parse_target(t)
    if info["kind"] == "af":
        p = os.path.join(RAW, t + ".pdb")
        if not os.path.exists(p):
            return None
        ats = read_pdb_atoms(p)
        lo, hi = info["start"], info["end"]
        return [a for a in ats if lo <= a[0] <= hi]

    p = os.path.join(RAW, t + ".cif")
    if not os.path.exists(p):
        return None
    ats = read_cif_atoms(p, info["chain"])
    if not ats:
        return None
    trace = os.path.join(CA, t + ".pdb")
    if not os.path.exists(trace):
        return ats                      # whole chain fallback
    tr = np.array([a[3] for a in read_pdb_atoms(trace)])
    if len(tr) == 0:
        return ats
    ca = [(a[0], a[3]) for a in ats if a[1] == "CA"]
    if not ca:
        return ats
    tree = cKDTree(np.array([c[1] for c in ca]))
    d, i = tree.query(tr, k=1)
    keep = {ca[j][0] for j, dist in zip(i, d) if dist < 0.5}
    if len(keep) < 0.5 * len(tr):
        return ats                      # matching failed; use whole chain
    return [a for a in ats if a[0] in keep]


def largest_cavity(xyz, elem, spacing=0.6, probe=1.4, bur_cut=0.88,
                   ray_max=15.0, ray_step=0.75):
    xyz = np.asarray(xyz, float)
    rad = np.array([VDW.get(e, 1.70) for e in elem])
    lo, hi = xyz.min(0) - 3, xyz.max(0) + 3
    axes = [np.arange(lo[i], hi[i], spacing) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    shape = grid.shape[:3]
    pts = grid.reshape(-1, 3)
    tree = cKDTree(xyz)
    d, i = tree.query(pts, k=1)
    clear = d - rad[i]
    free = clear > probe
    occ = (clear < 0).reshape(shape)
    cand = np.where(free)[0]
    if cand.size == 0:
        return 0.0
    cp = pts[cand]
    steps = np.arange(1.0, ray_max, ray_step)
    hits_ = np.zeros(len(cp))
    shp = np.array(shape)
    for dv in DIRS:
        blocked = np.zeros(len(cp), bool)
        for s in steps:
            q = cp + dv * s
            idx = ((q - lo) / spacing).astype(int)
            ok = np.all((idx >= 0) & (idx < shp), axis=1)
            b = np.zeros(len(cp), bool)
            v = idx[ok]
            b[ok] = occ[v[:, 0], v[:, 1], v[:, 2]]
            blocked |= b
        hits_ += blocked
    mask = np.zeros(len(pts), bool)
    mask[cand[(hits_ / len(DIRS)) >= bur_cut]] = True
    lab, n = ndimage.label(mask.reshape(shape))
    if n == 0:
        return 0.0
    sizes = ndimage.sum(mask.reshape(shape), lab, range(1, n + 1))
    return float(sizes.max() * spacing ** 3)


rows = []
_stream = [None, None]   # (handle, DictWriter): stream results as they are made
for n, h in enumerate(hits, 1):
    ats = domain_atoms(h)
    if not ats or len(ats) < 200:
        print(f"  [{n}/{len(hits)}] {h['target']:<34} skipped "
              f"({0 if not ats else len(ats)} atoms)", flush=True)
        continue
    # save the trimmed domain for inspection
    with open(os.path.join(DOM, h["target"] + ".pdb"), "w") as fh:
        for k, (rn, an, e, c) in enumerate(ats, 1):
            fh.write(f"ATOM  {k:5d} {an:<4}{'UNK':>4} A{rn:4d}    "
                     f"{c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}  1.00  0.00"
                     f"          {e:>2}\n")
    vol = largest_cavity([a[3] for a in ats], [a[2] for a in ats])
    r = dict(h); r["cavity_A3"] = round(vol, 1); r["n_atoms"] = len(ats)
    r["n_res"] = len({a[0] for a in ats})
    rows.append(r)
    # Stream to disk immediately. The original wrote the CSV only after the whole
    # loop, so a run killed at 50 of 266 -- by a walltime, a preemption, or an
    # interactive session ending -- lost every measurement it had already made.
    if _stream[0] is None:
        _stream[0] = open(os.path.join(OUT, "homolog_cavities_partial.csv"),
                          "w", newline="")
        _stream[1] = csv.DictWriter(_stream[0], fieldnames=list(r.keys()))
        _stream[1].writeheader()
    _stream[1].writerow(r)
    _stream[0].flush()
    os.fsync(_stream[0].fileno())
    print(f"  [{n}/{len(hits)}] {h['target']:<34} {r['n_res']:>4} res  "
          f"cavity={vol:7.1f} A3   "
          + " ".join(f"{p}:{h['aa'+str(p)]}" for p in POSITIONS), flush=True)

if _stream[0] is not None:
    _stream[0].close()

if not rows:
    print("no domains measured"); sys.exit(1)

csvout = os.path.join(OUT, "homolog_cavities_full.csv")
with open(csvout, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f"\nwrote {csvout} ({len(rows)} domains)")

vols = np.array([r["cavity_A3"] for r in rows])
print(f"\ncavity volume distribution over {len(rows)} helix-grip domains:")
for q in (0, 25, 50, 75, 90, 100):
    print(f"   p{q:<3} {np.percentile(vols, q):8.1f} A3")
big = float(np.percentile(vols, 75))

print("\n=== mean cavity volume by aligned identity at each PYR1 wall position ===")
summary = {}
for p in POSITIONS:
    groups = {}
    for r in rows:
        groups.setdefault(r[f"aa{p}"], []).append(r["cavity_A3"])
    items = [(aa, float(np.mean(v)), len(v),
              sum(1 for x in v if x >= big))
             for aa, v in groups.items() if aa != "-" and len(v) >= 3]
    items.sort(key=lambda x: -x[1])
    summary[str(p)] = items
    print(f"\n  PYR1 {p} (WT {WT[p]}):")
    for aa, mu, k, nbig in items[:8]:
        flag = "  <-- WT" if aa == WT[p] else ""
        print(f"     {aa}  mean={mu:8.1f} A3   n={k:<3} n_top-quartile={nbig}{flag}")

json.dump(summary, open(os.path.join(OUT, "position_signature.json"), "w"),
          indent=2)
print(f"\ntop-quartile threshold: {big:.1f} A3")
