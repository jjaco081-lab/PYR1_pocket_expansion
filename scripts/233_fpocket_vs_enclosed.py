#!/usr/bin/env python
r"""
233_fpocket_vs_enclosed.py -- measure every homolog pocket TWICE, with an
established tool and with ours, on the SAME physical site.

WHY. Jannis: "No since many do not have ligands this cannot be consistent."
Correct, and it killed my proposal to merge fpocket pockets by ligand occupancy —
most homologs are apo. He also asked whether we had run fpocket correctly,
because splitting ABA into two roughly equal pockets with no visible bottleneck
looked wrong. Both questions resolved:

 * fpocket's two "pockets on ABA" are 2.50 A apart at their closest alpha-sphere
   centres, against a default single-linkage clustering distance of 1.73 A -- so
   they are split by a THRESHOLD, not a bottleneck.
 * but they are not two halves of one pocket. Their linings say what they are:
       pocket with 16/18 sd03 lining residues, chain A only   = THE ABA POCKET
       pocket with 6 gate + 4 latch residues, spans A+B        = THE GATE CLEFT
   In the apo monomer the gate cleft is WIDE OPEN (842.7 A^3) and fpocket ranked
   it first by volume; with HAB1 bound it collapses to 326.7 while the ligand
   pocket is unchanged (359.2 -> 368.1). fpocket was right; my RANKING was wrong.

THE DEFINITIONAL GAP, which is the point of this script:
    fpocket finds SURFACE CLEFTS  -- alpha-sphere hulls that extend past the
                                     buried void
    our code finds ENCLOSED CAVITIES -- the burial criterion discards open cleft
                                     entirely
They answer different questions. For a receptor that must CLOSE around its
ligand, enclosed volume is the relevant one; for "is there a big pocket here at
all", fpocket's is standard and literature-comparable. Reporting both makes
their agreement a quality check, and their DISAGREEMENT flags the structures
where the pocket is an open groove rather than a sealed chamber -- which is
exactly the property that disqualifies a graft donor for a closed-state readout.

HOW THE SAME SITE IS IDENTIFIED WITHOUT A LIGAND AND WITHOUT AN ALIGNMENT:
each fpocket pocket is paired to one of our chambers by SPATIAL OVERLAP (fraction
of the pocket's alpha-sphere centres falling within 3 A of a chamber voxel). The
reported pair is our domain-lined chamber and whichever fpocket pocket overlaps
it most. No ligand, no residue numbering, no superposition.

⚠ Volumes are NOT interchangeable. fpocket's alpha-sphere volume for PYR1's ABA
pocket is 359.2 A^3 where our enclosed volume is 164.4. A ratio near 2 is normal;
the interesting cases are the EXTREMES, and those are what this prints.
"""
import glob, json, os, re, subprocess, sys, tempfile, shutil
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m151 = import_module("151_cavity_bottleneck")
m221 = import_module("221_cavity_domain_trim")
m229 = import_module("229_cavity_visualise")
FPOCKET = "/bigdata/cutlerlab/jjaco081/conda_envs/fpocket_env/bin/fpocket"
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
OUT = os.path.join(ROOT, "results", "homolog_cavities", "fpocket_vs_enclosed.json")
WORK = os.path.join(ROOT, "results", "fpocket", "runs")
OVERLAP = 3.0
DOMAIN_FRAC = 2.0 / 3.0


def run_one(target):
    rec = dict(target=target)
    path = None
    for ext in (".pdb", ".cif"):
        p = os.path.join(RAW, target + ext)
        if os.path.exists(p):
            path = p; break
    if path is None:
        rec["error"] = "no file"; return rec
    rng = m221.domain_range(target)
    rec["kind"] = "AlphaFold" if target.startswith("af_") else "experimental"
    rec["domain"] = list(rng) if rng else None
    try:
        # ---- our atoms (first model, whitelist, no H) ----
        ch = m221.read_trimmed(path, None)
        at = [a for v in ch.values() for a in v]
        if len(at) < 300:
            rec["error"] = f"{len(at)} atoms"; return rec
        X = np.array([x for _, x, _ in at])
        el = [e for e, _, _ in at]
        dom = np.array([rng is None or (rng[0] <= r <= rng[1])
                        for _, _, r in at], bool)
        plddt = None
        if path.endswith(".pdb"):
            b = [float(l[60:66]) for l in open(path)
                 if l.startswith("ATOM")
                 and (l[76:78].strip() or l[12:16].strip()[0]).upper() not in ("H", "D")]
            if len(b) == len(at) and 0 <= min(b) and max(b) <= 100.5:
                plddt = np.array(b)

        # ---- our enclosed chambers ----
        cs = m229.chambers_with_voxels(X, el, dom)
        dl = [c for c in cs if c["frac"] >= DOMAIN_FRAC]
        ours = dl[0] if dl else None
        rec["enclosed"] = round(ours["vol"], 1) if ours else 0.0
        rec["enclosed_frac"] = round(ours["frac"], 3) if ours else None
        rec["n_chamber"] = len(cs)

        # ---- fpocket on the same atoms ----
        d = tempfile.mkdtemp(dir=WORK)
        pdb = os.path.join(d, "s.pdb")
        with open(pdb, "w") as fh:
            for i, (e, x, r) in enumerate(at, 1):
                fh.write(f"ATOM  {i%99999:5d}  CA  ALA A{r%9999:4d}    "
                         f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00"
                         f"          {e:>2s}\n")
            fh.write("END\n")
        # ⚠ fpocket needs real atom names to type atoms. Feed the original file
        # when it is a PDB; only fall back to the rewrite for CIF inputs.
        src = path if path.endswith(".pdb") else pdb
        if path.endswith(".pdb"):
            shutil.copyfile(path, pdb)
            src = pdb
        subprocess.run([FPOCKET, "-f", src], cwd=d, timeout=600,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        base = os.path.splitext(os.path.basename(src))[0]
        pqr = os.path.join(d, base + "_out", base + "_pockets.pqr")
        info = os.path.join(d, base + "_out", base + "_info.txt")
        pk, vols = {}, {}
        if os.path.exists(pqr):
            tmp = {}
            for l in open(pqr):
                if l.startswith("ATOM"):
                    tmp.setdefault(int(l[22:26]), []).append(
                        [float(l[30:38]), float(l[38:46]), float(l[46:54])])
            pk = {k: np.array(v) for k, v in tmp.items()}
        if os.path.exists(info):
            cur = None
            for l in open(info):
                m = re.match(r"Pocket (\d+)", l)
                if m: cur = int(m.group(1))
                if "Volume :" in l and cur: vols[cur] = float(l.split(":")[1])
        rec["n_fpocket"] = len(pk)

        # ---- pair by spatial overlap with OUR chamber ----
        if ours is not None and pk:
            cav = ours["pts"]
            best, bestf = None, 0.0
            for k, S in pk.items():
                dmin = np.linalg.norm(S[:, None] - cav[None], axis=-1).min(1)
                f = float((dmin < OVERLAP).mean())
                if f > bestf:
                    best, bestf = k, f
            rec["fpocket_id"] = best
            rec["fpocket_vol"] = round(vols.get(best, 0.0), 1)
            rec["overlap"] = round(bestf, 3)
            rec["fpocket_top_vol"] = round(max(vols.values()), 1) if vols else 0.0
            if plddt is not None:
                nb = np.linalg.norm(X[:, None] - cav[None][:, ::7], axis=-1).min(1) \
                     if cav.shape[0] > 7 else np.linalg.norm(X[:, None] - cav[None], axis=-1).min(1)
                near = nb < 5.0
                rec["lining_plddt"] = round(float(plddt[near].mean()), 1) if near.any() else None
        shutil.rmtree(d, ignore_errors=True)
    except Exception as e:                                       # noqa: BLE001
        rec["error"] = f"{type(e).__name__}: {e}"
    return rec


def main():
    os.makedirs(WORK, exist_ok=True)
    files = sorted(glob.glob(f"{RAW}/*.cif") + glob.glob(f"{RAW}/*.pdb"))
    targets = [os.path.basename(f).rsplit(".", 1)[0] for f in files]
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 48
    print(f"{len(targets)} structures, {procs} procs", flush=True)
    with mp.Pool(procs) as p:
        recs = p.map(run_one, targets)
    json.dump(recs, open(OUT, "w"), indent=1)          # PERSIST FIRST
    print(f"wrote {OUT}")
    ok = [r for r in recs if "error" not in r and r.get("fpocket_vol")
          and r.get("enclosed")]
    print(f"{len(ok)} of {len(recs)} have both measures\n")
    import statistics as st
    ratio = [(r["fpocket_vol"] / r["enclosed"], r) for r in ok]
    rs = [x for x, _ in ratio]
    print(f"fpocket / enclosed ratio: median {st.median(rs):.2f}  "
          f"IQR {np.percentile(rs,25):.2f}-{np.percentile(rs,75):.2f}")
    print(f"  (PYR1 reference: 359.2 / 164.4 = 2.18)\n")
    print("=" * 84)
    print("LARGEST DISAGREEMENTS -- fpocket >> enclosed  (OPEN GROOVE, not a sealed chamber)")
    print("=" * 84)
    print(f"{'target':<42}{'fpocket':>9}{'enclosed':>10}{'ratio':>7}{'ovl':>6}{'pLDDT':>7}")
    for x, r in sorted(ratio, reverse=True)[:10]:
        print(f"{r['target'][:40]:<42}{r['fpocket_vol']:>9.1f}{r['enclosed']:>10.1f}"
              f"{x:>7.1f}{r.get('overlap',0):>6.2f}{str(r.get('lining_plddt','-')):>7}")
    print("\n" + "=" * 84)
    print("fpocket << enclosed  (we find a sealed void fpocket does not call a pocket)")
    print("=" * 84)
    print(f"{'target':<42}{'fpocket':>9}{'enclosed':>10}{'ratio':>7}{'ovl':>6}{'pLDDT':>7}")
    for x, r in sorted(ratio)[:10]:
        print(f"{r['target'][:40]:<42}{r['fpocket_vol']:>9.1f}{r['enclosed']:>10.1f}"
              f"{x:>7.1f}{r.get('overlap',0):>6.2f}{str(r.get('lining_plddt','-')):>7}")
    print("\n" + "=" * 84)
    print("BEST AGREEMENT (ratio nearest PYR1's 2.18) -- the trustworthy ones")
    print("=" * 84)
    for x, r in sorted(ratio, key=lambda t: abs(t[0] - 2.18))[:8]:
        print(f"{r['target'][:40]:<42}{r['fpocket_vol']:>9.1f}{r['enclosed']:>10.1f}"
              f"{x:>7.2f}{r.get('overlap',0):>6.2f}{str(r.get('lining_plddt','-')):>7}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
