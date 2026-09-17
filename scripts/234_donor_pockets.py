#!/usr/bin/env python
r"""
234_donor_pockets.py -- rank DONOR pockets, on the criterion that actually
applies.

⚠ THIS CORRECTS THE PREMISE OF 233. Jannis: "We do not need the pocket to close
on the receptor as long as it does after we graft the PYR1 mechanism onto it. We
are mainly just looking for good donor pockets right now to expand the pocket."

233 paired each fpocket pocket to OUR ENCLOSED CHAMBER and treated a large
open groove as disqualifying. Both are wrong for donor screening:
 * the donor does not have to close -- PYR1's gate, latch and HAB1 interface are
   what close the chimera, and they are being grafted ON;
 * pairing by overlap with an ENCLOSED chamber biases selection toward pockets
   that are already sealed, which is exactly backwards.

So the donor metric is the fpocket VOLUME of the pocket that belongs to the
SRPBCC domain, and our enclosed volume is demoted to a descriptor.

SELECTION, with no enclosed-volume bias and no ligand:
    the largest fpocket pocket whose alpha spheres are predominantly (>= 2/3)
    lined by CATH-DOMAIN residues.
The domain restriction survives from v4 because it fixed a real error -- 26 % of
the survey was reporting a chamber belonging to a different domain -- but nothing
about enclosure enters the choice.

THE BAR, from PYR1 measured the same way:
    ABA pocket   apo 359.2 A^3 / +HAB1 368.1    44 alpha spheres
                 solvent access 0.47   hydrophobic density 30.2   apolar 0.73
    gate cleft   apo 842.7                      93 alpha spheres
                 solvent access 0.49   hydrophobic density 48.3   apolar 0.87
⚠ Solvent access barely distinguishes the ABA pocket from the open gate cleft
(0.47 vs 0.49), so it is NOT a useful buriedness filter here. Volume and
hydrophobic density do separate them (359 vs 843, 30.2 vs 48.3), so both are
reported and the gate cleft is the worked example of what a large-but-wrong
pocket looks like.

A donor is interesting if its domain pocket is BIGGER than PYR1's 359 A^3, since
the whole point is to expand. Reported with enclosed volume and lining pLDDT
alongside so an inflated or badly-predicted one is visible, not hidden.
"""
import glob, json, os, re, shutil, subprocess, sys, tempfile
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
m221 = __import__("importlib").import_module("221_cavity_domain_trim")
m229 = __import__("importlib").import_module("229_cavity_visualise")
FPOCKET = "/bigdata/cutlerlab/jjaco081/conda_envs/fpocket_env/bin/fpocket"
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
OUT = os.path.join(ROOT, "results", "homolog_cavities", "donor_pockets.json")
WORK = os.path.join(ROOT, "results", "fpocket", "runs2")
DOMFRAC = 2.0 / 3.0
PYR1_VOL = 359.2
KEYS = ["Volume", "Druggability Score", "Mean alp. sph. solvent access",
        "Mean local hydrophobic density", "Apolar alpha sphere proportion",
        "Number of Alpha Spheres"]


def parse_info(f):
    out, cur = {}, None
    for l in open(f):
        m = re.match(r"Pocket (\d+)", l)
        if m:
            cur = int(m.group(1)); out[cur] = {}
        elif cur and ":" in l:
            k, v = l.rsplit(":", 1)
            try: out[cur][k.strip()] = float(v)
            except ValueError: pass
    return out


def one(target):
    rec = dict(target=target,
               kind="AlphaFold" if target.startswith("af_") else "experimental")
    path = None
    for ext in (".pdb", ".cif"):
        p = os.path.join(RAW, target + ext)
        if os.path.exists(p): path = p; break
    if path is None:
        rec["error"] = "no file"; return rec
    rng = m221.domain_range(target)
    rec["domain"] = list(rng) if rng else None
    d = None
    try:
        ch = m221.read_trimmed(path, None)
        at = [a for v in ch.values() for a in v]
        if len(at) < 300:
            rec["error"] = f"{len(at)} atoms"; return rec
        X = np.array([x for _, x, _ in at])
        dom = np.array([rng is None or (rng[0] <= r <= rng[1])
                        for _, _, r in at], bool)
        plddt = None
        if path.endswith(".pdb"):
            b = [float(l[60:66]) for l in open(path) if l.startswith("ATOM")
                 and (l[76:78].strip() or l[12:16].strip()[0]).upper() not in ("H", "D")]
            if len(b) == len(at) and 0 <= min(b) and max(b) <= 100.5:
                plddt = np.array(b)

        d = tempfile.mkdtemp(dir=WORK)
        src = os.path.join(d, "s" + os.path.splitext(path)[1])
        shutil.copyfile(path, src)
        subprocess.run([FPOCKET, "-f", src], cwd=d, timeout=900,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        base = os.path.splitext(os.path.basename(src))[0]
        od = os.path.join(d, base + "_out")
        pqr, info = os.path.join(od, base + "_pockets.pqr"), os.path.join(od, base + "_info.txt")
        if not (os.path.exists(pqr) and os.path.exists(info)):
            rec["error"] = "fpocket produced nothing"; return rec
        tmp = {}
        for l in open(pqr):
            if l.startswith("ATOM"):
                tmp.setdefault(int(l[22:26]), []).append(
                    [float(l[30:38]), float(l[38:46]), float(l[46:54])])
        pk = {k: np.array(v) for k, v in tmp.items()}
        desc = parse_info(info)
        rec["n_fpocket"] = len(pk)

        # ---- domain-lined fraction of each pocket, NO enclosure anywhere ----
        cands = []
        for k, S in pk.items():
            dmin = np.linalg.norm(X[:, None] - S[None], axis=-1).min(1)
            near = dmin < 5.0
            if not near.any(): continue
            f = float(dom[near].sum() / near.sum())
            cands.append((desc.get(k, {}).get("Volume", 0.0), f, k, near))
        rec["n_domain_pockets"] = sum(1 for v, f, k, _ in cands if f >= DOMFRAC)
        dl = [c for c in cands if c[1] >= DOMFRAC]
        if not dl:
            rec["donor_vol"] = 0.0; return rec
        vol, frac, pid, near = max(dl, key=lambda c: c[0])
        rec.update(donor_vol=round(vol, 1), donor_frac=round(frac, 3),
                   donor_id=pid)
        for kk in KEYS:
            rec[kk.replace(" ", "_").replace(".", "")] = desc.get(pid, {}).get(kk)
        if plddt is not None:
            rec["lining_plddt"] = round(float(plddt[near].mean()), 1)
        rec["largest_any_vol"] = round(max(v for v, _, _, _ in cands), 1)
    except Exception as e:                                       # noqa: BLE001
        rec["error"] = f"{type(e).__name__}: {e}"
    finally:
        if d: shutil.rmtree(d, ignore_errors=True)
    return rec


def main():
    os.makedirs(WORK, exist_ok=True)
    files = sorted(glob.glob(f"{RAW}/*.cif") + glob.glob(f"{RAW}/*.pdb"))
    tg = [os.path.basename(f).rsplit(".", 1)[0] for f in files]
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 48
    print(f"{len(tg)} structures, {procs} procs", flush=True)
    with mp.Pool(procs) as p:
        recs = p.map(one, tg)
    json.dump(recs, open(OUT, "w"), indent=1)        # PERSIST FIRST
    print(f"wrote {OUT}\n")
    ok = [r for r in recs if "error" not in r and r.get("donor_vol")]
    print(f"{len(ok)} of {len(recs)} have a domain-lined fpocket pocket\n")
    v5 = {r["target"]: r for r in json.load(open(os.path.join(
        ROOT, "results", "homolog_cavities", "cavity_v5.json")))["structures"]
        if "error" not in r}
    import statistics as st
    for lbl in ("experimental", "AlphaFold"):
        s = [r["donor_vol"] for r in ok if r["kind"] == lbl]
        if s:
            print(f"{lbl:<13} n={len(s):>3}  median {st.median(s):7.1f}  "
                  f"> PYR1's {PYR1_VOL}: {sum(1 for x in s if x > PYR1_VOL):>3}  "
                  f"> 2x PYR1: {sum(1 for x in s if x > 2*PYR1_VOL):>3}")
    print(f"\n{'='*104}")
    print(f"TOP DONOR POCKETS  (fpocket volume of the domain-lined pocket; "
          f"PYR1's ABA pocket = {PYR1_VOL} A^3)")
    print(f"{'='*104}")
    print(f"{'target':<40}{'kind':<5}{'fpocket':>9}{'x PYR1':>8}{'encl':>8}"
          f"{'aSph':>6}{'solv':>6}{'hydroph':>9}{'pLDDT':>7}")
    for r in sorted(ok, key=lambda r: -r["donor_vol"])[:20]:
        e = v5.get(r["target"], {}).get("main_median", 0.0)
        print(f"{r['target'][:38]:<40}{r['kind'][:4]:<5}{r['donor_vol']:>9.1f}"
              f"{r['donor_vol']/PYR1_VOL:>8.2f}{e:>8.1f}"
              f"{(r.get('Number_of_Alpha_Spheres') or 0):>6.0f}"
              f"{(r.get('Mean_alp_sph_solvent_access') or 0):>6.2f}"
              f"{(r.get('Mean_local_hydrophobic_density') or 0):>9.1f}"
              f"{str(r.get('lining_plddt','-')):>7}")
    print("\n⚠ PYR1's own GATE CLEFT measures 842.7 A^3 with hydrophobic density")
    print("  48.3 -- a large pocket that is NOT the ligand site. Any donor whose")
    print("  hydrophobic density is far above PYR1's pocket value of 30.2 may be")
    print("  the same kind of thing and should be looked at, not just ranked.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
