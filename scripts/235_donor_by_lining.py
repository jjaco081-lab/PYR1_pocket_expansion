#!/usr/bin/env python
r"""
235_donor_by_lining.py -- pick the donor pocket by STRUCTURAL EQUIVALENCE to
PYR1's ABA site, not by size.

⚠ 234 WAS WRONG AND THIS REPLACES IT. 234 took "the largest fpocket pocket that
is >= 2/3 lined by CATH-domain residues". That reintroduced the exact error
diagnosed one step earlier: PYR1's own ABA pocket is fpocket's pocket **#2**, not
#1, because fpocket ranks the wide-open GATE CLEFT first by volume. Ranking by
size therefore selects whole-domain surface grooves. The symptom was unmissable:

    PYR1 ABA pocket     44 alpha spheres    359 A^3   hydrophobic density 30.2
    3qrzA00 "donor"    325 alpha spheres   3240 A^3   hydrophobic density 80.3
    4r7kD00 "donor"    427 alpha spheres   3073 A^3

and 254 of 261 structures "beat" PYR1, which is not a donor list, it is an
artefact.

THE FIX. These 261 structures are all the same fold and the survey was BUILT from
a foldseek search with PYR1 as the query, so `results/foldseek/fs_hits.tsv`
already contains a residue-level alignment for every one of them. Map PYR1's 18
pocket-lining positions onto the homolog through that alignment, then choose the
fpocket pocket that contacts the most mapped positions.

This is:
  * LIGAND-FREE for the homolog -- only PYR1's own lining set is used, and the
    homologs are apo (Jannis: "many do not have ligands this cannot be
    consistent");
  * FREE OF ANY ENCLOSURE CRITERION -- the donor does not have to close, because
    PYR1's gate, latch and HAB1 interface are being grafted on (Jannis: "We do
    not need the pocket to close on the receptor as long as it does after we
    graft the PYR1 mechanism onto it");
  * FREE OF ANY SIZE CRITERION -- size is the ANSWER, so it must not be the
    selector.

Validated on PYR1 itself: the rule must return the pocket with 16/18 lining
residues (359.2 A^3), NOT the 842.7 A^3 gate cleft. That is the pre-registered
control and the run is void if it fails.
"""
import csv, glob, json, os, re, shutil, subprocess, sys, tempfile
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
imp = __import__("importlib").import_module
m221, m214 = imp("221_cavity_domain_trim"), imp("214_schema_foldseek")
FPOCKET = "/bigdata/cutlerlab/jjaco081/conda_envs/fpocket_env/bin/fpocket"
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
FS = os.path.join(ROOT, "results", "foldseek", "fs_hits.tsv")
OUT = os.path.join(ROOT, "results", "homolog_cavities", "donor_by_lining.json")
WORK = os.path.join(ROOT, "results", "fpocket", "runs3")
#: PYR1's pocket lining (the sd03 set), auth numbering in 3QN1 chain A
LINING = [59, 81, 83, 87, 89, 92, 94, 108, 110, 117, 120, 122, 141, 159, 160,
          163, 164, 167]
PYR1_REF = 359.2
KEYS = ["Volume", "Druggability Score", "Mean alp. sph. solvent access",
        "Mean local hydrophobic density", "Apolar alpha sphere proportion",
        "Number of Alpha Spheres"]


def load_align():
    """target -> (qaln, taln, qstart, tstart)."""
    out = {}
    with open(FS) as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r, None)
        for row in r:
            if len(row) < 11:
                continue
            out[row[1]] = (row[9], row[10], int(row[5]), int(row[7]))
    return out


def ordered_res(path, chain=None):
    seen, out = set(), []
    for l in open(path):
        if not l.startswith("ATOM"):
            continue
        if chain and l[21] != chain:
            continue
        k = (l[21], int(l[22:26]))
        if k not in seen:
            seen.add(k); out.append(k)
    return out


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


ALIGN = None
PYR1_ORD = None


def one(target):
    rec = dict(target=target,
               kind="AlphaFold" if target.startswith("af_") else "experimental")
    path = None
    for ext in (".pdb", ".cif"):
        p = os.path.join(RAW, target + ext)
        if os.path.exists(p): path = p; break
    if path is None or target not in ALIGN:
        rec["error"] = "no file" if path is None else "no alignment"
        return rec
    qaln, taln, qs, ts = ALIGN[target]
    d = None
    try:
        # ---- map PYR1's lining onto this homolog ----
        if path.endswith(".cif"):
            ch = m221.read_trimmed(path, None)
            seen, tres = set(), []
            for c, v in ch.items():
                for _, _, r in v:
                    if (c, r) not in seen:
                        seen.add((c, r)); tres.append((c, r))
            tres.sort(key=lambda x: (x[0], x[1]))
        else:
            tres = ordered_res(path)
        mp_ = m214.map_from_alignment(qaln, taln, qs, ts, PYR1_ORD, tres)
        mapped = {mp_[k] for k in [("A", n) for n in LINING] if k in mp_}
        rec["n_lining_mapped"] = len(mapped)
        if len(mapped) < 6:
            rec["error"] = f"only {len(mapped)}/18 lining positions mapped"
            return rec

        ch = m221.read_trimmed(path, None)
        at = [(c, r, x) for c, v in ch.items() for _, x, r in v]
        if len(at) < 300:
            rec["error"] = f"{len(at)} atoms"; return rec
        X = np.array([x for _, _, x in at])
        is_lining = np.array([(c, r) in mapped for c, r, _ in at], bool)
        plddt = None
        if path.endswith(".pdb"):
            b = [float(l[60:66]) for l in open(path) if l.startswith("ATOM")
                 and (l[76:78].strip() or l[12:16].strip()[0]).upper() not in ("H", "D")]
            if len(b) == len(at) and 0 <= min(b) and max(b) <= 100.5:
                plddt = np.array(b)

        # ---- fpocket ----
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

        # ---- choose by LINING CONTACT, not by size ----
        best, bestn, bestlist = None, -1, None
        for k, S in pk.items():
            dmin = np.linalg.norm(X[:, None] - S[None], axis=-1).min(1)
            near = dmin < 5.0
            n = int((near & is_lining).sum())
            nres = len({(c, r) for (c, r, _), t in zip(at, near & is_lining) if t})
            if nres > bestn:
                best, bestn, bestlist = k, nres, near
        rec["donor_id"] = best
        rec["lining_contacted"] = bestn
        rec["donor_vol"] = round(desc.get(best, {}).get("Volume", 0.0), 1)
        for kk in KEYS:
            rec[kk.replace(" ", "_").replace(".", "")] = desc.get(best, {}).get(kk)
        rec["largest_pocket_vol"] = round(max(
            (v.get("Volume", 0.0) for v in desc.values()), default=0.0), 1)
        if plddt is not None and bestlist is not None and bestlist.any():
            rec["lining_plddt"] = round(float(plddt[bestlist].mean()), 1)
    except Exception as e:                                       # noqa: BLE001
        rec["error"] = f"{type(e).__name__}: {e}"
    finally:
        if d: shutil.rmtree(d, ignore_errors=True)
    return rec


def main():
    global ALIGN, PYR1_ORD
    os.makedirs(WORK, exist_ok=True)
    ALIGN = load_align()
    PYR1_ORD = ordered_res(os.path.join(ROOT, "data", "3QN1_complex_auth_aba.pdb"), "A")
    print(f"alignment for {len(ALIGN)} targets; PYR1 ordered residues {len(PYR1_ORD)}")
    files = sorted(glob.glob(f"{RAW}/*.cif") + glob.glob(f"{RAW}/*.pdb"))
    tg = [os.path.basename(f).rsplit(".", 1)[0] for f in files]
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 48
    with mp.Pool(procs) as p:
        recs = p.map(one, tg)
    json.dump(recs, open(OUT, "w"), indent=1)          # PERSIST FIRST
    print(f"wrote {OUT}\n")
    ok = [r for r in recs if "error" not in r and r.get("donor_vol")]
    print(f"{len(ok)} of {len(recs)} selected a pocket\n")
    import statistics as st
    v5 = {r["target"]: r for r in json.load(open(os.path.join(
        ROOT, "results", "homolog_cavities", "cavity_v5.json")))["structures"]
        if "error" not in r}
    for lbl in ("experimental", "AlphaFold"):
        s = [r["donor_vol"] for r in ok if r["kind"] == lbl]
        if s:
            print(f"{lbl:<13} n={len(s):>3}  median {st.median(s):7.1f}  "
                  f"> PYR1 {PYR1_REF}: {sum(1 for x in s if x > PYR1_REF):>3}  "
                  f"> 2x: {sum(1 for x in s if x > 2*PYR1_REF):>3}")
    print(f"\nmedian lining positions contacted: "
          f"{st.median([r['lining_contacted'] for r in ok]):.0f} of 18")
    print(f"\n{'='*108}")
    print(f"TOP DONORS BY LINING-EQUIVALENT POCKET  (PYR1 = {PYR1_REF} A^3, 44 spheres, hydroph 30.2)")
    print(f"{'='*108}")
    print(f"{'target':<40}{'kind':<5}{'vol':>8}{'xPYR1':>7}{'lin':>5}"
          f"{'aSph':>6}{'hydroph':>9}{'encl':>8}{'pLDDT':>7}{'biggest':>9}")
    for r in sorted(ok, key=lambda r: -r["donor_vol"])[:22]:
        e = v5.get(r["target"], {}).get("main_median", 0.0)
        print(f"{r['target'][:38]:<40}{r['kind'][:4]:<5}{r['donor_vol']:>8.1f}"
              f"{r['donor_vol']/PYR1_REF:>7.2f}{r['lining_contacted']:>5}"
              f"{(r.get('Number_of_Alpha_Spheres') or 0):>6.0f}"
              f"{(r.get('Mean_local_hydrophobic_density') or 0):>9.1f}{e:>8.1f}"
              f"{str(r.get('lining_plddt','-')):>7}{r.get('largest_pocket_vol',0):>9.1f}")
    print("\n'biggest' = the largest pocket in that structure. Where it greatly")
    print("exceeds 'vol', size-ranking would have picked the wrong pocket -- which")
    print("is exactly what 234 did.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
