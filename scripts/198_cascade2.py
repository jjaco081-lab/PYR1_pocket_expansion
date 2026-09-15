#!/usr/bin/env python
r"""
198_cascade2.py -- the CORRECTED structural cascade, v2.

Three defects in v1 (192/195) were exposed when Jannis eye-checked arm0b
d000/d002/d009 and got the OPPOSITE ranking to my filters. All three are fixed
here, and all three were the same underlying mistake: measuring a proxy instead
of the thing.

  1. CONNECTIVITY IS NOW FILTER ZERO.
     v1 gated on the NUMBER of chain breaks. d000 failed on "3 breaks" that were
     the seams of a 16-residue orphan floating 59 A away -- droppable, and its
     169-residue main body is intact. d009 passed with "2 breaks" that actually
     split the protein into 112 + 77 residue halves 21 A apart, with scaffolded
     motif on BOTH sides. Break count cannot tell those apart; connectivity can.

  2. CAVITY AND ENVELOPE ARE MEASURED ON THE LARGEST CONNECTED PIECE ONLY.
     d009's cavity was 80% in the gap BETWEEN its two halves: 61 A^3 with both
     pieces, 1 A^3 with the main piece alone, and a 1354 A^3 "envelope" that was
     measuring the gap. This is the same class of error as the disconnected-void
     bug Jannis caught in 145 -- fixed there for volumes, never for the CHAIN.

  3. HAB1 CONTACT IS WHOLE-CHAIN, GATED ON THE LONGEST CONTIGUOUS RUN.
     v1's F4 counted contacts in the LEADING 25 residues only, because that was
     batch 1-3's defect. It reported 0 for both d002 and d009 while d009 had an
     18-residue segment wrapping HAB1 on an inter-motif linker. The appendage is
     what the eye tracks, so the run length is what gets gated.

Calibration of the run-length threshold is honest but thin: d000 (longest run 2)
reads clean to Jannis, d002 (5) and d009 (18) read as appendages. n = 3. The raw
counts are therefore written to the CSV so the cut can be moved WITHOUT
recomputing anything.

⚠ F4 is vacuous on wild-type PYR1 by construction -- PYR1 has no designed
residues, so it has no non-motif segment to wrap HAB1 with. PYR1's own extensive
HAB1 interface must not be counted against it, which is why the filter is
defined over DESIGNED residues only.

Usage:
  python 198_cascade2.py --arms results/rfd3/arm1slack_s3001 [...] --procs 64
  python 198_cascade2.py --calibrate          # wild-type PYR1 hard gate
"""
import argparse, csv, glob, json, os, sys
import multiprocessing as mp
import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m143 = import_module("143_pocket_backbone_shape")
m151 = import_module("151_cavity_bottleneck")
m192 = import_module("192_filter_cascade")
m195 = import_module("195_arm_compare")

CONTACT = 4.5
BREAK = 4.5
PYR1_ENVELOPE, PYR1_RMAX, PYR1_RG, PYR1_CHAMBER = 1184.0, 3.21, 15.0, 119.0
#: orphan fragments up to this many residues, at a terminus, may be truncated
TRUNCATABLE = 25
#: longest contiguous run of DESIGNED residues touching HAB1 that is allowed
MAX_APPENDAGE = 4

FIELDS = ["arm", "design", "n_pieces", "n_res", "n_main", "n_orphan",
          "orphan_termini_only", "motif_split", "cavity", "r_max", "envelope",
          "ncomp", "rg", "dev_max", "nonmotif_hab1", "appendage", "verdict",
          "fail_at"]


def pieces_of(ca):
    """[[resnums]] -- split the chain wherever consecutive CAs exceed BREAK."""
    res = sorted(ca)
    out = [[res[0]]]
    for i in range(len(res) - 1):
        if res[i + 1] == res[i] + 1 and \
                np.linalg.norm(ca[res[i + 1]] - ca[res[i]]) > BREAK:
            out.append([])
        out[-1].append(res[i + 1])
    return out


def run_lengths(touching, motif_out, res):
    """longest contiguous run of DESIGNED (non-motif) residues touching HAB1."""
    best = cur = 0
    for r in res:
        if r in touching and r not in motif_out:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def measure_one(args):
    arm, cif = args
    name = os.path.basename(os.path.dirname(cif))
    if not name.startswith("d"):
        name = "m" + os.path.basename(cif).split("_model_")[1].split(".")[0]
    row = dict.fromkeys(FIELDS, "")
    row.update(arm=os.path.basename(arm), design=name)
    try:
        js = cif.replace(".cif.gz", ".json").replace(".cif", ".json")
        at = m195.read_cif_gz(cif)
        imap = json.load(open(js)).get("diffused_index_map", {})
        dc = m195.design_chain(at, imap)
        dat = [a for a in at if a[0] == dc]
        hab1 = np.array([a[4] for a in at if a[0] != dc])
        ca, by = {}, {}
        for a in dat:
            by.setdefault(a[1], []).append(a[4])
            if a[2] == "CA":
                ca[a[1]] = a[4]
        res = sorted(ca)
        motif_out = {int(v[1:]) for v in imap.values() if v[0] == dc}

        # ---- F0 connectivity ------------------------------------------------
        ps = pieces_of(ca)
        main = max(ps, key=len)
        orphan = [r for r in res if r not in set(main)]
        term_only = all(r < main[0] or r > main[-1] for r in orphan)
        msplit = bool(motif_out - set(main))
        row.update(n_pieces=len(ps), n_res=len(res), n_main=len(main),
                   n_orphan=len(orphan), orphan_termini_only=int(term_only),
                   motif_split=int(msplit))
        if msplit:
            row.update(verdict="REJECT", fail_at="F0 motif split across pieces")
            return row
        if orphan and not (term_only and len(orphan) <= TRUNCATABLE):
            row.update(verdict="REJECT", fail_at="F0 internal/large orphan")
            return row

        # ---- everything downstream uses the MAIN PIECE ONLY -----------------
        keep = set(main)
        mat = [(a[1], None, a[2], a[3], a[4]) for a in dat if a[1] in keep]
        r, wall, pts = m192.cavity_and_wall(mat)
        env = m192.envelope(mat, wall) if wall else None
        mca = np.array([ca[x] for x in main])
        rg = float(np.sqrt(((mca - mca.mean(0)) ** 2).sum(1).mean()))
        row.update(cavity=round((r or {}).get("main", 0.0), 1),
                   r_max=round((r or {}).get("r_max", 0.0), 2),
                   envelope=round((env or {}).get("envelope", 0.0), 1),
                   ncomp=(env or {}).get("ncomp", 0), rg=round(rg, 1))

        # ---- motif identity + deviation -------------------------------------
        comp = {a[1]: a[5] for a in dat}
        ref = m195.ref_motif_ca()
        P, Q, bad = [], [], []
        for ch, lo, hi, seq in m195.MOTIF:
            for k, resi in enumerate(range(lo, hi + 1)):
                key = f"{ch}{resi}"
                if key not in imap:
                    continue                      # segment absent from this arm
                oi = int(imap[key][1:])
                if oi not in keep:
                    continue
                if m195.AA3.get(comp.get(oi, "?"), "?") != seq[k]:
                    bad.append(key); continue
                rk = ref.get((ch, resi))
                if rk is not None and oi in ca:
                    P.append(ca[oi]); Q.append(rk[0])
        if bad:
            row.update(verdict="REJECT", fail_at=f"motif identity ({len(bad)})")
            return row
        dev = m195.kabsch_rmsd(np.array(P), np.array(Q))[0] if len(P) > 3 else None
        row["dev_max"] = "" if dev is None else round(dev, 2)

        # ---- F4 whole-chain appendage, DESIGNED residues only ---------------
        touch = set()
        if len(hab1):
            for rr in main:
                xyz = np.array(by[rr])
                if np.linalg.norm(hab1[None, :, :] - xyz[:, None, :],
                                  axis=2).min() < CONTACT:
                    touch.add(rr)
        app = run_lengths(touch, motif_out, main)
        row.update(nonmotif_hab1=len(touch - motif_out), appendage=app)

        # ---- verdict ---------------------------------------------------------
        for ok, why in [
                (r is not None and r.get("main", 0) > 0, "F1 no enclosed cavity"),
                ((env or {}).get("envelope", 0) > PYR1_ENVELOPE, "F2 envelope <= PYR1"),
                (dev is not None and dev < 1.0, "F3 motif deviation"),
                (app <= MAX_APPENDAGE, "F4 HAB1 appendage"),
                (rg <= PYR1_RG * 1.5, "F6 Rg")]:
            if not ok:
                row.update(verdict="fail", fail_at=why)
                return row
        row.update(verdict="PASS", fail_at="")
        return row
    except Exception as e:                                   # noqa: BLE001
        row.update(verdict="ERROR", fail_at=f"{type(e).__name__}: {e}"[:120])
        return row


def calibrate():
    """wild-type PYR1 through v2. A filter that rejects PYR1 is miscalibrated."""
    at = [(a[0], a[1], a[2], a[3], a[4]) for a in m143.read_pdb(
        os.path.join(ROOT, "data", "stage1", "wt_aba.pdb"), want_atom=False)
        if a[1] in m143.THREE]
    shp = {r["name"]: r for r in json.load(open(os.path.join(
        ROOT, "results", "pocket_shape", "pocket_shape.json")))}
    wall = list(shp["PYR1"]["lining"])
    ca = {a[0]: a[4] for a in at if a[2] == "CA"}
    ps = pieces_of(ca)
    r, _, _ = m192.cavity_and_wall(at)
    env = m192.envelope(at, wall)
    mca = np.array(list(ca.values()))
    rg = float(np.sqrt(((mca - mca.mean(0)) ** 2).sum(1).mean()))
    checks = [("F0 one connected piece", len(ps) == 1, f"{len(ps)} piece(s)"),
              ("F1 enclosed cavity", r and r.get("main", 0) > 0,
               f"{(r or {}).get('main', 0):.0f} A^3"),
              ("F2 envelope >= PYR1", env["envelope"] >= PYR1_ENVELOPE * 0.999,
               f"{env['envelope']:.0f} vs {PYR1_ENVELOPE:.0f}"),
              ("F4 appendage (vacuous)", True, "no designed residues"),
              ("F6 Rg within 1.5x", rg <= PYR1_RG * 1.5, f"Rg {rg:.1f}")]
    ok = True
    print("HARD GATE -- wild-type PYR1 through cascade v2")
    for nm, passed, detail in checks:
        print(f"   {'PASS' if passed else 'FAIL':<5}{nm:<26} {detail}")
        ok &= bool(passed)
    print(f"\n   {'GATE PASSES' if ok else 'GATE FAILS -- RECALIBRATE'}")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="*", default=[])
    ap.add_argument("--procs", type=int, default=32)
    ap.add_argument("--out", default="results/rfd3/cascade2.csv")
    ap.add_argument("--calibrate", action="store_true")
    a = ap.parse_args()
    if a.calibrate:
        return calibrate()

    outp = os.path.join(ROOT, a.out)
    done = set()
    if os.path.exists(outp):                      # resumable
        with open(outp) as fh:
            for row in csv.DictReader(fh):
                done.add((row["arm"], row["design"]))
    jobs = []
    for arm in a.arms:
        d = os.path.join(ROOT, arm)
        cifs = sorted(glob.glob(os.path.join(d, "d*", "*.cif.gz"))) or \
            sorted(glob.glob(os.path.join(d, "*.cif.gz")))
        for c in cifs:
            nm = os.path.basename(os.path.dirname(c))
            if not nm.startswith("d"):
                nm = "m" + os.path.basename(c).split("_model_")[1].split(".")[0]
            if (os.path.basename(d), nm) not in done:
                jobs.append((d, c))
    print(f"{len(jobs)} designs to measure ({len(done)} already done), "
          f"{a.procs} procs", flush=True)
    new = not os.path.exists(outp)
    with open(outp, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if new:
            w.writeheader()
        with mp.Pool(a.procs) as pool:
            for i, row in enumerate(pool.imap_unordered(measure_one, jobs, 1), 1):
                w.writerow(row); fh.flush()
                if row["verdict"] in ("PASS", "ERROR") or i % 25 == 0:
                    print(f"  [{i}/{len(jobs)}] {row['arm']}/{row['design']} "
                          f"{row['verdict']} {row['fail_at']}", flush=True)
    print(f"wrote {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
