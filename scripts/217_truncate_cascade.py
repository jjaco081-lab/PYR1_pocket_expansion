#!/usr/bin/env python
r"""
217_truncate_cascade.py -- delete the HAB1-wrapping TERMINAL appendage from every
design, then re-measure. This is a cascade STEP, not a rescue applied by hand.

JANNIS, twice now: "It has a really long terminal disordered section that wraps
around HAB1. This can probably just be truncated." and "Often the core looks
great and it just has a strange disorder appendage reaching around HAB1."

He is right and the measurement agrees. On arm3b3rasa all 40 F4 failures had the
appendage in the LEAD -- none internal, none trailing -- and deleting it rescued
35/40 with the cavity surviving 40/40 and the envelope essentially unchanged
(median 1253 -> 1238 A^3). It also fixed Rg: 21.8 -> 15.5, 25.9 -> 16.6, so the
"designs are loose, Rg 22-27 vs PYR1's 15" reported since 111 was THE APPENDAGE,
not the fold.

⚠ WHY THIS MATTERS RIGHT NOW. The current arm table is misleading. 6helix has
95 % intact chains and 59 % with a cavity -- the best structural quality of
anything run -- but passes the cascade 0.1 % of the time, while 4rep passes 11 %
with far worse structures. That inversion is an artifact of measuring BEFORE
truncation: long-lead arms wrap HAB1 and fail F4, short-lead arms do not. No arm
comparison is valid until every arm is truncated the same way.

THE RULE, deliberately minimal. Only TERMINAL non-motif residues are removed --
the lead before the first motif residue and the tail after the last. An internal
appendage on an inter-motif linker CANNOT be truncated without breaking the
chain, so those designs are reported as unrescuable rather than quietly trimmed.
From each terminus we delete the SHORTEST piece that removes every HAB1 contact
in that segment, not the whole segment, so nothing is discarded that was not
implicated.

⚠ Truncation is re-scored on the LARGEST CONNECTED PIECE, and the full filter set
is re-applied -- including Rg, which 203 did NOT re-check. That omission made
"35/40 rescued" an understatement rather than an overstatement, but it was still
an untested filter.
"""
import argparse, csv, glob, json, os, sys
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m192 = import_module("192_filter_cascade")
m220 = import_module("220_core_truncate")      # THE truncation rule
m195 = import_module("195_arm_compare")
m198 = import_module("198_cascade2")

CONTACT = 4.5
PYR1_ENVELOPE, PYR1_RG = 1184.0, 15.0
MAX_APPENDAGE = 4
FIELDS = ["arm", "design", "n_main", "n_cut_lead", "n_cut_tail", "cuttable",
          "internal_run",
          "cav_before", "cav_after", "env_before", "env_after",
          "rg_before", "rg_after", "app_before", "app_after",
          "verdict_before", "verdict_after"]


def measure(dat, keep):
    mat = [(a[1], None, a[2], a[3], a[4]) for a in dat if a[1] in keep]
    if len(mat) < 50:
        return None
    r, wall, _ = m192.cavity_and_wall(mat)
    env = m192.envelope(mat, wall) if wall else None
    ca = np.array([a[4] for a in dat if a[1] in keep and a[2] == "CA"])
    rg = float(np.sqrt(((ca - ca.mean(0)) ** 2).sum(1).mean())) if len(ca) else 0.0
    return dict(cav=(r or {}).get("main", 0.0),
                env=(env or {}).get("envelope", 0.0), rg=rg)


def verdict(m, app):
    if m is None:
        return "fail"
    if m["cav"] <= 0:
        return "F1 no cavity"
    if m["env"] <= PYR1_ENVELOPE:
        return "F2 envelope"
    if app > MAX_APPENDAGE:
        return "F4 appendage"
    if m["rg"] > PYR1_RG * 1.5:
        return "F6 Rg"
    return "PASS"


def one(job):
    arm, cif = job
    name = os.path.basename(os.path.dirname(cif))
    row = dict.fromkeys(FIELDS, "")
    row.update(arm=os.path.basename(arm), design=name)
    try:
        js = cif.replace(".cif.gz", ".json")
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
        motif = {int(v[1:]) for v in imap.values() if v[0] == dc}
        main = max(m198.pieces_of(ca), key=len)
        keep0 = set(main)
        if not (motif & keep0):
            row.update(verdict_before="F0 motif not in main piece",
                       verdict_after="unrescuable")
            return row
        lo, hi = min(motif & keep0), max(motif & keep0)
        touch = set()
        for rr in main:
            xyz = np.array(by[rr])
            if len(hab1) and np.linalg.norm(
                    hab1[None, :, :] - xyz[:, None, :], axis=2).min() < CONTACT:
                touch.add(rr)

        def longest_run(sel):
            best = cur = 0
            for k in sorted(sel):
                if k in touch and k not in motif:
                    cur += 1; best = max(best, cur)
                else:
                    cur = 0
            return best

        app0 = longest_run(keep0)
        m0 = measure(dat, keep0)
        # CORE-GROWTH truncation (220): keep what packs against the motif core.
        # Replaces the HAB1-contact rule, which Jannis showed fails when an
        # appendage does not reach HAB1, and the local-packing rule, which fails
        # when the appendage packs against ITSELF.
        keep1, n_lead, n_tail, run = m220.truncate(by, keep0, motif)
        cut = keep0 - keep1
        internal = run >= m220.MIN_INTERNAL_RUN
        m1 = measure(dat, keep1)
        app1 = longest_run(keep1)
        row.update(n_main=len(main), n_cut_lead=n_lead, n_cut_tail=n_tail,
                   cuttable=int(not internal), internal_run=run,
                   cav_before=round(m0["cav"], 1) if m0 else "",
                   cav_after=round(m1["cav"], 1) if m1 else "",
                   env_before=round(m0["env"], 1) if m0 else "",
                   env_after=round(m1["env"], 1) if m1 else "",
                   rg_before=round(m0["rg"], 1) if m0 else "",
                   rg_after=round(m1["rg"], 1) if m1 else "",
                   app_before=app0, app_after=app1,
                   verdict_before=verdict(m0, app0),
                   verdict_after=verdict(m1, app1))
    except Exception as e:                                       # noqa: BLE001
        row.update(verdict_before=f"ERROR:{type(e).__name__}")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=64)
    ap.add_argument("--out", default="results/rfd3/truncated.csv")
    a = ap.parse_args()
    arms = sorted(glob.glob(os.path.join(ROOT, "results", "rfd3", "arm*")))
    jobs = []
    for d in arms:
        for c in sorted(glob.glob(os.path.join(d, "d*", "*.cif.gz"))):
            jobs.append((d, c))
    print(f"{len(jobs)} designs across {len(arms)} arms, {a.procs} procs", flush=True)
    outp = os.path.join(ROOT, a.out)
    with open(outp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        with mp.Pool(a.procs) as pool:
            for i, row in enumerate(pool.imap_unordered(one, jobs, 4), 1):
                w.writerow(row); fh.flush()
                if i % 200 == 0:
                    print(f"  [{i}/{len(jobs)}]", flush=True)
    print(f"wrote {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
