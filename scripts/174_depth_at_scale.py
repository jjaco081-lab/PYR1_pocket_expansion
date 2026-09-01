#!/usr/bin/env python
r"""
174_depth_at_scale.py -- the §98c depth/gate-contact lead, at n=362 instead of 6.

§98c found that the two compounds with ZERO responders make far MORE gate+latch
contacts (median 20) than the four with responders (median 6), and sit closer to
the mouth (4.8-6.6 A vs 7.5-8.2). The suggested mechanism is obstruction: a
ligand in the gate's path prevents the closure the assay reads.

At n=6, with both dead compounds carrying predicted poses AND clashing, that is
a confounded lead. §84's tractability set is the same question at scale: 181
ligands PYR1 yielded a sensor for, against 181 property-matched ligands it did
not, all co-folded with the SAME wild-type receptor and all passing the fold
assertions of §93.

⚠ §93a ALREADY tested a depth measure here and found nothing: ligand centroid to
the gate/latch mouth, 8.07 A for hits vs 8.00 for non-hits, AUC 0.508. This runs
the CONTACT-COUNT form instead, which is what §98c actually separated on, and
which is not the same quantity -- a large ligand can have a distant centroid and
still reach into the gate.

⚠ Contact count scales with ligand size, and the two dead compounds in §98c were
the largest tested. Heavy-atom count is therefore reported alongside, and the
contact count is also normalised per heavy atom, so "big ligands touch more
things" cannot masquerade as a mechanism.

⚠ Labels here are LIGAND-level (some variant bound it) against §98c's
variant-level labels, and the receptor is wild type rather than K59R. This tests
the same geometric idea, not the identical quantity.
"""
import glob, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from importlib import import_module
G = import_module("164_pose_geometry")

T = os.path.join(ROOT, "data", "tractability")
GATE = [85, 86, 87, 88, 89]
LATCH = [115, 116, 117]


def auc(v, lab):
    """Rank-sum form: O(n log n), not the O(n^2) pairwise loop. Identical value,
    verified against the pairwise count in 95's three-way check."""
    v = np.asarray(v, float); lab = np.asarray(lab, bool)
    n1, n2 = int(lab.sum()), int((~lab).sum())
    if not n1 or not n2:
        return float("nan")
    order = np.argsort(v)
    ranks = np.empty(len(v)); ranks[order] = np.arange(1, len(v) + 1)
    # average ranks over ties so the value matches the 0.5-for-tie convention
    for u in np.unique(v):
        m = v == u
        if m.sum() > 1:
            ranks[m] = ranks[m].mean()
    return float(1 - (ranks[lab].sum() - n1 * (n1 + 1) / 2) / (n1 * n2))


def perm_p(v, lab, B=2000):
    a = auc(v, lab); rng = np.random.default_rng(0); lab = np.asarray(lab)
    null = np.array([auc(v, rng.permutation(lab)) for _ in range(B)])
    return a, float(((null - .5) >= abs(a - .5)).mean() * 2)


def main():
    fbb, aba, _ = G.frame()
    core = [r for r in fbb if 6 <= r <= 180 and r not in GATE + LATCH]
    prot = {}
    for l in open(os.path.join(ROOT, "data", "stage1", "wt_aba.pdb")):
        if l.startswith("ATOM") and l[76:78].strip() != "H":
            prot.setdefault(int(l[22:26]), []).append(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    prot = {k: np.array(v) for k, v in prot.items()}
    gate = np.vstack([prot[r] for r in GATE if r in prot])
    latch = np.vstack([prot[r] for r in LATCH if r in prot])
    mouth = np.vstack([prot[r] for r in GATE + LATCH if r in prot]).mean(0)
    key = json.load(open(os.path.join(T, "key_matched.json")))

    rows = []
    for lig, meta in key.items():
        g = sorted(glob.glob(os.path.join(T, "out", f"boltz_results_{lig}",
                                          "predictions", "*", "*.cif")))
        if not g:
            continue
        bb, L, _, pl = G.read(g[0])
        if len(L) == 0 or pl < 80:
            continue
        P, Q = G.paired(fbb, bb, core)
        if len(P) < 400:
            continue
        R, mp, mq = G.kabsch(P, Q)
        Lf = (L - mq) @ R + mp
        ng = int((np.linalg.norm(Lf[:, None] - gate[None], axis=-1).min(1) < 4.5).sum())
        nl = int((np.linalg.norm(Lf[:, None] - latch[None], axis=-1).min(1) < 4.5).sum())
        rows.append(dict(lig=lig, hit=bool(meta["label"]), heavy=meta["heavy"],
                         gate=ng, latch=nl, gl=ng + nl,
                         gl_per_atom=(ng + nl) / max(len(Lf), 1),
                         depth=float(np.linalg.norm(Lf.mean(0) - mouth))))
    lab = [r["hit"] for r in rows]
    print(f"{len(rows)} ligands ({sum(lab)} hits, {len(lab)-sum(lab)} matched non-hits)")
    print(f"\n{'quantity':<24}{'hits':>9}{'non-hits':>11}{'AUC':>8}{'perm p':>9}"
          f"   direction tested")
    tests = [("gate+latch contacts", "gl", False, "FEWER contacts = hit (98c)"),
             ("contacts per atom", "gl_per_atom", False, "same, size-normalised"),
             ("depth to mouth", "depth", False, "DEEPER = hit (98c)"),
             ("gate contacts only", "gate", False, "fewer = hit"),
             ("latch contacts only", "latch", False, "fewer = hit"),
             ("heavy atoms", "heavy", True, "smaller = hit (control)")]
    out = {}
    for name, k, lo_is_hit, desc in tests:
        v = np.array([r[k] for r in rows], float)
        a, p = perm_p(v if lo_is_hit else -v, lab)
        h = np.median(v[np.array(lab)]); n = np.median(v[~np.array(lab)])
        print(f"{name:<24}{h:>9.2f}{n:>11.2f}{a:>8.3f}{p:>9.4f}   {desc}")
        out[name] = dict(auc=a, p=p, hits=float(h), nonhits=float(n))
    print("\n  AUC > 0.5 supports the stated direction. §98c predicts >0.5 for the")
    print("  first three rows: hits should make FEWER gate/latch contacts and sit")
    print("  DEEPER than non-hits.")
    json.dump({"rows": rows, "tests": out},
              open(os.path.join(ROOT, "results", "tractability", "depth_at_scale.json"), "w"),
              indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
