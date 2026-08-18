#!/usr/bin/env python
"""
60_loop_dynamics_figure.py -- figures for the gate/latch loop-dynamics result,
over all four WT systems (12 trajectories, 15 protomer-level analysis units).

Six panels:
  A  gate state coordinate S = d(open) - d(closed) vs time, every replicate.
     The zero line is the open/closed watershed. This is the "does it hold its
     state" panel.
  B  two-reference scatter: d(closed) against d(open), post-equilibration frames.
     Open frames land bottom-right, closed frames top-left. The crystal references
     themselves are marked, so the axes have a physical scale rather than an
     arbitrary one.
  C  gate-latch minimum heavy-atom distance, distribution per replicate.
  D  per-residue backbone RMSF with the gate, latch and Lb7a5 loops shaded.
  E  latch state coordinate vs time -- the second pre-registered loop, shown in
     its own panel now that there is an apo-closed protomer to test it on.
  F  gate drift: mean S over the first 50 ns after the discard against the last
     50 ns. Points below the diagonal moved toward OPEN. This is the panel that
     answers "does a closed protomer open once the ligand is gone".

Replicates are drawn individually, never averaged into a single line, because the
inferential unit here is the replicate and hiding that inside a mean would
misrepresent how much evidence there actually is.

COLOUR SCHEME. Blues are open-like states, oranges closed-like, purple the
ternary; the DARK shade of each pair is the monomer system and the LIGHT shade is
the corresponding protomer of the S3 dimer. Blue/orange is the safest pair for
red-green colour blindness, and the light/dark split carries the monomer/dimer
distinction without needing a third hue.

Run with the pyr1_docking env python, after 59.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LD = os.path.join(ROOT, "data", "loop_dynamics")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

MAP = json.load(open(os.path.join(LD, "residue_map.json")))
SUM = json.load(open(os.path.join(LD, "summary.json")))
# per-unit, because S4's numbering is shifted by one relative to S1/S2 and S3's
# second protomer does not start at prmtop residue 1 at all
SEQ2NAT = {s: {int(k): v for k, v in m["seq_to_native"].items()}
           for s, m in MAP["systems"].items()}
DISCARD = SUM["discard_ns"]
PSF = SUM["ps_per_frame"]
DF = int(DISCARD * 1000 / PSF)

# unit -> (legend label, colour, short tick label)
SYSTEMS = {
    "S1_apo_open":      ("open, apo (S1, 3K3K A monomer)",      "#1F4E79", "open"),
    "S3_dimer_openA":   ("open, apo (S3, dimer protomer A)",    "#6BAED6", "dimerA"),
    "S2_holo_closed":   ("closed +ABA (S2, 3QN1 A monomer)",    "#A63603", "closed"),
    "S3_dimer_closedB": ("closed, APO (S3, dimer protomer B)",  "#FD8D3C", "apo-closed"),
    "S4_ternary":       ("ternary +ABA +HAB1 (S4)",             "#6A51A3", "ternary"),
}
# drawing order puts the apo-closed protomer last so it is never buried
ORDER = ["S1_apo_open", "S3_dimer_openA", "S2_holo_closed", "S4_ternary",
         "S3_dimer_closedB"]
HERO = "S3_dimer_closedB"
LOOPS = {"gate": (MAP["gate_native"], "#8c564b"),
         "latch": (MAP["latch_native"], "#9467bd"),
         "Lb7a5": (MAP["lb7a5_native"], "#2ca02c")}
REPS = (0, 1, 2)


def dat(s, r, name, col=1):
    p = os.path.join(LD, s, f"rep{r}", f"{name}.dat")
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return None
    return np.loadtxt(p, comments="#", usecols=(col,), ndmin=1).astype(float)


fig, axes = plt.subplots(3, 2, figsize=(16, 15))
axA, axB, axC, axD, axE, axF = axes.ravel()


def state_coord(s, r, loop):
    o, c = dat(s, r, f"{loop}_to_open"), dat(s, r, f"{loop}_to_closed")
    return None if o is None or c is None else o - c


# ---------------- A / E: state coordinate vs time ----------------
for ax, loop in ((axA, "gate"), (axE, "latch")):
    for s in ORDER:
        lab, col, _ = SYSTEMS[s]
        first = True
        for r in REPS:
            S = state_coord(s, r, loop)
            if S is None:
                continue
            t = np.arange(len(S)) * PSF / 1000.0
            ax.plot(t, S, color=col, lw=0.9 if s == HERO else 0.5,
                    alpha=0.9 if s == HERO else 0.6,
                    label=lab if first else None, zorder=3 if s == HERO else 2)
            first = False
    ax.axhline(0, color="k", lw=1.2, ls="--")
    ax.axvline(DISCARD, color="grey", lw=1, ls=":")
    ax.set_xlabel("time (ns)")
    ax.set_ylabel(f"S = {loop} RMSD to open  -  to closed  ($\\AA$)")

axA.text(0.99, 0.51, "open / closed watershed", transform=axA.transAxes,
         fontsize=7.5, color="k", ha="right", va="bottom")
axA.set_title("A. Gate: does each state hold?  S<0 open-like, S>0 closed-like")
axA.legend(fontsize=8, loc="center left", framealpha=0.95)
axE.set_title("E. Latch: the second pre-registered loop")

# ---------------- B: two-reference scatter ----------------
for s in ORDER:
    lab, col, _ = SYSTEMS[s]
    for r in REPS:
        go, gc = dat(s, r, "gate_to_open"), dat(s, r, "gate_to_closed")
        if go is None:
            continue
        axB.scatter(go[DF::20], gc[DF::20], s=1.2, alpha=0.10, color=col,
                    rasterized=True, zorder=3 if s == HERO else 2)
ref = MAP["ref_open_vs_closed"]["gate_bb_rmsd"]
axB.scatter([0], [ref], marker="*", s=260, color="#1F4E79", edgecolor="k", zorder=5)
axB.scatter([ref], [0], marker="*", s=260, color="#A63603", edgecolor="k", zorder=5)
lim = max(axB.get_xlim()[1], axB.get_ylim()[1])
axB.plot([0, lim], [0, lim], color="k", lw=1, ls="--")
# The diagonal is the equidistance line. Saying so on the figure matters: a single
# RMSD cannot distinguish "converted to the other state" from "fell apart in some
# third direction", and the whole point of two references is that it can.
axB.text(0.30 * lim, 0.86 * lim, "closer to OPEN", fontsize=9, style="italic",
         color="#1F4E79", ha="center")
axB.text(0.78 * lim, 0.40 * lim, "closer to CLOSED", fontsize=9, style="italic",
         color="#A63603", ha="center")
axB.set_xlabel("gate backbone RMSD to OPEN reference ($\\AA$)")
axB.set_ylabel("gate backbone RMSD to CLOSED reference ($\\AA$)")
axB.set_title("B. Two-reference projection (post-equilibration)")
# Proxy handles: the real ones need markerscale to make a 1-pt scatter visible,
# which also inflates the reference stars into blobs that cover the data.
axB.legend(handles=[Line2D([], [], marker="o", ls="", ms=6,
                           color=SYSTEMS[s][1], label=SYSTEMS[s][0]) for s in ORDER]
                   + [Line2D([], [], marker="*", ls="", ms=13, color="#1F4E79",
                             mec="k", label="open crystal (3K3K A)"),
                      Line2D([], [], marker="*", ls="", ms=13, color="#A63603",
                             mec="k", label="closed crystal (3QN1 A)")],
           fontsize=7.5, loc="lower right", framealpha=0.9)

# ---------------- C: gate-latch min distance ----------------
pos, labels, colors = [], [], []
for s in ORDER:
    _, col, short = SYSTEMS[s]
    for r in REPS:
        d = dat(s, r, "gate_latch_contacts", 3)
        if d is None:
            continue
        pos.append(d[DF:])
        labels.append(f"{short}\nrep{r}")
        colors.append(col)
if pos:
    bp = axC.violinplot(pos, showmeans=True, showextrema=False)
    for body, c in zip(bp["bodies"], colors):
        body.set_facecolor(c)
        body.set_alpha(0.65)
    axC.set_xticks(range(1, len(labels) + 1))
    axC.set_xticklabels(labels, fontsize=6.5, rotation=35, ha="right")
axC.axhline(4.5, color="grey", ls=":", lw=1)
axC.set_ylabel("gate-latch minimum heavy-atom distance ($\\AA$)")
axC.set_title("C. Gate-latch closure, per replicate")

# ---------------- D: per-residue RMSF ----------------
for s in ORDER:
    lab, col, _ = SYSTEMS[s]
    first = True
    for r in REPS:
        p = os.path.join(LD, s, f"rep{r}", "rmsf_bb_byres.dat")
        if not os.path.exists(p):
            continue
        a = np.loadtxt(p, comments="#", ndmin=2)
        nat = [SEQ2NAT[s].get(int(x), np.nan) for x in a[:, 0]]
        axD.plot(nat, a[:, 1], color=col, lw=1.1 if s == HERO else 0.8,
                 alpha=0.9 if s == HERO else 0.65, label=lab if first else None,
                 zorder=3 if s == HERO else 2)
        first = False
# Both termini fluctuate 6-14 A -- they are unstructured tails, not a result --
# and letting them set the y-axis squashes the loops this panel exists to compare.
# Clipped rather than trimmed, and said so, so nothing is silently hidden.
axD.set_ylim(0, 6)
for lname, (nats, lcol) in LOOPS.items():
    axD.axvspan(min(nats) - 0.5, max(nats) + 0.5, color=lcol, alpha=0.16)
    axD.text(np.mean(nats), 5.6, lname, ha="center",
             fontsize=9, color=lcol, fontweight="bold")
axD.text(0.99, 0.02, "y clipped at 6 $\\AA$; both termini exceed it",
         transform=axD.transAxes, ha="right", fontsize=7.5, style="italic",
         color="grey")
axD.set_xlabel("PYR1 residue (native numbering)")
axD.set_ylabel("backbone RMSF ($\\AA$)")
axD.set_title("D. Per-residue flexibility")
axD.legend(fontsize=7.5, loc="upper left")

# ---------------- F: gate drift, first 50 ns vs last 50 ns ----------------
# The apo-closed question in one panel. A point ON the diagonal did not move; a
# point BELOW it drifted toward the open basin over 300 ns. Drawn from the same
# pre-registered state coordinate as panel A -- no new observable.
WIN = int(50 * 1000 / PSF)
for s in ORDER:
    lab, col, _ = SYSTEMS[s]
    first = True
    for r in REPS:
        S = state_coord(s, r, "gate")
        if S is None:
            continue
        seg = S[DF:]
        if len(seg) < 2 * WIN:
            continue
        axF.scatter(seg[:WIN].mean(), seg[-WIN:].mean(), s=110, color=col,
                    edgecolor="k", linewidth=0.6, zorder=4 if s == HERO else 3,
                    label=lab if first else None)
        first = False
lo = min(axF.get_xlim()[0], axF.get_ylim()[0])
hi = max(axF.get_xlim()[1], axF.get_ylim()[1])
axF.plot([lo, hi], [lo, hi], color="k", lw=1, ls="--")
axF.axhline(0, color="grey", lw=0.8, ls=":")
axF.axvline(0, color="grey", lw=0.8, ls=":")
axF.fill_between([lo, hi], [lo, lo], [lo, hi], color="#1F4E79", alpha=0.07)
axF.text(0.97, 0.06, "below the diagonal = drifted toward OPEN",
         transform=axF.transAxes, ha="right", fontsize=8, style="italic",
         color="#1F4E79")
axF.set_xlabel("mean gate S, first 50 ns after discard ($\\AA$)")
axF.set_ylabel("mean gate S, last 50 ns ($\\AA$)")
axF.set_title("F. Did the gate move over 300 ns?")
axF.legend(fontsize=7.5, loc="upper left", framealpha=0.9)

n = SUM.get("n_replicates", {})
counts = ", ".join(f"{SYSTEMS[u][2]} x{n[st]}"
                   for u, st in (("S1_apo_open", "open"),
                                 ("S3_dimer_openA", "dimer_openA"),
                                 ("S2_holo_closed", "closed"),
                                 ("S3_dimer_closedB", "dimer_closedB"),
                                 ("S4_ternary", "ternary")) if n.get(st))
fig.suptitle(
    f"PYR1 gate/latch loop dynamics -- WT MD baseline, 300 ns x  {counts}   "
    f"(first {DISCARD} ns discarded; {SUM.get('core_residues', '?')}-residue core fit)",
    fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = os.path.join(FIG, "loop_dynamics.png")
fig.savefig(out, dpi=160)
print(f"wrote {out}")
