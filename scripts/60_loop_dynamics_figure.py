#!/usr/bin/env python
"""
60_loop_dynamics_figure.py -- figures for the gate/latch loop-dynamics result.

Four panels:
  A  state coordinate S = d(open) - d(closed) vs time, every replicate.
     The zero line is the open/closed watershed. This is the "does it hold its
     state" panel.
  B  two-reference scatter: d(closed) against d(open), post-equilibration frames.
     Open frames land bottom-right, closed frames top-left. The crystal references
     themselves are marked, so the axes have a physical scale rather than an
     arbitrary one.
  C  gate-latch minimum heavy-atom distance, distribution per replicate.
  D  per-residue backbone RMSF with the gate, latch and Lb7a5 loops shaded.

Replicates are drawn individually, never averaged into a single line, because the
inferential unit here is the replicate (n=3 vs n=3) and hiding that inside a mean
would misrepresent how much evidence there actually is.

Run with the pyr1_docking env python, after 59.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LD = os.path.join(ROOT, "data", "loop_dynamics")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

MAP = json.load(open(os.path.join(LD, "residue_map.json")))
SUM = json.load(open(os.path.join(LD, "summary.json")))
# per-system, because S4's numbering is shifted by one relative to S1/S2
SEQ2NAT = {s: {int(k): v for k, v in m["seq_to_native"].items()}
           for s, m in MAP["systems"].items()}
DISCARD = SUM["discard_ns"]
PSF = SUM["ps_per_frame"]
DF = int(DISCARD * 1000 / PSF)

SYSTEMS = {"S1_apo_open": ("open (apo, 3K3K A)", "#1f77b4"),
           "S2_holo_closed": ("closed (+ABA, 3QN1 A)", "#d62728"),
           "S4_ternary": ("ternary (+ABA +HAB1)", "#7f4fa8")}
LOOPS = {"gate": (MAP["gate_native"], "#8c564b"),
         "latch": (MAP["latch_native"], "#9467bd"),
         "Lb7a5": (MAP["lb7a5_native"], "#2ca02c")}


def dat(s, r, name, col=1):
    p = os.path.join(LD, s, f"rep{r}", f"{name}.dat")
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return None
    return np.loadtxt(p, comments="#", usecols=(col,), ndmin=1).astype(float)


fig, axes = plt.subplots(2, 2, figsize=(15, 10))
axA, axB, axC, axD = axes.ravel()

# ---------------- A: state coordinate vs time ----------------
for s, (lab, col) in SYSTEMS.items():
    first = True
    for r in (0, 1, 2):
        go, gc = dat(s, r, "gate_to_open"), dat(s, r, "gate_to_closed")
        if go is None:
            continue
        S = go - gc
        t = np.arange(len(S)) * PSF / 1000.0
        axA.plot(t, S, color=col, lw=0.5, alpha=0.65,
                 label=lab if first else None)
        first = False
axA.axhline(0, color="k", lw=1.2, ls="--")
axA.text(0.99, 0.5, "open / closed watershed", transform=axA.transAxes,
         fontsize=7.5, color="k", ha="right", va="bottom")
axA.axvline(DISCARD, color="grey", lw=1, ls=":")
axA.set_xlabel("time (ns)")
axA.set_ylabel("S = gate RMSD to open  -  to closed  ($\\AA$)")
axA.set_title("A. Does each state hold?  S<0 open-like, S>0 closed-like")
axA.legend(fontsize=9, loc="best")

# ---------------- B: two-reference scatter ----------------
for s, (lab, col) in SYSTEMS.items():
    first = True
    for r in (0, 1, 2):
        go, gc = dat(s, r, "gate_to_open"), dat(s, r, "gate_to_closed")
        if go is None:
            continue
        axB.scatter(go[DF::10], gc[DF::10], s=1, alpha=0.10, color=col, rasterized=True)
        first = False
ref = MAP["ref_open_vs_closed"]["gate_bb_rmsd"]
axB.scatter([0], [ref], marker="*", s=260, color="#1f77b4", edgecolor="k", zorder=5)
axB.scatter([ref], [0], marker="*", s=260, color="#d62728", edgecolor="k", zorder=5)
lim = max(axB.get_xlim()[1], axB.get_ylim()[1])
axB.plot([0, lim], [0, lim], color="k", lw=1, ls="--")
# The diagonal is the equidistance line. Saying so on the figure matters: a single
# RMSD cannot distinguish "converted to the other state" from "fell apart in some
# third direction", and the whole point of two references is that it can.
axB.text(0.30 * lim, 0.86 * lim, "closer to OPEN", fontsize=9, style="italic",
         color="#1f77b4", ha="center", rotation=0)
# placed below the diagonal on the right, where nothing is plotted and the
# lower-right legend cannot clip it
axB.text(0.75 * lim, 0.42 * lim, "closer to CLOSED", fontsize=9, style="italic",
         color="#d62728", ha="center")
axB.set_xlabel("gate backbone RMSD to OPEN reference ($\\AA$)")
axB.set_ylabel("gate backbone RMSD to CLOSED reference ($\\AA$)")
axB.set_title("B. Two-reference projection (post-equilibration)")
# Build the legend from proxy handles. Using the real ones needs markerscale to
# make the 1-pt scatter visible, which also inflates the reference stars into
# giant blobs that cover the data.
from matplotlib.lines import Line2D  # noqa: E402

axB.legend(handles=[
    Line2D([], [], marker="o", ls="", ms=6, color=SYSTEMS["S1_apo_open"][1],
           label=SYSTEMS["S1_apo_open"][0]),
    Line2D([], [], marker="o", ls="", ms=6, color=SYSTEMS["S2_holo_closed"][1],
           label=SYSTEMS["S2_holo_closed"][0]),
    Line2D([], [], marker="o", ls="", ms=6, color=SYSTEMS["S4_ternary"][1],
           label=SYSTEMS["S4_ternary"][0]),
    Line2D([], [], marker="*", ls="", ms=13, color="#1f77b4", mec="k",
           label="open crystal (3K3K A)"),
    Line2D([], [], marker="*", ls="", ms=13, color="#d62728", mec="k",
           label="closed crystal (3QN1 A)"),
], fontsize=8, loc="lower right", framealpha=0.9)

# ---------------- C: gate-latch min distance ----------------
pos, labels, colors = [], [], []
for s, (lab, col) in SYSTEMS.items():
    for r in (0, 1, 2):
        d = dat(s, r, "gate_latch_contacts", 3)
        if d is None:
            continue
        pos.append(d[DF:])
        labels.append({'S1_apo_open': 'open', 'S2_holo_closed': 'closed',
                       'S4_ternary': 'ternary'}[s] + f'\nrep{r}')
        colors.append(col)
if pos:
    bp = axC.violinplot(pos, showmeans=True, showextrema=False)
    for body, c in zip(bp["bodies"], colors):
        body.set_facecolor(c)
        body.set_alpha(0.6)
    axC.set_xticks(range(1, len(labels) + 1))
    axC.set_xticklabels(labels, fontsize=8)
axC.axhline(4.5, color="grey", ls=":", lw=1)
axC.set_ylabel("gate-latch minimum heavy-atom distance ($\\AA$)")
axC.set_title("C. Gate-latch closure, per replicate")

# ---------------- D: per-residue RMSF ----------------
for s, (lab, col) in SYSTEMS.items():
    first = True
    for r in (0, 1, 2):
        p = os.path.join(LD, s, f"rep{r}", "rmsf_bb_byres.dat")
        if not os.path.exists(p):
            continue
        a = np.loadtxt(p, comments="#", ndmin=2)
        nat = [SEQ2NAT[s].get(int(x), np.nan) for x in a[:, 0]]
        axD.plot(nat, a[:, 1], color=col, lw=0.9, alpha=0.7,
                 label=lab if first else None)
        first = False
for lname, (nats, lcol) in LOOPS.items():
    axD.axvspan(min(nats) - 0.5, max(nats) + 0.5, color=lcol, alpha=0.16)
    axD.text(np.mean(nats), axD.get_ylim()[1] * 0.94, lname, ha="center",
             fontsize=9, color=lcol, fontweight="bold")
axD.set_xlabel("PYR1 residue (native numbering)")
axD.set_ylabel("backbone RMSF ($\\AA$)")
axD.set_title("D. Per-residue flexibility")
axD.legend(fontsize=9, loc="upper left")

fig.suptitle(
    f"PYR1 gate/latch loop dynamics -- WT MD baseline, 300 ns x "
    f"{SUM['n_open']} open + {SUM['n_closed']} closed (+ ternary) replicates "
    f"(first {DISCARD} ns discarded)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(FIG, "loop_dynamics.png")
fig.savefig(out, dpi=160)
print(f"wrote {out}")
