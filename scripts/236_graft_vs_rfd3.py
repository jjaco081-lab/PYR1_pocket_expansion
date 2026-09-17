#!/usr/bin/env python
r"""
236_graft_vs_rfd3.py -- ONE figure for ONE decision: graft a pocket from another
helix-grip protein, or keep building with RFd3?

Jannis: "I would like a graph showing why or why not we are moving on with
grafting a pocket from another helix grip motif onto the PYR1 mechanisms and HAB1
binding interface compared to just using RFd3."

BOTH AXES USE THE SAME MEASUREMENT so the comparison is fair: the calibrated
enclosed volume (chamber probe 1.4 A, validated against ABA -- 18 of 19 ABA atoms
inside, PYR1 = 164.4 A^3 apo / 182.3 with HAB1). fpocket is NOT used: its
alpha-sphere volumes include open surface and are 2.18x larger on the same PYR1
pocket, so they cannot be put on the same axis as a design's cavity.

WHAT THE FIGURE HAS TO SHOW, to be decision-relevant:
  * donors: pocket volume against SEQUENCE IDENTITY to PYR1, because identity is
    what decides whether a graft is survivable. A huge pocket at 8 % identity is
    not a donor, it is a different protein.
  * RFd3: the volumes actually achieved by designs that PASS the cascade after
    truncation, drawn as a band across the same y-axis.
  * PYR1 itself as the line to beat.
The question the reader should be able to answer at a glance is: IS THERE A DONOR
THAT IS BOTH BIGGER THAN PYR1 AND CLOSE ENOUGH TO GRAFT -- and if not, does RFd3
reach those volumes anyway?

Figure style per Jannis: no gridlines, capitalised axis titles, reference labels
to the RIGHT, crystal vs predicted by marker SHAPE, minimal small text, no grey
caption block. Numbers are restated in the terminal so we agree on the analysis.
"""
import csv, json, os, sys
import statistics as st
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "figures", "graft_vs_rfd3.png")
PYR1_APO, PYR1_HAB1 = 164.4, 182.3
RED, BLUE, GREY = "#D96C6C", "#3D7EBF", "#BFC4CC"


def main():
    w = {r["target"]: r for r in csv.DictReader(
        open(os.path.join(ROOT, "results", "foldseek", "wall_position_identity.csv")))}
    v5 = [r for r in json.load(open(os.path.join(
        ROOT, "results", "homolog_cavities", "cavity_v5.json")))["structures"]
        if "error" not in r and r["target"] in w]
    don = [(float(w[r["target"]]["fident"]) * 100, r["main_median"],
            r["kind"] == "experimental", r["target"]) for r in v5]
    rf = [float(x["cav_after"] or 0) for x in csv.DictReader(
        open(os.path.join(ROOT, "results", "rfd3", "truncated.csv")))
        if x["verdict_after"] == "PASS" and float(x["cav_after"] or 0) > 0]

    print(f"DONORS  n={len(don)}   RFd3 PASS designs with a cavity  n={len(rf)}")
    print(f"PYR1 reference: {PYR1_APO} apo / {PYR1_HAB1} with HAB1\n")
    rf_med, rf_max = st.median(rf), max(rf)
    rf_p90 = float(np.percentile(rf, 90))
    print(f"RFd3 cavity: median {rf_med:.1f}  p90 {rf_p90:.1f}  max {rf_max:.1f}")
    big = [d for d in don if d[1] > PYR1_APO]
    print(f"donors above PYR1's {PYR1_APO}: {len(big)} of {len(don)}")
    print(f"   their identity: median {st.median([d[0] for d in big]):.1f} %")
    # the decision cells
    for idcut in (20, 25, 30):
        close = [d for d in don if d[0] >= idcut]
        closebig = [d for d in close if d[1] > PYR1_APO]
        beat = [d for d in close if d[1] > rf_max]
        print(f"\n   identity >= {idcut}%:  {len(close):>3} donors, "
              f"{len(closebig):>2} bigger than PYR1, "
              f"{len(beat):>2} bigger than the BEST RFd3 design ({rf_max:.0f})")
        for d in sorted(closebig, key=lambda x: -x[1])[:4]:
            print(f"        {d[3][:38]:<40} {d[1]:7.1f} A^3 at {d[0]:.1f} % id")

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    ax.axhspan(min(rf), rf_max, color=BLUE, alpha=0.13, lw=0, zorder=0)
    ax.axhline(rf_max, color=BLUE, lw=1.4, zorder=1)
    ax.axhline(PYR1_APO, color="black", lw=1.2, ls=(0, (5, 4)), zorder=2)
    for expt, mk, lbl in ((True, "o", "Crystal structure"),
                          (False, "^", "AlphaFold model")):
        s = [d for d in don if d[2] == expt]
        ax.scatter([d[0] for d in s], [d[1] for d in s], marker=mk, s=27,
                   facecolor="none" if not expt else RED,
                   edgecolor=RED, linewidth=1.0, alpha=0.75, zorder=3, label=lbl)
    ax.set_xlabel("Sequence identity to PYR1 (%)")
    ax.set_ylabel("Pocket volume (Å³)")
    ax.set_title("Grafted donors vs RFdiffusion3 designs, measured the same way",
                 fontsize=11, pad=12)
    xr = ax.get_xlim()[1]
    ax.text(xr, rf_max, f"  best RFd3 design  {rf_max:.0f}", va="center",
            ha="left", color=BLUE, fontsize=9)
    ax.text(xr, PYR1_APO, f"  PYR1  {PYR1_APO:.0f}", va="center", ha="left",
            fontsize=9)
    # ⚠ the band label collided with the PYR1 line; put it low in the band
    ax.text(xr, min(rf) + 0.18 * (rf_max - min(rf)),
            "  RFd3 range\n  (655 designs)",
            va="center", ha="left", color=BLUE, fontsize=8.5)
    # the empty quadrant is the whole point -- name it
    ax.annotate("no donor here:\nbig pocket AND graftable",
                xy=(34, 430), ha="center", va="center", fontsize=9.5,
                color="0.35")
    ax.plot([24, 24], [PYR1_APO, 660], color="0.75", lw=1.0, ls=":", zorder=1)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(False)
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    fig.subplots_adjust(right=0.78)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
