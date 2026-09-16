#!/usr/bin/env python
r"""215_junction_figure.py -- the chimera-feasibility figure.

Replaces the cavity-vs-identity scatter, which Jannis correctly called out as
stating the obvious (close relatives are similar and have similar pockets).
The specific, non-obvious claim is about JUNCTIONS: where can the two chains be
spliced, and how much backbone has to be rebuilt.

Junction cost = total CA deviation, after superposition on the REAL foldseek
alignment, at the optimal splice point for each of the 8 crossovers required to
keep PYR1's HAB1 interface, gate, latch and tunnel residues (flank +/-2).
Every position in every gap is scanned, so the cost is the minimum achievable.

THE RESULT IS A NULL ON THE OBVIOUS AXIS: Spearman(cavity, junction cost) =
-0.131, p = 0.42. Large pockets are NOT harder to splice. Api g 1 gives 2.1x
PYR1's cavity with 7 of 8 clean junctions.

House style: no caption block, no gridlines, Capitalised axes, reference label
right of the line, marker shape for provenance, 20 % smaller canvas.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RED, BLUE, GREY, INK, DGREY = "#D96C6C", "#3D7EBF", "#BFC4CC", "#22262B", "#8A8F98"
d = json.load(open(f"{ROOT}/results/schema_chimera/foldseek_flank2.json"))
NAME = {"2PCS": "CoxG", "2BK0": "Api g 1", "6AWV": "Ara h 8",
        "3OQU": "PYL9", "4DSB": "PYL3", "3QRZ": "PYL5"}
BIG = {"2BK0", "6AWV", "2PCS"}

cav = np.array([r["cavity"] for r in d])
jc = np.array([r["junction_cost"] for r in d])

# 20 % smaller than the 9.5 x 6.0 canvas
fig, ax = plt.subplots(figsize=(7.6, 4.8), dpi=200)
ax.scatter(cav, jc, s=34, c=GREY, edgecolors="white", linewidths=0.5, zorder=2)
ax.axvline(164, color=BLUE, lw=2.0, ls=(0, (5, 3)), zorder=1)
ax.text(172, ax.get_ylim()[1] * 0.96, "PYR1", fontsize=12, color=BLUE,
        ha="left", va="top", fontweight="bold")

OFF = {"2BK0": (14, 10), "6AWV": (14, -6), "2PCS": (-14, 8),
       "3OQU": (16, 16), "4DSB": (-14, -18), "3QRZ": (14, -12)}
# ⚠ 3QRZ appears TWICE -- chains A and B of the same crystal, measured at 104
# and 210 A^3. Same protein, and the 2x spread between its own chains is a real
# caveat on every single-chain cavity number here. Label it once, at the larger.
seen = set()
for r in sorted(d, key=lambda x: -x["cavity"]):
    p = r["pdb"]
    if p not in NAME or p in seen:
        continue
    seen.add(p)
    col = RED if p in BIG else DGREY
    ax.scatter([r["cavity"]], [r["junction_cost"]], s=110, c=col,
               edgecolors="white", linewidths=1.5, zorder=5)
    ax.annotate(f"{NAME[p]}  ({p})", (r["cavity"], r["junction_cost"]),
                textcoords="offset points", xytext=OFF[p], fontsize=10,
                color=INK, va="center",
                ha="right" if OFF[p][0] < 0 else "left",
                arrowprops=dict(arrowstyle="-", color="#CBD0D6", lw=0.8,
                                shrinkA=0, shrinkB=6))

ax.set_xlabel("Cavity Volume  (Å$^3$)", fontsize=12, color=INK)
ax.set_ylabel("Junction Cost  (Å)", fontsize=12, color=INK)
ax.set_title("Bigger pockets are not harder to splice", fontsize=14,
             color=INK, pad=12, loc="left")
ax.grid(False)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_color("#C8CDD4")
ax.tick_params(colors=INK, labelsize=11)
fig.tight_layout()
out = f"{ROOT}/results/figures/chimera_junctions.png"
fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
print(f"wrote {out}")
from scipy.stats import spearmanr
print(f"  Spearman(cavity, junction cost) = {spearmanr(cav, jc)[0]:+.3f}, "
      f"p = {spearmanr(cav, jc)[1]:.2f}  (n = {len(d)})")
