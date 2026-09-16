#!/usr/bin/env python
r"""212_graft_figure.py -- the categorical-tension figure (Jannis house style).

STYLE RULES (Jannis, 2026-09-16), applied here and to be reused:
  * no grey subtitle / caption blocks -- the key numbers go in the response text
    and on the slide, not in small type inside the figure
  * no gridlines when individual points are labelled
  * axis titles Capitalised
  * reference lines labelled to the RIGHT of the line
  * minimal text overall
  * crystal vs predicted distinguished by MARKER SHAPE (circle vs triangle)

⚠ THE SET IS NOT ALL CRYSTAL STRUCTURES. 102 of the 147 are AlphaFold models and
only 45 are experimental PDB entries, so provenance is encoded by MARKER FILL
(filled = experimental, open = predicted) rather than claimed in the title.
Cavity volumes on predicted structures depend on modelled side chains, so the
experimental points carry the argument.
"""
import csv, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RED, BLUE, GREY, INK = "#D96C6C", "#3D7EBF", "#BFC4CC", "#22262B"
DGREY = "#8A8F98"

rows = list(csv.DictReader(open(f"{ROOT}/results/homolog_cavities/homolog_cavities.csv")))
cav = np.array([float(r["cavity_A3"]) for r in rows])
fid = np.array([100 * float(r["fident"]) for r in rows])
isaf = np.array([r["target"].startswith("af_") for r in rows])
geo = {g["pdb"]: g for g in json.load(open(f"{ROOT}/results/graft_candidates/graft_geometry.json"))}
PYR1_CAV = 164.0

#: labelled points. large donors in red; PYR/PYL relatives in ink.
#: Protein NAMES, read from each structure's own _entity.pdbx_description --
#: far more useful to a protein audience than a PDB code. Note what the large
#: donors actually are: two are FOOD ALLERGENS of the Bet v 1 family and one is
#: a hypothetical protein. None is a receptor.
NAME = {"2PCS": "CoxG", "2NS9": "APE2225", "2BK0": "Api g 1", "6AWV": "Ara h 8",
        "3OQU": "PYL9", "3QRZ": "PYL5", "4DSB": "PYL3"}
#: ⚠ 2PCS's PDB entity description is only 'Conserved protein'; CoxG is the
#: identification from the Mo-CODH literature (README 1440-1488), not the file.
#: Volume and RMSD are NOT printed on the points: the y-axis already encodes
#: volume, so a numeric label would be a duplicate encoding (Tufte eraser test).
DONOR = {"2PCS": (16, 6), "2NS9": (16, 6), "2BK0": (16, 10), "6AWV": (18, -22)}
RELATIVE = {"3OQU": (12, 14), "3QRZ": (12, -16), "4DSB": (-16, -24)}

fig, ax = plt.subplots(figsize=(9.5, 6.0), dpi=200)
# shape, not just fill: survives greyscale printing and colour-blind viewing
ax.scatter(fid[isaf], cav[isaf], s=34, marker="^", facecolors="none",
           edgecolors=GREY, linewidths=1.1, zorder=2)
ax.scatter(fid[~isaf], cav[~isaf], s=34, marker="o", c=GREY,
           edgecolors="white", linewidths=0.5, zorder=2)

ax.axhline(PYR1_CAV, color=BLUE, lw=2.0, ls=(0, (5, 3)), zorder=3)
ax.text(60.8, PYR1_CAV, "PYR1", fontsize=12, color=BLUE,
        ha="left", va="center", fontweight="bold")

for pdb, off in DONOR.items():
    g = geo.get(pdb)
    if not g:
        continue
    x, y = 100 * g["fident"], g["cavity"]
    ax.scatter([x], [y], s=120, c=RED, edgecolors="white", linewidths=1.5, zorder=5)
    ax.annotate(f"{NAME[pdb]}  ({pdb})",
                (x, y), textcoords="offset points", xytext=off, fontsize=10,
                color=INK, va="center",
                arrowprops=dict(arrowstyle="-", color="#CBD0D6", lw=0.8,
                                shrinkA=0, shrinkB=6))
seen = set()
for r, c_, f_ in zip(rows, cav, fid):
    pdb = r["target"][:4].upper()
    if pdb in RELATIVE and pdb not in seen and not r["target"].startswith("af_"):
        seen.add(pdb)
        ax.scatter([f_], [c_], s=90, c=DGREY, edgecolors="white",
                   linewidths=1.4, zorder=5)
        ax.annotate(f"{NAME[pdb]}  ({pdb})", (f_, c_), textcoords="offset points",
                    xytext=RELATIVE[pdb], fontsize=10, color=INK, va="center",
                    ha="right" if RELATIVE[pdb][0] < 0 else "left",
                    arrowprops=dict(arrowstyle="-", color="#CBD0D6", lw=0.8,
                                    shrinkA=0, shrinkB=6))

ax.set_xlabel("Sequence Identity to PYR1  (%)", fontsize=12, color=INK)
ax.set_ylabel("Cavity Volume  (Å$^3$)", fontsize=12, color=INK)
ax.set_title("Cavity volume across 147 helix-grip homologs",
             fontsize=14, color=INK, pad=14, loc="left")
ax.set_xlim(0, 60); ax.set_ylim(0, 700)
ax.grid(False)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_color("#C8CDD4")
ax.spines["bottom"].set_bounds(0, 60)
ax.spines["left"].set_bounds(0, 700)
ax.tick_params(colors=INK, labelsize=11)

h = [plt.Line2D([], [], marker="o", ls="", mfc=GREY, mec="white", ms=8),
     plt.Line2D([], [], marker="^", ls="", mfc="none", mec=GREY, ms=8)]
ax.legend(h, [f"Experimental  (n = {int((~isaf).sum())})",
              f"AlphaFold  (n = {int(isaf.sum())})"],
          loc="upper right", frameon=False, fontsize=10, labelcolor=INK,
          handletextpad=0.6)

fig.tight_layout()
out = f"{ROOT}/results/figures/graft_categorical.png"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
print(f"wrote {out}")
hi = fid > 25
print(f"  above 25% identity: n={int(hi.sum())}, largest cavity {cav[hi].max():.0f} Å³")
print(f"  above 300 Å³: n={int((cav>300).sum())}, median identity {np.median(fid[cav>300]):.0f}%")
