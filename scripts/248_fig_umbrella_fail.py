#!/usr/bin/env python
r"""
248_fig_umbrella_fail.py -- SLIDE 7 FIGURE: the umbrella PMF is set by the pull
direction, not by the physics.

THE TEST. Umbrella sampling along P88 CA - R116 CA was run as a CALIBRATION: it
had to reproduce a known answer before being trusted on a designed pocket. The
pre-registered pass required THREE things at once -- holo favours closed
(positive), apo favours open (negative), and ddG(holo - apo) positive, i.e. ABA
stabilising the closed state.

WHAT HAPPENED. 62 windows, 1.24 us of sampling, WHAM converged. And:
  * holo +5.02 and apo +6.92 -- BOTH favour closed, so the apo criterion fails
  * ddG = -1.90 where a pass needs positive
  * and the diagnostic that explains it: seeding the same windows from the
    OPEN side vs the CLOSED side gives answers that differ by 3.22 kcal/mol on
    average and FLIP SIGN (+4.61 vs -3.29 over the same interval)

A converged PMF is seed-independent. The half-split drift -- the honest internal
error -- is 0.22 (holo) and 1.05 (apo), so the forward/reverse disagreement is
roughly 3x to 15x larger than the run's own error bar. That is the whole point of
the figure: the disagreement is not noise, it dwarfs the noise.

⚠ AND IT IS NOT A SAMPLING-LENGTH PROBLEM. A 15x miss with a sign flip is not
closed by extending 20 ns windows. It is protocol. My own 127 repair destroyed
the two-directional seeding that made the diagnosis possible, which is why the
hysteresis had to be recomputed before the arm could be closed.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "figures", "umbrella_hysteresis.png")
RED, BLUE, GREY = "#D96C6C", "#3D7EBF", "#BFC4CC"
#: from results/us_pmf_rev.log
ARMS = [("WT + ABA", 4.61, -3.29, 3.22, 6.57, 0.22),
        ("WT apo",   3.14, -2.40, 2.74, 6.49, 1.05)]


def main():
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.6),
                                  gridspec_kw=dict(width_ratios=[1.35, 1.0],
                                                   wspace=0.32))
    y = np.arange(len(ARMS))[::-1]
    for i, (lbl, fwd, rev, mean, mx, drift) in zip(y, ARMS):
        ax.plot([rev, fwd], [i, i], color=GREY, lw=2.4, zorder=1,
                solid_capstyle="round")
        ax.scatter([fwd], [i], s=95, color=BLUE, zorder=3,
                   label="Forward (seeded from CLOSED)" if i == y[0] else None)
        ax.scatter([rev], [i], s=95, color=RED, zorder=3,
                   label="Reverse (seeded from OPEN)" if i == y[0] else None)
        ax.text((fwd + rev) / 2, i + 0.17, f"{abs(fwd-rev):.1f} kcal/mol apart",
                ha="center", va="bottom", fontsize=9, color="0.35")
    ax.axvline(0, color="black", lw=1.0, ls=(0, (4, 3)), zorder=0)
    ax.set_yticks(y); ax.set_yticklabels([a[0] for a in ARMS], fontsize=10)
    ax.set_ylim(-0.6, len(ARMS) - 0.25)
    ax.set_xlabel("ΔG over the same interval (kcal/mol)")
    ax.set_title("The same windows, seeded two ways, disagree in SIGN",
                 fontsize=11, loc="left", pad=26)
    ax.legend(frameon=False, fontsize=9, ncol=2, loc="lower left",
              bbox_to_anchor=(0.0, 1.0))

    # right panel: disagreement vs the run's OWN error bar
    w = 0.34
    x = np.arange(len(ARMS))
    d = [a[3] for a in ARMS]; e = [a[5] for a in ARMS]
    ax2.bar(x - w/2, d, w, color=RED, label="Forward vs reverse")
    ax2.bar(x + w/2, e, w, color=GREY, label="Half-split drift (own error)")
    for xi, (dv, ev) in enumerate(zip(d, e)):
        ax2.text(xi, max(dv, ev) + 0.12, f"{dv/ev:.0f}× larger",
                 ha="center", fontsize=9, color="0.35")
    ax2.set_xticks(x); ax2.set_xticklabels([a[0] for a in ARMS], fontsize=10)
    ax2.set_ylabel("kcal/mol")
    ax2.set_title("The disagreement dwarfs the noise", fontsize=11,
                  loc="left", pad=40)
    # ⚠ two-entry legend needs its own row or it lands on the title
    ax2.legend(frameon=False, fontsize=9, ncol=2, loc="lower left",
               bbox_to_anchor=(0.0, 1.0), columnspacing=1.2, handlelength=1.4)
    ax2.set_ylim(0, max(max(d), max(e)) * 1.25)
    for a in (ax, ax2):
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        a.grid(False)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print("PRE-REGISTERED PASS required all three:")
    print("   holo favours CLOSED (positive)      got +5.02   PASS")
    print("   apo  favours OPEN   (negative)      got +6.92   FAIL")
    print("   ddG(holo-apo) POSITIVE              got -1.90   FAIL")
    print("\nhysteresis, same windows seeded two ways:")
    for lbl, fwd, rev, mean, mx, drift in ARMS:
        print(f"   {lbl:<10} forward {fwd:+.2f}  reverse {rev:+.2f}  "
              f"mean |diff| {mean:.2f}  max {mx:.2f}  own drift {drift:.2f}"
              f"  -> {mean/drift:.0f}x its own error bar")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
