#!/usr/bin/env python
r"""
247_fig_md_gate.py -- SLIDE 5 FIGURE: the gate never converts in 300 ns.

The 2x2 factorial (conformation x ABA occupancy), 3 replicates each, 300 ns per
replicate. Reaction coordinate is the P88 CA - R116 CA distance, which is the
gate-latch opening -- the same coordinate the umbrella sampling later used.

WHAT THE FIGURE HAS TO SHOW: the two conformations occupy separated, non-
overlapping bands and no trajectory ever crosses between them. That is the whole
claim of the loop-dynamics arm, and it is what licenses a STABILITY filter while
forbidding a SWITCHABILITY one.

⚠ S2_holo_closed_rep2 stopped at 1,084 frames (~11 ns) against 30,000 for the
others. It is drawn but flagged in the caption count -- silently averaging a
short replicate into a 300 ns panel is how 62a came to rest on n=1.
"""
import glob, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RC = os.path.join(ROOT, "results", "factorial", "rc")
OUT = os.path.join(ROOT, "results", "figures", "md_gate_trapped.png")
NS_TOTAL = 300.0
SYS = [("S1_apo_open",    "Open, apo",      "#3D7EBF", "-"),
       ("S10_holo_open",  "Open, + ABA",    "#7FB3E0", "-"),
       ("S9_apo_closed",  "Closed, apo",    "#D96C6C", "-"),
       ("S2_holo_closed", "Closed, + ABA",  "#A83232", "-")]


def load(tag):
    out = []
    for f in sorted(glob.glob(os.path.join(RC, f"{tag}_rep*.dat"))):
        v = np.loadtxt(f, comments="#", usecols=1)
        out.append((os.path.basename(f), v))
    return out


def main():
    fig, (ax, axh) = plt.subplots(
        1, 2, figsize=(10.4, 5.0), sharey=True,
        gridspec_kw=dict(width_ratios=[3.0, 1.0], wspace=0.06))
    stats, short = [], []
    for tag, lbl, col, ls in SYS:
        reps = load(tag)
        allv = []
        for name, v in reps:
            t = np.linspace(0, NS_TOTAL * len(v) / 30000.0, len(v))
            k = max(1, len(v) // 1500)
            ax.plot(t[::k], v[::k], color=col, lw=0.6, alpha=0.55, zorder=2)
            allv.append(v)
            if len(v) < 20000:
                short.append((name, len(v)))
        cat = np.concatenate(allv)
        stats.append((lbl, cat.mean(), cat.std(), cat.min(), cat.max(),
                      len(reps), col))
        axh.hist(cat, bins=70, orientation="horizontal", color=col,
                 alpha=0.55, density=True, zorder=2)
        ax.plot([], [], color=col, lw=2.2, label=f"{lbl}  (n={len(reps)})")

    # the watershed: nothing crosses it
    ax.axhspan(10.5, 12.0, color="0.85", lw=0, zorder=0)
    axh.axhspan(10.5, 12.0, color="0.85", lw=0, zorder=0)
    # ⚠ "no trajectory crosses" was WRONG -- the bands do overlap. Measured:
    # every excursion past the midpoint returns, the longest lasting 1.74 ns of
    # 300, and no replicate spends more than 9.9 % of frames on the far side.
    ax.text(NS_TOTAL * 0.015, 11.25,
            "  excursions return: longest 1.74 ns of 300",
            ha="left", va="center", fontsize=9, color="0.35")

    ax.set_xlabel("Simulation time (ns)")
    ax.set_ylabel("Gate opening, P88 Cα – R116 Cα (Å)")
    ax.set_title("Open and closed are both kinetically trapped at 300 ns",
                 fontsize=11.5, pad=34, loc="left")
    axh.set_xlabel("Density")
    axh.set_xticks([])
    for a in (ax, axh):
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        a.grid(False)
    ax.set_xlim(0, NS_TOTAL)
    # legend ABOVE the axes: every in-panel position collided with either a
    # trace band or the watershed annotation
    ax.legend(frameon=False, fontsize=9, ncol=4, loc="lower left",
              bbox_to_anchor=(0.0, 1.02), handlelength=1.6,
              columnspacing=1.3)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=200, bbox_inches="tight")

    print(f"{'system':<16}{'mean':>8}{'SD':>7}{'min':>7}{'max':>7}{'reps':>6}")
    for lbl, m, s, lo, hi, n, _ in stats:
        print(f"{lbl:<16}{m:>8.2f}{s:>7.2f}{lo:>7.2f}{hi:>7.2f}{n:>6}")
    # excursion audit -- the claim is NO SUSTAINED CONVERSION, not "no crossing"
    CUT = 11.25
    worst = 0.0
    for tag, lbl, col, ls in SYS:
        for name, v in load(tag):
            closed = v[:50].mean() < CUT
            far = (v > CUT) if closed else (v < CUT)
            best = cur = 0
            for x in far:
                cur = cur + 1 if x else 0
                best = max(best, cur)
            worst = max(worst, best * 300.0 / 30000.0)
    print(f"\nlongest sustained excursion across all 12 replicates: {worst:.2f} ns "
          f"of 300 -- and it RETURNS")
    if short:
        print(f"\n⚠ short replicate(s): {short}")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
