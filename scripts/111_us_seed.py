#!/usr/bin/env python
r"""
111_us_seed.py -- seed umbrella-sampling windows from the existing trajectories.

THE CALIBRATION THIS SERVES (README 59b)
§24 established that UNBIASED MD cannot cross the gate barrier: neither state
converts once in 1.8 us. So the open/closed free-energy difference -- the quantity
that separates a signalling switch from a constitutively closed protein (§58a) --
has never been measurable here. Umbrella sampling returns it as a PMF.

Before that can be trusted on a designed pocket it has to reproduce a case whose
answer is known:

    WT + ABA   should favour CLOSED
    WT apo     should favour OPEN

That is the §13e discipline (calibrate on a known anchor first) applied to a new
method, and it is a genuine test: if the PMFs do not separate in the right
direction, the scheme is wrong and no designed-pocket number from it means
anything.

THE REACTION COORDINATE, chosen by measurement rather than assumption
Every gate-latch CA pair was scored on how far it moves between the open and
closed crystal references. The best is P88 CA - R116 CA:

    crystal closed  7.64 A        crystal open  17.52 A
    S2 holo-closed  6.06 +/- 0.75     S1 apo-open   16.90 +/- 1.13
    S9 apo-closed   6.08 +/- 0.43     S10 holo-open 16.3-16.6

Complete separation, no overlap between the basins over 5 x 300 ns. Note S9 and S2
sit at the SAME value (6.08 vs 6.06), i.e. the closed state is geometrically
identical with and without ligand -- which is exactly why a stability measurement
cannot tell them apart and a free-energy one is needed.

A distance is used rather than an RMSD because Amber restrains distances natively
through &rst; an RMSD coordinate would need a plugin and adds a failure mode.

SEEDING. The existing factorial already samples nearly the whole range, so windows
are seeded from REAL equilibrated frames rather than from a steered pull, which
avoids dragging non-equilibrium structures into the windows. Holo lacks only
10.0-11.0 A and apo 8.5-12.5 A; those are seeded from the nearest available frame
and the restraint closes the gap during equilibration.
"""
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "data", "md191")
OUT = os.path.join(ROOT, "data", "umbrella")
CPPTRAJ = "cpptraj"

WINDOWS = [round(x, 1) for x in np.arange(5.0, 20.5, 0.5)]
SYSTEMS = {
    "holo": {"top": f"{MD}/S2_holo_closed/system.prmtop",
             "trajs": [("S2_holo_closed", "rep0"), ("S10_holo_open", "rep0"),
                       ("S10_holo_open", "rep2")]},
    "apo":  {"top": f"{MD}/S9_apo_closed/system.prmtop",
             "trajs": [("S9_apo_closed", "rep0"), ("S1_apo_open", "rep1")]},
}


def coord(sysname, rep):
    f = f"/tmp/rc_{sysname}_{rep}.dat"
    if not os.path.exists(f):
        raise SystemExit(f"{f} missing -- measure the coordinate first")
    return np.array([float(l.split()[1]) for l in open(f) if not l.startswith("#")])


def main():
    os.makedirs(OUT, exist_ok=True)
    for tag, cfg in SYSTEMS.items():
        # NOTE the two arms use DIFFERENT topologies (holo carries A8S), so a seed
        # frame may only be taken from a trajectory built on the same topology.
        series = []
        for s, r in cfg["trajs"]:
            v = coord(s, r)
            for i, x in enumerate(v):
                if i >= 3000:                       # drop 30 ns equilibration
                    series.append((x, s, r, i + 1)) # cpptraj frames are 1-indexed
        series = np.array(series, dtype=object)
        vals = np.array([s[0] for s in series], dtype=float)
        print(f"\n=== {tag}: {len(series)} candidate frames, "
              f"range {vals.min():.1f}-{vals.max():.1f} A")
        for w in WINDOWS:
            d = os.path.join(OUT, tag, f"w{w:04.1f}")
            os.makedirs(d, exist_ok=True)
            rst = os.path.join(d, "seed.rst7")
            if os.path.exists(rst):
                print(f"  {tag} w{w:04.1f}: seed exists")
                continue
            k = int(np.argmin(np.abs(vals - w)))
            x, s, r, frame = series[k]
            gap = abs(float(x) - w)
            inp = os.path.join(d, "seed.in")
            with open(inp, "w") as fh:
                fh.write(f"parm {MD}/{s}/system.prmtop\n")
                fh.write(f"trajin {MD}/{s}/{r}/prod.nc {frame} {frame} 1\n")
                fh.write(f"trajout {rst} restart\ngo\nquit\n")
            subprocess.run([CPPTRAJ, "-i", inp], capture_output=True, text=True)
            ok = os.path.exists(rst) or os.path.exists(rst + f".{frame}")
            # cpptraj may append the frame number for single-frame restart writes
            if not os.path.exists(rst) and os.path.exists(rst + f".{frame}"):
                os.rename(rst + f".{frame}", rst)
            flag = "" if gap < 0.25 else f"  <-- seed {gap:.2f} A away, restraint closes it"
            print(f"  {tag} w{w:04.1f}: {s}/{r} frame {frame} at {float(x):.2f} A"
                  f"{' OK' if os.path.exists(rst) else ' FAILED'}{flag}")


if __name__ == "__main__":
    main()
