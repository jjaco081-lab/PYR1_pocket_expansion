#!/usr/bin/env python
"""81b_frames.py -- how many production frames actually exist on disk.

Resumability is decided by MEASUREMENT, not by a step counter we maintain. A
counter parsed from .out files drifts the moment a segment dies mid-write, and
README 24e already cost 11.6 ns of duplicated frames to exactly that class of
bookkeeping error: the restart cadence and the frame cadence disagreed, so an
interrupted run re-simulated from the last restart while the trajectory kept the
frames written past it.

Counting frames in the NetCDF headers cannot drift, because it asks the files.
"""
import glob
import re
import subprocess
import sys


def count(f):
    try:
        import netCDF4
        return len(netCDF4.Dataset(f).dimensions["frame"])
    except Exception:
        pass
    try:
        o = subprocess.run(["ncdump", "-h", f], capture_output=True, text=True).stdout
        m = re.search(r"frame = UNLIMITED ; // \((\d+) currently\)", o)
        return int(m.group(1)) if m else 0
    except Exception:
        return 0


if __name__ == "__main__":
    pat = sys.argv[1] if len(sys.argv) > 1 else "prod_seg*.nc"
    print(sum(count(f) for f in sorted(glob.glob(pat))))
