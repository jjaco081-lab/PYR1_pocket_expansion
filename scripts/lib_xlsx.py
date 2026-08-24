#!/usr/bin/env python
r"""
lib_xlsx.py -- minimal dependency-free .xlsx reader.

`openpyxl` is not installed in any conda env on this cluster, and 62 already
depends on it, so the Tian datasets could not be re-read. An .xlsx file is a zip
of XML, so the sheet can be parsed with the standard library alone. This reads
values only -- no formulas, no formatting -- which is all these datasets need.

Two traps this handles that a naive reader gets wrong:
  * strings live in a shared-strings table and cells reference them by INDEX
    (t="s"), so a cell's literal <v> is an integer pointing elsewhere;
  * empty cells are OMITTED from the XML rather than written blank, so rows must
    be reconstructed by the cell's column REFERENCE (A1, B1, ...) and not by the
    order in which cells happen to appear. Ignoring that silently shifts every
    column after the first gap.
"""
import re
import zipfile
import xml.etree.ElementTree as ET

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def _col_index(ref):
    """'A'->0, 'B'->1, ... 'AA'->26."""
    m = re.match(r"([A-Z]+)", ref)
    n = 0
    for ch in m.group(1):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def sheet_rows(path, sheet=0):
    """Yield each row as a list of strings ('' for empty cells)."""
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall(f"{NS}si"):
                shared.append("".join(t.text or "" for t in si.iter(f"{NS}t")))
        names = sorted(n for n in z.namelist()
                       if re.match(r"xl/worksheets/sheet\d+\.xml$", n))
        root = ET.fromstring(z.read(names[sheet]))
        out = []
        for row in root.iter(f"{NS}row"):
            cells = {}
            for c in row.findall(f"{NS}c"):
                ref = c.get("r") or ""
                if not ref:
                    continue
                i = _col_index(ref)
                t = c.get("t")
                if t == "inlineStr":
                    v = "".join(x.text or "" for x in c.iter(f"{NS}t"))
                else:
                    ve = c.find(f"{NS}v")
                    v = ve.text if ve is not None else ""
                    if t == "s" and v not in (None, ""):
                        v = shared[int(v)]
                cells[i] = (v or "").strip()
            if cells:
                out.append([cells.get(i, "") for i in range(max(cells) + 1)])
    return out


def table(path, sheet=0):
    """First row is the header; returns (header, list-of-dicts)."""
    rows = sheet_rows(path, sheet)
    if not rows:
        return [], []
    hdr = rows[0]
    recs = []
    for r in rows[1:]:
        if not any(x for x in r):
            continue
        r = r + [""] * (len(hdr) - len(r))
        recs.append(dict(zip(hdr, r)))
    return hdr, recs


if __name__ == "__main__":
    import sys
    h, recs = table(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    print(f"{len(recs)} rows, {len(h)} columns")
    print("header:", h)
    for r in recs[:3]:
        print({k: v for k, v in list(r.items())[:10]})
