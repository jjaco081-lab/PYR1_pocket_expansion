#!/usr/bin/env python
r"""
239_build_chimera.py -- write the actual PYR1/donor chimera, superimposed, with
ONLY the residues each parent contributes.

Jannis: "give me a correctly superimposed structure for the first two with each
only showing the parts that would make the protein."

WHAT IS IN THE FILE. For each donor, `scan_junctions` has already assigned every
aligned position a parent: PYR1 keeps the gate (85-89), the latch (115-117), the
tunnel positions (81, 83, 87, 89, 92, 117, 159), the HAB1 interface blocks
(58-65, 81-92, 111-121, 146-168) plus a 2-residue flank, and slides each boundary
to the position of least backbone deviation. The donor supplies everything else,
which is where the bigger pocket comes from. This writes exactly that:

    chain A = the PYR1-derived residues, in PYR1's own frame
    chain B = the DONOR-derived residues, transformed INTO PYR1's frame
    chain J = one pseudo-atom per junction, B-factor = its CA deviation in A

so the two halves are already in register and nothing that would be discarded is
drawn. ⚠ This is a SPLICE MAP, not a built model: no loop closure, no repacking,
no minimisation. Residues at the junctions are left exactly where the two parents
put them, so any gap you see is a real gap that inpainting would have to close.
That is the point -- Jannis: "I think 7/8 might be fine because we could use
inpainting or some other tools for some cases, although not all."

⚠ THE TRANSFORM. 214's `kabsch(P, Q)` returns R, mean(P), mean(Q) such that
PYR1 -> donor is (X - Pm) @ R + Qm. Here the donor is brought to PYR1 instead,
which is the inverse: (Y - Qm) @ R.T + Pm. Getting this backwards would still
produce a plausible-looking file, so it is asserted below: after transforming,
the mapped donor CAs must sit within the alignment RMSD of their PYR1 partners.
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
imp = __import__("importlib").import_module
m214, m237 = imp("214_schema_foldseek"), imp("237_junctions_all")
m240 = imp("240_junction_greedy")
RAW = os.path.join(ROOT, "results", "homolog_cavities", "raw")
OUT = os.path.join(ROOT, "results", "chimera")
AA3 = m214.AA3
ONE2THREE = {v: k for k, v in AA3.items()}


def atoms_of(path, chain, want):
    """[(resnum, atomname, element, resname, xyz)] for residues in `want`."""
    out = []
    if path.endswith(".cif"):
        cols, rows, inl = {}, [], False
        for line in open(path):
            if line.startswith("_atom_site."):
                cols[line.strip().split(".")[1]] = len(cols); inl = True; continue
            if inl:
                if line.startswith(("#", "loop_", "_")):
                    if rows: break
                    continue
                x = line.split()
                if len(x) >= len(cols): rows.append(x)
        g = lambda r, k: r[cols[k]]                              # noqa: E731
        for r in rows:
            if g(r, "group_PDB") != "ATOM" or g(r, "type_symbol") == "H":
                continue
            if chain and g(r, "auth_asym_id") != chain:
                continue
            try: n = int(g(r, "auth_seq_id"))
            except ValueError: continue
            if n not in want: continue
            out.append((n, g(r, "label_atom_id"), g(r, "type_symbol"),
                        g(r, "label_comp_id"),
                        np.array([float(g(r, "Cartn_x")), float(g(r, "Cartn_y")),
                                  float(g(r, "Cartn_z"))])))
    else:
        for l in open(path):
            if l.startswith("ENDMDL"): break
            if not l.startswith("ATOM"): continue
            if chain and l[21] != chain: continue
            e = (l[76:78].strip() or l[12:16].strip()[0]).upper()
            if e == "H": continue
            try: n = int(l[22:26])
            except ValueError: continue
            if n not in want: continue
            out.append((n, l[12:16].strip(), e, l[17:20].strip(),
                        np.array([float(l[30:38]), float(l[38:46]),
                                  float(l[46:54])])))
    return out


def w(fh, i, name, resn, ch, num, xyz, el, b=0.0):
    """⚠ PDB columns are fixed-width and the ALTLOC column (17) is easy to drop.
    Without it resName and chainID each shift one left, the chain field comes out
    BLANK, and every downstream chain selection silently matches nothing -- which
    is exactly what happened on the first build here.
    ATOM(6) serial(5) space(1) name(4) altLoc(1) resName(3) space(1) chainID(1)
    resSeq(4) = 26 characters before the insertion code."""
    nm = f" {name[:3]}" if len(name) < 4 else name[:4]
    fh.write(f"ATOM  {i%99999:5d} {nm:<4s} {resn:>3s} {ch}{num%9999:4d}    "
             f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00{b:6.2f}"
             f"          {el:>2s}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+", required=True)
    ap.add_argument("--flank", type=int, default=2)
    #: greedy placement keeps the donor everything past the FIRST
    #: acceptable junction instead of the globally best one -- 214's
    #: minimiser gave PYR1 126 of 179 residues when only 70 are required.
    ap.add_argument("--greedy", action="store_true")
    ap.add_argument("--cut", type=float, default=1.5)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    pca, pheavy, pseq = m214.read_cif(os.path.join(ROOT, "data", "3QN1.cif"), "A")
    pres = sorted(pca)
    aln = {}
    for line in open(f"{ROOT}/results/foldseek/fs_hits.tsv"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 11 and f[1] not in aln:
            aln[f[1]] = f
    K0 = m214.keep_set(a.flank)

    for t in a.targets:
        if t not in aln:
            print(f"!! {t}: no alignment"); continue
        path = None
        for ext in (".cif", ".pdb"):
            p = os.path.join(RAW, t + ext)
            if os.path.exists(p): path = p; break
        if path is None:
            print(f"!! {t}: no structure"); continue
        dchain = (t[4] if (path.endswith(".cif") and len(t) > 4) else None)
        dca, dheavy, dseq = m237.read_any(path, dchain)
        dres = sorted(dca)
        f = aln[t]
        m = m214.map_from_alignment(f[9], f[10], int(f[5]), int(f[7]), pres, dres)
        m = {k: v for k, v in m.items() if k in pca and v in dca}
        ks = sorted(m)
        P = np.array([pca[k] for k in ks]); Q = np.array([dca[m[k]] for k in ks])
        R, pm, qm = m214.kabsch(P, Q)
        # donor -> PYR1 frame is the INVERSE of 214's PYR1 -> donor
        gx = lambda Y: (Y - qm) @ R.T + pm                        # noqa: E731
        back = np.array([gx(dca[m[k]][None, :])[0] for k in ks])
        rms = float(np.sqrt(np.mean(np.sum((back - P) ** 2, axis=1))))
        fwd = np.array([(pca[k] - pm) @ R + qm for k in ks])
        rms_f = float(np.sqrt(np.mean(np.sum((fwd - Q) ** 2, axis=1))))
        assert abs(rms - rms_f) < 1e-3, (
            f"{t}: inverse transform is wrong ({rms:.2f} vs {rms_f:.2f})")
        dev = {k: float(np.linalg.norm(back[i] - P[i])) for i, k in enumerate(ks)}
        K = K0 & set(ks)
        parent, J = (m240.scan_greedy(ks, K, dev, a.cut) if a.greedy
                     else m214.scan_junctions(ks, K, dev))
        keep_p = {k for k in ks if parent[k] == "PYR1"}
        keep_d = {m[k] for k in ks if parent[k] == "DON"}
        pa = atoms_of(os.path.join(ROOT, "data", "3QN1.cif"), "A", keep_p)
        da = atoms_of(path, dchain, keep_d)
        out = os.path.join(OUT, f"chimera_{t}{'_greedy' if a.greedy else ''}.pdb")
        with open(out, "w") as fh:
            fh.write(f"REMARK  PYR1 / {t} splice map, flank {a.flank}\n")
            fh.write(f"REMARK  chain A = PYR1-derived ({len(keep_p)} res): gate, "
                     f"latch, tunnel, HAB1 interface\n")
            fh.write(f"REMARK  chain B = donor-derived ({len(keep_d)} res): the "
                     f"pocket-forming backbone\n")
            fh.write(f"REMARK  chain J = junctions, B-factor = CA deviation (A)\n")
            fh.write(f"REMARK  alignment RMSD {rms:.2f} A over {len(ks)} residues\n")
            fh.write("REMARK  NOT a built model: no loop closure, no repacking\n")
            i = 0
            for n, nm, el, rn, x in pa:
                i += 1; w(fh, i, nm, rn, "A", n, x, el)
            for n, nm, el, rn, x in da:
                i += 1; w(fh, i, nm, rn, "B", n, gx(x[None, :])[0], el)
            for jn, (pos, d) in enumerate(J, 1):
                i += 1
                w(fh, i, "O", "HOH", "J", jn, pca[pos], "O", b=round(d, 2))
            fh.write("END\n")
        clean = sum(1 for _, d in J if d < 1.5)
        print(f"{t:<12} PYR1 {len(keep_p):>3} res + donor {len(keep_d):>3} res   "
              f"RMSD {rms:.2f}   junctions {clean}/{len(J)} clean   -> {out}")
        for jn, (pos, d) in enumerate(J, 1):
            mark = "" if d < 1.5 else "   <-- needs inpainting"
            print(f"     junction {jn} at PYR1 {pos:>4}  deviation {d:5.2f} A{mark}")
        # a pml that shows only what survives
        pml = os.path.join(OUT, f"chimera_{t}{'_greedy' if a.greedy else ''}.pml")
        with open(pml, "w") as fh:
            fh.write(f"load {os.path.basename(out)}, chim\n")
            fh.write("hide everything\nshow cartoon, chim and chain A+B\n")
            fh.write("color skyblue, chim and chain A\n")
            fh.write("color firebrick, chim and chain B\n")
            fh.write("show spheres, chim and chain J\n")
            fh.write("color yellow, chim and chain J\nset sphere_scale, 0.6, chim and chain J\n")
            fh.write("set cartoon_transparency, 0.1\nbg_color white\n")
            fh.write("show sticks, chim and chain A+B and not name N+C+O\n")
            fh.write("set stick_radius, 0.12\n")
            fh.write("orient chim and chain A+B\n")
            fh.write('print "BLUE = PYR1 (gate, latch, tunnel, HAB1 interface)"\n')
            fh.write('print "RED  = donor (the pocket)"\n')
            fh.write('print "YELLOW spheres = junctions; B-factor is the CA deviation"\n')
    return 0


if __name__ == "__main__":
    sys.exit(main())
