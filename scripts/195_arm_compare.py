#!/usr/bin/env python
r"""
195_arm_compare.py -- compare RFd3 Stage-0 arms on the five selection criteria
fixed in the plan, one row per design.

Plan, Stage 0: "Select the winning arm on: fraction with an enclosed cavity,
largest-connected envelope vs PYR1's 1184, motif CA deviation (batch 1 achieved
0.1 A), chain breaks, and HAB1-contacting residues in the leading segment
(currently 46)."

Three things this does that 190 did not:

  * reads gzipped CIFs from the per-design directories 193_rfd3_gen.py writes,
    rather than a flat batch directory of plain .cif;
  * takes the motif correspondence from RFd3's OWN `diffused_index_map`
    (input residue -> output residue) instead of re-deriving it by counting,
    so the CA deviation cannot be silently computed against the wrong residues;
  * ASSERTS THE MOTIF RESIDUE IDENTITIES. The map is used to look up each output
    residue and its three-letter code is checked against the PYR1 motif sequence.
    A number-only check would pass even if the whole segment were off by one --
    which is exactly the -2 auth/label offset that nearly scaffolded the wrong
    residues in the first place.

Usage:
  python 195_arm_compare.py results/rfd3/arm0a_s2026 results/rfd3/arm0b_s2027
"""
import glob, gzip, json, os, sys
import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from importlib import import_module                                 # noqa: E402
m143 = import_module("143_pocket_backbone_shape")
m151 = import_module("151_cavity_bottleneck")
m192 = import_module("192_filter_cascade")

#: the scaffolded motif, from 193_rfd3_gen.py. One-letter, in contig order.
MOTIF = [("A", 34, 40, "HAQRIHA"), ("A", 58, 65, "YKHFIKSC"),
         ("A", 81, 92, "VIVISGLPANTS"), ("A", 111, 121, "IGGEHRLTNYK"),
         ("A", 146, 168, "DMPEGNSEDDTRMFADTVVKLNL")]
ONE = {v: k for k, v in m143.THREE.items()} if isinstance(
    getattr(m143, "THREE", {}), dict) else {}
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
BB = {"N", "CA", "C", "O", "OXT"}
LEAD = 25          # "leading segment" of the design chain, per the plan
PYR1_ENVELOPE, PYR1_RMAX, PYR1_RG = 1184.0, 3.21, 15.0


def read_cif_gz(path):
    """[(asym, seq, atom, elem, xyz, comp)] -- heavy atoms only."""
    op = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    cols, rows, inl = {}, [], False
    with op as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols[line.strip().split(".")[1]] = len(cols); inl = True; continue
            if inl:
                if line.startswith(("#", "loop_", "_")):
                    if rows:
                        break
                    continue
                x = line.split()
                if len(x) >= len(cols):
                    rows.append(x)
    out = []
    for r in rows:
        if r[cols["type_symbol"]] == "H":
            continue
        out.append((r[cols["label_asym_id"]], int(r[cols["label_seq_id"]]),
                    r[cols["label_atom_id"]], r[cols["type_symbol"]],
                    np.array([float(r[cols[k]]) for k in
                              ("Cartn_x", "Cartn_y", "Cartn_z")]),
                    r[cols["label_comp_id"]]))
    return out


def kabsch_rmsd(P, Q):
    """max per-atom deviation and RMSD after optimal superposition of P onto Q."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1, 1, d]) @ Wt
    dev = np.linalg.norm(Pc @ R - Qc, axis=1)
    return float(dev.max()), float(np.sqrt((dev ** 2).mean()))


def ref_motif_ca():
    """{('A',34): (xyz, 'H'), ...} from the auth-numbered complex.

    m143.read_pdb drops the chain ID, and the reference is a two-chain complex,
    so this reads it here. Chain-blind lookup would silently mix HAB1 numbering
    into the PYR1 motif.
    """
    ca = {}
    for line in open(os.path.join(ROOT, "data", "3QN1_complex_auth.pdb")):
        if line[:4] != "ATOM" or line[12:16].strip() != "CA":
            continue
        if line[16] not in (" ", "A"):
            continue
        ca[(line[21], int(line[22:26]))] = (
            np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
            AA3.get(line[17:20].strip(), "?"))
    # assert the reference itself is what the contig claims, before using it
    for ch, lo, hi, seq in MOTIF:
        got = "".join(ca.get((ch, r), (None, "-"))[1] for r in range(lo, hi + 1))
        assert got == seq, f"reference {ch}{lo}-{hi} is {got}, contig says {seq}"
    return ca


def design_chain(atoms, imap):
    """The diffused chain, taken from RFd3's own map and then ASSERTED.

    ⚠ The first version of this took the chain with the most residues, and it
    was WRONG on 6 of 20 designs: HAB1 block B283-461 is 179 residues, and any
    design shorter than that won it. Those rows silently measured HAB1's
    cavities and HAB1's radius of gyration. The motif-identity check is what
    exposed it -- the six bad rows were exactly the ones reporting 179 residues
    and a scrambled motif sequence.

    `diffused_index_map` covers HAB1 as well as PYR1, so the design chain is the
    one the PYR1 motif keys map into, not the one with the most entries. The
    result is then cross-checked against the HAB1 block sequences: if the chosen
    chain reproduces a context block verbatim, we picked context, not design.
    """
    ch = {imap[f"A{r}"][0] for c, lo, hi, _ in MOTIF for r in (lo, hi)
          if f"A{r}" in imap}
    assert len(ch) == 1, f"PYR1 motif maps into {ch or 'nothing'}, expected one chain"
    dc = ch.pop()
    seq = {}
    for a in atoms:
        if a[0] == dc:
            seq[a[1]] = AA3.get(a[5], "?")
    s = "".join(seq[k] for k in sorted(seq))
    for blk in hab1_blocks():
        assert s != blk, f"chain {dc} IS an HAB1 context block, not the design"
    return dc


_HAB1 = None


def hab1_blocks():
    """one-letter sequences of the four fixed HAB1 blocks, from the reference."""
    global _HAB1
    if _HAB1 is None:
        seq = {}
        for line in open(os.path.join(ROOT, "data", "3QN1_complex_auth.pdb")):
            if line[:4] == "ATOM" and line[12:16].strip() == "CA" and line[21] == "B":
                seq[int(line[22:26])] = AA3.get(line[17:20].strip(), "?")
        _HAB1 = ["".join(seq[r] for r in range(lo, hi + 1) if r in seq)
                 for lo, hi in ((185, 221), (232, 270), (283, 461), (466, 505))]
    return _HAB1


def chain_breaks(ca_by_seq):
    seqs = sorted(ca_by_seq)
    nb = 0
    for i in range(len(seqs) - 1):
        if seqs[i + 1] == seqs[i] + 1:
            if np.linalg.norm(ca_by_seq[seqs[i + 1]] - ca_by_seq[seqs[i]]) > 4.5:
                nb += 1
    return nb


def measure(cifp, name):
    js = cifp.replace(".cif.gz", ".json").replace(".cif", ".json")
    if not os.path.exists(js):
        return None
    at = read_cif_gz(cifp)
    imap = json.load(open(js)).get("diffused_index_map", {})
    dc = design_chain(at, imap)
    dat = [a for a in at if a[0] == dc]
    other = np.array([a[4] for a in at if a[0] != dc])

    # --- motif: identity assertion, then CA deviation ------------------------
    comp = {}
    ca = {}
    for a in dat:
        comp[a[1]] = a[5]
        if a[2] == "CA":
            ca[a[1]] = a[4]
    ref = ref_motif_ca()
    P, Q, bad = [], [], []
    for ch, lo, hi, seq in MOTIF:
        for k, resi in enumerate(range(lo, hi + 1)):
            key = f"{ch}{resi}"
            if key not in imap:
                bad.append(f"{key} unmapped"); continue
            oi = int(imap[key][1:])
            got = AA3.get(comp.get(oi, "???"), "?")
            if got != seq[k]:
                bad.append(f"{key}->{imap[key]} is {got} want {seq[k]}")
                continue
            rk = ref.get((ch, resi))
            if rk is None or oi not in ca:
                continue
            P.append(ca[oi]); Q.append(rk[0])
    dev = kabsch_rmsd(np.array(P), np.array(Q)) if len(P) > 3 else (None, None)

    # --- leading segment HAB1 contacts --------------------------------------
    # ⚠ imap covers HAB1's fixed blocks too, and those start at output index 1,
    # so taking min() over ALL values put the first motif residue at 1 and made
    # the leading segment empty on every design. Restrict to the PYR1 motif.
    motif_out = sorted(int(v[1:]) for k, v in imap.items() if v[0] == dc)
    lead_end = min(motif_out) - 1 if motif_out else 0
    lead = [r for r in sorted(ca) if r <= min(lead_end, LEAD)]
    lcon = 0
    if len(other):
        for r in lead:
            xyz = np.array([a[4] for a in dat if a[1] == r])
            if len(xyz) and np.linalg.norm(
                    other[None, :, :] - xyz[:, None, :], axis=2).min() < 4.5:
                lcon += 1

    # --- cavity, envelope, Rg, spanning -------------------------------------
    tup = [(a[1], None, a[2], a[3], a[4]) for a in dat]
    r, wall, pts = m192.cavity_and_wall(tup)
    env = m192.envelope(tup, wall) if wall else None
    # FAST=1 skips the F5 merge test, which recomputes the whole cavity grid
    # once per wall residue (~24-30 per design) and dominates the runtime on any
    # design that HAS a cavity. The arm comparison does not turn on F5.
    span = ([] if os.environ.get("FAST") else
            m192.spanning_residue(tup, pts, wall)) if wall else []
    cav = np.array(list(ca.values()))
    rg = float(np.sqrt(((cav - cav.mean(0)) ** 2).sum(1).mean()))
    return dict(name=name, nres=len(ca), rg=rg,
                main=(r or {}).get("main", 0.0), rmax=(r or {}).get("r_max", 0.0),
                env=(env or {}).get("envelope", 0.0),
                ncomp=(env or {}).get("ncomp", 0),
                dw=(env or {}).get("dw", 0.0),
                dev_max=dev[0], dev_rms=dev[1], nmotif=len(P),
                breaks=chain_breaks(ca), lead_n=len(lead), lead_hab1=lcon,
                span=len(span), bad=bad)


def emit(m):
    dv = "n/a" if m["dev_max"] is None else f"{m['dev_max']:.2f}"
    flag = f"  <- +{m['env']-PYR1_ENVELOPE:.0f}" if m["env"] > PYR1_ENVELOPE else ""
    print(f"  {m['name']:<7}{m['nres']:>5}{m['rg']:>6.1f}{m['main']:>8.0f}"
          f"{m['rmax']:>7.2f}{m['env']:>10.0f}{m['ncomp']:>6}{dv:>9}"
          f"{m['breaks']:>5}{m['lead_hab1']:>4}/{m['lead_n']:<5}"
          f"{m['span']:>6}{flag}", flush=True)
    if m["bad"]:
        print(f"       ⚠ MOTIF IDENTITY: {'; '.join(m['bad'][:4])}"
              f"{' ...' if len(m['bad'])>4 else ''}", flush=True)


def main():
    arms = sys.argv[1:]
    if not arms:
        print(__doc__); return 1
    print(f"  PYR1 reference: envelope {PYR1_ENVELOPE:.0f} A^3, r_max {PYR1_RMAX}, "
          f"Rg {PYR1_RG}\n  batch 3 (best so far): 7/10 cavity, all 7 > 1184, best 1353\n")
    summary = {}
    for arm in arms:
        rows = []
        cifs = sorted(glob.glob(os.path.join(ROOT, arm, "d*", "*.cif.gz")))
        if not cifs:          # flat batch directories (batch01-03)
            cifs = sorted(glob.glob(os.path.join(ROOT, arm, "*.cif.gz")))
        print(f"== {arm}   n={len(cifs)}")
        print(f"  {'design':<7}{'res':>5}{'Rg':>6}{'cavity':>8}{'r_max':>7}"
              f"{'ENVELOPE':>10}{'ncomp':>6}{'motifDev':>9}{'brk':>5}"
              f"{'leadHAB1':>10}{'span':>6}", flush=True)
        for c in cifs:
            nm = os.path.basename(os.path.dirname(c))
            if not nm.startswith("d"):
                nm = "m" + os.path.basename(c).split("_model_")[1].split(".")[0]
            m = measure(c, nm)
            if m:
                rows.append(m)
                emit(m)      # stream per design: an arm takes ~10 min and a
                             # buffered table gives no progress signal at all
        cav = [m for m in rows if m["main"] > 0]
        big = [m for m in cav if m["env"] > PYR1_ENVELOPE]
        summary[arm] = dict(
            n=len(rows), cavity=len(cav), over_pyr1=len(big),
            best_env=max([m["env"] for m in rows], default=0.0),
            med_dev=float(np.median([m["dev_max"] for m in rows
                                     if m["dev_max"] is not None]) or 0),
            med_lead=float(np.median([m["lead_hab1"] for m in rows])),
            breaks=sum(m["breaks"] for m in rows),
            motif_fail=sum(1 for m in rows if m["bad"]))
        s = summary[arm]
        print(f"  -> cavity {s['cavity']}/{s['n']}, over PYR1 {s['over_pyr1']}/{s['n']}, "
              f"best {s['best_env']:.0f}, median motif dev {s['med_dev']:.2f} A, "
              f"median leading-HAB1 {s['med_lead']:.0f}, breaks {s['breaks']}, "
              f"motif-identity failures {s['motif_fail']}\n")
    out = os.path.join(ROOT, "results", "rfd3", "arm_compare.json")
    json.dump(summary, open(out, "w"), indent=1)
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
