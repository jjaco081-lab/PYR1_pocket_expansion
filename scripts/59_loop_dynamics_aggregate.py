#!/usr/bin/env python
"""
59_loop_dynamics_aggregate.py -- does open stay open, does closed stay closed, and
is there a separation we can filter pocket-expanded designs on?

Consumes the cpptraj output from 58b and answers three questions in order. The
order matters: question 3 is only worth asking if 1 and 2 come out clean.

  1. VALIDITY. Did the systems stay themselves? Core backbone equilibrated, ABA
     still in the pocket in the holo systems. If ABA left, "closed" means nothing
     and the system is void.
  2. STATE STABILITY. Over 300 ns does each system sit on its own side of the
     open/closed axis, and how many times does it cross?
  3. SEPARABILITY. Is any pre-registered observable actually usable as a filter,
     and with what error bar?

THE STATISTICS TRAP THIS SCRIPT IS BUILT AROUND
-----------------------------------------------
There are ~30,000 frames per replicate but they are NOT 30,000 independent
samples. Loop conformations are correlated over nanoseconds, so a t-test across
pooled frames would report absurd significance for any difference at all. This
project has already been bitten once by a discrimination statistic computed the
wrong way (README 15/pilot results).

So:
  - the INFERENTIAL unit is the REPLICATE (n=3 open, n=3 closed, n=1 ternary),
    never the frame
  - the integrated autocorrelation time is measured and an effective sample size
    reported alongside every frame-pooled number
  - frame-pooled quantities (AUC, occupancy fractions) are reported as
    DESCRIPTIVE only, and labelled as such
  - with n=3 vs n=3 no p-value is worth printing; effect size and the
    replicate-level range are what get reported
  - S4 (ternary) has ONE replicate and carries HAB1 + Mn2+, so it is reported
    beside the open/closed contrast, never pooled into it

EQUILIBRATION DISCARD
---------------------
Chosen from the CORE backbone RMSD only -- never from the loop observables being
tested -- so the discard window cannot be tuned to flatter the answer. One window
is chosen for all replicates (the max over replicates, capped so the shortest
replicate keeps usable data), and section 5 shows the conclusion is unchanged
anywhere from 0 to 150 ns of discard.

Run with the pyr1_docking env python.
"""
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LD = os.path.join(ROOT, "data", "loop_dynamics")
MAP = json.load(open(os.path.join(LD, "residue_map.json")))

def natmap(sysname):
    """per-system native->sequential; S4 is shifted by one, see script 58"""
    return {int(k): v for k, v in MAP["systems"][sysname]["native_to_seq"].items()}

# S4 is deliberately its OWN state, never pooled with S2. It is closed too, but it
# also carries HAB1 and Mn2+, so pooling it into "closed" would mix the effect of
# the partner protein into a WT open-vs-closed contrast. The primary comparison
# stays S1 vs S2; S4 is reported alongside.
SYSTEMS = {"S1_apo_open": "open", "S2_holo_closed": "closed",
           "S4_ternary": "ternary"}
PRIMARY = ("open", "closed")
REPS = [0, 1, 2]
PS_PER_FRAME = 10.0     # ntwx=5000 steps * 2 fs
GATE, LATCH, LB7A5 = MAP["gate_native"], MAP["latch_native"], MAP["lb7a5_native"]
CORE_NATIVE = MAP["systems"]["_reference"]["core_native"]


def load_dat(path, col=1):
    """cpptraj .dat -> 1-D float array, skipping the '#' header"""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    v = np.loadtxt(path, comments="#", usecols=(col,), ndmin=1)
    return v.astype(float)


def autocorr_time(x, maxlag=None):
    """
    Integrated autocorrelation time in FRAMES, by the initial-positive-sequence
    estimator: sum the normalised autocorrelation until it first goes <= 0.
    tau=1 means independent frames; N_eff = N / (2*tau).
    """
    x = np.asarray(x, float)
    x = x - x.mean()
    n = len(x)
    if n < 100 or x.std() == 0:
        return np.nan
    maxlag = maxlag or min(n // 4, 20000)
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:maxlag].real
    ac /= ac[0]
    tau = 1.0
    for k in range(1, len(ac)):
        if ac[k] <= 0:
            break
        tau += 2.0 * ac[k]
    return float(tau)


def plateau_frame(rmsd, ngrid=60):
    """
    Equilibration detection by Chodera's rule (J Chem Theory Comput 2016, the
    method behind pymbar.timeseries.detectEquilibration): choose the discard t0
    that MAXIMISES the effective sample size of what remains,
        N_eff(t0) = (N - t0) / (2 * tau(t0)).
    Throwing away data shortens the series but shortens tau faster while the
    system is still relaxing, so the maximum lands at the end of the transient.

    Chosen deliberately over a hand-set threshold. The first version of this
    function used a "within 1 SD of the final-third mean" rule and demanded a
    270 ns discard from a 300 ns run -- an artefact of slow drift in RMSD-to-frame-1,
    not a real equilibration time. A criterion with a free tolerance parameter is
    exactly the kind of knob that can be turned until the answer looks good.

    Searched over the first half only; a system needing more than that has not
    equilibrated and should be reported as such, not quietly trimmed.
    """
    n = len(rmsd)
    best_t0, best_neff = 0, -np.inf
    for t0 in np.linspace(0, n // 2, ngrid, dtype=int):
        seg = rmsd[t0:]
        if len(seg) < 200:
            continue
        tau = autocorr_time(seg)
        if not np.isfinite(tau) or tau <= 0:
            continue
        neff = (n - t0) / (2.0 * tau)
        if neff > best_neff:
            best_neff, best_t0 = neff, int(t0)
    return best_t0


# ------------------------------------------------------------------ load
print("=" * 78)
print("LOADING")
print("=" * 78)
data = {}
for sysname in SYSTEMS:
    for rep in REPS:
        d = os.path.join(LD, sysname, f"rep{rep}")
        go = load_dat(os.path.join(d, "gate_to_open.dat"))
        if go is None:
            print(f"  {sysname} rep{rep}: absent -- skipped")
            continue
        rec = {
            "gate_open": go,
            "gate_closed": load_dat(os.path.join(d, "gate_to_closed.dat")),
            "latch_open": load_dat(os.path.join(d, "latch_to_open.dat")),
            "latch_closed": load_dat(os.path.join(d, "latch_to_closed.dat")),
            "core_first": load_dat(os.path.join(d, "core_rmsd_first.dat")),
            "gl_mindist": load_dat(os.path.join(d, "gate_latch_contacts.dat"), 3),
            "aba": load_dat(os.path.join(d, "aba_rmsd.dat")),
            "rmsf": np.loadtxt(os.path.join(d, "rmsf_bb_byres.dat"),
                               comments="#", ndmin=2),
        }
        n = len(go)
        data[(sysname, rep)] = rec
        print(f"  {sysname} rep{rep}: {n} frames = {n * PS_PER_FRAME / 1000:.0f} ns")

if not data:
    sys.exit("no cpptraj output found -- run 58b first")

# ------------------------------------------------ 1. validity
print()
print("=" * 78)
print("1. VALIDITY")
print("=" * 78)
discards = []
for (s, r), d in sorted(data.items()):
    p = plateau_frame(d["core_first"])
    discards.append(p)
    print(f"  {s} rep{r}: core RMSD plateau at {p * PS_PER_FRAME / 1000:6.1f} ns"
          f"   (final level {d['core_first'][int(2*len(d['core_first'])/3):].mean():.2f} A)")
DISCARD = int(np.ceil(max(discards) * PS_PER_FRAME / 1000 / 10.0) * 10)
# One window for every replicate so the systems stay comparable. Cap it: the
# shortest replicate must keep a usable amount of data, and a discard that eats
# most of the run means the answer is being manufactured by the trimming.
SHORTEST_NS = min(len(d["gate_open"]) for d in data.values()) * PS_PER_FRAME / 1000
CAP = int(min(100, SHORTEST_NS / 3))
if DISCARD > CAP:
    print(f"\n  !! detector asked for {DISCARD} ns; capping at {CAP} ns "
          f"(shortest replicate is {SHORTEST_NS:.0f} ns)")
    DISCARD = CAP
DF = int(DISCARD * 1000 / PS_PER_FRAME)
print(f"\n  -> discarding the first {DISCARD} ns of every replicate")
print(f"     (chosen from CORE backbone RMSD only, never from the loop "
      f"observables under test)")

print("\n  ABA validity control (S2 only; heavy-atom RMSD vs frame 1, protein-fit):")
aba_ok = True
for (s, r), d in sorted(data.items()):
    if d["aba"] is None:
        continue
    a = d["aba"][DF:]
    if len(a) == 0:
        print(f"    {s} rep{r}: no frames left after discard -- skipped")
        continue
    print(f"    {s} rep{r}: mean {a.mean():5.2f} A, max {a.max():5.2f} A, "
          f"final {a[-1]:5.2f} A")
    if a.mean() > 5.0:
        aba_ok = False
print(f"    -> ABA {'stays in the pocket' if aba_ok else 'HAS MOVED -- S2 is suspect'}")

# ------------------------------------------------ 2. state stability
print()
print("=" * 78)
print("2. STATE STABILITY  (state coordinate S = d_to_open - d_to_closed)")
print("=" * 78)
print("  S < 0 means closer to the OPEN reference; S > 0 means closer to CLOSED.")
print("  A system that holds its state keeps S on one side for the whole run.\n")

rows = []
for (s, r), d in sorted(data.items()):
    for loop, (o, c) in (("gate", ("gate_open", "gate_closed")),
                         ("latch", ("latch_open", "latch_closed"))):
        do, dc = d[o][DF:], d[c][DF:]
        S = do - dc
        tau = autocorr_time(S)
        neff = len(S) / (2 * tau) if tau and tau > 0 else np.nan
        cross = int(np.sum(np.diff(np.sign(S)) != 0))
        rows.append(dict(sys=s, rep=r, loop=loop,
                         d_open=do.mean(), d_closed=dc.mean(),
                         S=S.mean(), S_sd=S.std(), frac_open=float((S < 0).mean()),
                         tau_ns=tau * PS_PER_FRAME / 1000 if tau else np.nan,
                         neff=neff, crossings=cross))

hdr = (f"  {'system':<16}{'rep':>3} {'loop':<6}{'d>open':>8}{'d>closed':>9}"
       f"{'S':>8}{'sd':>6}{'%open':>7}{'tau/ns':>8}{'Neff':>7}")
print(hdr)
print("  " + "-" * (len(hdr) - 2))
for x in rows:
    print(f"  {x['sys']:<16}{x['rep']:>3} {x['loop']:<6}{x['d_open']:8.2f}"
          f"{x['d_closed']:9.2f}{x['S']:8.2f}{x['S_sd']:6.2f}"
          f"{100*x['frac_open']:6.0f}%{x['tau_ns']:8.1f}{x['neff']:7.0f}")
print("\n  'sd' is the within-replicate spread of S. Compare it with the "
      "open-vs-closed\n  separation in section 3: separation >> sd is what makes "
      "a filter possible.")
print("  tau is the integrated autocorrelation time and Neff = N/(2*tau) the")
print("  effective number of INDEPENDENT samples in that replicate -- single")
print("  digits in several cases, which is why frames are not the unit of "
      "inference.")

# ------------------------------------------------ 3. separability
print()
print("=" * 78)
print("3. SEPARABILITY  (replicate-level; n is replicates, NOT frames)")
print("=" * 78)


def auc(a, b):
    """P(random draw from b > random draw from a); 0.5 = no separation."""
    a, b = np.asarray(a), np.asarray(b)
    allv = np.concatenate([a, b])
    order = allv.argsort()
    ranks = np.empty(len(allv), float)
    ranks[order] = np.arange(1, len(allv) + 1)
    # average ties
    _, inv, cnt = np.unique(allv, return_inverse=True, return_counts=True)
    sums = np.zeros(len(cnt))
    np.add.at(sums, inv, ranks)
    ranks = (sums / cnt)[inv]
    ra = ranks[:len(a)].sum()
    return float((ra - len(a) * (len(a) + 1) / 2) / (len(a) * len(b)))


OBS = {
    "gate S (d_open-d_closed)": lambda d: d["gate_open"][DF:] - d["gate_closed"][DF:],
    "latch S":                  lambda d: d["latch_open"][DF:] - d["latch_closed"][DF:],
    "gate-latch min dist":      lambda d: d["gl_mindist"][DF:],
    "gate RMSD to closed":      lambda d: d["gate_closed"][DF:],
    "latch RMSD to closed":     lambda d: d["latch_closed"][DF:],
}

summary = {}
print(f"  {'observable':<26}{'open reps (mean/rep)':<30}{'closed (S2)':<20}"
      f"{'S4':<8}{'gap':>7}{'AUC':>6}")
print("  " + "-" * 97)
for name, fn in OBS.items():
    per = {st: [] for st in SYSTEMS.values()}
    pool = {st: [] for st in SYSTEMS.values()}
    for (s, r), d in sorted(data.items()):
        v = fn(d)
        per[SYSTEMS[s]].append(v.mean())
        pool[SYSTEMS[s]].append(v)
    if not per["open"] or not per["closed"]:
        continue
    a = np.concatenate(pool["open"])
    b = np.concatenate(pool["closed"])
    A = auc(a, b)
    lo_o, hi_o = min(per["open"]), max(per["open"])
    lo_c, hi_c = min(per["closed"]), max(per["closed"])
    # gap between the two replicate-level ranges; negative means they overlap
    gap = (lo_c - hi_o) if lo_c > hi_o else (lo_o - hi_c if lo_o > hi_c else
                                            -(min(hi_o, hi_c) - max(lo_o, lo_c)))
    summary[name] = dict(open_reps=per["open"], closed_reps=per["closed"],
                         ternary_reps=per["ternary"], auc=A, gap=float(gap))
    o_txt = " ".join(f"{v:.2f}" for v in per["open"])
    c_txt = " ".join(f"{v:.2f}" for v in per["closed"])
    t_txt = " ".join(f"{v:.2f}" for v in per["ternary"]) or "-"
    print(f"  {name:<26}{o_txt:<30}{c_txt:<20}{t_txt:<8}{gap:7.2f}{A:6.2f}")

print("\n  'gap' = distance between the open and closed REPLICATE-MEAN ranges.")
print("  Positive = the replicate means do not overlap at all (a usable filter).")
print("  Negative = they overlap by that much (not usable on its own).")
print("  AUC is frame-pooled and DESCRIPTIVE only -- frames are autocorrelated,")
print("  so it flatters separation; trust the replicate-level gap instead.")

# ------------------------------------------------ RMSF
print()
print("=" * 78)
print("4. RMSF by loop (backbone, fit to average structure)")
print("=" * 78)


def loop_rmsf(rec, nats, sysname):
    seq = [natmap(sysname)[n] for n in nats]
    arr = rec["rmsf"]
    idx = {int(row[0]): row[1] for row in arr}
    return float(np.mean([idx[s] for s in seq if s in idx]))


print(f"  {'system':<16}{'rep':>3}{'gate':>8}{'latch':>8}{'Lb7a5':>8}{'core':>8}")
rmsf_rows = []
for (s, r), d in sorted(data.items()):
    g, l, lb = (loop_rmsf(d, GATE, s), loop_rmsf(d, LATCH, s),
                loop_rmsf(d, LB7A5, s))
    core = loop_rmsf(d, CORE_NATIVE, s)
    rmsf_rows.append(dict(sys=s, rep=r, gate=g, latch=l, lb7a5=lb, core=core))
    print(f"  {s:<16}{r:>3}{g:8.2f}{l:8.2f}{lb:8.2f}{core:8.2f}")

for loop in ("gate", "latch", "lb7a5"):
    o = [x[loop] for x in rmsf_rows if SYSTEMS[x["sys"]] == "open"]
    c = [x[loop] for x in rmsf_rows if SYSTEMS[x["sys"]] == "closed"]
    if o and c:
        print(f"    {loop:<6} open {np.mean(o):.2f} (range {min(o):.2f}-{max(o):.2f})"
              f"   closed {np.mean(c):.2f} (range {min(c):.2f}-{max(c):.2f})"
              f"   ratio {np.mean(o)/np.mean(c):.2f}x")

# ------------------------------------------------ 5. discard sensitivity
print()
print("=" * 78)
print("5. SENSITIVITY TO THE EQUILIBRATION DISCARD")
print("=" * 78)
print("  The detector wanted different windows for different replicates, and one")
print("  common window had to be capped. If the separation is real it should not")
print("  care about that choice. Gate S, replicate-mean ranges:\n")
print(f"  {'discard':>8}  {'open reps':<24}{'closed reps':<24}{'gap':>7}")
for dns in (0, 25, 50, 100, 150):
    o, c = [], []
    for (s, rr), d in sorted(data.items()):
        f0 = int(dns * 1000 / PS_PER_FRAME)
        seg_o, seg_c = d["gate_open"][f0:], d["gate_closed"][f0:]
        if len(seg_o) < 1000:
            continue
        st = SYSTEMS[s]
        if st not in PRIMARY:      # keep the ternary out of the open/closed contrast
            continue
        (o if st == "open" else c).append(float((seg_o - seg_c).mean()))
    if not o or not c:
        print(f"  {dns:>6} ns  (a replicate is too short at this discard)")
        continue
    gap = min(c) - max(o)
    print(f"  {dns:>6} ns  {' '.join(f'{v:.2f}' for v in o):<24}"
          f"{' '.join(f'{v:.2f}' for v in c):<24}{gap:7.2f}")
print("\n  A gap that stays large and positive across every window means the")
print("  open/closed split is a property of the trajectories, not of the trimming.")

out = dict(discard_ns=DISCARD, ps_per_frame=PS_PER_FRAME,
           per_replicate=rows, separability=summary, rmsf=rmsf_rows,
           aba_ok=bool(aba_ok),
           n_open=sum(1 for k in data if SYSTEMS[k[0]] == "open"),
           n_closed=sum(1 for k in data if SYSTEMS[k[0]] == "closed"))
json.dump(out, open(os.path.join(LD, "summary.json"), "w"), indent=1, default=float)
print(f"\nwrote {LD}/summary.json")
