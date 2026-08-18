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

# EVERY UNIT IS ITS OWN STATE. Nothing is pooled that differs in more than one
# way, because the whole point of README 24e was that S1 and S2 differ in BOTH
# conformation and ligand occupancy.
#
#   open           S1, 3K3K chain A monomer, apo + open
#   closed         S2, 3QN1 chain A monomer, ABA + closed
#   dimer_openA    S3 chain A, apo + open, but inside a dimer
#   dimer_closedB  S3 chain B, apo + CLOSED  <-- the missing cell of 19b
#   ternary        S4, ABA + closed + HAB1 + Mn2+
#
# S4 is closed too, but pooling it into "closed" would mix the partner protein
# into a WT open-vs-closed contrast. S3's two protomers share a box and a
# thermostat, so they are not independent of each other -- they are reported as
# their own states and never averaged together or into S1/S2.
#
# THE CELL THAT MATTERS. dimer_closedB is a closed, ligand-shaped protomer
# simulated with an EMPTY pocket: 3K3K chain B was ABA-bound in the crystal
# (README 28g) and S3 was built protein-only. It was not designed as the
# apo-closed control -- it is one by accident -- so its confounds are stated
# rather than glossed: it sits in a dimer against an open protomer, and it is on
# the OLD build protocol. It is still the only apo-closed trajectory that exists.
SYSTEMS = {"S1_apo_open": "open", "S2_holo_closed": "closed",
           "S3_dimer_openA": "dimer_openA",
           "S3_dimer_closedB": "dimer_closedB",
           "S4_ternary": "ternary"}
# The headline open-vs-closed contrast stays S1 vs S2 -- the same pair README 24
# reported -- so adding systems cannot quietly redefine the published result.
PRIMARY = ("open", "closed")
STATE_ORDER = ["open", "closed", "dimer_openA", "dimer_closedB", "ternary"]
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

# ---- do the NAMES match the MEASUREMENTS? ----
# S3's two protomers are told apart only by a chain letter, and the label
# "S3_dimer_closedB" is an assertion about which one is closed. A label that
# encodes a claim will eventually be trusted instead of the claim, so the claim is
# checked here against the trajectory: every unit whose name says open must sit on
# the open side of the watershed, and vice versa. A swapped chain letter anywhere
# upstream stops the script instead of producing a mislabelled figure.
print("\n  Name-vs-measurement check (gate state coordinate S, post-discard):")
EXPECT_SIDE = {"open": -1, "dimer_openA": -1,
               "closed": +1, "dimer_closedB": +1, "ternary": +1}
for (sname, r), d in sorted(data.items()):
    S = (d["gate_open"][DF:] - d["gate_closed"][DF:]).mean()
    want = EXPECT_SIDE[SYSTEMS[sname]]
    side = "open-like" if S < 0 else "closed-like"
    assert np.sign(S) == want, (
        f"{sname} rep{r}: name says {'open' if want < 0 else 'closed'} but mean "
        f"gate S = {S:+.2f} ({side}). The chain letter or the state table is wrong.")
    print(f"    {sname:<18} rep{r}  S = {S:+5.2f}  {side}  -- agrees with its name")
# and the two S3 protomers must be on OPPOSITE sides, which is the whole reason
# S3 is analysed as two units rather than one
sA = np.mean([(data[("S3_dimer_openA", r)]["gate_open"][DF:]
               - data[("S3_dimer_openA", r)]["gate_closed"][DF:]).mean()
              for r in REPS if ("S3_dimer_openA", r) in data])
sB = np.mean([(data[("S3_dimer_closedB", r)]["gate_open"][DF:]
               - data[("S3_dimer_closedB", r)]["gate_closed"][DF:]).mean()
              for r in REPS if ("S3_dimer_closedB", r) in data])
assert sA < 0 < sB, "the S3 protomers are not on opposite sides of the watershed"
print(f"    -> S3 is a MIXED dimer in the trajectory too, not just the crystal: "
      f"protomer A {sA:+.2f}, protomer B {sB:+.2f}")

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
print("  Each cell is the REPLICATE MEANS for that state, one number per replicate.")
print("  'gap' and 'AUC' are the S1-vs-S2 contrast only (the README 24 pair).\n")
hdr3 = f"  {'observable':<26}" + "".join(f"{st:<22}" for st in STATE_ORDER) + f"{'gap':>7}{'AUC':>6}"
print(hdr3)
print("  " + "-" * (len(hdr3) - 2))
for name, fn in OBS.items():
    per = {st: [] for st in STATE_ORDER}
    pool = {st: [] for st in STATE_ORDER}
    for s, r in sorted(data):
        v = fn(data[(s, r)])
        per[SYSTEMS[s]].append(v.mean())
        pool[SYSTEMS[s]].append(v)
    if not per["open"] or not per["closed"]:
        continue
    A = auc(np.concatenate(pool["open"]), np.concatenate(pool["closed"]))
    lo_o, hi_o = min(per["open"]), max(per["open"])
    lo_c, hi_c = min(per["closed"]), max(per["closed"])
    # gap between the two replicate-level ranges; negative means they overlap
    gap = (lo_c - hi_o) if lo_c > hi_o else (lo_o - hi_c if lo_o > hi_c else
                                            -(min(hi_o, hi_c) - max(lo_o, lo_c)))
    summary[name] = dict({st + "_reps": per[st] for st in STATE_ORDER},
                         auc=A, gap=float(gap))
    cells = "".join((" ".join(f"{v:.2f}" for v in per[st]) or "-").ljust(22)
                    for st in STATE_ORDER)
    print(f"  {name:<26}{cells}{gap:7.2f}{A:6.2f}")

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

print()
for loop in ("gate", "latch", "lb7a5"):
    line = f"    {loop:<6}"
    for st in STATE_ORDER:
        v = [x[loop] for x in rmsf_rows if SYSTEMS[x["sys"]] == st]
        line += f"  {st} {np.mean(v):.2f} ({min(v):.2f}-{max(v):.2f})" if v else ""
    print(line)
o = [x["gate"] for x in rmsf_rows if SYSTEMS[x["sys"]] == "open"]
c = [x["gate"] for x in rmsf_rows if SYSTEMS[x["sys"]] == "closed"]
if o and c:
    print(f"\n    gate open/closed ratio {np.mean(o)/np.mean(c):.2f}x")

# ------------------------------------------------ 5. the apo-closed question
print()
print("=" * 78)
print("5. DOES A CLOSED PROTOMER OPEN WHEN THE LIGAND IS ABSENT?")
print("=" * 78)
print("""  README 24 could not ask this: every closed replicate there had ABA in it, so
  conformation and occupancy were confounded (24e). S3 chain B answers it by
  accident -- 3K3K chain B is closed and was ABA-bound in the crystal, and S3 was
  built protein-only, so it is a closed, ligand-shaped protomer around an empty
  pocket.

  Read DRIFT, not level. The question is not "is S positive" (it starts positive
  by construction) but "does it move toward the open basin over 300 ns".
  first/last are the mean of the first and last 50 ns after the discard.
""")
WIN = int(50 * 1000 / PS_PER_FRAME)
drift_rows = []
for loop, (ok, ck) in (("gate", ("gate_open", "gate_closed")),
                       ("latch", ("latch_open", "latch_closed"))):
    print(f"  {loop.upper()}")
    print(f"    {'unit':<18}{'rep':>3}{'S first':>9}{'S last':>9}{'drift':>8}"
          f"{'min':>8}{'max':>8}{'cross':>7}")
    for (sname, r), d in sorted(data.items()):
        S = d[ok][DF:] - d[ck][DF:]
        if len(S) < 2 * WIN:
            continue
        a, b = float(S[:WIN].mean()), float(S[-WIN:].mean())
        cross = int(np.sum(np.diff(np.sign(S)) != 0))
        drift_rows.append(dict(unit=sname, rep=r, loop=loop, S_first=a, S_last=b,
                               drift=b - a, S_min=float(S.min()),
                               S_max=float(S.max()), crossings=cross))
        print(f"    {sname:<18}{r:>3}{a:9.2f}{b:9.2f}{b-a:8.2f}"
              f"{S.min():8.2f}{S.max():8.2f}{cross:7d}")
    print()

# The comparison that decides it: the apo-closed protomer against the two closed
# systems that DO have a ligand. If removing ABA released the gate, apo-closed
# should drift negative while holo-closed and the ternary do not.
gd = [x for x in drift_rows if x["loop"] == "gate"]
by_state = {}
for x in gd:
    by_state.setdefault(SYSTEMS[x["unit"]], []).append(x["drift"])
print("  Gate drift over 300 ns, by state (negative = moving toward OPEN):")
for st in STATE_ORDER:
    if st not in by_state:
        continue
    v = by_state[st]
    print(f"    {st:<16} {np.mean(v):+6.2f} A   (per replicate: "
          f"{' '.join(f'{x:+.2f}' for x in v)})")
n_cross = sum(x["crossings"] for x in gd
              if SYSTEMS[x["unit"]] == "dimer_closedB")
print(f"\n  apo-closed protomer: {n_cross} crossings of the open/closed watershed "
      f"in {len([x for x in gd if SYSTEMS[x['unit']] == 'dimer_closedB'])} x 300 ns")
print("""
  A null here is WEAK EVIDENCE, and must be reported as such. Gate opening after
  ligand loss is a barrier crossing, and README 24 already showed no crossing is
  sampled in 1.8 us aggregate from either basin. "It stayed closed" is therefore
  consistent both with the gate being ligand-independent AND with the run simply
  being too short. A DRIFT toward open is the informative outcome; its absence is
  not the converse.""")

# ------------------------------------------------ 5. discard sensitivity
print()
print("=" * 78)
print("6. SENSITIVITY TO THE EQUILIBRATION DISCARD")
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
           drift=drift_rows, aba_ok=bool(aba_ok),
           core_residues=len(CORE_NATIVE),
           n_replicates={st: sum(1 for k in data if SYSTEMS[k[0]] == st)
                         for st in STATE_ORDER},
           n_open=sum(1 for k in data if SYSTEMS[k[0]] == "open"),
           n_closed=sum(1 for k in data if SYSTEMS[k[0]] == "closed"))
json.dump(out, open(os.path.join(LD, "summary.json"), "w"), indent=1, default=float)
print(f"\nwrote {LD}/summary.json")
