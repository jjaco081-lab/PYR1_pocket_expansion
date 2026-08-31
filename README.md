# PYR1 Binding-Pocket Expansion

Enlarging the ligand-binding cavity of *Arabidopsis thaliana* PYR1 so the
PYR1–HAB1 chemically-induced-dimerisation biosensor can accept larger and more
chemically diverse ligands, **without** altering the gate and latch loops that
drive the conformational switch.

Project started 2026-08-06. Working directory
`/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion`.

---

## 0. AI-assisted work: declaration of contributions

**Status of this work.** Everything in this repository is a **preliminary
computational viability assessment** — a pilot whose purpose is to decide whether
enlarging the PYR1 pocket is worth pursuing experimentally, and to identify which
designs would be worth building. **No wet-lab validation has been performed.**
Every quantity reported here is computational, and the substantive claims are
hypotheses supported by modelling, not experimental results.

**Model and dates.** All AI-assisted work to date was performed with
**Claude Opus 5**, **2026-08-06 through 2026-08-11**. If a different model
contributes later, it is recorded as a separate dated block here rather than
folded into this one.

### Division of contributions

**Jannis Jacobs (project lead) — direction, judgement, and domain expertise**

- Conceived the project and set its objective, and **re-specified the success
  criterion repeatedly** as the work exposed ambiguity in it: the deliverable is
  new *library positions* for a Tian-style oligo-pool screen rather than a larger
  cavity as such; a change at an already-sampled position still counts if the
  substitution or the combination was untested; and the goal is a **portfolio** of
  scaffolds rather than a single ranked winner. These are the definitions the
  entire analysis is built around.
- Supplied the biological and experimental expertise the modelling cannot supply:
  PYR1–HAB1 biology, Y2H screening with FOA counter-selection, Tian library design
  and protocol, which supplementary datasets carry focused second-round libraries,
  and the colony-picking practice that makes non-recovery of a sequence
  uninformative — a caveat that invalidates an entire class of benchmark design.
- Selected and supplied the literature (Tian 2025; Leonard 2026; Dorosh 2013;
  Melcher 2010; Peterson 2010; the BoltzMol-1 preprint) and read it independently.
- **Corrected AI errors on points that changed conclusions**, including: that HAB1
  contacts ABA through a bridging water (overturning a stated claim that it never
  touches the ligand); that the mandipropamid complexes 4WVO/8EY0 likely do *not*
  use that water (verified — ABA sensing through it is real but **not conserved**,
  which overturned a proposed design filter — §23a);
  that the earlier LigandMPNN benchmark is confounded by pose and should not be
  treated as evidence (§23b); and that a maximally open pocket remains a valid
  shrink-to-fit starting point, reversing an AI over-correction.
- Proposed analyses that became part of the work: the floppy-ligand burial survey
  (§23d), using Tian actives that modestly enlarge lobe 1 as a small-amplitude
  calibration, and PYR1\*/HAB1\* as a validated negative control pair.
- Made all operational decisions: partition and scheduling policy, which arms to
  run, when to wait for results rather than act on partial data, and the
  preemption strategy for the MD campaign.

**Claude Opus 5 — implementation, computation, and drafting**

- Wrote every script in `scripts/` (45 files) and the analysis they implement:
  the cavity/volume methods, the three Rosetta arms, the MD preparation and run
  system, the novelty and variance-decomposition analyses, the Foldseek and
  homolog surveys, and the benchmark tooling.
- Executed and managed all computation, including SLURM job construction,
  scheduling, monitoring, and recovery.
- Performed the structural measurements reported here, including the W385 water
  analysis across 3QN1 / 4WVO / 8EY0 (§23a).
- Read the supplied literature and extracted the specific transferable methods and
  the points where they do not transfer.
- Wrote and maintains this README, the §22 change log, and the memory of
  superseded reasoning.
- Performed the statistical analysis, including identifying several of its own
  errors — the incorrect discrimination statistic (§14a), the sequence
  position-mapping bug (§23e), and two latent MD restart bugs (§23c).

### Why this is stated in this much detail

The scientific direction, the success criteria, the domain constraints, and the
corrections that changed conclusions are the project lead's. The implementation,
computation, and drafting are AI-generated and should be read as such: fast and
broad, but with a documented error rate.

That error rate is the reason for §22 and for version control. Several AI-produced
errors were caught during this work — some by the project lead, some by subsequent
analysis — and each is recorded in place rather than silently corrected. **There is
no guarantee that all have been found.** Anyone building on this should treat
§22 as a live list, re-derive load-bearing numbers before relying on them, and
weight the experimental claims of §12–§16 accordingly.

---

## 1. Motivation and the design constraint

The PYR1–HAB1 system is unusually tolerant of ligand reprogramming — a handful
of pocket substitutions retargets it to new small molecules — but it has a
**size ceiling**. Ligands beyond roughly the ABA volume class are not
accommodated, which excludes a large part of chemical space from this sensor
architecture.

The constraint that makes this hard is that **PYR1–HAB1 is a ratchet, not a
molecular glue**. Contact analysis of PDB 3QN1 (`scripts/01_pocket_contacts.py`,
output `logs/01_pocket_contacts.txt`) shows 21 PYR1 residues within 5.0 Å of
ABA and exactly **one** HAB1 residue — Trp385, at 4.67 Å. HAB1 contributes
essentially nothing to the binding site; it reads the **closed-state surface**
of PYR1 and wedges W385 into the gate–latch groove. Ligand chemistry is
invisible to HAB1 except through PYR1's conformation.

Two consequences define the whole project:

1. **The pocket and the readout are structurally decoupled but conformationally
   coupled.** Cavity geometry is a free variable; the closed-state epitope is
   not.
2. **The gate (β3–β4, S85–A89) and latch (β5–β6, H115–L117) are off-limits.**
   They must still close, and must close *only* when ligand is bound. A variant
   that binds a large ligand but no longer switches — or that switches
   constitutively — is a failed sensor even with a perfect pocket.

This is why the RayGun generative approach was abandoned: it preferentially
inserted residues into the gate and latch loops, precisely the regions that
must be conserved.

### Why not simply insert residues?

Considered and deferred, for reasons worth recording:

- **Mid-strand single insertions are register-breaking.** A +1 insertion in a
  β-strand flips side-chain in/out alternation and shears the backbone H-bond
  register. The minimum register-preserving unit on a strand is **+2**, and
  because the PYR1 sheet is curved, paired strands likely need to move
  together. Single-position insertion scanning would therefore mostly report
  false negatives.
- **Loop grafting from large-cavity homologs** (START domains such as STARD1 /
  CERT / MLN64 bind cholesterol and ceramide in the same helix-grip fold at
  >1000 Å³) remains attractive precisely because donor geometry is native and
  folds more reliably than *de novo* insertions. Script `08` is the first step
  toward that: identify whether large-cavity homologs share a signature.
- **RFdiffusion motif-scaffolding was rejected.** The gate and latch are not a
  motif, they are a two-state switch; diffusion models sample one state. A
  scaffold displaying closed-state geometry carries no guarantee it can open,
  and presents a novel surface HAB1 has never seen. (For future reference, the
  ligand-aware entry points are **RFdiffusionAA** (Krishna *et al.*, *Science*
  2024) and **RFdiffusion2**, both installed under `HAB1_shrinkinator/`.)

The pilot documented here therefore attacks the **wall** rather than the
backbone length: which side chains occlude the cavity, whether removing them
actually yields a larger pocket after repacking, and whether the switch
survives.

---

## 2. Key structural finding — the pocket wall is a buttressed cluster

Naive inspection identifies **K59** and **F108** as the restricting residues,
and the volume scan agrees. But the contact network
(`scripts/01_pocket_contacts.py`) shows they are held by a second-shell
scaffold:

| Interaction | Distance |
|---|---|
| E94 OE1 – R79 NE | 2.76 Å |
| E94 OE1 – R79 NH2 | 2.93 Å |
| F108 ring centroid – R79 guanidinium (cation–π) | 3.55 Å |
| K59 NZ – ABA carboxylate O12 | 2.85 Å |
| K59 NZ – F108 ring centroid | 4.39 Å |

**F108's rotamer is pinned by R79, which is locked by the E94 salt bridge.**

Minimum **heavy-atom side-chain** distances to ABA (not Cα — Cα distances run
6.6–11.7 Å and are not informative here):

| Residue | closest atom pair | min heavy-atom dist | Cα–ligand dist |
|---|---|---|---|
| K59 | NZ···ABA O12 | 2.85 Å | 6.62 Å |
| F108 | CZ···ABA O12 | 3.48 Å | 8.38 Å |
| Y120 | OH···ABA C2 | 3.60 Å | 7.73 Å |
| E141 | OE2···ABA O12 | 3.99 Å | 8.83 Å |
| **E94** | OE2···ABA C15 | **4.33 Å** | 8.40 Å |
| **R79** | NH2···ABA O11 | **5.81 Å** | 11.71 Å |

**⚠️ Distance-to-ligand is the WRONG metric for "lines the pocket."** ABA does
not fill the PYR1 cavity, so a residue can face the cavity surface directly and
still be >5 Å from the ligand — it simply lines a region the ligand does not
occupy. Ranking by ligand distance mislabels exactly those residues, which are
the interesting ones for pocket expansion.

`scripts/04_cavity_lining.py` measures lining directly: cavity grid points
within 4.5 Å of each residue's heavy atoms, plus a line-of-sight test.

| Residue | lining pts (WT) | lining pts (expanded) | dist to cavity surface (WT) | dist to ABA |
|---|---|---|---|---|
| K59 | 545 | 3712 | 2.96 Å | 2.85 Å |
| Y120 | 628 | 986 | 2.92 Å | 3.60 Å |
| E94 | 408 | 2876 | 2.93 Å | 4.33 Å |
| F108 | 244 | 5987 | 3.06 Å | 3.48 Å |
| E141 | 57 | 484 | 2.93 Å | 3.99 Å |
| **R79** | **4** | **3374** | **4.16 Å** | 5.81 Å |

Line-of-sight, marching the straight segment from each atom to the nearest ABA
atom and testing for intervening van der Waals spheres:

| Segment | Length | Result |
|---|---|---|
| R79 **NH2** → ABA | 5.81 Å | **CLEAR — no intervening PYR1 atom** |
| R79 **NE** → ABA | 7.28 Å | **CLEAR** |
| R79 NH1 → ABA | 7.64 Å | blocked by K59 NZ |
| E94 OE2 → ABA | 4.33 Å | CLEAR |

So the correct description of the cluster is:

- **E94 is first shell** (4.33 Å): it lines the pocket directly and occludes
  volume in its own right, *in addition to* anchoring R79.
- **R79 is a cavity-RIM residue on the far side of the pocket from ABA** — not
  a remote second-shell residue. Its guanidinium touches the edge of the WT
  cavity (4 lining points at 4.16 Å) with **open, unobstructed space between
  NH2 and the ligand**, and it becomes one of the dominant lining residues of
  the expanded cavity (3374 points, 0.27 Å). It simultaneously stacks on F108
  (3.55 Å cation–π) and is locked by E94.

  *(Note: the WT lining count of 4 is an artefact of this project's deliberately
  conservative cavity criterion — 0.88 buriedness cut, single connected
  component. The ChimeraX cavity-finder surface, which is more permissive,
  shows R79 clearly contacting the pocket. The line-of-sight test above is
  criterion-independent and is the more reliable statement.)*

This reframes the R79A/E94A result (+84.8 Å³, §3) and the project as a whole:
R79 and E94 are not merely holding a wall from behind — **they line a
sub-pocket that ABA does not occupy**, on the far side of the cavity. The
expanded-cavity column shows the consequence: truncating the quad recruits
**I48, V49, R50, V81, I62, S122** as new lining residues, i.e. it opens a
contiguous region that is sealed off in WT. That is precisely the "satellite
lobe separated by a wall" this project set out to find, and is consistent with
the hypothesis that the E94–R79 salt bridge partitions a larger ancestral
pocket.

---

## 3. Rigid-backbone truncation scan (upper bound)

`scripts/02_rigid_truncation_scan.py` → `results/02_rigid_truncation.csv`,
`logs/02_rigid_truncation.txt`

Each first-shell side chain is deleted beyond CB *in silico* with nothing else
moved, and the enclosed cavity re-measured. The delta is an **upper bound** on
the volume that residue withholds.

**WT closed-state cavity: 166.4 Å³**

| Variant | Cavity (Å³) | Δ (Å³) |
|---|---|---|
| F108→A | 239.9 | **+73.5** |
| K59→A | 238.1 | **+71.8** |
| E141→A | 235.8 | +69.4 |
| Y120→A | 212.4 | +46.0 |
| L117→A | 198.4 | +32.0 |
| E94→A | 184.8 | +18.4 |
| I110→A | 183.0 | +16.6 |
| **F61→A** | 149.9 | **−16.5** |
| **L87→A** | 150.2 | **−16.1** |
| **F159→A** | 151.9 | **−14.5** |

Combinations:

| Combination | Cavity (Å³) | Δ (Å³) |
|---|---|---|
| K59A / F108A | 294.2 | +127.9 |
| **R79A / E94A** | 251.1 | **+84.8** |
| F108A / R79A | 265.1 | +98.8 |
| F108A / R79A / E94A | 308.5 | +142.1 |
| **K59A / F108A / R79A / E94A** | **375.5** | **+209.1** |

Two results matter:

- **R79A/E94A alone unlocks +84.8 Å³ — more than either K59A or F108A —
  despite neither residue contacting the ligand.** This is direct evidence for
  the buttress model: the second shell, not the first, is holding the wall.
- The full quad reaches **375.5 Å³, a 2.3-fold increase** over WT.

**Negative deltas are informative, not noise.** F61, L87 and F159 *seal* the
cavity from bulk solvent; deleting them causes the connected component to leak
and buriedness to fall below cut-off. They must be left intact — removing them
opens the pocket to solvent rather than enlarging it.

⚠️ These are rigid-backbone upper bounds with no repacking. The post-relax
numbers (§5) are the real result; the gap between them is side-chain infilling.

---

## 4. Structural homology across the helix-grip/SRPBCC superfamily

`scripts/05_foldseek_search.sh` → `results/foldseek/fs_hits.tsv` (2053 hits),
parsed by `scripts/03_parse_foldseek.py`.

PYR1 chain A searched against foldseek **CATH50**; 266 hits at
alnTM ≥ 0.5 constitute the helix-grip core set (CATH superfamily 3.30.530.20).
Residue identity at each wall position, mapped through the foldseek alignment:

| PYR1 position | WT | Superfamily distribution (n ≈ 230–246 ungapped) |
|---|---|---|
| 59 | K | **D 15%**, K 10%, L 9%, S 8%, P 7%, A 7%, N 6%, V 6% |
| 79 | R | **R 19%**, F 10%, V 10%, I 9%, L 8%, A 8%, Y 7% |
| 94 | E | **E 23%**, S 11%, G 7%, Q 7%, A 7%, L 6% |
| 108 | F | **Y 24%**, F 14%, I 9%, V 9%, W 7%, L 6% |
| 120 | Y | **Y 19%**, L 16%, V 13%, F 10%, S 8% |
| 141 | E | **V 15%**, E 11%, L 11%, W 11%, Y 8% |

Interpretation:

- **K59 is not conserved.** Lys is only 10%; the superfamily consensus is Asp,
  and hydrophobics (L/V/I/A ≈ 28%) are common. K59 is a **PYR/PYL-specific
  adaptation for the ABA carboxylate** — and if the goal is larger, chemically
  different ligands, that anchor is not needed. This is the safest large-gain
  position in the panel.
- **F108 runs the other way.** Aromatic character dominates (Y 24% + F 14% +
  W 7% ≈ 45%), and **Tyr is more common than Phe — the superfamily consensus is
  *larger* than what PYR1 has.** Combined with the R79 cation–π stacking, F108
  is load-bearing. Shrinking it goes against the grain, which is why F108Y is
  included as an against-the-grain control.
- **R79, E94, E141 are variable across the superfamily but invariant across
  close PYR/PYL homologs** (alnTM > 0.89 all show K59/R79/E94/Y120/E141). That
  pattern — variable in the fold, fixed in the family — is the signature of
  **family-specific functional constraint rather than fold constraint**. The
  fold should tolerate losing them; the *switch* might not. Testing exactly
  that is the purpose of Arm 2.
- A distinct distant subfamily (alnTM 0.79–0.87) carries a **hydrophobic-59 /
  Y108 / W141** signature, which motivates the `K59L_F108Y_E141W` variant.

---

## 4b. Prior art: what Tian *et al.* already sampled

`scripts/11_tian_library_overlap.py` and `scripts/12_tian_combinations.py`
→ `logs/11_tian_overlap.txt`, `logs/12_tian_combinations.txt`

Source: supplementary datasets of **Tian *et al.*, PNAS,
doi 10.1073/pnas.2519924122**, stored at
`/bigdata/cutlerlab/jjaco081/mutation_prediction_benchmark/pnas.2519924122.sd*.xlsx`.
Dataset S3 (`sd04`) is the **library design** — per-position allowed amino
acids. Datasets S6/S7/S8 (`sd07`/`sd08`/`sd09`) are **419 characterised
clones** with per-position recovered residues.

Six libraries (TSM, DSM-Hao, Coumarin, PFAS, TNTv1, TNTv2) randomise **18
pocket positions**: 59, 81, 83, 87, 89, 92, 94, 108, 110, 117, 120, 122, 141,
159, 160, 163, 164, 167.

### Overlap with our six wall positions

| Position | In their libraries? | Recovered in hits | Volume-opening subs recovered |
|---|---|---|---|
| K59 | yes — all 6 libraries | 262 clones (N, M, Q, D, A, L, R, T) | A D N Q T |
| **R79** | **NEVER randomised** | **0 clones** | — |
| E94 | yes — 4 libraries | 189 clones (G, D, A, S, T) | A D G N S T |
| F108 | yes — 3 libraries | 42 clones | H L Q (3 clones only) |
| Y120 | yes — all 6 libraries | 205 clones | A D E G I K L M N Q S T V |
| E141 | yes — all 6 libraries | 130 clones | A D G S T |

### Three genuine gaps this project occupies

1. **R79 is completely unexplored.** It appears in no library design and in no
   characterised clone. Given §2 shows it is a cavity-rim residue with clear
   line of sight to the ligand, and §3 shows R79A/E94A unlocks +84.8 Å³, this
   is the single most under-sampled position in the pocket.

2. **F108 was only ever made LARGER.** Of 42 clones mutating F108, **39 (93%)
   are F108W** — bigger than Phe (228 vs 190 Å³). Only 3 clones carry a smaller
   residue (Q, L, H). Smaller substitutions *were* permitted by the DSM-Hao and
   TSM library designs but were essentially never recovered, consistent with
   F108 being load-bearing (§2, §4) — or simply not useful for their ligands.

3. **K59 and F108 were never simultaneously reduced.** 29 clones mutate both;
   **zero** reduce both in volume. Every one is the K59{A,Q} + F108W pattern —
   shrink one side, enlarge the other. The `K59A_F108A` variant in our panel
   (rigid +127.9 Å³, post-relax **+129.3 Å³**) is therefore untested territory.

### Two results that de-risk this project

- **The E94–R79 salt bridge is dispensable for function.** E94 is mutated in
  189 characterised clones, and **132 of those go to residues that cannot
  salt-bridge R79** (G 66, A 39, S 22, T 4, L 1). These are functional sensor
  hits, so breaking the bridge from the E94 side demonstrably preserves the
  switch. That is direct experimental evidence bearing on the main risk in
  Arm 2.
- **Large wall-volume release is tolerated.** 334 of 419 clones show a net
  side-chain volume *decrease* across the wall positions, the largest being
  **−267 Å³** (`G_31.3`: K59M/E94S/Y120G/E141G). Notably the entire top-12 by
  volume released are **PFAS binders** — long-chain perfluorinated compounds —
  i.e. the pocket was hollowed to accommodate an elongated ligand.

**Framing:** their libraries reprogram specificity *within* the existing cavity
envelope by side-chain swapping, and much of the naive "remove bulk" space is
already explored. What remains genuinely novel is the **R79 axis** and its
**combination with F108 reduction** — the moves that open the sealed
sub-pocket identified in §2.

### Was R79 ever in a PYR1 library? A survey of the engineering literature

| Study | Library scope | R79 included? |
|---|---|---|
| **Mosquna 2011** PNAS ([10.1073/pnas.1112838108](https://doi.org/10.1073/pnas.1112838108)) | site-saturation at **39 conserved residues** contacting ABA (LIG) or PP2C (PPI); all 741 single substitutions. Full list in `data/mosquna2011_ssm_residues.tsv`: 55, 59, 60, 61, 62, 63, 81, 83, 84, 85, 86, 87, 88, 89, 92, 94, 108, 110, 115, 116, 117, 120, 122, 141, 148, 150, 151, 154, 155, 156, 158, 159, 160, 162, 163, 164, 166, 167, 170 | **no** |
| **Park 2015** Nature ([10.1038/nature14123](https://doi.org/10.1038/nature14123)) | pocket library → PYR1^MANDI | **no** |
| **Beltrán 2022** Nat Biotechnol ([10.1038/s41587-022-01364-5](https://doi.org/10.1038/s41587-022-01364-5)) | **19 DSM positions**: 59, 62, 81, 83, 87, 89, 92, 94, 108, 110, 117, 120, 122, 141, 159, 160, 163, 164, 167 | **no** |
| **Park 2024** Nat Chem Biol ([10.1038/s41589-023-01447-7](https://doi.org/10.1038/s41589-023-01447-7)) | NNK single-site saturation across the **entire coding sequence** | **formally yes** (whole-ORF), but R79 appears nowhere in the results |
| **Tian 2025** PNAS (this project's §4b) | 18 positions across 6 libraries | **no** |

### R79 was selected by the rule and removed by hand — the reason is on record

The Beltrán supplementary methods give the library-design filter verbatim:

> We examined PYR1 bound to ABA (PDB ID 3K90, chain A) and selected all
> residues either:
> - with any atom within 5 Å of ABA, or
> - with any atom within 6 Å of ABA and a Cα–Cβ vector not directed away from ABA.
>
> After manual curation of the resulting list of positions, we removed the
> following residues:
> - **H60**: Faces away from ABA and involved in PYR1 dimerization
> - **F61**: Outward facing and involved in PYR1 dimerization
> - **R79**: *Appears to have only second-shell effects on ligand binding*
> - **P88**: Conserved in the gate loop
> - **T91**: Outward facing
> - **H115**: Conserved in the latch loop

`scripts/13_beltran_criterion.py` reimplements this filter on 3K90 chain A
(→ `logs/13_beltran_criterion.txt`). It reproduces their set **exactly**:

| | n | residues |
|---|---|---|
| rule 1 (≤ 5.0 Å) | 22 | 59, 60, 61, 81, 83, 87, 88, 89, 91, 92, 94, 108, 110, 115, 117, 120, 141, 159, 160, 163, 164, 167 |
| rule 2 (≤ 6.0 Å, Cα–Cβ toward ABA) | 3 | **79**, 116, 122 |
| **total** | **25** | matches their stated "25 residues … close contact to ABA" |

**R79 passes the automated filter** — 5.89 Å from ABA in 3K90, with a Cα–Cβ
cosine of **+0.294, i.e. pointing toward the ligand**. It entered the candidate
set on geometry and was then removed by hand.

So R79 was neither overlooked nor excluded by a distance threshold. It was
excluded by a **stated judgement**: *"Appears to have only second-shell effects
on ligand binding."* That is exactly the claim this project's measurements
contradict (§2):

- R79 **NH2 and NE have clear line of sight to ABA** — no intervening PYR1
  atom over 5.81 Å (3QN1).
- R79 borders the WT cavity and becomes one of its dominant lining residues
  once the wall is removed (4 → 3374 lining grid points).
- R79A/E94A unlocks **+84.8 Å³**, more than either K59A or F108A alone.

Checked independently, five of the six removals hold up well — and this project
had already excluded four of them for its own reasons:

| Residue | min dist to ABA | Cα–Cβ orientation | independent check |
|---|---|---|---|
| T91 | 4.34 Å | **away** (−0.968) | rigid truncation ΔV = **+0.0 Å³** — "outward facing" confirmed |
| F61 | 3.83 Å | toward (+0.903) | ΔV = **−16.5 Å³** — a *sealing* residue; excluded here too |
| P88 | 3.15 Å | toward (+0.317) | gate loop — excluded here too |
| H115 | 4.19 Å | toward (+0.848) | latch loop — excluded here too |
| H60 | 4.76 Å | toward (+0.019) | dimerisation interface; not a cavity residue |
| **R79** | **5.89 Å** | **toward (+0.294)** | **the one removal this project's data contradicts** |

*(Two minor discrepancies against their final 19: position 62 appears in their
list but not in this reproduction, and 116 here but not theirs. Their curation
evidently adjusted those two as well. Neither affects the R79 conclusion, which
turns on rule 2 alone.)*

### Mosquna 2011 excluded R79 by the hard 5 Å rule — verified

Mosquna *et al.* state their criterion explicitly:

> *"Thirty-nine residues within 5 Å of ABA, four water-molecules that contact
> ABA, or the PP2C HAB1 (Fig. 1A) were identified using PYR1 structure
> coordinates, and all possible 741 single amino acid substitutions at these
> sites were constructed by site-directed mutagenesis."*

So the rule is a union: **(a)** within 5 Å of ABA, **(b)** within 5 Å of one of
the ordered waters contacting ABA, or **(c)** contacts HAB1. Condition (b) is a
genuine second chance for R79 — a bridging water could have brought it in.

`scripts/14_mosquna_criterion.py` (→ `logs/14_mosquna_criterion.txt`) tests
this on both crystal forms:

| | R79 → ABA | (a) ≤ 5 Å of ABA | (b) ≤ 5 Å of an ABA-contacting water |
|---|---|---|---|
| 3K90 (their coordinates) | 5.89 Å | **fail** | **fail** |
| 3QN1 (ternary complex) | 5.81 Å | **fail** | **fail** |

The water rule works as described and is not vacuous — it is what brings **116
and 122** into the set (both present in their Table S1) without either being
within 5 Å of ABA. But **no ABA-contacting water bridges to R79**, so condition
(b) does not rescue it.

**R79 was therefore legitimately excluded by Mosquna's stated criterion** — no
oversight, no inconsistency. The exclusion traces entirely to the hard 5 Å
ligand-distance cutoff.

### The two exclusions together

| Study | Rule | R79 outcome |
|---|---|---|
| **Mosquna 2011** | ≤ 5 Å of ABA / ABA-water / HAB1 | **excluded by the rule** (5.89 Å, no bridging water) |
| **Beltrán 2022** | ≤ 5 Å of ABA, **or** ≤ 6 Å with Cα–Cβ toward ABA | **selected by the rule**, then removed by hand as "second-shell" |

Both exclusions are defensible on their own terms, and they fail in the same
place: each is anchored on **distance to ABA**. Mosquna's cutoff is simply too
tight to reach R79; Beltrán's is wide enough to catch it, at which point the
manual "second-shell" judgement removed it. The result is that across the whole
PYR1 engineering literature, R79 has been randomised only incidentally, in
Park 2024's whole-ORF single-site saturation, and **never in a combinatorial
pocket library**.

**Bottom line:** excluding R79 rests on a reasonable but untested structural
intuition — that a residue not directly touching ABA cannot matter for ligand
binding. That intuition is sound for *reprogramming specificity within* the
existing cavity, which is what every library above was built to do. It fails
for *enlarging* the cavity, because ABA does not fill the pocket, so
"second shell relative to the ligand" and "second shell relative to the cavity"
are different things (§2). **R79 is the former but not the latter** — and that
distinction is the opening this project is built on.

### ⚠️ Correction: F108A is already a validated, functional variant

The §4b claim that "F108 was only ever made larger" holds **only within the
Tian datasets**. It is not true of the field. Sequence comparison of **PDB
4WVO** ("An engineered PYR1 mandipropamid receptor in complex with
mandipropamid and HAB1") against WT PYR1 gives, directly from the structures:

> **PYR1^MANDI = K59R / V81I / F108A / F159L**

So **F108A is present in the flagship engineered receptor**, in a
crystallographically characterised, functional CID module. This is important
for two reasons:

1. **F108A demonstrably preserves the ratchet.** The main risk in Arm 2 —
   that removing the F108 buttress breaks the switch — is already answered
   experimentally for the single mutant.
2. **K59 and F108 have been co-mutated in a working receptor** (K59R + F108A),
   though K59R is a size-neutral swap rather than a volume-opening one. The
   untested move remains reducing **both** simultaneously (`K59A_F108A`,
   post-relax **+129.3 Å³**), and doing so with **R79** released.

4WVO is also the correct calibration anchor for §8 step 2: any scoring
pipeline used here must rank PYR1^MANDI + mandipropamid above WT PYR1 +
mandipropamid before its numbers on novel pockets can be trusted.

---

## 4c. Cavity SHAPE, carotenoids, and the ligand panel

`scripts/15_carotenoid_start.py`, `scripts/16_cavity_shape.py`,
`scripts/17_ligand_ladder.py` → corresponding files in `logs/`

### The START fold does bind intact carotenoids

**ABA is itself an apocarotenoid** — the NCED cleavage product of
9-*cis*-violaxanthin/neoxanthin — so PYR1 already binds a carotenoid *fragment*.
The ancestral hypothesis does not require the fold to gain a new capability,
only to have held the uncleaved precursor.

And the fold demonstrably does: ***Bombyx mori* carotenoid-binding protein
(BmCBP)** is a genuine START-domain protein (STARD3/MLN64-orthologous) that
binds lutein *in vivo* and produces the yellow-cocoon phenotype. Structures:
**7ZVR** (zeaxanthin complex), 7ZTQ apo, 7ZTR/7ZTU/7ZVQ apo point mutants,
8AAQ. Reference: *Structural basis for the carotenoid binding and transport
function of a START domain*, Structure (2022), PMID 36356587.

BmCBP also positions its Ω1 loop over the cavity via an **R173–D279 salt
bridge** — an arginine salt bridge controlling a loop above the binding site,
architecturally analogous to PYR1's R79–E94 (§2). That is a second independent
instance of the motif in this fold.

> ⚠️ **Zeaxanthin is only half-modelled in 7ZVR.** It is C₄₀H₅₆O₂ = 42 heavy
> atoms, but only 20 are deposited (C1–C20 + O3 — one β-ionone ring and the
> polyene to ~C14/C20) at occupancies 0.60–0.65. The distal half is disordered.
> Cavity volumes seeded on the modelled fragment are therefore a LOWER BOUND.
> **This is also a positive result for design: a C40 carotenoid does not
> require a fully enclosed 28 Å tunnel — a buried ionone ring plus polyene with
> a solvent-exposed distal tail is the actual binding mode in the only
> START-fold carotenoid complex solved.**

Volumes on this project's grid criterion: PYR1 166.4 Å³, BmCBP holo (7ZVR)
238.5 Å³, BmCBP apo (7ZTQ) 367.4 Å³.

### Merging the two lobes elongates the cavity — measured

Volume alone does not decide whether a polyene fits; the **longest internal
straight segment** does. `16_cavity_shape.py` computes, for each truncation
set, the longest straight path through the cavity holding ≥1.6 Å clearance,
plus PCA extents and the narrowest cross-section along that axis:

| variant | volume Å³ | **span Å** | PCA extents Å | r_med | **r_min** | max R |
|---|---|---|---|---|---|---|
| WT | 167.0 | **10.6** | 9.7 / 7.3 / 5.3 | 4.99 | 1.25 | 3.22 |
| K59A | 244.8 | 16.0 | 18.6 / 7.0 / 5.7 | 4.43 | 3.20 | 3.22 |
| F108A | 245.6 | 12.2 | 19.6 / 7.2 / 5.8 | 5.02 | 3.08 | 3.22 |
| K59A/F108A | 303.6 | 16.0 | 18.0 / 6.9 / 5.9 | 4.95 | 3.36 | 3.22 |
| R79A/E94A | 256.5 | 15.5 | 19.1 / 8.0 / 5.5 | 6.04 | 2.07 | 3.38 |
| F108A/R79A/E94A | 315.8 | 15.5 | 18.5 / 7.8 / 5.9 | 6.16 | 2.07 | 3.38 |
| **K59A/F108A/R79A/E94A** | **391.1** | **17.9** | 16.7 / 7.5 / 6.1 | 5.30 | 2.09 | 3.95 |
| quad + Y120A | 441.5 | 17.9 | 15.8 / 9.4 / 6.3 | 5.95 | 2.09 | 3.95 |
| quad + Y120A/E141A | **478.8** | 17.9 | 15.5 / 10.1 / 6.9 | 7.43 | 2.09 | 3.95 |

**Merging the lobes does elongate the pocket:** span 10.6 → 17.9 Å (+69 %),
principal long axis 9.7 → 16.7 Å, aspect ratio ~1.8 → ~2.7. The cavity stops
being a ball and becomes a channel. So the earlier "PYR1 is globular, carotenoids
need a tunnel" objection was drawn from **WT** geometry and does not survive
contact with the expanded geometry.

*Metric calibration:* WT span is 10.6 Å and the bound ABA conformer measures
10.45 Å end-to-end — ABA fits WT with ~0.2 Å to spare. The metric is therefore
well scaled, not arbitrarily strict.

**The bottleneck, not the length, is now the limiting feature.** Every
R79A/E94A-containing set has `r_min` ≈ 2.1 Å, i.e. the merged channel is pinched
mid-way — narrower than the ~3.3 Å half-width of a xanthophyll. Widening that
neck should become an explicit design objective alongside volume; the residues
lining it are the ones recruited on merging (I48, V49, R50, V81, I62, S122, §2).

### Ligand panel (`17_ligand_ladder.py`)

MMFF-minimised geometries measured the same way as the cavity:

| ligand | MW | heavy atoms | length Å | half-width Å | rot. bonds | fits quad span (17.9 Å) |
|---|---|---|---|---|---|---|
| β-ionone | 192.3 | 14 | 7.6 | 2.67 | 2 | ✅ |
| GR24 (strigolactone) | 248.3 | 18 | 9.2 | 2.45 | 2 | ✅ |
| ABA *(reference)* | 264.3 | 19 | 11.4 | 2.70 | 3 | ✅ |
| retinal | 284.4 | 21 | 14.1 | 2.69 | 5 | ✅ |
| retinoic acid | 300.4 | 22 | 14.5 | 2.67 | 5 | ✅ |
| mandipropamid *(4WVO anchor)* | 442.9 | 31 | 15.5 | 5.88 | 10 | ✅ |
| crocetin | 328.4 | 24 | 20.7 | 2.21 | 8 | ✗ (−2.8 Å) |
| β-apo-8′-carotenal | 416.6 | 31 | 24.0 | 2.81 | 9 | ✗ (−6.1 Å) |
| β-carotene (C40) | 536.9 | 40 | 27.3 | 3.92 | 10 | ✗ (−9.4 Å) |
| zeaxanthin (C40) | 542.8 | 40 | 27.8 | 3.35 | 9 | ✗ (−9.9 Å) |

**C40 carotenoids are included in the panel deliberately**, despite failing the
closed-state span test, for three reasons:

1. **The span test assumes full enclosure, which BmCBP shows is not required.**
   Half of zeaxanthin is disordered and unmodelled in 7ZVR — the physiological
   mode is a buried ionone ring plus polyene with a protruding tail. PYR1 needs
   to bury enough of the ligand to trigger gate closure, not all of it.
2. **These are rigid-backbone, closed-state numbers.** The protocol constrains
   Cα positions by design (§5) and cannot open a channel; real backbone
   relaxation and loop motion will lengthen the achievable span.
3. **The 4Å-scale deficit is the right size of target.** Going from 10.6 to
   17.9 Å by side-chain truncation alone already covers the apocarotenoid range;
   closing the remaining ~10 Å is exactly what the insertion/grafting arm (§1)
   is for.

Tiering for screening: **tier 1** (β-ionone → retinoic acid) is reachable now
and is where hits should be expected; **tier 2** (crocetin, β-apo-8′-carotenal)
is the stretch target for the merged-and-widened cavity; **tier 3** (zeaxanthin,
β-carotene) tests the ancestral hypothesis directly and should be scored with a
partially-exposed-tail binding mode in mind, not rejected on enclosed volume.

### ⚠️ 4d. The bottleneck IS the gate and latch — a structural conflict

`scripts/18_exits_and_bottleneck.py` → `logs/18_exits_bottleneck.txt`

Two results here materially change the outlook for oversized ligands, and both
argue for caution about tier 3.

**(a) The expanded cavity is essentially sealed — there is no tail exit.**
Searching for "mouths" (free, low-buriedness grid points bordering the cavity)
in the K59A/F108A/R79A/E94A cavity found **3 points** — i.e. no opening of any
meaningful size. Merging the lobes produces a larger *enclosed* chamber, not a
channel to the surface. A ligand too long to fit has nowhere to put its tail.

**(b) The constriction sits at the gate/latch end of the channel.**
Profiling cross-section along the 17.9 Å span axis, the three narrowest slices
(radius 1.97–2.64 Å, at t = +15.5 to +18.1 Å, i.e. one END of the span) are
lined by:

| rank | residue | weight | note |
|---|---|---|---|
| 1 | **F159** | 64 | C-terminal helix; a *sealing* residue (§3, ΔV = −14.5 Å³) |
| 2 | **H115** | 62 | **LATCH** |
| 3 | **A89** | 61 | **GATE** |
| 4 | **L117** | 55 | **LATCH** |
| 5 | **P88** | 51 | **GATE** |
| 6 | **R116** | 33 | latch-adjacent |
| 7 | L87 | 21 | gate-adjacent; sealing residue (ΔV = −16.1 Å³) |
| 8 | V83 | 20 | β-strand |

**Six of the eight bottleneck residues are gate, latch, or sealing residues.**
Widening the neck therefore means mutating exactly the elements the readout
depends on — the region this project declared off-limits in §1. This is a
genuine structural conflict, not a tractable engineering step:

- the merged channel runs from the newly opened satellite lobe at one end **to
  the gate/latch at the other**;
- an elongated ligand long enough to span it would press its far end **into the
  gate–latch groove**, which is precisely where HAB1 W385 must wedge;
- so the exposed-tail binding mode that works for BmCBP is not available to
  PYR1 without sacrificing the ratchet.

*Caveat:* both measurements are static, rigid-backbone, and use this project's
deliberately conservative buriedness criterion (0.88), which under-reports
marginal openings by construction. A ligand could in principle force a
transient opening that this analysis cannot see. The finding is strong enough
to set priorities, not strong enough to close the question — testing it is
exactly what tier 3 is for.

**Revised recommendation on tier 3.** Keep the C40 carotenoids in the panel,
but reframe their purpose: they are a **diagnostic**, not a target. Score them
in Arm 2 (ternary complex) specifically to measure whether the ligand
(i) occupies the gate–latch groove, (ii) displaces or competes with W385, and
(iii) degrades `dG_separated`. Expect them to fail; the value is in learning
*how* they fail and confirming the geometric argument above with energetics.

**And this reopens the deepening/insertion arm as the real route.** If the only
constriction is at the readout end, then the correct way to lengthen the pocket
is the opposite direction — **extend the cavity away from the gate/latch**, at
the distal β-hairpins and the helix N-cap (§1, "deepen"). That adds span
without touching the readout and would also create a tail exit on the face
*opposite* HAB1. The bottleneck result therefore promotes the insertion/grafting
arm from "future work" to the primary path for anything larger than tier 2.

### 4e. There is NO inter-lobe bottleneck — the limit is LENGTH, and it is backbone-bound

`scripts/19_interlobe_neck.py` → `logs/19_interlobe_neck.txt`

Profiling cross-sectional radius in fine bins along the whole span of the
merged (quad) cavity gives a flat, wide channel with a single narrow tail:

| region | t (Å) | radius (Å) |
|---|---|---|
| main chamber | +0.4 → +16.0 | **4.0 – 6.3** (no constriction) |
| gate/latch tail | +16.5 → +19.4 | 1.4 – 2.7 |

**The two lobes have fully merged into one continuous chamber.** The
"inter-lobe neck" detector returns t = +9.0 Å at radius **5.02 Å**, which is not
a constriction at all — it is the shallow local minimum of a broad plateau. All
twelve residues lining it are either already truncated in the quad (K59, F108,
E94, R79) or freely mutable (I110, N167, Y120, I62, S92, V163, V81, S122); none
are gate, latch, or sealing.

**Width is not the problem.** The chamber's radius (4–6 Å) already exceeds a
xanthophyll's ~3.35 Å half-width along its entire length. A C40 polyene is thin;
the merged pocket is wide enough for it.

**Length is the problem, and side-chain truncation has saturated:**

| variant | volume Å³ | **span Å** | neck radius Å | gap to C40 |
|---|---|---|---|---|
| quad | 391.8 | **18.6** | 5.02 | **+9.2 Å** |
| quad + neck(I110,N167,Y120,I62) | 513.9 | 16.9 | 4.48 | +10.9 Å |
| quad + Y120A/E141A | 480.1 | 17.9 | 6.21 | +9.9 Å |
| quad + Y120A/E141A + neck | **550.1** | 17.0 | 5.69 | +10.8 Å |

Volume keeps climbing (391 → **550 Å³**, 3.3× WT) but **span does not increase —
it slightly decreases.** Removing more side chains makes the chamber rounder and
fatter, not longer. The cavity's ends are bounded by **backbone**, not by side
chains, so mutation alone cannot extend it further.

### Verdict on complete C40 enclosure

| | length |
|---|---|
| merged cavity span (best achieved) | **18.6 Å** |
| β-carotene | 27.3 Å |
| zeaxanthin | 27.8 Å |
| **shortfall** | **~9–9.2 Å (the cavity is ~67 % of the required length)** |

**A C40 carotenoid cannot be fully enclosed by any side-chain-only design.**
The deficit is ~9 Å of *length*, and it is backbone-limited. Closing it requires
roughly 2–3 β-strand register units (+2 residues each on paired strands, §1) or
a hairpin graft from a large-cavity homolog — i.e. **the insertion/grafting arm
is not optional for tier 3, it is the only route.**

This sharpens the ancestral hypothesis into something testable rather than
decorative: if an ancestral SRPBCC member bound an intact carotenoid, it would
have needed a **longer** cavity, so modern PYR1 would represent a *shortening*.
Arm 3 can test this directly — do the large-cavity helix-grip homologs achieve
their volume by longer β-hairpins/elongated helix-grip geometry, and if so, at
which structural elements? That indel map is precisely the graft donor set.

> **Note on BmCBP provenance:** *Bombyx mori* is the domestic **silkworm** — a
> lepidopteran moth, not a spider. Its carotenoid-binding mechanism (Ω1 loop
> pinned by the R173–D279 salt bridge, half the ligand disordered) is one
> solution; an ancestral plant-lineage carotenoid binder need not have used the
> same one, and the fully-enclosed mode hypothesised here would in fact be a
> *distinct* mechanism.

---

## 5. The pilot: three arms

### Arm 1 — does the cavity stay open? (`scripts/06_rosetta_cavity_scan.py`)

Rigid truncation is an upper bound. The real question is whether an enlarged
cavity survives repacking or is simply refilled by second-shell side chains.
That is a **repacking** question, which is why this uses Rosetta FastRelax.

> **Why not ESMFold/AF2?** Structure predictors are close to blind to point
> substitutions: they return near-identical backbones at high confidence for
> every member of this panel, including variants that would collapse. Worse,
> MSA-based predictors actively "correct" toward the family consensus. They
> cannot answer this question. (This also applies to any future insertion
> scanning — use single-sequence mode there.)

- Input `data/pyr1_A.pdb` — PYR1 chain A, **apo**, closed conformation.
  Apo is deliberate: the ligand is removed so the cavity is free to close.
- Mutate → FastRelax (`ref2015_cst`, 3 repeats) with CA coordinate constraints
  to starting coordinates (`-relax:constrain_relax_to_start_coords`).
  Constraints hold the closed backbone so the measurement isolates side-chain
  infilling; side chains are unconstrained and free to collapse inward — that
  collapse is the signal.
- Re-score with plain `ref2015`; measure cavity with `lib_cavity` seeded at the
  crystallographic ABA centroid; record CA-RMSD to input as a control that the
  backbone was in fact held.
- WT is relaxed through the identical protocol and is the reference for both
  ΔREU and ΔVolume, so the comparison is internally controlled.

### Arm 2 — does the switch survive? (`scripts/07_rosetta_ratchet.py`)

Because R79/E94/K59/E141 are PYR/PYL-invariant, the dominant risk is not a
collapsed pocket but a **broken ratchet**. Run on the full ternary complex:

- Input `data/complex_AB_ABA.pdb` — PYR1 (A) + HAB1 (B) + ABA (X).
- Same mutate → FastRelax protocol, then `InterfaceAnalyzerMover` across the
  `A_B` jump with `pack_separated=True`.
- Readouts: `dG_separated` (interface energy — the ratchet proxy), `dSASA`,
  shape complementarity `sc_value`, gate and latch heavy-atom RMSD to input,
  and the **HAB1 W385 → ABA minimum distance** (the wedge that locks the closed
  state).

Decision table:

| Arm 1 | Arm 2 | Interpretation |
|---|---|---|
| cavity opens | dG_separated ≈ WT | **candidate** — bigger pocket, intact switch |
| cavity opens | dG_separated much worse | ratchet broken; pocket gain is useless |
| cavity opens | dG_separated much better | possible **constitutive** binder (ligand-independent false positive) — check the apo arm before trusting |

### Arm 3 — is there a big-pocket signature? (`scripts/08_homolog_cavities.py`)

Superfamily cavity volume spans ~400 Å³ (PYR/PYL) to >1000 Å³ (START domains)
**in the same fold** — an existence proof that this scaffold tolerates 2–3×
the PYR1 cavity. For the top homologs by alnTM, this script extracts structures
from the foldseek CATH50 database, measures the largest enclosed cavity with an
unbiased grid search (same `lib_cavity` criterion), and tests per position
whether residue identity separates large- from small-cavity domains.

If large-cavity members share an identity signature at 59/79/94/108/120/141,
**that signature is the design prescription — read off evolution rather than
guessed.**

Note: Arm 3 seeds the cavity search unbiasedly (largest component anywhere in
the domain) whereas Arms 1–2 seed at the ligand centroid. Volumes are
comparable *within* Arm 3 but only qualitatively comparable to Arm 1.

### Variant panel (`scripts/variants.py`, 41 including WT)

1. **Single-substitution scan** — 6 positions (59, 79, 94, 108, 120, 141) ×
   {A, G, S, V, L} = 30 variants.
2. **Superfamily-guided singles** — `F108Y` (consensus is *larger* than WT),
   `K59D` (consensus residue at 59), `E141W` (distant-subfamily signature).
3. **Wall-breaks** — `R79A_E94A`, `F108A_R79A`, `K59A_F108A`,
   `F108A_R79A_E94A`, `K59A_F108A_R79A_E94A`. These test directly whether F108
   needs its buttress removed to matter.
4. **Superfamily-signature combinations** — `K59L_F108Y`, `K59L_F108Y_E141W`.

Deliberately excluded: gate (S85–A89), latch (H115–L117), and the sealing
residues F61 / F159 / L87.

---

## 6. Methods detail for reproducibility

### Software versions

| Tool | Version / identifier |
|---|---|
| PyRosetta | **2026.06+release.1a56185c2592611dec4c9c75ddc9468cd2227c1f** (2026-01-30), `PyRosetta4.conda.ubuntu.cxx11thread.serialization.Ubuntu.python311.Release` |
| Rosetta (for `molfile_to_params.py` only) | 2023.45, module `rosetta/2023.45`, `$ROSETTA3=/opt/linux/rocky/8.x/x86_64/pkgs/rosetta/2023.45` |
| Rosetta score function | `ref2015` (scoring), `ref2015_cst` (relax) |
| foldseek | commit **941cd33ff0771cd2e3f144e3293e22a2b87e9fda**, binary `/bigdata/cutlerlab/jjaco081/tools/foldseek/foldseek/bin/foldseek` |
| foldseek CATH50 database | downloaded **2024-01-25** |
| Python (analysis) | 3.12.13 — env `/bigdata/cutlerlab/jjaco081/conda_envs/esmfold2` |
| Python (PyRosetta) | 3.11.15 — env `/bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis` |
| NumPy / SciPy / Biopython | 2.5.0 / 1.18.0 / 1.87 |
| RDKit | 2026.03.3 |

### Structural input

PDB **3QN1** — *Crystal structure of the PYR1 Abscisic Acid receptor in complex
with the HAB1 type 2C phosphatase catalytic domain*. Chain A = PYR1
(residues 1–181, 179 resolved), chain B = HAB1 (185–505, 295 resolved),
ligand **A8S** = (+)-abscisic acid, plus 3 Mn²⁺ and 465 waters (all discarded).
Alternate locations: only altloc `' '` or `'A'` retained. Hydrogens discarded
and re-added by Rosetta.

### Ligand preparation — protonation state (important)

The RCSB *ideal* SDF for A8S is the **neutral carboxylic acid**. In 3QN1 the
ABA carboxylate salt-bridges K59 NZ at 2.85 Å, so the bound species is the
**deprotonated carboxylate, net charge −1**. Because K59 is one of the mutated
positions, this charge assignment materially affects the K59 energetics.
`scripts/00_setup_inputs.py` therefore:

1. reads crystal heavy-atom coordinates from 3QN1;
2. assigns bond orders from the ideal-SDF template
   (`RDKit AssignBondOrdersFromTemplate`);
3. deprotonates the carboxylic acid (SMARTS `[CX3](=O)[OX2H1,OX1H0-]`,
   asserted unique) and sets formal charge −1;
4. adds hydrogens with coordinates (`Chem.AddHs(addCoords=True)`);
5. generates Rosetta params with `molfile_to_params.py -n A8S`.

`molfile_to_params.py` emits `A8S_0001.pdb` carrying both Rosetta-consistent
atom names and the input coordinates; that file is used directly as the ligand,
which sidesteps any PDB↔Rosetta atom-name mismatch. **Verified: 19 heavy atoms
in, 19 out, centroid shift 0.0000 Å, max per-atom deviation 0.0000 Å** — the
crystallographic pose is preserved exactly.

### Cavity volume algorithm (`scripts/lib_cavity.py`)

1. Cubic grid, 0.5 Å spacing (0.6 Å in Arm 3), over a ±14 Å box centred on the
   crystallographic ABA centroid.
2. A grid point is **free** if its distance to the nearest heavy atom centre
   exceeds that atom's van der Waals radius + a **1.4 Å** solvent probe
   (C 1.70, N 1.55, O 1.52, S 1.80, P 1.80 Å).
3. A free point is **buried** if ≥ **0.88** of **26** evenly distributed ray
   directions strike protein within 15 Å (0.75 Å march step). This separates
   interior cavity from bulk solvent and surface grooves.
4. Free ∧ buried points are connected-component labelled
   (`scipy.ndimage.label`); the component containing — or nearest to — the seed
   point is reported. Volume = *n* voxels × spacing³.

⚠️ **Absolute volumes are conservative relative to fpocket/CASTp** because of
the strict buriedness cut-off and single-component restriction. They are valid
for **relative comparison between structures processed identically** and should
not be quoted against literature values computed with another tool.

### Known caveats

- **ΔREU is not a benchmarked ΔΔG.** It is the difference of independently
  relaxed total scores (no `cartesian_ddg`, no reference-state correction).
  Use it to flag grossly destabilising variants, not to rank near-neutral ones.
- **`dG_separated` is not an affinity.** It is a Rosetta interface energy on a
  relaxed model — a rank-order triage signal only.
- **Backbone motion is restrained.** CA coordinate constraints mean the pilot
  measures side-chain infilling, not backbone remodelling. Genuinely large
  cavity gains may require backbone changes this protocol cannot sample; those
  belong to the insertion/grafting arm.
- **Both arms model the CLOSED state only.** Neither tests whether the apo
  variant can still *open*. That requires the apo-open state or MD, and is the
  main gap in the current design.

### Cluster execution notes

- **`/scratch` is node-local xfs and invisible to compute nodes.** A SLURM job
  with any path under `/scratch` fails at 00:00:00 with no log written. All
  job files live under `/bigdata` (GPFS).
- **The login node enforces a 1 GB per-user memory cap.** Both foldseek
  (~1.15 GB peak) and FastRelax exceed it and are silently OOM-killed
  (foldseek reports only `Kmer matching step died`). Everything must go through
  SLURM. PyRosetta *import and init* fit under the cap; FastRelax does not.
- Jobs run on partition **`cutlerlab`**.

---

## 7. Directory layout

```
PYR1_pocket_expansion/
├── README.md                      this file
├── data/
│   ├── 3QN1.cif                   source structure (RCSB)
│   ├── A8S_ideal.sdf              RCSB ideal ligand (neutral acid)
│   ├── aba_deprot.sdf             deprotonated carboxylate, crystal pose
│   ├── aba_xtal.pdb               ABA heavy atoms, crystal pose
│   ├── aba_params/A8S.params      Rosetta ligand definition (charge −1)
│   ├── aba_params/A8S_0001.pdb    Rosetta-named ligand, crystal coordinates
│   ├── pyr1_A.pdb                 PYR1 chain A, apo, closed  [Arm 1 input]
│   ├── hab1_B.pdb                 HAB1 chain B
│   ├── complex_AB.pdb             PYR1 + HAB1
│   └── complex_AB_ABA.pdb         ternary complex            [Arm 2 input]
├── scripts/
│   ├── lib_cavity.py              grid cavity-volume library (shared)
│   ├── variants.py                the 41-member panel + per-position rationale
│   ├── 00_setup_inputs.py         structure/ligand preparation
│   ├── 01_pocket_contacts.py      pocket definition + buttress network
│   ├── 02_rigid_truncation_scan.py  rigid-backbone upper bounds
│   ├── 03_parse_foldseek.py       superfamily conservation at wall positions
│   ├── 05_foldseek_search.sh      SLURM: foldseek vs CATH50
│   ├── 06_rosetta_cavity_scan.py  ARM 1
│   ├── 07_rosetta_ratchet.py      ARM 2
│   └── 08_homolog_cavities.py     ARM 3
├── results/
├── logs/
```

### Reproducing from scratch

```bash
P=/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
PY_ANALYSIS=/bigdata/cutlerlab/jjaco081/conda_envs/esmfold2/bin/python
PY_ROSETTA=/bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python

cd $P/data && curl -O https://files.rcsb.org/download/3QN1.cif \
  && curl -o A8S_ideal.sdf https://files.rcsb.org/ligands/download/A8S_ideal.sdf

$PY_ANALYSIS $P/scripts/00_setup_inputs.py
$PY_ANALYSIS $P/scripts/01_pocket_contacts.py
$PY_ANALYSIS $P/scripts/02_rigid_truncation_scan.py
sbatch        $P/scripts/05_foldseek_search.sh          # must be SLURM
$PY_ANALYSIS $P/scripts/03_parse_foldseek.py
sbatch        $P/scripts/09_submit_pilot.sh             # Arms 1+2, array
$PY_ANALYSIS $P/scripts/08_homolog_cavities.py          # Arm 3
```

---

## 5b. Pilot results — Arm 1 succeeded, Arm 2 is blind

Jobs **27260956** (pilot array 0–40) and **27261151** (Arm 3) both completed
2026-08-06/07: 41/41 JSONs in each arm, 147 homolog cavities in Arm 3.
Aggregated by **`scripts/20_aggregate_pilot.py`** →
`results/20_pilot_summary.csv` and `results/20_pilot_summary.md`.

WT reference (mean of 3 relax replicates): cavity **174.4 Å³**, apo score
**−510.58 REU**, `dG_separated` **−67.30**, `dSASA` **1797 Å²**, `sc` **0.754**,
gate RMSD 0.18 Å, latch RMSD 1.23 Å, W385–ABA **4.60 Å** (crystal 4.67 Å).

### ⚠️ Read this before reading the table: Arm 2 does not discriminate

The aggregator computes, for every metric, the spread **between** variant means
against the mean **within**-variant replicate spread. A metric only carries
information if that ratio exceeds ~2.

| metric | arm | within sd | between sd | ratio | |
|---|---|---|---|---|---|
| `cavity_A3` | 1 | 5.47 | 47.11 | **8.61** | informative |
| `score_ref2015` | 1 | 1.53 | 4.24 | **2.78** | informative |
| `dG_separated` | 2 | 2.78 | 2.04 | 0.73 | ❌ noise ≥ signal |
| `dSASA` | 2 | 11.85 | 9.89 | 0.83 | ❌ |
| `sc_value` | 2 | 0.006 | 0.004 | 0.69 | ❌ |
| `gate_rmsd` | 2 | 0.004 | 0.002 | 0.65 | ❌ |
| `latch_rmsd` | 2 | 0.229 | 0.212 | 0.92 | ❌ |
| `w385_aba_min_dist` | 2 | 0.025 | 0.046 | 1.84 | ❌ |

**Every Arm-2 metric fails.** The raw table appears to say "0 of 41 variants
break the readout" — that is *not* what happened. The assay cannot tell. Two
causes, both real and both fixable:

1. **nstruct = 3 is far too few for a Rosetta interface energy.** The replicate
   sd of `dG_separated` is ~2.8 REU, larger than the entire between-variant
   range. The apparent "improvements" of −5 REU in the table are noise.
2. **`-relax:constrain_relax_to_start_coords` pins every Cα to the input
   complex.** Gate and latch Cα RMSD are therefore bounded *by construction* —
   within-variant sd 0.004 Å. The near-constant gate RMSD of 0.18 Å across all
   41 variants is measuring the restraint, not the protein. This was the right
   setting for Arm 1 (it stops the backbone drifting so cavity volumes stay
   comparable) and the wrong one for Arm 2, whose entire question is whether the
   backbone moves.

**Fix before Arm 2 is used for any decision:** raise nstruct to ≥20, and either
drop the Cα constraints on the gate/latch or replace the RMSD metric with one
the constraints do not bound (e.g. gate–latch closure distance, or ΔΔG of
binding from independently relaxed bound/unbound ensembles). Arm 1 is unaffected
and its ranking stands on its own.

→ **Both fixes implemented in §5c.** `07_rosetta_ratchet.py` is superseded by
`21_rosetta_ratchet_v2.py`; the v1 results in `results/ratchet/` are retained
only as the record of the failed assay and must not be cited.

### 5c. Arm 2 v2 — the corrected readout assay

`scripts/21_rosetta_ratchet_v2.py`, submitted by `scripts/22_submit_ratchet_v2.sh`,
merged by `scripts/23_merge_ratchet_v2.py` → `results/ratchet_v2/`.

**Fix 1 — sampling.** nstruct 3 → **20**. One FastRelax replicate on the
474-residue ternary complex is ~50 min, so 20 serial replicates per variant
would be ~17 h. The work is chunked instead: **205 array tasks = 41 variants ×
5 chunks × 4 replicates = 820 relax runs**, ~3.5 h per task, on `cutlerlab`.
Task *i* → variant *i* mod 41, chunk *i* ÷ 41, replicates 4·chunk … 4·chunk+3.

**Fix 2 — stop the constraints from measuring themselves.** The global flag
`-relax:constrain_relax_to_start_coords` is **removed**. Coordinate constraints
are applied explicitly via `CoordinateConstraintGenerator` (Cα-only,
`set_sd(0.5)`, harmonic) to every residue **except** the gate (85–89) and latch
(115–117) of chain A, selected with `NotResidueSelector(ResidueIndexSelector(…))`.
`ref2015_cst` already carries `coordinate_constraint` at weight 1.0, verified.

The rest of PYR1 and all of HAB1 stay pinned — so variants remain comparable and
the pose cannot drift wholesale — while the gate and latch backbone is free to
respond to the mutation and be measured. This is the controlled form of the
experiment: everything fixed except the thing being asked about.

**New metrics the constraints do not bound** (the primary readout):

| metric | what it asks |
|---|---|
| `gate_bb_rmsd`, `latch_bb_rmsd` | backbone (N, Cα, C, O) displacement — now free |
| `gate_latch_min_dist` | closest gate↔latch heavy-atom approach; the gate closing *onto* the latch **is** the ratchet, so a variant that opens this has lost the switch |
| `gate_latch_com_dist` | centroid version of the same, less noise-prone |
| `lig_gate_min_dist` | does ABA still contact the loop it is meant to hold shut? |

Retained for continuity: `dG_separated`, `dSASA`, `sc_value`, `w385_aba_min_dist`,
`total_score`, and the v1 all-heavy-atom `gate_rmsd`/`latch_rmsd`.

`23_merge_ratchet_v2.py` re-runs the same variance decomposition that condemned
v1, so the assay is re-validated before its numbers are used. **If the ratios are
still < 2 after this, the problem is neither sampling nor constraints and the
metric itself must change.**

### Arm 1 — cavity volume after full repacking (informative, ratio 8.6)

23 of 40 variants open the cavity by ≥25 Å³ *after* FastRelax lets the
neighbourhood collapse. Top of the ranking:

| variant | cavity Å³ | ΔCav | ΔREU |
|---|---|---|---|
| K59A/F108A/R79A/E94A | 382.4 | **+208.0** | +11.74 |
| F108A/R79A/E94A | 313.7 | +139.3 | +10.11 |
| **K59A/F108A** | **304.6** | **+130.2** | **+2.98** |
| F108A/R79A | 262.9 | +88.5 | +7.59 |
| R79A/E94A | 255.4 | +81.0 | +10.29 |
| K59G | 253.9 | +79.5 | +4.95 |
| F108G | 249.2 | +74.8 | +3.67 |
| K59A | 242.4 | +68.0 | **−2.40** |
| F108A | 240.3 | +65.9 | +2.38 |

Three findings worth keeping:

- **The cavities survive repacking.** The rigid scan (§3) predicted the quad at
  375.5 Å³; the relaxed quad gives 382.4 Å³. The wall does not collapse back in.
  This was the central risk of the whole approach and it did not materialise.
- **K59A/F108A is the standout.** +130 Å³ for **+2.98 REU** — an order of
  magnitude cheaper per Å³ than the quad (+208 Å³ for +11.74 REU), and K59A
  alone is *stabilising* (−2.40 REU). Both positions are individually validated
  in real receptors (F108A in PYR1^MANDI/4WVO, §4b).
- **R79 shows strict epistasis with E94.** R79A alone gives **−2.8 Å³** — no
  gain at all — but R79A/E94A gives **+81.0 Å³**. Deleting the arginine without
  its salt-bridge partner lets E94 relax into the vacated space and reseal the
  wall. The buttress must be broken from *both* sides or not at all. This is
  invisible to the rigid scan, which reported R79A/E94A as simply additive, and
  it is a direct argument for pairing these positions in any library rather than
  randomising them independently.

The bulky substitutions behave as expected: E141W (−30.5 Å³), E94L (−28.5),
F108Y (−12.9) all *shrink* the pocket, confirming the volume metric responds to
side-chain size in the right direction.

---

## 7b. RESUME HERE — state as of 2026-08-07
> ⚠️ **SUPERSEDED by §21.** Kept as the historical record of the pilot handoff.

### Jobs — all complete

| job | what | outcome |
|---|---|---|
| **27260956** | pilot array 0–40, Arms 1+2, nstruct=3 | ✅ 41/41 both arms, ~2.5 h/task |
| **27261151** | Arm 3 compute stage | ✅ 147 homolog cavities, 32 min |

**Re-running anything:** every step is idempotent. `06`/`07` take `--index`
(0–40, panel order in `variants.py`) and overwrite their own JSON. Arm 3's
fetch stage is cached — it skips files already in `results/homolog_cavities/raw`.

### Open threads, in priority order

1. **Re-run Arm 2 properly** (§5b). It is currently uninformative: nstruct ≥20
   and a metric the Cα constraints do not bound. Until then **no variant has
   been shown safe for the readout** — the Arm-1 ranking is a cavity ranking
   only. 41 tasks × ~20 relax replicates on the 474-residue ternary complex is
   the main compute ask; ~25 min per replicate on `cutlerlab`.
2. **The §4d conflict is the biggest open scientific question.** The bottleneck of the
   merged channel is the gate/latch itself, and the expanded cavity has no tail
   exit. This makes tier-3 C40 ligands a diagnostic rather than a target, and
   promotes the **deepening/insertion arm** (extend away from the gate/latch, at
   the distal β-hairpins) to the primary route for anything larger than tier 2.
   Nothing has been built for that arm yet.
3. **Calibrate before trusting any ligand score.** Anchor is PDB **4WVO**
   (PYR1^MANDI = K59R/V81I/F108A/F159L). The pipeline must rank
   PYR1^MANDI + mandipropamid above WT + mandipropamid first.
4. **No open-state test exists.** Both arms model the closed state only; nothing
   confirms an expanded variant can still *open*. Needs MD (GROMACS 2024.3 or
   Amber 22, both available) apo and holo on the top candidates.
5. **Arm 3's signature output** (`position_signature.json`) has not been read —
   it is the input to any grafting design.

### Things established that should not be re-derived

- Distance-to-ligand is the wrong metric for "lines the pocket"; use
  `04_cavity_lining.py` (§2). R79 is a cavity-rim residue with clear line of
  sight to ABA, not a remote second-shell residue.
- R79 has never been in a combinatorial PYR1 pocket library; Beltrán's own
  geometric filter selected it and manual curation removed it (§4b).
- F108A is already validated in PYR1^MANDI (4WVO) — reducing F108 preserves the
  ratchet.
- `/scratch` is node-local; SLURM work must live on `/bigdata`. Login node has a
  1 GB cap that kills foldseek and FastRelax.

---

## 8. Status
> ⚠️ **SUPERSEDED — this section and §7b describe the state on 2026-08-07, before
> the arms read out. See §21 for the current status.** Both are kept because the
> §5b "Arm 2 is blind" line below is the claim that §14a/§14b later corrected, and
> the record of that correction matters.

| Step | Status |
|---|---|
| Inputs prepared, ligand params (charge −1) verified | ✅ |
| Pocket contacts + buttress network (§2) | ✅ |
| Cavity-lining + line-of-sight analysis (§2) | ✅ |
| Rigid-backbone truncation scan (§3) | ✅ |
| foldseek CATH50 search + conservation (§4) | ✅ |
| Prior-art audit: Tian, Mosquna, Park×2, Beltrán (§4b) | ✅ |
| Carotenoid/START survey + cavity shape + ligand ladder (§4c) | ✅ |
| Exit-vector and bottleneck analysis (§4d) | ✅ |
| Arm 1 — cavity scan, 41 variants | ✅ 41/41, discriminating (ratio 8.6) |
| Arm 2 — ratchet scan, 41 variants | ⚠️ 41/41 ran, but **blind** — must re-run (§5b) |
| Arm 3 — homolog cavity signature | ✅ 147 full-atom homolog cavities |
| Aggregation of pilot results (`20_aggregate_pilot.py`, §5b) | ✅ |
| Insertion/grafting arm | ❌ not started (promoted to primary path by §4d) |
| Open-state / MD validation | ❌ not started |

### Next steps

1. **Re-run Arm 2** with nstruct ≥20 and an unbounded readout metric (§5b).
   Nothing about readout safety can be concluded until this is done.
2. **Act on §4d.** The gate/latch *is* the bottleneck and the expanded cavity is
   sealed, so build the deepening/insertion arm — extend the cavity away from
   the readout at the distal β-hairpins, using +2 register-preserving units on
   paired strands and hairpin grafts from large-cavity homologs (Arm 3 output).
3. **Calibrate on 4WVO** before trusting any affinity or docking number.
4. **Add the open-state test** — short MD, apo and holo, on the top candidates.
5. **Remember the screen.** With IGGYPOP and the yeast ternary readout, the
   in-silico job is **enrichment, not prediction** — shift hit rate from ~0.1%
   to ~5% in a 10³–10⁴ library, with counterselection against constitutive
   HAB1 binding.

---

## 9. Six design objections, tested (2026-08-07)

Raised by the user against the §5 arm design. Each was checked against data
rather than argued; three overturned or qualified a claim in this README.

### 9a. ⚠️ CORRECTION — HAB1 W385 *does* sense the ligand, through a water

§1 and §5 stated HAB1 "never touches the ligand" and reads only the closed-state
surface. **That is wrong as written.** Measured on 3QN1 (`data/3QN1.cif`):

| contact | distance |
|---|---|
| W385 NE1 → **water A197** | 3.04 Å |
| water A197 → ABA **O10** | 2.72 Å |
| W385 → ABA, direct heavy-atom | 4.67 Å |

A single ordered water bridges the HAB1 tryptophan to the ABA carboxylate — a
canonical water-mediated hydrogen bond. In other published receptor/ligand
complexes the ligand itself occupies that water's position and contacts W385
directly. So the ratchet framing survives (HAB1 contributes no *pocket*
residues, §2), but the stronger claim that HAB1 is blind to ligand chemistry
does not: **the ligand's polar head is read by HAB1 at one remove.**

Design consequences, both live:
- A larger ligand that displaces water A197 changes an interface hydrogen bond,
  not just pocket packing. `w385_aba_min_dist` cannot see this — Rosetta relax
  here carries no explicit waters. **This is an unmodelled term in Arm 2.**
- Whether W385A abolishes binding through lost ligand sensing or through the
  ~100 Å³ interface cavity left by Trp→Ala is **not resolved by our data**, and
  the two make different predictions for oversized ligands.

### 9b. R79/E94 epistasis is symmetric, and K59 is independent

Neither single mutation opens anything; only the pair does.

| variant | ΔCav Å³ | additive expectation | **epistasis** |
|---|---|---|---|
| R79A | −2.8 | — | — |
| E94A | +7.8 | — | — |
| **R79A/E94A** | **+81.0** | +5.0 | **+76.0** |
| F108A/R79A | +88.5 | +63.1 | +25.4 |
| F108A/R79A/E94A | +139.2 | +70.8 | +68.4 |
| quad | +208.0 | +138.8 | +69.2 |
| **K59A/F108A** | **+130.2** | +133.8 | **−3.6** |

The effect is **not directional** — it is symmetric synergy. Whichever partner
is deleted first, the other relaxes in and reseals the wall. This splits the six
positions into two mechanistically distinct groups:

- **R79 / E94 / F108** — a coupled buttress. Must be mutated together; single-
  position libraries at these sites will read as dead.
- **K59** — energetically independent (epistasis −3.6 Å³, i.e. additive within
  noise). Can be varied on its own.

### 9c. ⚠️ Arm 2's free set is ~10× too tight — measured against the real stroke

`scripts/24_state_comparison.py`. 3K3K chain A is **apo** (open protomer), chain
B carries ABA; 3QN1 chain A is the closed ternary state. Superposed on core Cα
excluding gate/latch and their ±3 flanks:

| region | apo→ternary (**the switch stroke**) | ABA-bound→ternary |
|---|---|---|
| gate 85–89 backbone | **5.23 Å** | 0.93 Å |
| latch 115–117 backbone | **5.34 Å** | 0.96 Å |
| gate ±3 flanks (Cα) | **3.98 Å** | 0.74 Å |
| latch ±3 flanks (Cα) | **3.52 Å** | 0.63 Å |
| all residues (Cα) | 2.18 Å | 1.05 Å |

Two conclusions:

1. **ABA binding does the closing; HAB1 binding adds almost nothing** (0.93 Å).
   This is direct experimental support for the ratchet model — and it means Arm 2,
   which starts from the closed ternary, is measuring *stability of the closed
   state*, not the switch.
2. **The flanks move 3.5–4.0 Å, and Arm 2 pins them.** Releasing only 85–89 and
   115–117 lets the loops bulge but not hinge; the v2 smoke produced 0.36/0.64 Å
   against a real stroke of 5.2/5.3 Å. For the stability question this is
   defensible; for "is the mechanism preserved" it is not.

**v3 free set should be 82–92 and 112–120** (gate/latch ±3), which is exactly the
span that moves experimentally. Not yet run.

### 9d. There is a real interface change beyond the gate/latch

Comparing the ABA-bound protomer to the ternary complex — the comparison that
already agrees to 1.05 Å overall — the HAB1-contacting surface deviates *more*
than the protein average:

| set | Cα RMSD |
|---|---|
| all residues | 1.05 Å |
| HAB1 interface, gate/latch excluded (12 res) | **1.71 Å** |
| largest movers | G150 **7.4 Å**, N151 4.7, E149 4.5, S152 2.8, E153 2.6 |

The excess is localised to the **149–156 loop**, which carries interface
residues P148, N151, D155, T156. This is a contiguous, HAB1-contacting element —
not scattered noise, which is what crystallographic variation would look like.
So there is a **third moving element** besides gate and latch, and no Arm 2
metric currently watches it. Caveat: two crystal forms, so packing cannot be
fully excluded; confirming needs more ternary structures.

### 9e. Why nstruct=20 for Arm 2 but 3 is fine for Arms 1 and 3

Required replicates scale with (noise ÷ effect)², so the two arms are not
comparable:

| | replicate sd | between-variant sd | SEM at n=3 | verdict |
|---|---|---|---|---|
| Arm 1 `cavity_A3` | 5.47 Å³ | 47.1 Å³ | 3.2 Å³ | SEM is ~2% of a 130 Å³ effect — **n=3 ample** |
| Arm 2 `dG_separated` | 2.78 REU | 2.04 REU | 1.6 REU | SEM is **78% of the entire signal range** |

For a two-sample comparison at 80% power, n ≈ 2(σ(z₀.₉₇₅+z₀.₈)/Δ)². With
σ = 2.78 REU: detecting Δ = 3 REU needs **n ≈ 14**; n = 20 resolves ~2.5 REU.
Arm 3 needs no replicates at all — it is a deterministic grid calculation on
fixed coordinates, with zero stochastic component.

### 9f. The homolog set is *not* PYR/PYL-only — and F108 is the fold-wide lever

Conservation here is **structural**, not a sequence MSA: foldseek against CATH50,
147 hits at alnTM ≥ 0.5, residue identity read off the structural alignment.

| property | value |
|---|---|
| sequence identity to PYR1 | median **13.4%** (range 6.4–58.4%) |
| >50% id (PYR/PYL-like) | **5** of 147 |
| ≤30% id (distant fold) | **138** of 147 |
| structural alnTM | 0.605–0.970 |
| domain length | 93–270 res (PYR1 ≈ 180) |
| cavity volume | 0–634 Å³ (PYR1 WT 174) |

So the comparison is genuinely fold-wide and the ABA-only concern does not
apply. Within the 138 distant homologs:

| position | small residue | large residue | interpretation |
|---|---|---|---|
| **108** | median **194 Å³** (n=35) | median **74 Å³** (n=103) | **2.6× — a real fold-wide size lever** |
| 59 | 126 Å³ (n=53) | 125 Å³ (n=85) | **no effect — K59 is PYR1-specific** |

And cavity volume tracks domain length (Spearman ρ = 0.381, p = 8×10⁻⁵):
93–140 res → 20 Å³; 140–170 → 135 Å³; 170–220 → 153 Å³; 220–300 → **185 Å³**.
Independent support for §4e: **big pockets in this fold come substantially from
more backbone, not just smaller side chains.**

Counterexample worth chasing: `2pcsA00` reaches **569.6 Å³ in only 152
residues** — a large cavity without extra length. Whatever it does differently
is the most direct template available.

### 9g. Proposed Arm 4 — graft the gate/latch onto a large-cavity scaffold

User's proposal, and the data supports it over enlarging PYR1 further. §4e showed
PYR1 is ~9 Å short for a C40 with both cavity ends backbone-bounded and side-chain
truncation saturated. Rather than fight that, start from a homolog that already
has the volume (top hits: 634, 570, 554, 455 Å³ vs PYR1's 174) and transplant the
gate and latch onto it, then filter as in Arms 1–2.

Prerequisites, none yet done: check that candidate scaffolds present β3–β4 and
β5–β6 loops in a compatible geometry; that the HAB1-contacting face (§9d — now
known to include the 149–156 loop, not just the loops) can be reconstructed; and
that the scaffold's cavity actually opens to where a gate would sit. Arm 3's
`position_signature.json` plus these cavity volumes are the input.

---

## 10. `2pcsA00` identified — it is a CoxG, a quinone-extracting SRPBCC protein

The largest-cavity compact homolog from Arm 3 (§9f) turns out to be a member of
a family whose biology is directly relevant to this project.

### What the structure is

| field | value |
|---|---|
| PDB | **2PCS**, 2.40 Å, X-ray, deposited 2007-03-30 |
| title | "Crystal structure of conserved protein from *Geobacillus kaustophilus*" |
| depositors | Bonanno, Gilmore, Bain, Wu, Romero, Smith, Wasserman, Sauder, Burley, Almo — **NYSGXRC** structural genomics |
| citation | **"To be published"** — never written up |
| UniProt | **Q5QL47**, gene **GKP20** (plasmid-borne), "hypothetical conserved protein" |
| Pfam / InterPro | **COXG (PF06240)** / CO_DH_gsu (IPR010419), START-like superfamily (IPR023393) |
| SUPFAM / Gene3D | Bet v1-like (SSF55961) / **3.30.530.20** — same CATH node as PYR1 |
| ligand | **UNL** — unknown ligand, unidentified density |

So there is no literature on this protein *specifically*: it is an unpublished
structural-genomics deposition of a hypothetical protein. But its **family** was
characterised in 2025, and that explains the unknown ligand.

### The unknown ligand is a lipid tail

Measured from the deposited coordinates:

- **30 atoms, all modelled as carbon**, zero heteroatoms
- **18.59 Å** end to end, consecutive spacings 1.87–2.76 Å (a chain, not a ring system)
- mean B-factor **74.6** — poorly ordered, consistent with a partly mobile tail
- lining is essentially **all hydrophobic**: Y74, L20, L76, F136, L26, I63, Y109,
  A89, L125, Y59, L132, I30, L48, V55, V65, V78, I50, V87, I133, L36 … with a
  single polar residue (E90) at 4.3 Å

A 30-carbon unbranched hydrophobic chain in a wholly hydrophobic tunnel. The 2007
depositors had no way to assign it; in light of the family's function it reads as
a **co-purified isoprenoid/lipid tail**.

### What the family does (Wilson *et al.*, *Nat Chem Biol* 2025)

CoxG is the accessory protein of [MoCu]-CO dehydrogenase. The 2025 work on
*Mycobacterium smegmatis* Mo-CODH shows CoxG is a **membrane-bound
menaquinone-binding protein** that:

- has an **N-terminal SRPBCC domain** — seven-stranded antiparallel β-sheet plus
  two α-helices enclosing a large hydrophobic cavity (**521.4 Å³** by their
  measurement) — joined by a flexible linker to a **C-terminal membrane-anchoring
  helix**;
- **extracts menaquinone out of the membrane**, selectively binding
  **DH-MQ9 / MQ9** (~34% occupancy from *M. smegmatis* membranes; ubiquinone-8
  binds poorly, so this is real selectivity);
- delivers it to Mo-CODH, which reduces it with electrons from atmospheric CO,
  after which CoxG **returns the menaquinol to the membrane**;
- opens through a face **flanked by five lysines**, a positively charged rim that
  engages phospholipid head groups during extraction.

Crystal structure of the soluble domain: **PDB 8UDS**.

### Cross-validation of our cavity method

Measured with `lib_cavity.py` on the identical criterion used throughout:

| structure | our measurement | published |
|---|---|---|
| CoxG-NT, *M. smegmatis* (**8UDS**) | **597.0 Å³** | 521.4 Å³ |
| **2PCS** chain A | **571.5 Å³** | — |
| PYR1 WT (3QN1 A) | 174 Å³ | — |
| PYR1 quad-truncated | 382 Å³ | — |

Agreement with an independently published volume to ~14% is a useful external
check on the grid criterion — the first this project has had.

### Why this matters for the design

1. **The fold demonstrably binds something far longer than a C40 carotenoid.**
   MQ9 carries a **C45 polyprenyl tail**. This is the strongest evidence yet that
   the helix-grip fold is not intrinsically size-limited.
2. **But it does not fully enclose it.** The docking in that paper places the
   redox-active head group *near the opening*, with the tail in the interior.
   That is precisely the partially-enclosed, tail-out arrangement flagged as a
   risk in §4c/§4d — and here it is the *functional* arrangement, because the
   head must be presented to Mo-CODH.
3. **Loading is mechanistic, not diffusive.** CoxG gets its ligand by membrane
   extraction using a C-terminal TM helix and a lysine-rimmed opening. PYR1 has
   neither, and a soluble cytoplasmic biosensor cannot easily acquire one.
4. **The direct conflict with our readout is now explicit.** CoxG's pocket works
   because it is *open*. PYR1's readout works because the gate and latch *close*.
   Any PYR1 variant that accommodates a C40–C45 lipid the way CoxG does would
   present an open face where HAB1 expects a closed one.
5. **A compact scaffold can still be large.** 2PCS reaches ~570 Å³ in **152
   residues**, fewer than PYR1's ~180. So the §9f length correlation (ρ = 0.381)
   is a tendency, not a requirement — cavity volume here comes from how the
   β-sheet is splayed, not from extra chain. That makes CoxG-like proteins the
   better grafting scaffolds for Arm 4 (§9g).

**Caveat on 2PCS specifically:** the COXG assignment is by Pfam, and *GKP20* is
plasmid-borne in a thermophilic *Bacillus* relative with no verified Mo-CODH
context. Treat it as a CoxG-*family* protein of unproven function. **8UDS is the
better template** — it is a characterised CoxG with a published mechanism.

---

## 11. What this project is actually for (framing, corrected 2026-08-07)

Recorded because two earlier framings in this README were wrong, and the wrong
one drives the wrong experiments.

### 11a. The objective

**Goal: expand the ligand pool the sensor system can detect.** PYR1 has a
demonstrated size limit, so enlarging the pocket is the route — but the pocket is
not the product. **The product is a new backbone exposing additional library
positions**, from which a Tian-style oligo-pool library screen can be run. The
user already has Tian's libraries (built by combining oligo pools); the natural
extension adds oligos that both enlarge the pocket and introduce substitutions at
positions those pools never reached.

**This is still a pilot.** It answers two questions only: is expanding the pocket
plausible, and do the computational experiments give evidence the closing
mechanism would survive in vivo. A library follows only if the answer is yes.

### 11b. What counts as a win — and a correction

Not a win: shrinking residues already in the Tian libraries and calling the
enlarged pocket a new substitution library. Those substitutions are **already
inside the existing oligo pools**; the result is the same positions relabelled,
with no new combinations and no new ligands.

A win: **new positions entering the pocket lining** —
1. opening the second lobe (R79/E94; recruits I48, V49, R50, V81, I62, S122),
2. insertions creating genuinely new pocket positions,
3. relaxing key residues so more residues line the pocket,
4. grafting a larger scaffold.

⚠️ **Correction to an over-correction.** An earlier draft of this section said
shrinking Tian-library residues "is not a win," full stop. That is wrong as
guidance. It is not a win *as a final deliverable*, but a maximally open pocket
is a legitimate and necessary **design starting point** — the shrink-to-fit
workflow (§13) needs a large cavity to design into. Hence §12b.

**Mosquna's 39-residue set is a weaker constraint than Tian's 18.** The user does
not have those libraries and the set extends well beyond the pocket. Treat Tian's
randomised positions as the binding prior-art constraint and Mosquna as a
"less intriguing if covered" flag.

### 11c. Prior-art status of the six wall positions

| position | Mosquna 2011 SSM (39) | Tian randomised (18) | recovered in characterised hits (419 clones) |
|---|---|---|---|
| K59 | ✅ | ✅ | 262 clones (N, M, Q, D, A, L, R, T) |
| **R79** | ❌ **absent** | ❌ **never randomised** | **never mutated in any hit** |
| E94 | ✅ | ✅ | 189 clones (G, D, A, S, T, L, N) |
| F108 | ✅ | ✅ | 42 clones (W×39, Q, L, H) |
| Y120 | ✅ | ✅ | 205 clones (L, M, A, S, T, F, G…) |
| E141 | ✅ | ✅ | 130 clones (Q, G, F, K, W, T, I…) |

**R79 is the only position the field has never sampled**, and §9b showed it is
inert alone and only functions with E94 — a plausible mechanistic reason
single-position scanning missed it.

⚠️ **Competing explanation, on the record.** "Never sampled" is a statement about
experimental design, not fitness. R79 is **invariant across PYR/PYLs** while
variable across the superfamily — the signature of family-specific functional
constraint. Our own numbers show a cost (R79A +7.55 REU, R79A/E94A +10.29), and
§12a shows R79 single mutants actively *seal* the second lobe. So R79 is
unexplored **and** plausibly constrained; the pilot must distinguish "overlooked"
from "avoided for a reason." The open-state test (§13c) is how.

---

## 12. Where the existing landscape actually fails, measured

### 12a. Ligand size and shape (`27_tian_ligand_size.py`, `28_tian_ligand_shape.py`)

RDKit (ETKDG + MMFF) over Tian's 3366-compound screening library (sd01), 194
hits. **MW is not used for any conclusion** — their largest hit by mass is Tetrac
(747.8 Da) but four iodines supply ~508 Da of that. Heavy-atom count and longest
dimension decide pocket fit.

**Size arrives as bulk, not length:**

| heavy atoms | n | median length | length/heavy | ROD / DISC / SPHERE |
|---|---|---|---|---|
| <20 | 1277 | 7.5 Å | 0.55 | 66% / 29% / 3% |
| 20–25 | 699 | 10.3 Å | 0.47 | 65% / 30% / 3% |
| 25–30 | 552 | 12.0 Å | 0.45 | 65% / 26% / 7% |
| 30–40 | 611 | 13.4 Å | 0.40 | 64% / 28% / 6% |
| ≥40 | 220 | 15.8 Å | **0.33** | 49% / 35% / 15% |

**Hit rate by size × shape — the table that decides tunnel vs chamber:**

| heavy atoms | ROD | DISC | SPHERE |
|---|---|---|---|
| 20–25 | 9.8% | 8.9% | 4.0% |
| 25–30 | 4.7% | 2.0% | 4.7% |
| 30–40 | 1.5% | 2.3% | 2.7% |
| **≥40** | **0/108** | **0/77** | **0/35** |

Two separable limits:
- **A volume wall at ~40 heavy atoms, independent of shape** — 0 hits from 220
  compounds, rods and spheres alike.
- **A length wall at ~15 Å.** Within every size bin, hit and miss *median*
  lengths are nearly identical, but the tails diverge: no hit anywhere exceeds
  **15.3 Å**, while misses reach **21.0 Å** at 30–40 heavy atoms.

⚠️ Earlier drafts framed the ceiling as length alone. It is both.

**The carotenoids were screened and failed.** Longest failures ≥30 heavy atoms:
**lutein (29.1 Å)**, Trypan Blue (28.8), **β-carotene (27.9)**, Evans Blue (27.8),
Elbasvir (27.2), Velpatasvir (26.0), Ledipasvir (25.8). So the C40 target is
empirically unreachable from the existing landscape — which supports building a
new backbone and simultaneously warns that C40 is the wrong *first* target.

**Target zone: 30–45 heavy atoms, 15–20 Å** — where hit rate collapses from ~5%
to 0. The merged-lobe cavity (span 18.6 Å, ~390 Å³ relaxed) sits in that band.
Full C40 needs backbone work (§4e).

⚠️ Do not over-read the span↔hit-length link: it is a correlation within one
receptor's libraries, and "hit" folds in solubility, permeability and expression.
It bounds the existing landscape well; it predicts a new scaffold's reach loosely.

### 12b. Second-lobe access and the maximal pocket (`29_second_lobe_access.py`)

ΔCav is the wrong ranking metric: a variant can hollow the existing chamber,
gain 100 Å³, and recruit no new lining residues. What matters is second-lobe
access. WT satellite lobe = 39.5 Å³, 2.95 Å from the K59/R79/E94/F108 cluster.

**Robust result — R79 single mutants uniquely SEAL the second lobe.** R79A/G/L/S/V
all give `lobe2 = 0.0, merged = NO`. Every other variant retains or increases it.
Removing the arginine alone lets E94 collapse in and close the region — independent
confirmation of the §9b epistasis from a different measurement.

**F108 is the gatekeeper, not R79/E94.** F108 mutants top the ranking
(49–55 Å³ lobe2); the quad reaches 59.9; but R79A/E94A gives only 41.9 vs WT 38.1
(**+3.8**). The pair opens +81 Å³ of *total* cavity while adding little to the
*satellite* region — this qualifies how §9b was described.

⚠️ **Metric caveat.** Relaxed WT already registers 38.1 Å³ in the satellite region
and reads as merged, while rigid WT on the same mask gives 0.0. Either FastRelax
alone breaches the wall, or the 1 Å mask dilation bleeds into the main cavity.
Per-variant absolute numbers are **provisional** until this is separated; the R79
sealing flip is robust because it is a qualitative reversal.

**Maximal pocket (rigid truncation, gate/latch never touched):**

| truncation set | n | cavity | lobe2 |
|---|---|---|---|
| Tian pocket set (sealing kept) | 14 | 545.6 Å³ | 54.5 |
| **Tian set + R79** | 15 | **602.6 Å³** | 61.4 |
| Tian set + R79/E94 + lobe2 lining | 19 | 608.0 Å³ | 57.6 |

**Truncate every Tian pocket position → 545.6 Å³. Add R79 → +57 Å³.** That is
volume the existing oligo pools cannot reach, quantified. 608 Å³ is the design
ceiling for shrink-to-fit, and it is in CoxG territory (571–597 Å³, §10).

---

## 13. Proposed workflow

### 13a. Screening strategy — use Tian's own failures as the probe set

**Do not select on ABA.** ABA already fits the small pocket; second-lobe mutations
can only hurt it there, so an ABA screen selects *against* the expansion.

Instead, filter Tian's library to the band the expanded cavity targets —
**28–50 heavy atoms, 15.3–22 Å**:

> **229 compounds. Zero hits. 217 rods, 12 discs.**

Screened against six libraries covering 18 randomised positions, and none worked.
They are commercial, already size/shape-characterised in
`results/28_tian_ligand_shape.csv`, and **the negative control is a published
screen** — a hit on any of them is unambiguous and directly comparable. Widening
to >28 heavy and >15.3 Å gives 289.

Easiest first targets (just past the 15.3 Å ceiling, not far beyond):
cholecalciferol (28 heavy, 15.3 Å), deoxycholic acid (15.5), estradiol
dipropionate (15.5), terazosin (15.6), sitagliptin (16.1), etravirine (16.8),
glasdegib (18.3).

Pooling: user will follow Tian's protocol; not a constraint on design.

### 13b. The failure mode counterselection does not catch

5-FOA on the Y2H removes constitutives cheaply. It does **nothing** for receptors
that never close — those are *silent*, indistinguishable from "none of these
chemicals binds." If much of an expanded library is conformationally dead, the
effective library collapses invisibly.

> **Hard requirement: every backbone entering a screen must retain a measurable
> response to at least one known ligand.** Without a per-backbone positive
> control, a negative screen is uninterpretable. ABA works as a diagnostic even
> though it is the wrong selection ligand; better is one of Tian's larger hits
> (28–36 heavy atoms), which tests closure nearer the target size.

**Division of labour:** predicting novel binding is unreliable; predicting whether
a mutant retains a **conformational switch** is tractable. Computation guarantees
switchability and pocket geometry; the chemical screen finds ligands.

### 13c. Arm 2b — the open state (highest priority, not yet built)

Both existing arms start from the closed ternary complex and ask whether it stays
put. Neither tests the **switch**: apo/open → ligand → closed → HAB1 binds.
Enlarging a pocket generally destabilises the closed state relative to open, so
the failure mode we are least equipped to see is the one expansion most plausibly
causes.

**3K3K chain A is an experimentally determined apo, open PYR1** (§9c) and is now
in `data/`. The same 41 variants can be relaxed on the open state and scored for
the closed−open gap. No new structural modelling required. **This should rank
above the grafting arm** — cheaper, and it gates everything downstream.

### 13d. Setting filter cutoffs empirically

Do not pick a threshold on predicted affinity. **Calibrate on Tian's 419
characterised clones** — experimentally confirmed functional receptors, 132 of
which broke the E94–R79 bridge and still worked. Run any filter over that set,
see where true positives fall, and set the threshold to retain ~90–95% of them.
That converts an arbitrary cutoff into a measured false-negative rate on real
data, with 4WVO/PYR1^MANDI as an independent anchor.

This is also the honest test of whether our metrics mean anything: if the
switchability filter rejects many clones that demonstrably work in yeast, the
filter is wrong and we learn it before committing to a library.

### 13e. Shrink-to-fit design — use LigandMPNN as a LIBRARY designer

Build a maximally open cavity (§12b), then redesign the lining to complement a
specific ligand. Correct shape for the problem, with two provisos.

1. **Pose first, and it is the hard part.** LigandMPNN needs a fixed backbone and
   a ligand pose; docking a novel large ligand into a de-designed cavity is where
   docking is least reliable, and a wrong pose yields confident nonsense. Anchor
   on ABA's carboxylate geometry where the ligand has a comparable polar head
   (PYR1 is head-**in** — K59 plus the W385 water, §9a), generate several poses,
   design against each, and keep only designs that agree.
2. **Design distributions, not sequences.** LigandMPNN emits individual sequences
   with many simultaneous substitutions; an 8-mutation design may be unreachable
   in an oligo pool. Instead constrain it to a defined position set, sample
   deeply, and take **per-position amino-acid distributions** — that output *is*
   an oligo-pool specification, and per-position entropy shows which positions do
   work versus which are tolerant.

Every design must pass the ratchet test (Arm 2 v3): a sequence designer optimises
binding and is indifferent to the mechanism.

**Calibrate before trusting any affinity number** — PDB 4WVO (PYR1^MANDI =
K59R/V81I/F108A/F159L) with mandipropamid. If the pipeline cannot rank that above
WT + mandipropamid, its numbers mean nothing on novel ligands. Boltz skills are
available in-session for structure/binding and small-molecule screening.

---

## 14. Pilot readout — all arms complete (jobs finished 2026-08-08, analysed 2026-08-10)

Jobs 27275686 (Arm 2 v2), 27286483 (Arm 2 v3), 27286484 (Arm 2b) all COMPLETED:
820/820, 820/820, 492/492 units, zero failures.

### 14a. ⚠️ METHODOLOGICAL CORRECTION — the discrimination test was wrong

Scripts 20, 23 and 33 judged an arm by

    sd(per-variant means) / mean(within-variant replicate sd)  >= 2

**That statistic is wrong.** Replicate sd is the spread of individual relax runs;
what limits separation of variants is the noise on each variant's *mean*, which
is sd/√n — smaller by 4.5× at n=20. The correct test is the one-way ANOVA ratio
`F = s²_between / (s²_within/n)`, with the true effect recovered as
`s²_true = s²_between − s²_within/n`. `scripts/35_variance_decomposition.py`
implements this and re-reports every arm.

| arm / metric | n | within sd | SEM | between sd | **TRUE sd** | F | p |
|---|---|---|---|---|---|---|---|
| Arm 1 `cavity_A3` | 3 | 6.92 | 4.00 | 47.11 | **46.9 Å³** | 139.0 | 7e-61 |
| Arm 1 `score_ref2015` | 3 | 1.71 | 0.99 | 4.24 | 4.13 REU | 18.5 | 3e-27 |
| Arm 2 v2 `dG_separated` | 20 | 2.57 | 0.57 | 0.76 | **0.50 REU** | 1.77 | 3e-03 |
| Arm 2 v3 `dG_separated` | 20 | 2.91 | 0.65 | 0.78 | **0.43 REU** | 1.44 | 4e-02 |
| Arm 2 v3 `gate_bb_rmsd` | 20 | 0.026 | 0.006 | 0.008 | 0.005 Å | 1.89 | 9e-04 |
| Arm 2 v3 `latch_bb_rmsd` | 20 | 0.043 | 0.010 | 0.029 | 0.027 Å | 8.88 | 2e-41 |
| **Arm 2b `dE_close`** | 12 | 2.19 | 0.63 | 1.90 | **1.79 REU** | **9.02** | 6e-37 |
| Arm 2b `score_closed` | 12 | 1.91 | 0.55 | 3.97 | 3.93 REU | 51.9 | 4e-143 |
| Arm 2b `score_open` | 12 | 1.14 | 0.33 | 3.53 | 3.51 REU | 114.6 | 9e-210 |

**Arm 2b was mislabelled "NOT DISCRIMINATING" by the old rule (ratio 0.87). Its
F is 9.0 and its true effect is 1.79 REU. It works.**

Judge by TRUE sd in physical units, not F: with n=20 a negligible effect becomes
significant while staying useless — `latch_bb_rmsd` has F=8.9 and a true spread
of **0.027 Å**.

### 14b. ⚠️ Arm 2 revised: the readout is INVARIANT, not "blind"

§5b called Arm 2 blind. With the corrected statistic at n=20 that is wrong, and
the right conclusion is more useful:

**Pocket mutation at these six positions does not measurably perturb the
closed-state PYR1–HAB1 interface.** True between-variant spread of
`dG_separated` is **0.43–0.50 REU against a −67 REU interface — under 1%.**
v2 and v3 agree (0.50 vs 0.43) from independent runs with *opposite* constraint
schemes, so this is robust.

This is a **positive result and exactly what the ratchet model predicts** (§1,
§9a): HAB1 contributes no pocket residues and reads the closed-state surface, so
changing the pocket lining leaves the interface alone.

Two consequences:
- **Retire Arm 2 as a ranking tool.** There is nothing to rank. Report it as an
  invariance result — pocket engineering does not break the interface *provided
  the closed state is reached*.
- **Freeing the flanks and the 149–156 loop (v3) changed nothing** versus v2. The
  §9c objection was correct in principle, and empirically it did not matter,
  because the closed state is insensitive either way.

### 14c. Arm 2b — the switch, and the finding that matters

Conformer integrity verified: the two states stayed 5.20–5.24 Å apart throughout
(starting stroke 5.21 Å), gate backbone RMSD 0.08–0.12 Å within each state. No
replicate interconverted.

**WT `dE_close` = −5.49 ± 3.02 REU** — WT intrinsically prefers closed, as it must.

**No variant is predicted dead.** All 41 retain a negative `dE_close`, from
−10.74 to −1.43. Even the quad (+208 Å³) sits at **−2.42**, retaining ~44% of
WT's closing margin.

> ### The key design result: cavity expansion does NOT predict switch cost.
>
> `ΔCavity` vs `ΔΔE_switch`: **r = +0.034, p = 0.84.** No relationship at all.

The intuition that a bigger pocket is necessarily harder to close is **not
supported**. The cost is set by *which* residues are changed, not by how much
volume is opened:

| variant | ΔCav Å³ | Δbeyond-wall Å³ | ΔΔE_switch | dE_close |
|---|---|---|---|---|
| **K59A/F108A** | **+130.2** | +35.2 | **+0.10** | −5.39 |
| **F108A/R79A** | **+88.5** | +32.9 | **+0.75** | −4.74 |
| F108A | +65.9 | +35.1 | +0.50 | −4.98 |
| R79A/E94A | +81.0 | +11.9 | +1.98 | −3.51 |
| quad | +208.0 | +83.9 | +3.07 | −2.42 |
| K59G | +79.5 | +4.6 | **+3.62** | −1.87 |
| Y120L | +19.6 | — | **+4.06** | −1.43 |

K59A/F108A buys 130 Å³ for essentially nothing; K59G buys 80 Å³ for 3.6 REU;
Y120L buys 20 Å³ for 4.1 REU. **Expansion can be nearly free if the right
residues are chosen** — which is precisely what a portfolio strategy needs.

**Best novel candidate: `F108A/R79A`** — +88.5 Å³, +32.9 Å³ beyond the wall,
ΔΔE_switch +0.75 REU (near-free), and it contains **R79, the position the field
has never sampled** (§11c).

**Constitutive risk flagged:** `K59V` (ΔΔE −5.25) and `K59D` (−4.71) favour the
closed state. Directly testable — Tian recovered **K59D in 19 characterised,
ligand-dependent clones**, so if those are not constitutive, this arm over-calls
in that direction. That is a calibration case, not yet a contradiction, since
those clones carry additional mutations.

### 14d. Second-lobe metric repaired — and F108 is the gatekeeper

`scripts/34_second_lobe_fixed.py` resolved the §12b ambiguity by diagnosis rather
than assumption:

| structure | components | ABA component | reaches satellite centroid? |
|---|---|---|---|
| WT **crystal** | 38 | 167.0 Å³ | **no** |
| WT **relaxed** | 46 | 178.0 Å³ | **YES** |

**Cause (a) holds: FastRelax alone breaches the wall.** The "sealed satellite
lobe" is a property of the crystal structure, not of the protein — side-chain
relaxation opens it even in WT. The metric now uses a dilation-free
**beyond-the-wall volume** (plane through the K59/R79/E94/F108 side-chain
centroid, normal from the ABA centroid toward the satellite centroid; sanity
check: 100% of the satellite and 0% of the ABA chamber lie beyond it).

**This settles the §12b open question: F108 is the gatekeeper, not R79/E94.**
Ranking on Δbeyond-wall: quad +83.9, F108A/R79A/E94A +47.6, E141G +38.9,
K59A/F108A +35.2, F108A +35.1, F108A/R79A +32.9 — every leading variant contains
F108 (or E141), while **R79A/E94A alone gives only +11.9**. ΔCavity and
Δbeyond-wall correlate at r = 0.750, so expansion does largely go past the wall,
but the *pair* is not what drives it.

The novelty argument for R79 therefore rests on **total volume** (+57 Å³ on top
of the full Tian truncation set, §12b) and on its **epistasis with E94**, not on
second-lobe access specifically. §9b's framing is corrected accordingly.

---

## 15. Novelty, done properly (`scripts/36_novelty_check.py`)

### 15a. Refined win condition (user, 2026-08-10)

Novelty is **not** "a position Tian never touched". A change at a randomised
position is still new if the *substitution* was outside their oligo pools, or was
allowed but never *recovered*, or occurs in a *combination* that never arose.
Their libraries were screened against a chemical library, so what has actually
been tested is a set of (sequence, chemistry) pairs. **What must be avoided is
re-testing combinations that were already tested** — not touching a position
that happened to appear in a library.

This corrects §11b, which treated "in the Tian libraries" as disqualifying and
was too strict.

### 15b. Four levels of novelty, from sd04 (design) and sd07/08/09 (419 clones)

| verdict | n | meaning |
|---|---|---|
| **NOVEL-POSITION** | 9 | contains a position never randomised — all R79 variants |
| NOVEL-SUBST | 0 | substitution outside the oligo pools |
| **NOVEL-UNRECOVERED** | 13 | pools allowed it, screen ran, it never came back |
| NOVEL-COMBO | 17 | substitutions individually recovered, never together |
| **RETREAD** | **1** | `E94G` only — already made and screened |

**39 of 40 variants are new at some level.**

### 15c. What the pools allowed versus what came back

| pos | WT | allowed by pools | recovered in hits | **allowed but NEVER recovered** |
|---|---|---|---|---|
| 59 | K | ADEFGHILMNQRSTVWY | ADLMNQRT | EFGHISVWY |
| **79** | **R** | **— never randomised —** | — | — |
| 94 | E | ADFGHIKLMNQRSTVWY | ADEGLNST | FHIKMQRVWY |
| **108** | **F** | ADEGHIKLMNQRSTVWY | **HLQW** | **ADEGIKMNRSTVY** |
| 120 | Y | ADEFGHIKLMNQRSTVW | ADFGIKLMNQSTV | EHRW |
| 141 | E | ADFGHIKLMNQRSTVWY | ADFGHIKLQSTW | MNRVY |

**F108 is the striking one: 13 of 17 allowed substitutions never came back, and
`F108A` is among them.** The pools could make it, 3366 compounds were screened,
and only H, L, Q, W were recovered — W dominating at 39/42. The *smaller*-F108
direction was reachable throughout and never selected. That sharpens the §4b note
("F108 was only ever made larger") into a precise, checkable claim.

"Allowed but never recovered" is the cheapest genuine novelty at a sampled
position: untested **in practice** rather than merely unexplored.

### 15d. Consequence for the lead candidate

`F108A/R79A` is now novel on **two independent axes** — R79 was never randomised,
and F108A was allowed but never recovered — and the pilot says it costs almost
nothing:

| | value |
|---|---|
| ΔCavity | **+88.5 Å³** |
| Δbeyond-wall | +32.9 Å³ |
| ΔΔE_switch | **+0.75 REU** (near-free; WT margin −5.49 → −4.74) |
| interface | invariant (§14b) |
| novelty | NOVEL-POSITION **and** NOVEL-UNRECOVERED |

`F108A` alone is the cheaper single: +65.9 Å³, +35.1 Å³ beyond the wall,
ΔΔE +0.50, and NOVEL-UNRECOVERED on its own.

---

## 16. Satellite lobe across all deposited chains — I was wrong twice

`scripts/37_satellite_across_structures.py`. Every independently determined PYR1
chain, superposed onto 3QN1 A and analysed on the identical grid criterion:

| structure | state | ABA cavity | beyond wall | satellite void | **connected?** |
|---|---|---|---|---|---|
| 3QN1 A | ABA + HAB1, closed ternary | 167.0 | 0.0 | 39.5 | **no** |
| 3K3K A | homodimer, APO protomer (open) | 206.4 | 0.0 | 39.1 | **no** |
| 3K3K B | homodimer, ABA-bound protomer | 182.5 | 0.0 | 0.0 | **no** |
| 3K90 A | ABA-bound | 210.0 | **42.8** | 0.0 | **no** |
| 3K90 B | apo (ACY/GOL) | 201.9 | 0.0 | 0.0 | **no** |
| 3K90 C | apo (GOL) | 179.4 | 0.0 | 0.0 | **no** |
| 3K90 D | ABA-bound | 171.0 | 0.0 | 0.0 | **no** |

**Not one of seven deposited chains has the satellite lobe connected to the ABA
site** — across three crystal structures, apo and holo, open and closed. The
discrete satellite void exists in only 2 of 7; in the other five there is no
enclosed void there at all.

### ⚠️ This reverses §14d

§14d concluded "the seal is a crystal artefact; FastRelax reveals the true
arrangement." **The evidence now points the other way.** Seven independent
crystallographic observations agree the lobe is sealed; one relax protocol says
otherwise. The parsimonious reading is that **the FastRelax merge is a Rosetta
rotamer artefact**, not a discovery.

Balanced statement of what is and is not known:
- Crystal structures at 2.0–2.5 Å may not resolve narrow channels, and the 0.88
  buriedness cut is strict — so "sealed" is a statement about this criterion.
- Conversely, crystal packing can restrain side chains, so a rotamer Rosetta
  finds may be genuinely accessible in solution.
- **3K90 chain A shows 42.8 Å³ of ABA-cavity volume beyond the wall plane**, so
  beyond-wall volume does occur crystallographically — just never continuity to
  the satellite centroid.

**Practical consequence.** Relaxed WT has a beyond-wall volume of 54.4 Å³ where
crystal 3QN1 A has 0.0, so the Arm-1 WT baseline is probably inflated. Variant
rankings should survive (one protocol throughout, same bias in every row), but
**absolute claims that a variant "merges the lobes" must be dropped.** Settling
this properly needs MD, not more static modelling.

---

## 17. Literature: three insights and one failure mode we cannot see

### 17a. Dorosh et al. 2013 (PLoS Comput Biol 9:e1003114) — MD of PYR1 activation

**Independent confirmation of §9d.** They identify loop **Lβ7α5 (P148–D155)** and
report that "the flexibility of loop Lβ7α5, which is only indirectly involved in
gate closure, is affected even more strongly than that of the gate upon
ABA-binding." That is the same element found in §9d from 3K3K-vs-3QN1 geometry
(largest movers G150 7.4 Å, N151 4.7, E149 4.5, S152 2.8, E153 2.6). A structural
comparison and an MD study converge on the same third moving element — freeing it
in Arm 2 v3 was correctly motivated even though it did not change discrimination.

**The W385 water network is bigger than §9a said.** They report **R116 and L117
(latch) form water-mediated hydrogen bonds with HAB1 W385**, which itself is
"inserted between gate and latch". So W385 is water-bridged to *both* the ligand
(§9a, water A197) and the latch. It is a hydration-network hub, not a single
contact — and our models carry no explicit waters at all.

**"Constitutive" is relative to a non-zero baseline.** They found apo-PYR1 makes
substantial PP2C interaction — up to **80% inhibition of ABI2 at 4:1
receptor:phosphatase** — validated experimentally. WT already has ligand-
independent activity, consistent with our WT `dE_close` = −5.49 REU (closed
favoured even apo). Constitutivity is therefore a matter of degree, which is
exactly why the cutoff must be calibrated (§13d) rather than set at zero.

**A mechanism our model omits entirely: dimerisation.** PYR1 is a homodimer;
H60 governs the oligomeric state (H60P gives mixed monomer/dimer with weak basal
activity), and dimerisation competes with PP2C binding. Arm 2b models a monomer
in two conformations and cannot see the dimer↔monomer equilibrium.

> **This plausibly explains the K59 over-call (§14c).** H60 and F61 are
> dimerisation residues (Beltrán excluded both on that basis), and **K59 sits
> immediately adjacent to them**. K59 substitutions may act on dimerisation
> rather than on closure — invisible to a monomer model, and a reason to distrust
> Arm 2b specifically at K59.

### 17b. Melcher 2010 / Peterson 2010 — not training data, but a third failure mode

Both concern pyrabactin agonist/antagonist selectivity. Their mutations
(PYL1 I137V, A190V, V193I; PYL2 A93F) alter **ligand-mode selectivity**, not
constitutivity or viability, so they are **not usable labels** for the §13d
calibration. Mosquna remains the only properly labelled set.

What they do supply is a failure mode neither our arms nor our cavity metrics can
detect: **non-productive binding.** In PYL2, pyrabactin binds in the *open*
conformation without stabilising gate closure — the pyridine sits 4.8–7.0 Å from
the gate closed, 11–13 Å open. **A ligand can occupy the pocket and fail to
trigger the switch.**

For an enlarged pocket this is a first-order risk: more room means more ways to
bind unproductively. Three consequences:
- Cavity volume, span and beyond-wall volume are all **necessary but not
  sufficient** — none of them can distinguish productive from non-productive.
- Arm 2b is apo and therefore blind to it by construction.
- Detecting it requires **ligand-bound** modelling: does the ligand stabilise the
  closed state relative to open, not merely fit.

**PYR1 I110 is the selectivity determinant** (≡ PYL1 I137, PYL2 V114) — I137V
converts PYL1 to pyrabactin-inhibited, V114I converts PYL2 to pyrabactin-
activated. Note that I110 **is in Tian's randomised set** and **lines the neck
region** identified in script 19. It is the position most likely to flip a novel
ligand between agonist and antagonist, and it deserves explicit attention in any
library built here.

Also noted: **PYL6 has high ligand-independent basal PP2C binding** — a natural
constitutive control if one is needed at the bench.

---

## 18. Dimerisation: is it worth keeping? (measured on 3K3K)

PYR1 dimer interface = chain A residues within 5 Å of chain B:
**M1, H60, F61, I62, K63, I84, S85, G86, L87, P88, A89, S152, D155, T156, M158,
F159, T162, V163, L166** (19 residues).

### The dimer is a competitive OFF state, not an incidental oligomer

Compare against the HAB1 interface (§9c: H60, F61, K63, I84, S85, G86, L87, P88,
A89, R116, L117, P148, N151, D155, T156, M158, F159, T162, L166).

> **15 of the 19 dimer-interface residues are also HAB1-interface residues.**

Dimerisation and PP2C binding use the *same* face — the gate/latch surface — so
they are mutually exclusive by construction. The apo dimer is a genuine
competitive off-state that suppresses ligand-independent signalling.

Consistent with Mosquna: **H60 alone supplies 9 of the 29 activating
(constitutive) substitutions**, the largest single class. Breaking the dimer
raises basal activity, exactly as this geometry predicts.

### Our six positions do not threaten it

| position | distance to partner protomer |
|---|---|
| K59 | 7.0 Å (sequence-adjacent to interface residues H60, F61, I62) |
| Y120 | 12.2 Å |
| E94 | 13.7 Å |
| E141 | 13.9 Å |
| F108 | 15.4 Å |
| R79 | 15.6 Å |

None is at the interface; five are 12–16 Å clear. **Pocket expansion as designed
does not endanger dimerisation** — but K59 is flanked by three interface
residues, an independent reason (alongside §17a) to distrust Arm 2b at K59
specifically.

### Verdict: keep it, monitor it, do not engineer it

- **Keep**: it is free background suppression, and lower basal means fewer
  FOA-removable constitutives and better dynamic range.
- **Cost**: the dimer competes with the productive complex, lowering apparent
  sensitivity; and **ABA sits only 5.6 Å from the partner protomer**, so a
  substantially enlarged pocket extending toward that face could clash — large
  ligands may require dimer dissociation to load.
- **Do not** introduce mutations at the 19 interface residues unless deliberately
  tuning basal activity.
- ⚠️ **A grafted CoxG-like scaffold (§9g) will not dimerise at all**, so it
  forfeits this background suppression and should be expected to show higher
  basal activity. That is a predictable, quantifiable cost of the grafting route.

---

## 19. WT MD baseline (job 27333712)

### 19a. Why

Three reasons, in order of weight:

1. **Candidate MD is uninterpretable without a WT reference under identical
   conditions.** Every later "does this variant keep its loop dynamics?" run has
   to be compared against something.
2. **It settles §16.** Seven deposited chains say the satellite lobe is sealed;
   one FastRelax says merged. That is a question about whether a side-chain
   rearrangement is thermally accessible — only dynamics can answer it.
3. **Every cavity number in this project is a point estimate from one
   structure.** MD turns cavity volume into a distribution, so variant
   differences can finally be tested against thermal fluctuation — the failure
   mode that has already bitten this project twice (§5b, §14a).

### 19b. Systems (`scripts/38_md_prep.py`, `scripts/38b_md_build.py`)

| system | contents | residues | solvated atoms | waters | Na⁺/Cl⁻ |
|---|---|---|---|---|---|
| S1_apo_open | 3K3K chain A, ligand-free protomer | 178 | 46,820 | 11,042 | 30/30 |
| S2_holo_closed | 3QN1 chain A + ABA | 178 | 43,991 | 10,325 | 28/28 |
| S3_apo_dimer | 3K3K A+B, ABA removed | 366 | 61,008 | 13,866 | 37/37 |
| S4_ternary | 3QN1 A + ABA + HAB1 B + Mn²⁺ | 474 | 111,680 | 26,159 | 71/71 |

All net-neutral. S1/S2 reuse the script-30 files, so MD and Arm 2b share inputs.

**Force field:** ff19SB protein / **OPC** water (the pairing ff19SB was fitted
for) / GAFF2 + AM1-BCC ligand / Li–Merz 12-6 ions. ABA is built from
`data/aba_deprot.sdf`, verified to be in the crystal frame (**0.000 Å per-atom
deviation** from `aba_xtal.pdb`, same atom order) with formal charge **−1**,
matching `00_setup_inputs.py`. Because its coordinates are already correct it is
loaded with `loadmol2` — no placement step. `parmchk2` reported **no ATTN flags**,
i.e. no parameters were guessed.

⚠️ **tleap gotcha:** `leaprc.water.opc` already loads `frcmod.ionslm_126_opc`,
which covers Mn²⁺. Adding `frcmod.ions234lm_126_opc` (TIP3P-era naming) makes
tleap exit with "Could not open file" — that file does not exist for OPC.

### 19c. Measured throughput, and the job shape it forces

Benchmarked on **A100-SXM4-80GB** (`scripts/39b_md_benchmark.sh`):

| system | atoms | **ns/day** | 300 ns |
|---|---|---|---|
| S1_apo_open | 46,820 | **293.5** | 1.0 d |
| S4_ternary | 111,680 | **166.9** | 1.8 d |

**The account allows only 4 concurrent GPUs**, so the campaign is a 12-task array
throttled `%4` — one task per (system, replicate). Packing 3 replicates into one
job would put S4 near 6 days plus equilibration; one replicate per task keeps
every job well inside the limit. Whole set ≈ **4 days wall**.

`--gres=gpu:a100:1` is explicit: the `gpu` partition also holds P100 and K80
(K80 down), and the newest cards (ada6000, blackwell6000) exist **only** under
the preemptible `preempt_gpu`, which is not used here by choice.

**Protocol:** min (restrained → free) → NVT heat 0→300 K / 200 ps → NPT eq
500 ps restrained → NPT eq 1 ns free → **NPT production 300 ns**, 2 fs, SHAKE,
PME, 10 Å cutoff, Langevin γ=2.0, MC barostat, frames every 10 ps. Production
auto-resumes from `prod.rst7`. HMR at 4 fs would roughly double throughput and is
the first optimisation if this must be repeated across many candidates.

### 19d. Pre-registered observables

Fixed **before** the runs finish, and identical to what candidates will be judged
on later, so nothing can be cherry-picked afterwards:

- gate (85–89) and latch (115–117) backbone RMSF and RMSD
- gate–latch minimum heavy-atom distance
- **cavity volume time series** via `lib_cavity.py` on frames — directly
  comparable to every static number in this README
- **satellite-lobe connectivity fraction per frame** — the §16 question
- beyond-the-wall volume distribution (§14d metric)
- W385–ABA water-bridge occupancy (§9a; S4 only)
- **Lβ7α5 (148–156) RMSF** — Dorosh 2013 reports this loop is *more*
  ABA-sensitive than the gate (§17a)

### 19e. Protocol details, ionic strength, reproducibility, ligand parameters

**Stages actually executed per replicate** (`scripts/39_md_run.sh`), six in order,
each skipped if its `.rst7` already exists so a resubmitted job resumes:

| stage | ensemble | length | restraints on solute heavy atoms |
|---|---|---|---|
| min1 | — | 5000 steps | 10 kcal/mol/Å² |
| min2 | — | 5000 steps | none |
| heat | NVT 0→300 K | 200 ps | 5.0 kcal/mol/Å² |
| eq1 | NPT 300 K | 500 ps | 1.0 kcal/mol/Å² |
| eq2 | NPT 300 K | 1 ns | none |
| **prod** | **NPT 300 K** | **300 ns** | none |

Total equilibration before production: **1.7 ns**.

⚠️ 1.7 ns is modest; 5–20 ns is common practice. Rather than fix a discard window
now, **determine it from the data**: find where backbone RMSD plateaus and discard
that much of production. With 300 ns available, discarding 10–20 ns is free. Do
this before computing any observable in §19d.

**Ionic strength** — measured from the built systems, not assumed:

| system | [Na⁺] | [Cl⁻] | ionic strength |
|---|---|---|---|
| S1_apo_open | 177 mM | 152 mM | **164 mM** |
| S2_holo_closed | 184 mM | 151 mM | **168 mM** |
| S3_apo_dimer | 193 mM | 149 mM | **171 mM** |
| S4_ternary | 156 mM | 151 mM | **154 mM** |

Physiological ionic strength is ~150 mM, so these are comparable. They sit
slightly above 0.15 M because neutralising counterions are added on top of the
0.15 M NaCl; Na⁺ exceeds Cl⁻ because all four solutes are net negative.

⚠️ **Caveat for the methods section:** intracellular composition is K⁺-dominated
(~140 mM K⁺, ~10 mM Na⁺), and PYR1/HAB1 are cytosolic/nuclear, so KCl would be
marginally more apt than NaCl. Ionic *strength* — the quantity that governs
electrostatic screening — is correct, and NaCl at ~150 mM is the standard
convention, so the systems were not rebuilt.

**Reproducibility.** Ensemble-level: yes. Trajectory-level: deliberately not.

- `ig=-1` draws a fresh random seed per run. This is intentional — it is what
  makes the three replicates *independent samples* rather than three copies.
- The seed pmemd actually used is written into each stage's `.out` file, so **any
  individual trajectory can be reproduced exactly after the fact.** Harvested
  seeds are tabulated below by `scripts/40_md_provenance.py`.
- `addionsrand` places ions randomly, so a rebuild would not give identical
  coordinates — but `system.prmtop` / `system.inpcrd` are archived, preserving
  the exact systems used.
- Everything else (build, protocol, observables) is scripted and version-pinned.

**Ligand force-field quality — GAFF2 penalty scores.** `parmchk2` reports penalty
scores, the GAFF2 analogue of CGenFF penalties. **Zero ATTN flags**, but 8
annotated parameters, the highest being:

```
hc-c3-cf-ce   same as hc-c3-c2-c2,    penalty score = 324.0
hc-c3-cf-cf   same as hc-c3-c2-c2,    penalty score = 324.0
c3-cf-c2-ha   Same as X-X-ca-ha,      penalty score =  47.1  (general term)
c -cf-ce-ha   Same as X-X-ca-ha,      penalty score =  46.8  (general term)
```

All are **dihedrals** — bonds, angles, charges and van der Waals are clean. The
324s substitute generic sp² carbon (`c2`) for conjugated sp² (`cf`) in methyl
torsions on ABA's conjugated diene.

Assessment: acceptable **for this study**, where ABA acts as a restraint on the
pocket and the observables are protein loop dynamics and cavity volume. It is
**not** acceptable for ligand-focused work — binding free energies, novel ligand
poses, or the productive-vs-non-productive binding question (§17b) — because
these torsions govern ABA's side-chain conformation, i.e. exactly what determines
pocket fit. Fix if needed: a QM torsion scan and refit of those two dihedrals.

Note that CGenFF would not obviously be better: conjugated polyenes with methyl
substituents are a known hard case for both force fields and typically draw high
CGenFF penalties too. The approximation is intrinsic to ABA's chemistry, not to
the choice of antechamber over CHARMM-GUI.

<!-- MD_PROVENANCE_START -->
**Run provenance** (generated by `scripts/40_md_provenance.py`) — 5 replicate directories, 5 with production complete.

| system | rep | stages complete | production | ns/day | GPU |
|---|---|---|---|---|---|
| S1_apo_open | rep0 | 6/6 | 301.7 ns | 293.7 | NVIDIA A100-SXM4-80GB |
| S1_apo_open | rep1 | 6/6 | 301.7 ns | 311.2 | NVIDIA A100-SXM4-80GB |
| S1_apo_open | rep2 | 6/6 | 301.7 ns | 311.0 | NVIDIA A100-SXM4-80GB |
| S2_holo_closed | rep0 | 6/6 | 301.7 ns | 308.8 | NVIDIA A100-SXM4-80GB |
| S4_ternary | rep0 | 6/6 | 32.6 ns | — | — |

**Random seeds** — `ig=-1` draws these at runtime; recording them is what makes an individual trajectory reproducible. To repeat a run exactly, set `ig=<seed>` in the corresponding stage input.

| system | rep | heat | eq1 | eq2 | prod |
|---|---|---|---|---|---|
| S1_apo_open | rep0 | 322767 | 652241 | 881386 | 870734 |
| S1_apo_open | rep1 | 576824 | 853792 | 7655 | 502762 |
| S1_apo_open | rep2 | 269852 | 827165 | 822703 | 172814 |
| S2_holo_closed | rep0 | 3239 | 295676 | 209900 | 48576 |
| S4_ternary | rep0 | 61382 | 960030 | 587401 | 201027 |
<!-- MD_PROVENANCE_END -->

---

## 20. Complete script inventory

Every script, what it does, which env runs it, and where its output lands.
Envs: **E** = `/bigdata/cutlerlab/jjaco081/conda_envs/esmfold2/bin/python`
(numpy/scipy/biopython/rdkit); **T** = `.../tier1_analysis/bin/python` (PyRosetta);
**D** = `.../dockenv/bin/python` (openpyxl + rdkit); **A** = `module load amber/22`.

### Setup and shared libraries

| script | env | purpose |
|---|---|---|
| `lib_cavity.py` | E | shared cavity engine: grid, solvent-excluded criterion, ray-marched buriedness (0.88), connected components. Used by everything that reports a volume |
| `variants.py` | — | the 41-variant panel (`PANEL`) and WT identities (`WT`) |
| `00_setup_inputs.py` | E | builds `pyr1_A.pdb`, `hab1_B.pdb`, `complex_AB_ABA.pdb`, ABA params. Sets the deprotonated carboxylate (charge −1) used everywhere downstream |

### Structural analysis of the pocket

| script | env | purpose | output |
|---|---|---|---|
| `01_pocket_contacts.py` | E | contact map and the buttress network (§2) | stdout |
| `02_rigid_truncation_scan.py` | E | rigid-backbone truncation scan, upper bounds (§3) | `results/02_rigid_truncation.csv` |
| `04_cavity_lining.py` | E | cavity-lining vs distance-to-ligand; line-of-sight test (§2) | `results/04_cavity_lining.csv` |
| `16_cavity_shape.py` | E | PCA extents, max internal span, cross-section (§4c) | stdout |
| `17_ligand_ladder.py` | E | ligand size ladder ABA→C40 (§4c) | stdout |
| `18_exits_and_bottleneck.py` | E | cavity mouths, tail-exit vectors (§4d) | stdout |
| `19_interlobe_neck.py` | E | inter-lobe neck profiling (§4e) | `logs/19_interlobe_neck.txt` |
| `24_state_comparison.py` | E | open (3K3K A) vs closed (3QN1 A) stroke; the 149–156 loop (§9c–d) | stdout |
| `29_second_lobe_access.py` | E | **superseded by 34** — first second-lobe metric, mask-dilation biased | `results/29_second_lobe_access.csv` |
| `34_second_lobe_fixed.py` | E | dilation-free beyond-the-wall volume (§14d) | `results/34_second_lobe_fixed.csv` |
| `37_satellite_across_structures.py` | E | satellite lobe across all 7 deposited chains (§16) | stdout |

### Homology and the superfamily

| script | env | purpose | output |
|---|---|---|---|
| `05_foldseek_search.sh` | — | foldseek vs CATH50 (**must** run under SLURM) | `results/foldseek/` |
| `03_parse_foldseek.py` | E | conservation at the six wall positions (§4) | stdout |
| `08_homolog_cavities.py` | E | Arm 3: fetch full-atom homologs, measure cavities | `results/homolog_cavities/` |
| `10_submit_arm3.sh` | — | SLURM wrapper for `08` compute stage | `logs/arm3.log` |
| `15_carotenoid_start.py` | E | carotenoid-binding START folds; BmCBP (§4c) | stdout |

### Prior art

| script | env | purpose | output |
|---|---|---|---|
| `11_tian_library_overlap.py` | **D** | Tian library design vs our positions (§11c) | stdout |
| `12_tian_combinations.py` | **D** | combination analysis of Tian clones | stdout |
| `13_beltran_criterion.py` | E | reproduces Beltrán's selection rule exactly (§4b) | stdout |
| `14_mosquna_criterion.py` | E | tests R79 against Mosquna's 5 Å + water rule (§4b) | stdout |
| `27_tian_ligand_size.py` | **D** | heavy atoms / length of the 3366-compound library (§12a) | `results/27_tian_ligand_size.csv` |
| `28_tian_ligand_shape.py` | **D** | PMI shape classes, hit rate by size × shape (§12a) | `results/28_tian_ligand_shape.csv` |
| `36_novelty_check.py` | **D** | four-level novelty classification of the panel (§15) | `results/36_novelty.csv` |

### The pilot arms

| script | env | purpose | output |
|---|---|---|---|
| `06_rosetta_cavity_scan.py` | T | **Arm 1** — does the cavity survive repacking | `results/cavity_scan/` |
| `07_rosetta_ratchet.py` | T | **Arm 2 v1** — superseded by 21; kept as the record of the failed assay | `results/ratchet/` |
| `21_rosetta_ratchet_v2.py` | T | **Arm 2 v2** — nstruct 20, gate/latch released (§5c) | `results/ratchet_v2/` |
| `25_rosetta_ratchet_v3.py` | T | **Arm 2 v3** — also frees gate/latch ±3 flanks and the 149–156 interface loop; adds `iloop_bb_rmsd`, `iloop_hab1_min_dist` | `results/ratchet_v3/` |
| `30_prepare_open_closed.py` | E | builds the matched 178-residue open/closed monomers for Arm 2b; verifies identical residue sets and reproduces the 5.2 Å stroke | `data/pyr1_{open,closed}_A.pdb` |
| `31_rosetta_switch.py` | T | **Arm 2b** — closed-vs-open preference, `dE_close` (§14c) | `results/switch/` |
| `09_submit_pilot.sh` | — | SLURM: Arms 1+2 v1, array 0–40 | `logs/pilot_*.log` |
| `22_submit_ratchet_v2.sh` | — | SLURM: Arm 2 v2, work-stealing, 6×8 cores | `logs/ratchet2_*.log` |
| `26_submit_ratchet_v3.sh` | — | SLURM: Arm 2 v3, same shape | `logs/ratchet3_*.log` |
| `32_submit_switch.sh` | — | SLURM: Arm 2b, **rep-major** ordering so partial completion still gives a full panel at low replication | `logs/switch_*.log` |

### Aggregation and statistics

| script | env | purpose | output |
|---|---|---|---|
| `20_aggregate_pilot.py` | E | joins Arm 1 × Arm 2 into the decision table | `results/20_pilot_summary.{csv,md}` |
| `23_merge_ratchet_v2.py` | E | merges chunked Arm 2 output; takes `--dir ratchet_v2\|ratchet_v3` | `results/ratchet_v*/<variant>.json` |
| `33_merge_switch.py` | E | merges Arm 2b, conformer-integrity check, two-sided verdict | `results/33_switch_summary.csv` |
| `35_variance_decomposition.py` | E | **corrected** discrimination test (ANOVA F, true sd) across all arms (§14a) | stdout |

⚠️ `20`, `23` and `33` still print the **old, incorrect** discrimination ratio.
It is retained so earlier conclusions can be re-read against it, but **judge any
arm by `35_variance_decomposition.py`**, not by those columns.

### Molecular dynamics

| script | env | purpose | output |
|---|---|---|---|
| `38_md_prep.py` | E | extracts and cleans the four MD systems | `data/md/S*/protein.pdb` |
| `38b_md_build.py` | **A** | ABA parameters (antechamber/parmchk2) + tleap solvation, two-pass for 0.15 M NaCl | `data/md/S*/system.prmtop` |
| `39b_md_benchmark.sh` | — | SLURM: measures real ns/day on A100 before committing the campaign | `logs/md_bench.log` |
| `39_md_run.sh` | — | SLURM: 12-task array `%4`, six stages, resumable | `data/md/S*/rep*/` |
| `40_md_provenance.py` | E | harvests seeds, throughput, stage completion; **rewrites the README provenance block** | README §19e |

### Scratch (not part of the pipeline)

`smoke.sh`, `smoke2.sh` — early PyRosetta smoke tests, superseded.

### Order to reproduce from scratch

```
00 → 01 → 02 → 04            pocket characterisation
05 (SLURM) → 03              superfamily conservation
08 + 10 (SLURM)              Arm 3 homolog cavities
09 (SLURM) → 20              Arms 1 + 2 v1
22, 26 (SLURM) → 23          Arm 2 v2 / v3
30 → 32 (SLURM) → 33         Arm 2b
35                           corrected statistics over everything
38 → 38b → 39b → 39 → 40     MD baseline
```

⚠ **That list covers §1–§20 only and stops at script 40.** Later work has its own
ordering blocks, because the repo now runs to script 70:

| scripts | what | where |
|---|---|---|
| 41–52 | MD runs, preemption handling, post-relax RMSD | §19, §23 |
| 53–55 | stage-1 Rosetta arm, favour-native sweep, K59 decomposition | §23f–j |
| 56–57 | 804_1 ternary MD | §23 |
| 58–60 | gate/latch loop dynamics | §24, **§29**, order in **§30e** |
| 62–66 | coumarin ground truth, non-cognate docking and builds | §25–§27 |
| 67–68 | 191-residue rebuild and its MD systems | §28, order in **§30e** |
| 69–70 | the open + ABA cell and the factorial run | §30, order in **§30e** |

---

## 21. Current status (2026-08-10) — ⚠ SUPERSEDED BY §39; kept for the record

### Complete

| step | result |
|---|---|
| Pocket characterisation, buttress network, truncation scan | §2–§3 |
| Superfamily conservation (foldseek/CATH50, 266 hits) | §4 |
| Prior-art audit: Tian, Mosquna, Park ×2, Beltrán | §4b, §11c |
| Carotenoid/START survey; CoxG identified as `2pcsA00` | §4c, §10 |
| Cavity shape, ligand ladder, exit vectors, inter-lobe neck | §4c–§4e |
| **Arm 1** — cavity scan, 41 variants ×3 | ✅ true sd 46.9 Å³, F=139 |
| **Arm 2 v1/v2/v3** — 41 variants ×3, ×20, ×20 | ✅ ran; **readout is invariant** (§14b) |
| **Arm 2b** — closed vs open, 41 variants ×12 | ✅ F=9.0, true sd 1.79 REU (§14c) |
| **Arm 3** — 147 full-atom homolog cavities | ✅ §9f |
| Corrected discrimination statistic | §14a, `35_variance_decomposition.py` |
| Ligand-size/shape boundary of the existing landscape | §12a |
| Novelty classification of the panel | §15 |
| Satellite lobe across all 7 deposited chains | §16 |
| Dimerisation interface analysis | §18 |
| MD systems built and benchmarked | §19b–c |

### Running

| job | what | state |
|---|---|---|
| **27333712** | WT MD, 12 tasks `%4`, 4 systems × 3 replicates × 300 ns | running since 2026-08-10 12:56; at ~4 h, tasks 0–3 (S1×3, S2 rep0) were 45–47/300 ns at 299–313 ns/day, matching the benchmark; 8 queued on `JobArrayTaskLimit`. ETA ~4 days total. **S4_ternary finishes last** (task order) |

When it finishes: run `scripts/40_md_provenance.py` (fills the §19e provenance
block with seeds and throughput), then compute the §19d observables.

### Headline results

1. **Expansion does not predict switch cost** (r = +0.034, p = 0.84). Cost is
   residue-specific: K59A/F108A buys +130 Å³ for +0.10 REU; Y120L buys +20 Å³ for
   +4.06. Expansion can be nearly free with the right residues (§14c).
2. **The readout is invariant to pocket mutation** — true spread of
   `dG_separated` is 0.43–0.50 REU against a −67 REU interface (§14b). Exactly
   what the ratchet model predicts.
3. **The existing landscape fails above ~40 heavy atoms and ~15.3 Å**, and
   lutein and β-carotene were screened and failed. Target zone 30–45 heavy atoms,
   15–20 Å (§12a).
4. **Lead candidate `F108A/R79A`** — +88.5 Å³, ΔΔE_switch +0.75, novel on two
   axes: R79 never randomised, F108A allowed but never recovered (§15d).
5. **Screen against Tian's 229 failures** (28–50 heavy atoms, 15.3–22 Å, zero
   hits), not ABA (§13a).

### Not started

| item | why it matters |
|---|---|
| **Cutoff calibration on Mosquna's 741 labelled singles** | the only properly labelled two-class set; needed before any filter rejects anything (§13d) |
| **Insertion / grafting arm** | the only route past the backbone-bounded length limit (§4e, §9g) |
| **Ligand-bound modelling** | non-productive binding is invisible to every current arm (§17b) |
| **4WVO calibration** | no affinity number is trustworthy until PYR1^MANDI ranks above WT (§13e) |
| **Portfolio/diversity selection** | ranking is the wrong instrument for a portfolio (§11a) |

### Open corrections on the record

Three claims in this README were made and then overturned by later data. All are
kept with their corrections, because the reasoning matters:

- §5b "Arm 2 is blind" → §14b: it is **invariant**, not blind; the statistic was
  wrong (§14a).
- §14d "the sealed satellite lobe is a crystal artefact" → §16: **7 of 7
  deposited chains are sealed**; the FastRelax merge is more likely the artefact.
- §9b "R79/E94 opens the second lobe" → §14d: **F108 is the gatekeeper**;
  R79A/E94A adds only +11.9 Å³ beyond the wall. R79's novelty case rests on total
  volume and epistasis instead.

---

## 22. Change log — what we believed, and what changed it

Kept so that superseded reasoning is recoverable and so an overwrite cannot
silently erase why a decision was made. Dates are anchored to SLURM submit/end
times and script mtimes, not recollection.

**Maintenance rule: append here whenever a claim in this README is changed or
reversed. Never delete the old claim — strike it through in place and add a row.**

### Timeline

| date | phase | anchor |
|---|---|---|
| 2026-08-06 13:36 | project start; `lib_cavity.py`, inputs, ABA protonation | file mtime |
| 2026-08-06 15:37 | pilot submitted (Arms 1+2 v1, job 27260956) | SLURM submit |
| 2026-08-06 18:19 | pilot complete | SLURM end |
| 2026-08-07 11:31 | first aggregation; Arm 2 v1 declared "blind" | mtime `20_` |
| 2026-08-07 13:47 | Arm 2 v2 submitted (27275686) | SLURM submit |
| 2026-08-07 14:31 | open-vs-closed state comparison (§9c–d) | mtime `24_` |
| 2026-08-07 16:07 | Tian ligand size/shape boundary (§12a) | mtime `27_`, `28_` |
| 2026-08-07 18:15 | Arm 2 v3 + Arm 2b submitted (27286483/4) | SLURM submit |
| 2026-08-08 04:03–11:50 | v2, v3, Arm 2b all complete | SLURM end |
| 2026-08-10 09:54 | **discrimination statistic found to be wrong** (§14a) | mtime `35_` |
| 2026-08-10 10:50 | satellite lobe checked across 7 chains (§16) | mtime `37_` |
| 2026-08-10 12:56 | WT MD campaign submitted (27333712) | SLURM submit |
| 2026-08-11 11:35 | MD tasks 0–3 complete (S1 ×3, S2 rep0), 300 ns each | SLURM end |
| 2026-08-11 12:17 | 4WVO / 8EY0 downloaded; W385 water re-measured (§23a) | file mtime |
| 2026-08-11 12:40 | ligand-blind oracle exposes the recovery metric (§23b) | analysis |
| 2026-08-11 13:05 | A100 nodes DRAINING; S4 moved to `preempt_gpu` (27386629) | SLURM submit |
| 2026-08-11 13:10 | LigandMPNN bias test submitted (27386657) | SLURM submit |
| 2026-08-11 13:5x | floppy-ligand PDB burial survey launched (§23d) | mtime `44_` |
| 2026-08-12 | stage-1 **LigandMPNN** arm complete (27412775): F108A/F159L recovered, K59R/V81I missed, V81I inverted | SLURM end |
| 2026-08-13 | stage-1 **Rosetta FastDesign** arm (27421344); null arm fails its own control at 29 % WT, K59 0 % | SLURM end |
| 2026-08-13 | favour-native sweep (27438220) fails at every weight; per-residue decomposition (27438909) finds **+10.1 REU** desolvation | SLURM end |
| 2026-08-13 | Sanger: 4/4 colonies are design 804_1 | wet lab |
| 2026-08-14 | **coumarin benchmark** reframes the task — positions saturated, substitutions carry the ligand information (§25) | mtime `62_` |
| 2026-08-14 | loop dynamics analysed over S1/S2/S4 (§24); canonical MD settings fixed; S6–S9 built, not submitted | mtime `58_`–`66_` |
| 2026-08-17 | **PYR1 rebuilt on the full 191-residue sequence** (§28); FastRelax packer trap found | commits `38583c3`, `1372987` |
| 2026-08-18 | **3K3K exposed as a MIXED dimer**; loop dynamics re-run over all 12 trajectories (27545237), apo-closed cell filled by accident (§29) | SLURM end |
| 2026-08-18 | conformation × occupancy **factorial built and launched** (27547457); the open+ABA cell had never existed (§30) | SLURM submit |
| 2026-08-18 | per-position admissibility fails (§31); **pairwise enumeration** gets 3 of 4 into the top 150 of 13,457 (§32); electrostatic screen (§33) | mtime `71_`–`73_` |
| 2026-08-19 | six papers read from their **Methods** (§34); Leonard's protocol needs a known weak hit before it can place a pose | file mtime |
| 2026-08-19 | hydrogen-bond-as-clash fix (§35); **coupled moves** (27560390) — sampling is not the limit (§36) | SLURM end |
| 2026-08-19 | **MM-GBSA makes the K59 flip** (27561665 / 27568930 / 27569546, §37); WIN pre-registered as held-out (§38) | SLURM end |
| 2026-08-20 | MM-GBSA **extended** to the stratified candidate set, 46 tasks, 250 ps ensembles (27676964) | SLURM submit |
| 2026-08-20 | homolog fetch completed (all 66 solved domains); **graft-direction** screen; cavity measurement for all 266 (27684153) | SLURM submit |
| 2026-08-21 | MM-GBSA extended run **completed** (46/46, ~9.6 h each); aggregated with the pool ranking and a convergence audit (75c, 75d) | `sacct` 27676964, `results/mmgbsa/` |
| 2026-08-21 | cavity measurement completed, **261/266** rows (27684153, 1 h 23 m) | `results/homolog_cavities/homolog_cavities_full.csv` |
| 2026-08-21 | measured WHY the implicit reference drifted: ABA slides 2-4 A, K59 salt bridge intact 100%, penalty is EGB not VDW (§41) | `scripts/82`, cpptraj |
| 2026-08-21 | **explicit-solvent rebuild** of both references, ff19SB/OPC/KCl, 10 ns x 3 seeds (27693167) + MM-GBSA rescore (27697744) | `data/mmgbsa_explicit/`, §42 |
| 2026-08-21 | MM-GBSA **retired** as the ranking method (kept as an option); non-cognate MD relaunched at 150 ns without the redundant S9 legs (27697950) | §43a |
| 2026-08-21 | **TI pilot built and launched**: V81I + K59R x {ABA, mandi, apo} x 12 lambda = 72 windows (27698144, 27698165) | `data/ti/`, §43b |
| 2026-08-22 | TI pilot **returns 72/72**; error bar found to be ~20x too small, §43d premise withdrawn, mandi legs found to start from a 0.62 A clash | `results/ti/`, §44 |
| 2026-08-22 | `.gitignore` given **global** extension rules after the 4th per-directory miss left a 7.9 GB `prod.nc` untracked | `.gitignore` |
| 2026-08-23 | **non-cognate MD complete** (9 x 150 ns, 1.35 us); quadruple TI built and 40/48 windows run | §45, `results/noncognate/`, `data/ti_quad/` |
| 2026-08-24 | **quadruple TI 48/48**: crystal arm -1.12 +/- 2.34 (no call); docked arm invalid, ligand collapses into ghost F108 | S46, `results/ti_quad/` |
| 2026-08-24 | crystal arm **extension launched** (27725295, 24 windows to 120 ns); project-wide reflection written | §47 |
| 2026-08-24 | sensor-chemistry benchmark built; PFAS/TNT frozen as a held-out prospective test | §48, `scripts/106_sensor_benchmark.py` |
| 2026-08-24 | charge identity shown ligand-dependent but pose-free; K59 retention test pre-registered on PFAS/TNT | §49 |
| 2026-08-24 | lookup baseline run and killed; ligand-blind frequency null quantified as the bar | §50, `scripts/107_lookup_baseline.py` |
| 2026-08-24 | ground truth reframed as positive-unlabeled; metric switched to hit-retention vs library size | §51 |
| 2026-08-24 | incumbent corrected to the two-round process; potency bias found; target set at 17 % of 1 uM clones at ~80K members | §52 |
| 2026-08-24 | library sizes corrected for substitution depth; round 2 reframed as narrow+deep, not smaller | §53 |
| 2026-08-24 | round-1 -> round-2 carryover measured; design problem reframed as vocabulary subset selection | §54 |
| 2026-08-24 | co-folding validation submitted (27727052, 4 runs incl. WT negative control); ML data audit | §55 |
| 2026-08-24 | donor graft pockets measured (2PCS 27 lining side chains vs PYR1 19-20); switch still unsimulable | §56 |
| 2026-08-24 | donor ligands checked (2PCS = UNL, 3TFZ = buffer); five routes around library x library ranked | §57 |
| 2026-08-24 | factorial found to be n=1 not n=3; core-mask atom mismatch found; under-filled pocket problem stated | §58 |
| 2026-08-24 | partial-occupancy question answered from the clones; enhanced-sampling and parallel-evolution designs specified | §59 |
| 2026-08-24 | umbrella-sampling calibration set up and seeded (62 windows, WT+ABA vs WT apo) | §60, `scripts/111-113` |

### Reversals and corrections

| # | date | we thought / did | what changed it | now |
|---|---|---|---|---|
| 1 | 08-06 | PYR1–HAB1 is a **molecular glue** | HAB1 contributes **zero** pocket residues; it reads the closed-state surface | **ratchet** (§1) |
| 2 | 08-06 | pocket residues = those near the ligand (**distance-to-ligand**) | user's ChimeraX cavity view: R79 clearly touches the cavity | **cavity-lining + line-of-sight** (`04_`, §2). E94 was mislabelled second-shell (it is first shell, 4.33 Å); R79 is a cavity-rim residue |
| 3 | 08-06 | "F108 was only ever made **larger**" | PDB **4WVO** = PYR1^MANDI carries **F108A** | true of the Tian datasets only (§4b) |
| 4 | 08-06 | Arm 3 cavities were all 0.0 Å³ | foldseek stores Cα only → `convert2pdb` emitted backbone traces | rebuilt with real full-atom fetch (§9f) |
| 5 | 08-06 | the bottleneck is at the **gate/latch** end | user: "the one **between the two lobes**" | there is **no inter-lobe bottleneck**; the limit is length (§4e) |
| 6 | 08-07 | **"HAB1 never touches the ligand"** | water **A197** bridges W385 NE1 (3.04 Å) to ABA O10 (2.72 Å) | HAB1 senses the ligand head at one remove (§9a) |
| 7 | 08-07 | **Arm 2 is statistically blind** | corrected statistic at nstruct 20 | Arm 2 is **invariant, not blind** — true spread 0.43–0.50 REU vs a −67 REU interface (§14b) |
| 8 | 08-07 | win = positions **outside the Tian libraries** | user: substitutions and *combinations* also count | four-level novelty; only **1 of 40** variants is a true retread (§15) |
| 9 | 08-07 | over-corrected: "shrinking Tian residues **is not a win**" | user: a maximally open pocket is a valid **design starting point** | both true — not a deliverable, but a needed substrate for shrink-to-fit (§11b) |
| 10 | 08-07 | rank variants on **ΔCavity** | user: the goal is a **portfolio of scaffolds**, not a winner | select for diversity; ΔCav is a means (§11a) |
| 11 | 08-07 | sized SLURM jobs by **walltime** | user: walltime is abundant (30 d); **concurrency** is the limit | job shape set by concurrent-job and GPU caps (§19c) |
| 12 | 08-07 | nodes were "**draining**" | `State=MIXED+PLANNED`; `sinfo -R` shows no drain reason | planned/backfill, not drained |
| 13 | 08-07 | array pending because **CPUs** were full | r11 had 20 idle CPUs but only 7.6 GB allocatable | **memory** was the binding constraint |
| 14 | 08-10 | discrimination = between-sd ÷ **replicate sd** | that denominator ignores √n | ANOVA `F = s²_between/(s²_within/n)`; flipped Arm 2b from "blind" to F=9.0 (§14a) |
| 15 | 08-10 | the sealed satellite lobe is a **crystal artefact** | **7 of 7** deposited chains are sealed | the **FastRelax merge** is the more likely artefact (§16) |
| 16 | 08-10 | **R79/E94 opens the second lobe** | dilation-free metric: R79A/E94A adds only **+11.9 Å³** beyond the wall | **F108 is the gatekeeper**; R79's case rests on total volume and epistasis (§14d) |
| 17 | 08-11 | because the W385 water senses ABA, requiring ligand contact with it is a sound design filter | measured in 3QN1 / 4WVO / 8EY0: the water is present in **all three**, with W385, Pro88 and Arg116 contacts conserved to 0.3 Å, but **mandipropamid never touches it** (5.1–5.2 Å vs ABA's 2.72 Å) | the water still senses ABA, but it has **two separable roles** and only the **gate–latch–HAB1 staple** is conserved; sensing is not. Leonard et al.'s H-bond constraint would have **excluded mandipropamid** — use it to rank, never to exclude (§23a) |
| 18 | 08-11 | LigandMPNN's 28.7% top-1 recovery on the coumarin set measures method quality | a **ligand-blind oracle** on the same labels reaches **79.1%**, and **all 33** MPNN hits fall at positions where the experiment is ligand-invariant | sequence recovery is the wrong metric — most of it is winnable without the ligand. The number is separately confounded by pose (§23b) |
| 19 | 08-11 | `39_md_run.sh` was safe to repoint at any partition | its resume path overwrites `prod_cont.nc`, and `prod.in` reruns a **full** 300 ns because `irest=1` continues the clock | harmless on non-preemptible `gpu`, destructive under preemption; `42_md_run_preempt.sh` numbers segments and computes the remainder from the restart clock (§23c) |
| 20 | 08-12 | **stage 1 showed LigandMPNN cannot recover PYR1^MANDI** (§23g) | both mandipropamid arms held **no ligand** — shifted HETATM columns put `3` in the altLoc field and ProDy's default `altloc='A'` dropped all 29 atoms; ABA survived only because `A8S` starts with `A` | §23g **retracted**. Corrected run: LigandMPNN *does* find **F108A (+0.111)** and **F159L (+0.071)**, misses K59R, and V81I inverts to **−0.841**. Validate inputs by parsing them with the consuming library and asserting on counts (§23h) |
| 21 | 08-12 | a protocol that gives consistent numbers across trajectories is working | `NeighborhoodResidueSelector` measures neighbour-atom distances, and a 29-atom ligand has one NBR atom, so an 8 Å shell caught **3 residues**; the frozen pocket made every trajectory identical and dG a tidy −12.29 | all-heavy-atom shell → 43 residues, dG −19.11, trajectories that differ. **Suspiciously low variance is a bug signature, not a quality signal** (§23i) |
| 22 | 08-13 | the ligand-swap null only had to supply a subtraction, so it did not need checking itself | the ABA null arm is WT + native ligand, so its own correct answer is **zero mutations** — it was at **4/15**, deleting a 2.85 Å K59 salt bridge in 100% of trajectories. Never scored, because a null is only ever read as a difference | **score the control's own correct answer, not only the contrast.** K59R and V81I are uninterpretable, not negative; stage 1 does not clear (§23j) |
| 23 | 08-13 | a params file summing to 0.000 net charge meant the formal charge was dropped | field 4 of a params `ATOM` line is the MM type — the literal `X` — and the charge is field 5. Summing field 4 returns 0.000 for **any** ligand. The real sums were −0.970 and +0.100, both correct | caught by the validator written to act on it, before any commit or compute; params regenerated **byte-identical**. Validation now uses field 5 with a rounding-aware tolerance. **A number shaped like a finding still needs its reader checked** (§23j) |
| 24 | 08-14 | one residue-numbering map would serve every MD system, since they are all PYR1 | S1/S2 come from script 30's 3K3K∩3QN1 intersection (178 res, gate = seq 82-86); S4 was built from 3QN1 alone and keeps residue 2 as ALA (179 res, gate = seq **83-87**). Reusing S1's core mask on S4 gave a 162-residue fit set against the references' 161 — cpptraj **set it up anyway** and returned gate RMSD ≈ **84 Å** | maps are rebuilt **per system** and verified by residue identity before any frame is read; the core set is asserted identical across systems. Also: address ligands by NAME (`:A8S`), since `n_res+1` is ABA in S2 but HAB1's first residue in S4 (§24f) |
| 25 | 08-14 | the gate–latch contact is the staple that holds the closed state | in the closed MONOMER the contact is intermittent — broken in **60 %** of frames in S2 rep0 — and it is the weakest pre-registered observable (replicate gap 0.90 Å). Adding HAB1 gives the tightest distribution of any system (3.72 Å) | the staple is clamped by the PARTNER, not held by PYR1 alone — consistent with the ratchet framing. Do not filter designs on gate–latch distance (§24d) |
| 26 | 08-14 | a stable closed trajectory would show a design's switch works | neither state converts even once in 1.8 μs of aggregate WT sampling; both are kinetically trapped at 300 ns | MD licenses a **stability** filter (gate RMSD to closed, threshold 3.5–4.0 Å, ~1 % error each way) and is **blind to switchability** — the more likely failure mode for an enlarged pocket (§24b) |
| 27 | 08-17 | 3QN1 breaks at C(1)–N(3) as well as 68→71, and 3K3K shares those gaps | audited both CIFs per chain: **3K3K is entirely gap-free** (A 1–183, B 2–184, zero broken bonds) and **3QN1 chain A misses only 69–70**. The 3.02 Å break is in OUR `data/md/S1_apo_open/protein.pdb`, which script 30 made by intersecting the two crystals; 3QN1 does model residue 2, as the P2A **mutation** | a defect in an intermediate of ours was attributed to the PDB. PYR1 needs **one** 2-residue loop, not two, and the open form needs none (§28a) |
| 28 | 08-17 | 3K3K is the apo dimer | **chain B has ABA bound** (canonical 20-residue pocket) and is **CLOSED** — gate 0.90 Å / latch 0.68 Å from the 3QN1 closed reference, against chain A's 5.24 / 5.08 Å | 3K3K is a **mixed** open-apo + closed-holo dimer. `S3_apo_dimer` (3 × 300 ns, already run) simulated its closed protomer with the ABA **stripped**. §24 is unaffected — script 58 covers S1/S2/S4 only (§28g) |
| 29 | 08-17 | a frozen backbone means the crystal structure was preserved, and a MoveMap freezes what it names | core preservation reported **0.000 Å** while two defects went through: A87–B116 at **0.36 Å** (crystal 3.09 Å) and **K59 collapsing 2.85 → 1.68 Å** into the empty ABA cavity — K59 was in an explicit "frozen" set and moved anyway | **`FastRelax.set_movemap()` restricts minimisation only**; repacking needs a TaskFactory or the whole pose is repacked. Backbone RMSD also watches N/CA/C only, so both were invisible. Fixed in `lib_rosetta.py`; 178/191 residues now keep their input rotamers, and the dimer assembles with **zero** clashes before any refinement (§28d) |
| 30 | 08-18 | the §24b closed-state stability filter reports something about the ligand | the apo-closed protomer (S3 chain B, closed and ligand-shaped, pocket empty) holds its state for **3 × 300 ns with zero crossings**, drift **−0.07 Å**, and has the **least mobile gate of all fifteen units** (RMSF 0.72 Å vs 1.39 with ABA) | the filter is **conformation-reporting, not ligand-reporting**: a good gate-RMSD-to-closed score carries no evidence about occupancy. Confounded by the dimer, which alone drops the open gate 3.05 → 1.07 Å, so `S9_apo_closed` (the apo-closed **monomer**) is now the highest-value unrun system (§29c) |
| 31 | 08-18 | the apo-closed cell could be filled by running the S9 built on 08-14 against the existing S1/S2 | the old-tree S9 is **KCl** (K⁺ 33/Cl⁻ 28) while S1/S2 are **NaCl** (Na⁺ 35/34) — read from the topologies, not the build logs | pairing them would confound ligand removal with a **cation swap** inside the one comparison the factorial exists to make. All four cells rebuilt on `md191`; job 27547457, 12 tasks (§30a) |
| 32 | 08-18 | per-position steric admissibility could generate ligand-specific menus, since geometry excludes rather than ranks | it admits **14.7 of 20** residues at a typical pocket position, ABA and mandipropamid menus overlap at **Jaccard 0.84**, and all four mandipropamid ground-truth substitutions are admissible **for ABA too** — recalled without being ligand-conditional, exactly like the stage-1 clash baseline | the per-position step is dead, but the diagnosis is useful: every miss (V83W, V164W, V164H, V83F, A89V, E141F) is a **grow** mutation, and the real sensors are **compensating shrink/grow pairs** (F159A+A160I, Y120G+A160G). A per-position filter judges each mutation in a context where its partner has not happened. Enumerate **pairs** against a re-evaluated pocket (§31c) |
| 33 | 08-18 | evaluating mutation pairs jointly would expose the compensation a per-position filter misses | with a **summed** (hence additive) ligand-overlap objective it exposes nothing: relief of every top pair equalled the sum of its singles, and the best residue at a position changed with its partner in **15 of 4187 contexts**. Adding a packing term that counts only NON-overlapping contacts recovers **F108A rank 5, F159L rank 75, V81I rank 150 of 13,457**, with shrink/grow enriched to 62 % vs 42 % for the ABA control | pairs work, but only once the objective can represent the grow half. The remaining miss, K59R, is blind **by category** — it does not clash, so a steric method cannot reach it; that is the same residue and the same reason that defeated both stage-1 methods (§32d) |
| 34 | 08-18 | the K59R miss is a scoring problem — ref2015's desolvation in §23j, category blindness in §32 | transplanting the **crystal** Arg59 from 4WVO scores **2 hydrogen bonds** (NE–O2 2.64 Å, NH1–O2 3.26 Å) with a worst genuine overlap of 0.56 Å, **below** the feasibility threshold — the scoring accepts it. Its χ3 ≈ **100°** is non-rotameric and the backbone-independent library never proposes it | at this stage the miss is **SAMPLING**, not scoring, and chi minimisation cannot fix it: relieving strain does not create a hydrogen bond. Also found: hydrogen bonds are being counted as steric clashes (donor–acceptor at 2.64 Å = 0.43 Å hard-sphere overlap), which reaches back into §31/§32 and likely explains §31e's R79/V83/H115 control failures (§33b–c) |
| 35 | 08-19 | Leonard et al.'s protocol is a general computational route to new PYR1 specificity | read from their Methods: the pose comes from **"dock to sequence"**, in which the WT sequence is first **mutated to match a known low-affinity binder found by Y2H screening**. It also hard-filters on an H-bond to the latch water, which §23a measured mandipropamid to miss by **5.07–5.19 Å** | the method **improves a known weak hit**; it cannot start from a ligand alone, and applied to mandipropamid it would reject the correct pose. Our benchmark is strictly harder and comparisons must say so (§34a) |
| 36 | 08-19 | the hydrogen-bond-as-clash bug explained §31e's R79/V83/H115 control failures | with the exemption applied those three fail at the **same** overlaps (0.65/0.69/1.12). The crystal R79 is now clash-free (0.306 → 0.000 Å) but no library **rotamer** reproduces it | §31e's original diagnosis — wrong chi basin — was right and my §33c speculation was wrong. The fix is still correct and matters at hard-sphere tolerance (WT-ok 12 → 16 of 26 at tol 0), but §31's and §32's verdicts are unchanged (§35) |
| 37 | 08-19 | the K59R miss might be a sampling limit that backbone flexibility would fix | **coupled moves** — the method Kortemme built for exactly this benchmark, sampling sequence + side chains + backbone + ligand pose — recovers F108A at **+0.74** but retains K59 in **0 % of 50 trajectories in the ABA arm**, where the cognate ligand makes K59 unambiguously correct. Identical to fixed-backbone FastDesign | **scoring, confirmed by a second independent sampler.** Adding flexibility cannot recover K59R; any method using ref2015's desolvation on buried charge inherits it. Coupled moves still helps everywhere else — WT retention 29 % → 42 % (§36a) |
| 38 | 08-19 | no scoring function available to us can make the K59 flip | **MM-GBSA can.** With mandipropamid, R at 59 beats Q and N (−1.53 vs +0.04, +2.53); with ABA, wild-type Lys beats every substitution (R least-badly at +7.06). Generalised Born instead of Lazaridis–Karplus, over 100-frame MD ensembles | the §23j desolvation diagnosis was not just correct but **actionable** — changing the solvation model does what no sampler could. ⚠ Ensembles are only 50 ps, the ligand-swap difference is dominated by damage to the ABA complex, and F108A comes out NULL because relaxation absorbs the clash it exists to relieve (§37a–b) |
| 39 | 08-20 | the MM-GBSA result at position 59 was enough to call it a method | it was a **3-way within-position** comparison among K/Q/N, on 4 mutations we already knew the answers to. Whether it RANKS correctly inside a 23-variant pool is untested, and the raw pairwise top-24 all contain position 108, so a naive extension would have tested 'which partner goes with F108X' | extended to the **stratified** candidate set (best pair per distinct position pair, 14 positions instead of 5) with 5x longer ensembles, job 27676964. A within-position win is a signal; a ranking inside a pool is a method (§39b) |
| 40 | 08-21 | MM-GBSA made the K59 flip (belief 38) | the flip was an artefact of an **unconverged reference**. The WT–ABA complex scores −32.74 at 50 ps and −25.00 at 250 ps — a **+7.75 kcal/mol** move, the largest of any run — and it is *still* drifting +2.62 within the 250 ps window. Every ddG in the ABA arm is measured against it. Recomputed at 250 ps the pre-registered test gives the **opposite** answer: selectivity −8.58 (PASS) → **+1.31 (FAIL)** | **belief 38 is retracted.** Sign-based MM-GBSA verdicts are void at these lengths. Rank-based ones survive (a constant shift cannot reorder an arm) and the four known mutations do land at ranks 1/3/4/18 of 22, p = 0.049 — but that is marginal, post-hoc, and bounded by a median per-variant drift of 1.05 kcal/mol, the same size as the spacing it is ranking on. **WIN stays sealed** (§40) |
| 41 | 08-21 | the 250 ps implicit numbers were the better-converged ones, so the ABA arm just needed longer runs | explicit solvent puts WT-ABA at **−32.46 ± 2.72**, which matches the **50 ps** value (−32.74), not the 250 ps one (−25.00). The 250 ps "convergence" was the ligand sliding 4 A out of a pose only explicit water can hold — drift AWAY from the answer, not toward it | **§40's verdict stands but its reasoning was wrong.** Implicit solvent is unreliable for this anion at *every* length tested; length was never the variable. Neither implicit number should be quoted |
| 42 | 08-21 | fixing the solvent would make MM-GBSA usable for ranking | explicit solvent fixes the **pose** (ABA max RMSD 4.88 → 2.5-3.5 A) but not the **precision**: between-seed sd is **2.72** kcal/mol for ABA, and one seed drifts **+10.09** within its own 10 ns. A ddG carries ±3.85 from one replicate each | the ranking must resolve 0.01-1.8 kcal/mol. Reaching ±0.5 needs **n=60 per variant** = 2760 runs = **1748 GPU-h (18 days on the 4-GPU cap)**. Pose was never the binding constraint — variance is (§42d) |
| 43 | 08-21 | the way to rescue a noisy ddG is a better solvent model or more replicates | the variance is a property of **subtracting two separately-computed absolute energies**. TI computes the same ddG as a *difference along a path*, so the ~43,000 unchanged atoms cancel by construction rather than being subtracted — and for SELECTIVITY the apo leg cancels exactly, leaving two legs instead of four | switched to TI at ~40 ns/leg vs MM-GBSA's ~1200 ns/variant. Not a precision upgrade bought with compute — a different estimator with a smaller variance by construction (§43b) |
| 44 | 08-22 | TI's tight error bar meant TI was precise | the bar was `sd/sqrt(n)` with **n=4004**, and n was wrong twice: pmemd prints every frame **twice** (bit-identical DV/DL, 0 of 2000 steps disagree) and the trailing AVERAGES / RMS banners were parsed as samples. Real n = **2000**, correlated to **tau = 75**, so `n_eff` falls to **13** in the worst K59R windows | the +/-0.06 on V81I was ~20x too small. 88 now reports **stat / conv / quad** separately (89_ti_reparse.py); quadrature is clean (0.00-0.01), V81I is converged, **K59R is not** and its second-half estimate moves further positive, not toward zero (§44a) |
| 45 | 08-22 | both V81I and K59R should individually shift selectivity toward mandipropamid, because both appear in 4WVO (the §43d premise) | **PYR1^MANDI is a FOUR-mutation set** (K59R/V81I/F108A/F159L) selected together; nothing requires a member to work alone in a WT background. Measured here, native V81 is **4.85 A from ABA** -- second shell, no contact (native V83, also a valine, is the one at 2.99 A). LigandMPNN had independently scored V81I at **-0.841**, anti-correlated with mandipropamid (S23h) | **§43d is withdrawn as a test of the method.** A converged TI, a sequence model and the structure all agree V81I *alone* is not mandi-favouring. The premise was mine and should have been challenged when written, not after it returned an unwelcome answer (§44b) |
| 46 | 08-22 | the TI mandipropamid legs represented a bound complex | 85 loads `3UZ.mol2` into the **WT** pocket with F108/F159 present and **no pose relaxation**: F108-ligand **0.62 A** heavy-atom, F159 1.38 A, V81 1.40 A, with 8 (V81I) and **18** (K59R) contacts under 2.0 A. Minimisation relieved it and nothing dissociates (ligand RMSD 1.1-3.9 A after CA superposition), but the resulting pose is validated against nothing | **every selectivity number rests on that leg.** The decisive run is the one §13e already requires and TI never had: the **quadruple** K59R/V81I/F108A/F159L, ABA vs mandipropamid, built from **4WVO's own coordinates** -- 24 windows, ~2.5 h (§44c, §44e) |
| 47 | 08-23 | if closed PYR1 held its state around a ligand it was never built for, the §24b filter would be conformation-only and disqualified | it is **not** purely conformational: the gate moves ~1.0 A more around imperatorin (+1.02) and flutamide (+1.08) than around ABA, run-level ranges non-overlapping. But **alpha-estradiol is indistinguishable from the cognate ligand** on every observable (ligand RMSD 1.54 vs 1.67, gate 1.69 vs 1.79, exact p 0.90) despite needing three mutations to work as a sensor | the filter carries ligand information but **passes a true negative**, which is the dangerous direction for a design filter. And the design cannot be significant: n=3 vs n=3 has a two-sided p floor of **0.10**, exactly where both "separate" results sit (§45)
| 48 | 08-24 | the quadruple TI would settle whether TI can rank selectivity | the CRYSTAL arm gives **-1.12 +/- 2.34** -- the right sign, but the error is twice the effect, so **NO CALL** on the S13e test; second-half-only moves to -0.10. Quadrature is clean (0.01) and the apo leg's 4.61 drift cancels exactly out of selectivity, as designed | still undecided, and now COSTED: stat 0.5 needs ~121 ns/window = 2.9 us = ~8.3 GPU-days over 24 windows (S46a) |
| 49 | 08-24 | a docked pose would just give a worse NUMBER than the crystal pose, letting us price "no structure" | the docked arm returned **-1381 kcal/mol** with the wrong CURVE SHAPE. Measured: the ligand collapses into the decoupled WT Phe108 ring -- **all 14 sub-2 A contacts at lambda=0.885 are F108 ring atoms**, min 1.22 A, while the crystal arm has zero in every window | it fails STRUCTURALLY, not noisily. The pose was docked into the F108A **cavity**, i.e. exactly the volume the transformation deletes, so **any TI from WT to a cavity-creating mutant with the ligand docked into that cavity is ill-conditioned in dual topology** -- which is the intended production workflow. The "cost of no structure" is not quantifiable this way (S46b) |
| 50 | 08-24 | §32's 3-of-4 recall meant a steric library would be about 75 % complete | measured over 691 + 78 + 45 characterised sensors: **steric is only 37.5-50 % of distinct substitutions** in every set, CHARGE is 26-30 %, polar 20-35 %. Position 59 carries **8-11 % of all substitutions** in every set and is **89 % charge chemistry** (185 of 208) | a steric library is not 75 % complete, it is **blind to a category**. But covering it is cheap: **K59 alone is 44 % of all charge chemistry and seven positions are 90 %**, so a position-typed menu fixes it at a cost smaller than the libraries Tian already builds (§48a-b) |
| 51 | 08-24 | choosing WHICH residue at a charge position needs the ligand's charged group located, hence a pose -- which §46b measured at 8.33 A error into a pocket that does not exist yet | holding the library fixed, 59-charge substitutions per clone run **0.05 (anionic ligands) -> 0.17 -> 0.31 -> 0.78 (cationic)**, a ~16x spread in the direction electrostatics demands, and acids relocate their charge chemistry to E94/V81/I110 instead. 17 of 19 acid ligands never touch 59 | **the circularity is binding for STERIC placement and largely absent for CHARGE** -- the discriminating feature is formal charge class, a 2D SMILES property needing no pose. Polar positions (120, 163) are NOT covered by this and plausibly do need geometry (§49a-b) |
| 52 | 08-24 | a chemical-similarity lookup over 194 already-screened ligands might match anything we build, so it had to be ruled out first | leave-one-ligand-out over 89 ligands: similarity beats a **ligand-blind frequency menu** by **+0.010** recall@20 (23 win / 13 loss / **53 tie**), and copying the single nearest ligand is far WORSE (0.20 vs 0.35). It helps only where a close neighbour exists (+0.045 at Tanimoto>=0.5) and **only 0.4 % of ligand pairs reach 0.5** (median 0.110) | **lookup is dead** -- nothing to look up for a novel ligand. But the null it exposed is the real prize: **freq recall@20 = 0.35, @40 = 0.52**, ligand-blind and free, and the bar every method here has never been measured against. oracle@20 = 0.99 also scopes the task to a 20-40 substitution menu (§50) |
| 53 | 08-24 | recall of observed substitutions was the right way to score a designed library | the hit sets are **positive-unlabeled**: the landscape is not fully tested, libraries are SAMPLED not enumerated, and selection removes variants for reasons unrelated to pocket binding (URA3 activity, constitutive interface binding, and non-uniform promiscuity filtering -- pan-PFAS cross-reactivity was a FEATURE there). 'Not observed' means unknown, not non-functional | **precision is unmeasurable, only recall**; the frequency baseline is **advantaged by construction** so failing to beat it is weak evidence; `oracle@20=0.99` was misread. Metric replaced by **fraction of ligands whose library contains >=1 hit, vs library SIZE** -- you need *a* sensor, not every sensor. New bar: **58 % at 10^5.9, 72 % at 10^7.3**, matching Tian's own focused-library sizes (§51) |
| 54 | 08-24 | hit retention vs library size was the metric that would show a designed library winning | every ground-truth hit came OUT of the existing libraries, so a designed menu is a SUBSET and can only lose hits -- ceiling 100 % at full size. And hit RATE is not recoverable at all: the data record characterised hits, not screening depth. The real incumbent is free (DSM glycerol stock) then a **focused ~1e5 round-2 library built from round-1 hit profiles** -- Tian coumarin **77,327**, TNT **506,229** -- which demonstrably works, rescuing 5 ligands round 1 missed | **retention can only measure shrinkage, never advantage.** The novel claim available is different and stronger: Tian's focused libraries NEED a round-1 screen, so designing one from chemistry alone **removes an experimental round**. Measurable headroom found: frequency ordering is biased to WEAK sensors -- at ~80,000 members it captures **17 % of 1 uM clones vs 24 % of 100 uM** (§52) |
| 55 | 08-24 | library size = product over positions of (allowed residues + 1) | the primaries are **substitution-DEPTH limited**: DSM-Hao is a **double**-substitution library and TSM a **triple**, confirmed in the clones (dsm mode 2, 193/266; tsm mode 3, 346/403). True sizes are **DSM-Hao 36,140** and **TSM 332,863**, not 3.8e21 and 7.7e15 -- wrong by up to **16 orders of magnitude**. Every library is 1e4-1e7 and the PRIMARIES are the smallest | **round 2 is not smaller, it is differently SHAPED** -- broad+shallow (18 positions, 2-3 deep) becomes narrow+DEEP (11-14 positions, 7-8 deep), reaching combinations round 1 cannot express at any screening depth. '>=10x size reduction' is the wrong axis; the task is **which positions are worth combining deeply**. Retention percentages survive (set containment); only the size axis was wrong (§53) |
| 56 | 08-24 | a secondary library is built by extrapolating from that ligand's own round-1 hit | measured on the 11 coumarin round-2 sensors: **own-ligand carryover is 0-64 %** (median ~23 %, and **4 of 11 had NO round-1 hit at all**), while **pooled round-1 across 194 ligands covers 75-100 %** (mostly 93-100). Round-1 explores 144 of 342 possible single substitutions and contains 75 % of all round-2 chemistry; the residual 8 sit at positions 65/71/74/109/124/134/178/184, OUTSIDE the original 18 | **an initial hit is neither necessary nor sufficient.** The design problem is **subset selection from a known ~144-substitution vocabulary plus a depth choice**, not extrapolation -- use the hit to WEIGHT the vocabulary, never to restrict it, since restricting would have failed for every ligand in the table (§54) |
| 57 | 08-24 | ~1150 labelled clones might be enough to bias an ML tool from ligand SMILES to sequence | the effective sample size is the number of independent **LIGANDS (~208)**, not clones -- clones sharing a ligand are repeats under one input. 69 ligands have exactly 1 clone; only 47 have >=5; median pairwise Tanimoto **0.106** | **generative SMILES->sequence is 2-4 orders of magnitude short**, and LigandMPNN fine-tuning inherits the pose problem AND starts from a model that INVERTS V81I (-0.841, §23h). What IS supported: a ~10^2-parameter conditional model over position x residue-class (§49a generalised), and -- unused so far -- a **tractability classifier** on 208 positives + **229 documented negatives** = 437 ligand-level labels (§55b) |
| 58 | 08-24 | the big-cavity graft donors offer PYR1's pocket with more room | measured for the first time from the DONOR side: cavity volume and lining side chains track at **r = +0.99**, so **2PCS lines its 570 A^3 cavity with 27 side chains against PYR1's 19-20** (~40 % more positions) at 9.8 % identity, and volume per lining residue rises from ~9-12 to 21 A^3 | they are **different architectures, not larger PYR1s**. The cavity/transplantability tension is **categorical, not gradual**: everything >300 A^3 sits at 3.45-4.03 A core RMSD and ~10 % identity, everything transplantable (<1.5 A) is PYR1-sized. Also **2NS9 and 2BK0 are APO**, so the previously-favoured 2BK0's 344 A^3 is a cavity-detection number being compared against a ligand-contact one (§56b-c) |
| 59 | 08-24 | grafting from a donor with a KNOWN ligand collapses the chemical dimension and avoids library x library | checked the het codes: **2PCS's ligand is literally `UNL`, "Unknown ligand"** (unmodelled density), 3TFZ's `CXS` is **CHES buffer**, and 2NS9/2BK0 are **apo**. Only **6AWV / (-)-epicatechin** has an identified physiological ligand -- and it is the weakest of the big four to graft (3.93 A coreRMSD, 9.8 % identity, **17/19** machinery coverage) | the premise holds for **one** candidate, which is also the hardest graft. And the product is a **DISCOVERY** cost that does not apply when a weak hit already exists -- then it is one library x one analyte, i.e. §54's vocabulary-subset problem. Grafting's remaining distinct use is as a **positive control that a big pocket can switch at all** (§57) |
| 60 | 08-24 | the §30 conformation x occupancy factorial was complete (4 cells x 3 reps x 300 ns) and merely unanalysed | I checked that `prod.nc` EXISTED, not how long it was. **5 of 12 reps reached 300 ns -- one per cell**; three array tasks failed outright and four more truncated at 6-65 ns. Separately the core superposition mask selected **656 atoms in the system vs 644 in the references** (both refs miss residues 2, 69-70, 182-191), so the fit silently did not happen and gate RMSDs came out at **39-41 A** for a five-residue loop | the 2x2 is **n=1 per cell** and cannot give the replicated S9-vs-S2 contrast it was built for. Same 'verify against the physics, not the file' failure as §44a's frame count. Corrected 161-residue common core written to `data/md191/core_mask.txt` (§58b) |
| 61 | 08-24 | a weak hit collapses the chemical dimension, so expanded-pocket designs can be tested against it | **the gate closes ONTO the ligand.** A weak hit is for a ligand that fits the CURRENT envelope; enlarge the cavity and that ligand no longer reaches the gate, giving either no closure (no signal) or ligand-independent closure (constitutive, which the counter-selection removes) | expansion must be paired with a **LARGER** ligand -- the 28-50 heavy-atom band, i.e. the 229 documented failures, which are the matched test set. And §24's kinetic trapping means we can measure the stability of a starting state but **never the open/closed free-energy difference**, which is the quantity separating signal from constitutive. **The question is one this class of method cannot answer** -- it needs Y2H (§58) |
| 62 | 08-24 | an enlarged cavity needs a larger ligand to fill it, or transduction breaks (§58) | measured over 691 clones: **94 release >=100 A^3 of side-chain volume, and 28 of those bind ligands of <=20 heavy atoms** -- ABA-sized. **Honokiol and Magnolol (isomers, 20 heavy) both give 1 uM sensors at -161 A^3**, close to doubling PYR1's 174 A^3 cavity; Carpropamid 1 uM at -144 | **the ligand does NOT have to fill the enlarged cavity.** Partial occupancy reaches the top potency class empirically, so the under-filling worry is materially weakened -- with the caveat that net side-chain volume is not cavity volume (the pocket may be RESHAPED rather than voided, and water fills the rest) (§59a) |
| 63 | 08-24 | the open/closed free-energy difference is simply out of reach here (§24, §58a) | that is a property of UNBIASED MD. Umbrella sampling returns it, and a coordinate exists: **P88 CA - R116 CA**, picked by scoring every gate-latch CA pair against the two crystals and then checked against 5 x 300 ns -- **closed basin 6.06-6.08 A, open basin 16.3-16.9 A, no overlap**. S9 and S2 sit at the SAME value, so the closed state is geometrically identical with and without ligand -- exactly why stability cannot separate them | 62 windows seeded from real equilibrated frames (not a steered pull), 1.24 us. Run as a **CALIBRATION on a known answer first** (holo must favour CLOSED, apo OPEN); the reported quantity is the holo-apo DIFFERENCE so coordinate error cancels, and overlap/drift/hysteresis print beside every number (§60) |
| 64 | 08-24 | the umbrella restraint atoms (P88 CA / R116 CA) are the same indices in every arm | **K59R adds atoms before residue 88**: the quad arms put them at **1408/1814**, the WT arms at **1403/1819**. 112 hardcoded the WT pair | a hardcoded pair would have restrained **the wrong atoms in every quad window, silently**. Indices are now derived per topology with the residue identity asserted. Also caught: transplanted rotamers appended AFTER residue 191 rather than in residue order (tleap `Atom .R<THR 191>.A<OXT 15> does not have a type`), and a 0-byte prmtop passing a `-f` guard that should have been `-s` (§60d) |
| 65 | 08-24 | the 2x2 factorial would show what the ligand contributes to holding the gate closed | analysed at n=1 on the P88-R116 coordinate: **closed/apo 6.08 A vs closed/+ABA 6.06 A -- a 0.02 A difference** -- and ABA does not close an open gate either (S10 stays 16.3-16.6 A over 2 reps). **ZERO transitions in any cell, either direction, with or without ligand** (hysteretic assignment; a naive midpoint threshold had reported 18 spurious crossings for S10, whose min is 11.01 A and never nears the closed basin at 6 A) | occupancy changes NOTHING measurable about conformation on this timescale. The §58a question is not unanswered by unbiased MD, it is **unanswerable** by it, and the 7 truncated reps would only have added error bars to a quantity carrying no ligand information. Sets §60's acceptance bar: the two closed states are GEOMETRICALLY identical, so any scheme separating them must do it on free energy (§61) |
| 66 | 08-24 | the closed state is identical with and without ligand, so MD sees nothing (§61) | that was true of the MEAN and false of the DISTRIBUTION. The two closed basins differ **threefold in width** (k_eff 3.26 apo vs 1.06 holo) and up to **26x in barrier-ward excursions** (>8 A: 0.19 % apo vs 5.00 % holo). The apo pocket collapses and STIFFENS, matching §29's 'most rigid gate of 15 units' | **MD is not blind -- it measures the basin's CURVATURE and gets ~0.33 kcal/mol of the answer** (harmonic entropy of the softer well). What it cannot measure is the basin's DEPTH relative to the other basin, where most of a switch's ddG lives; depth and curvature are independent. The ratio saturates: 0 of 27,000 frames open in BOTH arms, so both return the same upper bound. Not force field, not entropy -- **ergodicity** (§62) |
| 67 | 08-25 | §53b settled that library SIZE is the wrong axis for round 2 | that was right about what round 2 IS and wrong about what it COULD be. The **exact** smallest library holding a known round-2 sensor for all 11 coumarin ligands is **9,216 with wild-type offered (15x) and 4,608 with one position forced (30x)** against Tian's 138,240; 10/11 costs 2,592 (53x). Its shape: **11 positions but 15 substitutions**, so 8 of 11 need exactly one residue | a >=10x reduction IS available, and it is entirely in **residue identity** -- positions are non-negotiable because every round-2 clone spans 7-8 of them. ⚠ A greedy oracle returned 31,104 (3.4x too pessimistic) AND was `PYTHONHASHSEED`-dependent; both fixed by exact branch and bound with sorted iteration (§63b) |
| 68 | 08-25 | round-1 frequency should be enough to design a secondary library | ligand-blind gets **1 of 11** ligands at Tian's budget. Class-weighting (ECFP4 Tanimoto to the target class) moves positions to **11/11 in the top 12** and the optimal residue to **rank 1 at 7/11 positions** (blind: 4/11), and the best design reaches **8-9 of 11 ligands, 8/11 under leave-one-ligand-out** -- still short of Tian's 11/11 at the same size | the first place in this project where knowing the ligand has PAID. But a **global weight-per-log-size greedy captured 1/11** (it deepens loud positions and shuts quiet ones) and a per-position profile plateaued at **4/11**; only a **coverage/co-occurrence** objective reaches 9/11, and that sits on a knife edge (exponent 1 or 3 -> 2/11). Lead, not method (§63c-d) |
| 69 | 08-25 | sd07 recall measures whether a designed coumarin library is good | sd07's clones were **DRAWN FROM** Tian's library, so it contains them by construction and 11/11 is automatic. A different library of equal quality scores badly because the sensors IT would have found were never screened | the metric measures **rediscovery of Tian's choices**, not library quality -- [[feedback_hits_are_not_optima]] applied to libraries. The symmetric metric is recovery of Tian's 23 substitutions from the same round-1 data: **17-18 of 23** (§63e) |
| 70 | 08-25 | the three substitutions §32 recovered clash and K59R does not | **V81I clashes LESS than K59R (0.12 vs 0.21 A) and was still recovered at rank 150.** The separating axis is VOLUME: F108A -101 A^3, V81I +27, F159L -23, **K59R +5**. K59R is **triply invisible** -- near-isosteric, charge-neutral (Lys +1 -> Arg +1), non-clashing; what it changes is H-bond GEOMETRY | the signature is not hydrophobicity, burial or contact count, it is **does the substitution move volume**. Three separate blindnesses is why §32, §33 and §36 all missed the same residue. The blind spot is **20-28 % of real sensor chemistry** (|dVol| < 25 A^3) (§64a-b) |
| 71 | 08-25 | WT residue volume predicts which direction a position mutates (rho +0.62, p 0.007) | scored against the shrink fraction implied by each position's **own sd04 menu**, mean enrichment is **+0 percentage points** (sd 22). Most of the correlation is the trivial null: a big residue shrinks because most amino acids are smaller. What survives is per-position -- **F108 -54 points, V163 -41, V83 -37** | no structural descriptor reaches significance against position usage on n = 18 (8 tested, best p = 0.12). Residue identity comes from **pooled/class frequency**, not from geometry (§64c) |
| 72 | 08-25 | a position whose substitutions are concentrated can be fixed to save library size | concentration GIVEN a substitution is not universality. **V163: 96 % W when mutated, and 95 % of round-2 sensors carry it -- forceable, worth 2x.** V83L (71 % when mutated) drops capture to 9/11 and F108W (66 %) to 8/11 | the statistic is **two numbers**: P(mutated \| class) AND P(residue \| mutated, class). V163 was **callable in advance** -- pooled round-1 rates it 9 % (bottom of the list) but class-restricted it is 22-25 % mutated and **9 of 9 own-class mutations are W**. All 30x of the forcing headroom is this one position. Pre-registered for PFAS/TNT (§64d-e) |
| 73 | 08-25 | sd03 is the round-1 input for every target class | **the 3,366-compound screening deck contains no PFAS at all**, and Tian ran a SEPARATE PFAS round-1 screen whose data sits inside sd09 under `mut_lib`: **DSM-Hao = 89 clones / 18 ligands (round 1)**, **PFOS_GenWT = 154 / 25 (round 2)**. My first prospective run used sd03, got an empty class set, and "passed" by having nothing to say | that verdict was **void**. Re-run on the right input the §64e rule is a genuine **true negative**: 89 class clones, fired nowhere, and nothing in round 2 was forceable (best carry-top K59M **0.44** vs the 0.80 bar). 7 of 25 round-2 targets had no round-1 hit, again supporting "an initial hit is not necessary" (§65a-b) |
| 74 | 08-25 | P(top residue \| mutated) in round 2 measures selection | **where a library offers exactly ONE substitution at a position it is 1.00 by construction.** Every coumarin position I quoted at 1.00 -- V83, F108, V163, N167 -- offered exactly one residue | the round-2 half of §64d was inflated. What survives and is not circular is **P(mutated)** -- did the winners keep wild-type when they could have. V163 still stands (0.95, and its single option makes it collapsible, 2x); A160 at 0.87 is wild-type-droppable but keeps 4 residues, worth only 1.25x. §63b's branch-and-bound was never affected -- it works from clone sets and forced V163 alone (§65c) |
| 75 | 08-25 | forcing is a general lever worth ~2x per position | **it pays in proportion to how SHALLOW the menu is.** Dropping wild-type at PFAS's three high-P(mutated) positions (E94 0.87, Y120 0.86, K59 0.83) takes 49.5M -> 27.5M = **1.8x only**, because they offer 3/5/8 residues; coumarin's V163 offered ONE and was worth 2x | forcing is worth (1+n)/n, so it is a **shallow-position lever**. Separately, round-1 P(mutated) predicts round-2 P(mutated) at **Spearman +0.74** over 17 positions -- but it missed K59 badly (round-1 rank 7 at 6 %, round-2 83 %) (§65d) |
| 76 | 08-25 | PFAS sensors spread out because the PFAS panel is chemically broader than the coumarins | **refuted**: mean pairwise ECFP4 Tanimoto within class is **coumarin 0.352 vs PFAS 0.339**, indistinguishable | whether a class admits a forceable position is **not** explained by the class's chemical spread, and nothing here predicts it in advance. The rule can say "no" honestly; it has said "yes" once, on one position, retrospectively (§65e) |
| 77 | 08-25 | Tian's libraries are all sized about as tightly as the coumarin one (15x off optimal) | **the PFAS library is over-hedged by three orders of magnitude**: the exact smallest library holding a round-2 sensor for all 25 PFAS is **14,400 (3,441x) with WT offered, 9,600 (5,161x) forced**, against 49,545,216 built. Its optimum uses **8 variable positions, not 13** -- V83, L87, A89, V163, N167 contribute nothing, every sensor keeps wild-type there | library width tracks **designer confidence**, which tracks class-relevant round-1 data: coumarin had 36 own-class clones and a clean position signal, PFAS had 89 whose best position reached only 0.61 and whose identities never converged. **For PFAS the prize is which positions to open (3,441x), not residue identity (1.5x)** -- the opposite of coumarin, and the half class-weighted round-1 already recovers (rho +0.74) (§65f) |
| 78 | 08-25 | the position half of the design problem is only demonstrated in-domain (coumarin) | on PFAS, the **out-of-domain** class, ranking positions by round-1 P(mutated) over 89 clones recovers **7 of 8 needed positions in its top 8 and 8 of 8 in its top 10**, almost monotonically (ranks 1-7 all needed) | strongest position-selection result in the project, and it is on the hard class. But a library from it reaches only **15/25 ligands at 933,120 (53x smaller than Tian's 49.5M)** because residue identity needs top-5 depth to get anywhere. **Positions recoverable, residues not -- the same boundary as coumarin, on different chemistry** (§65g) |
| 79 | 08-25 | with the menu fixed by round-1 data, a scorer only has to RANK combinations -- an easier and pose-free job | **no separation.** 670 variants, protein-only ref2015 with a paired wild-type control repacked in the identical shell: **REAL vs LIBRARY AUC 0.542 (p 0.28), REAL vs WILD AUC 0.443 (p 0.14)**; excluding the clash tail makes it worse (0.333). Real sensors are if anything MORE strained than random library members | a sensor works by **remodelling** the cavity, and remodelling costs packing energy -- the property that makes a combination good is the one ref2015 penalises. Protein-only stability is mildly ANTI-correlated, not merely uninformative. Ruled out: cheap repack scoring as a library filter. Not ruled out: FastRelax (~100x cost, poor prior) and ligand-aware scores (still blocked by the pose problem). ⚠ The first version scored against the RAW wild-type and gave every variant -130 to -165 REU -- input strain being repacked away (§66a) |
| 80 | 08-25 | the START/SRPBCC fold has a useful body of ligand-bound structures to learn from | **it has 11.** Of 261 foldseek homologs, 195 are AlphaFold models and 66 experimental; 11 of those 66 carry a buried pocket ligand, **2 of them ABA** and 2 more buffer additives (HEZ, CXS) -- so **~7 non-ABA biological ligands** | against Tian's 208 ligands with sensor data the structural record adds nothing on the chemistry axis. ⚠ `data/survey_cache/` is the WRONG directory (built by scripts 44-46 for the floppy-ligand survey; 97 % ligand-bound, topped by bacteriochlorophyll and detergent), and any compositional count returns 100 % because glycans and modified residues are HETATM. Census must be geometric (§66b) |
| 81 | 08-25 | the pocket and the HAB1 interface are structurally separate because HAB1 has no pocket residues | **they share the C-terminal helix.** alpha3 (153-180) carries **5 library positions (159, 160, 163, 164, 167) and 6 HAB1-contacting residues (155, 156, 158, 159, 162, 166)** on OPPOSITE FACES, and **F159 is both** -- 3.2 A from the ligand and 3.2 A from HAB1, and the most-mutated position in round 1 (317 of 1747). 4 of 18 library positions touch HAB1: L87, A89, L117, F159 | the ratchet framing survives (HAB1 still contacts no pocket residue) but a mutation at 160/163/164/167 repacks the core of the helix presenting the interface -- **an untested route to dead or constitutive sensors**. ⚠ cpptraj DSSP failed silently (no amide H); argmax over an all-zero row called PYR1 20 strands and ZERO helices (§66c) |
| 82 | 08-25 | §66a settled that protein-only Rosetta cannot rank combinations | **fixed-backbone scoring INVERTED the signal.** Cartesian FastRelax on the identical 670 variants gives **REAL vs LIBRARY AUC 0.284 (p 1.7e-8)** and **REAL vs WILD 0.117 (p 2.1e-23)**, medians REAL -1.1 < LIBRARY 3.3 < WILD 8.6 -- against repack's 0.542/0.443 with REAL apparently worst | a sensor needs backbone motion to accommodate its side chains; denying it charges the sensor for strain it never carries. §66a's "real sensors are more strained" was an artefact. ⚠ Two silent bugs on the way: torsion-space relax moved "frozen" residues **0.98 A** (lever arm; Cartesian gives 0.000), and the AUC label was backwards for two runs (§67a) |
| 83 | 08-25 | the relax score is a black box that has to be run per variant | it is **89 % additive** (R^2 0.892 over 175 substitutions) and its per-substitution coefficients are **orthogonal to round-1 frequency (Spearman +0.03)**. Predicting which of Tian's 23 coumarin substitutions round 2 actually uses: **Rosetta +0.46 (p 0.029)**, class-weighted r1 +0.13, **pooled r1 -0.27 (wrong direction)** | round-1 frequency badly over-weights F159 and V81Y -- **F159I is round-1's most common substitution (66) and appears in ONE round-2 sensor**; V81Y (61) is used 12 times and is the residue relax penalises hardest. The additive fit is also the practical route: score a few hundred combinations, apply the coefficients to a whole library for free (§67b) |
| 84 | 08-25 | if relax ranks residues that well it should improve the library design | **it makes it worse: 0/11 ligands vs class-frequency's 8/11.** Rosetta selects for STABILITY, and the residues that create a new binding site are not the stable ones | **relax is a filter over combinations inside a menu someone else chose, not a menu generator.** Worth **1.74x at 90 % sensor retention** (57.3 % of the library kept), and much better at rejecting non-library residues (31 % of WILD kept). Division of labour: positions from round-1 rate, residues from round-1 class frequency, COMBINATIONS from relax (§67c-e) |
| 85 | 08-25 | the relax filter is a general viability screen | its value is a **DEPTH effect**. Separation widens monotonically with substitution count -- AUC **0.359 (3-5 subs) -> 0.219 (8-10)** -- and Spearman(n_sub, ddG) is **-0.12 in REAL vs +0.14 LIBRARY / +0.29 WILD**: random stacking accumulates strain, real sensors do not. The deepest real sensors are the MOST stable group | stability only becomes limiting at round-2 depth (5-8 substitutions), never at round-1 depth (2-3). That is also why round-1 frequency cannot encode it (§67b) -- it is measured in a regime where the constraint does not bind. Real sensors are mutually COMPENSATING combinations, which is what relax sees and a per-substitution rule cannot (§68b) |
| 86 | 08-25 | §67 validated the relax filter for library design | **it validated it in the FAVOURABLE regime and it is expected to fail on our actual goal.** ref2015 penalises cavities, so it prefers GROW (mean beta **+0.43** vs **+1.75** for SHRINK; Spearman(beta, dVol) -0.25). Coumarins are 11-15 heavy atoms vs ABA's 19, so their sensors GROW the lining (+23.1 A^3 per substitution) -- the same direction ref2015 wants. Across 125 sd03 ligands, Spearman(ligand size, mean dVol) = **-0.30**: <=15 atoms **+2.4**, 16-20 **-2.3**, 21-28 **-9.0**, >=29 **-15.1 A^3** | **for ligands larger than ABA real sensors SHRINK the lining, which is the move ref2015 penalises most.** Left to choose a menu it picks K59D/A160L/N167Y/S122Q -- none in any round-2 sensor. Fix the OBJECTIVE not the method: score **abs(cavity volume - target ligand volume)**, which is still pose-free (volume is a 2D property of the SMILES), and keep ddG as a constraint rather than an objective (§68c-e) |
| 87 | 08-25 | §68d: ddG works on coumarins only because they are SMALLER than ABA, so it should fail on pocket expansion | **refuted.** On PFAS (target 224.6 A^3 vs PYR1's 165.9) ddG gives **AUC 0.334 raw and 0.206 once matched on substitution count** -- BETTER than coumarin's 0.251. The premises were right (mean beta +0.43 grow vs +1.75 shrink; Spearman(ligand size, dVol) -0.30); the inference was wrong | **the MENU already encodes the direction.** Tian's PFAS menu is built from pocket-opening substitutions, so LIBRARY members open the cavity too (+31.6 A^3 vs REAL's +39.2, against WILD's +1.2) and ref2015's grow-bias never gets to express itself. Within a directionally-correct menu ddG ranks **combinatorial compatibility**, not cavity size -- §68b survives, §68c only bites when the score must CHOOSE a menu (§69c) |
| 88 | 08-25 | scoring |cavity volume - target ligand volume| would beat ddG (§68e) | **dead in both arms**: coumarin AUC **0.540**, PFAS **0.461**. It cannot reproduce even the arm where ddG works. Its only significant result is REAL vs WILD on PFAS (0.301), i.e. it distinguishes a pocket-opening library from an arbitrary one | it re-measures the direction the menu has already fixed, so it adds nothing the wet lab needs. ⚠ Also: REAL sensors OVERSHOOT -- coumarin sensors shrink the cavity 35 A^3 below a target only 3 A^3 away, consistent with the ligand not having to fill the pocket (§59a) (§69a) |
| 89 | 08-25 | random draws matched to the real substitution-count range are an adequate null | the PFAS draws were **uniform over 2-6 while the real distribution is skewed high (5.4 vs 3.9 mean)**, and since separation WIDENS with depth (§68b) that **understated** the effect: PFAS ddG AUC 0.334 raw -> **0.206 stratified** | match the null on the confound, or stratify and re-weight; a range match is not a distribution match. Corrected filter payoff is **~2x at 90 % sensor retention in BOTH classes** (coumarin 2.10x, PFAS 2.05x), up from §67d's unstratified 1.74x (§69b, §69d) |
| 90 | 08-25 | the umbrella was 24/62 complete and just needed resubmitting | **38 windows had NEVER RUN.** The four md191 systems were solvated independently (S2 46,570 / S10 51,588 / S9 46,483 / S1 51,585 atoms); 111 seeded the open-basin windows from the OPEN systems while 112 runs every window under the CLOSED topology, so all 38 died instantly with `natom mismatch`. Complete windows span **5.0-10.5 A -- the closed basin only** | **no PMF was ever computable** -- §60's quantity is G(open) - G(closed) and the open basin is at 16.3-16.9 A. My guard in 111 checked holo-vs-apo and missed closed-vs-open WITHIN an arm, because I treated "topology" as which molecules are present when the thing that must match is the ATOM COUNT. Also: they failed in 23 s against a 2 h window and I read the zero as "not started". Fix = adiabatic pulling from the completed 10.5 A window, 0.5 A x 200 ps steps, 3.8 ns/arm (§70) |
| 91 | 08-26 | the coumarin headroom is 15x and PFAS is over-hedged 3,441x (§63b, §65f) | **both were measured against too easy a target.** Those oracles accept ANY known hit, including ones **50x weaker** than the best available for that ligand. Requiring the BEST measured sensor per ligand: coumarin exact minimum goes **9,216 -> 55,296**, i.e. Tian's 138,240 is only **2x** above optimal, not 15x | Tian's library is close to optimally sized; the apparent waste was my metric, not their design. And our best designed library drops from **8/11 ligands to 2/11** on the best-hit target, catching sensors ~4.5x weaker than available. Use the best-hit target from here (§71) |
| 92 | 08-26 | the best-hit correction (§71) deflates the headroom everywhere | **only for coumarin.** Requiring the best measured sensor per ligand: coumarin **15x -> 2x** (Tian sized it near-optimally), but **PFAS 3,441x -> 510x** -- still enormously over-hedged | §65f's claim that library width tracks DESIGNER CONFIDENCE survives the correction. Coumarin had 36 own-class round-1 clones and a clean position signal; PFAS had 89 with a weak signal against a chemically novel class. ⚠ Also: **sd09 encodes potency in a `min_conc (µM)` column from a finer RETEST (9 levels)**, which disagrees with its 3-point primary screen on **98 of 154 rows** (always lower). sd07 has no such column and uses the explicit ladder instead (§71e) |
| 93 | 08-26 | Tian's PFAS library holds 49,545,216 variants (product of its 13 positions) | **the arithmetic is right and the MODEL is wrong -- it is depth-capped at ~6 substitutions, so 2,193,926 (4.4 %).** Tell: **93 of 154 sensors carry exactly 6 substitutions and one carries 7**, while a full-combinatorial member would average 9.2. Coumarin by contrast expects 7.0 and observes 6.8 over a 3-10 range -- genuinely uncapped | headroom overstated **18x (any-hit: 3,441x -> 192x)** and **12x (best-hit: 510x -> 43x)** once both sides are costed under the same cap. Coumarin stands at **2.5x**. Same error as §53's DSM-Hao 3.79e21. **Standing check: before quoting a library size, compare the expected substitution count of a random member against the observed distribution in that library's own hits, and look for a hard ceiling** (§73) |
| 94 | 08-26 | our library-design method can be seeded from a single round-1 hit | **it cannot, for the residue step.** Recovering Tian's coumarin design: full class-weighted evidence puts the optimal residue at rank 1 for **7/11** positions, ligand-blind pooled for **4/11**, and **one own-ligand hit + pooled prior for 4.1/11** over 200 draws -- i.e. a single hit adds essentially NOTHING over ligand-blind | **positions survive on one hit (9.9/11), residues do not.** One clone carries 2-3 substitutions; an 11-position profile cannot be built from it. So a single hit may suffice where the prize is POSITIONS (PFAS, §73) and not where it is IDENTITY (coumarin). Fix = **SSM on the single hit** to make the profile experimentally, exactly as Baker's hcy129 -> hcy129.1 did (§74c) |
| 95 | 08-26 | PFAS has 43x headroom because its library includes residues that turned out useless | **no -- only 4 of 45 offered substitutions are never used (9 %), and for coumarin it is 0 of 23.** Both libraries are well utilised. The driver is that size is EXPONENTIAL in positions: PFAS offers **3.46 substitutions/position where its winners need 1.92**, over 13 positions (ratio 1.53), against coumarin's **2.09 -> 1.64** over 11 (ratio 1.17); 1.53^13 ~ 244 vs 1.17^11 ~ 6 | hedging by ~1.5 residues per position is invisible locally and enormous globally. The lever with the most to gain is therefore **trimming per-position depth**, which is what a residue-ranking step does (§74b) |
| 96 | 08-26 | with a fixed pocket and ~1,150 labelled clones, a model could output a few candidate SEQUENCES per ligand | **the target does not exist.** Sensors for the SAME ligand share only **Jaccard 0.45**, there is exactly **1 identical pair** among all within-ligand pairs, and one ligand's own sensors span **10-16 distinct substitutions over 8-11 positions**. Each ligand has a MANIFOLD of solutions, not an answer | the correct output is a **library, not a shortlist** -- a property of the biology, not of our methods. Constructively, the **core IS predictable**: V163W is in every sensor for **10 of 11** ligands. Per-ligand LOLO, our designed library contains a hit for **9/11 at 69,120**, against an oracle of 32 for one known sensor but **10^3-10^4 to cover the solution union** -- so the realistic remaining prize is **10-50x, not 2,000x** (§75a-b) |
| 97 | 08-26 | the Baker NTF2 result argues that ligand->sequence design is out of reach | **their hard problem is the one we do not have.** They built >10,000 novel backbones and docked into pockets that did not exist: 54,500 oligos, **0.36 %** hit rate, cortisol **1 hit from 630**. We change side chains in a backbone that already folds, binds and transduces, and read out by GROWTH SELECTION | **screening 10^5 PYR1 variants is one flask; 10^5 novel backbones is a campaign.** Library size is therefore NOT our scarce resource (Y2H handles 10^6-10^7), so shrinking 69,120 to 5,000 buys little. ⚠ The real risk is whether a NOVEL chemotype has any solution in the 144-substitution vocabulary at all -- which is what the sealed prospective test probes (§75c-d) |
| 98 | 08-26 | our ligand-present failures could be sampling, scoring or both -- we never separated them | **it is SCORING.** Handed PYR1^MANDI's sequence and the crystal mandipropamid pose, FastRelax converges to the SAME conformation from the crystal rotamers and from a scrambled repack -- **max |A-B| = 0.035 A over 18 pocket positions** -- and that shared minimum sits **0.61 A from the crystal overall, 1.08 A at the four mutated positions** | sampling reaches the score's minimum reliably; the minimum is simply wrong. ⚠ Backbone frame mismatch does NOT explain it (Spearman(bb dev, sc RMSD) = **+0.13**), though K59's own backbone deviates 1.15 A so its number is partly frame. **No ligand-present run before this one tested structural recovery** -- §23j/§36/§37-42/§43-46 all asked for ranking or design (§76a-b) |
| 99 | 08-26 | a ~1 A side-chain RMSD at the mutated positions is a modest modelling error | **it destroys the interaction.** ARG59-ligand contacts measured within one structure, so frame error cancels: crystal **NE-O2 2.64 A + NH1-O2 3.26 A (bidentate)** becomes **NH1-O2 2.78 A, NE-O2 3.96 A (monodentate)** after relax -- two H-bonds under 3.5 A become one | this refines §33b rather than contradicting it: the library never PROPOSES crystal Arg59's non-rotameric chi3, AND the score does not KEEP it when handed it. Same interaction ref2015 overcharged by +10.1 REU in §23j. **An iterative dock/relax/mutate loop converges on ref2015's fixed point, so more iterations cannot help -- only a better score can** (§76c-d) |
| 100 | 08-26 | docking the ligand and relaxing it with the protein should beat the protein-only score | **it is at chance.** REAL vs LIBRARY goes **0.284 -> 0.471** pooled and **0.250 -> 0.504** n_sub-matched; only REAL vs WILD survives (0.117 -> 0.332). 670 variants, 40/40 chunks, zero failures | **the noise is 3-5x the signal**: sd of dG_bind among REAL sensors FOR THE SAME LIGAND is 1.4-5.2 REU (range to 16) against a between-set median difference of ~1 REU. dG_bind also tracks ligand SIZE (Spearman -0.39). Within-ligand AUC recovers to 0.417, still far off protein-only. Both pre-registered causes are demonstrated -- pose noise (one smina run per variant) AND §76's scoring defect -- and this test cannot apportion them. **The protein-only filter works BECAUSE it never touches the ligand** (§77) |
| 101 | 08-26 | relaxing before scoring is harmless bookkeeping | **it destroys the discrimination.** On crystal poses the mandipropamid cross-over is correct at raw (**-1987 REU**) and repack (**-1487**) and WRONG after relax (**+0.75**). WT+mandipropamid starts at **+1994 REU** -- F108 is 0.62 A from a ligand atom -- and relax drives it to **-37.5**, indistinguishable from the real sensor | relaxation lets the wild-type pocket make room for a ligand it cannot accommodate; §37b's "relaxation absorbs the clash it exists to relieve" now shown for the whole cross-over. The ABA arm is correct at all three protocols (§78a) |
| 102 | 08-26 | getting the cross-over right means the score sees complementarity | **it sees one clash.** Single-mutant decomposition: **F108A alone is 99.8 % of the quad's repack signal (-1478.81 of -1481.50)**; F159L is -0.02 after repack, V81I +0.82, and **K59R -3.47 (0.23 %)** | the "success" is a steric collision that §32a already detected geometrically with no energy function, naming F108 (2.78 A) and F159 (1.42 A) unprompted. An unrelaxed score is a clash detector, and we had a cheaper one. K59R stays invisible, though at the right sign and above §69's 1.34 REU replicate noise (§78b) |
| 103 | 08-26 | the cross-over result might be specific to PYR1/mandipropamid | **it reproduces in PYL2.** A 2x2 inside one system (7MWN WIN sensor, 3KDI wild-type + ABA, same numbering) is **correct in all six cells** -- WIN -46.00/-2.44/-1.44, ABA -35.83/-0.67/-2.83 -- and each column is again carried by ONE substitution: **Q64K is 98.6 % of the WIN signal, V166I is all of the ABA signal** | ⚠ **3KDJ is PYL1 + ABI1, not PYL2** (27 % identity, numbered 31-209); using it would have mutated R64/I165/W166 and returned a confident cross-over on the wrong protein. 7MWN's deposited record independently confirms `K64Q, F165A, V166I`. K64 is visible where PYR1's K59 was not **because it CLASHES with WIN** (K59 sits 2.86 A from mandipropamid and does not) (§79a) |
| 104 | 08-26 | ref2015 balances several terms when it discriminates a sensor | **it is ONE term.** Per-term WT-minus-sensor: raw **fa_rep 1995.12 = 100.4 %** of the total, repack **1484.69 = 100.2 %**; fa_atr, fa_sol, fa_elec and hbond_sc contribute nothing and **fa_atr/fa_elec point the WRONG way**. After relax fa_rep falls 1995 -> **1.76** and the residual cancels to **-0.04 REU** | **FastRelax does not overpack, it OVER-RELIEVES** -- WT+mandipropamid goes +1994 -> -37.5, a structure that cannot exist. The score is a **clash detector with a working sign and no usable magnitude**: right when a collision exists, silent when none does, and empty once the structure is physical. Same boundary §32 reached geometrically, now confirmed term by term on 2 receptors and 4 cells (§79b-c) |
| 105 | 08-27 | with real structures in hand, a structural or confidence metric will rank sensor quality | **none does.** 590 Boltz-2 sensor structures over 172 ligands, correlated with measured min_conc: within-ligand rho is **-0.030 (confidence), -0.084 (ligand_iPTM), -0.029 (contacts), +0.089 (H-bonds), +0.061 (burial)** -- chance is 21.5/43 ligands and every descriptor sits on it, with H-bonds and burial pointing the WRONG way. The only strong pooled correlate is **ligand size (-0.246)**, a confound | ⚠ but the task is harder than §69's: these are all WORKING sensors, so this ranks potency AMONG POSITIVES rather than separating positives from random variants, and the label is only 3 levels. Poses are fine (0.49-0.66 A on 4WVO) -- **confidence does not know about binding** (WT+mandipropamid, which cannot bind, scores ligand_iPTM 0.97); seed spread separates 4.9x better (0.35 vs 1.70 A). Mapping recovered from the CIFs themselves: **637 of 718 uniquely resolved** (§80) |
| 106 | 08-27 | with all 62 windows at 20 ns the umbrella would deliver its calibration | **it FAILS 2 of 3 pre-registered criteria.** Holo favours closed (+5.02 kcal/mol) ✅ but **apo also favours closed (+6.92)** ❌ and **ddG = -1.90 is the wrong sign** ❌. The easy explanations are excluded: no window has <5 occupied bins, coordinate drift is 0.02-0.03 A, and **ABA stays bound (52-82 contacts even at 20 A)** | ⚠ **my own §70 repair is the prime suspect**: it replaced every open-basin seed in BOTH arms with adiabatic pulling OUTWARD FROM CLOSED, so the two-directional seeding §60 built in specifically to expose hysteresis is gone, and one-directional pulling would inflate the open free energy in both arms -- the observed pattern. §113's docstring promises a hysteresis check the code never implemented. The dimer mis-specification (apo-open is a DIMER state, we ran a monomer) does NOT rescue it: ddG should still be positive. **No designed-pocket number may be quoted** (§81) |
| 107 | 08-27 | Boltz-2's binder-vs-decoy classifier is the strongest untried idea (EF ~18x in BoltzMol-1) | **it calls 90 % of real PYR1 sensors NON-BINDERS.** Median `affinity_probability_binary` over 637 experimentally confirmed sensors is **0.303**; only **9.6 %** clear p > 0.5. The archive already held the affinity output (2,870 files) so this cost **no GPU** | there IS a signal but at the wrong LEVEL: "called binder" is monotone in potency (20.0 / 11.3 / 5.8 % at 1 / 10 / 100 uM, pooled AUC 0.634) yet **within-ligand rho is -0.028** (20/43, chance 21.5) and `affinity_pred_value` goes **+0.327 pooled -> -0.000 within-ligand**. The signal is BETWEEN ligands, not between sensors. **Boltz-2 varies the LIGAND against a fixed protein; we vary the PROTEIN** -- a direction its training never constrains. Protein-only fa_rep (§69) remains the only method above chance (§82) |
| 108 | 08-27 | a ligand class where Boltz ranks variants would be worth finding, and a corrected permutation test can find it | **calibration yes, ranking no.** Boltz is more confident on **aromatic (rho +0.387, p 4e-7)** and **H-bond-donating (+0.262)** ligands -- drug-like ChEMBL-dense space, NOT non-polar. But for within-ligand RANKING the best subgroup (high-HBD, median rho -0.528, n=9) **passed a multiplicity-corrected permutation at p=0.013 and still failed** | **capsaicin and zucapsaicin have IDENTICAL ECFP4 fingerprints (Tanimoto 1.00) and give rho +0.488 vs +0.183**; capsaicin vs nonivamide differ by **1.056**. Per-ligand rho is sampling noise and the 9 "independent" ligands are ~6 chemotypes, violating the permutation's exchangeability. **The check with teeth is not a p-value but whether near-duplicate inputs agree** -- and the screening library supplies those free (§83) |
| 109 | 08-27 | a hit/non-hit ligand benchmark just needs non-hits sampled from the screen | **a random draw is rigged, and by POLARITY not size.** With random negatives **TPSA alone separates at AUC 0.298** (HBD 0.307, HBA 0.328, rotb 0.364, cLogP 0.632) while heavy atoms manage only 0.412 -- PYR1's hits skew greasy and aromatic. Any model would "succeed" on a descriptor | matched negatives by greedy nearest-neighbour on 7 z-scored descriptors bring the **worst |AUC-0.5| to 0.021**, so the benchmark cannot be won on size or polarity. Both classes co-folded against the SAME wild-type PYR1 (a hit is a ligand some VARIANT bound, so an evolved-pocket/wild-type split would be a fatal asymmetry). 181 vs 181, 362 runs, built not yet run (§84) |

### Bugs caught before they cost anything

| date | bug | how it surfaced |
|---|---|---|
| 08-06 | `/scratch` is node-local → SLURM job died at 00:00:00 with no log | job failed instantly |
| 08-06 | foldseek OOM-killed by the 1 GB login-node cap | bare `Killed` message |
| 08-06 | `VRT` residue has no CA → FastRelax RMSD crash | filtered with `is_protein()` |
| 08-07 | chained `sed` renamed the path before the filename rule matched, so `26_submit_ratchet_v3.sh` called a **nonexistent** `21_rosetta_ratchet_v3.py` | caught while writing `32_`; queued job cancelled and resubmitted **before it ran** |
| 08-10 | `loadamberparams frcmod.ions234lm_126_opc` — TIP3P-era naming, no such file for OPC | tleap exited on all four MD systems |
| 08-17 | protomer B relaxed in isolation drove R116/M158 side chains into protomer A — heavy atoms at **0.36 Å** | inter-chain clash check in `67c`, before any MD was built |
| 08-17 | the closed form was relaxed with an **empty pocket**, collapsing K59 into ABA's site (2.85 → 1.68 Å) | pocket-contact assertion against the crystal in `67c`; backbone checks had passed at 0.000 Å |
| 08-17 | FastRelax seeded only through numpy — identical runs gave −195.4 and −202.4 REU | noticed re-running the build twice; fixed with `-run:constant_seed -run:jran`, now byte-identical |

### Standing methodological lessons

1. **Judge an effect by its size in physical units, not by significance.**
   `latch_bb_rmsd` has F = 8.9 and a true spread of **0.027 Å**.
2. **Compare between-group spread against the SEM, not the replicate sd** (§14a).
3. **A null result is only a result if the assay can discriminate** — always run
   the variance decomposition before concluding "no effect".
4. **Static modelling could not settle the satellite lobe.** Seven crystals
   versus one relax protocol is why the MD baseline exists (§19a).

---

## 23. Session of 2026-08-11 — water staple, benchmark hygiene, preemption

### 23a. The W385 water does sense ABA — but that sensing is not conserved

Leonard et al. (Nat Commun 17, 1234, 2026) anchored their opioid-biosensor poses
by requiring the ligand to hydrogen bond the conserved transduction water. Before
adopting that constraint, it was tested against the deposited mandipropamid
complexes. Superposing on the phosphatase chain B:

| structure | res. | water | W385 NE1 | nearest ligand | Pro88 O | Arg116 N |
|---|---|---|---|---|---|---|
| 3QN1 (ABA) | 1.80 Å | A197 | 3.04 | **2.72** (O10) | 2.71 | 2.92 |
| 4WVO (mandipropamid) | 2.25 Å | A328 | 2.96 | 5.07 (CAY) | 2.52 | 3.21 |
| 8EY0 (mandipropamid) | 2.40 Å | A314 | 3.04 | 5.19 (OAX) | 2.85 | 3.05 |

The water is present in **all three**. Its contacts to HAB1 **W385**, the PYR1
**gate** (Pro88 carbonyl) and the PYR1 **latch** (Arg116 amide) are conserved to
within 0.3 Å. Only the ligand contact varies: ABA engages at 2.72 Å,
mandipropamid sits 5.1–5.2 Å away.

The literature's account of ABA is correct and stands: in the WT structure ABA is
sensed by HAB1 through this water-mediated hydrogen bond. What the mandipropamid
structures add is that the water has **two separable roles**, and only one of them
is conserved.

Its **structural** role — a three-way staple locking gate + latch + HAB1 in the
closed state — is present in all three complexes. Its **ligand-sensing** role is
not: mandipropamid retains the staple intact while never engaging the water.

So ligand sensing through this water is one viable strategy, not a requirement of
the mechanism. **Requiring the contact as a hard filter would exclude
mandipropamid** — the most successful engineered PYR1 agonist there is,
orthogonalised as PYR1\*/HAB1\* in 8EY0 and validated in planta. Use it to rank,
never to exclude.

Two consequences. The staple is load-bearing in a way Arm 2 cannot see, since
Rosetta relax as configured carries no explicit waters — S4_ternary MD is the only
system that contains it. And 8EY0 supplies a rare validated true-negative protein
pair: PYR1\* does not bind wild-type HAB1.

*Wording that must not drift:* the water **does** sense ABA in the WT structure.
The finding is that this function is secondary and demonstrably not conserved —
not that it is "not a sensing water". A Boltz-2 supplementary figure
in Tian et al. suggests TETRA, liothyronine, WIN and mezlocillin may occupy the
water site; worth an orthogonal AF3 check, but co-folding models agreeing tells us
they share inductive biases, not that they are right. Neither resolves ordered
water. PDB ID trap: the orthogonal structure is **8EY0** (digit zero); `8EYO` is
human mitochondrial malic enzyme 3. Mandipropamid is ligand **3UZ**.

### 23b. Sequence recovery is the wrong benchmark metric

Prior work in `mutation_prediction_benchmark/` docked 11 coumarins into WT PYR1
and redesigned 19 library positions with LigandMPNN, scoring top-1 agreement
against the Tian wetlab sequences. Two baselines were missing:

| method | top-1 recovery |
|---|---|
| LigandMPNN (ligand-aware, structure-based) | 33/115 = **28.7%** |
| predict the wild-type residue | 0/115 = 0.0% |
| **ligand-blind oracle** (per position, pooled consensus of the experimental answers) | 91/115 = **79.1%** |

The oracle is built from the test labels and is therefore not deployable — it is a
**ceiling**, showing what any ligand-blind method could reach. It beats the
ligand-aware model by 50 points. Worse, **all 33** LigandMPNN hits fall at
positions where the experiment is ligand-invariant (108: 11/11, 81: 9/11,
167: 7/7, 160: 5/11); at the four positions where reality varies with ligand it
scores 5/43.

Ligand identity therefore moves the answer at only ~21% of position-ligand cells
in this panel, and the model captures none of it. Ranking methods by sequence
recovery would rank them by fidelity to a ligand-independent consensus.

**Proposed replacement — ligand-conditional discrimination.** Score every hit
sequence against every ligand in the panel and ask whether the cognate ligand
ranks first. A ligand-blind method scores 0.5 AUROC by construction, and negatives
come free: a sequence selected for Imperatorin, scored against Psoralen, is
verifiably expressed, folded, and present in the library. This matters because
non-recovery carries no information — 13 of 17 *allowed* F108 substitutions never
came back from Tian's screens (§15), and not every colony was sequenced.

**Two confounds recorded, both real.** (1) Poses were docked into the *wild-type*
pocket and then redesigned, so a wrong pose packs the wrong wall; the 28.7% is
ambiguous between model quality and pose quality. (2) The coumarin sequences come
from a **second-round focused library** built on first-round hits, so they share
ancestry and are not independent draws — some of the 79.1% consensus is phylogeny,
not chemistry. Focused sets must be split as units. The oracle result itself is
untouched by (1), since it never uses the docking.

`scripts/43_ligandmpnn_bias_test.py` tests whether the errors are a conservative
prior (model proposes V163I where experiment demands V163W) by biasing W/F/Y at
0/1/2/3. Primary metric is **overall** recovery, not position 163 — biasing toward
W would make 163 succeed regardless of truth. Secondary split: 13 ligand-invariant
vs 6 ligand-variable positions. Only a pose fix can lift the variable ones.

### 23c. MD moved to `preempt_gpu`, and two latent bugs in `39_md_run.sh`

Tasks 0–3 finished cleanly (S1 ×3, S2 rep0; 300 ns each, ~23–24 h, rc=0). Backbone
RMSD plateaus by ~40 ns at 1.7–2.4 Å; rep1 starts strained at 3.1 Å and relaxes.
Then all three A100 nodes went `DRAINING`, pushing tasks 4–11 to an estimated start
of 2026-08-18.

S4_ternary was moved to `preempt_gpu` (ada6000; `-A preempt`) as job 27386629.
Array tasks 9/10/11 of 27333712 write the same directories and were **held, not
cancelled**, so they keep queue position — `scontrol release` restores them,
`scancel` retires them. Two `pmemd` processes on one `prod.rst7` would corrupt it.

`39_md_run.sh` could not simply be repointed. Both bugs are invisible on a
non-preemptible partition:

1. **Lost segments** — the resume path always wrote `prod_cont.nc` with `-O`, so a
   second preemption overwrites the first resumption's trajectory.
2. **Runaway length** — `prod.in` has a fixed 300 ns `nstlim` and restarts use
   `irest=1`, which *continues* the clock, so resuming ran a further 300 ns.

`42_md_run_preempt.sh` numbers segments `prod_cont_NNN.nc` and computes the
remainder from the restart file's own clock (`ncdump -v time`; a finished replicate
reads **301700 ps**, since the clock includes 1.7 ns of equilibration). It keeps
`prod_backup.rst7` and validates `prod.rst7` before use, bounding loss to ~1 ns.

`scripts/41_md_view.sh` builds stripped, imaged, CA-aligned trajectories in
`data/md/view/` (~102 MB per replicate at 100 ps spacing) for visual QC.

### 23d. Do receptors actually engulf floppy ligands? (`scripts/44_floppy_ligand_survey.py`)

The entropy cost of ordering a 15–22 Å rod inside a channel is not computable by
anything in this pipeline, but the question has an empirical proxy: across the PDB,
do receptors bury large flexible ligands whole, or grip a head and expose a tail?

Per ligand, in one representative complex: buried fraction
(1 − SASA_complex/SASA_free) and an **axial burial profile** — atoms projected onto
their first principal axis, split into three equal-length bins.

Bins are oriented so `head` is the more buried end, which makes the head−tail gap
positive **by construction**. A permutation null (40 shuffles of atoms between bins,
preserving bin sizes and the orientation step) measures that bias directly; only
the **excess** gap over null is evidence of real asymmetry.

Limits recorded: only protein chains are kept, so burial is by protein and not by
cofactors; no resolution filter; crystal-packing neighbours unmodelled; and
selection is biased toward ligands that crystallised at all, which under-samples
exactly the floppy chemistry at issue. That bias runs **toward** finding engulfment,
so a tails-exposed result is the conservative one.

If tails are systematically exposed, the design objective changes from "enclose the
ligand" to "grip a head group and tolerate an exit vector" — a materially different
target than the one §13 currently proposes.

**Result (360 ligands, `results/44_floppy_survey.csv`):**

| heavy atoms | n | buried | head | tail | gap | null | **excess** |
|---|---|---|---|---|---|---|---|
| <20 | 33 | 0.89 | 0.97 | 0.74 | 0.17 | 0.12 | **0.08** |
| 20–30 | 104 | 0.83 | 0.93 | 0.74 | 0.16 | 0.10 | **0.05** |
| 30–40 | 91 | 0.82 | 0.92 | 0.72 | 0.17 | 0.11 | **0.07** |
| 40–55 | 78 | 0.79 | 0.89 | 0.63 | 0.14 | 0.12 | **0.05** |
| >55 | 54 | 0.72 | 0.84 | 0.59 | 0.15 | 0.11 | **0.06** |

`corr(n_heavy, buried) = -0.26`; `corr(n_heavy, excess gap) = **+0.07**`.

**Receptors do engulf large ligands.** Burial declines only modestly across a
threefold size range, and even the largest bin is still 72% buried. Crucially the
head-tail asymmetry *in excess of the null* is small (0.05-0.08) and **flat** — it
does not grow with size, which is the signature head-gripping would have produced.

So the §13 objective stands: enclose the ligand, rather than redesign around an
exit vector. The conformational-entropy objection is **weakened but not
eliminated** — nature evidently does bury large flexible ligands routinely, but
that says nothing about the *cost*, only that it is payable. The sample is biased
toward chemistry that crystallised, and that bias runs toward engulfment, so this
is the conservative reading rather than an optimistic one.

Not yet done: subsetting by elongation to isolate rods specifically (the probe set
is 217/229 rods), which is the population this project actually cares about.

### 23e. LigandMPNN bias test — position-mapping bug, and an inconclusive result

Job 27386657 ran 44 designs (11 ligands x bias 0/1/2/3 on W/F/Y). The first
scoring pass returned a flat **0.9% at every bias level**, which was a bug, not a
finding: LigandMPNN emits the **full 181-residue sequence**, and `consensus()` was
reading characters 0-18 instead of the designed positions. Those characters are
outside the design set, so they never changed. Fixed by mapping residue number to
sequence index from the pose PDB's CA order (position 59 -> index 57).

Corrected, and still flat:

| bias | overall | invariant | variable |
|---|---|---|---|
| 0 | 7/115 = 6.1% | 6/51 = 11.8% | 1/64 = 1.6% |
| 1 | 7/115 = 6.1% | 6/51 = 11.8% | 1/64 = 1.6% |
| 2 | 7/115 = 6.1% | 6/51 = 11.8% | 1/64 = 1.6% |
| 3 | 7/115 = 6.1% | 6/51 = 11.8% | 1/64 = 1.6% |

The bias flag *does* work — predictions move (pos 120 L/F -> F; pos 163 picks up
F) — but they move to the **wrong** bulky residue. Experiment wants W at both 108
and 163; the model gives Y and F and will not give W at any bias level tested.

**This test is inconclusive as designed, and the reason must be recorded.** The
bias-0 arm was supposed to reproduce the repo's 28.7% baseline and instead scores
6.1%, so this configuration (temperature 0.1, `--ligand_mpnn_use_side_chain_context 1`,
`sorted()[0]` pose) is **not** the configuration that produced the earlier numbers.
The comparison *across bias levels* is internally valid — same seed, same pose,
bias the only variable — but it cannot be attributed back to the earlier study.

What can be said: within this configuration, biasing toward bulk does not improve
agreement in any stratum. That is weak evidence against the conservative-prior
hypothesis and weak evidence for the pose being the problem, consistent with the
model preferring the wrong bulky residue rather than being uniformly timid. Before
this becomes a claim, the original run configuration has to be recovered from
`mutation_prediction_benchmark/scripts/` so bias 0 reproduces 28.7%.

**Scheduling arrangement as of 2026-08-11 16:0x.** `preempt_gpu` turned out to be
capped at **one GPU** for this account (`preempt` + `gpu` QOS =
`cpu=48,gres/gpu=1`), so it reaches newer hardware but not more of it — three S4
replicates there would run strictly serially at ~29 h each. ada6000 does run S4
faster than the A100 benchmark (~250 vs 167 ns/day, no preemption in the first
90 min), so the split adopted is **hybrid**:

| replicate | where | state |
|---|---|---|
| S4 rep0 | `preempt_gpu` (gpu10, ada6000), job 27386629_0 | running |
| S4 rep1, rep2 | `gpu`, tasks 27333712_10/11 | released, queued for A100s |
| S4 rep0 duplicate | `gpu`, task 27333712_9 | **held** — preempt owns rep0 |
| S2 rep1/2, S3 ×3 | `gpu`, tasks 4–8 | queued |

Task 9 must stay held for as long as the preempt job owns `S4_ternary/rep0`; two
`pmemd` processes on one `prod.rst7` corrupt it. Note `squeue` collapses an array
into a single row and shows the union of reasons, so a partially-held array reads
as fully held — use `scontrol show job <id>_<task>` to see per-task state.

### 23f. Stage 1 — retrospective recovery of PYR1^MANDI (`scripts/47_stage1_inputs.py`)

Can a sequence designer recover a known engineered receptor from WT plus a
perfectly placed ligand? Ground truth from 4WVO vs 3QN1 over the 174 residues
resolved in both: **K59R, V81I, F108A, F159L** — four substitutions, all inside
Tian's randomised positions.

**Library choice is not free.** All four are reachable only in **DSM-Hao** and
**TSM**; under the Coumarin alphabet F108A is not offered at all (position 108
allows W only), so that arm would be unwinnable by construction.

| library | reachable |
|---|---|
| Coumarin | 1/4 | 
| PFAS | 1/4 |
| TNTv1 / TNTv2 | 2/4 |
| **DSM-Hao / TSM** | **4/4** |

This also reframes §23b: F108W appears in 11/11 coumarin sequences partly because
**W was the only option the library provided**, not purely because it was preferred.

**Inputs** in `data/stage1/` — `wt_mandi.pdb` (true 4WVO pose superposed into the
WT frame, 0.48 Å CA RMSD over 174 CA), `polygly_mandi.pdb` (15 designable
positions truncated to Gly), `wt_aba.pdb`. Designable set is Tian's 18 randomised
positions minus 87/89 (gate) and 117 (latch); none of the four true mutations sits
there, so the exclusion is free. `INCLUDE_GATE_LATCH` runs the alternative arm as a
diagnostic — a method proposing gate mutations is proposing to break transduction.

**Two baselines, both of which constrain the design.**

*Trivial clash baseline.* Ranking WT side chains by steric overlap with the ligand
gives 108 (0.62 Å), 159 (1.98), 59 (2.86), 83 (3.14) — **3 of 4 true positions**,
missing only V81I, which barely clashes (3.28 Å, one atom within 4 Å) and is
evidently a packing optimisation rather than clash relief. So **position
identification is nearly free; identity selection is the real task**, and any
method must beat 3/4 on positions to have demonstrated anything.

*Chance level.* Uniform draws from the DSM-Hao alphabets across 15 positions:

| N sequences | E[true subs recovered] | P(all 4) |
|---|---|---|
| 1 | 0.31 | 0.000 |
| 10 | **2.11** | 0.068 |
| 50 | 3.84 | **0.845** |

Recall-at-N is therefore nearly uninformative at generous N — 50 sequences recover
all four 85% of the time at random.

**But that uniform null is the wrong reference, and this is the load-bearing
methodological point.** Sampled sequences are **not independent** — designers
mode-collapse, so 10 samples may carry 2–3 effective draws, making a null built on
10 independent draws too harsh. Rather than guess a correction factor, the primary
null is the **ligand-swap control**: run the identical protocol on `wt_aba.pdb`,
changing only the ligand. It inherits whatever convergence the method has, so no
independence assumption is needed. Scoring is three numbers per position —
probability mass on the true residue, the same from the ABA run, and their
difference. **Only the difference is evidence.** The uniform column is retained,
relabelled as an upper bound on chance under independence.

Effective sample size (distinct sequences, mean pairwise Hamming distance,
per-position entropy) is reported per run: if 50 sequences carry 3 effective draws,
that is itself a finding about the method.

**Noted for later — dual-library design.** Mirror Tian's two-stage approach: round 1
wide positions with few substitutions per sequence, round 2 narrow positions with
more substitutions each. Fits the probability-based scoring, since round 1 needs
per-position marginals and round 2 needs joint combinations.

### 23g. ~~Stage 1 result — LigandMPNN does not recover PYR1^MANDI~~ (job 27392339) — **RETRACTED 2026-08-12**

> **⚠ RETRACTED. Do not cite any number in this section.** The two mandipropamid
> arms contained **no ligand**. `47_stage1_inputs.py` wrote the ligand's PDB
> records one column left of specification, which put the `3` of `3UZ` into the
> altLoc field; ProDy — LigandMPNN's parser — keeps only altLoc `' '` or `'A'` by
> default and silently discarded all 29 atoms. ABA survived solely because its
> code `A8S` happens to begin with `A`, so **the ligand-swap null was the only
> arm that ever held a ligand**, and `delta = P_mandi − P_aba` measured *apo minus
> ABA-holo*. The corrected run is **§23h**, and it reverses the conclusion.
>
> The section is kept in full rather than deleted. The failure mode — an input
> defect that every downstream step tolerated, yielding a complete and internally
> consistent result — is the most transferable thing stage 1 has produced.

Six arms, 50 sequences each at T=0.2, scored as §23f specifies.

**Primary — probability mass on the true residue, against the ABA ligand-swap null**

| mutation | P(wt_mandi) | P(polygly) | P(wt_aba) | **delta** |
|---|---|---|---|---|
| K59R | 0.165 | 0.165 | 0.060 | **+0.105** |
| V81I | 0.970 | 0.970 | 0.996 | **−0.026** |
| F108A | 0.001 | 0.001 | 0.000 | +0.001 |
| F159L | 0.000 | 0.000 | 0.000 | +0.000 |

Identical under the DSM-Hao and unrestricted alphabets to three decimals.

**This is a clear negative.** Only **K59R** carries ligand-conditional signal
(+0.105). **V81I is negative** — the ABA control prefers it *more* than the
mandipropamid arm, so its apparent recovery is a ligand-independent preference,
exactly what the null was built to expose. And **F108A and F159L are assigned
essentially zero probability** in every arm. F108A is this project's own lead
prediction from cavity volume and novelty analysis (§14d, §15); LigandMPNN gives
it p = 0.001.

**Secondary — recall-at-N, which the null defeats**

| arm | recalled | which |
|---|---|---|
| wt_mandi (both alphabets) | 3/4 | K59R, V81I, F159L |
| polygly_mandi (both) | 3/4 | K59R, V81I, F159L |
| **wt_aba (the null)** | **2/4** | K59R, V81I |

Uniform chance at N=50 is **3.84/4**. So 3/4 is *below* the independent-draw
chance level, and the ligand-swap null reaches 2/4 on its own. Recall-at-N would
have looked like a 75% success; it is not one. This is the second time in this
project that a plausible-looking recovery number dissolved against a proper null
(cf. the ligand-blind oracle, §23b).

**Effective sample size** — 50 sequences give 49 distinct, mean pairwise Hamming
5.51 over 15 positions, mean per-position entropy 0.87 bits. So the samples are not
collapsed, and the concern that motivated the ligand-swap null did not bite here;
the null was still the right reference, since it is what showed V81I to be an
artefact. Note the ABA arm is measurably tighter (Hamming 3.96, entropy 0.61):
LigandMPNN is more certain about wild-type PYR1 with its native ligand, which is
the sane direction.

**A flaw in my arm design, recorded rather than hidden.** `wt_mandi` and
`polygly_mandi` returned *identical* numbers because LigandMPNN masks the side
chains of designable positions regardless — `--ligand_mpnn_use_side_chain_context`
supplies context only from NON-designed residues. Poly-Gly truncation touched only
the 15 designable positions, i.e. exactly the atoms the model already ignores, so
the two inputs are the same from the model's point of view. **The poly-Gly arm is
not a valid separate condition for LigandMPNN.** It remains meaningful for
Rosetta and for docking, where side chains are real, and the poly-Gly structure is
still the right input for the admissibility test of §23. To make it a genuine arm
for a masked model, the truncation would have to extend to non-designed residues
lining the pocket.

**Interpretation.** With a *perfect* pose — the crystallographic ligand placement
from the engineered receptor itself — LigandMPNN recovers one of four substitutions
above its own ligand-swap control. Since stage 2 removes pose knowledge, and stage
1 is the ceiling, the cascade should not proceed to stage 2 on LigandMPNN alone.
Per the gate condition, that is the outcome to act on rather than a reason to stop:
add **Rosetta FastDesign** to stage 1, which fails differently from a learned model
and is the Leonard baseline, before concluding anything about the architecture.

Consistent with §23b: LigandMPNN reproduces natural-looking pockets, and the
mutations engineering actually needs — F108A, F159L, both removing bulk to make
room — are ones it will not propose.

### 23h. The apo bug, and the corrected stage-1 LigandMPNN result (job 27412775) — 2026-08-12

#### What broke

`47_stage1_inputs.py` assembled its outputs by writing protein records with
Biopython's `PDBIO` and then appending the ligand with a hand-rolled f-string.
That f-string placed `resName` in columns 17–19 instead of 18–20, shifting every
field from `altLoc` onward one column left:

| field | columns | should be | was read as |
|---|---|---|---|
| altLoc | 17 | `' '` | **`'3'`** (first char of `3UZ`) |
| resName | 18–20 | `3UZ` | `UZ` |
| chainID | 22 | `A` | `' '` |
| x, y, z | 31–54 | correct | **still correct** |

The coordinates survived because the values are short enough that the displaced
8-character windows still contained them — `'  1.579 '` parses to `1.579`. That
is precisely why nothing caught it. The files opened correctly in PyMOL, the
ligand sat in the pocket, and every script that touched them succeeded.

ProDy — which is what LigandMPNN parses with — defaults to `altloc='A'`, keeping
only records whose altLoc is `' '` or `'A'`. Mandipropamid's `'3'` matched
neither, so ProDy returned a protein-only structure:

```
wt_mandi       total=1414  protein=1414  hetero=0     ← 29 atoms discarded
polygly_mandi  total=1352  protein=1352  hetero=0     ← 29 atoms discarded
wt_aba         total=1433  protein=1414  hetero=19    ← kept
```

Re-parsing with `altloc='all'` recovers all 29, confirming the mechanism.

**ABA survived by coincidence, which is the worst possible outcome.** Its CCD
code `A8S` put a literal `A` into the altLoc column — the one character ProDy
accepts. Had both ligands been dropped, all six arms would have been apo, the
deltas would have been ~0, and the result would have looked broken. Instead the
*null* was the only arm with a ligand, so `delta = P_mandi − P_aba` was
`apo − ABA-holo`: a sign-inverted quantity that still produced a plausible table.

#### Fixes

- `het_line()` writes strict PDB columns and asserts each record's altLoc,
  resName, chainID and coordinate fields after writing.
- `write()` re-parses every file **with ProDy at its default settings** and
  refuses to proceed unless all ligand atoms are visible. Column asserts alone
  would not have sufficed — the file was readable, just readable as something
  else. Only a parse with the consuming library catches that.
- `48_stage1_run.py` runs `preflight()` before any compute, checking ligand atom
  **counts** per arm (29/29, 29/29, 19/19) and exiting otherwise.
- Neither fix uses `altloc='all'`. That would restore the atoms while leaving
  the malformed file in place, and hide the next occurrence.

**Rule for the project:** *validate inputs by parsing them with the library that
will consume them, and assert on counts.* Visual inspection and shape checks both
pass here. Deleted results are in `results/stage1_RETRACTED_apo_bug/`.

#### The corrected result — a partial recovery, not a failure

Same six arms, N=50, T=0.2, seed 37; only the inputs changed.

| mutation | P(wt_mandi) | P(polygly) | P(wt_aba) | **delta** | retracted delta |
|---|---|---|---|---|---|
| K59R | 0.000 | 0.000 | 0.060 | **−0.060** | +0.105 |
| V81I | 0.155 | 0.155 | 0.996 | **−0.841** | −0.026 |
| F108A | 0.111 | 0.111 | 0.000 | **+0.111** | +0.001 |
| F159L | 0.071 | 0.071 | 0.000 | **+0.071** | +0.000 |

Identical to three decimals under both alphabets — the DSM-Hao restriction never
binds, because everything the model wants at these positions is already in the
library.

**The conclusion inverts on the two positions that matter most.** With the ligand
actually present, LigandMPNN finds **F108A (+0.111)** and **F159L (+0.071)**,
both with *zero* mass in the ABA arm — unambiguously ligand-conditional. The
retracted section claimed these were the model's blind spot and read that as
confirmation of §23b. That reading was an artefact: an apo pocket has no reason
to open itself up, so of course truncation carried no weight. **F108A is this
project's own lead prediction (§14d, §15), and the method does propose it.**

**V81I is now strongly negative (−0.841), and this is the null doing its job.**
With ABA the model puts I at position 81 with probability 0.996 — a mutation away
from wild type, in the arm whose correct answer is *zero* mutations. Mandipropamid
*suppresses* it to 0.155. So V81I is not merely ligand-independent; it is
anti-correlated with the ligand that actually requires it. Any protocol scoring
raw recovery would bank V81I as a hit in both arms.

**K59R is missed (−0.060).** The model gives R at 59 zero mass with
mandipropamid and 0.060 with ABA.

**What the two recovered mutations have in common.** From the trivial baseline in
§23f, ranked by steric overlap with the ligand: F108 (min dist **0.62 Å**, 23
atoms within 4 Å) and F159 (**1.98 Å**, 12 atoms). K59 (2.86 Å) and V81 (3.28 Å)
barely touch it. **LigandMPNN recovers exactly the two positions where the clash
is severe, and misses both where it is not.** The clash baseline already
identifies 3 of 4 *positions* for free. So the honest statement is: the model
supplies ligand-conditional *identities* at positions that steric overlap alone
would have flagged, and adds nothing where the required change is not driven by
overt clash. That is a real capability and a real limit, and it is a much
narrower claim than either §23g or a naive reading of "2/4 recovered".

**Recall-at-N is unchanged as a metric and still uninformative.** Mandipropamid
arms 3/4 (V81I, F108A, F159L), ABA null 2/4 (K59R, V81I), uniform chance 3.84/4.
Note the mandi arms "recall" V81I whose delta is −0.841 — recall-at-N counts the
single most anti-ligand proposal in the study as a success. §23b's retirement of
sequence recovery stands, and is now doubly supported.

**Effective sample size.** 38/50 distinct, mean Hamming 3.25, entropy 0.54 bits —
*more* collapsed than the apo run (49/50, 5.51, 0.87), as expected once a ligand
constrains the pocket. The ABA arm stays looser (41/50, 3.96, 0.61).

**The poly-Gly finding from §23g survives.** `wt_mandi` and `polygly_mandi` remain
identical to three decimals, confirming this is a property of how LigandMPNN masks
designable side chains and not a consequence of the bug.

#### Effect on the gate

§23f's gate said a LigandMPNN failure requires a physics-based second method
before any claim about the architecture. That still holds, but the reason has
changed: this is now a **partial success**, and the open question is whether
Rosetta finds K59R and V81I — the two positions where clash is weak and an
explicit energy function, with a charged ABA carboxylate at K59, might see what a
learned model does not. Built as §23i.

### 23i. Stage 1, Rosetta FastDesign arm — built and running (job 27421344) — 2026-08-12

**Status at time of writing: 300 trajectories in flight, results not yet merged.**
Everything below the "how to finish" heading is protocol, not outcome.

> **Completed 2026-08-13.** Results are in **§23j**, together with the null-arm
> failure that constrains how much they can be asked to carry. The pre-registered
> table at the end of this section is applied there verbatim, including the part
> that goes against the arm.

#### Why this arm exists

§23f's gate: one method's result is not a statement about the architecture. With
§23h showing LigandMPNN recovers F108A and F159L but misses K59R and V81I, the
sharp question is whether a physics-based method finds the two it missed. Those
are the two positions where steric clash is *weak* (K59 2.86 Å, V81 3.28 Å from
the ligand, versus F108 0.62 Å and F159 1.98 Å) — exactly where a learned model
has least to go on and an explicit energy function has most. K59 in particular
salt-bridges the ABA carboxylate, which is an electrostatic fact Rosetta scores
directly. Rosetta is also Leonard et al.'s baseline.

#### Pipeline

| script | does |
|---|---|
| `49_stage1_ligand_params.py` | crystal coordinates + CCD bond orders + explicit H + aromatic perception → `molfile_to_params.py` → `.params` |
| `50_stage1_rosetta.py` | one arm-block: FastRelax with a design-enabled task factory |
| `50b_submit_stage1_rosetta.sh` | 30 array tasks × 10 trajectories = 6 arms × 50 |
| `51_stage1_rosetta_merge.py` | ligand-swap delta, recall, effective sample size, proposed substitutions, interface energies |

Ligand topology cannot come from a PDB — Rosetta needs bond orders — and it must
not come from the CCD's *ideal* conformer, because the **bound** conformer is the
object of interest (§23d, scripts 45/46). The two are married: coordinates from
the stage-1 file, chemistry from the CCD SMILES, verified one-to-one at
**0.0000 Å** for both ligands.

**ABA is modelled as the carboxylATE, not the CCD's neutral acid.** At assay pH
it is charged, and in 3QN1 that charge salt-bridges K59 — the very position under
test. `--aba-neutral` builds the acid for the sensitivity check, which is worth
running before any claim about K59.

#### Protocol choices, and what they cost

- **Design** at the 15 non-gate/latch Tian positions; WT identity always allowed,
  so a position whose WT residue is outside the library alphabet is not forced to
  mutate.
- **Repack** an 8 Å all-heavy-atom shell around the ligand (43 residues);
  **freeze** the remaining 121. This is cheaper *and* better: side chains 30 Å
  away flipping between trajectories inject variance into sequence frequencies
  that has nothing to do with the ligand. Both arms use the identical rule, so
  the ligand-swap null is unaffected.
- **Backbone and ligand jump fixed.** Stage 1's premise is a perfect pose. If the
  ligand could translate, a failure could be blamed on drift and a success could
  come from relocating the ligand somewhere WT already accommodates.
- **Chi minimisation only where the packer acts**, so minimisation cannot quietly
  relax the 121 frozen residues.
- **Poly-Gly is a genuine arm here**, unlike for LigandMPNN (§23g): Rosetta packs
  explicit rotamers, so truncation really does remove the clash signal.

**Limit, stated plainly:** Rosetta optimises a fixed backbone against a fixed
pose, so it can only relieve the clash it is handed — the same perfect-pose
assumption LigandMPNN got. That is what keeps the comparison fair; neither method
is being asked stage 2's harder question.

#### Three sizing errors, recorded because each was silent

1. **Whole-protein repacking** did not finish one trajectory in 28 min under
   `-ex1 -ex2`. Fixed by the pack shell.
2. **`NeighborhoodResidueSelector` at 8 Å selected 3 residues.** It measures
   between *neighbour atoms*, and `molfile_to_params` gives a 29-atom ligand a
   single NBR atom, so the sphere was drawn from one point on a 12 Å molecule.
   The pocket froze solid, every trajectory returned an identical sequence, and
   the interface energy came out a tidy **−12.29**. With an all-heavy-atom shell:
   43 residues, **−19.11**, and trajectories that actually differ. A protocol bug
   that produces *more* consistent numbers is the dangerous kind.
   I first misdiagnosed this as a seeding failure; both PyRosetta reseed calls
   were verified reproducible, and the real cause was the frozen environment.
3. **All 30 tasks wrote one shared input path**, leaving a 16 MB file. Every task
   still loaded a valid pose — the heavy-atom-count and one-to-one pose asserts
   passed in all 30 logs — but that was timing luck. Input paths are now
   per-unit.

#### How to finish this after the session ends

```bash
cd /bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion
squeue -j 27421344 -h -r -o "%T" | sort | uniq -c        # expect 30 COMPLETED
python scripts/51_stage1_rosetta_merge.py                # writes the tables
```
Results land in `results/stage1_rosetta/<arm>__b<block>.json`, written
incrementally after every trajectory, so a killed task still contributes.
Six arms must be present: `{wt_mandi, polygly_mandi, wt_aba} × {dsm_hao, free}`.

**Read the output in this order, and do not skip to recall:**

1. `delta = freq_mandi(true) − freq_aba(true)`. Only delta is evidence.
2. **K59R and V81I specifically** — the LigandMPNN misses. These decide whether
   the two methods fail the same way (a claim about the problem) or differently
   (a claim about the methods).
3. Effective sample size. If distinct-sequence count is low, the frequencies are
   worth less than N suggests.
4. Proposed substitutions at *all* 15 positions, including the 11 with no ground
   truth. This is the project's actual deliverable — candidate library positions.
   Discount anything appearing at the same rate in the ABA arm.
5. Recall-at-N last, and only as an upper bound. It counted V81I — delta −0.841,
   the most anti-ligand proposal in the study — as a success in §23h.

**Pre-registered reading of the result**, so it is not chosen after the fact:

| Rosetta finds | means |
|---|---|
| K59R and/or V81I | methods are complementary; stage 1 half-works; **proceed to stage 2** with both |
| only F108A/F159L | both methods recover only what the clash baseline gives free; stage 1 does not clear |
| nothing above null | physics and learning agree the pose alone is insufficient; stop and reconsider the architecture |

Then run `--aba-neutral` as the protonation sensitivity check before writing any
K59 conclusion.

---

### 23j. Stage 1 Rosetta result: F108A only, and a null arm that fails its own control (jobs 27421344, 27438220, 27438909) — 2026-08-13

All 30 array tasks of job 27421344 completed: 6 arms × 50 trajectories.

#### The primary table

`delta = freq_mandi(true) − freq_aba(true)`. Identical under both alphabets, so
only `dsm_hao` is shown; `free` differs by ≤0.02 on every row.

| mut | f(wt_mandi) | f(wt_aba) | delta | reading |
|---|---|---|---|---|
| **F108A** | 1.000 | 0.000 | **+1.000** | recovered, ligand-conditional |
| F159L | 0.000 | 0.000 | 0.000 | position found, identity wrong — proposes **Asp/Glu** |
| K59R | 0.000 | 0.000 | 0.000 | uninterpretable — see the null-arm failure below |
| V81I | 0.000 | 0.980 | **−0.980** | **inverted** |

**V81I inverting replicates across methods.** LigandMPNN gave −0.841 (§23h),
Rosetta −0.980. Two unrelated methods independently prefer Ile at 81 *in the arm
whose correct answer is zero mutations*. This is no longer attributable to one
scoring function, and it is a second, independent strike against recall-at-N,
which scores V81I as a success in both arms.

**F159 is a near-miss of a specific and informative kind.** Rosetta puts position
159 under 100% ligand-conditional pressure — the ABA arm holds Phe at 100%, the
mandipropamid arm mutates it in 100% of trajectories — and then proposes Asp
(dsm_hao) or Glu (free), i.e. burying a carboxylate against a chlorophenyl ring.
The position signal is right and the chemistry is not.

#### The null arm fails its own WT-recovery control

The ABA arm is WT protein with its native ligand. Its correct answer is **zero
mutations at all 15 designable positions** — it is the WT-recovery control, and
nobody had scored it as one. It keeps WT at **4 of 15 (29%)**:

```
 59 wt=K  keeps-WT   0%   N82%, I18%   <- deletes the ABA carboxylate salt bridge
 81 wt=V  keeps-WT   0%   I100%
 83 wt=V  keeps-WT   0%   L100%
 92 wt=S  keeps-WT   0%   A100%
141 wt=E  keeps-WT   0%   K74%, Q18%
108 wt=F  keeps-WT  92%   <- quiet
159 wt=F  keeps-WT 100%   <- quiet
```

The consequence is precise and it is not recoverable by argument: **at any
position the null also mutates, "no ligand-conditional signal" cannot be
distinguished from "the protocol cannot hold a native contact."** K59 and V81 are
both such positions. F108 and F159 are readable only because the null happens to
be quiet at exactly those two — which are, again, the two severe-clash positions.

This is the third instance in this project of the same failure shape, and the
general form is now recorded: **score the control's own correct answer, not only
the contrast.** A null used solely as a subtraction can be broken indefinitely
without showing it.

#### Ruling out the cheap explanations

Two candidate causes were checked and eliminated before blaming the scorefunction.

**Pose.** K59 NZ sits **2.85 Å** from ABA carboxylate O4 and **3.02 Å** from O3,
measured on the pose Rosetta actually scores. That is a textbook bidentate salt
bridge. The geometry is not the problem.

**Ionisation.** `A8S_anion.params` partial charges sum to **−0.970** over 38
atoms — the intended −1, to within the 2-decimal rounding of the charge column.
`3UZ` sums to +0.100 over 51 atoms, correct for a neutral molecule.

> An intermediate reading of "exactly 0.000" for both, which briefly looked like a
> dropped formal charge, was an off-by-one in the reader: field 4 of a params
> `ATOM` line is the MM type — the literal `X` for every ligand atom — and the
> charge is field 5. Summing field 4 returns 0.000 for *any* ligand, which is
> precisely the shape of a real finding. `49_stage1_ligand_params.py` now
> validates the generated params with the correct field, against a tolerance that
> scales as `0.005·n_atoms + 0.02` because per-atom rounding accumulates. The
> regenerated params are **byte-identical** to the ones job 27421344 used, so no
> result depends on this.

#### The favor-native calibration, and why it failed (job 27438220)

`ref2015`'s reference energies are fit for soluble monomer design, not for holding
a native complex, so `50_stage1_rosetta.py` gained `--favor-native`
(`FavorNativeResidue` — which lives in `protocols.protein_interface_design`, *not*
`protocols.simple_moves`, in PyRosetta 2026.06). The weight was swept **on the
null arm only**, against a criterion fixed in the submit script before any output
existed: smallest weight keeping WT ≥70% and K59 ≥50% in both alphabets,
disqualified if it freezes the packer.

| weight | WT frac | WT pos | K59 Lys | distinct seqs (of 10) |
|---|---|---|---|---|
| 0 | 29% | 4/15 | **0%** | 16 (of 50) |
| 0.25 | 45% | 7/15 | **0%** | 7 |
| 0.5 | 59% | 10/15 | **0%** | 10 |
| 1.0 | 65% | 10/15 | **0%** | 7 |
| 1.5 | 85% | 13/15 | **0%** | 3 |

Global WT recovery responds smoothly to the bonus, so the term works. **K59 stays
at 0% at every weight**, including 1.5, where 13 of 15 positions are otherwise
held and the packer is approaching frozen (3 distinct sequences of 10). No weight
passes. Per the pre-registration, the sweep was **not** widened.

#### Why no weight could have worked (job 27438909)

Forcing each candidate at position 59, repacking only the 8 Å shell, scoring with
`ref2015` — per-residue decomposition, relative to WT Lys:

| aa | Δ pose | fa_elec | fa_sol | fa_atr | hbond_sc |
|---|---|---|---|---|---|
| **K** | 0.000 | **−4.019** | **+10.118** | −9.162 | **−0.532** |
| N | −2.939 | −1.660 | +6.742 | −7.406 | 0 |
| **I** | **−4.211** | −0.505 | **+3.524** | −8.425 | 0 |
| R | **+2.274** | −2.114 | +8.925 | −9.718 | 0 |
| Q | +0.701 | −1.748 | +7.061 | −8.281 | 0 |
| M | +0.515 | −0.734 | +4.741 | −8.657 | 0 |

**Rosetta sees the salt bridge and rewards it correctly.** Lys has the best
`fa_elec` of any candidate by 1.9 REU and is the only one earning `hbond_sc`. It
then pays a **+10.1 REU Lazaridis-Karplus desolvation penalty**, 6.6 REU more than
Ile, for burying the ammonium. The electrostatic reward is real and is
outweighed by more than two to one. This is the documented ref2015 buried-salt-
bridge pathology, here quantified on our own system.

Two things follow.

**The sweep's failure was predictable from this number.** Ile beats Lys by 4.211
REU. A favor-native bonus must exceed that to hold K59, and 1.5 already left the
packer at 3 distinct sequences of 10. The weight needed to rescue the salt bridge
is deep inside the regime where the benchmark is frozen and therefore vacuous.
The two independent jobs agree quantitatively.

**K59R is doubly unreachable.** Arg scores **+2.274 worse than Lys**, i.e. 6.5 REU
worse than the best option. Even a protocol that held Lys at 59 would never
*propose* Arg. The Rosetta arm cannot recover K59R for reasons that have nothing
to do with mandipropamid.

#### Applying the pre-registered table

Rosetta found F108A, and F159's position but not its identity. It did not find
K59R or V81I. §23i's table, applied verbatim:

> | only F108A/F159L | both methods recover only what the clash baseline gives free; **stage 1 does not clear** |

**Stage 1 does not clear.** Both methods recover exactly the severe-clash
positions (F108 0.62 Å, F159 1.98 Å) and neither recovers the weak-clash ones
(K59 2.86 Å, V81 3.28 Å). Ranking WT side chains by steric overlap already
supplies 3 of the 4 positions at no computational cost, and neither method beats
that baseline on *positions*. What LigandMPNN adds is ligand-conditional
*identities* where clash already flags the position; what Rosetta adds is one
such identity, F108A, at 100%.

The `--aba-neutral` sensitivity check §23i called for is **moot** and was not run:
neutral ABA weakens the very salt bridge whose burial cost is the problem, so it
can only deepen the effect. It is not a test that can change this reading.

#### What this licenses, and what it does not

The honest scope of the negative result:

- It is a statement about **stage 1's gate**, not about Rosetta generally. The
  arm's failures are concentrated at a *charged* contact and, per the
  decomposition, are a property of `ref2015`'s solvation model.
- It is therefore **directly actionable for the deliverable**: for library design
  around charged ligands or charged pocket contacts, this protocol's proposals
  should not be trusted. Positions 92, 122 and 141 — which the null mutates at
  100% — are exactly the ones the delta correctly nulls out, which is the machinery
  working as designed.
- It does **not** license widening the sweep, upweighting `fa_elec`, or
  constraining the salt bridge. Each would be tuning the control until it agrees;
  the constraint version would additionally be telling the method the answer.

#### Files

| path | holds |
|---|---|
| `results/stage1_rosetta/<arm>__b<block>.json` | 30 files, 6 arms × 50 trajectories |
| `results/stage1_rosetta/stage1_rosetta_summary.json` | merged |
| `results/stage1_rosetta/k59_decomposition.json` | the per-residue table above |
| `results/stage1_rosetta_fnsweep/w{0.25,0.5,1.0,1.5}/` | the calibration sweep |
| `scripts/52…55b` | sweep, picker, six-arm re-runner (unused), decomposition |

`54_submit_stage1_rosetta_fn.sh` is retained but **was never run**: no weight
passed the criterion, so there was nothing to re-run at. It refuses to launch
without an explicit `FN=`.

#### Two defects fixed in `51_stage1_rosetta_merge.py`

Both surfaced by re-running the script rather than by reading its output.

1. **WT identities were inferred** from whichever trajectory happened to be
   unmutated at a position, so a position mutated in *all 50* printed `?` —
   rendering as missing data when it meant the opposite, a position under total
   selection pressure. Now read from the input PDB and cross-checked against the
   ground-truth table by `assert`.
2. **The summary JSON was written into the directory the script globs**, so the
   second run ingested its own output and died on a missing `arm` key.

The merge now also prints `mut%`/`aba%` per position and a position-level
hypergeometric test. That test is what exposed the null-arm problem: 3 of 15
positions are ligand-conditional, 2 of the 4 true ones fall inside, p = 0.145.

## 24. Gate/latch loop dynamics: both states are kinetically trapped at 300 ns

`scripts/58_loop_dynamics_prep.py` → `58b_loop_dynamics_run.sh` (job 27444932/27444988,
cpptraj) → `59_loop_dynamics_aggregate.py` → `60_loop_dynamics_figure.py`.
Figure: `figures/loop_dynamics.png`. Summary: `data/loop_dynamics/summary.json`.

Computed **only** the §19d pre-registered observables, which were fixed before these
runs finished, so nothing here was chosen after seeing the data.

### 24a. The headline

**Open stays open, closed stays closed, and neither ever converts.** Across
6 × 300 ns of WT monomer (3 apo-open + 3 holo-closed, all now complete) plus one
300 ns ternary replicate, the gate never crosses the open/closed watershed in a sustained way:

| system | gate S = d(open) − d(closed) | % of frames open-like | gate RMSD to closed |
|---|---|---|---|
| S1 apo open, 3 reps | −3.22, −3.10, −3.63 | 99, 98, 100 % | 6.69, 5.65, 7.75 Å |
| S2 holo closed, 3 reps | +2.48, +3.83, +4.16 | 0, 0, 0 % | 2.36, 1.62, 1.39 Å |
| S4 ternary (+HAB1), 1 rep | **+4.66** | 0 % | **1.13 Å** |

The **replicate-mean gap is 5.57 Å on the gate and 6.27 Å on the latch**, against a
within-replicate spread (sd of S) of only 0.3–1.2 Å. That is a separation of roughly
5σ between states, and it is **insensitive to the equilibration discard** — the gap
is 5.52–5.66 Å at every window from 0 to 150 ns (§59 section 5).

### 24b. What this does and does not license

**It licenses a stability filter.** A single threshold on *gate backbone RMSD to the
closed reference* separates the states with a wide margin:

| threshold | open frames called closed | closed frames called open |
|---|---|---|
| 3.0 Å | 0.003 % | 4.80 % |
| **3.5 Å** | **0.077 %** | **1.11 %** |
| 4.0 Å | 0.93 % | 0.19 % |

Replicate means never come close to overlapping (open 5.65–7.75 Å, closed 1.13–2.36 Å).

**It does not license a switchability filter, and this is the important caveat.**
Because neither state ever converts, MD on this timescale can only tell you whether a
design *holds* the state you built it in. It cannot tell you whether a design *can
close* — that transition is not sampled even once in 1.8 μs of aggregate WT sampling.
So this filter detects the failure mode "the closed state falls apart" and is blind
to "cannot reach the closed state", which is the more likely failure for an enlarged
pocket. Do not read a stable closed trajectory as evidence of a working switch.

### 24c. Flexibility: a monotonic ladder, and one loop that does not follow it

Backbone RMSF, fit to the average structure:

| loop | open | closed | ternary | open/closed |
|---|---|---|---|---|
| gate 85–89 | 3.05 Å (2.98–3.14) | 1.39 Å (1.28–1.51) | **0.76 Å** | 2.19× |
| Lβ7α5 148–156 | 3.38 Å (2.99–3.68) | 1.60 Å (1.39–1.75) | 1.28 Å | 2.11× |
| latch 115–117 | 1.09 Å (1.01–1.20) | 1.44 Å (1.18–1.91) | 0.86 Å | **0.76×** |
| core (161 res) | 0.80–1.05 Å | 0.71–0.94 Å | 0.75 Å | — |

Three things follow.

1. **Ligand then partner each roughly halve the gate's motion** — 3.05 → 1.39 → 0.76 Å.
   The ordering is clean and monotonic, which is what makes RMSF usable as a graded
   readout rather than a binary one.
2. **Dorosh 2013 (§17a) is confirmed, and then some.** Lβ7α5 is not merely
   ABA-sensitive, it is *more mobile than the gate itself* in the open state
   (3.38 vs 3.05 Å). It was included only because Dorosh reported it; it turns out to
   be the single most dynamic element of the open receptor.
3. **The latch does not discriminate on RMSF and must not be used as a filter.**
   Its ratio inverts (0.76×) and the ranges overlap outright (open 1.01–1.20,
   closed 1.18–1.91). The latch discriminates well on *position* (S gap 6.27 Å) and
   not at all on *flexibility*. Two observables on the same three residues, opposite
   verdicts.

### 24d. The gate–latch "staple" is intermittent without HAB1

Minimum gate–latch heavy-atom distance, and the fraction of frames with no contact
at the 4.5 Å cutoff:

| system | mean min-dist | frames > 4.5 Å |
|---|---|---|
| open rep0/1/2 | 6.02 / 6.64 / 7.17 Å | 62 / 98 / 87 % |
| closed rep0 | 5.12 Å | **60 %** |
| closed rep1 | 4.00 Å | 13 % |
| closed rep2 | 3.78 Å | 3 % |
| **ternary (+HAB1)** | **3.72 Å** | **0.8 %** |

In the closed monomer the staple is **not a persistent contact** — one replicate has it
broken 60 % of the time — and it is the weakest of the pre-registered observables
(replicate gap 0.90 Å, and the only one whose frame distributions visibly overlap,
panel C). Adding HAB1 gives the tightest and narrowest distribution of any system.
This is consistent with the ratchet framing in §23: the staple is not something PYR1
holds shut by itself, it is something the partner protein clamps.

⚠️ The ternary arm is **n = 1** and confounds partner binding with everything else
HAB1 brings (Mn²⁺, a large interface). It is reported beside the open/closed
contrast, never pooled into it. S4 rep1/rep2 are queued.

### 24e. Two confounds and one artefact, stated plainly

- **Conformation is confounded with ligand occupancy.** S1 is apo *and* open; S2 is
  holo *and* closed. There is no apo-closed system in the §19b set, so "the closed
  state is rigid" and "ABA rigidifies" cannot be separated here. This does not damage
  the filtering application — candidates will be judged holo against these same
  references — but it forbids the mechanistic claim.
- **Effective sample sizes are small.** Integrated autocorrelation times are
  3.4–44.9 ns, so each 300 ns replicate carries only **3–35 independent samples**.
  Frame counts are not evidence. Everything inferential here is at replicate level
  (n = 3 vs 3); the frame-pooled AUCs are labelled descriptive and are not the basis
  of any claim.
- **S4's trajectory contains 11.6 ns of duplicated frames.** Its 24 preemption
  segments sum to 31,162 frames where 300 ns at 10 ps/frame gives 30,000, because
  `ntwr` writes the restart every 1 ns while `ntwx` writes a frame every 10 ps, so
  frames after the last restart get re-simulated. ~3.7 % of S4 frames are duplicated
  or superseded, slightly over-weighting the moments before each preemption. Harmless
  at distribution level, wrong for anything time-resolved — deduplicate on the NetCDF
  `time` variable first. Fix for future preempt runs: set `ntwr` = `ntwx`.

### 24f. Numbering: the map is per system, and that mattered

The gate is native 85–89 but **sequential 82–86 in S1/S2 and 83–87 in S4**. S1/S2 come
from script 30, which intersected 3K3K and 3QN1 and dropped residue 2; S4 was built
from 3QN1 alone and keeps it (as ALA, the P2A of §23). A mask of `:85-89` selects the
wrong residues in every system, and S1's mask is wrong for S4.

Two failures caught by asserting rather than assuming:

1. Reusing S1's core mask on S4 gave a **162-residue** fit set against the references'
   **161**. cpptraj set the superposition up anyway and returned **gate RMSD ≈ 84 Å**
   instead of refusing. Now the core is restricted to residues common to all systems
   and asserted identical across them.
2. Addressing ABA as `n_residues + 1` gives 179 in S2 (correct) and 180 in S4 —
   which is **HAB1's first residue**, not the ligand, since A8S is residue 478 there.
   Ligands are now addressed by residue name, `:A8S`.

`58_loop_dynamics_prep.py` rebuilds the map from each system's own `protein.pdb` and
verifies every loop by identity (gate = SER GLY LEU PRO ALA, latch = HIS ARG LEU)
before any frame is read; the map is saved to `data/loop_dynamics/residue_map.json`.

⚠️ Also note an unrelated hazard found here: **an unquoted bash heredoc
command-substitutes backticks**, so a comment reading `` `nofit` `` inside the
generated cpptraj input made bash execute `nofit` and silently delete the word.
Harmless this time; it would not always be.

## 25. Could Tian's round-1 screen have been skipped? Ground truth and headroom

`scripts/62_coumarin_ground_truth.py` → `data/coumarin_benchmark/tian_library_truth.json`.

The question (user, 2026-08-14): Tian mined their primary screen for the coumarin
cluster, derived a profile from the hits, and built a focused 77,327-member library
that yielded sensors for three ligands the original screen missed. Could the sites
have been predicted instead, so round 1 need not be run? It does not have to be
perfect — it has to narrow the window.

### 25a. The ground truth was already on disk

`pnas.2519924122.sd04.xlsx` holds the DESIGN of **every** library — per library, per
position, the exact residues offered — including the focused Coumarin, PFAS, TNTv1
and TNTv2 libraries. ⚠️ Earlier in the session I claimed this file was the DSM/TSM
design only and that the coumarin design was missing. **That was wrong**; the whole
truth set is present. (The column is `Amino Acids allowed for mutation`; an earlier
read printed it as `Amino Acids allowed ` only because the display was clipped at
20 characters — a reminder not to key on a truncated header.)

**The paper's "11 of the 19 binding-pocket sites" is exactly:**

```
[59, 81, 83, 108, 120, 122, 159, 160, 163, 164, 167]      24 substitutions
K59 AQ · V81 IY · V83 L · F108 W · Y120 AM · S122 EGN
F159 HIV · A160 GIMV · V163 W · V164 EFS · N167 D
```

Independently, substructure-matching the coumarin core against the 692 characterised
clones of the primary screen returns **exactly 8 ligands** (Imperatorin, Osthole,
Isopsoralen, Methoxsalen, citropten, Bergapten, Psoralen, Scopoletin), all in
`chem_cluster` 26, over 36 clones — matching the paper's "sensors for eight". Those
36 clones touch **15** positions, and Tian's chosen 11 are a strict subset; the four
dropped (A89, I110, L117, E141) are among the lowest-frequency, and A89/L117 sit in
the gate and latch.

### 25b. Positions are the wrong prediction target — measured, not assumed

The three focused libraries barely differ in *which* positions they randomise:

| | shared | Jaccard |
|---|---|---|
| Coumarin vs PFAS | 9/15 | 0.60 |
| Coumarin vs TNTv2 | 10/15 | 0.67 |
| PFAS vs TNTv2 | 12/15 | 0.80 |

Nine positions are shared by all three — **59, 83, 120, 122, 159, 160, 163, 164,
167** — and exactly one is unique to each library (F108 coumarin, E94 PFAS,
L117 TNTv2).

**So naming the shared nine scores 9/11 = 82 % recall on coumarin using no ligand
information whatsoever.** A position-prediction benchmark is saturated before it
starts. That is the same defect that sank stage 1 (§23j), where a trivial clash
ranking supplied 3 of 4 positions free, and it is why headroom is now measured
before a task is fixed rather than after.

### 25c. The substitutions are almost perfectly ligand-specific

At the nine shared positions:

| pos | WT | Coumarin | PFAS | TNTv2 | 3-way J |
|---|---|---|---|---|---|
| 59 | K | AQ | LMN | DMNRT | 0.00 |
| 83 | V | L | I | I | 0.00 |
| 120 | Y | AM | AGILMNST | F | 0.00 |
| 122 | S | EGN | LRW | E | 0.00 |
| 159 | F | HIV | HIL | M | 0.00 |
| 160 | A | GIMV | ILMNVW | V | 0.14 |
| 163 | V | W | S | GMW | 0.00 |
| 164 | V | EFS | GKLMNSW | DKLN | 0.00 |
| 167 | N | D | G | QV | 0.00 |

**Mean 3-way Jaccard = 0.02.** At 8 of the 9 shared positions the three libraries
share *no allowed residue at all*.

The pocket positions are a fixed lining set Tian re-randomises every time; **all of
the ligand information is in which residue goes there.** The benchmark must predict
substitution menus, not positions.

### 25d. What the window is worth

Combinatorial design-space size (WT plus the allowed residues at each position):

| library | positions | members |
|---|---|---|
| DSM-Hao | 18 | 3.8 × 10²¹ |
| TSM | 18 | 7.7 × 10¹⁵ |
| **Coumarin** | 11 | **1.4 × 10⁵** |
| TNTv2 | 14 | 1.2 × 10⁶ |
| PFAS | 13 | 5.0 × 10⁷ |

Tian's coumarin library as built was 77,327 members. The focusing step is worth
roughly **sixteen orders of magnitude** — it converts an unscreenable space into one
yeast transformation. That ratio, not position recall, is the thing a method has to
earn.

## 26. PRE-REGISTRATION: can a method propose Tian's focused library without round 1?

Written **before** any method is run, so the reading cannot be chosen afterwards.
Ground truth and headroom are in §25; this section fixes the task, the metric, the
baselines and the verdict table.

### 26a. Task

For a ligand class **L ∈ {Coumarin, PFAS, TNTv2}**, output for each pocket position
a **ranked list of amino acids to offer**. Positions are GIVEN, not predicted — §25b
showed they are near-ligand-independent and that naming the shared nine already
scores 82 % recall with no ligand information. Predicting them would measure nothing.

Scored on the **9 positions shared by all three focused libraries**
(59, 83, 120, 122, 159, 160, 163, 164, 167). The class-unique positions (F108
coumarin, E94 PFAS, L117 TNTv2) are reported separately as a bonus, never folded
into the primary score — one position cannot carry a claim.

**Allowed inputs:** WT PYR1 structure, the ligand structures for L, any
general-purpose model, and **Tian screen data for classes other than L**. That last
one is deliberate: a real prospective user would have prior screens on other
chemotypes.

**Forbidden inputs:** any round-1 screen hit for class L, and L's focused library
design. Leakage here voids the result, so the input manifest is written to disk
before the run.

### 26b. Ground truth, and its limitation stated up front

Truth = Tian's focused-library menu for L (`data/coumarin_benchmark/tian_library_truth.json`).

⚠️ **This is what Tian chose to build, not the set of substitutions that work.** It
carries their judgment as well as their data — e.g. A89 was dropped from the coumarin
library despite the same 4/36 clone frequency as Y120 and N167, which were kept.
A method proposing a *better* library would score badly. So a second, function-facing
score is computed alongside: recall against the substitutions actually observed in
the round-2 characterised sensors (sd07/sd08/sd09). Neither is decisive alone;
disagreement between them is itself reportable.

### 26c. Primary metric — library size at fixed recall

Take the top-k residues per position, build the combinatorial library
size = Π(kᵢ + 1), and report the size needed to reach 50 %, 75 % and 100 % recall of
the true menu. This is the metric because it is the thing the user asked for — a
narrowed window — measured in the unit that decides whether a screen is possible.

Reference points, already computed:

| design | members | recall of coumarin truth |
|---|---|---|
| DSM-Hao (do nothing) | 3.8 × 10²¹ | 100 % by construction |
| Tian's Coumarin library | 1.4 × 10⁵ | 100 % |

The gap between those two rows is the entire prize.

### 26d. Baselines that must be beaten

1. **Generic hotspot.** Rank residues at each position by frequency across the 692
   characterised clones **excluding class L**. This is the analogue of stage 1's
   ligand-swap null — and unlike that one it cannot collapse, because it is built
   from real data at every position.
2. **Chemistry-only.** Rank by similarity to WT (BLOSUM) with no ligand term.
3. **DSM-Hao**, as the do-nothing upper bound on size.

### 26e. The discrimination check — a method that ignores the ligand must fail

Compute the **3-way Jaccard of the method's own predicted menus** across the three
classes. Truth is **0.02** (§25c). A method emitting near-identical menus is not
using ligand information, however good its recall looks, because the shared menu
alone can score.

This check is computed from the method's outputs, so it cannot fail silently the way
the ABA null did in §23j.

### 26f. Pre-registered verdict table

| result | means |
|---|---|
| ≥50 % recall at ≤10⁶ members **and** predicted Jaccard ≤0.2 | the method could have replaced round 1 — **apply it prospectively to a new scaffold** |
| ≥50 % recall but predicted Jaccard >0.5 | recall is coming from the generic menu, not the ligand; **not usable**, report and stop |
| <25 % recall at 10⁶ members | does not narrow the window enough; **stop** |
| anything else | report the numbers, make no claim |

10⁶ is set as the ceiling because Tian's own libraries were 7.7 × 10⁴ – 5 × 10⁵
members and a yeast transformation comfortably covers 10⁶.

### 26g. Order of reading

1. the discrimination check — if it fails, nothing else is worth reading
2. library size at 50 % recall, against the generic-hotspot baseline
3. the function-facing score (round-2 sensors), and whether it agrees with the
   design-facing one
4. the class-unique positions, last and as a bonus only

## 27. PLANNED MD: is the loop-dynamics signal actually ligand-dependent?

§24 established that closed PYR1 holds its state and open holds its, with a 5.6 Å
replicate-mean gap. It did **not** establish that this has anything to do with the
ligand — §24e flagged the confound explicitly, because S1 is apo *and* open while S2
is holo *and* closed, and the §19b set contains no apo-closed system.

If the closed state holds equally well with a ligand PYR1 was never built for, then
the §24b stability filter reports on **conformation only** and is blind to ligand
identity — which would disqualify it as a design filter, since every candidate will
be judged holo with a non-native ligand. That is the point of these runs.

### 27a. Design — what to run, and what not to

The proposal on the table was 3 structures × 3 replicates × open and closed = 18
runs. Most of that spends GPU-days on cells that cannot discriminate. The cheaper
design below closes the confound and answers the ligand question:

| arm | system | replicates | status |
|---|---|---|---|
| 1 | closed + ABA | 3 | **have** (S2) |
| 2 | open, apo | 3 | **have** (S1) |
| 3 | **closed, apo** | 3 | NEW — breaks the §24e confound |
| 4 | **closed + non-cognate ligand** | 3 | NEW — the ligand test |
| 5 | open + non-cognate, sanity | 1 | NEW — user's own suggestion |

**7 new replicates ≈ 7 GPU-days, ≈2 days wall at 4 concurrent GPUs** — against ~18
for the full factorial, for the same inferential content.

Arm 3 is the highest value per GPU-day and is worth running even if nothing else is:
it is the missing cell that makes the *existing* 6 replicates interpretable. Without
it, "closed is rigid" and "ABA rigidifies" cannot be separated at all.

Arm 5 tests the user's point that a docked-complex PYR1 should relax into the S1 apo
ensemble. One replicate suffices: it is a check, and §24's own two-reference
projection is the readout. If it does not converge onto S1, that is a finding about
the docked backbone and arm 4 must be re-read in that light.

### 27b. Ligand choice

A **coumarin** (from the eight in §25a), docked into the closed 3QN1 pocket. Reasons:
WT PYR1 demonstrably does not respond to coumarins — that is precisely why Tian
needed 11 mutated positions — they are small enough to place without the steric
clash that disqualifies mandipropamid, and it ties this control to the §25/§26
thread rather than opening a third one.

### 27c. Pre-registered reading — and the asymmetry that must not be over-read

Observables are §19d's, unchanged, so this is directly comparable to §24.

| result | means |
|---|---|
| closed+coumarin behaves like closed+ABA | the filter is **conformation-only** and cannot rank ligands; use it to reject collapse, never to claim binding |
| closed+coumarin drifts toward apo-closed | the filter **is** ligand-sensitive; a real design filter |
| ligand leaves the pocket | strong evidence of poor binding — the most informative outcome |

⚠️ **Absence of release is not evidence of binding.** Residence times are typically
microseconds to milliseconds; 300 ns cannot sample unbinding. Release *if observed*
is informative; not observing it says nothing. This asymmetry is registered here so a
null cannot later be read as a positive.

⚠️ A docked pose is a hypothesis. If arm 4 differs from arm 1, the cause could be the
pose rather than the ligand. Mitigated by running the three replicates from
**different docked poses** rather than three seeds of one pose, so pose sensitivity
is measured rather than assumed.

### 27d. Second run: a pocket-expansion variant

`F108A_R79A_E94A` — already defined in `scripts/variants.py` and already measured in
the pilot:

| variant | cavity Å³ | Δ vs WT |
|---|---|---|
| WT | 174.4 | — |
| F108A | 240.3 | +65.9 |
| R79A_E94A | 255.4 | +81.0 |
| **F108A_R79A_E94A** | **313.7** | **+139.3** (1.8× WT) |
| K59A_F108A_R79A_E94A | 382.4 | +208.0 (2.2× WT) |

This is the variant that removes the gatekeeper **and** both halves of the salt
bridge buttressing it. 3 replicates × 300 ns in the closed state, apo, read on the
§24 observables.

The question is **not** whether the pocket is bigger — that is measured. It is
whether a pocket enlarged by 80 % still holds the closed state, i.e. whether the
ratchet survives losing its buttress. §24b's threshold (gate RMSD-to-closed
3.5–4.0 Å) is the pre-registered filter, and this is its first real test.

⚠️ Recall from the pilot that **expansion does not predict switch cost (r = 0.03)**,
so cavity volume cannot stand in for this measurement — it has to be run.

---

## 28. Rebuilding on the complete 191-residue protein (2026-08-17)

Every MD system up to now inherited whatever the crystal happened to model. That is
the defect this section closes, and closing it also closed a smaller one: an earlier
claim in this project about *what* the crystals are missing was itself wrong.

### 28a. What the crystals actually contain

Audited directly from the CIFs, per chain, with residue identity asserted (12/12
landmark checks pass on each) — `scripts/lib_resnumber.py`:

| file / chain | modelled | internal gaps | broken C–N bonds |
|---|---|---|---|
| 3K3K A (open) | 1–183 | **none** | **none** |
| 3K3K B | 2–184 | **none** | **none** |
| 3QN1 A (closed) | 1–181 | **69–70 only** | E68→F71, **4.64 Å** |
| 3QN1 B (HAB1) | 185–505 | 222–231, 271–282, 462–465 | 11.62, **15.93**, 6.36 Å |

⚠️ **This corrects §27's note**, which recorded a `C(1)–N(3) = 3.02 Å` break in 3QN1
and said 3K3K "had the same gaps". Neither is true. That 3.02 Å break belongs to our
own derived `data/md/S1_apo_open/protein.pdb`, which script 30 produced by
*intersecting* the two crystals — a defect in an intermediate of ours, attributed to
the PDB. And 3QN1 does model residue 2; it is **ALA**, the P2A substitution, which is
a mutation rather than a gap.

Practical consequence: **PYR1 needs one 2-residue loop, not two**, and the open form
needs no internal loop building at all.

The deposition's own `_pdbx_poly_seq_scheme` was used to fix numbering rather than
inferring offsets: 3QN1 chain A entity positions 3–193 → author residues 1–191;
3K3K entity 21–211 → 1–191. Both deposited sequences agree with UniProt O49686 at all
191 positions except residue 2, and re-deriving the sequence this way reproduces
`pyr1_sequence.py`'s canonical string exactly (md5 `5ff74f5d34`) — two independent
routes to the same 191-mer.

### 28b. Four build decisions, each settled by measurement rather than judgement

**1. Alternate conformations are chosen explicitly.** 3K3K carries 0.5-occupancy
alternates at ≥8 positions **including R116, a latch residue** (also S3, S31, S32,
Y58, N70, R95, S109); 3QN1 carries none. Letting each parser's default win would make
the open and closed systems differ by an unrecorded coin-flip at exactly the positions
this project measures. Altloc A is selected explicitly and the count is logged
(67 atom positions in chain A, 48 in chain B).

**2. The 69–70 loop is grafted from 3K3K, not invented.** Superposing the flanking
backbone gives **0.85 Å** — and it is still 0.85 Å using only residues 68 and 71, so
the graft is determined by the flanks rather than by the modeller. The loop also sits
**19.1 Å from the gate and 20.3 Å from ABA**, outside the switching machinery. Real
coordinates from the same protein beat a de novo loop.

**3. Both tails are rebuilt from 181.** 3K3K models two extra C-terminal residues, but
their B-factors are **142 and 151 against a core mean of 49** — 2.9× and 3.1×, present
but barely ordered. Keeping them would give the open form a 2-residue head start built
from real-but-meaningless density while the closed form got 10 modelled residues,
reintroducing the open/closed protocol asymmetry this rebuild exists to remove. Both
forms are truncated at 181 and get an identically generated 182–191.

**4. P2 is restored.** φ(A2) in 3QN1 is **−67.4°**, inside proline's −63 ± 15° window,
so this is a ring closure and not a backbone rebuild.

### 28c. The tail is a sample, not a prediction

182–191 reads **SGDGSGSQVT** — three glycines and three serines, a natural GS-linker,
disordered in every structure ever solved of this protein. No starting conformation
for it is meaningful, so the builder samples 12 seeded conformers and keeps the
best-scoring, rather than dictating an extended chain.

That choice is also what keeps it affordable. A fully extended 10-mer projects ~35 Å
into solvent; the sampled coil conformers cost **+4 % (open) and +9 % (closed)** in
periodic box volume, which is paid on every one of the 300 ns.

**Why model it at all** — the decision was the user's, and the geometry supports it.
The C-terminal Cα is **31.7 Å from the gate and 26.1 Å from ABA**, so a 10-residue
tail cannot realistically reach the §24 observables. But it is **16.3 Å from HAB1** and
**18.2 Å from the partner protomer**, which it reaches easily. If those residues touch
the partner in reality, capping at 181 would delete a real interaction rather than
avoid a nuisance one. The protein in the Y2H assay is full-length; capping simulates a
molecule that does not exist.

⚠️ **Pre-registered**: all core-fit and RMSF masks must EXCLUDE 182–191. The tail is
present so the system is physically honest, not so its motion can be interpreted.

### 28d. The pipeline

| script | role |
|---|---|
| `lib_resnumber.py` | canonical sequence + md5 guard, explicit altloc selection, identity assertions, peptide-bond audit, `verify_build`, `residue_map.json` |
| `lib_structqc.py` | core-preservation, clash, bond-geometry and solvent-cost checks |
| `67_build_pyr1.py` | stage A (crystal → gap-free 1–181, Biopython) → stage B (PyRosetta: restore P2, prepend missing N-terminus, sample the tail, relax only what was rebuilt) → stage C (verify) |
| `67c_assemble_complexes.py` | assemble the open dimer, assert the rebuilt pocket still matches the crystal ABA |
| `68_build_md_systems.sh` | drive `lib_solvate.sh` over the rebuilt structures into `data/md191/` |

**⚠ The trap that cost the most: `FastRelax.set_movemap()` does not restrict
repacking.** It governs minimisation degrees of freedom only. FastRelax repacks through
a **TaskFactory**, and given none it repacks *every residue in the pose* no matter what
the MoveMap says. Both structural defects below came from this, and both walked past a
backbone check reporting 0.000 A. `lib_rosetta.restrict_packing()` now builds a
TaskFactory (`RestrictToRepacking` + `OperateOnResidueSubset(PreventRepackingRLT(), …,
flip=True)`), asserts the resulting task packs exactly the intended set and designs
nothing, and is passed with `relax.set_task_factory()`. Result: **178 of 191 residues
keep their input rotamers** in the open build, 168 of 191 in the closed one.

**Reproducibility.** FastRelax is Monte Carlo. The first version seeded only numpy and
gave **−195.4 and −202.4 REU on two identical runs**; Rosetta's own RNG is now pinned
with `-run:constant_seed -run:jran <seed>`, and two runs then produce **byte-identical
ATOM records**. Seed and conformer choice are written to `build_report.json`.

**Verification, in-pipeline rather than by hand** — every build must pass:

- residue numbering exactly 1–191, sequence identical to canonical (md5 `5ff74f5d34`)
- 13 landmark residues asserted by identity, including the P2 tripwire
- **zero** peptide-bond breaks (measured as C–N > 1.5 Å, not as a numbering gap)
- **zero** held-fixed residues displaced by >0.01 Å — this is what proves the MoveMap
  actually froze the crystal core rather than merely being told to. Observed: max
  **0.000 Å** across all three chains.
- zero steric clashes below 2.2 Å between non-adjacent residues

Built: `pyr1_open_191`, `pyr1_closed_191`, `pyr1_open_191_protomerB`. C–N bonds fall
in 1.308–1.368 Å (ideal ≈1.33).

### 28e. `data/md191/`, not `data/md/`

The old tree holds **13 completed 300 ns replicates** (S1–S3 ×3, S4 ×3, S5 ×1) and is
the evidence behind §24. The two generations are **not comparable** — different protein
(191 aa vs a gapped 178–179), different salt (KCl vs NaCl), different histidine
assignment, and for S5 a different ion ordering. Keep both, label both, never pool.

### 28f. What is NOT done

**HAB1, and therefore S4.** 3QN1 chain B needs three de novo internal loops, the worst
**12 residues spanning 15.93 Å**. There is no local template: the only other HAB1
coordinates in this project (`PYRI_HABI_ABA.pdb`) are a copy of 3QN1 with the same
gaps. All three loops sit **32–50 Å from PYR1, the gate, ABA and the Mn site**, so
their conformation cannot affect any observable — they are being built to remove six
artificial charged termini, not because the loops matter. That makes a de novo build
acceptable here in a way it would not be at an interface.

⚠️ Note also that 3QN1 chain B is a **175–511 construct**, not full-length HAB1 (511 aa
with a 174-residue N-terminal extension). Unlike PYR1, no choice here yields "the whole
protein", so the construct boundary should be **capped** (ACE/NME) rather than
extended — a different decision from PYR1's tail, for a different reason.

### 28g. 3K3K is a mixed dimer, and `S3_apo_dimer` is misnamed

Found while writing the ligand-shell guard: **3K3K is not an apo structure.**

| 3K3K protomer | ligand | gate RMSD to closed (3QN1/A) | latch | state |
|---|---|---|---|---|
| chain A | none | 5.24 Å | 5.08 Å | **open, apo** |
| chain B | **A8S 1001** | **0.90 Å** | **0.68 Å** | **closed, ABA-bound** |

ABA in chain B is lined by the canonical 20-residue pocket — 59, 61, 83, 87, 88, 89,
91, 92, 94, 108, 110, 115 and the rest — the same set that lines it in 3QN1. Measured
by superposing each protomer's core (152 residues, gate and latch excluded) onto 3QN1
chain A.

So the crystal is a **half-occupied dimer**: one apo-open protomer and one holo-closed
one.

⚠️ **`S3_apo_dimer` — 3 × 300 ns, already complete — was built protein-only from both
chains.** Its chain B is therefore a closed, ligand-shaped protomer simulated with an
empty pocket: the same artefact §28d's `ligand_shell()` was written to prevent, except
here it ran for 900 ns. Do not read S3 chain B as an apo-open control, and do not call
the system apo.

**§24 is unaffected.** `scripts/58_loop_dynamics_prep.py`'s `SYSTEMS` map covers S1, S2
and S4 only — S3 was never in the loop-dynamics analysis.

Two consequences for the rebuild:

- **S1_apo_open from 3K3K chain A remains correct** — genuinely apo and genuinely open.
- **A rebuilt dimer must keep chain B's ABA**, and its coordinates must come from 3K3K.
  `data/md/A8S.mol2` carries the **3QN1** pose and is in the wrong frame — the same trap
  recorded for S5 in §23.

There is also a benefit: 3K3K chain B is an experimentally determined closed protomer
that is *not* 3QN1, so the closed state no longer rests on a single crystal.

### 28h. Facts are asserted from the structure now, not recalled

`S3_apo_dimer` was named for a property that is false. That name then acted as a
summary and outranked four separate records saying otherwise — including one that
contained its own refutation, `"S3_apo_dimer (3K3K A+B, ABA removed)"`. The
information was never missing; the label was simply trusted instead of it.

The fix is not to write the fact down again. It is `data/structure_provenance.json`,
which records every source structure **per chain** — modelled range, internal gaps,
peptide breaks, conformational state, and *which ligands the chain actually carries* —
and `lib_resnumber.assert_provenance()`, which script 67 runs before it touches
anything. If a structure disagrees with its record, the build stops.

Verified with a negative control: editing the record to claim 3K3K chain B is apo
raises `ligands in contact with this chain are ['A8S'], expected []`. The mistake that
cost 900 ns can no longer be made silently.

**Two rules that follow.** Never name a thing after a property that could be wrong —
prefer `S3_dimer_3K3K` to `S3_apo_dimer`. And check ligand occupancy and conformational
state **per chain, never per file**: a dimer's protomers need not be occupied alike, and
this one is not.

### 28i. First three systems built on `data/md191/`

| system | source | ligand | atoms | waters | KCl |
|---|---|---|---|---|---|
| `S1_apo_open` | 3K3K chain A | none (genuinely apo) | 51,585 | 12,201 | 33 |
| `S2_holo_closed` | 3QN1 chain A | A8S, crystal pose | 46,570 | 10,934 | 30 |
| `S9_apo_closed` | 3QN1 chain A | none — **by design**, the missing cell of §19b | 46,483 | 10,921 | 30 |

All neutral to <0.001 e. Verified with cpptraj rather than the build log: 191
non-solvent residues each (192 for S2, with A8S).

⚠️ **S3 is deliberately not built yet.** 3K3K chain B is closed and ABA-bound (§28g), so
a protein-only dimer would repeat the first generation's mistake. It needs an ABA in
**3K3K's** frame first — `data/md/A8S.mol2` holds the 3QN1 pose.

**Histidine tautomers differ between open and closed at exactly one position, H60.**
`reduce` assigns HIE in the open protomer and HID in both closed systems; the latch
H115 is HIE everywhere and H127 is HID everywhere. H60 is the dimerisation residue, and
the difference is `reduce` responding to genuinely different local H-bond networks
(3K3K chain A sits in a dimer, 3QN1 chain A does not). It is recorded rather than
forced: overriding a geometry-based assignment to manufacture consistency would be
worse than the confound. Flagged so it is not rediscovered as a surprise in an
open-vs-closed comparison.

---

## 29. All twelve WT trajectories, including the apo-closed protomer (2026-08-18)

§24 analysed three systems. This re-runs the same pre-registered observables over
**all four**, and the addition is the one §24e named as missing.

### 29a. What changed, and why the §24 numbers moved slightly

`S3_apo_dimer` was never analysed — script 58's system table covered S1/S2/S4
only. It is a **mixed** dimer (§28g): chain A apo-open, chain B **closed and
ABA-bound in the crystal**, built protein-only. Its chain B is therefore a closed,
ligand-shaped protomer simulated around an **empty pocket** — the apo-closed cell,
arrived at by accident rather than design.

| | before | now |
|---|---|---|
| trajectories | 7 | **12** (S1×3, S2×3, S3×3, S4×3) |
| analysis units | 7 | **15** — S3 contributes two, one per protomer |
| ternary | n=1 | **n=3** |
| core superposition set | 161 residues | **160** |

The core shrank by one because S3 chain B is modelled from native residue 2 and
cannot supply residue 1. Every unit is refitted on the intersection rather than
letting per-system cores drift apart, so the S1/S2/S4 numbers here are **not
bit-identical** to §24's. They are not materially different either: the gate gap
is 5.47 Å against §24's 5.57, and stays 5.47–5.61 at every discard from 0 to
150 ns.

**The ternary at n=3 confirms the n=1 result.** Gate S = +4.76 / +4.89 / +3.92
(§24 reported +4.66 from one replicate); gate RMSF 0.81 Å (§24: 0.76). ABA stays
in the pocket in all six ligand-bearing replicates (mean 1.21–2.50 Å).

Two new safety checks in `58_loop_dynamics_prep.py`: landmarks are re-verified
against each **prmtop's own `RESIDUE_LABEL` block**, not only the `protein.pdb`
tleap was fed, and the per-protomer prmtop range drives the RMSF mask — `:1-N`
would have reported protomer A's fluctuations as protomer B's, since chain B
occupies prmtop residues **184–366**.

### 29b. Removing the ligand does not open the gate — and does not loosen it

| state | gate S (per replicate) | sd | drift over 300 ns | crossings |
|---|---|---|---|---|
| open, apo (S1 monomer) | −3.21 −3.09 −3.64 | 0.42–1.15 | −0.63 Å | 0–166 |
| open, apo (S3 protomer A) | −3.90 −4.14 −3.85 | ~0.4 | +0.10 Å | 0 |
| closed +ABA (S2 monomer) | +2.42 +3.81 +4.14 | 0.65–0.81 | −0.61 Å | 0–30 |
| **closed, APO (S3 protomer B)** | **+4.76 +4.73 +4.74** | **0.13** | **−0.07 Å** | **0** |
| ternary (S4) | +4.76 +4.89 +3.92 | — | +0.12 Å | 0 |

**Zero crossings of the watershed in 3 × 300 ns**, and the drift is −0.07 Å —
indistinguishable from not moving. The apo-closed protomer is not merely stable;
it has the **tightest** gate in the entire set (sd 0.13 Å against 0.65–0.81 for the
holo monomer) and the tightest gate–latch staple (3.03–3.05 Å, tighter than the
ternary's 3.72–3.88).

Backbone RMSF ladder, gate:

```
open monomer (S1)        3.05
closed +ABA (S2)         1.39
open protomer   (S3 A)   1.07
ternary (S4)             0.81
closed apo protomer (S3 B) 0.72   <- least mobile gate of all
```

### 29c. What this licenses, and what it does not

**The dimer is a large confound and must be stated first.** Dimerisation alone
takes the open gate from 3.05 Å (S1 monomer) to 1.07 Å (S3 protomer A) — a 2.9×
drop with no change of conformation and no ligand. So the apo-closed protomer's
rigidity cannot be attributed to ligand removal; it is measured inside a clamp.

The **controlled** comparison is within the single S3 box, where both protomers
share a thermostat, a barostat, a build protocol and an absence of ligand, and
differ only in conformation: **open 1.07 Å vs closed 0.72 Å, a 1.5× ratio.** The
S1-vs-S2 ratio, where conformation *and* occupancy both differ, is 2.20×. Read
crudely, conformation supplies most of the §24 RMSF ladder and the ligand the
remainder.

1. **The §24b stability filter is conformation-reporting, not ligand-reporting.**
   §27 was designed to test whether the closed state holds equally around a
   non-cognate ligand. This is the stronger version of that test — the closed
   state holds with **no ligand at all**, for 900 ns, more rigidly than with ABA.
   A design that scores well on gate-RMSD-to-closed is being scored for holding a
   conformation, and that number carries **no evidence about occupancy**.
2. **A null here is weak evidence and is reported as such.** Gate opening after
   ligand loss is a barrier crossing; §24 already established that no crossing is
   sampled in 1.8 μs from either basin. "It stayed closed" is consistent both with
   the gate being ligand-independent and with 300 ns being too short. The
   informative outcome would have been a drift toward open; its absence is not the
   converse. See §29d.
3. **S9_apo_closed became more valuable, not less.** It is the apo-closed
   **monomer** — the same cell without the dimer clamp — and it is the only way to
   tell "closed holds without a ligand" from "the dimer holds it". It is built
   (`data/md191/S9_apo_closed`) and unsubmitted.

### 29d. Pre-registered reading, for the ligand-dependence runs

Fixed here, before the enantiomer and non-cognate systems are built, so §27's
observables cannot be reselected once their answers are known:

- **Do not read gate opening as the ligand-dependence signal.** Nothing opens on
  this timescale in any of 15 units. The comparison that *does* respond to ligand
  is the graded **RMSF ladder**, and it must be read against a same-box control.
- **Compare like with like.** A monomer contrast is only interpretable against
  monomers; the 2.9× dimerisation effect measured here is larger than the ligand
  effect being looked for.
- **Switchability still requires biased sampling** — a PMF along S — and no
  unbiased 300 ns run of any composition will supply it.

### 29e. Artefacts

| what | where |
|---|---|
| figure (6 panels, 15 units) | `figures/loop_dynamics.png` |
| full numeric report | `data/loop_dynamics/aggregate_report.txt` |
| machine-readable summary | `data/loop_dynamics/summary.json` |
| per-unit masks and maps | `data/loop_dynamics/residue_map.json` |
| per-unit cpptraj output | `data/loop_dynamics/<unit>/rep<N>/` |
| SLURM | job 27545237, 15 tasks, `cutlerlab`, ~2 min each |

Panels A/B/C/D are §24's, extended to five states; **E** adds the latch state
coordinate in its own panel; **F** plots mean gate S over the first 50 ns after
the discard against the last 50 ns, so "did anything move" is one glance.

### 29f. Indexing audit, and why two runs are both "open + apo"

Asked for directly (2026-08-18), because the latch band in panel D looks
misaligned with the nearest RMSF peak.

**The indexing is correct in all five units, verified three ways.** cpptraj writes
the *topology* residue number in `rmsf_bb_byres.dat`, so protomer B's rows really
do run 184–366; every row maps to a native residue with none left over; and taking
each unit's claimed sequential index straight into the **prmtop's own
`RESIDUE_LABEL`** returns SER-GLY-LEU-PRO-ALA at the gate and HIS-ARG-LEU at the
latch in every case:

| unit | gate native→seq | latch native→seq | latch labels |
|---|---|---|---|
| S1 | 85→82 | 115→112 | HIS ARG LEU |
| S2 | 85→82 | 115→112 | HIS ARG LEU |
| S3 protomer A | 85→85 | 115→115 | HIS ARG LEU |
| S3 protomer B | 85→**267** | 115→**297** | HIS ARG LEU |
| S4 | 85→83 | 115→113 | HIS ARG LEU |

**The apparent misalignment is real and is not an offset.** The local maximum in
that region sits at native **114**, one residue *before* the latch, and the latch
lies on its shoulder — in S2 the profile reads 2.02 (113) → **2.52 (114)** → 1.71
(115) → 1.49 (116) → 1.10 (117). The much larger peak that draws the eye is at
**131–134**, a different loop that is not one of the three pre-registered ones and
carries no shading. Panel D now labels each band with its residue range, tints the
three-residue latch more heavily, and carries a 5-residue minor grid so alignment
can be read off the axis rather than trusted.

**Two units are legitimately "open + apo".** S1 and S3 protomer A are the same
crystal chain — 3K3K chain A — once as an isolated monomer and once inside the
dimer. They are two runs of one state, and the 2.9× difference in gate RMSF
between them (3.05 vs 1.07 Å) is the measured cost of dimerisation, not a
duplicated label. The closed states are the ones that differ in occupancy: **S2 is
closed + ABA, S3 protomer B is closed + APO**, and that contrast is the whole
point of §29b. Legend labels now name conformation and occupancy separately.

**The labels are now asserted, not trusted** (§28h, `feedback` on names that carry
claims). Script 59 checks every unit's mean gate S against the side its name
claims and stops on disagreement, and separately asserts the two S3 protomers sit
on *opposite* sides. They do — protomer A −3.98, protomer B +4.74 — so S3 is a
mixed dimer in the trajectory, not only in the crystal.

---

## 30. The conformation × occupancy factorial (2026-08-18)

§19b set up a 2×2 and only the diagonal was ever run. §29 showed why that matters:
the closed state holds its conformation with no ligand at all, so the §24b
stability filter reports conformation rather than occupancy — but the only
apo-closed data was a protomer inside a dimer, and dimerisation alone changes gate
RMSF by 2.9×.

| | apo | + ABA |
|---|---|---|
| **open** | S1 ✅ 3 × 300 ns (old protocol) | **S10 — never existed** |
| **closed** | S9 built 08-14, never run | S2 ✅ 3 × 300 ns (old protocol) |

Both "right" forms were done and neither "wrong" form had a single step of MD.

### 30a. Why all four are re-run, including the two we already have

The old-tree S9 was built three days after S1/S2 and is **KCl**; S1 and S2 are
**NaCl** (verified from the topologies: S1 Na⁺ 35/Cl⁻ 30, S2 Na⁺ 34/Cl⁻ 28, S9 K⁺
33/Cl⁻ 28). Pairing them would confound ligand removal with a cation swap inside
the one comparison the factorial exists to make — and Na⁺ binds carboxylates more
strongly than K⁺ (§ canonical settings). A factorial whose cells differ in the
build is not a factorial, so all four cells run on `data/md191/`.

The `data/md/` trajectories are **not** discarded: they remain the evidence behind
§24 and §29 on their own protocol, and are never pooled with md191.

### 30b. Building S10 — open backbone, ABA present

`scripts/69_build_holo_open.py`. ABA is **transplanted, not docked**: superposed
from the closed structure by the rigid core (656 backbone atoms over 164 residues,
gate/latch/Lβ7α5 and the modelled 182–191 tail excluded, fit RMSD 1.07 Å). Docking
would answer a weaker question and add pose uncertainty to the one arm whose point
is that the ligand starts in its known correct position.

The open pocket holds ABA **half-formed**, which is the expected result:

- **retained** (11): 59, 83, 92, 94, 108, 110, 115, 120, 141, 159, 167
- **lost** (8): 61, **87, 88, 89** (the gate itself), 91, 117, 163, 164
- **gained** (3): 85, 156, 160

The transplant left one tight contact — ABA O2 against **H115 CE1**, the latch
histidine, at 2.30 Å. Repaired by repacking **only** the two residues within 3.0 Å
(59, 115) with a **PackRotamersMover**, not a FastRelax, so the backbone cannot
move by construction; then a **chi-only** minimisation, because a discrete rotamer
library cannot make a sub-ångström adjustment. Result: closest contact **2.72 Å,
and it is now O4···K59 NZ** — the carboxylate salt bridge itself, at an ideal
distance. ref2015 +79.1 → −36.3 REU.

⚠ **This system's entire value is that its backbone is OPEN**, so any relaxation
toward closed would pre-bias the experiment. Gate and latch backbones are asserted
**bit-identical** to `pyr1_open_191.pdb` (max 0.0000 Å) and all 189 non-repacked
residues are frozen all-atom (max 0.0000 Å).

### 30c. Pre-registered reading — fixed before the runs finish

- **S10 has three outcomes and only two are informative.** The gate **closes** (the
  ligand drives closure and MD can see it), or ABA **leaves** (the open pocket does
  not retain it). *Nothing happening is the uninformative case* — 300 ns against a
  barrier §24 showed is never crossed in 1.8 μs — and must not be written up as
  "ABA does not close the gate".
- **S9 vs S2 is the ligand-removal contrast without the dimer clamp** that
  confounds §29's S3 protomer B.
- **Observables are §19d's, unchanged.** Nothing new is added, so the answer cannot
  be shopped for. Analysis goes through scripts 58–60 with `data/md191` swapped in.

### 30d. Verification and submission

`69b_verify_md191.sh` re-derives each cell's claims from `system.prmtop` and
`system.inpcrd` — never the directory name (§28h) — via a solvent-stripped frame:

| system | protein | A8S | gate→open | gate→closed | verdict | K59 NZ→ABA |
|---|---|---|---|---|---|---|
| S1_apo_open | 191 (1–191) | 0 | 0.00 | 5.52 | open ✓ | — |
| S2_holo_closed | 191 | 1 | 5.52 | 0.00 | closed ✓ | **2.85 Å**, 19 lining |
| S9_apo_closed | 191 | 0 | 5.52 | 0.00 | closed ✓ | — |
| S10_holo_open | 191 | 1 | 0.00 | 5.52 | open ✓ | **2.72 Å**, 14 lining |

The K59 distance is the check that a **wrong-frame ligand** cannot survive: it
would pass residue counts, conformation and clash tests while sitting in bulk
solvent.

**Job 27547457**, 12 tasks (4 systems × 3 replicates), `70_md_run_factorial.sh`.
Split across partitions because `preempt_gpu` caps this account at **1 concurrent
GPU** — tasks 1/3/6/9 (one per system, so every cell gets an early replicate) moved
to `gpu` at 4 concurrent, the rest queued on `preempt_gpu`. ~300 ns/day, so ≈3 days
rather than ≈12.

⚠ Two bugs caught while building, both of the class this project keeps hitting:
mol2 **column 6 is the GAFF2 type** (`ca`, `ho` — lower case), so testing it
against `"H"` classed every hydrogen as heavy; and PDB residue name occupies
**columns 18–20**, so an atom name written from column 13 with no altLoc slot
shifted `A8S` left and Rosetta reported `Unrecognized residue: 8S`. A third was
caught by the verifier itself: cpptraj writes the ligand as `ATOM`, not `HETATM`,
so a hetflag test counted A8S as a 192nd protein residue and measured it against
itself — reporting a 0.00 Å "clash" and a 0.00 Å K59 distance that passed.

### 30e. Reproducing §28–§30

Written after being asked directly (2026-08-18) whether the README carries enough
to reproduce these runs. It did carry the **reasoning**; it did not carry the
**run order or the environments**, and the master list under §19 stops at script
40 while the repo now runs to 70. Both gaps are closed here.

**Environments** — these are not interchangeable. The base miniconda python has a
broken `Bio`, and PyRosetta exists in only one of them.

| tag | interpreter | holds |
|---|---|---|
| **T** | `/bigdata/cutlerlab/jjaco081/conda_envs/tier1_analysis/bin/python` | python 3.11.15, PyRosetta4 `.python311.Release` 2026.06 |
| **D** | `/bigdata/cutlerlab/jjaco081/conda_envs/pyr1_docking/bin/python` | python 3.10.20, numpy 2.2.6, biopython 1.87 |
| **A** | `module load amber/22` | AmberTools, cpptraj V6.4.4, tleap, pdb4amber |
| **G** | `module load amber/22_mpi_cuda` | `pmemd.cuda` (ships **no** tleap — that is why builds load `amber/22` instead) |

**Order**

```
# --- 28: rebuild PYR1 on the complete 191-residue sequence ---------------
T  scripts/67_build_pyr1.py            # seed 20260817, -> data/structures_191/
T  scripts/67c_assemble_complexes.py   # dimer + pocket assertion
A  bash scripts/68_build_md_systems.sh # -> data/md191/  (loads amber/22 itself)

# --- 30: the open + ABA cell --------------------------------------------
T  scripts/69_build_holo_open.py       # seed 20260818, -> pyr1_open_holo_191.pdb
A  bash scripts/68_build_md_systems.sh S10_holo_open
A+D bash scripts/69b_verify_md191.sh   # all four cells; MUST pass before submitting
G  sbatch scripts/70_md_run_factorial.sh   # 12 tasks, job 27547457

# --- 29: loop dynamics over finished trajectories ------------------------
D  scripts/58_loop_dynamics_prep.py    # per-protomer masks, verified vs prmtop
A  sbatch scripts/58b_loop_dynamics_run.sh  # 15 units, cutlerlab, job 27545237
D  scripts/59_loop_dynamics_aggregate.py    # -> summary.json, aggregate_report.txt
D  scripts/60_loop_dynamics_figure.py       # -> figures/loop_dynamics.png
```

To point 58–60 at the new trajectories instead of `data/md/`, change `MD` in 58
and 59 and `MD`/`OUT` in 58b. The masks are re-derived and re-asserted from each
system's own topology, so nothing else needs editing — that is what §29a's
`RESIDUE_LABEL` check is for.

**Determinism.** Both PyRosetta stages are seeded through Rosetta's own RNG
(`-run:constant_seed -run:jran`), not numpy — seeding numpy alone gave −195.4 vs
−202.4 REU on identical runs (§28, change-log row 27). Re-running 67 or 69 with
the recorded seed reproduces **byte-identical ATOM records**. Seeds live in
`data/structures_191/build_report.json` and `s10_build_report.json`. MD seeds are
drawn per stage and appended to `seeds.txt` in each replicate directory, so the
trajectories are reproducible only in distribution, which is the intent.

**Inputs not generated by these scripts**

| file | origin |
|---|---|
| `3K3K`, `3QN1` mmCIF | RCSB; per-chain facts asserted from `data/structure_provenance.json` (§28h) |
| `data/md/A8S.mol2`, `.frcmod` | `antechamber -i data/aba_deprot.sdf -fi sdf -c bcc -nc -1 -at gaff2 -rn A8S -s 2 -dr no`, then `parmchk2`; AM1 with `qmcharge=-1`. **Carries the 3QN1 crystal pose** — the frame matters (§30b) |
| `data/structures_191/A8S_open_frame.mol2` | written by script 69: same atoms, types and AM1-BCC charges, coordinates only rotated |
| `data/aba_params/A8S.params` | Rosetta ligand params, used solely for the 69 repack |
| `pyr1_open_A.pdb`, `pyr1_closed_A.pdb` | script 30's references, 178 residues — still the loop-dynamics references, and NOT the same numbering as `md191` |

**⚠ GPU selection.** Script 70 requests `--gres=gpu:1` — *any* card — and excludes
by node, because neither "newest" nor a single named type is right.

| nodes | cards | feature | status for amber/22 |
|---|---|---|---|
| gpu01–03 | k80 | `gpu_legacy` | too OLD — below the build's compute-capability floor |
| gpu05 | p100 | `gpu_prev` | runs, but a 2016 Pascal card — excluded for **throughput**, not compatibility |
| gpu06–08 | a100 | `gpu_latest,gpu_highmem` | ✅ verified |
| gpu09/10/12 | ada6000 | `gpu_latest` | ✅ verified |
| gpu11 | h100 | `gpu_latest,gpu_highmem` | untested; left in (2 cards) |
| gpu13–14 | blackwell6000 | `gpu_latest,gpu_highmem` | **FAILS in seconds** |

`pmemd.cuda` from `amber/22_mpi_cuda` has **no kernel for the Blackwell cards**:
tasks 2 and 4 of job 27547457 landed on `gpu14` and died in 10 and 3 seconds with
`cudaMemcpyToSymbol: SetSim copy to cSim failed invalid device symbol`. They were
resubmitted as job 27547617.

⚠ **`--constraint=gpu_latest` is the obvious fix and it is wrong.** The cluster
tags gpu13–14 as `gpu_latest`, so asking for the newest hardware selects exactly
the cards that fail. The feature names describe the silicon, not what this AMBER
build was compiled for. Exclude by node:
`--exclude=gpu01,gpu02,gpu03,gpu05,gpu13,gpu14`.

**What a reproducer should check rather than trust:** `69b_verify_md191.sh` and the
assertions inside 58 are the reproduction test. If the four cells verify and the
15 unit maps pass their `RESIDUE_LABEL` check, the pipeline is wired correctly;
if they do not, nothing downstream is worth reading.

---

## 31. Headroom check: steric admissibility does not generate menus (2026-08-18)

### 31a. The idea and why it was worth testing

Stage 1 failed because we asked a score function to **rank** residues, and
ref2015's ranking is wrong in this pocket for a quantified reason (§23j: +10.1 REU
to bury K59's ammonium). Geometry does not rank — it **excludes**. So: at each
pocket position, which residues can be placed at all with the ligand present? A
hard yes/no, ligand-conditional by construction, needing no null arm.

Measured **before** building anything on it, because §25 established that pocket
positions are near ligand-independent and the stage-1 clash baseline already gave
3 of 4 ground-truth positions free. `scripts/71_admissibility.py`, results in
`data/admissibility/admissibility.json`. Rotamers from Rosetta's
backbone-independent library used as a catalogue of observed geometry; accept/
reject by hard-sphere overlap on explicit Bondi radii. **ref2015 is never called**;
the only Rosetta energy used anywhere is `fa_rep` alone, to settle chi.

### 31b. Verdict: it does not clear

26 positions within 5 Å of either ligand, WT PYR1 backbone, ABA (`A8S`, 19 heavy)
vs mandipropamid (`3UZ`, 29 heavy).

| | ABA | mandipropamid |
|---|---|---|
| mean admissible, of 20 | **14.7** | **14.1** |
| mean per-position Jaccard(ABA, mandi) | colspan | **0.84** |

**Too permissive.** Three quarters of all residues are placeable at a typical
position. On the six positions that actually carry cannabinoid sensors
(59, 81, 108, 120, 159, 160) the menus give 7.1 × 10⁶ combinations against 6.4 ×
10⁷ unrestricted — a **9× reduction**, where Tian's coumarin focused library is
1.4 × 10⁵. Fifty times too big.

**And ligand-blind where it matters.** All four mandipropamid ground-truth
substitutions — K59R, V81I, F108A, F159L — are admissible, but **all four are
admissible for ABA too**. They are recalled without being ligand-conditional,
which is exactly how the stage-1 clash baseline behaved. Recall against the 45
Beltrán cannabinoid sensors is **93%** of in-pocket substitution occurrences — high
recall with no selectivity is not a filter.

### 31c. Why it fails, which is the useful part

The misses are all one kind: **V83W, V83F, V164W, V164H, A89V, E141F** — every one
a *grow* mutation. Admissibility on a fixed WT backbone rejects them because in the
unmutated pocket there is no room. Real sensors use them anyway.

Look at what the sensors actually are: **F159A + A160I** (WIN), **Y120G + A160G**
(the JWH series, 27 and 23 occurrences). These are **compensating shrink/grow
pairs** — the shrink is what creates room for the grow. A per-position filter
evaluates each mutation in a context where its partner has not happened yet, so it
is structurally incapable of seeing the dominant empirical pattern.

⚠ **This kills the per-position menu step, not the combinatorial idea.** The
prerequisite that looked obvious — build menus, then combine — is the part that
breaks. Pairs have to be enumerated directly against a re-evaluated pocket.

### 31d. Two other things the data says

**Position frequency across the 45 cannabinoid sensors** is even more saturated
than coumarin: **A160 in 38, Y120 in 33, F159 in 24** of 45. Substitutions are
reused too — A160G 27×, Y120G 23×. Within a chemical family the menu is largely
conserved; the ligand-specific signal sits at **K59** (Q/N/S/T/A/R), **H115Q**,
**V83** (W/L/F), **E141** (F/M), **V164** (H/W), **V81** (I/M).

**Nine sensor substitutions are outside the 5 Å pocket entirely** — E4G (5
sensors), Y23H, D26G, Y58F, E102K, M158I, Q169R, D184G, V190A. `PYR1^WIN` gets its
last order of magnitude (50 → 10 nM) from adding **E4G**. No pocket-based method
of any kind can reach these.

### 31e. Confidence, and the one residual defect

Controls: WT admissible at its own position **23/26**; GLY admissible everywhere;
TRP excluded in 33 of 52 position-ligand pairs. Two control failures were diagnosed
and fixed rather than absorbed into the tolerance — a peptide **bond** scored as a
2.05 Å clash (`PRO88:NV — LEU87:C`, d = 1.35 Å), and single-start chi minimisation
landing in the wrong basin (now 5 restarts).

**R79, V83 and H115 still fail.** They make no crystal contact worse than 0.4 Å, so
this is placement, not the pocket. It does not threaten the conclusion: a failure of
this kind wrongly **excludes** a residue, so it makes admissibility look *more*
selective than it is. The saturation cannot be an artefact of it — only understated.
⚠ V83 is the exception worth remembering: it fails its own control **and** its
grow-mutations V83W/V83F are among the misses, so V83's numbers specifically are
not trustworthy.

The verdict is stable across the whole tolerance sweep: Jaccard 0.83–0.95 from 0.0
to 1.0 Å, and mean |admissible| never drops below 8.5 at a tolerance where the WT
control already fails 14 of 26.

---

## 32. Pairwise enumeration against a jointly re-evaluated pocket (2026-08-18)

§31 killed the per-position menu step and gave the reason: real sensors are
compensating shrink/grow **pairs**, and a per-position filter judges each mutation
in a context where its partner has not happened. So evaluate both at once.
`scripts/72_pairwise.py`, results in `data/pairwise/pairwise.json`.

### 32a. The signal, measured before building anything on it

| ligand in WT PYR1 | max lig-protein overlap | clashing atoms | summed overlap |
|---|---|---|---|
| ABA (cognate) | **0.22 Å** | 0 of 19 | 0.29 Å |
| mandipropamid | **2.78 Å** | 7 of 29 | 31.07 Å |

The cognate ligand fits its own pocket and the non-cognate one does not, which
makes **ABA a built-in negative control**: no mutation should improve it. The clash
already names **PHE108 (2.78 Å)** and **PHE159 (1.42 Å)** unprompted — 2 of the 4
mandipropamid ground-truth mutations. K59R and V81I do **not** clash (2.86 and
3.28 Å away), so they cannot be reached by clash relief at all.

### 32b. Two objectives, and two mistakes made getting there

`lig_ov` = **summed** ligand-protein van der Waals overlap. Summed rather than
maxed so it is **additive over protein atoms**, which splits a double mutant into
fixed environment + rotamer i + rotamer j exactly, turning 25×20 × 25×20 into a
precompute plus a tiny per-pair search.

⚠ **Mistake 1 — an additive objective cannot show synergy.** Ranking on `lig_ov`
alone, every top pair had relief exactly equal to the sum of its singles, and the
best residue at position i changed with its partner in **15 of 4187 contexts
(0.4%)**, all trivial A→G. The decomposition that made it fast removed the very
thing being looked for. The grow half of a shrink/grow pair does not reduce clash —
it restores **packing** — so it is invisible to a clash-relief objective by
construction.

⚠ **Mistake 2 — counting contacts within a flat 4.5 Å rewards clashes.** The first
packing term put `83W+163Y` at the top of the Pareto front with 335 "contacts" and
a relief of **−71.7 Å**: it made the overlap catastrophically worse and scored as
the best-packed pair on the board. A contact is favourable only if the atoms touch
**without interpenetrating**, so `contacts` now counts pairs with overlap in
[−1.0, +0.25] Å.

### 32c. What it delivers once both objectives are right

Ranking by relief subject to *not losing packing* (contacts ≥ WT):

| ground truth | best packing-preserving pair | relief | rank of 13,457 |
|---|---|---|---|
| **F108A** | 81Q+108A | 24.16 | **5** (top 0.04 %) |
| **F159L** | 115H+159L | 4.39 | **75** (top 0.6 %) |
| **V81I** | 81I+83T | 0.46 | **150** (top 1.1 %) |
| K59R | 59R+81Q | −0.04 | 3907 (top 29 %) |

**3 of 4 in the top 150 of 13,457 — a 90× narrowing at 3/4 recall.** Full 4/4
recall needs 3907 pairs, i.e. K59R alone costs 26× the library.

**Shrink/grow enrichment is real and ligand-specific**: 123 of the top 200 mandi
pairs change volume in opposite directions (62 %) against **85 of 200 (42 %) for
ABA**, where chance is ~50 %. ⚠ No p-value is quoted — the top-200 pairs share
positions and are not independent, the same trap as §14a.

And the top mandi pairs are **81+108 shrink/grow combinations** — `81K+108V`,
`81Q+108A`, `81Q+108S`. The true answer is V81I + F108A. **The method finds the
right pair of positions and the wrong identity at 81**, which is §25's
positions-are-easy/substitutions-are-hard result reappearing at the pair level.

### 32d. Where the boundary actually is

K59R is not a failure of implementation, it is a failure of *category*. K59 does
not clash; its contribution is the salt bridge to the ligand, which is
electrostatic. A steric method cannot see it **by construction** — and this is the
same residue, for the same underlying reason, that defeated both LigandMPNN and
Rosetta FastDesign in §23j, where ref2015 saw the salt bridge and then overcharged
its desolvation by +10.1 REU.

So the honest scope: **this recovers sterically-driven substitutions and is blind
to electrostatically-driven ones.** That is a principled boundary, and it says what
a second stage must supply — electrostatic complementarity at non-clashing
positions — rather than leaving "it did not work" undiagnosed.

The ABA negative control behaves correctly throughout: best relief 0.27 Å against a
WT overlap of 0.29 Å, i.e. nothing to fix and nothing invented.

---

## 33. Electrostatic complementarity at non-clashing positions (2026-08-18)

§32d specified this stage: steric enumeration recovers sterically-driven
substitutions and is blind to electrostatically-driven ones, so K59R needs a
complementarity screen. `scripts/73_electrostatic.py`.

### 33a. The decisive test, and its result

The two ligands differ exactly where it matters — **ABA is an anion (−0.97)**
presenting a carboxylate; **mandipropamid is neutral (+0.10)** presenting an amide.
ABA's answer at 59 is WT Lys; mandipropamid's is Arg. One method, two ligands, one
flip.

| | rank 1 | LYS rank | ARG rank |
|---|---|---|---|
| ABA | **K** ✅ | 1 | 2 |
| mandipropamid | K ❌ | 1 | 3 |

**Half a pass.** With ABA it puts Lys first *for the right reason* — a bidentate
salt bridge, NZ–O4 **2.83 Å** and NZ–O3 **3.31 Å**. That is worth stating plainly:
stage 1's null arm retained K59 in **0 %** of trajectories even with the cognate
ligand (§23j), and this retains it at rank 1. The pathology that defeated ref2015
is genuinely absent here.

With mandipropamid it does not flip. **K59R was missed.**

### 33b. Why — sampling, not scoring, and that is testable

Transplanting the crystal Arg59 from 4WVO onto our model (core superposition, 170
Cα, 0.47 Å) and scoring it under *this script's own criteria*:

- **2 hydrogen bonds**: NE–O2 **2.64 Å**, NH1–O2 **3.26 Å**
- worst genuine steric overlap **0.56 Å** (against Tyr58) — **below** the 0.75 Å
  feasibility threshold
- its χ angles are −78 / 176 / **100** / 69; **χ3 ≈ 100° is non-rotameric**

So the conformation is acceptable to the scoring and would have been ranked highly.
It was never proposed. The backbone-independent rotamer library does not contain
it, and §31's chi-minimisation fix cannot help: minimising against `fa_rep`
relieves strain, it does not *create* a hydrogen bond.

**This is the first time the K59 miss has been localised to sampling.** In §23j it
was ref2015's desolvation; in §32 it was category (no clash). Here the physics is
right and the search is the limit — a different and much more tractable problem.

### 33c. Two bugs found, one of which reaches backwards

1. **The van der Waals table has no hydrogen radius**, so `VDW.get("H", 1.7)`
   silently gave hydrogens a carbon-sized 1.70 Å instead of 1.20 Å, inflating every
   H contact by 0.5 Å. Confined to diagnostics here — the script's steric test is
   heavy-atom only — but it is in `lib`-shaped code and will bite elsewhere.
2. **Hydrogen bonds are counted as steric clashes.** A donor–acceptor pair at
   2.64 Å has a hard-sphere overlap of 1.55 + 1.52 − 2.64 = **0.43 Å**. Polar–polar
   contacts are *supposed* to interpenetrate. ⚠ This reaches back into **§31 and
   §32**, whose feasibility tests use the same hard-sphere rule, and is the likely
   explanation for §31e's three unresolved control failures — **R79** (salt bridge
   to E94), **V83** and **H115** are exactly the polar/charged positions such a rule
   would wrongly exclude. Not yet confirmed; recorded so it is not rediscovered.

### 33d. The screen fails its own controls, in the predicted direction

Ranking every non-clashing position:

| control | ABA | mandipropamid |
|---|---|---|
| WT retained as top choice | **1/27 (4 %)** | **4/30 (13 %)** |
| ARG or LYS picked | **16/27 (59 %)** | 8/30 (27 %) |
| top choice differs between ligands | colspan | 18/27 (67 %) |

With no desolvation term, everything wants to be charged against an anionic ligand
— 59 % Lys/Arg for ABA against 27 % for the neutral one. This is the **predicted
mirror image** of ref2015's failure, written into the script's docstring before the
run, and the Arg-everywhere control exists precisely to catch it. It did.

So: **not usable as a ranker.** 4 % WT retention is disqualifying, and the 67 %
ligand discrimination cannot be credited while most of it is monopole attraction.
What the screen *does* do correctly is evaluate a supplied conformation — which is
the role §33b shows it should have: a **filter over candidates**, with sampling
supplied by something else.

### 33e. Three silent input failures, all caught by assertions

None of these raised an error on its own, and all three would have produced
confident numbers:

1. **The params were never used.** Both stage-1 params declare `NAME LIG` while
   `data/stage1/wt_*.pdb` contain `A8S`/`3UZ`. Rosetta silently built the ligand
   from its PDB-components dictionary — `pdb_A8S`, **neutral, charge sum 0.000**
   where the params sum to −0.970.
2. **Renaming did not fix it**: those PDBs carry RCSB atom names (`CAC CAB OAI…`)
   and the params carry generation order (`C12 N1 C11…`). They were never a matched
   pair. The matched inputs are stage 1's own Rosetta-written
   `results/stage1_rosetta/_input_wt_*.pdb`, ligand `LIG`, 38 and 51 atoms.
3. **Two ligands cannot share a name in one process.** Both params are `LIG`, and
   PyRosetta's residue type set is global, so a second `init` does not replace the
   first type — mandipropamid loaded with **ABA's charges** (−0.970 instead of
   +0.100). Fixed by one subprocess per ligand.

Caught only because the ligand net charge is asserted against the params file.
⚠ **§32 is unaffected**: it used coordinates and element-derived radii only, never
a charge, and the ligand coordinates are identical either way — confirmed, the
mandipropamid pose in our WT model is the 4WVO pose to **0.01 Å**.

---

## 34. What the field actually does — six papers, read from Methods (2026-08-19)

PDFs and extracted text in `data/papers/`. Read from the Methods sections, not
abstracts. Figures were not read.

### 34a. Leonard et al. 2026 — the only comparable computational protocol

| stage | what they did |
|---|---|
| conformers | **TREMD**: 27 replicas, **300–450 K**, **100 ns**, exchange every 2 ps. Boltzmann occupancy; keep "high occupancy" **>4 %**, discard **<1 %**. Code: `ajfriedman22/SM_ConfGen` |
| **pose** | **"dock to sequence"** — the WT sequence is **mutated to match a known low-affinity binding sequence found by Y2H library screening**, then conformers are SVD-superposed onto the ABA ketone in **3QN1**, treating the nitazene nitro as the required H-bond acceptor |
| perturb | PyRosetta Rigid Body Perturbation |
| filter 1 | no backbone clash, fast-grid search against a **poly-glycine shave** |
| filter 2 | an N or O close enough to the **latch water** to hydrogen bond |
| refine | PyRosetta Packer relax/repack, then full-protein clash check |
| score | standard PyRosetta scorefunction, **< −300 REU** "generally considered realistic" |
| design | Rosetta **FastDesign**. LigandMPNN comparison at **temperature 0.2**, **0.1 Å Gaussian noise**, `ligandmpnn_v_32_010_25.pt` |
| build | oligo pools in 3 cassettes, Golden Gate, Y2H |

**No ΔΔG, no shape complementarity, no SASA filter.** Clash + water H-bond + one
total-score cutoff. Simpler than assumed.

⚠ **Two things this confirms.** First, the user's read is correct *verbatim*: the
protocol is **seeded by a weak library hit** — it needs a known binding sequence
before it can place a pose. It improves a known weak binder; it does not find
specificity from nothing. Second, the **latch-water H-bond is a hard filter** in
their pose generation, and §23a measured that mandipropamid sits **5.07–5.19 Å**
from that water in 4WVO/8EY0. Their protocol, applied to mandipropamid, would
reject the correct pose.

### 34b. Rationalizing Diverse Binding Mechanisms (ACS Chem Biol 2024)

Same collaboration (Whitehead / Shirts / Cutler) and the mechanistic foundation
under §34a. Three things land directly on our work:

1. **"It is the H-bonding capability at position R59 that is essential for ligand
   recognition."** Independent support for treating 59 as the electrostatic
   position (§33).
2. **R79 hydrogen bonds the main-chain carbonyl of F52 at >90 % occupancy in both
   apo and holo, and does not contact the ligand.** ⚠ This independently explains
   §31e's unresolved **R79** control failure — its side chain is locked in a tight
   intra-protein H-bond, and §33c showed our hard-sphere test scores exactly that
   as a clash. Two routes, same conclusion.
3. **3 × 300 ns MD of PYR1^mandi and PYR1^WIN with ligand + HAB1 gave only ONE
   ligand conformer**, matching crystals 4WVO and 7MWN; TREMD puts that bound
   conformer **< 1 k_BT** above the solution minimum.

### 34c. An et al. 2024 — and the honest part

Rotamers from RDKit, Rosetta-relaxed with constraint, ddG-scored, clustered at
**1.5 Å RMSD**, lowest-energy per cluster. RIFgen/RIFdock onto **9,703
pseudocycles**; pockets pre-packed with large hydrophobics (V/L/I/F/Y). Metrics:
`CMS` (contact molecular surface), `dsasa`, `ddg_norepack`, `atomic_depth`
(**6–7 Å**), `hole_around_lig`, `total_hb_to_lig`. AF2 gate: **pLDDT > 90**,
**Cα-RMSD < 3 Å**, **PAE < 5**.

⚠ **Their cutoffs were "selected by manual inspection of docks at different
ranks."** Stated plainly in the supplement. The field does not have principled
thresholds either — which is worth remembering before treating any published
number as a standard.

### 34d. Tian 2025 and Park 2023 — the library baseline

Tian's targeted libraries are built from **sequence profiles of first-round hits**:
coumarin **77,327** members, TNT **506,229**. That biasing *works experimentally* —
it isolated sensors for 4-methylumbelliferone, 7-methoxycoumarin, TNT, DNT and
2ADNT, all missed in round 1. ⚠ Note the asymmetry: the same information has not
yet helped us **computationally** on the PFAS / TNT / coumarin controls.
Park 2023 is library-driven throughout (~400,000 yeast colonies; a ~12,000-clone
library), confirming that outside Leonard this field is screening, not design.

### 34e. What we should change

1. **Adopt a solution-conformer filter for prospective ligands** — TREMD, keep
   conformers within ~1 k_BT of the minimum (§34b validates the criterion). Our
   retrospective work uses crystal poses, which is correct, but a new ligand has no
   crystal pose.
2. **Do not adopt the latch-water H-bond as a hard filter** (§23a, now reinforced).
3. **Fix the H-bond-as-clash bug** (§33c) — R79 is now confirmed as a genuine
   intra-protein H-bond, not a modelling artefact.
4. **Our benchmark is harder than Leonard's by construction.** They start from a
   weak hit; we start from the ligand. Any comparison must say so.

---

## 35. The hydrogen-bond-as-clash fix (2026-08-19)

`scripts/lib_sterics.py` — one clash rule, now shared by 71, 72 and 73 so they
cannot drift apart.

**The bug.** A hard-sphere test treats every close approach as a clash, and a
hydrogen bond *is* a close approach: donor–acceptor at 2.6–3.2 Å against a Bondi
radius sum of 3.07 Å registers as ~0.4 Å of "overlap". The rule now exempts a pair
when one atom can donate and the other accept, treating it as a clash only below
**2.5 Å** — shorter than any real hydrogen bond, and below the 2.64 Å of the 4WVO
Arg59 contact. Acceptor–acceptor and donor–donor pairs are **not** exempt, which is
why donor and acceptor are tracked separately rather than collapsed to "polar".

Unit-tested: donor/acceptor at 2.8 Å → exempt; acceptor/acceptor at 2.8 Å → 0.27 Å
clash; donor/acceptor at 2.2 Å → 0.87 Å clash. On the real structure, R79's crystal
side chain goes from **0.306 Å overlap to 0.000**.

**What it changed — and what it did not.** It matters exactly where it should, at
hard-sphere tolerance, and is invisible at the tolerances §31 actually used, because
a 0.5 Å tolerance already absorbs a typical H-bond's 0.43 Å:

| tol (Å) | WT-ok before | WT-ok after |
|---|---|---|
| 0.00 | 12/26 | **16/26** |
| 0.25 | 21/26 | 21/26 |
| 0.50 | 23/26 | 23/26 |

§31's verdict is unchanged: mean |admissible| 14.7 / 14.1, Jaccard 0.84, still
saturated and largely ligand-blind. §32's ranking barely moves — F108A rank 5
(unchanged), V81I 150 → 147, F159L 75 → 81 — which is a useful stability check on
that result rather than a revision of it.

⚠ **CORRECTION to §33c.** I wrote that this bug was the "likely explanation" for
§31e's R79/V83/H115 control failures. **It is not.** With the exemption applied,
those three still fail at the same overlaps (0.65 / 0.69 / 1.12). The crystal R79 is
now clash-free, but no *library rotamer* reproduces it — so §31e's original
diagnosis, wrong chi basin, was right and the later speculation was wrong. The
literature agreement in §34b (R79 H-bonds F52's carbonyl at >90 %) confirmed the
chemistry, not my explanation of the control failure.

## 36. Coupled moves: sampling is not what defeats K59

`scripts/74_coupled_moves.py` + `74b`. Stage 1 re-run with the sampler the
literature prescribes (§34): Ollikainen & Kortemme's **coupled moves**, which
samples sequence, side chains, backbone **and the ligand's pose** together and gave
a 5.75× improvement on exactly this task. Same designable positions, same two arms,
same ligand-swap scoring rule as §23j — **only the sampler changed**.

Job 27560390, 100 tasks on `cutlerlab`, 50 trajectories per arm × 1000 trials,
`ligand_mode` on, backrub backbone mover, ref2015, all rc=0.

| ground truth | f(mandi) | f(ABA) | delta | verdict |
|---|---|---|---|---|
| **F108A** | **0.74** | 0.00 | **+0.74** | recovered, cleanly ligand-conditional |
| V81I | 0.08 | 0.28 | −0.20 | **inverted**, as in both stage-1 methods |
| K59R | 0.00 | 0.00 | 0.00 | missed |
| F159L | 0.00 | 0.00 | 0.00 | missed |

### 36a. The pre-registered reading, applied

§36 was registered before the run with two competing diagnoses and a rule for
telling them apart: *K59R recovered ⇒ sampling was the limit; K59R still missed ⇒
the scoring diagnosis stands and no sampler will fix it.*

**K59R was missed, and the decisive number is the control.** In the ABA arm — with
the **cognate** ligand, where K59 is unambiguously correct — coupled moves retains
K59 in **0 % of 50 trajectories**. Identical to fixed-backbone FastDesign, and to
every favour-native weight up to 1.5.

That cannot be a sampling failure. Lysine is in every rotamer library, its native
conformation is in the input pose, and `IncludeCurrent` was on. Backbone and ligand
flexibility were added and it changed nothing. **The +10.1 REU Lazaridis–Karplus
desolvation penalty for burying the ammonium (§23j) is confirmed as the cause, now
by a second and methodologically independent sampler.**

Coupled moves *is* better overall — WT retention rises from **29 % to 42 %** — so
the flexibility helps everywhere except the one position whose problem is energetic.

### 36b. Where each K59R miss now sits

Four methods, four localisations, and they are no longer contradictory:

| method | why it missed K59R |
|---|---|
| LigandMPNN (§23h) | learned prior; no mass on a non-clashing substitution |
| Rosetta FastDesign (§23j) | **scoring** — buried-charge desolvation |
| steric pairs (§32) | **category** — K59 does not clash at all |
| electrostatic screen (§33b) | **sampling** — crystal Arg is non-rotameric, χ3 ≈ 100° |
| **coupled moves (§36)** | **scoring, confirmed** — better sampling, same 0 % retention |

§33b's sampling diagnosis was about *my* hand-rolled rotamer enumeration, and stands
for that script. It was never a claim about the design stage, and §36 now closes the
design stage: the limit there is the energy function, not the search.

⚠ **What this rules out.** Adding flexibility — more rotamers, backbone moves,
ligand moves — will not recover K59R. Any method that scores buried charge with
ref2015's desolvation term inherits the same failure. That is a constraint on every
future stage of this pipeline, and it is now measured twice rather than argued once.

---

## 37. MM-GBSA rescoring: the K59 flip, at last (2026-08-19) — ⚠ RETRACTED BY §40

> **⚠ The headline result in this section did not replicate.** The 50 ps
> ensembles used here scored the WT–ABA reference ~7.75 kcal/mol too favourably;
> at 250 ps the pre-registered test gives the opposite answer. See **§40**.
> Kept unedited for the record — the reasoning was sound, the ensemble was not.

§36 established with two independent samplers that K59R fails for a reason inside
ref2015 — the +10.1 REU desolvation penalty for burying the ammonium. So the
solvation model was replaced: generalised Born instead of the pairwise
Lazaridis–Karplus approximation. `scripts/75_mmgbsa_build.py`, `75b`, `75c`, `75e`.

Rosetta builds the structures and **never scores them**; every energy below is Amber
(ff14SB + GAFF2/AM1-BCC, igb=8, 0.15 M salt).

### 37a. The pre-registered criterion is met

ddG relative to WT within each arm, so receptor, ligand and force field all cancel:

| variant | ddG mandi | ddG ABA | difference | reading |
|---|---|---|---|---|
| **K59R** * | **−1.53** | **+7.06** | **−8.58** | **prefers mandipropamid** |
| K59Q | +0.04 | +13.73 | −13.69 | prefers mandipropamid |
| K59N | +2.53 | +8.93 | −6.40 | prefers mandipropamid |
| V81I * | −0.19 | +5.84 | −6.03 | prefers mandipropamid |
| F108A * | −1.69 | +2.23 | −3.92 | NULL, within the spread |
| F159L * | −1.19 | +6.96 | −8.15 | prefers mandipropamid |

\* a real 4WVO mutation. n = 100 frames per arm; per-frame sd 2.2–3.4 kcal/mol,
propagated in quadrature and **not** divided by √n, because 100 frames from one
50 ps run are strongly correlated.

**Within position 59 the ordering is right in both directions**, which is the whole
test and the thing no previous method managed:

* with **mandipropamid**: R (−1.53) is better than Q (+0.04) and N (+2.53) — 4WVO is
  K59R.
* with **ABA**: wild-type Lys beats every substitution, R least-badly (+7.06) and Q
  worst (+13.73) — the cognate receptor should prefer its own residue.

ref2015 could not produce that flip at any favour-native weight, and coupled moves
retained K59 in 0 % of trajectories. Changing the solvation model did it.

### 37b. Three caveats that must travel with the number

1. **The ensembles are 50 ps.** GBn2 runs at ~380 steps/min on one core for this
   2,900-atom system — measured, and pmemd is no faster — so a converged run was not
   affordable. 50 ps is enough to give frames that genuinely differ; it is **not**
   enough to call the averages converged.
2. **The `difference` column is dominated by the ABA arm.** Mutations are near-neutral
   for mandipropamid (−1.69 to +2.53) and strongly unfavourable for ABA (+2.23 to
   +13.73). The ligand-swap rule is therefore mostly detecting *damage to the cognate
   complex*, which is a real and expected signal but not the same as demonstrated
   complementarity for the new ligand.
3. **F108A comes out NULL**, though it relieves the largest clash (2.78 Å, §32a). The
   explanation is mechanical: minimisation and MD relieve that clash in the WT complex
   too, so by the time an energy is computed the strain the mutation exists to remove
   has already been absorbed. Structure-relaxed rescoring is blind to it — the same
   shape of error as §32b's additive objective hiding compensation.

### 37c. ⚠ The first attempt's error bars were fabricated by my own design

The first run built 8 structures per variant by repacking with 8 different Rosetta
seeds and reported the spread as the uncertainty. **All 8 mandipropamid repacks came
out byte-identical**, and the ABA arm collapsed to 2 distinct structures after
minimisation. A 7-residue packing problem has one optimum and annealing finds it
every time, so the seeds sampled nothing and every `sd = 0.00` was an artefact of
zero diversity rather than a measure of precision. The verdicts read off that table
were meaningless and are discarded; 75e replaces them with real MD ensembles, and
75c now refuses to report a row whose sd is exactly zero with n > 1.

Two smaller traps, both fixed in place: `igb=8` requires **mbondi3** radii or sander
will not start at all, and `ante-MMPBSA.py` writes no complex topology when its input
is already unsolvated — the complex prmtop is the input itself.

---

## 38. PRE-REGISTRATION: WIN 55,212-2 as a held-out test (2026-08-19)

> **Status 2026-08-21: still sealed, NOT run.** Its gating condition — MM-GBSA
> recovering the known ABA/mandipropamid answers — was not met (§40d). Protocol
> below is unchanged and has received no WIN-specific tuning.

Written **before** the extended MM-GBSA run, and before anything about WIN has been
computed, so that WIN remains a genuine held-out test rather than a third training
set. Every choice below is fixed here; if any of them changes later, the change and
its reason get recorded and the test is reported as no longer blind.

### 38a. Why WIN, and why now

Everything in §31–§37 was developed against **two** ligands: ABA (cognate) and
mandipropamid (non-cognate). Each round of method-fixing used those two to decide
what to change. That is exactly the setup in which a method quietly overfits, and
the user raised it before the extension was run rather than after.

WIN 55,212-2 is the right third case because it is **the same positions with
different answers**:

| position | mandipropamid (4WVO) | WIN (Beltrán PYR1^WIN) |
|---|---|---|
| 59 | **R** | **Q** |
| 159 | **L** | **A** |
| 160 | — | **I** |
| 81 | I | — |
| 108 | A | — |

So a method that has learned "position 59 wants Arg" fails; one that has learned to
read the ligand should say Arg for mandipropamid and Gln for WIN. K59Q is also
neutral, so the buried-charge pathology that dominated §23j is absent — a different
part of the scoring function is being tested.

All four high-sensitivity WIN sensors agree on the pocket set (**K59Q, F159A,
A160I**; the best adds only surface mutations E4G/Y23H/D26G), so the ground truth is
unambiguous.

### 38b. What will be run, fixed now

1. **Protocol frozen at whatever §39 settles on for ABA/mandipropamid.** No
   WIN-specific tuning of ensemble length, force field, solvation model, candidate
   set or scoring rule. If WIN needs a different setting to work, that is a finding
   about the method, not a setting to adopt.
2. **The null is a ligand swap, as always**: WIN vs ABA, scored as
   `ddG(WIN) − ddG(ABA)`.
3. **Receptor**: PYR1. If a WIN-bound pose is needed it comes from **7MWN**, which is
   **PYL2**, not PYR1 — 88 % identity over the 25 ABA-proximal positions, globally
   alignable to 0.55 Å. ⚠ Beltrán transposed mutations PYR1→PYL2 to get that crystal,
   so using it means transposing **back**, and the numbering must be re-derived by
   residue identity, never assumed.
4. **Success**: at position 59, Q ranks above R **and** above wild-type Lys for WIN,
   while for mandipropamid R still ranks above Q. Both directions, or it is not a
   demonstration of ligand-reading.
5. **Partial**: the right positions (59, 159, 160) surface without the right
   identities. That is §25's positions-are-easy result again and is reported as such.
6. **Failure**: the method gives WIN the same answer it gives mandipropamid. That is
   overfitting to the two-ligand training pair, and it retires the approach as a
   ligand-discriminating method regardless of how well it does on ABA/mandipropamid.

### 38c. One thing that would invalidate the test

If the WIN pose has to be built by a procedure that used mandipropamid or ABA
information to place it, the test is contaminated. The pose must come from 7MWN
directly, or from a poly-glycine dock that never sees the other two ligands.

---

## 39. Current status (2026-08-20) — ⚠ superseded by §40 for the MM-GBSA arm

§21 was written on 2026-08-10 and still listed the first MD campaign as running.
This replaces it. **Written before compaction so the state is recoverable without
the conversation.**

### 39a. Running right now

| job | what | state at 2026-08-20 | ETA |
|---|---|---|---|
| **27676964** | MM-GBSA extended — 23 variants × 2 arms, 250 ps ensembles, 46 tasks, `cutlerlab` | 37,500 / 125,000 steps; 0/46 scored | overnight |
| **27684153** | cavity measurement, all 266 helix-grip hits, `cutlerlab` | 63 / 266 streamed to disk | ~1 h |
| **27547457** | conformation × occupancy factorial MD, GPU | **5 / 12** replicates complete | ~1.5 days |

Every one is a SLURM **batch** job, deliberately — an interactive allocation kills
its processes on exit, `nohup` included.

### 39b. What each is for, and what to do when it lands

**MM-GBSA extended** is the test that separates a method from a signal. §37 showed
MM-GBSA wins a 3-way *within-position* comparison at residue 59 — which we half
expected. The open question is whether the four known mutations rank highly inside a
**23-variant pool**. Candidates are the **stratified** top of §32's pairwise ranking:
the best pair per distinct *position pair*, because the raw top-24 all contain
position 108 and would have made the test "which partner goes with F108X". The
ground-truth pairs sit at ranks 18 / 25 / 49 of the full list.
→ aggregate with `scripts/75c_mmgbsa_aggregate.py`.

**Cavity measurement** adds the **21 solved domains that were never measured** — Arm 3
used a stricter alnTM cut and stopped at 147. Results stream to
`homolog_cavities_partial.csv` and land in `homolog_cavities_full.csv`, deliberately
*not* the original 147-row file.
→ then re-run `scripts/79_graft_geometry.py` over the full set; **3H3Q** (CERT START
domain + ceramide) and **4QDC** (a steroid) finally get volumes.

**Factorial MD** fills the 2×2 that §19b set up and only ever ran on the diagonal.
→ analyse with scripts 58–60 pointed at `data/md191`.
⚠ `preempt_gpu` caps this account at **1 concurrent GPU**; `gpu` allows 4. Moving
pending tasks across cut the estimate from 6.4 days to ~1.5. Do it again if they
pile up: `scontrol update jobid=<id> partition=gpu account=cutlerlab`.

### 39c. The next scientific step is already fixed

**§38: WIN 55,212-2 as a held-out test.** Written *before* the extension so it stays
blind. Same positions, different answers — 59 wants **Q** not R, 159 wants **A** not
L, and 160I is a position never scored. No WIN-specific tuning; success requires
**both** directions.

### 39d. Parked, with reasons

- **Graft direction** (`METHODS_REVIEW.md` appendix) — lower priority by decision,
  because the literature success rate is poor. Best candidate **2BK0**, the celery
  allergen Api g 1, but core RMSD is 3.4–4.0 Å at 10–17 % identity: a rebuild with a
  template, not a transplant.
- **Pose-sensitivity sweep** — perturb the mandipropamid pose, re-run the pairwise
  screen, plot top-150 recall against pose RMSD. Decides whether a pose-search
  programme is worth building at all.
- **S6–S8 non-cognate MD** and the **ABA-enantiomer / phaseic-acid ladder** — built or
  specified, not run.

### 39e. Where the project actually stands

The K59R problem has been missed by five methods for four distinct reasons, and only
the fifth solved it. That progression is the substance of §23–§37:

| method | why it missed K59R |
|---|---|
| LigandMPNN | learned prior; no mass on a non-clashing substitution |
| Rosetta FastDesign | **scoring** — buried-charge desolvation, +10.1 REU |
| steric pairwise | **category** — K59 does not clash at all |
| electrostatic screen | **sampling** — the crystal Arg is non-rotameric, χ3 ≈ 100° |
| coupled moves | **scoring, confirmed** — better sampling, same 0 % retention |
| **MM-GBSA** | **solved it** — generalised Born instead of Lazaridis–Karplus |

What is *not* yet established is whether that generalises. Two things test it: the
extended pool now running, and the WIN held-out set.

---

## 40. The extended MM-GBSA run refutes §37 (2026-08-21)

Job 27676964 finished 46/46 tasks clean (~9.6 h each, 250 ps GB ensembles,
100 frames). It was built to answer one question — *does MM-GBSA rank the four
known mandipropamid mutations highly inside a 23-variant pool, or was §37 a
3-way within-position win on mutations we already knew the answer to?*

It answered a different and more important question first.

### 40a. The 50 ps result does not replicate

The same seven variants exist at both lengths. Six of the fourteen numbers moved
by less than 1.7 kcal/mol. One moved by **+7.75**:

| variant | arm | 50 ps | 250 ps | shift |
|---|---|---|---|---|
| **WT** | **ABA** | **−32.74** | **−25.00** | **+7.75** |
| WT | mandi | −36.33 | −38.00 | −1.67 |
| K59R | ABA | −25.69 | −26.71 | −1.02 |
| K59R | mandi | −37.85 | −38.39 | −0.54 |

The one that moved is the **reference**. Every ddG in the ABA arm subtracts
WT–ABA, so a 7.75 kcal/mol error in that single number propagates into all of
them at full strength. Recomputing the pre-registered §38 test:

| | ddG mandi | ddG ABA | selectivity | verdict |
|---|---|---|---|---|
| 50 ps | −1.53 | **+7.06** | **−8.58** | PASS — prefers mandipropamid |
| 250 ps | −0.40 | **−1.71** | **+1.31** | **FAIL — prefers ABA** |

The entire §37 "K59 flip" was the WT–ABA complex being scored before it had
relaxed. At 50 ps it still sat near its minimised starting structure, which
flattered it by ~8 kcal/mol; the mutant complexes, already perturbed, had no
such advantage. The apparent selectivity was the reference falling away from the
mutants, not K59R moving toward mandipropamid.

**This is exactly the failure the run was designed to expose**, and it is the
third time on this project that a result has turned out to be a property of the
protocol rather than the protein (cf. §31 params-never-loaded, §37's own
byte-identical "ensembles"). The pattern is consistent enough to be a standing
rule: *a number that has never been recomputed at a different setting is not yet
a measurement.*

### 40b. 250 ps is not converged either

Splitting each ensemble into halves (`scripts/75d_mmgbsa_convergence.py`):

- median |drift| **1.05 kcal/mol**; **23 of 46** runs drift more than 1.0
- WT–ABA drifts **+2.62**, in the *same direction* as its 50→250 ps move

So WT–ABA has not settled at 250 ps — it is still climbing. Extending to 1 ns
would very likely move it again.

Note what the block standard errors in §40c do *not* say. They are within-run
scatter about the current mean, so they stay near 0.3–1.3 kcal/mol while that
mean marches. **An error bar cannot see a drift it is centred on.** Only the
half-split can.

### 40c. What survives: ranks, not signs

The two conclusions are not equally damaged, and separating them is the whole
value of the run.

A uniform reference error adds the same constant to every ddG in an arm. That
**flips sign-based verdicts** — which is what killed §37 — but **cannot reorder
them**. Rank-based conclusions are therefore invariant to the specific failure
above. On selectivity rank (`results/mmgbsa/ranking.txt`):

| rank | variant | ddG mandi | ddG ABA | selectivity |
|---|---|---|---|---|
| 1 | **F159L** * | −0.69 | −0.70 | +0.01 |
| 3 | **K59R** * | −0.40 | −1.71 | +1.31 |
| 4 | **V81I** * | +1.00 | −0.36 | +1.36 |
| … | | | | |
| 18 | **F108A** * | +2.46 | −5.59 | +8.05 |

\* = a substitution in the 4WVO mandipropamid sensor

Three of four land at ranks 1, 3, 4 of 22; mean rank 6.5 against 11.5 expected
by chance, exact one-sided permutation **p = 0.049**.

That is a real signal and it is *not* nothing — but it is one point below the
conventional threshold, on n = 4, chosen post hoc after the pre-registered test
had already failed. It is also bounded by §40b: the spacing between ranks 1–5 is
0.01–1.8 kcal/mol, and the median per-variant drift is 1.05. **The method is
ranking on differences the same size as its own irreproducibility.**

F108A at rank 18 is the same NULL as §37b, now worse, and for the known reason —
relaxation absorbs the clash the mutation exists to relieve.

### 40d. WIN stays sealed

The §38 pre-registration gates the WIN 55,212-2 held-out test on MM-GBSA first
recovering the known ABA/mandipropamid answers. **It did not.** The pre-registered
test fails outright at 250 ps, and the surviving rank-based signal is marginal and
post-hoc.

Running WIN now would spend the one genuinely held-out ligand on a method whose
own positive controls are unresolved — and whichever way it came out, the result
would be uninterpretable: a hit would be unattributable, a miss unattributable.
**Held out means held out.** §38 is unchanged and no WIN-specific tuning has been
applied to anything.

### 40e. What would actually settle it

In priority order, cheapest first:

1. **Longer ensembles on the reference alone.** WT–ABA and WT–mandi at 1–2 ns,
   nothing else. If WT–ABA stops moving, the arm becomes usable; 2 tasks, not 46.
2. **Independent replicates, not longer single runs.** 3 × 250 ps from different
   velocity seeds gives a between-run error bar, which is the honest one for a
   quantity this drifty. Within-run block SE has now twice been the thing that
   hid a problem.
3. **Test the rank claim where it is cheap.** The ranking is the part that
   survived; it can be checked against Beltrán's cannabinoid sensors — real
   substitutions, already in hand (§30), never used for tuning — without touching
   WIN.

Until (1) or (2) lands, MM-GBSA is a **ranking heuristic under audit**, not the
scoring solution §37 claimed. §37's headline is retracted; its diagnosis of
ref2015 desolvation (§23j, §36) is untouched — that rested on Rosetta numbers,
not these.

---

## 41. Why the implicit reference drifted — measured, not assumed (2026-08-21)

§40 established *that* WT–ABA moved 7.75 kcal/mol. This is *why*.

**The energy decomposition localises it.** Of the +7.75, **EGB (polar solvation)
contributes +7.65**; VDWAALS contributes +1.51 and ESURF +0.15. Mandipropamid's
EGB moved only +1.60.

**But the ligand did move**, and an inference I drew from the decomposition alone
was wrong. Measured (protein-fitted, `scripts/82`):

| | 1st half | 2nd half | max | final |
|---|---|---|---|---|
| ABA, implicit | 1.09 Å | **3.18 Å** | 4.88 Å | 4.13 Å |
| mandipropamid, implicit | 1.45 Å | 1.70 Å | 2.45 Å | 1.72 Å |

**The salt bridge never breaks.** K59 NZ ↔ nearest carboxylate oxygen: 2.77 Å at
frame 1, 2.83/2.81 Å by half, **100% intact** below 4 Å throughout. The lysine
side chain *tracks* the ligand — ABA, K59 and the gate migrate together as a unit.
This is a coordinated relaxation of the whole site, not a release event.

**Mechanism.** The build contained **zero explicit waters** (verified: the input
PDB had none and `leap` did no solvation). ABA's crystallographic pose is held
partly by ordered bridging waters; GB replaces each with a featureless continuum,
and a continuum cannot donate a directional hydrogen bond. So the crystal pose is
a minimum of *reality* but not of *the Hamiltonian we integrated*, and the system
slid to the model's own minimum ~2 Å away.

That inverts the intuition: **WT–ABA had the furthest to fall precisely because it
was the only system positioned by real physics.** Every other run — 44 mutants and
WT/mandi — started from a repacked or docked pose the model already liked.

**Why EGB and not VDW.** VDW measures *contact*, which a 2 Å slide inside an
enclosing pocket preserves. EGB measures *burial depth of charge* and goes as
q²/R_eff. ABA carries a localised formal −1 (`A8S.mol2` = −1.001); mandipropamid
is neutral (`3UZ.mol2` = −0.003) with its EGB spread over many small dipoles.
**Both effects compound: ABA moves more, and its energy is more sensitive to
moving.**

---

## 42. Explicit solvent: the pose is fixed, the precision is not (2026-08-21)

Both references rebuilt with `lib_solvate.sh` (ff19SB/OPC, truncated octahedron,
12 Å buffer, 0.15 M **KCl**), 10 ns × 3 velocity seeds, scored with the *same*
igb=8 settings on water-stripped frames. Exactly one variable moved.

**No crystallographic waters, deliberately** (user, 2026-08-21). Placing ABA's
ordered waters while mandipropamid and any future test ligand got only bulk
solvent would hand the cognate ligand a structural advantage no test ligand can
have — rebuilding the very asymmetry §41 diagnosed.

### 42a. The pose holds

| system | lig 1st | lig 2nd | lig max | K59 salt bridge |
|---|---|---|---|---|
| aba_WT_s0 | 1.01 | 1.31 | 2.49 | 89% |
| aba_WT_s1 | 1.76 | 2.79 | 3.52 | 100% |
| aba_WT_s2 | 1.90 | 1.93 | 3.18 | 60% |
| *aba, implicit* | *1.09* | *3.18* | *4.88* | *100%* |

Max excursion drops from 4.88 Å to 2.5–3.5 Å. Note the salt bridge is now
sometimes **water-mediated** (60–100%) rather than permanently direct — implicit
solvent had over-stabilised the direct contact, exactly as the missing-water
diagnosis predicts.

### 42b. It reverses §40's reasoning, not its verdict

| arm | 50 ps | 250 ps | **explicit** | closer to |
|---|---|---|---|---|
| ABA | −32.74 | −25.00 | **−32.46 ± 2.72** | **50 ps** |
| mandi | −36.33 | −38.00 | **−46.66 ± 1.29** | 250 ps |

§40 assumed the 250 ps number was the better-converged one. **It was not.** The
250 ps "convergence" was the ligand drifting 4 Å out of a pose only explicit water
can hold — motion *away* from the answer. Length was never the variable; implicit
solvent is unreliable for this anion at every length tested.

**The retraction still stands.** Not because 250 ps was right, but because no
implicit number here is quotable, and the *mutants* have not been re-run in
explicit solvent. Nothing licenses reinstating the K59 flip.

### 42c. The score still does not converge

| system | 1st half | 2nd half | drift |
|---|---|---|---|
| aba_WT_s0 | −35.33 | −33.06 | +2.27 |
| **aba_WT_s1** | −38.92 | −28.83 | **+10.09** |
| aba_WT_s2 | −27.54 | −31.11 | −3.57 |
| mandi_WT_s0/1/2 | | | +0.52, −2.63, +0.01 |

Mandipropamid is stable. ABA is not — one seed moves 10 kcal/mol inside its own
10 ns, and the drifts now point in *inconsistent directions*, which is scatter
rather than the systematic march seen implicitly.

**Between-seed sd: ABA 2.72, mandipropamid 1.29 kcal/mol.** This is the honest
error bar, and the first one this project has had — neither the 8 byte-identical
repacks (sd 0.00) nor the 250 ps block SE (0.3–1.3) could have produced it.

### 42d. What it would cost to make this work

A ddG needs two means, so one replicate each carries **±3.85 kcal/mol**. The
ranking must resolve **0.01–1.8**.

| target precision on a ddG | replicates per variant |
|---|---|
| ±1.0 kcal/mol | 15 |
| **±0.5 kcal/mol** | **60** |
| ±0.25 kcal/mol | 238 |

At n=60: 2760 runs × 38 min = **1748 GPU-hours ≈ 18 days** on the 4-GPU cap, for
one ligand pair.

**The pose was never the binding constraint — the variance is.** Explicit solvent
fixed the thing §41 diagnosed and the method is still not precise enough, which is
a cleaner negative result than either implicit run could have given.

---

## 43. Retiring MM-GBSA, and what replaced it (2026-08-21)

§42d priced MM-GBSA honestly: n=60 replicates per variant per ligand to reach
±0.5 kcal/mol. For a screen of many substitutions × many ligands — which is the
actual goal — that does not scale. **Retired as the ranking method, explicitly not
ruled out**: its rank-based use survives (§40c) and it stays available as a coarse
filter.

### 43a. Non-cognate MD, relaunched (27697950)

Nine systems already built and never run: Imperatorin, Flutamide, α-Estradiol —
matched to ABA at 19–20 heavy atoms and MW 270–276, spanning furanocoumarin,
nitroaromatic anilide and steroid — each from three independently docked poses.

Two changes:

- **Dropped the S9_apo_closed legs.** That cell has since been answered by three
  300 ns replicates in the md191 tree (§29), and re-running it here would have
  used the *older* `data/md` build (NaCl, 178 residues), so it would not have been
  comparable to the answer already in hand.
- **300 → 150 ns.** The original length was sized to catch **gate opening**, but
  §24 and §29 showed open and closed are both kinetically trapped and apo-closed
  does not open either — opening is not observable for *any* ligand on this
  timescale. The ligand-discriminating readout is **pose stability**, which
  separated ABA from mandipropamid within 10 ns in §41–42. 150 ns is 15× that.

### 43b. TI: a different estimator, not a bigger one

The instinct after §42 is "better solvent model, or more replicates." Both treat
the variance as something to be overwhelmed. It isn't — it is a property of the
estimator. MM-GBSA forms a ddG by subtracting **two separately-computed absolute
energies**, so every error in the ~43,000 unchanged atoms enters twice and cancels
only to the extent that two independent averages happen to agree. That is exactly
how one unconverged WT–ABA reference poisoned all 22 ddGs (§40a).

TI computes the same quantity as a **difference along a path**. The mutation is
run in both legs of a thermodynamic cycle, so the unchanged atoms cancel *by
construction*:

```
   WT·L  ──ΔG₁ (mutate in complex)──>  MUT·L
    │                                    │
 ΔG_bind(WT)                       ΔG_bind(MUT)
    │                                    │
   WT + L ──ΔG₂ (mutate in apo)───>  MUT + L

   ΔΔG_bind = ΔG₁ − ΔG₂
```

And for **selectivity the apo leg cancels exactly**:

```
ΔΔG(mandi) − ΔΔG(ABA) = ΔG₁ᵐᵃⁿᵈⁱ − ΔG₁ᴬᴮᴬ
```

so ranking a mutation between two ligands needs **two legs, not four** — the
selectivity is the more trustworthy number by four error sources.

**Cost**, corrected from an earlier claim in this project that FEP was "far
costlier": ~40 ns per leg, against MM-GBSA's ~1200 ns per variant at usable
precision. TI is **cheaper and rigorous**; the real cost is setup, not compute.

### 43c. Setup decisions, each measured rather than assumed

| decision | why |
|---|---|
| **Full sidechain dual topology**, not a common core | 11 of the 15 atoms Val and Ile nominally share change partial charge, and their CG1 hydrogens do not sit at matching coordinates (methyl vs methylene). A common-core mask would have been quietly wrong |
| **Identity asserted at build** | native 81 → seq **79**, 59 → seq 59. Sequential 81 is **also** a valine (native V83), two positions away |
| **tiMerge tolerance set from a measurement** | It refused at its 0.01 Å default. Rather than raise it until it passed, `86` audits first: exactly **one rebuilt hydrogen at 0.024 Å, zero heavy atoms displaced**, in all six legs. Heavy-atom displacement is a hard failure |
| **Masks read from tiMerge, never retyped**, then range-compressed with an atom-count assertion | Amber echoes mdin at fixed width; a long namelist string is an avoidable place for a mask to be clipped |
| **`noshakemask` over the TI region** | Those bonds are being alchemically changed and must not be constrained to an endpoint geometry |
| **12-point Gauss–Legendre** | Nodes never touch λ = 0 or 1, so the softcore endpoint singularity is avoided by construction rather than extrapolated away |
| **One window tested before the other 71** | Caught that `ifsc=1` requires `ntmin=2` — an 8-second failure instead of 72 |

Measured throughput: **271 ns/day under TI** vs 465 plain (~40% alchemical
overhead), ~24 min per window, ~7.2 h for all 72 at the 4-GPU cap.

### 43d. The pre-registered test — ⚠ WITHDRAWN BY §44b

**V81I and K59R are both substitutions in the 4WVO mandipropamid sensor, so both
selectivities should come out negative.** MM-GBSA put them at **+1.36** and
**+1.31** — both the wrong sign, at ±3.85.

If TI recovers the sign, it replaces MM-GBSA for ranking and scales to the real
screen. If it fails, **read the per-window drift the aggregator prints before
blaming the method** — that diagnostic is precisely what MM-GBSA never had, and
is why an error bar there misled twice (§40b, §42c).

Aggregate with `python3 scripts/88_ti_aggregate.py`; it tolerates partial data and
reports N/12 per leg.


---

## 44. The TI pilot returns, and the pre-registered test does not apply (2026-08-22)

All 72 windows COMPLETED (27698144, 27698165), 12/12 in every leg, no empty
outputs. The headline is easy to state and misleading to stop at:

| | selectivity = ddG(mandi) − ddG(ABA) | verdict |
|---|---|---|
| **V81I** | **+4.20 ± 0.64** | prefers ABA |
| **K59R** | **+3.34 ± 1.77** | prefers ABA |

Both positive — the same sign MM-GBSA got (§40c: +1.36, +1.31), now with error
bars ~6× tighter. Read naively this says TI failed the §43d test exactly as
MM-GBSA did, only more confidently. That reading is wrong on two counts, and
both were found by checking the setup before accepting the number.

### 44a. The first error bar was wrong by ~20×

`88_ti_aggregate.py` originally reported ±0.06 on the V81I selectivity. It
computed `sd/sqrt(n)` with `n = 4004`, and **both halves of that were wrong**:

- pmemd prints every timestep **twice** — one energy block per TI region — with a
  bit-identical `DV/DL` in each (checked: 0 of 2000 steps disagree). The regex in
  87 collected both copies.
- It also swept up the trailing `A V E R A G E S` and `R M S  F L U C T U A T I O N S`
  banners, whose `DV/DL` fields are a mean and a standard deviation, not samples.
  Those are the 4 extras in `4004 = 2000×2 + 4`.

So the real frame count is **2000**, and those frames are heavily correlated:
integrated autocorrelation times reach **τ = 75**, driving `n_eff` as low as
**13** in the worst K59R windows. This is the *third* time in this project a
tight within-run error bar has sat on a mean that was still moving (§40b, §42c).
`scripts/89_ti_reparse.py` now re-extracts from `prod.out` with de-duplication
and correlation correction, and 88 reports three error terms separately —
**stat**, **conv** (first-half vs second-half drift), **quad** — because folding
a drift into a statistical error is precisely how the previous two misreads
happened.

Quadrature is **not** the problem: spline-vs-Gauss-Legendre disagreement is
0.00–0.01 kcal/mol in every leg. 12 nodes resolve this integrand fine.

Convergence is only a problem for K59R. V81I legs drift 0.04–0.60 kcal/mol;
K59R legs drift 0.27–1.45 with `n_eff` 13–25, and the second-half-only estimate
moves K59R **further** positive (+3.34 → +4.20), not toward zero. **K59R is not
converged enough to call. V81I is — and it is confidently the "wrong" sign.**

### 44b. The pre-registered test does not constrain single mutants

§43d asserted: *"V81I and K59R are both substitutions in the 4WVO mandipropamid
sensor, so both selectivities should come out negative."* That inference does
not hold. **PYR1^MANDI = K59R / V81I / F108A / F159L** — a *four*-mutation set
selected together. Nothing about it requires each member, alone in a WT
background, to shift selectivity toward mandipropamid; epistasis is the norm in
evolved multi-mutant receptors.

Two independent pieces of evidence already in this repo point the same way:

- **V81 is second shell.** Measured here: native V81 is **4.85 Å** from ABA and
  makes no direct contact. Native V83 — also a valine — is at **2.99 Å**. §23f
  had already noted V81I "barely clashes."
- **LigandMPNN independently called V81I inverted** at **−0.841** (§23h),
  i.e. actively anti-correlated with mandipropamid, by a completely different
  method.

So a well-converged TI, a sequence model, and the structure all agree that V81I
*alone* is not mandipropamid-favouring. That is evidence the **premise** was
wrong, not the estimator. I wrote that premise, and it should have been
challenged when it was written rather than after it returned an inconvenient
answer.

### 44c. The mandipropamid legs start from a physically impossible pose

The more serious defect. `85_ti_build.sh` builds the mandi legs by loading
`3UZ.mol2` into the **WT** protein — F108 and F159 still present — with no pose
relaxation. Measured against the unmutated copy:

| contact | distance |
|---|---|
| **F108 – mandipropamid** | **0.62 Å** |
| F159 – mandipropamid | 1.38 Å |
| V81 – mandipropamid | 1.40 Å |
| heavy-atom contacts < 2.0 Å | 8 (V81I) / **18** (K59R) |

A 0.62 Å heavy-atom separation is two atoms essentially superimposed.
Minimisation relieved it — final poses are ligand-RMSD 1.1–3.9 Å from start and
stay bound (58–101 contacts within 4 Å) — but "relaxed into whatever local
minimum the WT pocket allowed" is not a validated pose, and **every selectivity
number depends entirely on that leg**.

*(Correction to my own first pass: I initially measured ligand displacement in
the raw frame and read 6–18 Å, i.e. dissociation. That was whole-box drift — the
same no-superposition mistake §24 already logged for gate/latch RMSD. After CA
superposition nothing dissociates in any leg.)*

### 44d. What the pilot did and did not establish

**Established.** The machinery works end to end: 72/72 windows, correct residues
(sequential 79 = native 81, confirmed against ligand contacts rather than by the
"is it a VAL" assertion in 85, which native 83 would also pass), correct ligands
(38 atoms/−1.001 ABA, 51/−0.003 mandipropamid), ligands stay bound, quadrature
adequate, ~21 min per window. Cost is ~40 ns per leg against MM-GBSA's ~1200 ns
per variant, and the apo leg cancels exactly in selectivity.

**Not established.** Whether TI can rank selectivity — because the test it was
given cannot answer that, and the arm it depended on started from a 0.62 Å clash.

Recorded deviation: these systems were **neutralised only** (10–11 K⁺, 0 Cl⁻),
not the canonical 0.15 M KCl of §27. It largely cancels in a difference, but it
is a deviation.

### 44e. The test that would actually decide it

§13e already specifies the right calibration and it has never been run for TI:
**PYR1^MANDI + mandipropamid must rank above WT + mandipropamid.** Run the
*quadruple* K59R/V81I/F108A/F159L, ABA vs mandipropamid — 2 legs × 12 windows =
**24 windows, ~2.5 h at 4 concurrent**, since the apo leg cancels for
selectivity. Critically, build the mandipropamid end from **4WVO's own
coordinates**, protein and ligand together, so the bound state is a real
crystallographic pose rather than a ligand dropped into a pocket that cannot
hold it.

- If the quadruple comes out **negative**: TI works, and §43d was simply the
  wrong test.
- If the quadruple also comes out **positive**: TI fails the project's own
  calibration anchor and is disqualified for ranking, exactly as MM-GBSA was.

Either way it is decisive, which §43d was not.

---

## 45. Non-cognate MD: the closed-state filter is ligand-sensitive, and still lets a true negative through (2026-08-23)

Nine 150 ns runs finished (27697950): Imperatorin, Flutamide, α-Estradiol —
matched to ABA at 19–20 heavy atoms and MW 270–276, spanning furanocoumarin,
nitroaromatic anilide and steroid — each from three **independently docked**
poses, against S2 (ABA) as reference on a byte-identical receptor. 1.35 µs total.

Numbering was verified by sequence in all twelve units before anything was
measured: these are 178-residue systems, so native gate 85–89 is sequential
**82–86** (`SGLPA`) and native latch 115–117 is **112–114** (`HRL`). Ligand RMSD
is measured **after** superposing on the protein core — without that it reports
whole-box tumbling, the mistake §24 and §44c each made once.

### 45a. Pose stability separates two of three, but not the third

Run-level means (Å), range across the three runs:

| ligand | ligand RMSD | gate RMSD | verdict vs ABA |
|---|---|---|---|
| **ABA** (cognate) | 1.67 [1.23–2.30] | 1.79 [1.48–2.21] | — |
| Imperatorin | 3.18 [2.23–4.09] | **2.81 [2.59–3.06]** | gate SEPARATE |
| Flutamide | 2.36 [1.58–3.16] | **2.86 [2.21–3.58]** | gate SEPARATE |
| **α-Estradiol** | **1.54 [1.23–2.00]** | **1.69 [1.34–1.89]** | **indistinguishable** |

The ordering is unchanged at every discard window from 0 to 75 ns, so it is not
an artefact of where equilibration was cut.

### 45b. What this does and does not license

**The filter is NOT purely conformation-reporting.** That was the disqualifying
outcome §27 was built to detect, and it did not happen: the gate moves ~1.0 Å
more around imperatorin and flutamide than around ABA, with non-overlapping
run-level ranges. So closed-state stability does carry *some* ligand information.

**But it passes α-estradiol**, which is a genuine non-cognate for WT PYR1 — it
needs F159V/V163W/V164G to work as a sensor. On every observable α-estradiol is
indistinguishable from the cognate ligand: ligand RMSD 1.54 vs 1.67, gate 1.69 vs
1.79, exact p = 0.90. For a design filter that is the dangerous direction of
error: it would wave a non-binder through.

**The design cannot reach conventional significance, by construction.** With
n=3 vs n=3 there are only C(6,3) = 20 splits, so the smallest attainable
two-sided p is **2/20 = 0.10**. Both "SEPARATE" results sit exactly at that
floor. "The ranges do not overlap" reads as strong and is p = 0.10 — worth
stating plainly, because this project has mis-stated a discrimination statistic
before (§15).

**Effective sampling is far smaller than it looks.** 15,000 frames per run, but
the integrated autocorrelation time of the ligand RMSD reaches **τ = 2059
frames (20.6 ns)**, giving n_eff as low as **3** in one run and ≤ 68 in all
twelve. Frame-pooled numbers are descriptive only.

Both pre-registered asymmetries (§27c) still bind: **a ligand that stays put has
not been shown to bind** — 150 ns cannot sample µs–ms unbinding, so only release
is informative and none was seen; and ABA's three runs are three *seeds of one
pose* while each non-cognate's are three *different poses*, so ABA's range is the
optimistic one and the contrast is, if anything, conservative against the
non-cognates.

### 45c. One observation worth a hypothesis, not a claim

α-Estradiol is both the ligand WT tolerates best here **and** the one with the
best evolved sensor in the Tian set — clone 124_10 at **1 µM**, against 10 µM for
the imperatorin and flutamide clones. If WT tolerance predicted evolvability that
would be directly useful for picking targets. On three ligands, one of which
carries the entire signal, it is an anecdote. Recorded so it can be tested, not
banked.

---

## 46. The quadruple TI: right sign, no call — and the docked arm fails structurally (2026-08-24)

48/48 windows. The four low-λ windows of both mandi arms had to be rescued with a
descending λ-ladder (§46c); everything else ran in parallel.

### 46a. The crystal arm gives the right sign and cannot resolve it

| leg | dG | stat | conv | quad | τ_max | n_eff min |
|---|---|---|---|---|---|---|
| aba | 14.88 | 1.37 | 1.46 | 0.01 | 175 | 14.2 |
| mandi_xtal | 13.75 | 1.07 | 0.58 | 0.01 | 245 | 10.2 |
| apo | 18.05 | 1.19 | **4.61** | 0.01 | 261 | 9.6 |
| mandi_dock | **−1381.09** | 1.14 | 4.84 | 0.13 | 69 | 36.0 |

> **selectivity, crystal pose = −1.12 ± 2.34** (stat 1.74 | conv 1.57)
> second-half-only −0.10

**Negative is the correct sign** — the quadruple is predicted to shift preference
toward mandipropamid, which is what PYR1^MANDI does. But the error bar is twice
the effect, so this is a **NO CALL** on the §13e test: not a pass, not a failure.
The second-half-only value (−0.10) moves toward zero, which argues against
reading the sign as meaningful on this much sampling.

Two things worked as designed. Quadrature error is 0.01–0.02, so 12 nodes resolve
even this integrand. And the **apo leg's 4.61 kcal/mol drift — the worst of the
four — cancels exactly out of the selectivity**, which is why selectivity carries
conv 1.57 while the absolute ddG_bind values carry ±5. That is the whole reason
selectivity was chosen as the readout (§43b).

Resolving it is affordable but not free: reaching stat = 0.5 needs ~12× more
sampling, i.e. **121 ns/window, 2.9 µs over 24 windows, ~8.3 GPU-days**. Note
that buys down `stat` only; `conv` shrinks only if the windows genuinely
equilibrate.

### 46b. The docked arm did not fail noisily — it failed structurally

−1381 kcal/mol is not a measurement. Its ⟨∂V/∂λ⟩ curve has the wrong *shape*:
where ABA, crystal and apo all rise to +200…+350 near λ ≈ 0.1 and cross zero
around λ ≈ 0.6, the docked arm plunges to **−2675**.

Measured cause, per window, ligand-to-WT-sidechain minimum heavy-atom distance:

| window | crystal arm | docked arm |
|---|---|---|
| λ=0.32 | 3.59 Å, 0 contacts | 2.04 Å, 0 |
| λ=0.56 | 3.22 Å, 0 | 1.72 Å, **11** |
| λ=0.79 | 3.22 Å, 0 | 1.40 Å, **15** |
| λ=0.88 | 3.44 Å, 0 | **1.22 Å, 14** |

**Every one of those 14 contacts is with the WT Phe108 ring** — CZ, CE1, CE2,
CD1, CD2, CG. The docked ligand has collapsed into the volume occupied by the
decoupled WT phenylalanine. At high λ that sidechain is a ghost, so sitting
inside it costs nothing; at intermediate λ the λ-derivative of that overlap is
enormous, and the integral is meaningless.

**This is not bad luck, and it is not really about docking quality.** The pose was
docked into a model of the *quad* pocket, which carries F108A — a cavity. Docking
fills cavities; that is what it is for. So the ligand is placed exactly in the
volume that the alchemical transformation removes. **Any TI from WT to a
cavity-creating mutant, with the ligand docked into that cavity, is
ill-conditioned in dual topology.** The intended production workflow — design an
enlarged pocket, dock a novel ligand into it, compute ΔΔG by TI — walks into this
by construction.

The crystal pose escapes it because real mandipropamid in 4WVO does not fill the
F108 cavity that way.

Consequence: **the "cost of not having a structure" cannot be quantified by this
experiment.** The docked arm does not degrade gracefully into a worse number; it
leaves the domain where the estimator is defined. The −1394 kcal/mol "gap"
printed by the aggregator is an artefact and must not be quoted.

### 46c. The λ-ladder, and what it did and did not fix

The four low-λ windows of both mandi arms died in the parallel run: at λ ≈ 0 the
WT sidechains are fully coupled and F108 sits 0.62 Å from the ligand, and
minimisation cannot escape that — measured, E = 4.5×10⁸ kcal/mol and |F|max =
8.4×10⁶, flat over 10,000 steps. Heating then died with `illegal memory access
… kNLSkinTest`. A second, backbone-only minimisation stage was tried first and
**was the wrong diagnosis**: there is no downhill path out of a 0.62 Å contact.

Seeding each window from the equilibrated structure of the window above
(lam04 → 03 → 02 → 01 → 00) fixed it for the crystal arm. It did **not** fix the
docked arm, and could not have: its seed, lam04, was already contaminated, so
windows 3 and 2 inherited the collapse (−1181, −294) and only by windows 1 and 0
had the ligand relaxed back out (+27, +10).

Also recorded: amber22 `pmemd.cuda` has no kernels for **h100 (gpu11)**, which
fails in seconds with `invalid device symbol` — the same class as the known
blackwell and k80 failures. Working set is a100 + ada6000 only.

---

## 47. Where this project actually stands (2026-08-24)

Written after 49 recorded belief changes, ~6 µs of MD, two retracted results and
one estimator switch. The question asked was: are we at marginal returns, and
should the goal change?

### 47a. The scoreboard

| method | asked to | outcome |
|---|---|---|
| Rosetta FastDesign / ratchet | recover PYR1^MANDI | F108A + F159L only; K59R and V81I missed — **stage 1 does not clear** (§23j) |
| LigandMPNN | same | F108A +0.111, F159L +0.071; K59R missed; V81I **inverted** at −0.841 (§23h) |
| Per-position steric admissibility | generate residue menus | **fails** — 14.7 of 20 residues admissible, ABA/mandi menus at Jaccard 0.84, i.e. not ligand-conditional (§31) |
| **Pairwise enumeration + packing term** | narrow the search | **works** — 3 of 4 ground truth in the top 150 of 13,457, a **90× narrowing**; shrink/grow enrichment 62 % vs 42 % (§32) |
| Electrostatic complementarity | reach K59R | K59 rank 1 for ABA, but K59R still missed — crystal χ3 ≈ 100° is non-rotameric (§33) |
| Coupled moves (Kortemme) | fix K59R by sampling | **0 % K59 retention even with the cognate ligand** → it is scoring, not sampling (§36) |
| MM-GBSA | rank ΔΔG | **retracted** — unconverged reference; ±3.85 needs n=60/variant ≈ 1748 GPU-h (§40–42) |
| MD closed-state filter (§24b) | rank ligands | ligand-sensitive, but **passes α-estradiol**, a true negative (§45) |
| TI | rank selectivity | crystal arm right sign, **no call** at ±2.34; docked arm structurally invalid (§46) |

One clear success. One weak-positive. Seven failures or no-calls.

### 47b. Three failure modes, each recurring

**1. K59R is invisible to every method, always for the same reason.** It does not
clash, so steric methods miss it *by category* (§32d). Its benefit is a bidentate
hydrogen bond whose crystal χ3 is non-rotameric, so rotamer libraries never
propose it (§33b). And ref2015 pays +10.1 REU of desolvation to bury the salt
bridge it correctly rewards (§23j). More sampling does not help — the method
built for exactly this benchmark retains K59 in 0 % of trajectories.

**2. Every affinity estimator's noise exceeds the spacing it must resolve.**
MM-GBSA: ±3.85 against a ranking spacing of 0.01–1.8. TI: ±2.34 against an effect
of ~1. This is not a bug to fix; it is what computing small differences of large
numbers in a flexible protein costs.

**3. Setup defects that return a plausible number instead of an error.** Shifted
PDB columns silently dropped a ligand and voided a whole arm (§23g); a params
parser summed the atom-*type* column and returned 0.000 charge for any ligand;
smina typed mandipropamid as dummy atoms and scored −9.7 for it; a regex
double-counted every frame and made an error bar 20× too small; docking into a
designed cavity broke the alchemical path. **Nearly every serious error in this
project produced a number, not a crash.** The assertion discipline is the reason
they were caught, and it should not be relaxed.

### 47c. Are we at marginal returns?

**For the affinity-prediction line: yes, and the arithmetic is not close.**
Resolving the selectivity of *one* variant pair, with a crystal structure in
hand, is now costed at **8.3 GPU-days** (§46a). A modest screen — 100
substitutions × 3 ligands — is ~2,500 GPU-days. That is not a screen. And §46b
showed the cheaper path is worse than slow: docking into an enlarged pocket puts
the ligand in the volume the alchemical transformation deletes, so the
dock-then-TI workflow is ill-conditioned *by construction*, not by bad luck.

**For the project: no — but only if the goal changes.** The goal was never "compute
ΔΔG". It is a sensor for a novel ligand. Physics-based ranking was one candidate
engine for that, and it has now been tested about as fairly as we can afford.

### 47d. The reframing this argues for

Every real PYR1 sensor in the literature — Park 2015, Tian, Beltran — came from
**directed evolution over a library**, not from a prediction. And the field's own
baseline is low: RbsB, 2 M variants → 1.2–1.5× improvement (§34). If even large
libraries are low-yield, the value computation can add is **making the library
smaller and richer**, not picking the winner.

We already have the one result that does this. §32's pairwise enumeration turned
13,457 candidates into a top-150 containing 3 of 4 ground-truth mutations. If a
Y2H library of 10³–10⁴ is screenable — and this lab runs Y2H — then a
computationally enriched 10³ library beats a random 10⁶ one.

So the proposal is to change the success metric from

> *rank mutations by predicted ΔΔG* — which nothing here achieves to the needed precision

to

> *design a maximally enriched ~10³-variant pocket library for a given ligand*,
> scored by **recall of known sensors in the top N**.

That metric is measurable today, on data already in hand, without a single GPU.

### 47e. The tests we are missing — in priority order

1. **The Beltran 45-sensor table is extracted and has never been used as a
   benchmark** (§31). Forty-five real sensors is a genuine test set, and top-N
   recall against it directly measures the reframed goal. Cheap, CPU-only, and
   the single highest-value thing not yet done.
2. **Tian's 229 failures have never been screened** (§8 step 5). Failures
   discriminate: a method that cannot separate 229 known failures from the
   successes is not usable, however well it recovers PYR1^MANDI.
3. **Nothing has ever been predicted prospectively.** Every result is
   retrospective recovery of a 4-mutation answer we already knew. §32's 3/4 in
   the top 150 could be overfitting to that one set — points 1 and 2 are exactly
   what would show it.
4. **§38's held-out test is gated on a retired method.** WIN 55,212-2 stays
   sealed pending MM-GBSA, which no longer exists as a ranking method. That gate
   should be rewritten against the library-recall benchmark, not quietly dropped.
5. **We never established the target.** What ΔΔG does a working sensor need? Without
   it, "±2.34 kcal/mol" has no pass/fail meaning. Beltran's table carries EC50s
   and can supply the number.
6. **No computational output has reached the bench.** The lab does Y2H. Nothing
   from this project has been handed over.

### 47f. What I would do next

1. **Finish the running extension** (27725295) — it is already paid for and
   settles whether TI is usable *at all* when a structure exists. That is worth
   knowing even if TI never becomes the screening engine.
2. **Benchmark §32's pairwise method on Beltran's 45 sensors and Tian's 229
   failures.** Days of CPU. This is the decisive test of the reframing, and it
   can invalidate it cheaply.
3. **Rewrite §38's gate** against that benchmark.
4. If the benchmark holds, **design a ~10³ library for one target ligand and hand
   it to Y2H.**

The honest summary: the physics arm has produced reliable *negative* knowledge —
open and closed are both kinetically trapped, the gate–latch staple is clamped by
HAB1, closed-state stability is a weak positive and not a gate, K59R is a scoring
failure and not a sampling one — and one reusable *positive* result, the pairwise
narrowing. The negative results are worth having and were expensive to get
honestly. But they say what will not work. The pairwise result is the only thing
so far that says what will, and it has never been tested outside the one case it
was built on. That test should come before any more GPU time.

---

## 48. What chemistry real sensors use, and why a steric library is not enough (2026-08-24)

§47 proposed reframing the goal to library enrichment, with §32's pairwise steric
method as the candidate engine. §32 recovers 3 of 4 of PYR1^MANDI, missing K59R.
The objection that motivated this section: if a whole **category** of substitution
is invisible to a steric objective, 75 % recall is not "good", it is blind to
chemistry that may be required — and a library built that way would exclude it by
construction.

Measured across every characterised sensor set we hold.

### 48a. Steric chemistry is a minority in every dataset

Distinct substitutions, by category:

| set | clones | distinct | steric | CHARGE | polar |
|---|---|---|---|---|---|
| sd03 (mixed screen) | 691 | 144 | 39.6 % | 25.7 % | 34.7 % |
| sd07 (coumarin) | 78 | 32 | 37.5 % | 28.1 % | 34.4 % |
| Beltran-45 (cannabinoid) | 45 | 40 | 50.0 % | 30.0 % | 20.0 % |

**A purely steric objective addresses at most half of the chemistry real sensors
use**, and the fraction is stable across three unrelated ligand classes. This is
the same lesson as §31 in a different form: **volume tells you the position, not
the residue.** At position 59 a volume objective can say "make it smaller"; it
cannot choose among Q/N/S/T/A/M/L, and it cannot propose R at all, because R is
*bigger*.

K59R is not a PYR1^MANDI quirk. Position 59 carries **8–11 % of all
substitutions in every set**, with 11 distinct residues in sd03 alone
(A, D, E, L, M, N, Q, R, S, T, V).

### 48b. Charge chemistry is concentrated, so covering it is cheap

Charge-changing substitutions are 417 of 2426 occurrences (17 %) across the dev
sets, and they are not spread evenly:

| position | WT | n | steric | CHARGE | polar | charge share |
|---|---|---|---|---|---|---|
| 59 | K | 208 | 23 | **185** | 0 | **89 %** |
| 122 | S | 131 | 50 | 46 | 35 | 35 % |
| 167 | N | 117 | 46 | 38 | 33 | 32 % |
| 159 | F | 364 | 240 | 5 | 119 | 1 % |
| 160 | A | 245 | 236 | 0 | 9 | 0 % |
| 120 | Y | 239 | 57 | 1 | 181 | 0 % |
| 163 | V | 133 | 7 | 0 | 126 | 0 % |

**K59 alone accounts for 44 % of all charge chemistry; seven positions account
for 90 %** (59, 122, 167, 94, 141, 81, 164).

So the pocket separates into position *types*:

- **steric positions** (159, 160, 83, 87, 89) — a volume objective is the right tool
- **charge positions** (59, 122, 167, 94, 141) — a volume objective is the wrong tool
- **polar positions** (120 at 76 % polar, 163 at 95 %) — neither

That argues for a **position-typed library**: run the steric method where steric
chemistry dominates, and simply *enumerate the charge/polar classes* where they
dominate. Covering ~90 % of the charge chemistry costs a full charge menu at
seven positions — well inside what real designs already spend, since Tian's own
TSM library offers 11 residues at position 59 and DSM-Hao offers 17.

### 48c. The experimenters already do this, and it validates the rule

Read from the library designs (sd04), position 59:

| library | allowed at 59 | charge-preserving option? |
|---|---|---|
| TSM | A D E L M N Q R S T V | R |
| DSM-Hao | 17 residues incl. H, R | H, R |
| TNTv1 | D M N R | R |
| TNTv2 | D M N R T | R |
| Coumarin | A Q | **none** |
| PFAS | L M N | **none** |

The TNT libraries deliberately span the charge spectrum at 59 — **R** (preserve
+1) and **D** (reverse to −1) in the same menu. Human designers spend their
positional budget exactly where a steric method is blind.

**This also kills a test I was about to pre-register.** The obvious prospective
test — "does PFAS, which carries a carboxylate, retain the K59 charge more than
neutral TNT?" — is confounded: the PFAS library *offers no charge-preserving
residue at 59 at all*, while TNTv1/v2 both offer R. The observed rates (sd07
coumarin 0.372 vs Beltran 0.733) reflect what each library permitted, not what
selection preferred. Any test at this position must condition on the library
design.

### 48d. PRE-REGISTRATION: the PFAS / TNT prospective test

`scripts/106_sensor_benchmark.py` declares a frozen split in its header:

- **DEV** — sd03, sd07, Beltran-45
- **HELD OUT** — **sd08 (TNT, 96 clones)** and **sd09 (PFAS, 245 clones)**

The holdout is split by **ligand chemistry**, not at random: it asks whether a
rule fitted on cannabinoids and coumarins transfers to a nitroaromatic and a
perfluorinated acid. That is the generalisation a design method must make, and
nothing in this project has ever been tested prospectively (§47e).

**Declared contamination.** Library *designs* (sd04) have been read, including the
allowed sets at position 59 and the randomised-position lists for PFAS and TNT —
they appear in the table above. The **recovered clones in sd08 and sd09 have not
been read**, and the script refuses to read them without `--reveal-holdout`.

**The test.** Build a per-position residue menu for TNT and for PFAS from the
ligand structure and WT PYR1 only, by two methods:
  1. **steric only** — §32's overlap-relief objective;
  2. **steric + position-typed charge/polar enumeration** — §48b's rule.
Then score **recall of the substitutions observed in functional sensors** in
sd08/sd09, and report library-menu size as the cost.

**Predicted, before looking:**
- steric-only recall will fall well short of 1.0, and the misses will be
  concentrated in the CHARGE and polar categories rather than spread evenly;
- the position-typed rule will raise recall materially at a menu cost no larger
  than the libraries Tian actually built;
- position 59 will be the single largest source of steric-only misses.

If steric-only recall is already high on the holdout, §48b is wrong and the
simpler method wins. That outcome is what the split exists to allow.

---

## 49. Charge identity is ligand-dependent — but predictable without a pose (2026-08-24)

§48 showed *where* charge chemistry lives. The open question was whether the
*identity* of the residue can be chosen, and the obvious worry was circularity:
picking a charge at position 59 seems to need the ligand's charged group located,
which needs a pose — and §46b measured docking at **8.33 Å** from the crystal
pose, into a pocket that does not exist yet because it is what we are designing.

### 49a. The measurement

Holding the library **fixed** (mut_lib = tsm, identical allowed residues for every
ligand), grouping the 194 screened ligands by 2D chemistry:

| ligand chemistry | ligands | clones | 59-charge substitutions per clone |
|---|---|---|---|
| **acid (anionic)** | 19 | 38 | **0.05** |
| neutral polar (HBD+HBA ≥ 4) | 66 | 211 | 0.17 |
| neutral apolar | 43 | 127 | 0.31 |
| **base (cationic)** | 11 | 27 | **0.78** |

A **~16× spread**, in the direction physical chemistry demands: an anionic ligand
benefits from K59's +1 (a salt bridge — precisely ABA's situation), a cationic
ligand is repelled by it and must remove it. Acidic ligands also relocate their
charge chemistry elsewhere entirely — E94 (8), V81 (5), I110 (4) — rather than to 59.

Robust in both directions rather than driven by one ligand: **17 of 19 acid
ligands never charge-mutate 59**, and the base effect spans 7 ligands (Riluzole
6/6, Cinacalcet 6/6, Guanabenz 2/2). Weaknesses to keep in view: n is small for
acid (38 clones) and base (27, of which Cinacalcet and Cinacalcet HCl are the same
compound, leaving ~6 distinct bases), and formal-charge class is assigned from
SMARTS, not from pKa at assay pH.

### 49b. Why this dissolves the circularity — for charge

The discriminating feature is the ligand's **formal charge class**, which is a
2D property of the SMILES. It needs no pose, no docking, and no model of the
designed pocket. **The circularity is binding for steric placement and largely
absent for charge**, which is convenient, because charge is exactly where the
steric methods were blind (§48a).

This does not extend to the polar positions (120 at 76 % polar, 163 at 95 %).
Hydrogen-bond geometry is directional and short-range, so those plausibly do need
a pose. They are left open here rather than assumed away.

### 49c. Options, rated

| option | cost | potential | success chance | verdict |
|---|---|---|---|---|
| **E. chemical-similarity lookup** — nearest of 194 screened ligands, reuse its menu | very low | high in-domain | high in-domain, unknown out | **do first** |
| **B. charge-class rule from 2D chemistry** (§49a) | low | high — sets identity, not just position | moderate–high | **do second** |
| A. ligand-agnostic fixed charge menu at the 7 positions | ~zero | medium (90 % of charge chemistry) | high | fallback; discards the 16 × signal |
| D. pose-*ensemble* consensus, P(position sees a polar atom) | moderate | moderate | unknown | only if B underperforms |
| C. pose-derived electrostatic complementarity | high | high if it worked | **low** — 8.33 Å docking error, hypothetical pocket | **do not** |

**E has never been done, and it gates the rest.** No trivial baseline has ever
been established in this project. If "find the chemically nearest ligand already
screened and copy its library" recalls as much as anything built from structure,
then every physics result here is competing with a lookup table — and that has to
be known before more effort is spent. It is the cheapest experiment available and
it can invalidate a great deal.

### 49d. PRE-REGISTRATION: the K59 charge-retention test on the holdout

Written before reading sd08/sd09.

**PFAS = perfluorooctanoic acid / PFOS — a strong acid, anionic at assay pH.**
**TNT = 2,4,6-trinitrotoluene — neutral, electron-poor.**

Both libraries randomise position 59, so *retention* is comparable even though the
substitution menus differ (PFAS offers L/M/N, TNTv1-2 offer D/M/N/R/T).

> **Prediction: the fraction of clones retaining a POSITIVE charge at position 59
> — wild-type Lys, or a mutation to Arg/His — will be HIGHER for PFAS than for
> TNT.**

Dev-set anchors: acid ligands charge-mutate 59 in 5.3 % of clones (2/38), neutral
ligands in 22.5 % (76/338). So the expectation is PFAS retention well above TNT's.

This is falsifiable in the obvious way: if PFAS ≤ TNT, §49a does not transfer to
new chemistry and the charge-class rule should not be built on it. It is also the
first genuinely prospective prediction this project has made (§47e).

---

## 50. The lookup baseline: chemical similarity adds almost nothing, but frequency sets a real bar (2026-08-24) — ⚠ METRIC SUPERSEDED BY §51

§49c argued this had to be run before anything structural, because no trivial
baseline had ever been established here. `scripts/107_lookup_baseline.py`,
leave-one-ligand-out over the 89 sd03 ligands with ≥3 clones.

The comparison that matters is not "does lookup work" but **does chemical
similarity beat a menu that never looks at the ligand at all** — common
substitutions are common, and a frequency menu will recall a lot for free.

### 50a. Result

Recall of a held-out ligand's observed substitutions, top-K menu:

| scorer | @5 | @10 | @20 | @40 |
|---|---|---|---|---|
| **freq** (ligand-blind) | 0.148 | 0.236 | **0.348** | **0.517** |
| **sim** (Tanimoto-weighted) | 0.152 | 0.233 | 0.358 | 0.554 |
| 1nn (nearest ligand's menu) | 0.113 | 0.172 | 0.200 | 0.201 |
| *oracle* (ceiling) | 0.567 | 0.855 | 0.991 | 1.000 |

Per-library (target and neighbours sharing allowed residues) gives the same
picture: tsm +0.023, dsm +0.024 at K=20.

**Similarity beats frequency by +0.010 on average.** Head-to-head at K=20:
23 wins, 13 losses, **53 ties** — the top-20 menus are mostly identical. And
**1nn is far worse than frequency** (0.20 vs 0.35), so copying a single
neighbour's library is actively bad.

### 50b. Why it fails — the reference set is chemically sparse

Similarity does help *when a close neighbour exists*:

| | n | sim | freq | delta |
|---|---|---|---|---|
| close neighbour (Tanimoto ≥ 0.5) | 24 | 0.405 | 0.360 | **+0.045** |
| far (< 0.5) | 65 | 0.340 | 0.343 | −0.003 |

But over the 89 ligands, pairwise Tanimoto has **median 0.110**, 90th percentile
0.186, and **only 0.4 % of pairs reach 0.5**. So the method is not wrong in
principle — there is simply nothing to look up. For a genuinely novel ligand,
which is the whole use case, a close neighbour will essentially never exist.

**Option E (§49c) is dead.** Lookup is not a competitive method and should not be
built on.

### 50c. What this bought: the bar, quantified for the first time

The frequency menu is ligand-blind, costs nothing, and is not weak:

> **recall@20 = 0.35, recall@40 = 0.52**

Any structural or rule-based method must beat that to be worth its compute. This
project has run nine methods (§47a) and has never once measured one against a
null. From here, nothing gets reported without it.

Two further things the baseline hands us:

- **The task is well-scoped.** oracle@20 = 0.99, so ~20 substitutions covers
  essentially all of any single ligand's observed chemistry. A menu of 20–40 is
  the right target size, not 200.
- **The frequency menu is 13 steric / 5 polar / 2 charge**, against a true
  population of roughly 40/35/26 (§48a). So the cheap baseline is *itself*
  skewed toward steric chemistry and under-represents charge — which is exactly
  the gap §48b and §49a propose to fill, and gives the position-typed rule a
  specific, measurable thing to beat rather than a vague improvement.

The top-20 ligand-blind menu, for reference: F159I/V/T/L/A/G, V81Y/L/R, Y120A/G,
V163W, V83L, A160L/I/V/M, K59D, L87M, V164F.

---

## 51. The ground truth is positive-unlabeled, and the metric has to change (2026-08-24) — ⚠ LIBRARY SIZES IN 51d ARE WRONG, SEE §53

A correction from Jannis, and it invalidates part of how §50 was framed.

### 51a. What the hit sets actually are

The characterised clones are **hits that were found**, not best sequences, for
three separate reasons:

1. **The landscape is not fully tested.** Each library randomises a subset of
   positions with a subset of residues, so a substitution absent from the data may
   simply never have been offered.
2. **Not every combination reaches screening.** Libraries are *sampled*, sized so
   that a hit is likely to be found if one exists — not enumerated. Absence from
   the hit set does not mean the variant failed; it may never have been made.
3. **Selection removes variants for reasons unrelated to pocket binding** —
   constitutive activity (traditionally an interface effect, not a pocket one), a
   mutation permitting URA3 activity, and sometimes promiscuity. Promiscuity
   filtering is *not* uniform: the PFAS sensors are marketed as pan-PFAS, i.e.
   cross-reactivity between congeners was a feature there, not a reject criterion.

So "not observed" means **unknown**, never "non-functional". This is a
positive-unlabeled sample.

### 51b. Three consequences, two of which §50 got wrong

- **Precision is not measurable.** Only recall. Any metric that penalises a
  proposed substitution for being absent from the hit set is scoring the screen,
  not the biology.
- **The frequency baseline is advantaged by construction.** It reproduces what
  people put into libraries and then found — drawn from the same design process
  that generated the labels. **Beating it is therefore strong evidence; failing to
  beat it is weak evidence of failure.** A structural method proposing something
  never offered scores zero even if it is right.
- **`oracle@20 = 0.99` was misread.** It means 20 substitutions cover the
  *observed* hits, not that 20 suffice.

### 51c. The metric that matches the actual goal

Substitution-level recall quietly asks *did we find the best answer*. The real
question is *can the library window be narrowed while still containing a hit*.
So the primary metric is now:

> **fraction of ligands for which the designed library contains ≥ 1 known hit
> clone, as a function of library size** — where a clone counts only if *every*
> one of its substitutions is in the menu.

Ligand-blind frequency menu, 89 ligands with ≥3 hit clones, 550 clones:

| K | library size | ligands with ≥1 hit | clone-level |
|---|---|---|---|
| 10 | 240 | 29 % | 8 % |
| **20** | **13,440** | **51 %** | 16 % |
| 30 | 737,280 | 58 % | 22 % |
| **40** | **21.0 M** | **72 %** | 33 % |
| 60 | 6.7 B | 88 % | 51 % |
| 144 (all) | 3.2 × 10¹⁹ | 100 % | 100 % |

Ligand-level is the right row to read, and clone-level is the wrong one, for
exactly Jannis's reason: **you need *a* sensor, not every sensor.**

### 51d. Grounded against the libraries actually built

| library | positions | size |
|---|---|---|
| Coumarin | 11 | 138,240 |
| TNTv2 | 14 | 1.24 M |
| TNTv1 | 16 | 4.67 M |
| PFAS | 13 | 49.5 M |
| TSM | 18 | 7.7 × 10¹⁵ |
| DSM-Hao | 18 | 3.8 × 10²¹ |

The **focused** libraries sit at 10⁵–10⁷ — deliberately inside Y2H's ~10⁶–10⁷
screening capacity. The general ones are astronomically larger and are sampled,
not enumerated. Which is itself the argument for this metric: even the
experimenters were never enumerating a landscape, they were sizing a window.

### 51e. The revised bar

> **A ligand-blind frequency library retains a known hit for 58 % of ligands at
> 10⁵·⁹ variants, and 72 % at 10⁷·³ — the same size range as Tian's own focused
> libraries.**

A ligand-aware method has to beat that *at matched library size*, or reach the
same retention in a smaller window. And because of §51b, all of these are
**lower bounds** — a library may well contain a working sensor that was never
found.

The §50 numbers (recall@20 = 0.35) are not wrong, but they answer a question we
no longer care about, and that section now carries a pointer here.

---

## 52. The right incumbent, and what "something novel" would actually mean (2026-08-24) — ⚠ SIZE ARITHMETIC CORRECTED IN §53

Jannis: a designed library is only worth building if it delivers **≥10× size
reduction at comparable hit rate**, or **the same size with more/better hits**.
Otherwise use the DSM library already sitting in glycerol stocks, which is free.

That is correct, and §51 was still benchmarking against the wrong thing.

### 52a. A flaw in the metric, stated plainly

Every ground-truth hit came **out of** the existing libraries. So any menu we
design is a *subset* and can only lose hits, never gain them — the ceiling is
100 % at full size. **The hit-retention metric can never show a designed library
beating the incumbent.** It can only measure how far the window can be narrowed
before hits are lost. That is worth knowing, but it is not evidence of advantage,
and §51 came close to reading it as though it were.

Related: **hit *rate* is not recoverable from this data at all.** The datasets
record characterised hits, not screening depth, so density-per-transformant — the
quantity that would actually decide the economics — cannot be computed here.

### 52b. The incumbent is a two-round process, and round 2 is already focused

| library | actual members |
|---|---|
| Tian coumarin (focused) | **77,327** |
| Tian TNT (focused) | **506,229** |
| Park 2023 | ~12,000 clones, ~400,000 colonies screened |

The design *spaces* quoted in §51d (10¹⁵–10²¹) are not what gets built. The real
process is:

1. screen the existing DSM/TSM stock — free, already made;
2. if that fails, build a **focused ~10⁵ library from sequence profiles of
   round-1 hits**.

And step 2 works: it isolated sensors for 4-methylumbelliferone,
7-methoxycoumarin, TNT, DNT and 2ADNT, **all missed in round 1** (§34d).

So the target is not "beat the DSM stock" — nothing beats free. It is round 2.

### 52c. What would genuinely be novel

**Tian's focused libraries require a round-1 screen to build the profile.** A
method that designs an equally focused (~10⁵) library from ligand chemistry and
structure alone — with no round-1 hits to learn from — **removes an entire
experimental round per ligand.** That is a concrete, economically meaningful
advance, and it is the honest framing of what computation can contribute here.
It is also the one claim not undermined by §52a, because it is not a claim about
beating the library; it is a claim about not needing the first screen.

The second unambiguous target already exists: **the 229 documented failures**
(§13, 28–50 heavy atoms, 15.3–22 Å, screened against all six libraries, zero
hits). A hit on any of those is directly comparable to a published negative.

### 52d. Option B is open, and now quantified

Frequency ordering is **biased toward weak sensors**. At matched potency, the
fraction of clones whose every substitution is in the menu:

| K | library size | 1 µM | 10 µM | 100 µM |
|---|---|---|---|---|
| 20 | 13,824 | 13 % | 10 % | 21 % |
| **25** | **82,944** ≈ Tian coumarin | **17 %** | 13 % | **24 %** |
| 30 | 663,552 | 20 % | 17 % | 29 % |
| 40 | 15.6 M | 33 % | 29 % | 38 % |

Potent clones are captured **less** often than weak ones, and not because they
are more complex — substitutions per clone is 2.53 / 2.57 / 2.44 across the three
potency classes. Their substitutions are simply rarer, and a frequency menu fits
the mode, which is weak hits (316 clones at 100 µM against 60 at 1 µM).

**So there is real headroom for "same size, better hits":**

> at ~80,000 members — the size Tian actually builds — a ligand-blind frequency
> library captures **17 % of the potent (1 µM) clones**. Beat that at matched
> size and the library is genuinely better tuned, not merely smaller.

This is the first target in this project that is both measurable and, if met,
clearly worth doing. It also has an obvious cheap baseline to clear first —
potency-weighted frequency (weight each substitution by 1/min_conc) — which must
be run before any structural method claims the credit.

---

## 53. Correction: the libraries are substitution-DEPTH limited, and round 2 is a change of shape (2026-08-24)

Jannis: DSM-Hao is a **double-substitution** library — a limited number of total
substitutions, combined two at a time. The secondary libraries *limit the total
substitutions but allow more per binder*, so they reach combinations the primary
library cannot express at all, despite offering fewer possibilities overall.

That is correct, and the library sizes in §51d and §52b were badly wrong.

### 53a. What I got wrong

I computed library size as the product over positions of (allowed residues + 1),
which assumes every position can vary simultaneously. The primaries are
depth-limited, so that is not the space they span. Confirmed directly from the
clone data:

| library | substitutions per clone |
|---|---|
| **dsm** (primary) | mode **2** — 193 of 266 clones |
| **tsm** (primary) | mode **3** — 346 of 403 |
| sd07 coumarin (secondary) | mode **7**, range 3–10 |
| sd08 TNT (secondary) | mode **7**, range 1–10 |

Recomputing under the real constraint:

| library | positions | depth | **true size** | what §51d claimed |
|---|---|---|---|---|
| **DSM-Hao** | 18 | **2** | **36,140** | 3.79 × 10²¹ |
| **TSM** | 18 | **3** | **332,863** | 7.70 × 10¹⁵ |
| Coumarin | 11 | all 11 | 138,240 | 138,240 |
| TNTv2 | 14 | all 14 | 1,244,160 | — |
| TNTv1 | 16 | all 16 | 4,665,600 | — |
| PFAS | 13 | all 13 | 49,545,216 | — |

Wrong by up to **16 orders of magnitude**. Every library is 10⁴–10⁷, and the
**primaries are the smallest of them**.

### 53b. What this changes

**Round 2 is not smaller. It is differently shaped.**

- round 1 — **broad and shallow**: 18 positions, 2–3 substitutions per variant
- round 2 — **narrow and deep**: 11–14 positions, 7–8 substitutions per variant

Round 1 asks *which positions and residues matter*. Round 2 asks *which of them
combine*. The 7–8-deep combinations that make the secondary sensors work are
**not reachable in the primary library at any screening depth**, because it
cannot express them.

So two earlier framings are now void:

- **"≥10× library-size reduction" is the wrong axis** (§52). Round 2 is
  the same size or larger than round 1. Size is not what round 2 buys.
- **The design task is not "make it smaller"** — it is *which few positions and
  residues are worth combining deeply*.

What survives from §52 is the claim that matters: **Tian's round-2 libraries are
built from round-1 hit profiles, so producing one without round 1 removes an
experimental round.** And §52d's potency target survives untouched, since it is
about which clones a menu captures, not about library size.

Also unaffected: the **hit-retention percentages** in §51c are set containment
and do not depend on library size. Only the size axis beside them was wrong.

### 53c. The task, as actually posed

> Build a pipeline that produces focused **secondary** libraries for a range of
> molecules.

Since every characterised hit so far comes from this limited context, continuing
in it is legitimate. And a method that recovers the known substitutions well
*and* proposes extra ones lets those extras be considered for the library too —
the proposal set is not capped by what has been tried.

**The retrospective form of the test**, which we can run now: take the DSM-Hao
**round-1** hits for PFAS / TNT / coumarin, design a secondary library from them,
and ask whether it recovers the known **round-2** hits — in a smaller library
than the one actually built.

⚠ **Bias control is the hard part here, not the method.** Round-1 and round-2 hits
for the same ligand share chemistry by construction, so groups must be held out
deliberately — by ligand and by chemical class — or the evaluation will report
memorisation. §48d's frozen PFAS/TNT split is the starting point but is not
sufficient on its own, because round-1 hits for those same ligands would leak.

### 53d. Standing note on prospective targets

Jannis has selected an initial target molecule and others with known weak initial
hits, and is **deliberately withholding their identities to prevent unintentional
biasing**. Do not ask for them, and do not attempt to infer them from the data.
Any method must be specified and frozen before those identities are revealed.

---

## 54. Round-1 hits do not predict round-2 hits — the vocabulary does (2026-08-24)

The obvious model of the two-round process is "extrapolate from this ligand's
weak round-1 hit". Measured on the coumarin family (sd07 is DEV; PFAS/TNT stay
sealed), that model is wrong.

### 54a. Per-ligand carryover is near zero; the pooled vocabulary is near complete

For each round-2 sensor, what fraction of its substitutions had already been seen
in round-1 hits **for that same ligand**, versus in the **pooled** round-1 data
across all 194 ligands:

| ligand | round-2 distinct | own round-1 | **pooled round-1** |
|---|---|---|---|
| Osthole | 14 | 64 % | 100 % |
| Imperatorin | 13 | 46 % | 100 % |
| Isopsoralen | 11 | 36 % | 100 % |
| Methoxsalen | 13 | 23 % | 100 % |
| Scopoletin | 16 | 6 % | 100 % |
| Psoralen | 17 | 6 % | 94 % |
| Bergapten | 16 | **0 %** | 94 % |
| Citropten | 13 | **0 %** | 100 % |
| 7-Methoxycoumarin | 15 | **0 %** | 93 % |
| 5,7-DH-4-methylcoumarin | 11 | **0 %** | 91 % |
| 4-Methylumbelliferone | 20 | **0 %** | 75 % |

**Four of eleven round-2 sensors had no round-1 hit for their ligand at all** —
they were round-1 *failures*, rescued by round 2. And even where a round-1 hit
existed, it typically supplied under a quarter of the chemistry the round-2
sensor used.

Meanwhile round-1 across all ligands explores **144 of the 342 possible single
substitutions** (42 %) at the 18 pocket positions, and **75 % of all round-2
chemistry already appears somewhere in it**. The residual 8 substitutions are at
positions **65, 71, 74, 109, 124, 134, 178, 184** — *outside* the original 18, so
round 2 also widened the position set.

### 54b. What the design problem actually is

Not extrapolation from a weak hit. **Subset selection from a known, small, largely
closed vocabulary**, plus a decision about depth:

> given ~144 substitutions that are known to do something in PYR1, choose ~11
> positions and their residues, to be combined **7–8 deep**, such that the
> resulting ≤10⁵ library contains a sensor for the target chemical class.

That is a far better-posed problem than anything attempted in §20–46, and it
explains why Tian's method works: their profiles are built from first-round hits
**pooled across the library**, not from one ligand's hits.

Two consequences for a prospective target:

- **An initial hit is not necessary** — 4 of 11 succeeded without one.
- **An initial hit is not sufficient** — own-ligand carryover is 0–64 %. It should
  be used to *weight* the vocabulary, never to restrict it. Restricting to the
  hit's own substitutions would have failed for every ligand in the table.

### 54c. The next test, and it is runnable now

Design a **coumarin-class secondary library** from the pooled round-1 vocabulary
plus class chemistry — **without touching sd07** — and ask:

1. does it contain the known round-2 hits, and for how many of the 11 ligands?
2. at what library size, against the **138,240** Tian actually built?
3. does it recover the four ligands that round 1 missed entirely — the cases where
   the method has to work without a hit to lean on?

Bias control: pooled round-1 is legitimate input, because it is what Tian had.
sd07 is not. PFAS and TNT stay sealed for the prospective test (§48d).

---

## 55. Co-folding validation, and whether the data can support an ML design tool (2026-08-24)

### 55a. The co-folding test (27727052, running)

§49c rated pose-derived methods do-not-attempt, on the measured basis that
docking into a not-yet-existing pocket lands 8.33 Å from the crystal pose (§46b).
**An initial hit changes that problem**: it supplies a sequence *known to bind*, so
the task becomes modelling a complex that exists rather than designing one. Boltz
2.2.1 is installed locally, so this costs nothing but GPU time.

Four runs, and the controls carry the weight:

| run | purpose |
|---|---|
| PYR1^MANDI + mandipropamid | the design scenario |
| PYR1^MANDI + HAB1 + mandipropamid | matches 4WVO, which is ternary |
| **WT + mandipropamid** | **negative control** |
| **WT + ABA** | **positive control** (3QN1) |

⚠ **Corrected (Jannis).** An earlier draft said the gate–latch only closes on
HAB1 binding. That is wrong: **ligand binding closes the gate; HAB1 then locks the
closed form together.** So the *binary* run can reach the closed state on its own,
which makes it the more informative of the two rather than a degraded version of
the ternary one — and it is also the configuration a design pipeline would
actually have.

The negative control carries the most weight: **WT is known not to bind
mandipropamid** — that is why the quadruple was evolved — so if WT scores as
confidently as the quadruple, the method has no discriminative power however good
the quadruple's pose looks.

Gate: if 4WVO is reproduced, structure becomes usable for the steric positions
(159/160/83/87/89) and the polar ones (120/163) that §49b left open. If not, the
stack stays pose-free and those positions fall back to frequency priors.

### 55b. Is there enough data to bias an ML tool from SMILES to sequence?

Counted rather than guessed:

| quantity | value |
|---|---|
| distinct **ligands** with ≥1 hit | **~208** (194 sd03 + 14 Beltran) |
| labelled clones | ~1150 |
| distinct substitutions (**output space**) | **144** over 18 positions |
| documented **ligand-level negatives** | **229** |
| ligands with ≥5 clones | 47 (69 have exactly 1) |
| median pairwise ligand Tanimoto | **0.106** |

**The effective sample size is ~208, not ~1150.** Clones sharing a ligand are not
independent draws for learning a ligand→sequence map; they are repeats under one
input. So:

**Generative SMILES → sequence: no.** Two to four orders of magnitude short, and
the chemical space is too sparse (median Tanimoto 0.106) for interpolation to
rescue it — the same sparsity that killed the lookup baseline in §50b.

**Fine-tuning LigandMPNN: not advisable, for two reasons.** It consumes a
*structure*, not a SMILES, so it inherits the pose problem this was meant to avoid.
And §23h measured the base model on exactly this task: it recovers F108A (+0.111)
and F159L (+0.071), misses K59R, and **inverts V81I to −0.841**. Fitting 208
ligands on top of that would overfit long before it corrected the inversion.

**What the data does support** is a small, heavily regularised conditional model —
position × coarse residue class (~18 × 5 = 90 cells) conditioned on a handful of
ligand descriptors (formal charge class, heavy-atom count, HBD/HBA, aromatic
rings). That is §49a generalised from one position to all of them, and 208 ligands
can carry it if the parameter count stays near 10².

**And one genuinely well-posed ML task is sitting unused.** The 229 failures are
*ligand-level* labels — no sensor found for that molecule in any of six libraries.
Combined with the 208 successes that is **437 ligands with binary labels**, enough
for a modest classifier on ligand descriptors answering:

> *is this molecule tractable for a PYR1 sensor at all?*

That is target triage, not design, and it is the one place where the label count
is adequate rather than marginal. It has never been attempted here, and it would
directly inform which molecules are worth a library.

---

## 56. Graft status, and the donor pockets measured for the first time (2026-08-24)

### 56a. Where the graft arm stands

`79_graft_geometry.py` screened 19 candidates for whether PYR1's 19-residue
machinery (four discontinuous segments, §9c) could be carried onto a large-cavity
SRPBCC relative. **Geometry only — no design was attempted, and the arm has been
parked at lower priority since.** Its conclusion stands: *rebuild, not transplant.*

### 56b. The tension is categorical, not gradual

Across all 19 candidates cavity vs core RMSD gives r = +0.10 — no trend. The
structure is a clean split:

| | cavity Å³ | core RMSD | fident |
|---|---|---|---|
| **big-cavity** (2PCS, 2NS9, 2BK0, 6AWV) | 319–570 | **3.45–4.03 Å** | **0.098–0.122** |
| **transplantable** (4DSB, 3OQU) | **162–204** | **1.03–1.46 Å** | 0.457–0.511 |

The transplantable ones are close PYR1 homologs (≈50 % identity) whose pockets are
**the same size as PYR1's 174 Å³** — no gain. The big-cavity ones are distant
relatives at ~10 % identity, and 3.5–4 Å of core RMSD is a rebuild. There is
nothing in between in this candidate set.

### 56c. The donor pockets, measured (`109_graft_pocket_landscape.py`)

79 measured coverage *from PYR1's side* — how many of PYR1's 19 machinery
residues have an equivalent position in the donor. It never asked what the
**donor's own** pocket looks like. Measured now, by the same criterion used for
PYR1 (side-chain heavy atom within 4.5 Å of the bound ligand):

| PDB | cavity Å³ | ligand | pocket residues | **side-chain lining** | Å³ per lining residue |
|---|---|---|---|---|---|
| **2PCS** (CoxG) | 570 | UNL, 30 heavy | 30 | **27** | 21 |
| 6AWV | 319 | 28E, 42 heavy | 26 | **21** | 15 |
| 3TFZ | 207 | CXS, 28 | 18 | 17 | 12 |
| 3OQU | 204 | ABA, 19 | 17 | 17 | 12 |
| 4DSB | 162 | ABA, 19 | 17 | 17 | 10 |
| **PYR1** | **174** | ABA, 19 | — | **19–20** | ~9 |
| 2NS9 / 2BK0 / 3QRZ | 455 / 344 / 210 | **apo** | not measurable | — | — |

**Cavity volume and the number of lining side chains track each other almost
exactly: r = +0.99 (n = 5).** So the big donors are not "PYR1's pocket with more
room" — they are pockets with **more positions**. 2PCS lines its cavity with **27
side chains against PYR1's 19–20**, roughly 40 % more design positions, in a
scaffold sharing 9.8 % sequence identity. Volume per lining residue also rises,
from ~9–12 Å³ in the PYR1-like pockets to 21 Å³ in 2PCS, so it is a more open
architecture as well as a larger one.

⚠ **The second and third largest candidates are apo.** 2NS9 (455 Å³) and 2BK0
(344 Å³) have no bound ligand, so their volumes come from cavity detection rather
than ligand contact and are not strictly comparable to PYR1's 174 Å³, and their
lining-residue counts cannot be obtained at all. **2BK0 was the candidate
previously flagged as most promising** (19/19 covered, gate RMSD 2.20 Å); that
recommendation rests on a cavity number of a different kind from the one it was
being compared against.

### 56d. No, we still cannot simulate the switch — and for grafts that is the binding constraint

Confirmed, and 79's own header says so: it screens what is *geometrically*
possible and "says nothing about whether the graft would still cycle between open
and closed."

§24 and §29 are why. Open and closed are **both kinetically trapped**: neither
state converts even once in **1.8 µs of aggregate WT sampling**, and apo-closed
holds for 3 × 300 ns with **zero crossings**. So the accessible observable is
*stability of a state*, never *transition between states*.

For the pocket-expansion arm that is a limitation. **For the graft arm it is
disqualifying at the validation step**, because a graft's entire value
proposition is that the transplanted machinery still cycles. What we can measure —
does the closed state hold — is satisfied equally well by a **constitutively
closed** protein, which is a failure mode, not a success. §45 sharpens this: the
closed-state filter is ligand-sensitive but passed α-estradiol, a true negative.

So a graft could be designed, folded, and shown stable, and none of that would
distinguish a working switch from a dead one. Any serious graft attempt needs the
switch assayed experimentally (Y2H reports exactly this), or an enhanced-sampling
method this project has not built.

---

## 57. How to test expanded pockets without running library × library (2026-08-24)

The concern: screening a pocket library against a chemical library is a product,
and characterising interactions afterwards is worse. Grafting from a donor whose
ligand is already known collapses one dimension — you know what to test.

That reasoning is sound. The premise mostly fails on the donors we have.

### 57a. The donors' "known ligands", checked

| donor | cavity Å³ | het code | what it actually is |
|---|---|---|---|
| **2PCS** | **570** | UNL | **"Unknown ligand"** — unmodelled density |
| 2NS9 | 455 | — | apo |
| 2BK0 | 344 | — | apo |
| **6AWV** | 319 | 28E | **(−)-epicatechin**, C15H14O6 — a real ligand |
| 3TFZ | 207 | CXS | **CHES buffer** — a crystallisation additive |

**Exactly one big-cavity donor has an identified physiological ligand.** The
570 Å³ champion's ligand was never assigned; 3TFZ's is buffer; two are apo. And
6AWV is the weakest of the big four on transplantability — core RMSD 3.93 Å,
9.8 % identity, and **17/19** machinery coverage, the only one of the four that
loses machinery residues.

So grafting buys the known-ligand advantage in one case, and that case is the
hardest graft.

### 57b. The product only exists if you are doing discovery

Library × library is a **discovery** cost — it applies when you are looking for
*new* pocket/ligand pairs with no starting point. It is not incurred when a weak
hit already exists for the target, because then the chemical dimension is already
collapsed to one compound and the problem is **one library × one analyte**, which
is exactly the round-2 optimisation of §54.

Note also that the product is smaller than it looks even for discovery: Y2H
selection is parallel over variants, so the unit is *one screen per compound*, not
one assay per (variant, compound) pair. Tian ran ~194 compounds this way. The cost
is N screens, and N is what wants reducing.

### 57c. Five ways to avoid the product, ranked by what they cost

1. **Start from a weak hit.** No chemical dimension at all. Available today for
   targets that already have one, and §54 says the design problem there is subset
   selection from a ~144-substitution vocabulary — the best-posed problem in this
   project. *Cost: one library per target.*
2. **Tractability triage** (§55b). 208 positives + 229 documented negatives = 437
   ligand-level labels, enough for a classifier that predicts whether a compound
   is worth screening at all. Cuts N before any bench work. *Cost: days of CPU,
   never attempted.*
3. **Pooled screening with deconvolution.** Screen against pools of 10–20
   compounds, deconvolute the survivors. Cuts N by the pool size. The standard
   risk is a promiscuous variant dominating the pool — worth noting that
   promiscuity is not always a defect here, since the PFAS sensors are marketed as
   pan-PFAS (§51a). *Cost: bench, standard practice.*
4. **Probe-set screening.** A small chemically diverse set (~20–50 compounds) run
   against the expanded pocket to answer *does this pocket bind anything and still
   switch* before committing to a specific analyte. Separates "is the pocket
   functional" from "does it bind X", which §56d says we cannot answer
   computationally. *Cost: one screen against a fixed set, reusable across designs.*
5. **Grafting.** Collapses the chemical dimension only where the donor ligand is
   real — one candidate (6AWV/epicatechin), which is also the hardest graft. Plus
   §56b (transplantable donors are PYR1-sized; big-cavity donors need a rebuild)
   and §56d (the switch cannot be validated computationally, and a constitutively
   closed graft passes every filter we have). *Cost: highest, and it carries the
   one risk we have no assay for.*

### 57d. What I would actually do

**If the target already has a weak hit, grafting is solving a problem you do not
have.** The chemical dimension is already collapsed; the remaining problem is
library design for a known analyte, which is §54's vocabulary-subset problem and
is testable retrospectively today on coumarin.

Grafting keeps one distinct use: as a **positive control that a large pocket can
be made switchable at all**. 6AWV/epicatechin is the only candidate where that
control would have an identified ligand to report against. That is a real
experiment, but it answers a feasibility question rather than producing a sensor,
and it costs a rebuild.

The step that removes the most future cost for the least effort is (2): a
tractability classifier on 437 labelled ligands. It is the only item here that
shrinks N *before* anything reaches the bench, and the data for it has been
sitting unused since §31.

---

## 58. The under-filled pocket problem, and a correction to §30's status (2026-08-24)

Jannis raised the objection that is prior to every scoring question in this
project:

> a weak hit found *without* an expanded pocket cannot test an expanded pocket,
> because expansion introduces residues the library never offered — and if the
> cavity is enlarged while the ligand stays the same size, the ligand no longer
> reaches the gate. Either the gate cannot close (no signal), or it closes
> anyway (constitutive), unless the cavity can hold open without a large energy
> cost.

This is correct and it is mechanistic, not logistical. **The gate closes onto the
ligand.** Enlarging the cavity without enlarging the ligand breaks the
transduction step, in one of two directions, and both are failure modes —
constitutive closure is precisely what the counter-selection removes.

### 58a. What we already know, and why it does not settle it

§29 measured exactly the relevant cell on the old tree at n = 3: **apo-closed
PYR1 does not open in 3 × 300 ns, with zero crossings, and has the most rigid gate
of the 15 units profiled.** So a *closed and empty* pocket is at least kinetically
stable, which cuts both ways — permissive for signalling (closure does not
require the cavity to be filled) and permissive for the constitutive failure too.

And §24 is the reason it cannot be settled here: **open and closed are both
kinetically trapped**, neither converting once in 1.8 µs of aggregate sampling. So
the accessible observable is the stability of whichever state you start in, never
the free-energy difference between them — and the free-energy difference is
exactly the quantity that separates "signal" from "constitutive".

**The question Jannis is asking is the one our methods are structurally unable to
answer.** Not underpowered — unable, on this class of method.

### 58b. Correction: the §30 factorial is n = 1, not n = 3

I reported the conformation × occupancy factorial as complete. It is not. I
checked that `prod.nc` existed, not how long it was:

| cell | rep0 | rep1 | rep2 |
|---|---|---|---|
| S1 open/apo | 65 ns | **300 ns** | 41 ns |
| S2 closed/+ABA | **300 ns** | 8 ns | 10 ns |
| S9 closed/apo | **300 ns** | 6 ns | 8 ns |
| S10 open/+ABA | **300 ns** | 6 ns | **300 ns** |

**5 of 12 reps reached 300 ns — one per cell.** Three array tasks failed outright
(§ timeline, 27547457_2/4/5) and four more were truncated. The 2 × 2 was designed
at n = 3 specifically so the S9-vs-S2 contrast would have replication; at n = 1 it
cannot deliver that, and its whole purpose was to remove the NaCl/KCl confound
from the §29 comparison.

A second defect, found while extracting: the core superposition mask selected
**656 atoms in the system but 644 in the references**, because both reference PDBs
are missing residues 2, 69–70 and 182–191. Mismatched counts meant the fit
silently did not happen and the gate RMSDs came out at 39–41 Å — impossible for a
five-residue loop. This is the trap script 58's own header documents ("restrict
the CORE to residues present in EVERY unit *and* in the references") and I did not
apply it. The corrected 161-residue common core is written to
`data/md191/core_mask.txt`.

### 58c. What actually follows for testing expanded pockets

1. **An expanded pocket must be paired with a larger ligand.** The expansion is
   not for the ligands that already fit — it is for the **28–50 heavy-atom band**
   (§13), which is why the 229 documented failures were selected in that band in
   the first place. They are the matched test set: same size range, published
   negatives, already characterised.
2. **The constitutive-activity risk has to be screened experimentally.** Y2H
   reports it directly — growth without ligand — and no filter we have
   distinguishes a switch from a constitutively closed protein (§45 passed
   α-estradiol; §56d makes the same point for grafts).
3. **A weak hit remains useful, but for a different purpose than I implied in
   §57.** It collapses the chemical dimension only for ligands that fit the
   current envelope. For an expanded-pocket test it does not transfer, exactly as
   Jannis says.

So the expanded-pocket experiment is necessarily *expanded library × large-ligand
set*, with N bounded by the 229 rather than open-ended — and with the
ligand-independent-activation control run alongside, because that is the failure
mode computation cannot see.

---

## 59. Three follow-ups: enhanced sampling, partial occupancy, and parallel evolution (2026-08-24)

### 59a. Can a ligand bind at the pocket edge and leave a void behind it? — YES, measured

§58 raised the risk that an enlarged cavity with a same-sized ligand cannot
transduce. The Tian clones answer it directly. Net side-chain volume change
against ligand size, over 691 clones with SMILES:

| ligand heavy atoms | n | median volume released |
|---|---|---|
| ≤18 | 222 | −4 Å³ |
| 18–22 | 225 | +21 Å³ |
| 22–26 | 149 | +26 Å³ |
| ≥26 | 95 | +51 Å³ |

Bigger ligands do get more hollowing (r = +0.25), as expected. **But the tail is
the point: 94 clones release ≥ 100 Å³, and 28 of those bind ligands of ≤ 20 heavy
atoms — the same size as ABA.** Several reach the best potency class:

| clone ligand | heavy atoms | ΔV | subs | min_conc |
|---|---|---|---|---|
| **Honokiol** | 20 | **−161 Å³** | 3 | **1 µM** |
| **Magnolol** | 20 | **−161 Å³** | 3 | **1 µM** |
| **Carpropamid** | 20 | −144 Å³ | 4 | **1 µM** |
| Carpropamid | 20 | −185 Å³ | 2 | 10 µM |
| Monobenzone | 15 | −183 Å³ | 3 | 100 µM |

PYR1's cavity is 174 Å³, so −161 Å³ of side-chain volume is close to doubling it —
and Honokiol and Magnolol (isomers, so a consistent pair rather than a fluke)
still give 1 µM sensors with a 20-heavy-atom ligand.

**So the ligand does not have to fill the enlarged cavity.** Empirically, an
ABA-sized ligand can occupy part of a substantially hollowed pocket and still
produce a top-potency sensor. That materially weakens the under-filling worry —
though note the caveat that net side-chain volume is not the same as cavity
volume: the pocket may be *reshaped* rather than left as a void, and water can
fill what remains. What the data establishes is that the **mechanism tolerates it**,
not what the void does.

### 59b. Enhanced sampling — yes, still available, with one specific design

§24's limitation is a property of **unbiased** MD: the barrier is not crossed in
300 ns. Umbrella sampling or metadynamics along a gate coordinate returns the PMF,
which is exactly the open↔closed free-energy difference §58a said we cannot reach.

The risk is the reaction coordinate. Gate closure is a loop rearrangement with
orthogonal slow modes — latch, ligand pose, side-chain repacking — so a single
RMSD coordinate will show hysteresis. That is the same failure class that broke
the TI: a path-dependent estimator whose path is contested.

Two things make it worth doing anyway:

- **Take a difference, not an absolute.** The quantity that matters is not the PMF
  of one pocket but ΔΔG between a **filled** and an **under-filled** cavity.
  Systematic coordinate error largely cancels in that difference — the same logic
  that made TI *selectivity* usable while absolute ddG was not (§43b, §46a).
- **A calibration case already exists.** WT + ABA closes, WT apo is open. Any
  scheme must reproduce that sign before being trusted on a designed pocket —
  the §13e discipline applied to a new method.

Cost is ~20–30 windows × 50–100 ns per system, so 1–3 µs — the same order as
already spent on the factorial, and unlike the factorial it targets the quantity
that decides the question.

### 59c. Parallel directed evolution — feasible, and it tests what computation cannot

Three ligands × five starting structures, selecting for responsiveness.

This is workable, and its main argument is §58a: **the constitutive-versus-signal
question is not answerable by any method in this project, and Y2H answers it
directly** — growth without ligand is the readout for the exact failure mode no
computational filter here can see.

Two design notes:

- **Do not build 15 libraries.** Five starting structures for one ligand will
  share most positions, so one library per ligand spanning all five costs three
  libraries rather than fifteen, and the lineages compete within a single
  selection.
- **Include a positive control.** If all lineages fail, "expansion does not work"
  and "the five starting structures were badly chosen" are indistinguishable. A
  ligand/pocket pair already known to work makes the negative interpretable —
  the same lesson as the ABA null arm in §23j, which was never scored on its own
  correct answer.

Precedent supports the shape of it: Tian's round 2 evolves from round-1 profiles
and rescued five ligands round 1 missed, and §54 found 4 of 11 round-2 sensors had
no round-1 hit at all — so starting points help but are not required.

---

## 60. Enhanced-sampling calibration: WT+ABA vs WT apo (2026-08-24)

§59b argued umbrella sampling is the way to reach the open↔closed free-energy
difference that §24's kinetic trapping puts out of reach — and that it must
reproduce a known answer before being trusted anywhere else. Set up, seeded, and
ready to submit.

### 60a. The reaction coordinate, chosen by measurement

Every gate–latch Cα pair was scored on how far it moves between the open and
closed crystal references. The winner is **P88 Cα – R116 Cα**, and it survives
contact with the unbiased trajectories:

| | P88–R116 (Å) |
|---|---|
| crystal closed / open | 7.64 / 17.52 |
| S2 holo-closed (300 ns) | **6.06 ± 0.75** (5.1–10.3) |
| S9 apo-closed (300 ns) | **6.08 ± 0.43** (5.0–9.4) |
| S1 apo-open (300 ns) | **16.90 ± 1.13** (11.8–20.7) |
| S10 holo-open (2 × 300 ns) | 16.3–16.6 (11.0–20.6) |

**Complete separation, no overlap between basins over 5 × 300 ns.** Note that S9
and S2 sit at the *same* value (6.08 vs 6.06) — the closed state is geometrically
identical with and without ligand, which is precisely why a stability measurement
cannot distinguish them and a free-energy one is required.

A distance is used rather than an RMSD because Amber restrains distances natively
through `&rst`; an RMSD coordinate needs a plugin and adds a failure mode.

### 60b. Setup

- **31 windows per arm**, 5.0–20.0 Å in 0.5 Å steps, k = 10 kcal/mol/Å². Thermal
  width √(kT/k) = 0.24 Å, so neighbours overlap at ~2σ.
- **62 windows seeded from real equilibrated frames**, not from a steered pull —
  the existing factorial already samples nearly the whole range. Only 4 windows
  (apo 10.0–11.5 Å) have seeds more than 0.25 Å off target, max 1.15 Å, which the
  restraint closes during a 2 ns restrained equilibration.
- **Atom indices asserted at runtime** in each arm's own topology (holo and apo use
  different prmtops) — a silent index shift would bias the PMF invisibly.
- 20 ns production per window ⇒ **1.24 µs total**.

### 60c. What is being tested, and how it can fail

> **PASS requires holo to favour CLOSED, apo to favour OPEN, and
> ΔΔG = [G_open−G_closed]_holo − [...]_apo to be POSITIVE.**

The headline is the **difference**, not either absolute PMF: systematic error from
the coordinate largely cancels in it, the same reason TI *selectivity* was usable
when its absolute ddG was not (§43b, §46a).

**The coordinate is the risk and it is stated up front.** Gate closure is a loop
rearrangement with orthogonal slow modes — latch, ligand pose, side-chain
repacking — so a single distance can show hysteresis. Two mitigations are built
in rather than hoped for:

- windows are seeded from **both** basins, so the low and high arms approach the
  barrier from opposite directions and disagreement in the overlap region is
  visible rather than averaged away;
- `113_us_pmf.py` prints **overlap, half-split drift and hysteresis** beside every
  number, the §44a discipline after a tight error bar sat on a moving mean twice.

If this calibration fails, the scheme is wrong and no designed-pocket PMF from it
should be quoted — which is the entire point of running it on a case whose answer
is already known.

### 60d. Extended to a receptor × ligand cross-over (2026-08-24)

Two arms test whether a gate can close. Four test whether the method can tell
*which* ligand closes *which* receptor:

| | WT PYR1 | PYR1^MANDI |
|---|---|---|
| apo | negative (open) | — |
| **+ ABA** | **POSITIVE** (closed) | negative |
| **+ mandipropamid** | negative | **POSITIVE** (closed) |

**The diagonals have to flip.** That is a selectivity test, and it is far stronger
than the two-arm version — a scheme that merely reports "a closed gate is stable"
passes the two-arm test and fails this one.

Built and verified (`114`/`115`), all at 0.15 M KCl on the md191 closed frame:

| arm | res 59 | ligand | P88–R116 | min prot–lig | contacts < 2 Å |
|---|---|---|---|---|---|
| quad_mandi | ARG | 51 atoms | 6.77 Å | 1.88 Å | 2 |
| quad_aba | ARG | 38 atoms | 6.77 Å | 2.15 Å | 0 |
| **wt_mandi** | LYS | 51 atoms | 6.77 Å | **0.27 Å** | **8** |

Two things the build surfaced:

**The restraint atom indices are not shared.** The quad arms put P88 Cα / R116 Cα
at **1408 / 1814**, the WT arms at **1403 / 1819** — K59R adds atoms before
residue 88. `112` originally hardcoded the WT pair; it now derives them from each
topology and asserts the residue identity, because a hardcoded pair would have
restrained the wrong atoms in every quad window without any error.

**wt_mandi starts at 0.27 Å.** That is the physically correct answer — WT cannot
accommodate mandipropamid, which is the entire reason the quadruple was evolved —
but it is also exactly the configuration whose minimisation §46c measured as
inescapable (E = 4.5 × 10⁸, |F|max = 8.4 × 10⁶, flat over 10,000 steps). It must
be relaxed with the backbone restrained and the ligand and side chains free, which
is the fix that worked in `102`'s ladder. If a low-coordinate wt_mandi window
still fails, the ladder remedy applies: seed it from the neighbouring window that
succeeded rather than from the clashed start.

Two build defects, both of the same family — a file that looks right and a parser
that disagrees: the transplanted rotamer atoms were appended after residue 191
instead of in residue order (tleap: `Atom .R<THR 191>.A<OXT 15> does not have a
type`), and a 0-byte prmtop left behind by the failed run passed a `-f` guard that
should have been `-s`.

---

## 61. The factorial, analysed at n=1 — and why it closes the argument for §60 (2026-08-24)

§58b established the 2×2 is n = 1 per cell rather than n = 3. Analysed anyway,
because a *qualitative* observation (does a cell stay in its basin?) does not need
replication, and the result is the strongest available argument for the umbrella
work.

Measured on the reaction coordinate P88 Cα – R116 Cα rather than gate RMSD to a
reference — the coordinate needs no reference and separates the basins better.
Basin assignment is **hysteretic** (closed only below 9 Å, open only above 14 Å)
so a brief excursion into the empty gap is not miscounted as a transition; a naive
midpoint threshold reported 18 spurious "crossings" for S10, whose minimum is
11.01 Å and which never approaches the closed basin at 6 Å.

| cell | rep | mean (Å) | sd | min | max | transitions |
|---|---|---|---|---|---|---|
| open / apo | rep1 | 16.90 | 1.13 | 11.81 | 20.72 | **0** |
| open / +ABA | rep0 | 16.63 | 1.37 | 11.01 | 20.59 | **0** |
| open / +ABA | rep2 | 16.30 | 1.34 | 10.96 | 19.99 | **0** |
| closed / apo | rep0 | **6.08** | 0.43 | 5.01 | 9.35 | **0** |
| closed / +ABA | rep0 | **6.06** | 0.75 | 5.08 | 10.30 | **0** |

### 61a. What it shows

**The closed state sits at the same coordinate with and without ligand — 6.08 vs
6.06 Å, a difference of 0.02 Å.** And ABA does not close an open gate either:
both S10 replicates stay at 16.3–16.6 Å.

⚠ **Corrected in §62: that is true of the MEAN and false of the distribution.**
The two closed basins differ threefold in width and up to 26-fold in
barrier-ward excursions. "Just as stable" was wrong.

**Zero transitions in any cell, in either direction, with or without ligand.**

### 61b. What it does not show, and why that is the point

This is not a null result about the biology. It is a restatement of §24's
trapping in the cleanest possible form: **occupancy changes nothing measurable
about the conformation on this timescale**, because the only accessible observable
is the stability of whichever state you start in.

Which means the §58a question — can an under-filled cavity hold closed, and at
what cost — is not merely unanswered by unbiased MD. It is *unanswerable* by it,
and no amount of additional replication would change that. The 7 truncated
replicates would have added error bars to a quantity that carries no information
about the ligand.

That is the argument for §60. The PMF is not a refinement of this measurement; it
is the only way to measure the thing at all. And §61's numbers set its
acceptance criterion concretely: the two closed states are *geometrically*
identical, so any scheme claiming to separate them must do so on free energy, and
must produce that separation from ensembles whose mean coordinates differ by
0.02 Å.

---

## 62. Why unbiased MD cannot see it — and the part where it can (2026-08-24)

Jannis asked whether the failure in §61 is the force field, or entropy, or
something else — and noted that even a small effect ought to show up.

It is none of those, and the question exposed an overstatement in §61 that needs
correcting first.

### 62a. Correction: MD does show a difference. I compared the wrong statistic.

§61 reported the two closed cells as identical because their means differ by
0.02 Å. Their **distributions** do not:

| | mean | sd | p99 | max | k_eff (kT/var) |
|---|---|---|---|---|---|
| closed / apo | 6.08 | **0.43** | 7.45 | 9.35 | **3.26** kcal/mol/Å² |
| closed / +ABA | 6.06 | **0.75** | 8.65 | 10.30 | **1.06** kcal/mol/Å² |

| frames above | apo | +ABA |
|---|---|---|
| 7.0 Å | 3.78 % | 11.42 % |
| 8.0 Å | 0.19 % | **5.00 %** |
| 9.0 Å | 0.00 % | 0.36 % |

**The holo closed basin is threefold softer and excurses toward the barrier up to
26× more often.** So "an empty closed pocket is just as stable" was wrong — it is
*more rigid*, which is exactly §29's finding that apo-closed has the most rigid
gate of the 15 units profiled. The empty pocket collapses and stiffens; with ABA
present the gate rests on the ligand and has room to breathe.

### 62b. What that difference is worth, and what it is not

Treating the basin as harmonic along this coordinate, the softer well carries more
vibrational entropy:

> ΔF = −kT/2 · ln(k_apo/k_holo) = **−0.33 kcal/mol** in favour of the holo closed state.

Real, in the right direction, and **small**. For comparison, a 10× shift in the
closed:open ratio is 1.37 kcal/mol and a 100× shift is 2.75.

So MD is not blind here. It measures the **curvature** of the basin it is sitting
in, and gets ~0.3 kcal/mol of the answer. What it cannot measure is the basin's
**depth relative to the other basin** — and that is where most of a switch's
ΔΔG lives. A ligand can shift a well's depth by several kcal/mol without changing
its width at all; depth and curvature are independent.

### 62c. Why the ratio is unmeasurable, concretely

The quantity that separates a switch from a constitutive binder is
P(closed)/P(open). Unbiased MD estimates it by counting frames:

| | frames open | estimate |
|---|---|---|
| apo | **0 of 27,000** | P(open) < 3.7 × 10⁻⁵ |
| +ABA | **0 of 27,000** | P(open) < 3.7 × 10⁻⁵ |

Both are **upper bounds, and the same upper bound**. A hundredfold difference
between them is invisible because neither numerator is ever non-zero. The
observable saturates at "never seen".

The timescale arithmetic says why, taking τ = τ₀ exp(ΔG‡/kT) with τ₀ ≈ 1 ns:

| barrier | mean transition time | vs one 300 ns replicate |
|---|---|---|
| 5 kcal/mol | 4.4 µs | 15× |
| 8 | 0.67 ms | 2.2 × 10³ × |
| **10** | **19 ms** | **6.4 × 10⁴ ×** |
| 12 | 0.55 s | 1.8 × 10⁶ × |

§24 saw zero transitions in 1.8 µs of aggregate sampling, which already puts the
barrier above ~5 kcal/mol. At a gate barrier of 10 kcal/mol — unremarkable for a
loop rearrangement with a latch — one transition needs ~19 ms, about **10⁵ times**
what we ran.

### 62d. So: not the force field, and not entropy either

- **Not the force field.** It would govern the *accuracy* of a PMF, but the
  problem here is sampling, not accuracy — and the force field is visibly
  producing a ligand-dependent signal (the threefold width difference).
- **Not entropy as such.** MD captures entropy within a basin it samples; the
  0.33 kcal/mol above *is* an entropy term, measured. What is missing is the
  entropy and enthalpy of the basin never visited, and the population ratio
  between them.
- **It is ergodicity.** The estimator needs both basins; the barrier prevents
  visiting both; so the estimator returns a bound rather than a value.

This is precisely what umbrella sampling fixes: it forces occupancy of the whole
coordinate and removes the bias afterwards, converting an unmeasurable ratio into
a PMF. §61's numbers set the bar — any scheme claiming to separate these two
states must do it on free energy, from ensembles whose means differ by 0.02 Å and
whose widths differ by threefold.

---

## 63. The coumarin secondary library, designed from round-1 alone (2026-08-25)

§54c specified this test and it is now run. `scripts/116_secondary_library.py`,
`scripts/118_force_positions.py`, results in `results/secondary_library/`.

Inputs are sd03 (692 round-1 clones over 194 ligands) and the target SMILES —
what Tian had before building the coumarin library. **sd07 is the answer sheet
and is never an input.** sd08/sd09 stay sealed (§48d).

### 63a. The vocabulary is closed — this really is subset selection

All **23 of 23** of Tian's coumarin substitutions already appear in pooled
round-1's 144. The only round-2 substitutions absent from round-1 are 8
singletons at positions *outside* the 18 — C65F, F71S, R74C, S109R, T124M,
R134L, M178I, D184Y — i.e. PCR carry-through, not designed content. Nothing has
to be invented. §54b's framing is confirmed exactly.

### 63b. The headroom, measured exactly rather than guessed

The smallest library containing at least one known round-2 sensor for **every**
one of the 11 ligands, by branch and bound over the clone choices:

| arithmetic | exact minimum | vs Tian's 138,240 |
|---|---|---|
| wild-type always offered (§51d) | **9,216** | **15×** |
| positions may be **forced** (WT dropped) | **4,608** | **30×** |

and at partial capture, forcing allowed: 10/11 at 2,592 (53×), 9/11 at 1,152
(120×), 8/11 at 384 (360×).

⚠ **A greedy version of this oracle was wrong twice.** It returned 31,104 against
the true 9,216 — understating the prize by 3.4× and nearly killing the exercise —
and because it iterated a *set* of ligand names, `PYTHONHASHSEED` made it return
different answers on different runs. Both are fixed: the search is exact and every
iteration order is sorted.

**So a ≥10× reduction is genuinely available.** §53b said size was the wrong axis
for round 2; that was right about what round 2 *is* and wrong about what it
*could be*.

The shape of the optimum matters more than its size: **11 positions but only 15
substitutions**, so 8 of 11 positions need exactly one residue. Positions are
non-negotiable — every round-2 clone spans 7–8 of them and is lost if one is shut
— and the entire prize is in **residue identity**.

### 63c. Positions are recoverable. Residues are the bottleneck.

Ranking the 18 positions by round-1 substitution weight:

| | true positions in top 11 | in top 12 |
|---|---|---|
| ligand-blind | 10/11 | 10/11 |
| class-weighted (max ECFP4 Tanimoto to the 11 targets, α = 4) | 10/11 | **11/11** |

Rank of the *optimal* residue inside its own position's profile:

| | ranked 1st | top-3 |
|---|---|---|
| ligand-blind | 4/11 | 6/11 |
| **class-weighted** | **7/11** | **9/11** |

This is the first place in the project where knowing the ligand has paid.

### 63d. What the designed libraries actually capture

Ligands captured / clones captured, re-optimising at each budget:

| design | 3,000 | 10,000 | 30,000 | 100,000 | 138,240 |
|---|---|---|---|---|---|
| blind/global | 0/0 | 0/0 | 1/1 | 1/1 | 1/1 |
| blind/profile | 0/0 | 0/0 | 0/0 | 0/0 | 0/0 |
| class/global | 2/3 | 2/3 | 2/3 | 2/3 | 2/3 |
| class/profile | 1/1 | 1/1 | 3/4 | **8/12** | 8/12 |
| class/cover | 1/1 | 1/1 | 2/3 | 6/9 | **9/14** |
| *Tian* | | | | | *11/69* |

**We do not beat the incumbent.** Best is 8–9 of 11 ligands at Tian's own budget.
Leave-one-ligand-out — deleting each target's own round-1 clones before designing
— gives **8/11**, and the three failures are not the four ligands that had no
round-1 hit at all, so this is not simple memorisation.

⚠ Two objectives were tried and discarded, and both failures were informative.
A **global weight-per-log-size greedy** captured 1 of 11: pooled counts are wildly
uneven (F159 alone carries 317 of 1747 round-1 substitutions), so it spent the
whole budget deepening loud positions while quiet ones stayed shut. A
**per-position profile** with forced breadth plateaued at 4 of 11: marginal
frequency picks the wrong residue too often. What works better is a **coverage**
objective over whole round-1 clones — i.e. residue **co-occurrence**, not marginal
frequency. That is a real signal, but its 9/11 sits on a knife edge (the same
objective at exponent 1 or 3 gives 2/11), so it is reported as a lead, not a
method.

### 63e. ⚠ The comparison with Tian is not symmetric, and that bounds the test

sd07's clones were **drawn from** Tian's library, so it contains them by
construction — 11/11 is automatic, not earned. A different library of equal
quality scores badly here because the sensors *it* would have found were never
screened. Hits are positive-unlabeled ([[feedback_hits_are_not_optima]]), so this
metric measures **rediscovery of Tian's choices**, not library quality.

The symmetric metric is recovery of Tian's 23 substitutions — method against
method, both working from the same round-1 data:

| design | 3,000 | 30,000 | 138,240 |
|---|---|---|---|
| blind/global | 9/23 | 13/23 | 13/23 |
| class/profile | 10/23 | 14/23 | **17/23** |
| class/cover | 9/23 | 14/23 | **18/23** |

**17–18 of 23**, from round-1 data alone. That is the defensible claim: we
reproduce most of an expert's library choice, and we cannot yet beat it.

---

## 64. What ties real sensor substitutions together (2026-08-25)

`scripts/117_hit_signature.py`, results in `results/hit_signature/`. Structures
are `data/stage1/wt_{aba,mandi}.pdb`, native numbering, **all 18 positions
asserted to carry their expected wild-type residue** before anything is computed
(§44d).

### 64a. The PYR1^MANDI four — and my first reading was wrong

The obvious story is "the three §32 recovered clash, K59R does not". It is false:

| substitution | clash | ΔVolume | charge | §32 rank of 13,457 |
|---|---|---|---|---|
| F108A | 2.78 Å | **−101 Å³** | 0 → 0 | 5 |
| F159L | 1.42 Å | **−23 Å³** | 0 → 0 | 75 |
| K59R | 0.21 Å | **+5 Å³** | +1 → +1 | **3907** |
| V81I | **0.12 Å** | **+27 Å³** | 0 → 0 | 150 |

**V81I clashes less than K59R and was still recovered.** Clash is not the
discriminator — it is what makes 108 and 159 *easy*, not what makes 59 *hard*.

The axis that separates them cleanly is **volume**. F108A and F159L relieve
overlap by shrinking; V81I is the **grow** half of a shrink/grow pair, recovered
by §32's packing term with no clash to relieve; K59R moves +5 Å³, which is
nothing.

**K59R is triply invisible.** Near-isosteric (+5 Å³), charge-neutral (Lys +1 →
Arg +1), and non-clashing. What it changes is hydrogen-bond *geometry* —
guanidinium is planar and bidentate where ammonium is not, and the crystal makes
two bonds with it (NE–O2 2.64 Å, NH1–O2 3.26 Å, §33a). A volume term, a
formal-charge term and a clash term are three separate blindnesses, which is
exactly why §32, §33 and §36's coupled moves all missed the same residue.

> **The signature is not hydrophobicity, not burial, not contact count. It is:
> does the substitution move volume?**

Percentiles across the 18 positions confirm the negative half — burial 6th–83rd,
contacts 39th–94th, polar-neighbour count 0th–89th. Nothing separates the four.

### 64b. How big the blind spot is

Fraction of substitutions that are near-isosteric (|ΔVolume| < 25 Å³, the band
K59R sits in):

| set | distinct | by occurrence |
|---|---|---|
| sd03 (194 ligands) | 32/113 (28 %) | 25.5 % |
| sd07 (coumarin) | 7/32 (22 %) | 11.9 % |
| Beltran-45 | 10/36 (28 %) | 19.9 % |

A quarter to a third. Not a corner case, not the majority — which is why §32
works at all and why it stops where it does.

### 64c. No structural descriptor reaches significance on n = 18

Spearman against how often real sensors use a position, over all 18:

| descriptor | ρ | perm p |
|---|---|---|
| max overlap with mandipropamid | 0.38 | 0.12 |
| Σ overlap | 0.38 | 0.12 |
| min distance to ligand | −0.36 | 0.15 |
| contacts, burial, polar neighbours, WT volume | 0.03–0.26 | 0.29–0.91 |

Eight descriptors on 18 points. Nothing is called a rule.

⚠ **And the one correlation that did look strong is mostly a null.** WT residue
volume vs % of substitutions that shrink gives ρ = +0.62, p = 0.007 — but a big
residue shrinks because most of the other amino acids are smaller, and the
libraries do not offer all 19 anyway. Scoring against the shrink fraction implied
by each position's **own sd04 menu**: mean enrichment **+0 percentage points**
(sd 22). Most of the signal was the null. What survives is per-position and large
at a few positions — **F108 −54 points** (grows far more than its menu implies),
**V163 −41**, **V83 −37** — but that is empirical, not derived from geometry.

### 64d. The design payoff: which residues can be *forced*

The design-relevant statistic is two numbers, not one: **P(mutated | class)** and
**P(residue | mutated, class)**. High concentration with a low mutation rate means
the residue is right *when used*, not that it should always be used.

A position can be forced — wild-type dropped, so it costs a factor of 1 instead of
2 — only when the second is near 1. Across the 11 coumarin positions, exactly one
qualifies, and it is the one the exact optimum forces:

| set | V163 mutated | identities |
|---|---|---|
| pooled round-1 | 59/692 (9 %) | W 47, M 7, H 3, Q 2 |
| class round-1 (Tanimoto ≥ 0.4) | 11/49 (22 %) | **W 10**, M 1 |
| own-class round-1 | 9/36 (25 %) | **W 9** |
| round-2 sensors | **74/78 (95 %)** | **W 74** |

Pooled round-1 rates V163 at 9 % — bottom of the list, which is why every
ligand-blind method here leaves it out. Restricted to the class the rate rises to
22–25 % **and the identity becomes unanimous**. This one was callable in advance.

Forcing it is worth the 6.5× between 30,000 and 4,608, and **all of the forcing
headroom comes from this single position** — every other position in the optimum
still needs wild-type, because different ligands' sensors disagree there.

⚠ Forcing the next two most concentrated positions **fails**: V163+V83L drops to
9/11 ligands, +F108W to 8/11. Concentration *given a substitution* is not the same
as universality, and 66 % is not 95 %.

Mutual information against a within-position permutation null makes all 11
positions ligand-dependent at p < 0.05 — but with 50–360 substitutions each that
test has power to call trivial effects, so the column that matters is P(top |
mutated), not p.

### 64e. PRE-REGISTRATION: the forcing rule on the sealed sets

Written before sd08/sd09 are read.

> For a target class, over round-1 clones whose ligand has Tanimoto ≥ 0.4 to the
> class, compute P(mutated) and P(top residue | mutated). **Predict that any
> position with P(top | mutated) ≥ 0.85 on n ≥ 8 class clones carries that residue
> in ≥ 80 % of that class's round-2 sensors, and may be forced.**

On the coumarin dev set the rule fires **exactly once** (V163W, 91 % on n = 11)
and is right (95 %). A rule that fires once on one class is not a method; PFAS and
TNT decide it. It fails if it fires on a position whose round-2 frequency is
< 80 %, or fires nowhere on either sealed class while a forceable position exists
in their round-2 data.

---

## 65. The forcing rule, run prospectively on PFAS (2026-08-25)

Jannis authorised unsealing sd09. `scripts/119_pfas_prospective.py`, results in
`results/secondary_library/pfas_prospective.*`.

**TNT is not run.** Its round-2 library was not built from hits on chemically
related ligands — none had hits — but from ligands sharing specific structural
features. The rule keys on 2D whole-molecule similarity, so its input does not
exist for TNT. sd08 stays sealed.

### 65a. ⚠ The first run used the wrong round-1 input and its verdict is void

I took round 1 to be sd03 — the 194 hits from the 3,366-compound Selleck/LATCA
deck — which contains **no PFAS at all** (0 compounds with "perfluoro" in the
name; the most fluorinated is perflubron at 17 F, not a hit). The class set came
out empty, the rule fired nowhere by default, and it "passed" by having nothing
to say.

Tian in fact ran a separate PFAS round-1 screen: an improved DSM-Hao library
against a panel of 103 PFAS. **That data is inside sd09 itself**, separated by
`mut_lib`:

| `mut_lib` | clones | ligands | role |
|---|---|---|---|
| `DSM-Hao` | 89 | 18 | **round 1** — the design input |
| `PFOS_GenWT` | 154 | 25 | **round 2** — the answer sheet |

Seven of the 25 round-2 targets had no round-1 hit — PFOA, PFHxS, PFNA,
perfluorodecanoic acid, and three fluorotelomer alcohols. Same shape as
coumarin's 4 of 11 (§54a), and it again supports "an initial hit is not
necessary".

### 65b. The pre-registered test, on the right data

Class set: 89 clones at Tanimoto ≥ 0.40 — well past the n ≥ 8 trigger, so the
rule had every opportunity to fire.

**It fired nowhere.** No position's round-1 residue identity reaches
P(top | mutated) ≥ 0.85. The most-mutated position, E94 (54 of 89 clones), splits
G 0.59 / others; Y120 splits A 0.23; A160 splits I 0.28.

**And nothing was forceable.** Best round-2 "carries one specific residue" is
K59→M at **0.44**, against the 0.80 bar. Nothing else exceeds 0.39.

> **Verdict: true negative.** The rule declined and there was genuinely nothing to
> fire on. Unlike the void first run, this one had the data to catch a false
> positive and did not produce one — but it still cannot demonstrate the rule
> *finds* things, because PFAS has nothing to find.

### 65c. ⚠ A correction to §64d's statistic — P(top | mutated) is partly circular

Checking why PFAS lacks consensus exposed a confound in my own measure. **Where a
library offers exactly one substitution at a position, round-2 P(top | mutated) is
1.00 by construction, not by selection.**

Every coumarin position I reported at 1.00 — V83, F108, V163, N167 — offered
exactly one residue. So the round-2 half of the coumarin claim was inflated.

What survives, and is *not* circular, is **P(mutated)**: whether the winners kept
wild-type when they could have. Re-reading §64d through that lens:

| coumarin position | options offered | P(mutated) | forceable? |
|---|---|---|---|
| **V163** | W (1) | **0.95** | **yes — collapses to one residue, 2×** |
| A160 | GIMV (4) | 0.87 | wild-type droppable, but 4 residues remain → 1.25× |
| V83 | L (1) | 0.74 | no |

The exact branch-and-bound in §63b was never affected — it forced V163 and only
V163, working from clone sets rather than from this statistic. It is the *rule's*
predictor that needed the correction, and the round-1 side is clean (DSM-Hao
offers 6–17 residues per position, so round-1 unanimity is real selection).

### 65d. What did transfer: P(mutated), at ρ = +0.74

Post-hoc, not part of the test. Round-1 mutation rate predicts round-2 mutation
rate across 17 positions at **Spearman +0.74**. Three positions drop wild-type in
round 2 at ≥ 80 %: **E94 (0.87), Y120 (0.86), K59 (0.83)**.

Round-1 called E94 (rank 1 of 18) and Y120 (rank 4). It missed **K59 badly** —
rank 7, mutated in 6 % of round-1 clones and 83 % of round-2 clones. So the
ordering transfers, the top of it does not.

**And the payoff is small here, for a structural reason worth keeping:** forcing
pays in proportion to how *shallow* the menu is. Dropping wild-type at all three
takes 49,545,216 → 27,525,120, a **1.8×**. A position offering one substitution
halves the library; one offering eight saves 11 %. Coumarin's V163 offered exactly
one residue; PFAS's K59/E94/Y120 offer 3, 5 and 8.

### 65e. The two classes differ in kind, and I could not explain why

- **coumarin** — one position converges on one residue → force it, 2×
- **PFAS** — three positions drop wild-type but spread across residues → 1.8×, and
  nothing to collapse to

⚠ **My proposed explanation is refuted.** I guessed PFAS sensors spread out
because the PFAS panel is chemically broader than the coumarin panel. Mean
pairwise ECFP4 Tanimoto within each class: **coumarin 0.352, PFAS 0.339** —
indistinguishable. Class breadth does not explain it.

What remains unexplained is therefore worth stating plainly: **whether a class
admits a forceable position appears to be a property of the receptor–class
interaction, not of the class's chemical spread, and nothing here predicts it in
advance.** The rule can say "no" honestly. It has said "yes" exactly once, on one
position, retrospectively. That is where it stands.

### 65f. The PFAS headroom is enormous — and it is positions, not forcing

Same exact branch-and-bound as §63b, over the 154 round-2 clones (16.9 M nodes,
exhaustive):

| arithmetic | exact minimum | vs Tian's 49,545,216 |
|---|---|---|
| wild-type always offered | **14,400** | **3,441×** |
| forcing allowed | **9,600** | **5,161×** |

Against coumarin's 15× / 30×, that is a different regime entirely. The optimum
menu is

```
K59:MN(forced)  E94:ADEGS  Y120:LMTVY  S122:LS  E141:EGQ  F159:FHIL  A160:ALMV  V164:VW
```

— **8 variable positions, not 13.** Five of Tian's positions (V83, L87, A89,
V163, N167) contribute nothing: every one of the 154 sensors keeps wild-type
there. Almost all of the 3,441× is **positions that did not need to be open**, and
only 1.5× of it is the K59 forcing.

So the two libraries were sized very differently against their own results.
Tian's coumarin library was near-optimally tight (15× off); the PFAS library was
**over-hedged by three orders of magnitude**. That tracks how much class-relevant
round-1 data each had: coumarin had 36 own-class round-1 clones with a clean
position signal, PFAS had 89 clones whose best position (E94) reached only 0.61
and whose identities never converged. Less confidence, wider library.

⚠ **Two things this is not.** It is not a claim that a 14,400-member library would
have *found* those 154 sensors — they were found by screening 49.5 M, and hits are
positive-unlabeled ([[feedback_hits_are_not_optima]]); a tighter library contains
them but was not where anyone was looking. And the optimum is fitted to the clones
actually recovered, so it is a bound on the prize, not a design.

But it does relocate where the prize is. For coumarin the whole 15× sat in residue
identity at fixed positions. For PFAS the 3,441× sits almost entirely in **which
positions to open at all** — and that is the half §63c already showed is
recoverable, with class-weighted round-1 putting 11/11 true coumarin positions in
its top 12, and round-1 P(mutated) predicting round-2 P(mutated) at ρ = +0.74
here.

### 65g. Position selection works; residue selection is still the bottleneck

Ranking the 18 positions by round-1 P(mutated) on the 89 PFAS round-1 clones,
against the 8 positions the exact optimum actually needs:

| | positions recovered |
|---|---|
| top 6 | 6/8 |
| **top 8** | **7/8** |
| top 10 | **8/8** |

The ranking is almost monotone in usefulness — ranks 1–7 are all needed, and the
only break is S122 at rank 10. **This is the strongest position-selection result
in the project**, and it is on the out-of-domain class.

But turning it into a library still fails on residue identity:

| positions × residues | size | ligands | clones |
|---|---|---|---|
| top 8 × top 2 | 6,561 | 4/25 | 4/154 |
| top 8 × top 3 | 65,536 | 7/25 | 12/154 |
| top 8 × top 4 | 312,500 | 14/25 | 21/154 |
| **top 8 × top 5** | **933,120** | **15/25** | 26/154 |
| *Tian* | *49,545,216* | *25/25* | *154/154* |

So a **53× smaller** library reaches 15 of 25 targets. That is a real trade rather
than a win, and the positive-unlabeled caveat cuts both ways here: those 154
sensors were found by screening 49.5 M, and a 933 K library would plausibly find
*different* sensors for some of the 10 it misses — but nothing here measures that.

The pattern is identical to coumarin (§63c): **positions recoverable, residues
not.** Two classes, two chemistries, same boundary.

---

## 66. Three questions: hybrid scoring, the START-fold ligand census, and pocket vs interface (2026-08-25)

### 66a. Can a scorer rank COMBINATIONS once the wet lab has fixed the menu? — no

`scripts/122_combination_viability.py`, `123_viability_aggregate.py`, job 27752872
(cutlerlab, CPU, no GPU allowance touched).

Every previous attempt asked a structure-based tool to pick the right *residue*,
and every one failed. §63–§65 showed round-1 data already narrows the menu. So the
tools no longer have to generate candidates — only to **rank combinations inside a
menu the wet lab chose**. And that can be done **pose-free**, which sidesteps the
circularity that killed everything else (§49b, §46b's 8.33 Å docking error): a
working sensor must at minimum fold and pack, and asking whether a combination of
5–8 mutations is structurally viable needs no ligand.

670 variants, protein-only ref2015, shell repacked, each against a **paired
wild-type control repacked in the identical shell**:

| set | n | median ΔΔG | p25 | p75 | catastrophic (> 200 REU) |
|---|---|---|---|---|---|
| REAL (sd07 sensors) | 70 | 66.7 | 31.8 | 125.7 | **22.9 %** |
| LIBRARY (random from Tian's 138,240) | 300 | 53.1 | 20.9 | 762.4 | 32.3 % |
| WILD (random from DSM-Hao at the same 11 positions) | 300 | 89.8 | 41.3 | 292.8 | 30.3 % |

| comparison | AUC | p |
|---|---|---|
| REAL vs WILD | **0.443** | 0.14 |
| REAL vs LIBRARY | **0.542** | 0.28 |

**No separation, from either null.** Excluding the catastrophic tail does not
rescue it: REAL vs LIBRARY AUC falls to **0.333** and REAL vs WILD to 0.550.

And the direction is the interesting part. Real sensors are *slightly more
strained* than random library members, which is mechanistically sensible: a sensor
works by **remodelling** the cavity, and remodelling costs packing energy. The
property that makes a combination good is the one ref2015 penalises. Protein-only
stability is not merely uninformative here — it is mildly anti-correlated.

The only real signal is the clash tail: 22.9 % of real sensors are catastrophic
against ~31 % of both nulls. A hard clash filter would cut ~30 % of the library
while losing ~23 % of the sensors. That is barely better than deleting at random.

⚠ **What this does and does not close.** Ruled out: cheap protein-only repack
scoring as a library filter. Not tested: a FastRelax-based version (~100× the
cost, and the prior is now poor), and any ligand-aware score — still blocked by
the pose problem, not by cost.

⚠ **A confound found and fixed mid-run.** The first version scored the repacked
mutant against the *raw* wild-type and returned ΔΔG of −130 to −165 REU for every
variant — that was the input structure's own strain being repacked away, plus a
bias toward variants with more mutations, which open larger shells. The paired
control removes both. The early "signal" that motivated the full run was two
extreme WILD variants and did not survive n = 670.

### 66b. START/SRPBCC-fold proteins with a known pocket ligand: **11**

`scripts/120_start_ligand_census.py`. Of the 261 homologs the foldseek search
returned, **195 are AlphaFold models** (no ligands by construction) and **66 are
experimental**. Of those 66, **11 have a buried pocket ligand**:

| | ligand | heavy atoms | | | ligand | heavy atoms |
|---|---|---|---|---|---|---|
| 2flhA00 | ZEA | 16 | | 4jhiA00 | EMU | 17 |
| 3h3qA00 | H13 | 27 | | 4n3eY00 | 2AN | 21 |
| 3oquB00 | **A8S (ABA)** | 19 | | 4qdcA02 | ASD | 21 |
| 3putA00 | HEZ | 8 | | 5nonC00 | 93H | 42 |
| 3tfzE00 | CXS | 14 | | 6awvC00 | 28E | 21 |
| 4dsbB00 | **A8S (ABA)** | 19 | | | | |

Two are ABA in PYL-family proteins, and HEZ (hexanediol) and CXS (a CHAPS-like
sulfonate) are buffer additives that happen to be buried. **So the fold offers
about seven non-ABA biological pocket ligands.**

⚠ **The obvious directory is the wrong one.** `data/survey_cache/` holds 478 CIFs
and returns 97 % ligand-bound with bacteriochlorophyll, chlorophyll, carotenoid
and lauryl maltoside at the top — it was built by scripts 44–46 for the
floppy-ligand survey and is full of light-harvesting complexes, not START
proteins. And a purely compositional count over any set returns 100 %, because
NAG/BMA/MAN glycans, TYS/CSO/TPO modified residues and ACE caps are all HETATM.
The census is geometric: ≥ 8 heavy atoms, ≥ 70 % of atoms contacting protein
within 4.5 Å, ≥ 20 contacts.

**Seven ligands is not a training set.** Against Tian's 208 ligands with sensor
data, the structural record for this fold adds essentially nothing on the chemistry
axis. Where it could still matter is geometry — the 42-heavy-atom 93H shows the
fold *can* hold something twice ABA's size.

### 66c. The pocket and the HAB1 interface share one helix, and F159 is both

`scripts/121_pocket_vs_interface.py`, on `data/complex_AB_ABA.pdb`, secondary
structure from PyRosetta DSSP.

| set | n | residues |
|---|---|---|
| HAB1-contacting (< 4.5 Å) | 19 | 60, 61, 63, 84–89, 116, 117, 148, 151, 155, 156, 158, **159**, 162, 166 |
| ABA-contacting (< 4.5 Å) | 19 | 59, 61, 83, 87, 88, 89, 91, 92, 94, 108, 110, 115, 117, 120, 141, **159**, 163, 164, 167 |
| **both** | 6 | 61, 87, 88, 89, 117, **159** |

**Four of the 18 library positions touch HAB1 directly: L87, A89, L117, F159.**
The first three are the gate (87–89) and latch (117) — the ratchet residues, which
is expected. **F159 is the one that should give pause**: it is the most-mutated
position in round 1 (317 of 1747 substitutions) and it sits 3.2 Å from *both* the
ligand and HAB1.

The structural picture is one element doing two jobs. The C-terminal helix
**α3 (153–180)** carries **five library positions (159, 160, 163, 164, 167)** and
**six HAB1-contacting residues (155, 156, 158, 159, 162, 166)**, interleaved —
and per-residue distances show they are on **opposite faces**:

```
D155  HAB1 3.8   T156  HAB1 4.1   M158  HAB1 2.6
F159  ABA 3.2 / HAB1 3.2   <- both
A160  ABA 4.6   V163  ABA 3.5   V164  ABA 4.4   N167  ABA 3.8
T162  HAB1 3.7   L166  HAB1 3.7
```

This is why HAB1 contains no pocket residues (the ratchet framing is intact) and
yet a third of the library sits on the helix that presents the interface. A
mutation at 160, 163, 164 or 167 never touches HAB1, but it repacks the core of
the helix that does — a route to a dead or constitutive sensor that nothing in
this project currently screens for, and a plausible contributor to the
constitutive clones seen in the Y2H work.

⚠ **cpptraj's DSSP failed silently** on this input — no backbone amide hydrogens,
so every residue came back Turn or all-zero. My first parser took argmax over an
all-zero row, which returns column 0, so 0.0 became "Extended" and PYR1 was
reported as **20 strands and zero helices**. A helix-grip fold with no helices
should have stopped me. Replaced with PyRosetta DSSP (53 H, 64 E, 62 L).

---

## 67. FastRelax reverses §66a — and the reversal is the result (2026-08-25)

`scripts/124_viability_relax.py`, `124b` (job 27803659, cutlerlab CPU),
aggregated with `123 --dir results/viability_relax`.

Identical experiment to §66a — same 670 variants, same three sets, same 8 Å shell,
same paired wild-type control — with **Cartesian FastRelax** (backbone and side
chains free inside the shell) in place of fixed-backbone repacking. 3 replicates
per structure, minimum taken.

### 67a. Opposite answer

| | REAL (70) | LIBRARY (300) | WILD (300) |
|---|---|---|---|
| **repack** median ΔΔG (§66a) | 66.7 | **53.1** | 89.8 |
| **relax** median ΔΔG | **−1.1** | 3.3 | 8.6 |

| comparison | repack AUC | **relax AUC** | relax p |
|---|---|---|---|
| REAL vs LIBRARY | 0.542 (REAL worse) | **0.284 (REAL better)** | **1.7 × 10⁻⁸** |
| REAL vs WILD | 0.443 | **0.117 (REAL better)** | **2.1 × 10⁻²³** |

Fixed-backbone scoring did not merely fail to see the signal — it **inverted** it.
§66a's conclusion that real sensors are "more strained" was an artefact of denying
them the backbone motion they need. With relaxation the ordering is clean and
monotone: REAL < LIBRARY < WILD.

Controls: n_sub 6.6/6.7/6.6 and shell 52.1/51.8/51.6 across sets, so the obvious
confound is matched. Replicate spread median 1.34 REU (mutant), 0.30 (wild-type),
against between-set gaps of 4.3–9.6 REU — **noise below signal**. Non-shell drift
**0.000 Å**.

⚠ **Two protocol bugs found on the way, both silent.** Torsion-space FastRelax with
the backbone free only inside the shell moved "frozen" residues by **0.98 Å** — a
phi/psi change at *i* rotates the whole chain after *i*, and the shells differ
between sets, so it would have been a systematic between-set artefact. Cartesian
fixes it. And the aggregator's AUC label read `> 0.5 means REAL scores BETTER`
when the statistic is P(REAL ranks *above* the null) on a lower-is-better
quantity — the prose in §66a was taken off the medians and is unaffected, but the
label was backwards for two runs.

### 67b. The score is 89 % additive, and orthogonal to round-1 frequency

An additive per-substitution model fitted on LIBRARY+WILD only (600 variants, 175
substitutions) explains **R² = 0.892**. Removing it moves REAL vs LIBRARY from
AUC 0.284 to **0.400** — most of the separation is per-substitution, some residual
epistasis survives.

The per-substitution coefficients are the useful object, because they are
**independent of what round 1 already told us**: Spearman(Rosetta β, round-1
pooled frequency) = **+0.03**.

Predicting which of Tian's 23 coumarin substitutions actually appear in round-2
sensors:

| predictor | Spearman |
|---|---|
| **Rosetta relax β** | **+0.46** (permutation p = 0.029) |
| round-1 class-weighted frequency | +0.13 |
| round-1 pooled frequency | **−0.27** |
| any weighted combination of the two | ≤ +0.46 |

Round-1 pooled frequency points the **wrong way**, and the reason is visible in the
table: **F159I is the single most common round-1 substitution (66 occurrences) and
appears in exactly one round-2 sensor**; F159V (64) likewise once; **V81Y (61) is
used 12 times and is the substitution Rosetta penalises hardest (+13.7 REU)**.
Meanwhile V81I (round-1 count 11) is used 39 times. Round 1 over-weights F159 and
V81Y badly; relax does not.

### 67c. ⚠ But it cannot CHOOSE the menu — only rank within one

Substituting the Rosetta score into §63d's menu-selection step makes the library
**worse, not better**:

| scorer | best library ≤ 200 K | ligands captured |
|---|---|---|
| class-weighted frequency (§63d) | 103,680 | **8/11** |
| class × Rosetta | 30,720–73,728 | **0/11** |
| Rosetta alone | — | nothing under budget |

Rosetta selects for **stability**, and the residues that create a new binding site
are not the stable ones. Used as a menu generator it picks residues that pack well
and do nothing. So the honest scope is narrow and specific:

> **Rosetta relax is a filter over combinations inside a menu someone else chose.
> It is not a menu generator.**

### 67d. What the filter is actually worth

| keep this fraction of known sensors | ΔΔG threshold | library retained | reduction |
|---|---|---|---|
| 95 % | ≤ 6.71 REU | 63.0 % | 1.59× |
| **90 %** | **≤ 5.32 REU** | **57.3 %** | **1.74×** |
| 80 % | ≤ 3.52 REU | 50.7 % | 1.97× |

**1.7× at 90 % sensor retention** — real, statistically solid, and modest. On WILD
the same threshold retains only 31 %, so it is much better at rejecting residues
the wet lab never selected than at pruning the ones it did.

⚠ LIBRARY is unlabelled, not a set of known failures, so "removes 43 % of the
library" is not "removes 43 % correctly" ([[feedback_hits_are_not_optima]]). The
retention figure for REAL is the reliable half.

### 67e. So the hybrid works, in one specific place

The division of labour the data supports:

| step | what does it | evidence |
|---|---|---|
| which positions to open | round-1 mutation rate | 7/8 in top 8 on PFAS (§65g), 11/11 in top 12 on coumarin (§63c) |
| which residues to allow | round-1 class-weighted frequency | §63c; **not** Rosetta (§67c) |
| which combinations to keep | **Cartesian FastRelax ΔΔG** | AUC 0.284, p 1.7 × 10⁻⁸; 1.7× at 90 % retention |

Cost: ~21 CPU-hours for 670 variants, i.e. ~2 CPU-minutes per variant. Filtering a
10⁵-member library outright is not affordable at that rate (≈ 3,800 CPU-days); the
additive model (R² = 0.89) is the practical route — score a few hundred sampled
combinations, fit per-substitution coefficients, and apply those to the full
library for free.

---

## 68. Why the relax filter works, and why it will not transfer to expansion (2026-08-25)

§67 established that Cartesian FastRelax separates real sensors from random
library members (AUC 0.284) but cannot choose a menu (0/11). Both halves have the
same explanation, and following it out produces a warning about our actual goal.

### 68a. What it measures: a necessary condition, not a sufficient one

Protein-only ΔΔG asks *does this combination of side chains pack into this
backbone*. A working sensor needs two things:

- **(A) it must fold and pack** — necessary, not sufficient
- **(B) the cavity must fit the new ligand** — needs the ligand

The score sees (A) and is blind to (B). That predicts exactly the observed pattern:
good at **rejecting** (things failing (A) are certainly not sensors), useless at
**choosing** (maximising (A) gives a well-packed protein with no site).

### 68b. Why (A) binds at all: it is a DEPTH effect

Round-1 clones carry 2–3 substitutions and round-2 clones 5–8. Stability is not
limiting at 2–3; at 7–8 it is, because strain accumulates. Prediction: the
separation should widen with depth. It does, monotonically:

| substitutions | n REAL | n LIBRARY | REAL median | LIBRARY median | AUC |
|---|---|---|---|---|---|
| 3–5 | 13 | 103 | −0.52 | +1.81 | 0.359 |
| 6 | 13 | 34 | −0.71 | +5.47 | 0.251 |
| 7 | 27 | 37 | +0.15 | +5.16 | 0.259 |
| **8–10** | 17 | 126 | **−2.73** | +4.20 | **0.219** |

And the direction inside each set says the same thing:
Spearman(n_sub, ΔΔG) = **−0.12 in REAL, +0.14 in LIBRARY, +0.29 in WILD.**
Random stacking accumulates strain; real sensors do not — the deepest real sensors
are the *most* stable group. Real sensors are mutually compensating combinations,
and that is precisely what relax can see and a per-substitution rule cannot.

**This also explains why round-1 frequency fails at it (§67b).** Round-1 frequency
is measured in the shallow regime where stability never binds, so it cannot encode
a constraint that only appears at depth.

### 68c. Why it cannot choose a menu: ref2015 fills cavities

A cavity is a packing defect and ref2015 penalises it. So the score systematically
prefers **growing** the lining:

| | mean β | n |
|---|---|---|
| GROW (ΔVol > +20 Å³) | **+0.43 REU** | 76 |
| SHRINK (ΔVol < −20 Å³) | **+1.75 REU** | 76 |

Spearman(β, ΔVolume) = **−0.25**. Left to choose, it picks K59D, A160L, N167Y,
S122Q, V163D — **none of which appears in a single round-2 sensor.** It designs a
well-packed protein with the pocket filled in.

### 68d. ⚠ The coumarin validation was in the FAVOURABLE regime

Coumarins are **11–15 heavy atoms; ABA is 19**. A smaller ligand needs the lining
to grow *inward*, which is the same direction ref2015 wants. Real coumarin sensors
average **+23.1 Å³ per substitution**, about +162 Å³ over a 7-substitution sensor.
So §67's success was partly luck of the target class.

Across 125 sd03 ligands with ≥ 4 substitutions,
Spearman(ligand heavy atoms, mean ΔVolume used) = **−0.30**:

| ligand size | ligands | mean ΔVolume per substitution |
|---|---|---|
| ≤ 15 (smaller than ABA) | 33 | **+2.4 Å³** |
| 16–20 (ABA-sized) | 36 | −2.3 |
| 21–28 | 46 | −9.0 |
| **≥ 29 (the §58c expansion band)** | 10 | **−15.1 Å³** |

> **For ligands larger than ABA, real sensors SHRINK the lining — and that is
> exactly the move ref2015 penalises most. The filter that works for coumarins is
> expected to work AGAINST us on pocket expansion.**

This is testable rather than speculative: rerun §67 with a large-ligand class
(sd03's ≥ 29-heavy-atom ligands, or Beltran's cannabinoids at 21–25) and the AUC
should move toward or past 0.5.

### 68e. The alternative the reasoning points to

The defect is not the sampling (Cartesian relax is fine) and not the backbone
(§67a fixed that). It is the **objective**: ΔΔG minimises cavity, and a sensor
needs a cavity of a *particular size*. So replace the target rather than the
method:

> score a combination by **|cavity volume − target ligand volume|**, not by ΔΔG.

That stays pose-free — it needs the ligand's *volume*, which is a 2D property of
the SMILES, not its pose, so it does not re-import the circularity of §49b. The
machinery exists (`lib_cavity.py`, §04/§16 measured PYR1's cavity at 174 Å³ and
donor pockets up to 570 Å³). ΔΔG is then kept as a **constraint**, not an
objective: reject the strained, then select on cavity match among what survives.

Retrospective test, runnable on the variants already scored: compute cavity volume
for all 670 and ask whether **cavity-match** separates REAL from LIBRARY for
coumarins *and* keeps separating them for a large-ligand class, where ΔΔG is
predicted to fail.

---

## 69. The cavity-match test: my mechanism was wrong, and the filter is better than I claimed (2026-08-25)

`scripts/125_cavity_match.py`, `126_cavity_aggregate.py`, jobs 27812914 (coumarin)
and 27812915 (PFAS), cutlerlab CPU. Predictions were fixed in 125's header before
the run; two of the three failed.

Two arms on opposite sides of the wild-type cavity (**165.9 Å³**, measured here,
consistent with §04):

| arm | target ligand volume | pocket must |
|---|---|---|
| coumarin | 162.7 Å³ | barely change |
| **PFAS** | **224.6 Å³** | **open by ~60 Å³** — the expansion regime |

### 69a. Verdict

| prediction | result | |
|---|---|---|
| 1. ΔΔG works on coumarin (AUC < 0.40) | **0.280** | HOLDS |
| 2. ΔΔG **fails** on PFAS (AUC ≥ 0.45) | **0.334** | **FAILS** |
| 3a. cavity match works on coumarin | 0.540 | **FAILS** |
| 3b. cavity match works on PFAS | 0.461 | **FAILS** |

**Cavity match is dead** — it cannot reproduce even the arm where ΔΔG already
works. And **§68d's central prediction is refuted**: ΔΔG did not fail in the
expansion regime.

### 69b. ⚠ A confound that was working AGAINST the result

PFAS REAL sensors carry 5.4 substitutions on average against LIBRARY's 3.9 — the
random draws were uniform over the observed range while the real distribution is
skewed high. Since §68b established that ΔΔG separation *widens* with depth, the
raw PFAS number was **understating** the effect. Stratifying by substitution count
and re-weighting:

| | raw AUC | **n_sub-matched AUC** |
|---|---|---|
| coumarin | 0.280 | **0.251** |
| **PFAS** | 0.334 | **0.206** |

**ΔΔG works *better* in the expansion regime than in the coumarin one.** And the
depth effect reappears inside PFAS exactly as §68b predicts: AUC 0.350 at 4
substitutions → 0.288 at 5 → **0.145 at 6**.

### 69c. Why §68d was wrong

The argument was: ref2015 prefers GROW substitutions; ligands larger than ABA need
SHRINK substitutions; therefore ΔΔG should oppose large-ligand sensors. The two
measured premises are still true — mean β **+0.43 grow vs +1.75 shrink**, and
Spearman(ligand size, ΔVolume used) = **−0.30**. The inference was wrong.

The direction check shows why:

| arm | REAL Δcavity | LIBRARY Δcavity | WILD Δcavity |
|---|---|---|---|
| coumarin (needs shrink) | **−35.4 Å³** | −3.0 | −4.4 |
| PFAS (needs grow) | **+39.2 Å³** | +31.6 | **+1.2** |

Real sensors move the cavity strongly in the required direction in **both**
arms — and so do random members of the same library (+31.6 for PFAS). **The menu
already encodes the direction.** Tian's PFAS menu is built from pocket-opening
substitutions, so every member of it opens the pocket, and ref2015's grow-bias
never gets the chance to express itself. Within a directionally-correct menu, what
ΔΔG ranks is **combinatorial compatibility**, not cavity size.

So §68b (the depth/compensation mechanism) survives and is the real explanation;
§68c's volume bias is correctly measured but only matters when the score is asked
to *choose* a menu (§67c), which is exactly where it already failed. The two
findings are consistent: **the menu supplies the direction, ΔΔG supplies the
compatibility.**

And that is also why cavity match adds nothing — it re-measures the direction the
menu has already fixed. Its one significant result is REAL vs **WILD** on PFAS
(AUC 0.301, p 5 × 10⁻¹²), i.e. it can tell a pocket-opening library from an
arbitrary one, which the wet lab does not need help with.

### 69d. The filter, restated on matched strata

| class | keep 95 % of sensors | keep 90 % | keep 80 % |
|---|---|---|---|
| coumarin | 1.93× | **2.10×** | 2.32× |
| PFAS | 1.54× | **2.05×** | 2.42× |

**~2× at 90 % sensor retention, in both regimes** — up from §67d's 1.74×, which was
computed without stratifying. Modest, consistent, and now demonstrated on the
class we actually care about.

⚠ LIBRARY is unlabelled, so these are lower bounds on the separation and the
reduction figure is what the filter *removes*, not what it removes correctly.

---

## 70. ⚠ The umbrella sampling has been half-dead since launch (2026-08-25)

Found while resubmitting the array. **38 of the 62 windows have never produced a
single frame**, and I had been reporting "24/62 complete, resubmit to finish" as
though the rest were merely unstarted.

### 70a. The failure

```
| ERROR:   natom mismatch in inpcrd/restrt and prmtop files!
```

The four md191 systems were solvated **independently** and have different atom
counts:

| system | atoms |
|---|---|
| S2_holo_closed | 46,570 |
| S10_holo_open | **51,588** |
| S9_apo_closed | 46,483 |
| S1_apo_open | **51,585** |

`111_us_seed.py` seeded the open-basin windows from the OPEN systems, while
`112_us_run.sh` runs every window in an arm under the **closed** system's topology.
The split is exact and total:

| windows | seeded from | result |
|---|---|---|
| 5.0–10.5 Å (24) | S2_holo_closed / S9_apo_closed | **20.0 ns each, complete** |
| **11.0–20.0 Å (38)** | S10_holo_open / S1_apo_open | **0.0 ns, dead on arrival** |

### 70b. Why it matters more than "38 jobs to rerun"

The completed windows span **5.0–10.5 Å — the closed basin only.** §60's headline
quantity is *G(open) − G(closed)*, and the open basin sits at 16.3–16.9 Å. **No PMF
was ever computable from what had run**, in either arm. The calibration was not
partially finished; it was structurally incapable of producing its result.

### 70c. My own guard was one level too shallow

`111`'s docstring says: *"the two arms use DIFFERENT topologies (holo carries A8S),
so a seed frame may only be taken from a trajectory built on the same topology."*
I checked holo-vs-apo and never considered closed-vs-open **within** an arm —
because I was thinking of topology as "what molecules are present" when the thing
that actually has to match is the atom count, which independent solvation changes.
This is the same class as §60d's restraint-atom indices, where K59R shifting atom
numbering would have silently restrained the wrong pair: **an identity that looks
obviously shared and is not.**

It also should have been caught immediately. The windows failed in 23 seconds
each; a 20 ns window takes 2 hours. I had a `have_ns` progress function and read
its output as "not started yet" rather than asking why nothing had started.

### 70d. The fix — adiabatic pulling, `scripts/127_us_reseed_open.sh`

The closed-system box can hold the open state, so no rebuild is needed:
clearance is **19.7–22.3 Å** against an open-state protein extent of 52.0–55.2 Å,
which is no larger than the closed states' 54.9/57.4 Å — the gate opens, the fold
does not expand.

So: start from window **10.5 Å**, which is complete, equilibrated and already in
the correct topology, and step the restraint centre outward in 0.5 Å increments of
200 ps, saving each step's final restart as the next window's seed. 19 steps per
arm = **3.8 ns**, against ~10 ns for a Jarzynski pull, and it produces one seed per
window directly.

Slow stepping is the point: `112`'s equilibration is 2 ns, which can settle a
structure already near its target but cannot drag one 10 Å — that is why this has
to be a separate stage rather than something the restraint absorbs.

Cost to finish: 3.8 ns × 2 arms of pulling, then 38 windows × 22 ns ≈ **840 ns**.

---

## 71. ⚠ The headroom numbers were measured against too easy a target (2026-08-26)

Jannis: the round-2 hits came out of Tian's library, so if that library was built
suboptimally we are asking a method to rediscover a mediocre answer — better to
target the **best** sensor per ligand. `scripts/128_best_hit_target.py`.

sd07 carries a full dose-response ladder (0.025–100 µM, 12 levels) and the
within-ligand spread reaches **50×** (4-methylumbelliferone 0.5 vs 25 µM), so "any
hit" and "best hit" are genuinely different targets.

### 71a. The correction

Exact smallest library, by branch and bound, containing for **every** one of the
11 ligands:

| target | WT offered | vs Tian | forcing allowed | vs Tian |
|---|---|---|---|---|
| **any** known hit (§63b) | 9,216 | **15×** | 4,608 | 30× |
| **the BEST** known hit | **55,296** | **2×** | 27,648 | 5× |

**Demanding the best sensor rather than any sensor costs 6× in library size, and
Tian's 138,240 is then only 2× above the exact optimum.**

So §63b's "the entire headroom is 15×" and §65f's "PFAS is over-hedged by 3,441×"
were both measured against a target that accepts a sensor **50× weaker** than the
one actually available. ⚠ §71e revises this: the coumarin claim does collapse, but
the PFAS one **survives at 510×**. On the harder and more meaningful target, **Tian's
coumarin library is close to optimally sized.** The apparent waste was mostly an
artefact of my metric, not of their design.

This cuts the other way from what the question anticipated: the worry was that
Tian's library might be suboptimal and we were chasing a bad answer. The data says
the opposite — the library is nearly right, and it was my *scoring* of it that was
too lenient.

### 71b. And our designed libraries catch weak hits

| designed library | size | ligands with a hit | ligands with the **best** hit | potency of what it caught |
|---|---|---|---|---|
| class/profile f=0.50 | 17,280 | 3/11 | **1/11** | 5.0× weaker than best |
| class/profile f=0.35 | 103,680 | **8/11** | **2/11** | 4.5× weaker than best |
| class/profile f=0.25 | 129,600 | 3/11 | 1/11 | 5.0× weaker than best |

The 8/11 headline from §63d collapses to **2/11** on the best-hit target, and what
it does capture is ~4.5× weaker than what was available. **The two metrics are not
interchangeable**, and the harder one is the one worth designing against from here.

### 71c. What this does and does not settle

It does **not** find the best possible *sequence* — that needs binding data for
sequences nobody made, and nothing in this project predicts affinity (§40–§46).
What it does is enumerate every possible *library* exactly against the strongest
label available. The ground truth is still positive-unlabeled
([[feedback_hits_are_not_optima]]); it is now anchored to the best measured sensor
rather than to an arbitrary one.

### 71d. Audit of the sd07 dose-response parse (2026-08-26)

The dose block is laid out for a human reader, not a parser, so the reading was
checked rather than assumed. `results/secondary_library/sd07_dose_audit.txt`.

**Four cell values, not two**: `+` (244), `-` (231), **`ND` (273)** and **blank
(188)**. Blank means *below* the tested window and `ND` means *above* it — each
compound was assayed over its own contiguous window, chosen to bracket its
sensors:

| compound | tested window (µM) | best limit |
|---|---|---|
| Methoxsalen | 0.025 – 1.0 | **≤ 0.025** |
| Imperatorin | 0.025 – 1.0 | 0.05 |
| Osthole | 0.05 – 1.0 | 0.1 |
| Citropten | 0.05 – 2.5 | 0.25 |
| Bergapten | 0.05 – 2.5 | 0.5 |
| 4-Methylumbelliferone | 0.25 – 25 | 0.5 |
| Isopsoralen | 0.25 – 5.0 | 1.0 |
| 5,7-DH-4-methylcoumarin | 0.25 – 25 | 1.0 |
| Psoralen | 0.25 – 25 | 2.5 |
| 7-Methoxycoumarin | 0.5 – 25 | 5.0 |
| Scopoletin | 2.5 – 100 | 10.0 |

**The parse survives the audit**: of 78 clones, **77 are cleanly bracketed** (a `-`
immediately below the first `+`), **0 are non-monotone** (no `-` ever follows a
`+`), and **0 failed to respond** anywhere in their window. So "lowest dose scored
`+`" is a real detection limit, not an artefact of where the window happened to
start.

⚠ **One clone is left-censored**: Methoxsalen **12A-1** is `+` at 0.025 µM, the
lowest dose anyone tested it at, so its true limit is **≤ 0.025** rather than
`= 0.025`. It is that ligand's best sensor, so the error understates its margin
over the runner-up (0.05 µM) rather than inflating it — the direction that cannot
create a false §71 conclusion.

⚠ **Windows differ by compound**, so a limit is only resolved to the ladder used
for that compound. Scopoletin's 10 µM and Methoxsalen's 0.025 µM are both "best",
400× apart, and neither was tested outside its own window.

Also checked: the off-position filter (§63a's 8 PCR singletons) **removes no
ligand's best sensor** — best-overall equals best-within-the-18 for all 11.

### 71e. PFAS on the best-hit target — the headroom survives there

⚠ **sd09 encodes potency differently from sd07, and the difference matters.**
sd07 carries the full explicit ladder (12 columns, `+`/`-`/`ND`/blank). sd09
carries only a **coarse three-point primary screen** (1 / 10 / 100 µM, `+`/`-`
only) **plus a separate `min_conc (µM)` column from a finer retest** with nine
levels (0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500 µM). The two disagree on **98 of
154 rows**, always with `min_conc` lower — 50 µM where the primary screen first
scores `+` at 100 — and 3 rows carry values outside the three columns entirely.
`min_conc` is the authoritative number; deriving potency from the dose columns
would have been wrong for two-thirds of the PFAS sensors.

Exact smallest library containing, for every one of the 25 PFAS ligands:

| target | exact minimum | vs Tian's 49,545,216 |
|---|---|---|
| **any** known hit | 14,400 | **3,441×** |
| **the BEST** known hit | **97,200** | **510×** ⚠ superseded by §73: **43×** |

So the correction behaves completely differently in the two classes:

| class | any-hit headroom | **best-hit headroom** | |
|---|---|---|---|
| coumarin | 15× | **2×** | Tian sized it near-optimally |
| **PFAS** | 3,441× | **510×** | still enormously over-hedged |

**§65f's claim survives the metric correction.** The coumarin library really was
close to right and my 15× was an artefact; the PFAS library really was
over-hedged, and remains so by ~500× even when the target is the best measured
sensor rather than any sensor. The distinction tracks exactly what §65f proposed —
class-relevant round-1 data. Coumarin had 36 own-class round-1 clones and a clean
position signal; PFAS had 89 clones whose best position reached 0.61 and whose
residue identities never converged, against a chemically novel class. Less
confidence, wider library — and for PFAS the widening was far beyond what its own
results needed.

---

## 73. ⚠ Correction: the PFAS library is depth-capped, and the headroom was overstated ~12–18× (2026-08-26)

Jannis flagged 49,545,216 as looking high. The arithmetic was right — it is the
13 positions' `(1+n)` factors, and equals the three assembly blocks crossed
(384 × 144 × 896). **The model was wrong.**

### 73a. The tell

If a library is full-combinatorial, a random member carries a predictable number
of substitutions. Comparing that to what the sensors actually carry:

| library | expected per random member | observed in sensors | |
|---|---|---|---|
| Coumarin | 7.0 | **6.8** (range 3–10) | consistent |
| PFAS | 9.2 | **5.4** (range 2–7) | **not consistent** |

And the PFAS shape is decisive: **93 of 154 sensors carry exactly 6 substitutions
and exactly one carries 7.** Selection gives smooth distributions; that is a wall.
The PFAS library's own variant count peaks at **9** substitutions, so under a
full-combinatorial model the modal sensor should have 9 — none does.

Coumarin's sensors peak at 7, and its library's variant count also peaks at 7
(34,789 members). That one really is full-combinatorial.

**So PFAS is depth-capped at ~6 substitutions**, making its true size
**2,193,926 — 4.4 % of the product, a 22.6× overstatement.**

This is the same error as §53, where DSM-Hao was taken as 3.79 × 10²¹ before it
turned out to be a *double*-substitution library. I reprinted that same absurd
number in §72's table without flagging it.

### 73b. The corrected comparison, like for like

The oracle minimises a full-combinatorial product, so its menus have to be
re-costed under the same ≤6 cap — otherwise a capped denominator is being divided
into an uncapped numerator.

| class | library as built | any-hit oracle | | best-hit oracle | |
|---|---|---|---|---|---|
| coumarin (uncapped) | 138,240 | 9,216 | 15× | 55,296 | **2.5×** |
| **PFAS (≤6 cap)** | **2,193,926** | **11,424** | **192×** | **51,008** | **43×** |

| PFAS headroom | reported | corrected | overstated by |
|---|---|---|---|
| any-hit | 3,441× | **192×** | 18× |
| best-hit | 510× | **43×** | 12× |

### 73c. What survives

The qualitative split holds: **coumarin is near-optimally sized (2.5×), PFAS is
genuinely over-hedged (43×)**, and §65f's explanation — library width tracks
class-relevant round-1 data — is unaffected. What does not survive is the
magnitude. "Over-hedged by three orders of magnitude" was wrong; it is a bit over
one order.

⚠ **A standing check this project keeps failing:** before quoting a library size,
confirm whether it is full-combinatorial or substitution-depth-limited. The
diagnostic is cheap and now written down — compare the expected substitution count
of a random member against the observed distribution in that library's own hits,
and look for a hard ceiling.

---

## 74. Lee/Pellock/Baker NTF2 paper; why PFAS has more headroom; and what our method actually needs (2026-08-26)

### 74a. The paper — peer-reviewed, and the affinity-maturation half is the usable part

*Small-molecule binding and sensing with a designed protein family*, Lee, Pellock,
Norn et al., **Nature Communications 2026, 17:4533** — so unlike §66's L-Caliby
this one is peer-reviewed and has real binding data.

>10,000 designed NTF2 scaffolds (hallucination + ProteinMPNN + Rosetta), RIFdock
to place six ligands, sequence design by native-guided Rosetta or LigandMPNN,
filtered on Rosetta ddG / H-bonds / contact molecular surface plus AF2. Then yeast
display.

**The screening burden is the number to take away:**

| target | oligos ordered | unique hits | ITC-characterised |
|---|---|---|---|
| **HCY (cortisol)** | **630** | **1** | 1 |
| ROC | 1,661 | 19 | 1 |
| APX | 9,024 | 8 | 4 |
| OHP | 7,573 | 117 | 4 |
| WRF | 16,276 | 46 | 2 |
| IRI | 19,390 | 8 | 2 |

~54,500 designs → 199 hits (**0.36 %**), best affinities high-nM to low-µM. Their
own honest note: the success rate "was lower than alternative methods utilizing the
four-helical bundles", and the designs that worked hugged native ketosteroid
isomerase structures — i.e. it worked where it stayed close to nature.

**What we can actually use is Fig. 4.** From **hcy129 — the single cortisol hit out
of 630** — they ran **site-saturation mutagenesis**, took the favourable mutations,
built a **combinatorial library of just those**, screened it, and got **hcy129.1 at
K_D 68 nM, a 31-fold improvement**. The same move on iri807 gave 1.5–5.5×.

That is our two-stage problem exactly — hit → per-position profile → focused
combinatorial library — executed from **one** starting hit. Which answers §74c.

### 74b. Why PFAS has ~43× headroom and coumarin only 2.5×

Not because the PFAS library wastes options — it barely does. Only **4 of 45**
offered substitutions never appear in any PFAS sensor (9 %), and for coumarin it is
**0 of 23**. Both libraries are well utilised at the substitution level.

The driver is that **library size is exponential in position count, so a modest
per-position excess compounds**:

| class | positions | offered | oracle needs | per position | ratio | ratio^positions | observed |
|---|---|---|---|---|---|---|---|
| coumarin | 11 | 23 | 18 | 2.09 → 1.64 | 1.17 | ~6× | **2.5×** |
| PFAS | 13 | **45** | 25 | **3.46 → 1.92** | **1.53** | ~244× | **43×** |

PFAS offers **1.5 extra substitutions per position over 13 positions**; coumarin
offers **0.45 over 11**. Neither looks extravagant per position — but 1.53¹³ is
244 while 1.17¹¹ is 6. (The estimate overshoots the measured 43× because the depth
cap and the uneven spread of substitutions across positions both damp it; the
mechanism and direction are right, the arithmetic is indicative.)

So the reason PFAS is loose is **not** bad judgement about which residues to
include — it is that hedging by ~1.5 residues per position is invisible locally
and enormous globally. That is also why our method has more to offer there: the
prize is trimming per-position depth, which is exactly what a residue-ranking step
does.

### 74c. What our method is, and it needs MORE than one hit

Stated plainly, the pipeline as it now stands:

| step | how | evidence |
|---|---|---|
| 1. which positions to open | rank the 18 pocket positions by round-1 mutation rate | §63c, §65g |
| 2. which residues per position | class-weighted round-1 frequency, keep within a fraction of the top | §63c |
| 3. which combinations to keep | Cartesian FastRelax ΔΔG, applied through its additive per-substitution fit | §67, §69 |

**Tested directly: can step 2 run off a single round-1 hit?** Recovering Tian's
coumarin design from different amounts of evidence:

| evidence | positions found | optimal residue ranked 1st |
|---|---|---|
| A. full class-weighted (49 clones, ≥0.4 Tanimoto) | 10/11 | **7/11** |
| B. pooled round-1, **ligand-blind** (692 clones) | 10/11 | 4/11 |
| C. **one own-ligand hit** + pooled prior (200 draws) | 9.9/11 | **4.1/11** |

> **A single hit adds essentially nothing over the ligand-blind baseline** (4.1 vs
> 4.0 residues). One clone carries 2–3 substitutions; an 11-position residue
> profile cannot be built from that.

**Positions survive on one hit — residues do not.** Step 1 draws on pooled
round-1 and is ligand-blind anyway, so it gives 10/11 regardless. Step 2 is the
part that needs a class, and one clone is not a class.

Two consequences:

- For a target like **PFAS**, where the prize is mostly position selection (§73),
  a single hit may be enough to be useful.
- For a target like **coumarin**, where the prize is residue identity, it is not —
  and the fix is the one §74a's paper demonstrates: **run SSM on the single hit**
  to generate the per-position profile experimentally, then apply steps 2–3 to
  that. That converts one sequence into exactly the input our method wants, at the
  cost of one extra round of wet-lab work.

---

## 75. Can we go from ligand to a few sequences? The target is the wrong shape (2026-08-26)

Jannis: we always modify the same pocket and have ~1,150 labelled clones, so could
an algorithm output a few sequences, or a small library, from a ligand alone? And
does §74a's Baker result argue against it?

### 75a. "A few sequences" is ill-posed, and that is a measurement not an opinion

There is no single answer per ligand to predict. Across the 78 coumarin sensors:

| | |
|---|---|
| mean within-ligand Jaccard between sensors | **0.45** |
| identical sensor pairs, out of all within-ligand pairs | **1** (one duplicate, Scopoletin) |
| distinct substitutions spanned by one ligand's own sensors | **10–16** |
| positions those sensors touch | **8–11** |
| substitutions present in **every** sensor for a ligand | **1–7** |

Every ligand has a *manifold* of solutions sharing under half their content. Asking
a model for "the" sequence is asking it to pick one winning ticket out of many
winning tickets. **The correct output format is a library, not a shortlist** — and
that is a property of the biology, not a limitation of our methods.

The constructive half: the **core is real and predictable**. **V163W appears in
every sensor for 10 of 11 ligands**; V83L in 6. Those are exactly the substitutions
§64d's forcing analysis identified from round-1 alone.

### 75b. What a per-ligand library actually costs, LOLO

The real use case is one target, not eleven. Deleting each ligand's own round-1
clones and designing from the rest:

| | median size |
|---|---|
| oracle, if you knew one specific working sensor | **32** (any hit) / 128 (best) |
| oracle, to cover a ligand's whole known solution union | **2¹⁰–2¹⁶ ≈ 10³–10⁴** |
| **our designed library, containing ≥1 hit** | **69,120** — works for **9/11 ligands** |
| our designed library, containing the **best** hit | 9.2 M — works for 5/11 |

So the gap to "knowing the answer" is ~2,000×, but the gap to the **realistic
floor** — covering the solution manifold rather than one point in it — is only
**~10–50×**.

### 75c. Why the Baker result does not transfer

Their hard problem was the one we do not have. They designed **>10,000 novel
backbones**, docked ligands into pockets that did not exist, and needed sequences
that both fold and bind — 54,500 oligos for a 0.36 % hit rate, and cortisol gave
**1 hit from 630**. Their own retrospective says the designs that worked hugged
native ketosteroid isomerase.

We change 11–18 side chains in a backbone that already folds, already binds, and
already transduces. Foldability is near-free (§67's ΔΔG separation is a
*within*-menu effect, not a folding rescue), the pose problem is bounded, and the
readout is a **growth selection** rather than yeast display plus FACS plus deep
sequencing.

Which flips the binding constraint. **Screening 10⁵ PYR1 variants is one flask;
screening 10⁵ novel backbones is a campaign.** So for us library size is not the
scarce resource — Y2H handles 10⁶–10⁷ — and shrinking 69,120 to 5,000 buys very
little. What matters is whether the library **contains a sensor at all**.

### 75d. So: yes to a program, no to a shortlist

| ask | verdict |
|---|---|
| ligand → a few sequences | **no** — the solution manifold has Jaccard 0.45 and no unique target |
| ligand → a ~10⁴–10⁵ library containing a sensor | **already ~there**: 9/11 ligands at 69,120, LOLO |
| ligand → a library containing the **best** sensor | **no** — 5/11 at 9.2 M |
| ligand → a small library **without any hit for that ligand** | positions yes, residues no (§74c) |

The honest position is that we are at **parity with the incumbent**, not ahead of
it, and the remaining prize on library size is ~10–50× rather than orders of
magnitude. The two levers are known and both are measured: residue ranking
(currently 7/11 at rank 1, §63c) and the ΔΔG combination filter (~2×, §69d).

⚠ **And the real risk is not library size at all.** Everything above is measured
on coumarins and PFAS, whose chemistry is inside the 144-substitution vocabulary
round 1 already explored (§63a). Nothing here shows that a genuinely novel
chemotype has *any* solution in that vocabulary. That — not shortlist length — is
what would sink a prospective target, and it is what §48d's sealed prospective
test was designed to probe.

---

## 76. FastRelax reproduces the pocket but breaks the salt bridge — it is SCORING (2026-08-26)

Jannis asked the question none of the ligand-present runs had asked: **does relax
reproduce the experimental pocket when handed the right sequence and ligand?** If
it does, our failures are in ranking, not modelling; if it does not, the problem is
the score and an iterative dock/relax/mutate loop would just converge on the wrong
answer faster. `scripts/133_relax_vs_crystal.py`.

Every previous ligand-present protocol asked the protocol to *rank* or *design* —
FastDesign (§23j), coupled moves (§36), MM-GBSA (§37–42), TI (§43–46). **None
tested structural recovery.** §67/§69's FastRelax was protein-only by construction.

### 76a. Design

PYR1^MANDI in the stage-1 frame (4WVO's four substitutions, its side chains
transplanted after a 696-atom backbone superposition at 0.53 Å), with the crystal
mandipropamid pose — `stage1/params/3UZ_0001.pdb` matches `wt_mandi.pdb`'s 3UZ to
**0.0000 Å**, so it *is* the crystal ligand in params-compatible naming.

| arm | start | asks |
|---|---|---|
| **A** | crystal rotamers | do they **stay**? → scoring |
| **B** | repacked from scratch | do they **reach** it? → sampling |

### 76b. Relax converges to one minimum, and it is not the crystal

| | mean over 18 pocket positions | the four mutated |
|---|---|---|
| A (from crystal, relaxed) | **0.61 Å** | **1.08 Å** |
| B (repacked, relaxed) | **0.61 Å** | **1.07 Å** |

**Maximum |A − B| across all 18 positions: 0.035 Å.** Relax finds the *same*
conformation whether started at the crystal or from a scrambled repack. B's start
was genuinely scrambled (L159 at 2.88 Å, N167 at 1.58 Å before relax) and still
converged to A.

So **sampling is not the bottleneck** for reaching the score's minimum — and that
minimum sits 0.6 Å from the crystal overall, 1.1 Å at the positions we care about.

⚠ **Control: backbone frame mismatch does not explain it.** Spearman(per-residue
backbone deviation, side-chain RMSD) = **+0.13**. L117 has 0.24 Å backbone
deviation and 1.50 Å side-chain error; I110 has 0.53 Å backbone and 0.26 Å side
chain. The one position where the caveat does bite is **K59 itself** (backbone
deviation 1.15 Å, the largest in the set), so its 1.22 Å is partly frame.

### 76c. The mechanism: bidentate becomes monodentate

The RMSD understates what is lost. ARG59 → ligand polar contacts, measured **within
one structure** so the frame mismatch cancels:

| | contacts |
|---|---|
| crystal (transplanted) | **NE–O2 2.64 Å, NH1–O2 3.26 Å** — bidentate |
| after FastRelax | NH1–O2 2.78 Å, **NE–O2 3.96 Å** — monodentate |

Two H-bonds under 3.5 Å become one. The crystal geometry reproduces §33a exactly
(2.64 / 3.26 Å), and **ref2015 walks away from it.**

### 76d. What this settles

> **It is the score, not the sampling.** Relax converges reliably to a well-defined
> minimum from any start; that minimum discards the bidentate guanidinium
> interaction the real sensor uses.

This refines rather than contradicts §33b. That section found crystal Arg59's
χ3 ≈ 100° is non-rotameric so the library never *proposes* it — true, and now we
know the deeper problem: **even handed the answer, the score does not keep it.**
Both failures point at the same interaction, and both are consistent with §23j,
where ref2015 saw the K59 salt bridge and overcharged its desolvation by +10.1 REU.

**Consequences for the iterative loop (§76's originating question):** iterating
dock → relax → mutate with ref2015 converges on ref2015's preferred pocket, which
is the wrong one at exactly the positions that distinguish a sensor. Adding more
iterations cannot fix a fixed point. What would help is the score.

**And a general score may not suffice.** §68c measured that ref2015 prefers GROW
over SHRINK (β +0.43 vs +1.75), and §69c that the required direction **flips with
ligand size** (Spearman(ligand heavy atoms, ΔVolume used) = −0.30). A single
global reweighting cannot satisfy both regimes, so the tractable form is
**class-conditioned** weights — binned by ligand size and polarity — rather than
one PYR1 function. The data supports ~208 independent ligands (§55b), which is
enough for a handful of parameters per bin and not enough for a general model.

---

## 77. Ligand-aware scoring is WORSE than protein-only (2026-08-26)

The gate test from §76's originating question, run to completion: 670 variants,
each with its cognate coumarin docked into the relaxed mutant pocket and protein +
ligand relaxed jointly, scored `dG_bind = E(complex) − E(protein) − E(ligand)`.
Job 27823020, 40/40 chunks, **zero failures**. `scripts/131_ligand_aware_score.py`.

### 77a. Adding the ligand destroys the discrimination

| comparison | protein-only (§69) | **ligand-aware** |
|---|---|---|
| REAL vs LIBRARY | **0.284** (p 2 × 10⁻⁸) | **0.471** (p 0.5) |
| REAL vs LIBRARY, n_sub-matched | **0.250** | **0.504** |
| REAL vs WILD | 0.117 (p 2 × 10⁻²³) | 0.332 (p 1 × 10⁻⁵) |

**Exactly chance on the comparison that matters.** The score still separates REAL
from WILD, so it is not pure noise — it can tell library residues from arbitrary
ones. It cannot tell a working combination from a random one inside the library,
which is the entire job.

### 77b. Why — the noise is several times the signal

| diagnostic | result |
|---|---|
| Spearman(ligand heavy atoms, dG_bind) | **−0.39** — partly measuring ligand size, not fit |
| within-ligand AUC (ligand identity removed) | **0.417** weighted mean — better than pooled, still far off 0.250 |
| sd of dG_bind among REAL sensors **for the same ligand** | **1.4 – 5.2 REU**, range to 16 REU |
| between-set difference in median dG_bind | **~1 REU** (REAL −20.64 vs LIBRARY −19.72) |

**The spread among sensors that all demonstrably work is three to five times the
difference between sensors and random variants.** No amount of careful reading
recovers a signal from that.

Per-ligand AUC swings from 0.25 (methoxsalen, good) to **0.78 (citropten,
backwards)** — consistent with a per-variant quantity dominated by pose noise
rather than by fit.

### 77c. The pre-registered ambiguity, and how far it can be resolved

§131's header committed in advance that a null here would be ambiguous between
"ligand-awareness adds nothing" and "the docking is too poor to tell". The
diagnostics narrow it but do not close it:

- **Docking noise is real and large** — one smina pose per variant, no replicates,
  and the within-ligand spread is exactly what that predicts.
- **§76's scoring defect is also real and independent** — handed the crystal
  complex, relax breaks the bidentate ARG59 salt bridge into a monodentate one.

So both candidate causes are demonstrated, and this experiment cannot apportion
them. The cheap next test would be best-of-N docking (say 5 poses per variant,
taking the lowest dG_bind), which would collapse the pose component and leave the
scoring component exposed. **§76 makes the prior poor**: even with a perfect pose,
the score walks away from the interaction that distinguishes the real sensor.

### 77d. Where this leaves the loop

An iterative dock → relax → mutate cycle is now doubly discouraged. §76 showed it
converges on ref2015's fixed point, which is the wrong conformation at the mutated
positions. §77 shows the ligand-aware objective it would optimise is **at chance**
on the discrimination task. Iterating a noisy score toward a wrong fixed point is
not a plan.

**What survives is the protein-only filter** (§67, §69): AUC 0.250 n_sub-matched,
~2× at 90 % sensor retention, in both the coumarin and PFAS regimes. It works
*because* it never touches the ligand, and so never incurs the pose noise.

---

## 78. The score gets the cross-over right on crystal poses — and relax destroys it (2026-08-26)

Jannis: does it work when we hand it the right pose, and what if we skip the relax?
`scripts/134_crystal_crossover.py`. A receptor × ligand cross-over so the diagonal
must win in **both** directions — getting one right is easy for the wrong reason
(a mutant pocket is emptier), getting both means complementarity rather than a size
artefact. All four cells built identically in the stage-1 frame, so none is a
crystal while another is a model.

### 78a. Right at raw and repack, wrong after relax

`dG_bind = E(complex) − E(protein) − E(ligand)`, same protocol per row:

| cell | raw | repack | relax |
|---|---|---|---|
| WT + mandipropamid | **+1994.0** | +1497.8 | −37.5 |
| PYR1^MANDI + mandipropamid | +6.9 | +11.3 | −36.8 |
| WT + ABA | −26.7 | −23.9 | −67.1 |
| PYR1^MANDI + ABA | +191.1 | −19.9 | −59.2 |

| selectivity | raw | repack | relax |
|---|---|---|---|
| mandipropamid (quad should win) | **−1987** ✅ | **−1487** ✅ | **+0.8 ❌** |
| ABA (WT should win) | **+218** ✅ | **+4.1** ✅ | +7.9 ✅ |

**Jannis's hypothesis is confirmed: relaxation destroys the mandipropamid signal.**
WT + mandipropamid starts at +1994 REU — F108 sits 0.62 Å from a ligand heavy atom
(§46c) — and relax drives it to **−37.5**, indistinguishable from the real sensor.
Relaxation lets the wild-type pocket make room for a ligand it cannot actually
accommodate. That is §37b's "F108A comes out NULL because relaxation absorbs the
clash it exists to relieve", now shown for the whole cross-over.

### 78b. ⚠ But the signal is one clash, not complementarity

Single-mutant decomposition against WT, same crystal pose:

| variant | Δ(raw) | Δ(repack) |
|---|---|---|
| K59R | +0.39 | **−3.47** |
| V81I | −19.91 | +0.82 |
| **F108A** | **−1751.08** | **−1478.81** |
| F159L | −216.54 | −0.02 |
| quad | −1987.14 | −1481.50 |

**F108A alone is 99.8 % of the quad's repack signal. K59R is 0.23 % of it.**
F159L, worth −217 raw, vanishes to −0.02 once side chains repack.

So the cross-over "success" is a single steric collision — the same thing §32a
detected geometrically without any energy function, which named F108 (2.78 Å) and
F159 (1.42 Å) unprompted. **The score is not seeing complementarity; it is seeing
that phenylalanine does not fit.** K59R, the substitution every method has missed,
is still invisible: −3.47 REU against F108's −1479, though it is at least the right
sign and above §69's 1.34 REU replicate noise.

### 78c. What this changes

- **Do not relax before scoring for discrimination.** The information lives in the
  unrelieved clash, and relax is very good at relieving it.
- **But an unrelaxed score is a clash detector**, and we already had one that is
  cheaper and needs no energy function at all (§32).
- The honest scope is unchanged from §68: shape gets recovered, chemistry does not.

⚠ Cosmetic defect in the decomposition output: the label prints the first letter of
the three-letter code, so K59**R** appears as "59A" (ARG). Values are correct.

---

## 79. The 2×2 in PYL2, and where the signal actually lives (2026-08-26)

Jannis supplied 7MWN (an engineered PYL2 WIN 55,212-2 sensor) and then a wild-type
PYL2 + ABA control, which together give a receptor × ligand cross-over **inside one
system with one numbering**. `scripts/135_win_crossover.py --system {win,aba}`.

⚠ **3KDJ was rejected and would have been a silent wrong-receptor error.** It is
**PYL1 + ABI1 + ABA**, not PYL2 — chain A is **27 % identical** to PYL2 and numbered
31–209. Applying K64Q/F165A/V166I there would have mutated R64/I165/W166 and
returned a confident cross-over on the wrong protein. **3KDI** is the right control:
wild-type PYL2 + ABA, 181 residues 7–187, **zero gaps, 100 % identity, same
numbering as 7MWN**. 7MWN's own deposited entity record independently confirms the
mutation set: `'K64Q, F165A, V166I'` against UNP PYL2_ARATH O80992.

### 79a. Six of six correct — and each on one substitution

| column | expected winner | raw | repack | relax |
|---|---|---|---|---|
| WIN (7MWN) | sensor | **−46.00** ✅ | **−2.44** ✅ | **−1.44** ✅ |
| ABA (3KDI) | wild-type | **−35.83** ✅ | **−0.67** ✅ | **−2.83** ✅ |

The diagonal is right everywhere, including after relax — better than PYR1's
mandipropamid arm, which flipped (§78a). But the decomposition repeats §78b's
pattern exactly, at a different position each time:

| WIN, revert one | raw | repack | | ABA, install one | raw | repack |
|---|---|---|---|---|---|---|
| **Q64K** | **+45.37** | +0.84 | | K64Q | −3.11 | −0.20 |
| A165F | +0.23 | +0.12 | | F165A | +0.84 | +0.82 |
| I166V | +0.39 | +1.48 | | **V166I** | **+38.11** | +0.05 |
| all three | +46.00 | +2.44 | | all three | +35.83 | +0.67 |

**Q64K is 98.6 % of the WIN signal; V166I is more than all of the ABA signal** (K64Q
points the wrong way there). One collision per column, and at repack the totals fall
to **−2.44 and −0.67 REU** — at or below §69's 1.34 REU replicate noise.

**Why K64 is visible where PYR1's K59 was not:** wild-type K59 sits **2.86 Å** from
mandipropamid and does not clash (§32a), so reverting it costs +0.39 REU. Wild-type
K64 **collides** with WIN's 32 heavy atoms, so reverting it costs +45 REU. The score
sees the lysine for the same reason it saw F108 — steric overlap, not chemistry.

### 79b. The score is one term

Per-term decomposition of the WT-minus-sensor difference (mandipropamid):

| protocol | **fa_rep** | fa_atr | fa_sol | fa_elec | hbond_sc | total |
|---|---|---|---|---|---|---|
| raw | **1995.12 (100.4 %)** | −5.90 | +4.39 | −5.78 | +0.29 | 1988.12 |
| repack | **1484.69 (100.2 %)** | −2.84 | +2.56 | −3.44 | +0.11 | 1481.07 |
| relax | **+1.76** | −2.13 | −0.99 | +1.99 | −0.84 | **−0.04** |

**Unrelaxed, the discrimination is 100 % steric repulsion.** Attraction, solvation,
electrostatics and hydrogen bonding contribute nothing measurable, and `fa_atr` and
`fa_elec` point the **wrong way** — they mildly favour the cell that cannot bind.

After relax `fa_rep` falls from 1995 to **1.76**, a thousandfold, and the residual
is competing terms of 1–2 REU that cancel to **−0.04**. The post-relax signal is not
a weaker version of the same thing; it is gone, replaced by term-level noise.

### 79c. The correction this forces on how we describe it

FastRelax does not *overpack*. It **over-relieves**: WT + mandipropamid goes
+1994 → −37.5 REU, landing indistinguishable from the real sensor — a structure that
does not exist, since wild-type PYR1 genuinely cannot bind mandipropamid. Because
the discriminating information lives entirely in unrelieved `fa_rep`, relaxing it
away destroys the measurement.

> **The score is a clash detector with a working sign and no usable magnitude.** It
> is right when a collision exists, silent when one does not, and the moment the
> structure is relaxed enough to be physical the collision — and the answer — is
> gone.

That is the same boundary §32 reached geometrically with no energy function at all,
now confirmed term by term, on two receptors, two ligand classes and four cells.

---

## 80. 718 Boltz-2 sensor structures: good poses, no potency signal (2026-08-27)

Jannis supplied `boltz_predictions_719sensors.zip` — Boltz-2 co-folded models of the
Tian sensors, June 2025, **one seed each**. This is the structural dataset §79
said we lacked, and it changes what can be tested.

### 80a. Recovering the mapping

The archive carries no input manifest, so sensor identity was recovered from the
CIFs themselves: the consensus over all 718 sequences reproduces wild-type PYR1
exactly (191 aa + an `SGDGSGSQVT` linker), each structure's substitutions follow by
diff, and the ligand is identified by heavy-atom count.

| | |
|---|---|
| structures extracted | **718** (1,436 files; MSAs left in the archive) |
| matched to an sd03/sd07/sd09 clone by substitution set | **701** |
| **uniquely resolved** (substitution set **+** ligand heavy-atom count) | **637** |
| distinct ligands | **185** |
| with a measured potency | **590** (1 µM: 55, 10 µM: 257, 100 µM: 278) |

⚠ 57 remained ambiguous (one substitution set used for several ligands of the same
size) and 24 had no size match. Both are dropped rather than guessed.

### 80b. The poses are good — where we can check them

From §55a's co-folding run, against crystals:

| system | ligand RMSD | ligand_iPTM | seed spread (5 models) |
|---|---|---|---|
| PYR1^MANDI + mandipropamid | **0.49–0.66 Å** | 0.985–0.988 | **0.35 Å** |
| WT PYR1 + ABA | 0.72–1.19 Å | 0.937–0.954 | 0.91 Å |
| **WT PYR1 + mandipropamid** (cannot bind) | 1.28–1.63 Å | 0.963–0.976 | **1.70 Å** |

Sub-Ångström on the real complex — far better than the smina poses whose noise
destroyed §77. **But confidence does not know about binding**: wild-type PYR1 +
mandipropamid, a complex that does not form, scores ligand_iPTM **0.97**. Seed
spread separates them 4.9× better than confidence does (0.35 vs 1.70 Å), matching
the binder project's finding that iPTM reproduces while the pose may not.

⚠ 4WVO and 3K3K are both in the PDB and so in training; the one case that is not
(WT + mandipropamid) looks like 4WVO's mode copied onto the wild-type sequence.
These are not clean tests.

### 80c. ⚠ Nothing structural predicts potency — 590 sensors, 172 ligands

Spearman against measured min_conc (lower = better sensor). "Within-ligand" removes
ligand identity by correlating only inside each ligand's own sensor set (43 ligands
with ≥ 4 sensors and potency variation, 328 sensors):

| metric | pooled ρ | within-ligand ρ | ligands in the right direction |
|---|---|---|---|
| Boltz confidence | −0.019 | −0.030 | 20/43 |
| ligand_iPTM | −0.112 | −0.084 | 23/43 |
| complex_pLDDT | +0.009 | −0.047 | 20/43 |
| protein–ligand contacts | −0.145 | −0.029 | 22/43 |
| polar contacts (H-bond proxy) | +0.024 | **+0.089** | 13/43 |
| buried fraction | +0.117 | **+0.061** | 13/43 |
| max / summed VdW overlap | −0.012 / −0.029 | −0.006 / −0.047 | 20/43, 22/43 |
| **ligand heavy atoms** | **−0.246** (p 2×10⁻⁹) | — | — |

**Chance is 21.5/43.** Every descriptor sits on it, and the two most chemically
meaningful — hydrogen bonds and burial — point the **wrong way**. The only strong
pooled correlate is **ligand size**, which is a confound, not a design signal, and
it explains most of the pooled numbers that look non-zero.

### 80d. What this does and does not close

It does **not** say the structures are wrong — poses validate at 0.5 Å where a
crystal exists. It says that **given good structures, neither co-folding confidence
nor simple structural descriptors rank sensor potency.**

Two caveats that matter for how hard to read it:

- **The task here is harder than the one our filter does.** These 590 are all
  *working sensors*; ranking potency *among positives* is a different and harder
  problem than separating positives from random library members, which §69 does at
  AUC 0.25. A null here does not overturn §69.
- **The label is coarse** — three levels (1/10/100 µM) from a round-1 screen, and
  many ligands have little potency variation among their sensors.

**The decisive experiment this dataset finally makes possible** is the one §79 said
we could not run: co-fold a matched set of *random library variants* — non-sensors —
and repeat §69's discrimination test with Boltz poses instead of smina poses. That
separates "the poses were too poor" from "the score cannot see it", which §77 could
not apportion. It needs GPU time and 300–600 co-folding runs.

---

## 81. The umbrella calibration FAILS — and my own repair is a prime suspect (2026-08-27)

All 62 windows reached 20.0 ns (1.24 µs total) after §70's reseed. `113_us_pmf.py`:

| arm | G(closed) | G(open) | G_open − G_closed | |
|---|---|---|---|---|
| **HOLO** (WT + ABA) | −1.70 | +3.32 | **+5.02** kcal/mol | CLOSED favoured ✅ |
| **APO** (WT) | −1.26 | +5.65 | **+6.92** kcal/mol | CLOSED favoured ❌ |
| | | | **ΔΔG = −1.90** | ❌ |

Pre-registered (§60c): holo favours closed ✅, **apo favours open ❌**, **ΔΔG
positive ❌**. **Two of three fail.**

### 81a. The failure is not technical — the easy explanations are all excluded

| check | result |
|---|---|
| window overlap | **no window has < 5 occupied bins** |
| coordinate equilibration | drift **0.02–0.03 Å**; windows sit on target within 0.05–0.09 Å |
| does ABA dissociate in the open windows? | **no** — 52–82 heavy-atom contacts < 4.5 Å even at 20 Å |
| holo convergence | half-split drift **0.22** kcal/mol |
| apo convergence | half-split drift **1.05** kcal/mol ⚠ |

So the restrained coordinate is sampled properly, the ligand stays bound, and the
windows overlap. The number is what it is.

### 81b. ⚠ The §70 repair destroyed the control that would have caught this

The original design (§60, `111_us_seed.py`) seeded windows from **both basins** —
closed-basin trajectories for the low windows, open-basin for the high ones —
specifically *"so the low and high arms approach the barrier from opposite
directions and disagreement in the overlap region is visible."*

§70's repair replaced every open-basin seed (11.0–20.0 Å, **both arms**) with
**adiabatic pulling outward from the closed basin**. It made the windows run — but
now **every open window descends from the closed state**, one-directionally, over
3.8 ns. The hysteresis control is gone, and a directional bias from pulling would
inflate the open-state free energy **in both arms** — which is exactly the
observed pattern.

⚠ `113_us_pmf.py`'s own docstring promises a hysteresis diagnostic that **the code
never implemented** and that the data can no longer support. I chose pulling
because it was the cheap fix; the expensive fix — re-solvating the open frames into
the matching topology — was the one that would have preserved the control.

### 81c. The other candidate, and why it does not rescue the result

The apo reference may simply be mis-specified. **Apo-open PYR1 is a dimer state**
(3K3K chain A, §28g), and we simulated a **monomer** — and §29b/§30 already found
the apo-closed monomer stable across 3 × 300 ns with the **most rigid gate of all
fifteen units**, with a "dimer confound" flagged at 2.9×. On that reading the PMF's
apo answer may be correct for the system actually simulated.

**But that does not save it.** Even granting a monomer that prefers closed, ABA
should make it prefer closed **more**, not less. ΔΔG = −1.90 kcal/mol says ligand
binding *destabilises* the closed state relative to open, which is backwards on any
reading of the mechanism. That failure is independent of the dimer question.

### 81d. Verdict

> **The scheme is not validated. No designed-pocket number may be quoted from it.**

Two repairs are needed before it could be, and they are separable:

1. **Restore both-basin seeding** — re-solvate the S10/S1 open frames into the
   S2/S9 topologies so the open windows descend from the open state, and implement
   the hysteresis check §113 promises. This tests §81b directly.
2. **Match the reference system** — run the apo arm as the dimer, or re-specify the
   pre-registration for a monomer, so the expectation matches what is simulated.

Until 1 is done the result cannot be interpreted, because the most likely
explanation for it is an artefact I introduced.

---

## 82. Boltz-2's binder classifier calls 90 % of real sensors non-binders (2026-08-27)

BoltzMol-1 (bioRxiv 2026.07.04.736485) builds hit discovery on Boltz-2, which
co-folds a complex and emits **a binder-vs-decoy classifier** plus a pIC50 head,
reporting **EF ≈ 18× at the top 0.5 %** on MF-PCBA. That classifier is exactly the
discrimination §69 measures, so it is the strongest untried idea we have.

**No GPU was needed.** The June 2025 archive already contains the affinity output —
2,870 files — so the classifier had been run on every sensor and never read.

### 82a. The result

`affinity_probability_binary` over the **637 experimentally confirmed sensors**,
every one isolated in a growth selection and validated by dose–response:

| | |
|---|---|
| median binder probability | **0.303** |
| called a binder at p > 0.5 | **9.6 %** (61/637) |
| at p > 0.3 | 50.9 % |
| at p > 0.2 | 81.6 % |

**Boltz-2 calls 90 % of real PYR1 sensors non-binders.**

### 82b. There is a signal, and it is at the wrong level

| min_conc | n | median p_bind | called binder |
|---|---|---|---|
| **1 µM** | 55 | 0.345 | **20.0 %** |
| 10 µM | 257 | 0.290 | 11.3 % |
| 100 µM | 278 | 0.301 | **5.8 %** |

The "called binder" rate is monotone in potency, and pooled AUC separating 1 µM
from 100 µM sensors is **0.634** — above chance.

⚠ **But within-ligand it vanishes**: ρ = **−0.028** (20/43 ligands in the right
direction, chance is 21.5). So the pooled signal is **between ligands**, not between
sensors — some ligands both yield better sensors and score higher. For design we
need discrimination *among variants for one ligand*, and there it is at chance.
`affinity_pred_value` behaves identically: pooled ρ = **+0.327** (p 2 × 10⁻¹⁵),
within-ligand ρ = **−0.000**.

### 82c. Why — the direction is untested by construction

BoltzMol-1 and Boltz-2's affinity training vary the **ligand** against a fixed
protein: that is what ChEMBL-style data contains, and what EF ≈ 18× measures. **We
vary the protein — 11–18 pocket mutations — against a fixed ligand.** Nothing in
the training distribution constrains that direction, and this is the first
measurement of it we know of.

Two things are separable in the failure and only one is fatal:

- **Calibration** — these sensors bind at 1–100 µM, weak by the standards of the
  affinity data Boltz-2 was trained on, so a low absolute probability is defensible.
- **Ranking** — being unable to order variants *within* one ligand is not, and that
  is the property a design filter needs.

### 82d. Where that leaves the scoreboard

| method | discriminates real sensors from random library variants? |
|---|---|
| protein-only Cartesian FastRelax ΔΔG (§69) | **yes** — AUC 0.250 n_sub-matched, ~2× at 90 % retention |
| ligand-aware Rosetta with docked pose (§77) | no — AUC 0.504 |
| cavity-volume match (§69) | no — AUC 0.540 |
| Boltz-2 confidence / iPTM (§80) | no — within-ligand ρ −0.03 to −0.08 |
| **Boltz-2 binder classifier (§82)** | **no** — within-ligand ρ −0.028 |

The protein-only clash filter remains the only thing that works, and §79 explains
why: it reads unrelieved `fa_rep`, which is the one quantity that actually differs
between a sensor and a random variant. Every richer method tested has been at
chance.

---

## 83. Is there a ligand class Boltz handles better? Calibration yes, ranking no (2026-08-27)

Jannis asked whether Boltz-2 does better on some chemical class — non-polar, say —
and whether it is even worth looking. `scripts/137_boltz_by_ligand_class.py`.
Descriptor list fixed before running; two questions kept separate.

### 83a. Calibration: yes, and not in the direction one would guess

Median `p_bind` per ligand vs eight descriptors, 171 ligands:

| descriptor | ρ | p |
|---|---|---|
| **aromatic rings** | **+0.387** | **4.4 × 10⁻⁷** |
| **H-bond donors** | **+0.262** | **6.4 × 10⁻⁴** |
| cLogP | +0.124 | 0.11 |
| rotatable bonds | −0.132 | 0.085 |
| heavy atoms, TPSA, HBA | −0.07 … +0.07 | ns |

**Boltz is most confident on aromatic, hydrogen-bond-donating ligands** — the
drug-like region where ChEMBL is densest — **not** on non-polar ones. That is a
real and interpretable property of the model. It is a statement about *calibration*,
not about whether the score works: it cannot pick a variant, because it is constant
across every variant for a given ligand.

### 83b. Ranking: a subgroup passed the permutation test and failed on inspection

Within-ligand ρ across the 43 ligands with ≥ 4 sensors and potency variation:
median **+0.000**, 20/43 negative against a chance rate of 21.5. Splitting into
tertiles, the best cell was **high H-bond donors: median ρ −0.528 (n = 9)**.

A permutation test that shuffles the 43 ρ values and re-takes the most extreme of
21 cells — so the multiple comparison is corrected — gave **p = 0.013**. It
survived.

**It should not have.** Four of the nine ligands in that cell are capsaicinoids:

| pair | Tanimoto | ρ A | ρ B | difference |
|---|---|---|---|---|
| capsaicin / zucapsaicin | **1.00** | +0.488 | +0.183 | 0.305 |
| capsaicin / dihydrocapsaicin | 0.71 | +0.488 | −0.528 | **1.016** |
| capsaicin / nonivamide | 0.68 | +0.488 | −0.569 | **1.056** |
| dihydrocapsaicin / nonivamide | 0.76 | −0.528 | −0.569 | 0.041 |

**Molecules with an identical ECFP4 fingerprint give ρ differing by 0.305, and
near-identical ones by up to 1.06.** If the subgroup effect were chemical they
would agree. They do not — so per-ligand ρ is dominated by sampling noise, and the
nine "independent" ligands are ~6 chemotypes, violating the exchangeability the
permutation assumed.

### 83c. Was it worth looking?

Yes, for two reasons, and the second is the more useful:

- §83a is a real finding about where Boltz is confident, and it predicts the model
  will be least useful exactly where our targets are unusual — PFAS being the
  obvious case.
- §83b is a worked example of a significance test being overturned by an **internal
  replicate**. The permutation corrected for multiplicity and still passed;
  what caught it was asking whether two versions of the same molecule agree.

> With 43 noisy per-ligand estimates, subgroup analysis will always produce a
> winner. The check that has teeth is not a p-value, it is whether near-duplicate
> inputs give the same answer.

That check is cheap here because the screening library contains salt forms and
close analogues by construction — capsaicin/zucapsaicin at Tanimoto 1.00 is a free
internal replicate, and it should be used routinely rather than only when a result
looks too good.

---

## 84. A property-matched tractability benchmark (built 2026-08-27, not yet run)

§82 showed Boltz-2 cannot rank pocket VARIANTS for a fixed ligand — a direction its
training never constrains. The ligand-level question runs **with** the training
direction: given a molecule, will PYR1 yield a sensor for it at all? sd01 supplies
the labels — **194 Hit? = Yes, 3,172 Hit? = No** — the first genuine negative set
this project has had. 181 hits and 2,531 non-hits survive SMILES parsing and
de-duplication.

### 84a. ⚠ A random negative set would be rigged, and by polarity not size

Jannis's constraint: the benchmark must not be winnable by "obviously too large for
the pocket", or a model would succeed for a reason `heavy_atoms` also delivers.
Measured, and the concern is real — though the giveaway is not the one expected:

| descriptor | AUC, **random** negatives | AUC, **matched** negatives |
|---|---|---|
| TPSA | **0.298** | 0.491 |
| H-bond donors | **0.307** | 0.499 |
| H-bond acceptors | **0.328** | 0.490 |
| rotatable bonds | **0.364** | 0.499 |
| cLogP | **0.632** | 0.517 |
| heavy atoms | 0.412 | 0.479 |
| formal charge | 0.492 | 0.500 |

With a random draw **TPSA alone reaches AUC 0.298**. Size is the weaker signal
(0.412); PYR1's hits skew greasy and aromatic, so *polarity* is what leaks.

Negatives are therefore matched to hits by greedy nearest-neighbour on z-scored
heavy atoms, cLogP, TPSA, HBD, HBA, charge and rotatable bonds, without
replacement. **Worst |AUC − 0.5| after matching: 0.021.** Median heavy atoms 20 vs
21; median cLogP 3.1 vs 3.1. The random set is kept on disk as the naive comparator,
to show how much easier the rigged version looks.

### 84b. Design decisions

- **Both classes get the same wild-type PYR1.** A hit is a ligand some *variant*
  bound, not one wild-type binds; giving hits an evolved pocket and non-hits
  wild-type would be a fatal asymmetry. The question becomes *does WT co-folding
  predict LIBRARY tractability* — a proxy one step from the label, but symmetric.
- **One MSA for all 362 runs** (`141`), since the receptor is identical throughout.
  Paying the public server 362 times would dominate cost and introduce variation in
  the one input that must be constant.
- Affinity head on, 5 diffusion samples, 5 affinity samples.

`139_tractability_set.py`, `140_tractability_yaml.py`, `141_tractability_msa.sh`,
`142_tractability_run.sh`. **362 runs, ~5 GPU-hours, ~1.5 h at 4 concurrent.**
Queued behind the umbrella hysteresis rerun by Jannis's instruction.

### 84c. What it can and cannot show

It can show whether Boltz-2 is useful **before** a campaign — is this target worth a
library — which is a real decision and one nothing in this project addresses. It
**cannot** rescue the design question: §82 settled that variant ranking is at
chance, and a ligand-level result would not change it.

---

## 85. The umbrella sampling FAILS its own pre-registered test — the arm is closed (2026-08-28)

§81 reported a failing PMF and named my own §70 repair as the prime suspect: that
repair rebuilt window 127 in a way that destroyed the two-directional seeding, so
every window was reached by pulling in one direction only. The fix was to reseed
inward from the equilibrated outermost window (`136`) and rerun (`138`), with the
acceptance criterion **fixed in advance in §81a**: a converged PMF is
seed-independent, so forward and reverse must agree to within the half-split
drift — **0.22 kcal/mol holo, 1.05 apo**.

All 38 reverse windows finished at 20 ns (job 27854551; 24 of the 62 array tasks
exited immediately because reverse seeds exist only for the 10.5–19.5 Å half of
the coordinate). `113_us_pmf.py` now implements the hysteresis check it had only
promised.

| | forward | reverse | disagreement | bar |
|---|---|---|---|---|
| holo, mean \|F−R\| after offset alignment | — | — | **3.22** | 0.22 |
| apo, mean \|F−R\| | — | — | **2.74** | 1.05 |
| holo G(16–19.5) − G(10.5–13) | **+4.61** | **−3.29** | sign flips | — |
| apo G(16–19.5) − G(10.5–13) | **+3.14** | **−2.40** | sign flips | — |

The magnitude misses by 15× and **the sign of the tilt reverses with the pull
direction**. Calibration fails on its own terms too: holo +5.02, apo +6.92,
ddG −1.90 where a pass needs positive. The pull direction is setting the answer.

**This is protocol, not physics, and it is not a sampling-length problem.** A
15× miss with a sign flip is not closed by extending 20 ns windows. §62 already
established that unbiased MD cannot reach the transition and §29 that both states
are kinetically trapped at 300 ns; this was the enhanced-sampling attempt to get
around that, and it does not work in this geometry with this coordinate.

**Do not rerun, extend, or re-window this scheme.** No open/closed ΔG from it is
quotable, and the designed-pocket switching-cost number it was built to produce
remains unavailable by simulation. That was the last route to it that did not
require new chemistry, so **switch cost stays a Y2H question** (§58).

---

## 86. The graft direction, closed — and three of my own metrics withdrawn (2026-08-28)

Jannis asked whether the pocket could be enlarged by adding residues to the β
sheet, and how the best transplantable donor's wall differs from PYR1's.
Scripts `143`–`152`, `results/pocket_shape/`.

### 86a. Enlargement is neither insertion nor displacement

| donor | cavity | wall residues | aligned to a PYR1 position | **inserted** |
|---|---|---|---|---|
| 2PCS | 570 Å³ | 45 | 43 | **2** |
| 2NS9 | 455 | 41 | 39 | **2** |
| 2BK0 | 344 | 39 | 36 | **3** |
| 6AWV | 319 | 47 | 44 | **3** |

No donor builds its pocket out of new backbone. Nor are the lobes further apart:
sheet-to-grip-helix separation at structurally equivalent positions is **14.31 Å
in PYR1 and 14.67 Å in 2PCS**, correlating with cavity volume at **r = −0.18**
across 19 relatives. Adding 2 aa per strand is the wrong lever, and in a
7-stranded sheet whose hairpins *are* the gate and latch it is also the most
expensive one.

### 86b. Three metrics of mine, withdrawn

1. **143's envelope was circular** — the C-β hull of residues the cavity itself
   selected, so a bigger cavity recruited more residues into its own hull. With
   a structurally defined wall (the same 24 aligned positions everywhere),
   cavity vs envelope falls from ~+0.99 to **r = +0.20**.
2. **143's poly-Gly column was meaningless** — 0–3 Å³ for five structures is the
   enclosed component *ceasing to exist* once buriedness drops below cut, not a
   small cavity. The leak flag missed it because a component that vanishes never
   touches the grid boundary.
3. **149's ligand-distance audit and 148's 9 Å adjacency filter were invalid in
   principle.** Jannis: *"the pocket goes beyond just ABA and selecting residues
   within a distance of ABA is not an effective way to quantify the pocket."*
   Correct — ABA is 19 heavy atoms in a cavity we are trying to enlarge, so
   scoring membership by proximity to it defines the pocket as the volume
   already occupied and excludes every expansion by construction. 148 discarded
   43 of 52 candidates on that basis. **150's `r_eq` was also wrong**: it is the
   equal-area radius of a cross-section of *grid points*, which are probe
   centres already satisfying clear > 1.4 Å, so a thin ribbon of them does not
   imply a thin passage.

### 86c. What a ligand-free measurement says

`151` decomposes each cavity by the **maximin of clearance along paths between
chambers** — the largest sphere that can be walked from one to the other. Exact
on the grid, no ligand, no straight-pocket assumption.

| structure | total | main chamber | r_max | satellites |
|---|---|---|---|---|
| 2PCS | 570 | **570** | 3.80 | none |
| 4DSB | 173 | 173 | 3.68 | none at r ≥ 2.4 |
| 2BK0 | 346 | 172 | 3.33 | 60, 103, 12 |
| 2NS9 | 462 | 157 | 3.06 | 114, 66, 127 |
| 3OQU | 204 | 146 | 3.23 | 58 @ r = 2.23 |
| 6AWV | 305 | 126 | 3.17 | 50 and 92 @ r = **1.43** |
| PYR1 | 164 | **119** | 3.21 | 19, 27 |

**A quarter of PYR1's own published volume is satellite, not chamber** — every
cross-structure volume comparison in §56/§79 was partly counting voids a ligand
cannot reach. Only 2PCS has a genuinely large single chamber; the other large
totals are fragmentation.

### 86d. Why the direction is closed anyway

I proposed 4DSB as a find on the strength of its 173 Å³ single chamber against
PYR1's 119. **Jannis rejected it, and was right:** on the cavity-independent
measures 4DSB has *less* room than PYR1 — envelope **1804 vs 1873**, poly-Gly
ceiling **1079 vs 1184** — its variable positions sit where PYR1's already are,
so it offers no new ligand chemistry, and "no satellites at r = 2.4" is not "no
bottleneck": the method only sees constrictions that *separate* chambers and is
blind to a waist within one. One metric said +45 %, two said negative, and I led
with the flattering one.

**The graft direction is closed.** It was already lower priority (§56); the
categorical tension is now explained rather than merely observed — the donors
with real chambers are at ~10 % identity and 3.5–4 Å core RMSD, and the ones
close enough to transplant have no more usable room than PYR1. Effort returns to
the non-grafting work: the library-enrichment reframing (§47c), the protein-only
ΔΔG filter (§69), and the tractability benchmark (§84), which is still built and
unrun.

### 86e. The direct test agrees, and refutes my own F108 claim

`152` ran the test 148's withdrawn filter was standing in for: truncate each
position and measure the change in **main-chamber** volume (§86c definition), with
PYR1's own 24 wall positions as a paired control. WT baseline: main chamber
122 Å³, total enclosed 164, r_max 3.21 Å.

| group | n | median ΔV_main | best | n > +20 Å³ |
|---|---|---|---|---|
| PYR1 wall | 24 | 0 | **+44** (Y120) | 3 |
| **borrowed (146)** | 50 | **0** | **+0** | **0** |

**Not one of the 50 borrowed positions opens any main-chamber volume.** The idea
is dead on the direct measurement, independently of Jannis's judgement call —
which it corroborates. The wall positions that do work are Y120 (+44), V163
(+34) and L117 (+22), all long known.

**And my §87 claim that F108 is the mouth of the lobe is refuted.** Scanning
F108 → I/L/V/A/G:

| variant | main chamber | ΔV_main | total enclosed |
|---|---|---|---|
| WT (F) | 122 | — | 164 |
| F108I | 125 | +3 | 172 |
| F108L | 125 | +4 | 173 |
| F108V | 124 | +3 | **242** |
| F108A | 124 | +3 | **242** |
| F108G | 123 | +1 | **248** |

Deleting F108 entirely adds **+84 Å³ of enclosed volume and +1 Å³ of chamber**.
The volume it releases stays behind a neck narrower than 2.4 Å, so no ligand atom
reaches it. §23's "F108 gatekeeps the second lobe" survives intact — this is what
gatekeeping looks like measured properly — but the stronger reading I offered in
§87, that 3OQU's isoleucine at this position opens a usable chamber, does not.
This is the same lesson as §58: **the second lobe is narrow, not merely blocked.**

---

## 87. The Beltran-45 benchmark: the first prospective result (2026-08-28)

§47e listed this as "the single highest-value thing not yet done", and §47f put
it before any more GPU time. `scripts/153_library_recall.py`.

**The task, and why it is the right one.** A library is a MENU — positions, each
with allowed residues, wild type always retained, size ∏(1 + nᵢ) per §53. A
sensor is CAPTURED when all of its substitutions lie inside the menu. The score
is recall at fixed size, never precision, because the ground truth is
positive-unlabeled (§51).

**Why Beltran-45 is a real test set.** Nothing in this project was built on it.
It uses 20 design positions and only 10 overlap Tian's 18; only 7 have ever been
scored here. A method tuned on Tian's vocabulary is being asked to work where it
has never looked — the transfer failure mode §47e point 3 warns about.

### 87a. Tian's frequency prior transfers

| library size | frequency (Tian sd03) | greedy-cover (sees the answers) | random null | p |
|---|---|---|---|---|
| 10³ | **33 %** | 51 % | 4 % ± 6 % | 0.001 |
| 10⁴ | **42 %** | 53 % | 7 % ± 9 % | 0.001 |
| 10⁵ | **64 %** | 76 % | 15 % ± 13 % | <0.001 |
| 10⁶ | **67 %** | 80 % | 27 % ± 17 % | 0.013 |

Substitution frequencies learned from Tian's coumarin/PFAS/etc. round-1 clones
recover **64 % of 45 real cannabinoid sensors at a library of 10⁵**, against a
random null of 15 %. This is the first result in the project that is
prospective rather than retrospective recovery of a known answer, and the first
that transfers across ligand class and across library design.

### 87b. And it saturates exactly at the vocabulary ceiling

Tian's **entire** 144-substitution vocabulary, taken as one menu, captures
**30/45 = 67 %**. The frequency ranking reaches 67 % at 10⁶ and stops — so
within Tian's vocabulary the frequency prior is already essentially optimal, and
**all remaining headroom is in POSITIONS, not identities.**

| ligand | captured by Tian's whole vocabulary |
|---|---|
| JWH-015 | 7/7 |
| WIN 55,212 | 6/9 |
| JWH-072 | 5/6 |
| **CBDA** | **0/3** |
| **∆9-THC** | **0/2** |
| **4F-MDMB-BUTINACA** | **0/1** |

Three ligands are **completely unreachable** no matter how the 144 substitutions
are re-ranked. This is the same conclusion the PFAS arm reached by a different
route (§71/§73): the headroom is positional.

### 87c. What is running

That makes the open question specific: can the one scorer that beats chance
(§69, protein-only Cartesian FastRelax ΔΔG) nominate substitutions at positions
where **no frequency prior exists at all**? `154` scores all 380 single
substitutions at Beltran's 20 positions — 10 of which this project has never
touched — by 124's method unchanged, with the paired same-shell wild-type
control that §66a's first version lacked. 20-task CPU array on `cutlerlab`,
gated behind an identity-assertion smoke test.

⚠ It remains a **stability** filter (§69, §77). It can say a substitution is
tolerated; it cannot say it makes a sensor. Its only job here is to supply a
ranking where frequency is silent.

### 87d. The false-negative side of the ΔΔG filter — and it is the side that matters

Jannis asked the question every previous report of this filter avoided: what is
the rate at which ΔΔG is *right* to eliminate something, and are there real
sensors carrying substitutions it would have thrown away?

§69 and §82d quoted **enrichment**. That is the wrong side of the ledger for
library design — a library you never build cannot be rescued by the fact that
what you did build was enriched. `155` reports recall first, with thresholds set
**within n_sub strata** so variant size cannot leak in.

| sensor retention | library kept | wild kept | **shrinkage** | enrichment |
|---|---|---|---|---|
| **100 %** | 51 % | 25 % | **1.9×** | 1.95× |
| 89 % | 49 % | 20 % | 2.0× | 1.80× |
| 76 % | 44 % | 15 % | 2.3× | 1.71× |
| 69 % | 40 % | 12 % | 2.5× | 1.73× |
| 53 % | 29 % | 7 % | 3.4× | 1.82× |

**To lose no real sensor at all, the library can only be halved.** Pushing to
3.4× costs 47 % of the sensors. And enrichment is **flat at 1.7–1.95× across the
entire range** — there is no operating point where this filter is especially
good. It is a uniform weak tilt, not a discriminator.

Set against §87a, where a substitution-frequency prior recovers 64 % of a
*held-out* sensor set at a library of 10⁵, a cost-free 1.9× is not a library
design method. It remains what §79 said it was: a clash detector, useful as a
final viability screen on an already-narrow set, not as the thing that does the
narrowing.

**Attribution is underpowered and one hypothesis of mine is refuted.** V81Y is
the most over-represented substitution among rejected sensors (+0.042 excess at
90 % retention, +0.055 at 80 %), and V164M and F159V appear only in eliminated
sensors — but on counts of 2–7, which is a lead and not a finding. I proposed
that ΔΔG, being a clash detector, should systematically reject sensors that GROW
side chains into the pocket, which is exactly what expansion needs. **It does
not: 69 of 70 real sensors net-grow**, so there is no contrast group, and the
partial correlation controlling n_sub is +0.19.

Whether a specific *residue* would be eliminated needs single-substitution ΔΔG,
which job 27925009 is computing for Beltran's 20 positions (§87c).

---

## 88. The mandipropamid benchmark was wrong, and the tuned score does not survive (2026-08-31)

### 88a. K59R was never a prediction target

Jannis read the mandipropamid paper's methods. The library was built **in the
ABA non-responsive PYR1(K59R) backbone**: K59R was isolated separately in
error-prone PCR screens against structurally dissimilar agrochemicals, and was
**forced into every library member**, not selected. Their library is stated
exactly — 475 variants, site-saturation at **25** pocket-lining residues (P55
F61 I62 V81 V83 L87 P88 A89 S92 E94 E141 F108 I110 H115 R116 L117 Y120 S122
M158 F159 A160 T162 V163 V164 N167). 25 × 19 = 475 on the nose, and 59 is not
among them.

**So §23j, §33, §36 and §47b were scoring a prediction nobody had to make.**
"K59R is invisible to every method" is true and irrelevant: the recorded verdict
*stage 1 does not clear* rested on it, and has to be revisited. `156` rebuilds
the real benchmark — all 475, scored in **both** backgrounds, asking whether
forcing K59R rescues V81I / F108A / F159L. Job 27977098.

Jannis's second correction: sensors were also screened for **cross-reactivity**,
so the best-affinity variant need not win. A low rank is evidence against the
score; a high rank is not proof of selection.

### 88b. The noise floor, measured for the first time

380 single substitutions at Beltran's 20 positions, 3 repacks each: **median
replicate spread 0.00 REU, 90th percentile 0.21.** The signal that must be
resolved — V81I + F108A + F159L are worth ~2.7 REU once F108A's clash is
relieved — is **13× the noise**. Whatever is wrong with this score, it is not
imprecision. That removes the main reason to doubt Goal 3 is reachable.

### 88c. Per-residue false negatives, and what causes them

If ΔΔG had chosen the library from those 380:

| keep top | sensor substitutions kept | recall | chance |
|---|---|---|---|
| 10 % | 9/40 | 22 % | 10 % |
| 25 % | 17/40 | 42 % | 25 % |
| 50 % | 21/40 | **52 %** | **50 %** |

Enrichment lives only at the very head and is **gone by the median**. AUC 0.606.
The worst false negatives name the mechanism: **Y120G at rank 356/380** (used by
4F-MDMB-BUTINACA, AB-PINACA, CBDA, JWH-007 — the most widely used substitution
in the set), Y120A 308, A160G 261, **F159G 258** (Beltran's round-1 DSM hit for
WIN 55,212).

| class | n | median ΔΔG |
|---|---|---|
| SHRINK (< −20 Å³) | 179 | **+2.40** |
| similar | 65 | +0.04 |
| GROW (> +20 Å³) | 136 | +0.54 |

**ref2015 penalises shrink substitutions specifically** — removing a side chain
leaves a void it scores as destabilising — and real sensors are shrink-biased
(median ΔV −16.3 Å³ against −3.7 for the rest). The filter penalises the class
it should favour. §89b's earlier test of this on 70 combinations was
uninformative because 69 of the 70 net-grow; at single-substitution resolution
there is finally contrast.

### 88d. A label-free correction, and it works

Within one volume class the raw score is already better than pooled — AUC 0.652
shrink-only, 0.683 grow-only, against 0.606 pooled. That is the signature of a
between-class offset, not a within-class failure. So a quadratic in ΔV is fitted
to ΔΔG using **only the 340 unlabelled substitutions** and subtracted; no sensor
label enters the fit.

| score | AUC | top 10 % | top 25 % | top 50 % |
|---|---|---|---|---|
| raw ref2015 ΔΔG | 0.606 | 22 % | 42 % | 52 % |
| **volume-corrected** | **0.666** | 25 % | 45 % | **70 %** |
| chance | 0.500 | 10 % | 25 % | 50 % |

**F159G moves from rank 258 to 41.** This is the first PYR1-tuned scoring change
in the project that improves anything, and it supports Jannis's proposition that
the function has to be tailored rather than generalised.

⚠ It is partial and not yet out-of-sample. The head of the list barely improves
(22 → 25 %), Y120G stays at 275, and the fit and the evaluation share the same
380 substitutions and the same structural context. The real test is the
Park-475 mandipropamid library — different ligand, different backbone,
different lab.

### 88e. The tractability array produced nothing, and why

All 362 runs failed. Every task logged `Error: '\x00'` and exited 0 with no
affinity output. Cause: **a single trailing NUL byte, the last byte of the
468,025-byte shared a3m**, which Boltz's parser rejects. One corrupt byte in the
one input deliberately shared across all 362 runs (§84b) cost the entire array.

This is the §47b failure mode again — *a setup defect that returns a plausible
state instead of an error*. Every task "COMPLETED". The MSA is repaired (1,820
sequences, NUL-free) and `142` now **asserts** the a3m is NUL-free and refuses
to run otherwise, rather than repairing it silently: a reader that tolerates the
corruption would hide the next one.

### 88f. ⚠ The volume correction does NOT transfer — 88d is downgraded

§88d reported a volume-corrected ΔΔG lifting AUC 0.606 → 0.666 on Beltran's 380
singles and moving F159G from rank 258 to 41. Tested out-of-sample on the
Park-475 mandipropamid library it **fails**, in two different ways.

**The quadratic form breaks on scale.** Beltran's singles span ΔΔG −5 to +22;
Park's span ±1,400 because the wild-type frame clashes with the mandi pose. Fitted
there, the quadratic is driven by clash outliers and its top hits become A160W,
S122W, A89W, P88W, T162W — *the largest residue at every position*, which is
what an over-fitted volume term does by construction.

**A scale-free version survives only weakly.** Converting ΔΔG to a within-volume-
class rank (shrink / similar / grow):

| | Beltran AUC | Beltran top 25 % | mandi V81I | mandi F159L |
|---|---|---|---|---|
| raw ΔΔG | 0.606 | 42 % | 239/456 | 128/456 |
| within-volume rank | 0.628 | **35 %** | 185/456 | **162/456** |

The Beltran AUC gain shrinks from +0.060 to +0.022, top-25 % recall gets *worse*,
and on mandi the two targets move in opposite directions. **The §88d result was
largely dataset-specific and should not be carried forward as a working
correction.** The shrink-penalty it was built on (§88c) is still real and
measured; the fix for it is not.

### 88g. What the Park-475 run actually showed, and the conditional follow-up

Two defects in `156`, both visible in its own output. It repacked a 6 Å shell
around the mutated position, so each of the 475 was scored in a *different*
shell and K59R — sitting outside almost every one — was frozen and **cancelled in
the paired difference**: it changed ΔΔG for 22 of 475 variants, only at
positions 55, 108, 141, 158, 164, 167. The A/B was mostly vacuous.

And one clash swamps everything. The mandi pose comes from the quadruple-mutant
crystal, so wild-type F108 clashes with it enormously: **the top 12 of 475 are
F108S/A/C/G/P/T/D/Q/N/E/V/M at −1358 to −1484 REU**, with F108A at rank 2, while
V81I ranks 260 and F159L 146. That is not a failure to find them — it is the
score correctly reporting that nothing else matters until F108 is relieved.

So the question is conditional, and `160` asks it with a **single fixed shell**
(union of all 25 pocket neighbourhoods, identical for every variant, so K59R can
act and cross-position ranks are comparable) across three backgrounds: K59R,
K59R+F108A, and K59R+F108A+F159L. If V81I and F159L rise once the dominant clash
is relieved, Goal 3 needs a **greedy sequential protocol**, not a better score.
Job 27977268.

---

## 89. Both closed-space enumerations return, and they disagree usefully (2026-08-31)

### 89a. The WIN DSM library: real but modest narrowing

All **35,863** allowed doubles of Beltran's DSM library (sd04's `DSM-Hao` menu,
18 mapped positions) scored on 7MWN reverted to wild-type PYL2, repack only,
fixed crystal WI5 pose. Replicate spread median 0.00, 90th percentile 0.00.

| round-1 WIN hit | rank of 35,863 | percentile |
|---|---|---|
| F159G+A160I | **1,547** | top 4.3 % |
| F159S+A160L | 4,805 | top 13.4 % |
| F159S+A160V | 6,632 | top 18.5 % |
| F159G+A160V | 8,276 | top 23.1 % |
| F159T+A160L | 17,532 | top 48.9 % |

Median rank **6,632 against a null median of 17,935**, permutation **p = 0.047**.
The best-of-five is not significant on its own (p = 0.20).

**Screen size to a first hit: 1,547 variants ranked, against 7,172 expected at
random — a 4.6× narrowing.** That is a real number for Goal 1, and the first time
this project has produced one on a closed, independently-defined library. It is
also modest: capturing all five still needs 49 % of the library.

**And it finds the right position for the wrong reason.** Position enrichment in
the top 1 %: V83 3.08×, **A160 2.63×**, S92 2.29×, K59 1.68× — but **F159 1.03×,
i.e. exactly chance**. The score sees that A160 must change and is blind to F159,
even though every one of the five hits pairs them. The top 10 doubles are
A160F/I/W and V83L combinations: it wants to *fill* the pocket at 160, which
half-matches the real answer (A160V/I/L) while missing F159→small entirely.

### 89b. The conditional mandi test fails: relieving the clash does not rescue V81I

`160`, fixed 64-residue shell, three backgrounds. The shell confirms F108's
dominance — background dG_bind **+1501.07 (K59R) → +17.16 (K59R+F108A)**.

| target | K59R | +F108A | +F108A+F159L |
|---|---|---|---|
| **V81I** | 224/437 | **222/437** | **221/437** |
| V81M | 400/437 | 21/437 | 22/437 |
| V163W | 426/437 | 426/437 | 426/437 |

**V81I does not move.** The greedy/sequential hypothesis of §88g is refuted: once
the dominant clash is relieved, the score still has no information about the real
substitution. Its top picks in every background are P55Y/P55M, N167G, N167P,
H115G, V81F, V81Y — and note it ranks V81F and V81Y at 7–8 while the real V81I
sits at 222, so it prefers large aromatics where the sensor uses a modest
V→I growth.

### 89c. The measurement that constrains Goal 3

The two runs bracket a genuine protocol tension.

| shell | replicate spread (90th pct) | cross-position ranks | background effects |
|---|---|---|---|
| per-position 6 Å (`156`) | **0.21 REU** | not comparable | cancel — K59R moved 22/475 |
| fixed 64-residue (`160`) | **5.01 REU** | comparable | act correctly |

The signal that must be resolved is ~2.7 REU for three substitutions, so ~1 REU
each. **The shell that makes variants comparable has a noise tail five times the
per-substitution signal; the shell that is quiet cannot compare them.** That, not
the energy function's accuracy, is the immediate obstacle — and unlike accuracy
it is addressable, by replicate averaging or by restricting the shell to the
ligand's first contact sphere.

### 89d. Two more silent-state failures in the tractability array

The repaired run skipped all 362 tasks. The NUL-poisoned first attempt had left
362 output directories containing only `processed/manifest.json`, and `142`
resumed on **directory existence** rather than on the artefact, so every task
logged "exists, skip" and exited 0. The guard now tests for
`predictions/*/affinity*.json` and deletes a stale directory before running.

That is the fourth defect in this project to return a plausible state instead of
an error, and the second in the same array (§88e). Resubmitted as 27978561.

---

## 90. What >80 % accuracy would require, and why it is not measurable today (2026-08-31)

Jannis asked whether to step back and work out what is needed to separate
plausible binders from variants *shown not to work*, at >80 % accuracy. The
question is the right one and it has a specific, checkable blocker.

**Accuracy needs true negatives, and we have almost none.** Every metric this
project has reported — recall at fixed library size (§87), hit retention (§89b),
rank of a known sensor (§89a) — is a RECALL metric, and that was not a
stylistic choice. §51 established the ground truth is positive-unlabeled: a
variant absent from a hit list may be a failure, or may simply never have been
screened. Inventorying every source we hold:

| source | level | positives | **tested negatives** |
|---|---|---|---|
| sd07 dose-response | sensor × ligand | 78 respond | **0** — all 78 rows respond somewhere in 0.025–100 µM |
| sd09 dose-response | sensor × ligand | 242 | **3** (`-` at 100 µM) |
| sd03 / sd04 clones | variant | 692 | 0 (hit lists only) |
| Beltran-45 | variant | 45 | 0 |
| **sd01** | **ligand** | **194** | **3,172** |

**Three variant-level true negatives exist in the entire corpus.** No
classification accuracy can be computed from that, at any threshold, by any
method. This is not a scoring problem and no amount of GPU time touches it.

### 90a. What that implies

1. **Accuracy is the wrong headline metric even where negatives exist.** At
   sd01's 194/3,172 balance a classifier that answers "no" every time scores
   94 %. Any target must be stated as balanced accuracy, AUC or MCC on a
   *matched* set — which is exactly why §84 property-matched its negatives and
   showed TPSA alone reaches AUC 0.298 on a random draw.
2. **The one measurable version of the question is running.** §84's tractability
   benchmark is the only true-negative test available: 181 hits against 2,531
   property-matched non-hits, worst residual descriptor leakage |AUC − 0.5| =
   0.021. It asks *given a ligand, will PYR1 yield a sensor* — the ligand-level
   form of Jannis's question. Job 27978561.
3. **At the variant level the honest answer is that the experiment has not been
   done.** The one library that was screened EXHAUSTIVELY is Park's 475 —
   every single mutant at 25 positions in the K59R backbone, assayed against
   mandipropamid. If that screen's per-variant outcome is recoverable from the
   supplementary, it converts the 437 variants already scored in §89b into a
   real classification benchmark with ~470 true negatives. **That table is the
   single highest-value missing input in the project**, and it is a literature
   retrieval, not a computation.

### 90b. What a credible >80 % claim would need

Stated in advance so it cannot be moved afterwards:

- a **held-out** set with real tested negatives, not unlabelled non-hits;
- **balanced accuracy or MCC**, with the trivial-majority baseline quoted beside it;
- negatives **property-matched** to positives, since §84 showed an unmatched draw
  is winnable by polarity alone;
- the **screening depth** stated, because a negative from a shallow screen is a
  weaker label than one from an exhaustive one;
- **cross-reactivity** treated as a separate axis (Jannis): a variant that binds
  the target but also three others is a chemistry failure, not a binding failure,
  and the two must not be pooled into one label.

Against those criteria nothing in this project currently qualifies, and the §89a
figure that comes closest — 4.6× narrowing to a first hit, p = 0.047 — is a
recall result on a positive-only set.

---

## 91. Park's Figure 5 decoded: 7,045 true negatives, and the first accuracy number (2026-08-31)

§90 named the missing input: the per-variant outcome of Park's exhaustively
screened 475-member pocket library. Jannis supplied the figure. It is an image,
so it was decoded programmatically, not transcribed.

### 91a. The encoding, and how it was verified

Each substitution cell is a **2 × 2 block of sub-cells, one per compound** —
Jannis's reading, and it is exactly right. Measured over the whole table, each
quadrant contains one compound and only that compound, with **zero
cross-contamination**:

| quadrant | compound | n |
|---|---|---|
| top-left | benzothiadiazole | 19 |
| top-right | mandipropamid | 17 |
| bottom-left | benoxacor | 9 |
| bottom-right | fludioxonil | 24 |

That is a strong self-check: a mis-registered grid would smear colours across
quadrants. Row and column geometry was anchored to the *label glyphs* (20 row
labels at 41.42 px pitch, columns fitted at 56.94 px) after an even-spacing
assumption drifted half a row by the bottom of the table. Detection separates
cleanly — 80 calls at coverage ≥ 0.19 and the next value is 0.06 — so the
threshold is not a tuning choice.

### 91b. The label set

- **475 variants** (25 residues × 19 substitutions), all in the **K59R**
  background, each against **15 agrochemicals at 100 µM** = **7,125 tested
  interactions**, matching the caption exactly.
- **8 of the 25 residues never conferred responsiveness** and are omitted from
  the figure: **55, 61, 62, 88, 110, 115, 116, 163** → 152 variants.
- **11 of the 15 compounds produced no responder at all.**

⚠ **Arithmetic corrected (Jannis).** An earlier version of this section added
"152 × 15 = 2,280" to "475 × 11 = 5,225" and reported 7,505 free negatives —
more than the 7,125 interactions that exist, which is the tell. The two sets
overlap in absent-residue × dead-compound, 152 × 11 = 1,672. The correct
statement is a partition:

| block | count | status |
|---|---|---|
| absent variants × dead compounds | 152 × 11 = 1,672 | all negative |
| absent variants × productive compounds | 152 × 4 = 608 | all negative |
| shown variants × dead compounds | 323 × 11 = 3,553 | all negative |
| **shown variants × productive compounds** | **323 × 4 = 1,292** | **the only block the figure decides** |
| total | **7,125** | ✓ |

So **5,833 interactions (81.9 %) are known negative without reading the figure
at all**, and only 1,292 cells had to be decoded. Of those, **80 are positive**
across 54 distinct variants and 1,212 negative.

**7,045 true negatives against 80 positives (1.12 % positive)** — 5,833 from
structure plus 1,212 read from the figure. This is the
first true-negative set in the project, and it is 2,348× larger than the three
variant-level negatives §90 could find in the entire published corpus.

### 91c. ⚠ V81I is not a mandipropamid hit

The decoded table says **V81I responds to benzothiadiazole and fludioxonil, not
mandipropamid**. Yet PYR1^MANDI is K59R+V81I+F108A+F159L. So V81I entered the
mandipropamid receptor through *combinatorial* mutagenesis of a hit isolated
against a different compound.

Every earlier scan that treated V81I as a mandipropamid answer was scoring the
wrong label — including §89b, which recorded V81I at rank 224/437 as a failure.
It was never a mandipropamid single-mutant positive. The genuine mandipropamid
singles are F108 (A,C,E,G,I,L,N,Q,S,T,V), F159 (A,C,I,L,M,T,V), A89W and S122G.

### 91d. The first accuracy measurement, against real negatives

Protein-only repack ΔΔG on all 475, scored against the decoded labels:

⚠ **MANDIPROPAMID ONLY.** An earlier version of this table also scored a
"responds to any compound" label. Jannis: *"we should only ever compare hits if
we use the right ligand."* That row is withdrawn and must not be reinstated. The
ligand modelled here is mandipropamid at 30 heavy atoms; benzothiadiazole is 13,
benoxacor 16, fludioxonil 18. A benoxacor responder scored against the
mandipropamid pose is being asked about a molecule twice the size of the one it
actually binds, so its rank carries no information. Negatives remain valid — a
variant tested against mandipropamid at 100 µM that did not respond is a true
negative for mandipropamid whatever else it binds.

| task | n | positives | **AUC** | best balanced acc. |
|---|---|---|---|---|
| mandipropamid, all 475 | 475 | 20 | **0.868** | 0.810 |
| mandipropamid, F108 excluded | 456 | 9 | **0.720** | 0.754 |
| trivial-majority baseline | — | — | 0.500 | 0.500 |

### 91e. It picks the substitution, not only the position

The project's standing belief (§63, §71) is that positions are recoverable and
substitutions are not. Restricted to the right ligand, that is too pessimistic.
Asking the within-position question — *given F108, which of its 19 substitutions
respond?* — removes the position signal entirely:

| position | responders | within-position AUC |
|---|---|---|
| **F108** | 11 of 19 | **0.784** |
| **F159** | 7 of 19 | **0.631** |

Ordered by ΔΔG, responders starred:

```
F108   S* A* C* G* P  T* D  Q* N* E* V* M  R  K  L* I* H  Y  W
F159   D  E  T* I* V* W  A* N  H  S  M* Q  C* G  P  L* K  R  Y
```

F108's failures are interpretable: P and D rank 5th and 7th and do not respond —
proline breaks the strand, aspartate buries a charge — while W, Y, H rank last
and correctly do not respond. F159 is weaker but still above chance.

The 20 positives sit at only 4 positions (F108 11, F159 7, S122 1, A89 1), so the
18 residues with no mandipropamid responder contribute negatives only. The
clearest false positives are the six F108 substitutions that score in the top 20
and respond to nothing at all: **F108P, D, M, R, K, H**.

**How AUC was computed.** AUC is the probability that a randomly chosen true hit
scores better than a randomly chosen true negative — 0.5 is chance, 1.0 perfect,
and it is threshold-free. For mandipropamid that is 20 × 455 = 9,100 ordered
pairs, each scored 1 if the hit has the lower ΔΔG, 0.5 if tied, 0 otherwise, and
averaged. Verified three ways that must agree: direct pairwise count 0.8682,
Mann-Whitney U/(n₁n₂) 0.8682, rank-sum formula 0.8682.

⚠ **The mandipropamid AUC is carried by the F108 block.** Median ΔΔG is −753.94
for the 20 hits against +0.70 for the 455 negatives, because 11 of the 20 are
F108X sitting near −1,400 REU. That is why the F108-excluded row is the
load-bearing one: 0.720 on the 9 remaining positives.

⚠ Balanced accuracy is at the best threshold chosen on the same data and is an
optimistic ceiling, not a held-out estimate; **AUC is the number to quote**.
And per Jannis, responsiveness is not affinity: cross-reactivity is a separate
axis this label set does not resolve.
