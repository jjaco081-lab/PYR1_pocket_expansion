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

## 21. Current status (2026-08-10) — supersedes §7b and §8

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

## 37. MM-GBSA rescoring: the K59 flip, at last (2026-08-19)

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
