# PYR1 Binding-Pocket Expansion

Enlarging the ligand-binding cavity of *Arabidopsis thaliana* PYR1 so the
PYR1–HAB1 chemically-induced-dimerisation biosensor can accept larger and more
chemically diverse ligands, **without** altering the gate and latch loops that
drive the conformational switch.

Project started 2026-08-06. Working directory
`/bigdata/cutlerlab/jjaco081/PYR1_pocket_expansion`.

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
*No MD output yet — jobs have not started writing. Re-run `scripts/40_md_provenance.py` once they do.*
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
| **27333712** | WT MD, 12 tasks `%4`, 4 systems × 3 replicates × 300 ns | queued/running, ~4 days |

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

### Bugs caught before they cost anything

| date | bug | how it surfaced |
|---|---|---|
| 08-06 | `/scratch` is node-local → SLURM job died at 00:00:00 with no log | job failed instantly |
| 08-06 | foldseek OOM-killed by the 1 GB login-node cap | bare `Killed` message |
| 08-06 | `VRT` residue has no CA → FastRelax RMSD crash | filtered with `is_protein()` |
| 08-07 | chained `sed` renamed the path before the filename rule matched, so `26_submit_ratchet_v3.sh` called a **nonexistent** `21_rosetta_ratchet_v3.py` | caught while writing `32_`; queued job cancelled and resubmitted **before it ran** |
| 08-10 | `loadamberparams frcmod.ions234lm_126_opc` — TIP3P-era naming, no such file for OPC | tleap exited on all four MD systems |

### Standing methodological lessons

1. **Judge an effect by its size in physical units, not by significance.**
   `latch_bb_rmsd` has F = 8.9 and a true spread of **0.027 Å**.
2. **Compare between-group spread against the SEM, not the replicate sd** (§14a).
3. **A null result is only a result if the assay can discriminate** — always run
   the variance decomposition before concluding "no effect".
4. **Static modelling could not settle the satellite lobe.** Seven crystals
   versus one relax protocol is why the MD baseline exists (§19a).
