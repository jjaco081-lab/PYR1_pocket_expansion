# Methods review — every computational method used in this project

One row per method: what it computes, how we used it, what it got right,
what it got wrong, and whether it is still in play.

> **A note on voice.** Throughout this document, 'I' and 'my' refer to Claude (Anthropic), the AI model that ran these analyses and wrote this file - not to the repository owner. Jannis Jacobs set the direction, chose the experiments, supplied the literature and caught several of the errors recorded here; the implementation, the mistakes written in the first person, and this summary are the model's. See README §0 for the full per-model attribution.

## Timeline

| 2026 | what happened |
|---|---|
| **Aug 6** | Project start. Pocket characterised; F108 identified as the wall sealing the satellite lobe. |
| **Aug 6-8** | Rosetta arms 1, 2, 2b: 41 variants relaxed in open and closed states. |
| **Aug 10** | The discrimination statistic is found to be WRONG. Arm 2b flips from 'blind' to F = 9.0. |
| **Aug 10** | WT MD campaign submitted - 4 systems x 3 replicates x 300 ns. |
| **Aug 11** | W385 water re-measured across 3QN1/4WVO/8EY0: it senses ABA but mandipropamid never touches it. |
| **Aug 12** | Stage 1, LigandMPNN arm: recovers F108A and F159L, misses K59R and V81I. |
| **Aug 13** | Stage 1, Rosetta FastDesign arm: same two recovered. The NULL arm fails its own control (K59 kept 0%). |
| **Aug 13** | Favour-native sweep fails at every weight. Per-residue decomposition finds the cause: +10.1 REU desolvation. |
| **Aug 14** | Coumarin benchmark reframes the task: positions are saturated, SUBSTITUTIONS carry the ligand information. |
| **Aug 17** | PYR1 rebuilt on the full 191-residue sequence. FastRelax's packer trap found and closed. |
| **Aug 18** | 3K3K exposed as a MIXED dimer; the apo-closed cell turns out to have been simulated by accident. |
| **Aug 18** | Per-position admissibility fails; pairwise enumeration gets 3 of 4 into the top 150 of 13,457. |
| **Aug 18** | Conformation x occupancy factorial built and launched - the open+ABA cell had never existed. |
| **Aug 19** | Six papers read from their Methods. Leonard's protocol needs a known weak hit before it can place a pose. |
| **Aug 19** | Coupled moves: better sampling, K59 still retained 0% of the time. Sampling is NOT the limit. |
| **Aug 19** | MM-GBSA rescoring makes the K59 flip. Changing the solvation model did what no sampler could. |
| **Aug 28** | Umbrella sampling FAILS its pre-registered hysteresis test by 15x with a sign flip. The arm is closed. |
| **Aug 28** | Graft direction closed: enlargement is neither insertion (<=3 inserted residues in any donor) nor displacement (sheet-helix r = -0.18). |
| **Aug 28** | Beltran-45, the first held-out benchmark: Tian's frequency prior transfers at 64 % recall, then saturates at the 67 % vocabulary ceiling. |
| **Aug 31** | Jannis: the mandipropamid library was built in the FORCED K59R backbone. Every 'K59R is invisible' record was scoring a prediction nobody had to make. |
| **Aug 31** | ddG replicate noise measured: 0.21 REU against a 2.7 REU signal. The score is precise; the shrink-penalty is what makes it wrong. |
| **Aug 31** | DSM doubles enumerated in full (35,863): 4.6x narrowing to a first hit, but A160 found at 2.63x and F159 at exactly chance. |
| **Aug 31** | TRUE NEGATIVES INVENTORIED: sd07 has 0, sd09 has 3, sd01 is ligand-level only. Accuracy cannot be measured on the data that exists. |

## Glossary

For readers less familiar with the structures, programs and notation used below.


**STRUCTURES (PDB entries)**

| term | what it is | why it matters here |
|---|---|---|
| `3QN1` | PYR1 + ABA + HAB1, closed ternary complex, 1.80 A. | The reference CLOSED state and the source of the ABA pose used throughout. |
| `3K3K` | PYR1 homodimer, 2 chains. | A MIXED dimer, not an apo one: chain A is apo and OPEN, chain B is CLOSED with ABA bound. The reference OPEN state comes from chain A. |
| `4WVO` | PYR1-MANDI + mandipropamid + HAB1, 2.25 A. | The engineered mandipropamid receptor. Its 4 mutations vs 3QN1 - K59R, V81I, F108A, F159L - are the ground truth every design method is tested against. |
| `8EY0` | Orthogonalised PYR1*/HAB1* + mandipropamid, 2.40 A. | Independent second mandipropamid structure; also a rare validated true-negative protein pair. |
| `7MWN` | PYL2-WIN + WIN 55,212-2. | Cannabinoid complex in PYL2, a close PYR1 relative (88% identity over the 25 ABA-proximal positions). |
| `3K90` | PYR1, 4 chains, mixed apo and ABA-bound. | Used in the satellite-lobe survey across all deposited chains. |

**LIGANDS**

| term | what it is | why it matters here |
|---|---|---|
| `ABA / A8S` | Abscisic acid; PDB chemical component A8S. Anionic at pH 7 (net -1). | The COGNATE ligand. Because WT PYR1 evolved for it, it serves as the ligand-swap NULL: a method that prefers a mutation for ABA too is not reading the ligand. |
| `Mandipropamid / 3UZ` | A fungicide, 29 heavy atoms, neutral. | The non-cognate test ligand. It CLASHES in the WT pocket (2.78 A, 7 of 29 atoms), which is the signal the steric methods exploit. |
| `WIN 55,212-2` | An aminoalkylindole synthetic cannabinoid. | Target of Beltran's best cannabinoid sensor, PYR1-WIN. |

**PROGRAMS AND FORCE FIELDS**

| term | what it is | why it matters here |
|---|---|---|
| `Rosetta / PyRosetta` | All-atom protein modelling suite. | Used for relax, design and structure building throughout. |
| `ref2015` | Rosetta's default all-atom energy function. | Its Lazaridis-Karplus desolvation term is the specific thing that made K59R unreachable - it pays +10.1 REU to bury a charged ammonium. |
| `FastRelax / FastDesign` | Rosetta protocols: repack + minimise, with (Design) or without (Relax) allowing identity changes. | TRAP: set_movemap() restricts MINIMISATION only. Repacking needs a TaskFactory or the whole pose is repacked. |
| `Coupled moves` | Rosetta sampler (Ollikainen & Kortemme 2015) coupling sequence, side-chain, backbone and ligand-pose moves. | The published fix for fixed-backbone specificity design. Tried here; did not rescue K59R. |
| `LigandMPNN` | Neural network designing sequence conditioned on backbone plus explicit ligand atoms. | The learned-model arm of stage 1. |
| `Amber / ff19SB / ff14SB` | Molecular dynamics engine and its protein force fields. | ff19SB (with OPC water) for explicit-solvent MD; ff14SB for the implicit-solvent MM-GBSA rescoring, which is what that method was validated with. |
| `GAFF2 / AM1-BCC` | General small-molecule force field and its charge model. | Used to parameterise ABA and mandipropamid for Amber. |
| `igb=8 / GBn2` | A generalised-Born implicit solvent model. | The alternative solvation treatment that finally produced the K59 flip. Requires mbondi3 radii or sander refuses to run. |
| `MM-GBSA` | Molecular Mechanics + Generalised Born Surface Area: dG_bind = G_complex - G_receptor - G_ligand, averaged over an ensemble. | An ensemble RESCORING method, not a design method - it ranks structures someone else generated. |

**METRICS**

| term | what it is | why it matters here |
|---|---|---|
| `REU` | Rosetta Energy Unit - Rosetta's arbitrary energy scale, loosely kcal/mol. | Comparable within one score function, not across. |
| `ddG (delta-delta-G)` | The change in a free-energy difference between two states or two variants. | Here usually mutant-minus-WT within one ligand arm, so receptor and force field cancel. |
| `RMSD / RMSF` | Root-mean-square deviation (from a reference) / fluctuation (about a mean). | RMSD says where a loop IS; RMSF says how much it MOVES. The gate needs both. |
| `Jaccard index` | Overlap of two sets: |A and B| / |A or B|. 1.0 = identical, 0 = disjoint. | Used to ask whether two ligands get the same residue menu. 0.84 means largely ligand-blind. |
| `kBT` | Thermal energy at the simulation temperature, ~0.6 kcal/mol at 300 K. | A natural unit for 'is this conformer accessible'. |
| `Y2H` | Yeast two-hybrid: growth or reporter readout of a protein-protein interaction. | How Tian and Beltran measure sensor activation - the ligand-dependent PYR1-HAB1 interaction drives the reporter. |

**STATISTICS**

| term | what it is | why it matters here |
|---|---|---|
| `s2 (s-squared)` | The VARIANCE - the square of the standard deviation. s2_between is the variance of the group means; s2_within is the variance within groups. | Written s2 here because these are variances, not standard deviations. The distinction is what the corrected F statistic turns on. |
| `F = s2_between / (s2_within / n)` | The one-way ANOVA F statistic. | The denominator s2_within/n is the SQUARED STANDARD ERROR OF THE MEAN. Dividing by n IS the sqrt(n) correction, expressed in variance units - so there is only one equation here, not two. The original error was using s_within itself as the yardstick. |
| `ICC(1)` | Intraclass correlation - the fraction of total variance that is between-subject rather than within-subject. | Used to ask whether a score reproduces across random seeds. |
| `Ligand-swap null` | Running the identical protocol with the COGNATE ligand and subtracting. | The project's standing scoring rule: only f(test) - f(cognate) is evidence. It inherits the method's own biases, so it needs no independence claim. |

## Status vocabulary

Fixed so it means the same thing in every row.

| Status | Meaning |
|---|---|
| **ACTIVE** | in use now, or its result is load-bearing for what comes next |
| **OPTION** | works; kept for a specific future purpose or to couple with another method |
| **PLANNED** | built or specified, not yet run |
| **RETIRED** | tested and set aside — the row records why, so it is not re-tried by accident |

Generated by `scripts/76_methods_review.py`. Spreadsheet: `data/methods_review.xlsx`.

## ACTIVE (17)

### Cavity chamber decomposition (bottleneck radius)

**Pocket - structure** · Aug 28 · CPU local · `143, 145, 150, 151, 152`

- **How it computes:** Decomposes an enclosed cavity into chambers by the MAXIMIN of clearance along paths between them - the largest sphere that can be walked from one to the other. Exact on the grid; makes no straight-pocket assumption.
- **How we used it:** Re-measured PYR1 and 19 relatives after the total-volume comparison was shown to be counting voids a ligand cannot reach.
- **Inputs:** PYR1 3QN1 frame + 10 large-cavity relatives
- **Key result:** A QUARTER of PYR1's own published cavity volume is satellite, not chamber: 164 total but 119 main. Only 2PCS has a genuinely large single chamber (570/570, r_max 3.80); 2NS9 462/157, 6AWV 305/126 with satellites behind r = 1.43 necks, 3QRZ 212/73. The large totals are mostly fragmentation.
- **Shortcomings & failures:** Only detects constrictions that SEPARATE chambers - a waist within one chamber is invisible to it, which is what Jannis saw in PyMOL for 4DSB. Volume split between chambers is nearest-seed Voronoi, approximate at boundaries; the bottleneck radii themselves are exact.
- **What we learned:** Report chamber volume, never total enclosed volume, when comparing pockets across structures. Two earlier metrics of mine were withdrawn on the way here: a circular C-beta envelope, and a poly-Gly column that reported 0-3 A^3 when the component ceased to exist.

### Library recall on Beltran-45 (held-out benchmark)

**Design - benchmark** · Aug 28 · CPU local · `153`

- **How it computes:** Scores a LIBRARY, not a variant: a menu of positions x allowed residues, size prod(1+n_i); a sensor is captured when all its substitutions lie inside. Recall at fixed size, never precision.
- **How we used it:** First genuinely held-out test in the project - 45 cannabinoid sensors nothing here was built on, using 20 design positions of which only 10 overlap Tian's.
- **Inputs:** Beltran 45 sensors; Tian sd03 as the frequency prior
- **Key result:** Tian's substitution frequency TRANSFERS: 64 % recall at a library of 1e5 against a random null of 15 % +- 13 % (p < 0.001), 33 % at 1e3 against 4 %. But it saturates at exactly the vocabulary ceiling - Tian's whole 144-substitution vocabulary captures 30/45 = 67 % - so all remaining headroom is POSITIONAL. CBDA 0/3, THC 0/2 and 4F-MDMB 0/1 are unreachable however the 144 are re-ranked.
- **Shortcomings & failures:** Beltran's sensors came from Beltran's libraries, so a position they never varied is invisible: absence is not evidence against a position.
- **What we learned:** The win is new POSITIONS, not re-ranking the known vocabulary. Any new library method must beat the frequency baseline, which is advantaged and hard to beat.

### Single-substitution ddG at Beltran's 20 positions

**Design - scoring** · Aug 28-31 · CPU cutlerlab · `154, 157, 158`

- **How it computes:** All 380 singles by protein-only Cartesian FastRelax ddG with a paired same-shell wild-type control; 3 replicates, minimum taken.
- **How we used it:** Asked the per-RESIDUE false-negative question: if ddG chose the library, which real substitutions would be discarded?
- **Inputs:** PYR1 3QN1 frame; Beltran's 20 design positions
- **Key result:** Replicate spread measured for the first time: median 0.00, 90th pct 0.21 REU, against a ~2.7 REU signal - the score is PRECISE, so imprecision is not the problem. Recall if ddG chose the library: 22 % at top-10 %, 42 % at top-25 %, and 52 % at top-50 % which is chance. AUC 0.606. ref2015 penalises SHRINK substitutions specifically (median +2.40 vs +0.54 grow) and real sensors are shrink-biased (median dV -16.3 vs -3.7), so it penalises the class it should favour.
- **Shortcomings & failures:** Y120G, the most widely used substitution in the whole set, ranks 356/380. A volume correction fitted label-free lifted AUC to 0.666 on this set but FAILED to transfer to mandipropamid and is withdrawn.
- **What we learned:** The filter is a clash detector suitable as a final viability screen on an already-narrow set, not as the thing that does the narrowing.

### Closed-space enumeration (Park-475, DSM doubles)

**Design - benchmark** · Aug 31 · CPU cutlerlab · `156, 159, 160`

- **How it computes:** Exhaustive scoring of a real published library in its own defined search space: Park's 475 site-saturation singles in the K59R backbone, and all 35,863 allowed doubles of Beltran's DSM library.
- **How we used it:** The sharpest available retrospective test - the search space is closed and stated in the papers' methods, so ranks are interpretable.
- **Inputs:** 4WVO mandi pose in the PYR1 frame; 7MWN reverted to wild-type PYL2 with WI5
- **Key result:** DSM doubles: the five round-1 WIN hits rank 1547, 4805, 6632, 8276, 17532 of 35,863; median 6632 vs a null median 17,935, permutation p = 0.047. Screen size to a first hit 1547 vs 7172 expected - a 4.6x narrowing, the first such number on a closed independently-defined library. Park-475: F108A ranks 2/475 but the whole top 12 is F108X at -1358 to -1484 REU, V81I 260 and F159L 146.
- **Shortcomings & failures:** It finds A160 (2.63x enriched in the top 1 %) and is BLIND to F159 (1.03x, exactly chance) although all five hits pair them. The conditional test refutes the greedy hypothesis: V81I stays at 224/222/221 across K59R, +F108A and +F108A+F159L backgrounds, and the score prefers V81F and V81Y at ranks 7-8.
- **What we learned:** Jannis found the benchmark itself was wrong: the mandi library was built in the FORCED K59R backbone, so every earlier record of 'K59R is invisible to every method' was scoring a prediction nobody had to make.

### Cavity characterisation (lib_cavity)

**Pocket - structure** · Aug 6-10 · CPU local · `00-04, lib_cavity.py`

- **How it computes:** Grid-based cavity detection; volume from probe-accessible grid points, with residues assigned by cavity lining + line of sight rather than distance to ligand.
- **How we used it:** Defined the pocket, measured expansion per variant, found the length limit and the lobe gatekeeper.
- **Inputs:** ABA (A8S); 3QN1, 3K3K, +7 deposited PYR1 chains
- **Key result:** Two findings about DIFFERENT states, which must not be compressed into one sentence. In WILD TYPE the satellite lobe is a sealed void, not connected to the ABA cavity, in 7 of 7 deposited chains - and F108 is the wall that seals it. In the QUAD mutant (K59/F108/E94/R79 truncated) the two lobes have fully merged into one continuous chamber with NO constriction: cross-sectional radius 4.0-6.3 A along the main span, and the apparent 'neck' at 5.02 A is the shallow minimum of a broad plateau. So the remaining limit is ligand LENGTH, not a bottleneck - the only narrow region is the gate/latch tail (1.4-2.7 A).
- **Shortcomings & failures:** Distance-to-ligand mislabelled E94 as second shell (it is first shell, 4.33 A).
- **What we learned:** Define pocket residues by cavity lining and line of sight, not a distance cutoff. R79 is a cavity-rim residue.

### Rosetta FastRelax cavity scan (Arm 1)

**Pocket - energetics** · Aug 6-7 · CPU cutlerlab · `09, 20; job 27260956`

- **How it computes:** All-atom relax with ref2015, then dCavity vs WT across a designed variant panel.
- **How we used it:** Ranked 41 pocket variants by how much cavity they open, x3 replicates.
- **Inputs:** ABA; 3QN1 closed, 3K3K open
- **Key result:** Worked as intended: true between-variant sd 46.9 A^3, F = 139. Produced the variant portfolio.
- **Shortcomings & failures:** Expansion does NOT predict switch cost (r = 0.03), so dCavity alone cannot rank designs.
- **What we learned:** Cavity opening is measurable and reproducible, but it is a means to a portfolio, not a selection criterion.

### Corrected discrimination statistic

**Pocket - statistics** · Aug 10 · CPU local · `35_variance_decomposition.py`

- **How it computes:** F = s2_between / (s2_within / n), where s2 is the VARIANCE (standard deviation squared). The denominator s2_within/n is the squared standard error of the mean, so dividing the within-group variance by n IS the sqrt(n) correction written in variance units - the two are the same equation, not two.
- **How we used it:** Re-derived every arm's discrimination after finding the original statistic wrong.
- **Inputs:** Re-analysis of existing Rosetta output; no new structures.
- **Key result:** Flipped Arm 2b from 'blind' to F = 9.0, correcting a central project claim.
- **Shortcomings & failures:** The original statistic used s_within itself as the yardstick - it compared between-group spread to the scatter of a SINGLE replicate rather than to the scatter of the replicate MEAN, which is sqrt(n) times smaller. The error had already propagated into several conclusions before it was found.
- **What we learned:** Compare between-group spread against the SEM, not the replicate sd. Run the variance decomposition BEFORE concluding 'no effect'.

### Per-residue energy decomposition

**Pocket - diagnosis** · Aug 13 · CPU cutlerlab · `55, 55b; job 27438909`

- **How it computes:** Breaks ref2015's total into per-term, per-residue contributions at one position.
- **How we used it:** Asked WHY K59 could not be retained.
- **Inputs:** ABA anion + K59 in 3QN1
- **Key result:** THE KEY DIAGNOSIS. ref2015 correctly sees the salt bridge (Lys best fa_elec by 1.9 REU, only candidate earning hbond_sc) then pays +10.118 REU of Lazaridis-Karplus desolvation to bury the ammonium. Ile beats Lys by 4.211 REU; Arg is 2.274 WORSE than Lys.
- **Shortcomings & failures:** Ruled out pose (K59 NZ 2.85 A from carboxylate) and ionisation (params sum -0.970 over 38 atoms) first, costing extra runs.
- **What we learned:** Diagnose the mechanism before proposing a fix. This one number explained three downstream failures and told us exactly which stage to replace.

### Pairwise enumeration

**Pocket - geometry** · Aug 18 · CPU local · `72_pairwise.py`

- **How it computes:** Enumerates mutation PAIRS against a jointly re-scored pocket. Objective = summed ligand-protein vdW overlap (additive, hence a cheap precompute) subject to not losing packing.
- **How we used it:** The version of the combinatorial idea the data supports, after per-position menus failed.
- **Inputs:** Mandipropamid (clashes in WT: 2.78 A, 7/29 atoms) vs ABA (fits: 0.22 A, 0/19) - a built-in negative control
- **Key result:** Best steric result so far: F108A rank 5, F159L rank 75, V81I rank 150 of 13,457 - 3 of 4 in the top 150, a 90x narrowing. Shrink/grow enriched to 62% of the top 200 vs 42% for ABA.
- **Shortcomings & failures:** K59R at rank 3907 - blind BY CATEGORY, since K59 does not clash at all. Picks the right POSITIONS (81+108) with the wrong identity at 81.
- **What we learned:** An ADDITIVE objective cannot show synergy. Counting contacts within a flat 4.5 A rewards clashes (83W+163Y topped the Pareto front with relief -71.7 A).

### Shared steric rule (lib_sterics)

**Pocket - geometry** · Aug 19 · CPU local · `lib_sterics.py`

- **How it computes:** One clash rule for all geometric screens: exempts donor-acceptor pairs down to 2.5 A, keeps acceptor-acceptor and donor-donor as clashes, includes a real hydrogen radius.
- **How we used it:** Fixed hydrogen bonds being scored as steric clashes across scripts 71/72/73.
- **Inputs:** All pocket screens
- **Key result:** Unit-tested three ways. R79's crystal side chain goes from 0.306 A 'overlap' to 0.000. WT-ok at hard-sphere tolerance improved 12 -> 16 of 26.
- **Shortcomings & failures:** Did NOT rescue the R79/V83/H115 control failures - I had claimed it would. Those are rotamer-basin failures and the earlier diagnosis was right.
- **What we learned:** A hard-sphere test scores every hydrogen bond as a ~0.4 A clash. The fix was invisible at the 0.5 A tolerance in use, which already absorbed it - a stability check, not a revision.

### MM-GBSA rescoring

**Pocket - energetics** · Aug 19 · CPU cutlerlab, 14 tasks · `75, 75b-e; jobs 27561665, 27568930, 27569546`

- **How it computes:** Replaces ref2015's pairwise Lazaridis-Karplus solvation with generalised Born: dG_bind = G_complex - G_receptor - G_ligand, averaged over an MD ensemble. ff14SB + GAFF2/AM1-BCC, igb=8, 0.15 M salt. Rosetta builds structures and never scores them.
- **How we used it:** Direct attack on the diagnosed cause of the K59 failure.
- **Inputs:** Mandipropamid vs ABA; 7 variants x 2 arms x 100-frame ensembles
- **Key result:** THE K59 FLIP. With mandipropamid R beats Q and N (-1.53 vs +0.04, +2.53); with ABA wild-type Lys beats every substitution (R least badly, +7.06). V81I and F159L also prefer mandipropamid.
- **Shortcomings & failures:** Ensembles only 50 ps (GBn2 runs ~380 steps/min here; pmemd no faster). The ligand-swap difference is dominated by DAMAGE TO THE ABA COMPLEX, not gain for mandipropamid. F108A comes out NULL because relaxation absorbs the clash it exists to relieve.
- **What we learned:** Changing the solvation model does what no sampler could - the desolvation diagnosis was right AND actionable. Also: my first error bars were fabricated (8 'independent' repacks came out byte-identical).

### Coumarin focused-library benchmark

**Pocket - benchmark** · Aug 14 · CPU local · `62_coumarin_ground_truth.py`

- **How it computes:** Extracts the exact per-position residue menus Tian et al. built for each focused library, then measures what a method must add over a ligand-blind baseline.
- **How we used it:** Asked whether round-1 screening could have been skipped computationally.
- **Inputs:** LIBRARY DESIGNS, not structures: Tian sd04 per-position residue menus for the Coumarin, PFAS and TNTv1/v2 focused libraries, plus 692 characterised sensor clones (sequences).
- **Key result:** Reframed the task. POSITIONS are near ligand-independent: Coumarin/PFAS/TNTv2 share 9 of their 11-14 positions, so naming the shared nine scores 82% recall with ZERO ligand information. SUBSTITUTIONS are ligand-specific: mean 3-way Jaccard 0.02.
- **Shortcomings & failures:** I initially claimed the coumarin library design was missing from the supplement. It was not - a column header truncated at 20 characters.
- **What we learned:** Predict MENUS, not positions. Focusing buys ~16 orders of magnitude (3.8e21 -> 1.4e5).

### Beltran cannabinoid sensor table

**Pocket - benchmark** · Aug 18-19 · CPU local · `data/beltran/win_sensors.json`

- **How it computes:** 45 experimentally validated sensors across 12 cannabinoids, with mutation sets and a graded dose-response (1e4 -> 10 nM).
- **How we used it:** A second, harder ground truth with an ORDER to recover, not just a set.
- **Inputs:** Sensor SEQUENCE table (no solved structures): 45 PYR1 variants vs WIN 55,212-2, the JWH series, D9-THC, CBDA, CP 47,497. WIN 55,212-2 is an aminoalkylindole synthetic cannabinoid, so it belongs in this set despite the name looking unlike the others.
- **Key result:** Extracted in full. A160 mutated in 38/45, Y120 in 33/45, F159 in 24/45. Ligand-specific signal sits at K59 (Q/N/S/T/A/R), H115Q, V83, E141, V164, V81.
- **Shortcomings & failures:** Even more position-saturated than coumarin; substitutions reused within a chemical family (A160G 27x, Y120G 23x). Nine substitutions sit OUTSIDE the 5 A pocket - PYR1-WIN gets its last order of magnitude from E4G.
- **What we learned:** A graded sensitivity ladder is a far stronger benchmark than a hit list. K59Q/K59N are neutral, so they should be reachable where K59R was not - a free sanity check.

### Literature protocol review

**Pocket - benchmark** · Aug 19 · local · `data/papers/`

- **How it computes:** Read six papers from their Methods sections, not abstracts, to establish what the field does and what the honest baseline is.
- **How we used it:** Calibrated our results against published practice before investing further.
- **Inputs:** Six papers (PDFs + extracted Methods text): Leonard 2026, An 2024, Tian 2025, Park 2023, ACS Chem Biol 2024, Nat Biotech 2026.
- **Key result:** Leonard's 'dock to sequence' MUTATES WT to a known weak-hit sequence before placing any pose - their method improves a known weak binder and cannot start from a ligand alone. An et al. picked cutoffs by MANUAL INSPECTION. RbsB precedent: 2M variants screened for 1.2-1.5x induction.
- **Shortcomings & failures:** Leonard hard-filters on an H-bond to the latch water, which mandipropamid misses by 5.07-5.19 A - their protocol would reject the correct pose.
- **What we learned:** The field's baseline is low, and our benchmark is strictly harder than Leonard's. The ACS Chem Biol paper independently confirmed R79 H-bonds F52's carbonyl at >90% occupancy, explaining one of our control failures from the other direction.

### WT MD baseline / loop dynamics

**MD** · Aug 10 - Aug 18 · GPU · `58-60; jobs 27333712, 27545237`

- **How it computes:** Explicit-solvent MD (ff19SB/OPC, 0.15 M KCl, 2 fs, 300 ns) with a two-reference projection: fit on the rigid core, then measure gate/latch RMSD to BOTH open and closed crystals.
- **How we used it:** Asked whether open and closed states are separable and whether MD can filter designs.
- **Inputs:** S1 apo-open (3K3K A), S2 holo-closed (3QN1 A + ABA), S3 dimer, S4 ternary (+HAB1, Mn2+); 12 trajectories
- **Key result:** Neither state ever converts in 1.8 us aggregate. Replicate-mean gap 5.57 A at the gate, stable at every discard from 0-150 ns. RMSF ladder: gate 3.05 open -> 1.39 closed -> 0.81 ternary.
- **Shortcomings & failures:** Licenses a STABILITY filter only (gate RMSD-to-closed, 3.5-4.0 A) and is BLIND to switchability - the likelier failure for an enlarged pocket. The latch fails as an RMSF filter (ratio inverts).
- **What we learned:** Never read a stable closed trajectory as a working switch. Inference is at replicate level: autocorrelation times 3.4-44.9 ns, so each 300 ns run holds only 3-35 independent samples.

### Apo-closed protomer analysis

**MD** · Aug 18 · GPU · `58-60; job 27545237`

- **How it computes:** Re-analysis of the S3 mixed dimer, whose chain B is a closed, ligand-shaped protomer simulated with the ABA stripped.
- **How we used it:** Filled the missing cell of the factorial using data already on disk.
- **Inputs:** 3K3K chain B (closed, was ABA-bound), 3 x 300 ns
- **Key result:** Apo-closed does NOT open: zero crossings, drift -0.07 A, and the MOST RIGID gate of all 15 units (RMSF 0.72 A vs 1.39 with ABA).
- **Shortcomings & failures:** Confounded by the dimer, which alone drops the open gate 3.05 -> 1.07 A - larger than the ligand effect being measured.
- **What we learned:** The stability filter is CONFORMATION-reporting, not ligand-reporting: a good gate-RMSD score carries no evidence about occupancy. Also: 3K3K is a MIXED dimer, and the name S3_apo_dimer had outranked four records saying so.

### 191-residue rebuild pipeline

**MD - structure** · Aug 17 · CPU local · `67, 67c, 68; lib_resnumber/structqc/rosetta`

- **How it computes:** Biopython graft + PyRosetta rebuild of the complete 191-aa sequence with no chain breaks, seeded and byte-reproducible, with identity and provenance assertions at every step.
- **How we used it:** Removed the gapped-receptor confound from every future MD system.
- **Inputs:** 3K3K A/B, 3QN1 A; ABA
- **Key result:** 7 verified structures; every chain 191 residues, 0 peptide-bond breaks, 0 clashes, ligand-lining side chains unmoved to 0.000 A.
- **Shortcomings & failures:** Cost two structural defects before the cause was found - a 0.36 A inter-chain clash and K59 collapsing 2.85 -> 1.68 A into an empty cavity.
- **What we learned:** FastRelax.set_movemap() restricts MINIMISATION only; repacking needs a TaskFactory or the whole pose is repacked. Backbone RMSD reported 0.000 A throughout - always check side chains and inter-chain contacts.

## ACTIVE (running) (1)

### Conformation x occupancy factorial

**MD** · Aug 18 - running · GPU · `69, 69b, 70; job 27547457` · **priority: now**

- **How it computes:** The 2x2 that README 19b set up and only ever filled on the diagonal: open/closed x apo/holo, all four cells on ONE build protocol.
- **How we used it:** Separates 'the gate holds its shape' from 'the ligand holds the gate'.
- **Inputs:** S1 open+apo, S2 closed+ABA, S9 closed+apo, S10 open+ABA (ABA transplanted, backbone bit-identical)
- **Key result:** S10 (open + ABA) never existed before and is now built and verified. 5 of 12 replicates complete, including both S10 reps 0 and 2.
- **Shortcomings & failures:** Not yet analysed. The old-tree S9 could not be reused: it is KCl where S1/S2 are NaCl, which would put a cation swap inside the key comparison.
- **What we learned:** Gate closure on ligand binding is the one arm where the barrier is plausibly downhill - pre-registered so 'nothing happened' cannot later be written up as 'ABA does not close the gate'.

## OPTION (7)

### Foldseek + CATH50 superfamily scan

**Pocket - structure** · Aug 6-7 · CPU cutlerlab · `03, 05, 08, 10` · **priority: later**

- **How it computes:** Structural (not sequence) homology search over CATH50; hits re-fetched full-atom and their cavities measured.
- **How we used it:** Asked whether an SRPBCC relative already has a bigger pocket to borrow from.
- **Inputs:** 266 structural hits; 147 full-atom cavities
- **Key result:** Found CoxG (2pcsA00 / 8UDS): ~570-600 A^3 pockets vs PYR1's 174, evolved to extract C45 menaquinone from membranes.
- **Shortcomings & failures:** CoxG's pocket is OPEN, which conflicts with PYR1's closed-state readout. Top-hit calls unstable to --max-seqs.
- **What we learned:** A large natural pocket in the same fold exists, but pocket size and allosteric readout are coupled - you cannot import one without the other.

### Rosetta open-vs-closed ddG (Arm 2 / 2b)

**Pocket - energetics** · Aug 7-10 · CPU cutlerlab · `22, 26, 32, 33; jobs 27275686, 27286483/4` · **priority: later**

- **How it computes:** Relax each variant in open and closed states; the ddG between them proxies the energetic cost of switching.
- **How we used it:** Tried to score whether an enlarged pocket still prefers the closed, HAB1-competent state.
- **Inputs:** ABA; 3QN1 vs 3K3K; 41 variants at nstruct 3, 20, 12
- **Key result:** Arm 2b reached F = 9.0 with a true spread of 1.79 REU.
- **Shortcomings & failures:** Arm 2 v1 was declared statistically blind - wrongly. The readout is INVARIANT, not blind: true spread 0.43-0.50 REU against a -67 REU interface.
- **What we learned:** A null is only a null if the assay can discriminate. Judge effects in physical units, not significance.

### Arm 3 homolog cavity survey

**Pocket - structure** · Aug 6-7 · CPU cutlerlab · `08, 10` · **priority: later**

- **How it computes:** Full-atom cavity measurement across foldseek homologs.
- **How we used it:** Asked how PYR1's pocket compares to its structural family.
- **Inputs:** 147 full-atom homolog structures
- **Key result:** Placed PYR1's 174 A^3 pocket in family context; supported the CoxG lead.
- **Shortcomings & failures:** First pass returned 0.0 A^3 for every cavity - foldseek stores Ca only, so convert2pdb emitted backbone traces. Silent and plausible-looking.
- **What we learned:** Validate a downloaded structure with the consuming tool before trusting a derived number. A zero can be a parsing artefact.

### LigandMPNN

**Pocket - design** · Aug 12 · GPU · `43, 46; job 27412775` · **priority: later**

- **How it computes:** Graph neural network designing sequence conditioned on backbone AND explicit ligand atoms; outputs per-position amino-acid probabilities.
- **How we used it:** Stage-1 retrospective test: can it recover the 4 known mandipropamid mutations from WT plus a perfectly placed ligand?
- **Inputs:** Structures: WT PYR1 with mandipropamid (ligand 3UZ, from 4WVO) vs ABA (ligand A8S, from 3QN1) as the ligand-swap null.
- **Key result:** Recovered F108A (+0.111) and F159L (+0.071) with zero mass in the ABA arm - genuinely ligand-conditional.
- **Shortcomings & failures:** Missed K59R and V81I. V81I INVERTED (-0.841): ABA prefers Ile at 81 with p=0.996 in the arm whose answer is no mutations. Recovers only severe-clash positions, which a trivial clash baseline already gives free.
- **What we learned:** Sequence recovery is the wrong metric - a ligand-blind oracle reaches 79.1% on the same labels. Score probability MASS against a ligand-swap null, never recall-at-N.

### Rosetta FastDesign

**Pocket - design** · Aug 13 · CPU cutlerlab · `50, 51; job 27421344` · **priority: later**

- **How it computes:** Monte Carlo sequence design with ref2015: repack and minimise while allowing identity changes at chosen positions.
- **How we used it:** The physics arm of the same stage-1 test, chosen because it fails differently from a learned model.
- **Inputs:** Mandipropamid vs ABA, 15 designable positions
- **Key result:** Recovered F108A at frequency 1.000 vs 0.000 in the null - the cleanest single recovery of the stage.
- **Shortcomings & failures:** F159L: right position, wrong identity (picks D/E). V81I inverted (-0.980). The NULL ARM FAILED ITS OWN CONTROL: with the cognate ligand it kept WT at 29% and K59 in 0% of trajectories.
- **What we learned:** If the null also mutates a position, 'no signal' is inseparable from 'the protocol cannot hold a native contact'. K59R and V81I were uninterpretable, not negative.

### Coupled moves (Kortemme)

**Pocket - design** · Aug 19 · CPU cutlerlab, 100 tasks · `74, 74b; job 27560390` · **priority: later**

- **How it computes:** Monte Carlo coupling sequence, side-chain, BACKBONE and ligand-pose moves - the published fix for fixed-backbone specificity design (5.75x on its own benchmark, 16/17 cases).
- **How we used it:** Re-ran stage 1 changing ONLY the sampler: same positions, same arms, same scoring rule.
- **Inputs:** Mandipropamid vs ABA; 50 trajectories/arm x 1000 trials, ligand_mode on
- **Key result:** F108A recovered cleanly (+0.74 against a 0.00 null). WT retention improved 29% -> 42%.
- **Shortcomings & failures:** K59R still missed. Decisively: in the ABA arm - cognate ligand, K59 unambiguously correct - K59 was retained in 0% of 50 trajectories, identical to fixed backbone.
- **What we learned:** SAMPLING IS NOT THE LIMIT at the design stage. Lysine is in every rotamer library and IncludeCurrent was on. Confirms the scoring diagnosis with a second independent sampler.

### Electrostatic / H-bond complementarity

**Pocket - geometry** · Aug 18 · CPU local · `73_electrostatic.py` · **priority: later**

- **How it computes:** Counts geometrically valid donor-acceptor pairs to the ligand plus a screened-Coulomb term. Deliberately models NO desolvation.
- **How we used it:** The stage pairwise enumeration specified: score the non-clashing positions electrostatics can reach.
- **Inputs:** ABA (anion, -0.97) vs mandipropamid (neutral, +0.10)
- **Key result:** Gets K59 right for ABA - Lys rank 1 via a bidentate salt bridge (NZ-O4 2.83 A, NZ-O3 3.31 A) where stage 1 kept it 0% of the time. The ref2015 pathology is genuinely absent.
- **Shortcomings & failures:** Missed K59R. Fails its own controls in the PREDICTED direction: WT top choice at only 4% of positions, Lys/Arg picked at 59% against an anionic ligand - with no desolvation, everything wants to be charged.
- **What we learned:** Not usable as a ranker; usable as a FILTER over supplied conformations. Its own miss was sampling: the crystal Arg59 (chi3 ~100 deg, non-rotameric) scores 2 H-bonds under its own criteria but is never proposed.

## PLANNED (2)

### Non-cognate ligand MD (S6-S8)

**MD** · Aug 14 (built) · GPU · `63-66` · **priority: P2 - after the factorial**

- **How it computes:** Same closed receptor, three ABA-sized ligands from three chemical classes, three docked poses each.
- **How we used it:** Designed to test whether the stability filter can rank LIGANDS or only conformations.
- **Inputs:** Imperatorin, Flutamide, alpha-Estradiol (19-20 heavy atoms, MW 270-276) in WT closed PYR1
- **Key result:** Built and pose-verified; three replicates per ligand start from three DIFFERENT poses, so pose sensitivity is measured rather than assumed.
- **Shortcomings & failures:** Never submitted. Partly superseded by the apo-closed result, which answered the same question more cheaply.
- **What we learned:** Absence of ligand release is NOT evidence of binding - residence times are us-ms and 300 ns cannot sample unbinding.

### ABA enantiomer / non-cognate controls

**MD** · planned · GPU · `-` · **priority: P2 - after the factorial**

- **How it computes:** Isosteric control: (-)-ABA is identical in formula, bulk and parameters, differing only in chirality, so any difference is chemistry rather than size.
- **How we used it:** Requested, to test whether gate behaviour is ligand-specific.
- **Inputs:** (-)-ABA plus an unrelated matched molecule; WT PYR1
- **Key result:** Not yet run.
- **Shortcomings & failures:** Must NOT be read as a gate-opening test: nothing opens on this timescale in any of 15 units. The graded RMSF ladder is the observable that responds to ligand.
- **What we learned:** Phaseic acid would add a third rung, turning a binary contrast into a graded ladder that is much harder to explain away.

## RETIRED (2)

### Favour-native weight sweep

**Pocket - calibration** · Aug 13 · CPU cutlerlab · `53, 54; job 27438220`

- **How it computes:** Adds a bonus for retaining wild-type identity; swept 0.25/0.5/1.0/1.5 on the NULL arm only, criterion fixed before submission.
- **How we used it:** Attempt to repair FastDesign's WT-retention failure.
- **Inputs:** ABA null arm
- **Key result:** WT recovery climbed 29 -> 45 -> 59 -> 65 -> 85%.
- **Shortcomings & failures:** K59 stayed at 0% AT EVERY WEIGHT, including 1.5 where the packer is nearly frozen (3 distinct sequences in 10). No weight passed.
- **What we learned:** The calibration failing IS the answer, not a setback. Do not widen a sweep until a control agrees - that is tuning until it says what you want.

### Per-position steric admissibility

**Pocket - geometry** · Aug 18 · CPU local · `71_admissibility.py`

- **How it computes:** For each position x each of 20 residues, place every library rotamer and ask whether ANY is clash-free against ligand and environment. Hard yes/no on Bondi radii, no score function.
- **How we used it:** Headroom check: does pure geometry narrow the residue choice per position?
- **Inputs:** ABA vs mandipropamid, 26 pocket positions, WT backbone
- **Key result:** Clean negative in an afternoon: mean 14.7/20 admissible, ABA-mandi Jaccard 0.84.
- **Shortcomings & failures:** DOES NOT CLEAR. Too permissive (9x narrowing where Tian's focused library is 50x smaller) and largely ligand-blind. All four ground-truth substitutions are admissible for BOTH ligands - recalled without being ligand-conditional.
- **What we learned:** YES, IT EXCLUDED REAL SENSOR MUTATIONS - but read the scope carefully. Against mandipropamid's OWN ground truth it excluded nothing: all four of K59R/V81I/F108A/F159L were admissible (the problem there was that they were admissible for ABA too). Against Beltran's cannabinoid sensors it wrongly excluded six real substitutions - 7% of in-pocket occurrences - and NONE of them come from a WIN sensor: V83W (JWH-007 x2, JWH-167), V83F (JWH-167), V164W (CBDA x2), V164H (JWH-015), A89V (JWH-193), E141F (CP 47,497). All NINE WIN sensors passed, and that is informative rather than lucky: every WIN mutation is F159->A/G/S/T (a shrink, always admissible), A160->I/L/V (a modest grow at a small residue), K59Q, or a surface position - exactly the substitution types a per-position filter can handle. The failures concentrate in the bulky aromatic GROWS (Trp, Phe, His at V83/V164) that need a partner shrink first. Caveat: all of these were scored in the ABA/mandipropamid pocket, not a cannabinoid-bound one, so it is not a clean test of them. Every miss was a GROW mutation, and real sensors are compensating shrink/grow PAIRS (F159A+A160I; Y120G+A160G in 27 and 23 of Beltran's 45 sensors) that a per-position filter structurally cannot see.
