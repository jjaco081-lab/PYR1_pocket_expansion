#!/usr/bin/env python
"""
76_methods_review.py -- one row per computational method used in this project.

Emits, from a single source of truth:
  * data/methods_review.xlsx   -- wide table, frozen header, colour-coded status
  * METHODS_REVIEW.md          -- the same content for GitHub

Status vocabulary, fixed so it means the same thing in every row:
  ACTIVE   -- in use now, or its result is load-bearing for what comes next
  OPTION   -- works, kept for a specific future purpose or to couple with another method
  PLANNED  -- built or specified, not yet run
  RETIRED  -- tested and set aside; the row records why, so it is not re-tried by accident
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COLS = ["Method", "Workstream", "What it computes (how)", "How we used it",
        "Inputs (structures / ligands / libraries)", "Dates (2026)", "Key result",
        "Shortcomings & failures", "What we learned", "Status", "Priority",
        "Where it ran", "Scripts / jobs"]

ROWS = [
["Cavity characterisation (lib_cavity)", "Pocket - structure",
 "Grid-based cavity detection; volume from probe-accessible grid points, with residues assigned by cavity lining + line of sight rather than distance to ligand.",
 "Defined the pocket, measured expansion per variant, found the length limit and the lobe gatekeeper.",
 "ABA (A8S); 3QN1, 3K3K, +7 deposited PYR1 chains", "Aug 6-10",
 "Two findings about DIFFERENT states, which must not be compressed into one sentence. In WILD TYPE the satellite lobe is a sealed void, not connected to the ABA cavity, in 7 of 7 deposited chains - and F108 is the wall that seals it. In the QUAD mutant (K59/F108/E94/R79 truncated) the two lobes have fully merged into one continuous chamber with NO constriction: cross-sectional radius 4.0-6.3 A along the main span, and the apparent 'neck' at 5.02 A is the shallow minimum of a broad plateau. So the remaining limit is ligand LENGTH, not a bottleneck - the only narrow region is the gate/latch tail (1.4-2.7 A).",
 "Distance-to-ligand mislabelled E94 as second shell (it is first shell, 4.33 A).",
 "Define pocket residues by cavity lining and line of sight, not a distance cutoff. R79 is a cavity-rim residue.",
 "ACTIVE", "-", "CPU local", "00-04, lib_cavity.py"],

["Foldseek + CATH50 superfamily scan", "Pocket - structure",
 "Structural (not sequence) homology search over CATH50; hits re-fetched full-atom and their cavities measured.",
 "Asked whether an SRPBCC relative already has a bigger pocket to borrow from.",
 "266 structural hits; 147 full-atom cavities", "Aug 6-7",
 "Found CoxG (2pcsA00 / 8UDS): ~570-600 A^3 pockets vs PYR1's 174, evolved to extract C45 menaquinone from membranes.",
 "CoxG's pocket is OPEN, which conflicts with PYR1's closed-state readout. Top-hit calls unstable to --max-seqs.",
 "A large natural pocket in the same fold exists, but pocket size and allosteric readout are coupled - you cannot import one without the other.",
 "OPTION", "later", "CPU cutlerlab", "03, 05, 08, 10"],

["Rosetta FastRelax cavity scan (Arm 1)", "Pocket - energetics",
 "All-atom relax with ref2015, then dCavity vs WT across a designed variant panel.",
 "Ranked 41 pocket variants by how much cavity they open, x3 replicates.",
 "ABA; 3QN1 closed, 3K3K open", "Aug 6-7",
 "Worked as intended: true between-variant sd 46.9 A^3, F = 139. Produced the variant portfolio.",
 "Expansion does NOT predict switch cost (r = 0.03), so dCavity alone cannot rank designs.",
 "Cavity opening is measurable and reproducible, but it is a means to a portfolio, not a selection criterion.",
 "ACTIVE", "-", "CPU cutlerlab", "09, 20; job 27260956"],

["Rosetta open-vs-closed ddG (Arm 2 / 2b)", "Pocket - energetics",
 "Relax each variant in open and closed states; the ddG between them proxies the energetic cost of switching.",
 "Tried to score whether an enlarged pocket still prefers the closed, HAB1-competent state.",
 "ABA; 3QN1 vs 3K3K; 41 variants at nstruct 3, 20, 12", "Aug 7-10",
 "Arm 2b reached F = 9.0 with a true spread of 1.79 REU.",
 "Arm 2 v1 was declared statistically blind - wrongly. The readout is INVARIANT, not blind: true spread 0.43-0.50 REU against a -67 REU interface.",
 "A null is only a null if the assay can discriminate. Judge effects in physical units, not significance.",
 "OPTION", "later", "CPU cutlerlab", "22, 26, 32, 33; jobs 27275686, 27286483/4"],

["Corrected discrimination statistic", "Pocket - statistics",
 "F = s2_between / (s2_within / n), where s2 is the VARIANCE (standard deviation squared). The denominator s2_within/n is the squared standard error of the mean, so dividing the within-group variance by n IS the sqrt(n) correction written in variance units - the two are the same equation, not two.",
 "Re-derived every arm's discrimination after finding the original statistic wrong.",
 "Re-analysis of existing Rosetta output; no new structures.", "Aug 10",
 "Flipped Arm 2b from 'blind' to F = 9.0, correcting a central project claim.",
 "The original statistic used s_within itself as the yardstick - it compared between-group spread to the scatter of a SINGLE replicate rather than to the scatter of the replicate MEAN, which is sqrt(n) times smaller. The error had already propagated into several conclusions before it was found.",
 "Compare between-group spread against the SEM, not the replicate sd. Run the variance decomposition BEFORE concluding 'no effect'.",
 "ACTIVE", "-", "CPU local", "35_variance_decomposition.py"],

["Arm 3 homolog cavity survey", "Pocket - structure",
 "Full-atom cavity measurement across foldseek homologs.",
 "Asked how PYR1's pocket compares to its structural family.",
 "147 full-atom homolog structures", "Aug 6-7",
 "Placed PYR1's 174 A^3 pocket in family context; supported the CoxG lead.",
 "First pass returned 0.0 A^3 for every cavity - foldseek stores Ca only, so convert2pdb emitted backbone traces. Silent and plausible-looking.",
 "Validate a downloaded structure with the consuming tool before trusting a derived number. A zero can be a parsing artefact.",
 "OPTION", "later", "CPU cutlerlab", "08, 10"],

["LigandMPNN", "Pocket - design",
 "Graph neural network designing sequence conditioned on backbone AND explicit ligand atoms; outputs per-position amino-acid probabilities.",
 "Stage-1 retrospective test: can it recover the 4 known mandipropamid mutations from WT plus a perfectly placed ligand?",
 "Structures: WT PYR1 with mandipropamid (ligand 3UZ, from 4WVO) vs ABA (ligand A8S, from 3QN1) as the ligand-swap null.", "Aug 12",
 "Recovered F108A (+0.111) and F159L (+0.071) with zero mass in the ABA arm - genuinely ligand-conditional.",
 "Missed K59R and V81I. V81I INVERTED (-0.841): ABA prefers Ile at 81 with p=0.996 in the arm whose answer is no mutations. Recovers only severe-clash positions, which a trivial clash baseline already gives free.",
 "Sequence recovery is the wrong metric - a ligand-blind oracle reaches 79.1% on the same labels. Score probability MASS against a ligand-swap null, never recall-at-N.",
 "OPTION", "later", "GPU", "43, 46; job 27412775"],

["Rosetta FastDesign", "Pocket - design",
 "Monte Carlo sequence design with ref2015: repack and minimise while allowing identity changes at chosen positions.",
 "The physics arm of the same stage-1 test, chosen because it fails differently from a learned model.",
 "Mandipropamid vs ABA, 15 designable positions", "Aug 13",
 "Recovered F108A at frequency 1.000 vs 0.000 in the null - the cleanest single recovery of the stage.",
 "F159L: right position, wrong identity (picks D/E). V81I inverted (-0.980). The NULL ARM FAILED ITS OWN CONTROL: with the cognate ligand it kept WT at 29% and K59 in 0% of trajectories.",
 "If the null also mutates a position, 'no signal' is inseparable from 'the protocol cannot hold a native contact'. K59R and V81I were uninterpretable, not negative.",
 "OPTION", "later", "CPU cutlerlab", "50, 51; job 27421344"],

["Favour-native weight sweep", "Pocket - calibration",
 "Adds a bonus for retaining wild-type identity; swept 0.25/0.5/1.0/1.5 on the NULL arm only, criterion fixed before submission.",
 "Attempt to repair FastDesign's WT-retention failure.",
 "ABA null arm", "Aug 13",
 "WT recovery climbed 29 -> 45 -> 59 -> 65 -> 85%.",
 "K59 stayed at 0% AT EVERY WEIGHT, including 1.5 where the packer is nearly frozen (3 distinct sequences in 10). No weight passed.",
 "The calibration failing IS the answer, not a setback. Do not widen a sweep until a control agrees - that is tuning until it says what you want.",
 "RETIRED", "-", "CPU cutlerlab", "53, 54; job 27438220"],

["Per-residue energy decomposition", "Pocket - diagnosis",
 "Breaks ref2015's total into per-term, per-residue contributions at one position.",
 "Asked WHY K59 could not be retained.",
 "ABA anion + K59 in 3QN1", "Aug 13",
 "THE KEY DIAGNOSIS. ref2015 correctly sees the salt bridge (Lys best fa_elec by 1.9 REU, only candidate earning hbond_sc) then pays +10.118 REU of Lazaridis-Karplus desolvation to bury the ammonium. Ile beats Lys by 4.211 REU; Arg is 2.274 WORSE than Lys.",
 "Ruled out pose (K59 NZ 2.85 A from carboxylate) and ionisation (params sum -0.970 over 38 atoms) first, costing extra runs.",
 "Diagnose the mechanism before proposing a fix. This one number explained three downstream failures and told us exactly which stage to replace.",
 "ACTIVE", "-", "CPU cutlerlab", "55, 55b; job 27438909"],

["Coupled moves (Kortemme)", "Pocket - design",
 "Monte Carlo coupling sequence, side-chain, BACKBONE and ligand-pose moves - the published fix for fixed-backbone specificity design (5.75x on its own benchmark, 16/17 cases).",
 "Re-ran stage 1 changing ONLY the sampler: same positions, same arms, same scoring rule.",
 "Mandipropamid vs ABA; 50 trajectories/arm x 1000 trials, ligand_mode on", "Aug 19",
 "F108A recovered cleanly (+0.74 against a 0.00 null). WT retention improved 29% -> 42%.",
 "K59R still missed. Decisively: in the ABA arm - cognate ligand, K59 unambiguously correct - K59 was retained in 0% of 50 trajectories, identical to fixed backbone.",
 "SAMPLING IS NOT THE LIMIT at the design stage. Lysine is in every rotamer library and IncludeCurrent was on. Confirms the scoring diagnosis with a second independent sampler.",
 "OPTION", "later", "CPU cutlerlab, 100 tasks", "74, 74b; job 27560390"],

["Per-position steric admissibility", "Pocket - geometry",
 "For each position x each of 20 residues, place every library rotamer and ask whether ANY is clash-free against ligand and environment. Hard yes/no on Bondi radii, no score function.",
 "Headroom check: does pure geometry narrow the residue choice per position?",
 "ABA vs mandipropamid, 26 pocket positions, WT backbone", "Aug 18",
 "Clean negative in an afternoon: mean 14.7/20 admissible, ABA-mandi Jaccard 0.84.",
 "DOES NOT CLEAR. Too permissive (9x narrowing where Tian's focused library is 50x smaller) and largely ligand-blind. All four ground-truth substitutions are admissible for BOTH ligands - recalled without being ligand-conditional.",
 "YES, IT EXCLUDED REAL SENSOR MUTATIONS - but read the scope carefully. Against mandipropamid's OWN ground truth it excluded nothing: all four of K59R/V81I/F108A/F159L were admissible (the problem there was that they were admissible for ABA too). Against Beltran's cannabinoid sensors it wrongly excluded six real substitutions - V83W, V83F, V164W, V164H, A89V, E141F - which is 7% of in-pocket occurrences. Caveat: those were scored in the ABA/mandipropamid pocket, not in a cannabinoid-bound one, so it is not a clean test of them. Every miss was a GROW mutation, and real sensors are compensating shrink/grow PAIRS (F159A+A160I; Y120G+A160G in 27 and 23 of Beltran's 45 sensors) that a per-position filter structurally cannot see.",
 "RETIRED", "-", "CPU local", "71_admissibility.py"],

["Pairwise enumeration", "Pocket - geometry",
 "Enumerates mutation PAIRS against a jointly re-scored pocket. Objective = summed ligand-protein vdW overlap (additive, hence a cheap precompute) subject to not losing packing.",
 "The version of the combinatorial idea the data supports, after per-position menus failed.",
 "Mandipropamid (clashes in WT: 2.78 A, 7/29 atoms) vs ABA (fits: 0.22 A, 0/19) - a built-in negative control", "Aug 18",
 "Best steric result so far: F108A rank 5, F159L rank 75, V81I rank 150 of 13,457 - 3 of 4 in the top 150, a 90x narrowing. Shrink/grow enriched to 62% of the top 200 vs 42% for ABA.",
 "K59R at rank 3907 - blind BY CATEGORY, since K59 does not clash at all. Picks the right POSITIONS (81+108) with the wrong identity at 81.",
 "An ADDITIVE objective cannot show synergy. Counting contacts within a flat 4.5 A rewards clashes (83W+163Y topped the Pareto front with relief -71.7 A).",
 "ACTIVE", "-", "CPU local", "72_pairwise.py"],

["Electrostatic / H-bond complementarity", "Pocket - geometry",
 "Counts geometrically valid donor-acceptor pairs to the ligand plus a screened-Coulomb term. Deliberately models NO desolvation.",
 "The stage pairwise enumeration specified: score the non-clashing positions electrostatics can reach.",
 "ABA (anion, -0.97) vs mandipropamid (neutral, +0.10)", "Aug 18",
 "Gets K59 right for ABA - Lys rank 1 via a bidentate salt bridge (NZ-O4 2.83 A, NZ-O3 3.31 A) where stage 1 kept it 0% of the time. The ref2015 pathology is genuinely absent.",
 "Missed K59R. Fails its own controls in the PREDICTED direction: WT top choice at only 4% of positions, Lys/Arg picked at 59% against an anionic ligand - with no desolvation, everything wants to be charged.",
 "Not usable as a ranker; usable as a FILTER over supplied conformations. Its own miss was sampling: the crystal Arg59 (chi3 ~100 deg, non-rotameric) scores 2 H-bonds under its own criteria but is never proposed.",
 "OPTION", "later", "CPU local", "73_electrostatic.py"],

["Shared steric rule (lib_sterics)", "Pocket - geometry",
 "One clash rule for all geometric screens: exempts donor-acceptor pairs down to 2.5 A, keeps acceptor-acceptor and donor-donor as clashes, includes a real hydrogen radius.",
 "Fixed hydrogen bonds being scored as steric clashes across scripts 71/72/73.",
 "All pocket screens", "Aug 19",
 "Unit-tested three ways. R79's crystal side chain goes from 0.306 A 'overlap' to 0.000. WT-ok at hard-sphere tolerance improved 12 -> 16 of 26.",
 "Did NOT rescue the R79/V83/H115 control failures - I had claimed it would. Those are rotamer-basin failures and the earlier diagnosis was right.",
 "A hard-sphere test scores every hydrogen bond as a ~0.4 A clash. The fix was invisible at the 0.5 A tolerance in use, which already absorbed it - a stability check, not a revision.",
 "ACTIVE", "-", "CPU local", "lib_sterics.py"],

["MM-GBSA rescoring", "Pocket - energetics",
 "Replaces ref2015's pairwise Lazaridis-Karplus solvation with generalised Born: dG_bind = G_complex - G_receptor - G_ligand, averaged over an MD ensemble. ff14SB + GAFF2/AM1-BCC, igb=8, 0.15 M salt. Rosetta builds structures and never scores them.",
 "Direct attack on the diagnosed cause of the K59 failure.",
 "Mandipropamid vs ABA; 7 variants x 2 arms x 100-frame ensembles", "Aug 19",
 "THE K59 FLIP. With mandipropamid R beats Q and N (-1.53 vs +0.04, +2.53); with ABA wild-type Lys beats every substitution (R least badly, +7.06). V81I and F159L also prefer mandipropamid.",
 "Ensembles only 50 ps (GBn2 runs ~380 steps/min here; pmemd no faster). The ligand-swap difference is dominated by DAMAGE TO THE ABA COMPLEX, not gain for mandipropamid. F108A comes out NULL because relaxation absorbs the clash it exists to relieve.",
 "Changing the solvation model does what no sampler could - the desolvation diagnosis was right AND actionable. Also: my first error bars were fabricated (8 'independent' repacks came out byte-identical).",
 "ACTIVE", "-", "CPU cutlerlab, 14 tasks", "75, 75b-e; jobs 27561665, 27568930, 27569546"],

["Coumarin focused-library benchmark", "Pocket - benchmark",
 "Extracts the exact per-position residue menus Tian et al. built for each focused library, then measures what a method must add over a ligand-blind baseline.",
 "Asked whether round-1 screening could have been skipped computationally.",
 "LIBRARY DESIGNS, not structures: Tian sd04 per-position residue menus for the Coumarin, PFAS and TNTv1/v2 focused libraries, plus 692 characterised sensor clones (sequences).", "Aug 14",
 "Reframed the task. POSITIONS are near ligand-independent: Coumarin/PFAS/TNTv2 share 9 of their 11-14 positions, so naming the shared nine scores 82% recall with ZERO ligand information. SUBSTITUTIONS are ligand-specific: mean 3-way Jaccard 0.02.",
 "I initially claimed the coumarin library design was missing from the supplement. It was not - a column header truncated at 20 characters.",
 "Predict MENUS, not positions. Focusing buys ~16 orders of magnitude (3.8e21 -> 1.4e5).",
 "ACTIVE", "-", "CPU local", "62_coumarin_ground_truth.py"],

["Beltran cannabinoid sensor table", "Pocket - benchmark",
 "45 experimentally validated sensors across 12 cannabinoids, with mutation sets and a graded dose-response (1e4 -> 10 nM).",
 "A second, harder ground truth with an ORDER to recover, not just a set.",
 "Sensor SEQUENCE table (no solved structures): 45 PYR1 variants vs WIN 55,212-2, the JWH series, D9-THC, CBDA, CP 47,497. WIN 55,212-2 is an aminoalkylindole synthetic cannabinoid, so it belongs in this set despite the name looking unlike the others.", "Aug 18-19",
 "Extracted in full. A160 mutated in 38/45, Y120 in 33/45, F159 in 24/45. Ligand-specific signal sits at K59 (Q/N/S/T/A/R), H115Q, V83, E141, V164, V81.",
 "Even more position-saturated than coumarin; substitutions reused within a chemical family (A160G 27x, Y120G 23x). Nine substitutions sit OUTSIDE the 5 A pocket - PYR1-WIN gets its last order of magnitude from E4G.",
 "A graded sensitivity ladder is a far stronger benchmark than a hit list. K59Q/K59N are neutral, so they should be reachable where K59R was not - a free sanity check.",
 "ACTIVE", "-", "CPU local", "data/beltran/win_sensors.json"],

["Literature protocol review", "Pocket - benchmark",
 "Read six papers from their Methods sections, not abstracts, to establish what the field does and what the honest baseline is.",
 "Calibrated our results against published practice before investing further.",
 "Six papers (PDFs + extracted Methods text): Leonard 2026, An 2024, Tian 2025, Park 2023, ACS Chem Biol 2024, Nat Biotech 2026.", "Aug 19",
 "Leonard's 'dock to sequence' MUTATES WT to a known weak-hit sequence before placing any pose - their method improves a known weak binder and cannot start from a ligand alone. An et al. picked cutoffs by MANUAL INSPECTION. RbsB precedent: 2M variants screened for 1.2-1.5x induction.",
 "Leonard hard-filters on an H-bond to the latch water, which mandipropamid misses by 5.07-5.19 A - their protocol would reject the correct pose.",
 "The field's baseline is low, and our benchmark is strictly harder than Leonard's. The ACS Chem Biol paper independently confirmed R79 H-bonds F52's carbonyl at >90% occupancy, explaining one of our control failures from the other direction.",
 "ACTIVE", "-", "local", "data/papers/"],

["WT MD baseline / loop dynamics", "MD",
 "Explicit-solvent MD (ff19SB/OPC, 0.15 M KCl, 2 fs, 300 ns) with a two-reference projection: fit on the rigid core, then measure gate/latch RMSD to BOTH open and closed crystals.",
 "Asked whether open and closed states are separable and whether MD can filter designs.",
 "S1 apo-open (3K3K A), S2 holo-closed (3QN1 A + ABA), S3 dimer, S4 ternary (+HAB1, Mn2+); 12 trajectories", "Aug 10 - Aug 18",
 "Neither state ever converts in 1.8 us aggregate. Replicate-mean gap 5.57 A at the gate, stable at every discard from 0-150 ns. RMSF ladder: gate 3.05 open -> 1.39 closed -> 0.81 ternary.",
 "Licenses a STABILITY filter only (gate RMSD-to-closed, 3.5-4.0 A) and is BLIND to switchability - the likelier failure for an enlarged pocket. The latch fails as an RMSF filter (ratio inverts).",
 "Never read a stable closed trajectory as a working switch. Inference is at replicate level: autocorrelation times 3.4-44.9 ns, so each 300 ns run holds only 3-35 independent samples.",
 "ACTIVE", "-", "GPU", "58-60; jobs 27333712, 27545237"],

["Apo-closed protomer analysis", "MD",
 "Re-analysis of the S3 mixed dimer, whose chain B is a closed, ligand-shaped protomer simulated with the ABA stripped.",
 "Filled the missing cell of the factorial using data already on disk.",
 "3K3K chain B (closed, was ABA-bound), 3 x 300 ns", "Aug 18",
 "Apo-closed does NOT open: zero crossings, drift -0.07 A, and the MOST RIGID gate of all 15 units (RMSF 0.72 A vs 1.39 with ABA).",
 "Confounded by the dimer, which alone drops the open gate 3.05 -> 1.07 A - larger than the ligand effect being measured.",
 "The stability filter is CONFORMATION-reporting, not ligand-reporting: a good gate-RMSD score carries no evidence about occupancy. Also: 3K3K is a MIXED dimer, and the name S3_apo_dimer had outranked four records saying so.",
 "ACTIVE", "-", "GPU", "58-60; job 27545237"],

["191-residue rebuild pipeline", "MD - structure",
 "Biopython graft + PyRosetta rebuild of the complete 191-aa sequence with no chain breaks, seeded and byte-reproducible, with identity and provenance assertions at every step.",
 "Removed the gapped-receptor confound from every future MD system.",
 "3K3K A/B, 3QN1 A; ABA", "Aug 17",
 "7 verified structures; every chain 191 residues, 0 peptide-bond breaks, 0 clashes, ligand-lining side chains unmoved to 0.000 A.",
 "Cost two structural defects before the cause was found - a 0.36 A inter-chain clash and K59 collapsing 2.85 -> 1.68 A into an empty cavity.",
 "FastRelax.set_movemap() restricts MINIMISATION only; repacking needs a TaskFactory or the whole pose is repacked. Backbone RMSD reported 0.000 A throughout - always check side chains and inter-chain contacts.",
 "ACTIVE", "-", "CPU local", "67, 67c, 68; lib_resnumber/structqc/rosetta"],

["Conformation x occupancy factorial", "MD",
 "The 2x2 that README 19b set up and only ever filled on the diagonal: open/closed x apo/holo, all four cells on ONE build protocol.",
 "Separates 'the gate holds its shape' from 'the ligand holds the gate'.",
 "S1 open+apo, S2 closed+ABA, S9 closed+apo, S10 open+ABA (ABA transplanted, backbone bit-identical)", "Aug 18 - running",
 "S10 (open + ABA) never existed before and is now built and verified. 5 of 12 replicates complete, including both S10 reps 0 and 2.",
 "Not yet analysed. The old-tree S9 could not be reused: it is KCl where S1/S2 are NaCl, which would put a cation swap inside the key comparison.",
 "Gate closure on ligand binding is the one arm where the barrier is plausibly downhill - pre-registered so 'nothing happened' cannot later be written up as 'ABA does not close the gate'.",
 "ACTIVE (running)", "now", "GPU", "69, 69b, 70; job 27547457"],

["Non-cognate ligand MD (S6-S8)", "MD",
 "Same closed receptor, three ABA-sized ligands from three chemical classes, three docked poses each.",
 "Designed to test whether the stability filter can rank LIGANDS or only conformations.",
 "Imperatorin, Flutamide, alpha-Estradiol (19-20 heavy atoms, MW 270-276) in WT closed PYR1", "Aug 14 (built)",
 "Built and pose-verified; three replicates per ligand start from three DIFFERENT poses, so pose sensitivity is measured rather than assumed.",
 "Never submitted. Partly superseded by the apo-closed result, which answered the same question more cheaply.",
 "Absence of ligand release is NOT evidence of binding - residence times are us-ms and 300 ns cannot sample unbinding.",
 "PLANNED", "P2 - after the factorial", "GPU", "63-66"],

["ABA enantiomer / non-cognate controls", "MD",
 "Isosteric control: (-)-ABA is identical in formula, bulk and parameters, differing only in chirality, so any difference is chemistry rather than size.",
 "Requested, to test whether gate behaviour is ligand-specific.",
 "(-)-ABA plus an unrelated matched molecule; WT PYR1", "planned",
 "Not yet run.",
 "Must NOT be read as a gate-opening test: nothing opens on this timescale in any of 15 units. The graded RMSF ladder is the observable that responds to ligand.",
 "Phaseic acid would add a third rung, turning a binary contrast into a graded ladder that is much harder to explain away.",
 "PLANNED", "P2 - after the factorial", "GPU", "-"]
]


ATTRIBUTION = (
    "Throughout this document, \"I\" and \"my\" refer to Claude (Anthropic), the AI "
    "model that ran these analyses and wrote this file - not to the repository owner. "
    "Jannis Jacobs set the direction, chose the experiments, supplied the literature "
    "and caught several of the errors recorded here; the implementation, the mistakes "
    "written in the first person, and this summary are the model's. See README \u00a70 "
    "for the full per-model attribution."
)

TIMELINE = [
    ("Aug 6",  "Project start. Pocket characterised; F108 identified as the wall sealing the satellite lobe."),
    ("Aug 6-8","Rosetta arms 1, 2, 2b: 41 variants relaxed in open and closed states."),
    ("Aug 10", "The discrimination statistic is found to be WRONG. Arm 2b flips from 'blind' to F = 9.0."),
    ("Aug 10", "WT MD campaign submitted - 4 systems x 3 replicates x 300 ns."),
    ("Aug 11", "W385 water re-measured across 3QN1/4WVO/8EY0: it senses ABA but mandipropamid never touches it."),
    ("Aug 12", "Stage 1, LigandMPNN arm: recovers F108A and F159L, misses K59R and V81I."),
    ("Aug 13", "Stage 1, Rosetta FastDesign arm: same two recovered. The NULL arm fails its own control (K59 kept 0%)."),
    ("Aug 13", "Favour-native sweep fails at every weight. Per-residue decomposition finds the cause: +10.1 REU desolvation."),
    ("Aug 14", "Coumarin benchmark reframes the task: positions are saturated, SUBSTITUTIONS carry the ligand information."),
    ("Aug 17", "PYR1 rebuilt on the full 191-residue sequence. FastRelax's packer trap found and closed."),
    ("Aug 18", "3K3K exposed as a MIXED dimer; the apo-closed cell turns out to have been simulated by accident."),
    ("Aug 18", "Per-position admissibility fails; pairwise enumeration gets 3 of 4 into the top 150 of 13,457."),
    ("Aug 18", "Conformation x occupancy factorial built and launched - the open+ABA cell had never existed."),
    ("Aug 19", "Six papers read from their Methods. Leonard's protocol needs a known weak hit before it can place a pose."),
    ("Aug 19", "Coupled moves: better sampling, K59 still retained 0% of the time. Sampling is NOT the limit."),
    ("Aug 19", "MM-GBSA rescoring makes the K59 flip. Changing the solvation model did what no sampler could."),
]

GLOSSARY = [
    ("STRUCTURES (PDB entries)", "", ""),
    ("3QN1", "PYR1 + ABA + HAB1, closed ternary complex, 1.80 A.",
     "The reference CLOSED state and the source of the ABA pose used throughout."),
    ("3K3K", "PYR1 homodimer, 2 chains.",
     "A MIXED dimer, not an apo one: chain A is apo and OPEN, chain B is CLOSED with ABA bound. The reference OPEN state comes from chain A."),
    ("4WVO", "PYR1-MANDI + mandipropamid + HAB1, 2.25 A.",
     "The engineered mandipropamid receptor. Its 4 mutations vs 3QN1 - K59R, V81I, F108A, F159L - are the ground truth every design method is tested against."),
    ("8EY0", "Orthogonalised PYR1*/HAB1* + mandipropamid, 2.40 A.",
     "Independent second mandipropamid structure; also a rare validated true-negative protein pair."),
    ("7MWN", "PYL2-WIN + WIN 55,212-2.",
     "Cannabinoid complex in PYL2, a close PYR1 relative (88% identity over the 25 ABA-proximal positions)."),
    ("3K90", "PYR1, 4 chains, mixed apo and ABA-bound.", "Used in the satellite-lobe survey across all deposited chains."),
    ("LIGANDS", "", ""),
    ("ABA / A8S", "Abscisic acid; PDB chemical component A8S. Anionic at pH 7 (net -1).",
     "The COGNATE ligand. Because WT PYR1 evolved for it, it serves as the ligand-swap NULL: a method that prefers a mutation for ABA too is not reading the ligand."),
    ("Mandipropamid / 3UZ", "A fungicide, 29 heavy atoms, neutral.",
     "The non-cognate test ligand. It CLASHES in the WT pocket (2.78 A, 7 of 29 atoms), which is the signal the steric methods exploit."),
    ("WIN 55,212-2", "An aminoalkylindole synthetic cannabinoid.", "Target of Beltran's best cannabinoid sensor, PYR1-WIN."),
    ("PROGRAMS AND FORCE FIELDS", "", ""),
    ("Rosetta / PyRosetta", "All-atom protein modelling suite.", "Used for relax, design and structure building throughout."),
    ("ref2015", "Rosetta's default all-atom energy function.",
     "Its Lazaridis-Karplus desolvation term is the specific thing that made K59R unreachable - it pays +10.1 REU to bury a charged ammonium."),
    ("FastRelax / FastDesign", "Rosetta protocols: repack + minimise, with (Design) or without (Relax) allowing identity changes.",
     "TRAP: set_movemap() restricts MINIMISATION only. Repacking needs a TaskFactory or the whole pose is repacked."),
    ("Coupled moves", "Rosetta sampler (Ollikainen & Kortemme 2015) coupling sequence, side-chain, backbone and ligand-pose moves.",
     "The published fix for fixed-backbone specificity design. Tried here; did not rescue K59R."),
    ("LigandMPNN", "Neural network designing sequence conditioned on backbone plus explicit ligand atoms.", "The learned-model arm of stage 1."),
    ("Amber / ff19SB / ff14SB", "Molecular dynamics engine and its protein force fields.",
     "ff19SB (with OPC water) for explicit-solvent MD; ff14SB for the implicit-solvent MM-GBSA rescoring, which is what that method was validated with."),
    ("GAFF2 / AM1-BCC", "General small-molecule force field and its charge model.", "Used to parameterise ABA and mandipropamid for Amber."),
    ("igb=8 / GBn2", "A generalised-Born implicit solvent model.",
     "The alternative solvation treatment that finally produced the K59 flip. Requires mbondi3 radii or sander refuses to run."),
    ("MM-GBSA", "Molecular Mechanics + Generalised Born Surface Area: dG_bind = G_complex - G_receptor - G_ligand, averaged over an ensemble.",
     "An ensemble RESCORING method, not a design method - it ranks structures someone else generated."),
    ("METRICS", "", ""),
    ("REU", "Rosetta Energy Unit - Rosetta's arbitrary energy scale, loosely kcal/mol.", "Comparable within one score function, not across."),
    ("ddG (delta-delta-G)", "The change in a free-energy difference between two states or two variants.",
     "Here usually mutant-minus-WT within one ligand arm, so receptor and force field cancel."),
    ("RMSD / RMSF", "Root-mean-square deviation (from a reference) / fluctuation (about a mean).",
     "RMSD says where a loop IS; RMSF says how much it MOVES. The gate needs both."),
    ("Jaccard index", "Overlap of two sets: |A and B| / |A or B|. 1.0 = identical, 0 = disjoint.",
     "Used to ask whether two ligands get the same residue menu. 0.84 means largely ligand-blind."),
    ("kBT", "Thermal energy at the simulation temperature, ~0.6 kcal/mol at 300 K.", "A natural unit for 'is this conformer accessible'."),
    ("Y2H", "Yeast two-hybrid: growth or reporter readout of a protein-protein interaction.",
     "How Tian and Beltran measure sensor activation - the ligand-dependent PYR1-HAB1 interaction drives the reporter."),
    ("STATISTICS", "", ""),
    ("s2 (s-squared)", "The VARIANCE - the square of the standard deviation. s2_between is the variance of the group means; s2_within is the variance within groups.",
     "Written s2 here because these are variances, not standard deviations. The distinction is what the corrected F statistic turns on."),
    ("F = s2_between / (s2_within / n)", "The one-way ANOVA F statistic.",
     "The denominator s2_within/n is the SQUARED STANDARD ERROR OF THE MEAN. Dividing by n IS the sqrt(n) correction, expressed in variance units - so there is only one equation here, not two. The original error was using s_within itself as the yardstick."),
    ("ICC(1)", "Intraclass correlation - the fraction of total variance that is between-subject rather than within-subject.",
     "Used to ask whether a score reproduces across random seeds."),
    ("Ligand-swap null", "Running the identical protocol with the COGNATE ligand and subtracting.",
     "The project's standing scoring rule: only f(test) - f(cognate) is evidence. It inherits the method's own biases, so it needs no independence claim."),
]

STATUS_FILL = {"ACTIVE": "C6E7C6", "ACTIVE (running)": "BFE3F5",
               "OPTION": "FFF0C2", "PLANNED": "E4DAF2", "RETIRED": "F2D4D4"}


def write_xlsx(path):
    from openpyxl import Workbook
    from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
    from openpyxl.utils import get_column_letter
    wb = Workbook()

    # sheet 1: how to read this, plus the timeline
    tl = wb.active
    tl.title = "Read me + timeline"
    tl.column_dimensions["A"].width = 16
    tl.column_dimensions["B"].width = 120
    tl.append(["Methods review", "Every computational method used in the PYR1 pocket-expansion project"])
    tl.cell(row=1, column=1).font = Font(bold=True, size=14)
    tl.append([])
    tl.append(["Attribution", ATTRIBUTION])
    tl.cell(row=3, column=1).font = Font(bold=True)
    tl.cell(row=3, column=2).alignment = Alignment(wrap_text=True, vertical="top")
    tl.row_dimensions[3].height = 76
    tl.append([])
    tl.append(["Status", "ACTIVE = in use now or load-bearing | OPTION = works, kept for a purpose | "
               "PLANNED = built or specified, not yet run | RETIRED = tested and set aside, reason recorded"])
    tl.cell(row=5, column=1).font = Font(bold=True)
    tl.cell(row=5, column=2).alignment = Alignment(wrap_text=True, vertical="top")
    tl.append([])
    tl.append(["TIMELINE", "2026"])
    tl.cell(row=7, column=1).font = Font(bold=True, color="FFFFFF")
    tl.cell(row=7, column=1).fill = PatternFill("solid", fgColor="3D4B5C")
    tl.cell(row=7, column=2).font = Font(bold=True, color="FFFFFF")
    tl.cell(row=7, column=2).fill = PatternFill("solid", fgColor="3D4B5C")
    for when, what in TIMELINE:
        tl.append([when, what])
        r = tl.max_row
        tl.cell(row=r, column=1).font = Font(bold=True)
        tl.cell(row=r, column=2).alignment = Alignment(wrap_text=True, vertical="top")

    ws = wb.create_sheet("Methods")
    thin = Side(style="thin", color="BBBBBB")
    ws.append(COLS)
    for c in range(1, len(COLS) + 1):
        cell = ws.cell(row=1, column=c)
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="3D4B5C")
        cell.alignment = Alignment(vertical="center", wrap_text=True)
    for r in ROWS:
        ws.append(r)
    widths = [30, 20, 62, 46, 46, 16, 62, 62, 62, 16, 22, 22, 34]
    for i, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(i)].width = w
    si = COLS.index("Status") + 1
    for row in range(2, len(ROWS) + 2):
        for c in range(1, len(COLS) + 1):
            cell = ws.cell(row=row, column=c)
            cell.alignment = Alignment(vertical="top", wrap_text=True)
            cell.border = Border(left=thin, right=thin, top=thin, bottom=thin)
        st = ws.cell(row=row, column=si)
        st.fill = PatternFill("solid", fgColor=STATUS_FILL.get(st.value, "FFFFFF"))
        st.font = Font(bold=True)
        ws.cell(row=row, column=1).font = Font(bold=True)
        ws.row_dimensions[row].height = 108
    ws.row_dimensions[1].height = 34
    ws.freeze_panes = "B2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(COLS))}{len(ROWS)+1}"
    # sheet 3: glossary
    gl = wb.create_sheet("Glossary")
    for i, w in enumerate([34, 62, 78], start=1):
        gl.column_dimensions[get_column_letter(i)].width = w
    gl.append(["Term", "What it is", "Why it matters here"])
    for c in range(1, 4):
        cell = gl.cell(row=1, column=c)
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="3D4B5C")
    for term, what, why in GLOSSARY:
        gl.append([term, what, why])
        r = gl.max_row
        if not what and not why:      # a section heading
            gl.cell(row=r, column=1).font = Font(bold=True, color="FFFFFF")
            for c in range(1, 4):
                gl.cell(row=r, column=c).fill = PatternFill("solid", fgColor="7A8899")
        else:
            gl.cell(row=r, column=1).font = Font(bold=True)
            for c in (2, 3):
                gl.cell(row=r, column=c).alignment = Alignment(wrap_text=True, vertical="top")
    gl.freeze_panes = "A2"

    wb.save(path)
    return path


def write_md(path):
    order = ["ACTIVE", "ACTIVE (running)", "OPTION", "PLANNED", "RETIRED"]
    out = ["# Methods review — every computational method used in this project", "",
           "One row per method: what it computes, how we used it, what it got right,",
           "what it got wrong, and whether it is still in play.", "",
           "> **A note on voice.** " + ATTRIBUTION.replace('"', "'"), "",
           "## Timeline", "",
           "| 2026 | what happened |", "|---|---|"]
    out += [f"| **{w}** | {t} |" for w, t in TIMELINE]
    out += ["",
           "## Glossary", "",
           "For readers less familiar with the structures, programs and notation used below.", ""]
    for term, what, why in GLOSSARY:
        if not what and not why:
            out += ["", f"**{term}**", "", "| term | what it is | why it matters here |", "|---|---|---|"]
        else:
            out += [f"| `{term}` | {what} | {why} |"]
    out += ["",
           "## Status vocabulary", "", "Fixed so it means the same thing in every row.", "",
           "| Status | Meaning |", "|---|---|",
           "| **ACTIVE** | in use now, or its result is load-bearing for what comes next |",
           "| **OPTION** | works; kept for a specific future purpose or to couple with another method |",
           "| **PLANNED** | built or specified, not yet run |",
           "| **RETIRED** | tested and set aside — the row records why, so it is not re-tried by accident |",
           "", f"Generated by `scripts/76_methods_review.py`. Spreadsheet: `data/methods_review.xlsx`.", ""]
    for st in order:
        rows = [r for r in ROWS if r[COLS.index("Status")] == st]
        if not rows:
            continue
        out += [f"## {st} ({len(rows)})", ""]
        for r in rows:
            d = dict(zip(COLS, r))
            out += [f"### {d['Method']}", "",
                    f"**{d['Workstream']}** · {d['Dates (2026)']} · {d['Where it ran']} · `{d['Scripts / jobs']}`"
                    + (f" · **priority: {d['Priority']}**" if d['Priority'] not in ('-', '') else ""), "",
                    f"- **How it computes:** {d['What it computes (how)']}",
                    f"- **How we used it:** {d['How we used it']}",
                    f"- **Inputs:** {d['Inputs (structures / ligands / libraries)']}",
                    f"- **Key result:** {d['Key result']}",
                    f"- **Shortcomings & failures:** {d['Shortcomings & failures']}",
                    f"- **What we learned:** {d['What we learned']}", ""]
    open(path, "w").write("\n".join(out))
    return path


if __name__ == "__main__":
    x = write_xlsx(os.path.join(ROOT, "data", "methods_review.xlsx"))
    m = write_md(os.path.join(ROOT, "METHODS_REVIEW.md"))
    from collections import Counter
    c = Counter(r[COLS.index("Status")] for r in ROWS)
    print(f"{len(ROWS)} methods, {len(COLS)} columns")
    for k, v in c.most_common():
        print(f"  {k:<16} {v}")
    print(f"\nwrote {x}\nwrote {m}")
