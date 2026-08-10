"""
variants.py -- the pilot mutant panel, shared by 06 (cavity) and 07 (ratchet).

Position rationale (PYR1 numbering, PDB 3QN1 chain A)
-----------------------------------------------------
K59   ABA carboxylate salt bridge, NZ-O12 2.85 A. Top-2 restrictor
      (+71.8 A3 rigid truncation). NOT conserved in the helix-grip/SRPBCC
      superfamily -- Lys only 10% of 235 ungapped homologs; consensus is Asp
      (15%), hydrophobics L/V/I/A ~28%. Read as a PYR/PYL-specific adaptation
      for the ABA carboxylate, therefore the safest large-gain position.
F108  Top restrictor (+73.5 A3). Ring centroid stacks on the R79 guanidinium
      at 3.55 A (cation-pi). Superfamily prefers BULKY AROMATIC here
      (Y 24%, F 14%, W 7%); Tyr is more common than Phe, i.e. the consensus is
      LARGER than PYR1. Likely load-bearing -- included with F108Y as an
      "against the grain" control.
R79   Minimum heavy-atom distance to ABA 5.81 A (NH2...O11) -- just OUTSIDE a
      5 A contact-radius pocket definition, so a first-shell scan does not list
      it, though it is still pocket-proximal. Salt-bridges E94 and pins the
      F108 rotamer (3.55 A cation-pi). The structural reason F108 may not relax
      when mutated alone.
E94   FIRST SHELL: OE2 is 4.33 A from ABA C15, so it occludes cavity volume
      directly, AND forms the other half of the R79 salt bridge
      (OE1-NE 2.76 A, OE1-NH2 2.93 A). Both roles are in play; the R79A/E94A
      double vs the two singles separates them.
      (All distances are side-chain heavy atom to ligand heavy atom. The
       corresponding CA-to-ligand distances are 6.6-11.7 A and are not used.)
Y120  Third-ranked restrictor (+46.0 A3), latch-adjacent -- watch for ratchet
      loss.
E141  Second-ranked restrictor (+69.4 A3), highly variable in the superfamily
      (V 15%, E 11%, L 11%, W 11%).

R79/E94/E141/K59 are all INVARIANT across close PYR/PYL homologs (alnTM>0.89)
but VARIABLE across the wider superfamily. That is the signature of
family-specific functional constraint rather than fold constraint: the fold
should tolerate losing them, but they may belong to the switch machinery.
Testing exactly that is why script 07 exists.

Positions deliberately EXCLUDED: gate (S85-A89), latch (H115-L117 backbone),
and the sealing residues F61 / F159 / L87, which gave NEGATIVE rigid-truncation
volumes because they breach the cavity to bulk solvent rather than enlarging it.
"""

WT = {59: "K", 79: "R", 94: "E", 108: "F", 120: "Y", 141: "E"}
SCAN_POSITIONS = [59, 79, 94, 108, 120, 141]
SCAN_SUBS = ["A", "G", "S", "V", "L"]


def build_panel():
    """Return list of (name, [(resnum, one_letter_aa), ...])."""
    panel = [("WT", [])]

    # 1. single-substitution scan: 6 positions x 5 substitutions = 30
    for p in SCAN_POSITIONS:
        for aa in SCAN_SUBS:
            if aa == WT[p]:
                continue
            panel.append((f"{WT[p]}{p}{aa}", [(p, aa)]))

    # 2. superfamily-guided singles
    panel += [
        ("F108Y", [(108, "Y")]),   # superfamily consensus is LARGER than WT
        ("K59D", [(59, "D")]),     # superfamily consensus residue at 59
        ("E141W", [(141, "W")]),   # distant-subfamily signature
    ]

    # 3. wall-breaks: test whether F108 needs its R79/E94 buttress removed
    panel += [
        ("R79A_E94A",             [(79, "A"), (94, "A")]),
        ("F108A_R79A",            [(108, "A"), (79, "A")]),
        ("K59A_F108A",            [(59, "A"), (108, "A")]),
        ("F108A_R79A_E94A",       [(108, "A"), (79, "A"), (94, "A")]),
        ("K59A_F108A_R79A_E94A",  [(59, "A"), (108, "A"), (79, "A"), (94, "A")]),
    ]

    # 4. superfamily-signature combination (hydrophobic-59 / Y108)
    panel += [
        ("K59L_F108Y", [(59, "L"), (108, "Y")]),
        ("K59L_F108Y_E141W", [(59, "L"), (108, "Y"), (141, "W")]),
    ]
    return panel


PANEL = build_panel()

if __name__ == "__main__":
    for i, (n, muts) in enumerate(PANEL):
        print(f"{i:3d}  {n:<24} {muts}")
    print(f"\n{len(PANEL)} variants (including WT)")
