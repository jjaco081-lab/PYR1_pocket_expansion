# af_A0A0R0IEX8_94_244_3.30.530.20: v2 318.6 -> v4 0.0; a 277.6 A^3 chamber exists but should NOT be in the domain
load af_A0A0R0IEX8_94_244_3_30_530_20_source.pdb, prot
load af_A0A0R0IEX8_94_244_3_30_530_20_cavity.pdb, cav
hide everything
select domain, prot and polymer and (resi 94-244)
select other, prot and polymer and not domain
select ligand, prot and not polymer and not resn HOH
show cartoon, domain
color grey80, domain
show lines, other
color palecyan, other
show sticks, ligand
color yellow, ligand
set stick_radius, 0.28, ligand
util.cnc ligand
show spheres, cav
set sphere_scale, 0.20, cav
set sphere_transparency, 0.55, cav
color firebrick, cav and chain A
color grey50, cav and chain B
color grey50, cav and chain C
color grey50, cav and chain D
color grey50, cav and chain E
color grey50, cav and chain F
color grey50, cav and chain G
color grey50, cav and chain H
set cartoon_transparency, 0.55
bg_color white
deselect
orient cav and chain A
zoom cav and chain A, 7
print "CHAMBERS (chain = rank, red = largest):"
print "  A    277.6 A^3   domain-lined fraction 0.00"
print "  B     57.2 A^3   domain-lined fraction 0.00"
print "  C     32.8 A^3   domain-lined fraction 0.00"
print "  D     30.7 A^3   domain-lined fraction 0.00"
print "  E     16.6 A^3   domain-lined fraction 0.00"
print "  F     12.7 A^3   domain-lined fraction 0.00"
print "grey cartoon = CATH domain; pale lines = rest of the file; YELLOW STICKS = ligand; red spheres = reported chamber"
print "the check: does the RED blob sit on the yellow ligand, inside the grey cartoon?"
