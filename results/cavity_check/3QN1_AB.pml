# 3QN1_AB: PYR1+HAB1: 182.3 A^3 as ONE chamber -- HAB1 seals it
load 3QN1_AB_source.pdb, prot
load 3QN1_AB_cavity.pdb, cav
hide everything
select domain, prot and polymer and (all)
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
print "  A    182.3 A^3   domain-lined fraction 1.00"
print "grey cartoon = CATH domain; pale lines = rest of the file; YELLOW STICKS = ligand; red spheres = reported chamber"
print "the check: does the RED blob sit on the yellow ligand, inside the grey cartoon?"
