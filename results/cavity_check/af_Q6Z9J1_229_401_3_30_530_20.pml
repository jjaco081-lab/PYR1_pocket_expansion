# af_Q6Z9J1_229_401_3.30.530.20: v2 216.6 -> v4 9.3, domain fraction 0.76
load af_Q6Z9J1_229_401_3_30_530_20_source.pdb, prot
load af_Q6Z9J1_229_401_3_30_530_20_cavity.pdb, cav
hide everything
select domain, prot and polymer and (resi 229-401)
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
print "  A    213.6 A^3   domain-lined fraction 0.53"
print "  B    186.2 A^3   domain-lined fraction 0.31"
print "  C     77.3 A^3   domain-lined fraction 0.00"
print "  D     56.8 A^3   domain-lined fraction 0.00"
print "  E      9.3 A^3   domain-lined fraction 0.76"
print "  F      9.3 A^3   domain-lined fraction 0.71"
print "  G      5.6 A^3   domain-lined fraction 0.76"
print "  H      3.7 A^3   domain-lined fraction 1.00"
print "grey cartoon = CATH domain; pale lines = rest of the file; YELLOW STICKS = ligand; red spheres = reported chamber"
print "the check: does the RED blob sit on the yellow ligand, inside the grey cartoon?"
