# 2pcsA00: CoxG donor 570.5 A^3 -- is it really one chamber?
load 2pcsA00_protein.pdb, prot
load 2pcsA00_cavity.pdb, cav
hide everything
show cartoon, prot and chain A
color grey80, prot and chain A
show lines, prot and chain B
color palecyan, prot and chain B
show spheres, cav
set sphere_scale, 0.22, cav
color firebrick, cav and chain A
color grey50, cav and chain B
color grey50, cav and chain C
color grey50, cav and chain D
color grey50, cav and chain E
color grey50, cav and chain F
color grey50, cav and chain G
color grey50, cav and chain H
set transparency, 0.6
bg_color white
zoom cav and chain A, 6
print "CHAMBERS (chain = rank, red = largest):"
print "  A    570.5 A^3   domain-lined fraction 1.00"
print "prot chain A = CATH domain (cartoon); chain B = rest of the file (lines)"
