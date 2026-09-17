# af_Q6Z9J1_229_401_3.30.530.20: v2 216.6 -> v4 9.3, domain fraction 0.76
load af_Q6Z9J1_229_401_3_30_530_20_protein.pdb, prot
load af_Q6Z9J1_229_401_3_30_530_20_cavity.pdb, cav
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
print "  A    213.6 A^3   domain-lined fraction 0.53"
print "  B    186.2 A^3   domain-lined fraction 0.31"
print "  C     77.3 A^3   domain-lined fraction 0.00"
print "  D     56.8 A^3   domain-lined fraction 0.00"
print "  E      9.3 A^3   domain-lined fraction 0.76"
print "  F      9.3 A^3   domain-lined fraction 0.71"
print "  G      5.6 A^3   domain-lined fraction 0.76"
print "  H      3.7 A^3   domain-lined fraction 1.00"
print "prot chain A = CATH domain (cartoon); chain B = rest of the file (lines)"
