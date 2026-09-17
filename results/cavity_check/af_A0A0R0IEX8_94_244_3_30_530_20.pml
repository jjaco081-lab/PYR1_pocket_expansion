# af_A0A0R0IEX8_94_244_3.30.530.20: v2 318.6 -> v4 0.0; a 277.6 A^3 chamber exists but should NOT be in the domain
load af_A0A0R0IEX8_94_244_3_30_530_20_protein.pdb, prot
load af_A0A0R0IEX8_94_244_3_30_530_20_cavity.pdb, cav
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
print "  A    277.6 A^3   domain-lined fraction 0.00"
print "  B     57.2 A^3   domain-lined fraction 0.00"
print "  C     32.8 A^3   domain-lined fraction 0.00"
print "  D     30.7 A^3   domain-lined fraction 0.00"
print "  E     16.6 A^3   domain-lined fraction 0.00"
print "  F     12.7 A^3   domain-lined fraction 0.00"
print "prot chain A = CATH domain (cartoon); chain B = rest of the file (lines)"
