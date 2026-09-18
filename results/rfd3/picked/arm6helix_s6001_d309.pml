# arm6helix_s6001/d309
# cavity 180.4 A^3 (1.10x PYR1), SS 0.88, orange 5 total / 1 longest run
load arm6helix_s6001_d309.pdb, arm6helix_s6001_d309
hide everything
show cartoon, arm6helix_s6001_d309 and chain A
spectrum b, grey70 orange red, arm6helix_s6001_d309 and chain A, 0, 2
show cartoon, arm6helix_s6001_d309 and not chain A
color palecyan, arm6helix_s6001_d309 and not chain A
set cartoon_transparency, 0.55, arm6helix_s6001_d309 and not chain A
show sticks, arm6helix_s6001_d309 and not polymer
color yellow, arm6helix_s6001_d309 and not polymer
bg_color white
orient arm6helix_s6001_d309 and chain A
print "GREY = core   ORANGE = kept but not core   RED = truncated"
print "pale cyan = HAB1 context (four chains)   yellow = ligand if present"
