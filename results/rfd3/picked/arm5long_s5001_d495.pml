# arm5long_s5001/d495
# cavity 230.9 A^3 (1.40x PYR1), SS 0.84, orange 7 total / 2 longest run
load arm5long_s5001_d495.pdb, arm5long_s5001_d495
hide everything
show cartoon, arm5long_s5001_d495 and chain A
spectrum b, grey70 orange red, arm5long_s5001_d495 and chain A, 0, 2
show cartoon, arm5long_s5001_d495 and not chain A
color palecyan, arm5long_s5001_d495 and not chain A
set cartoon_transparency, 0.55, arm5long_s5001_d495 and not chain A
show sticks, arm5long_s5001_d495 and not polymer
color yellow, arm5long_s5001_d495 and not polymer
bg_color white
orient arm5long_s5001_d495 and chain A
print "GREY = core   ORANGE = kept but not core   RED = truncated"
print "pale cyan = HAB1 context (four chains)   yellow = ligand if present"
