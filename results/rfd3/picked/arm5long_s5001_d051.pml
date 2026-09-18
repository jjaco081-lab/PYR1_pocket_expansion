# arm5long_s5001/d051
# cavity 176.5 A^3 (1.07x PYR1), SS 0.86, orange 5 total / 3 longest run
load arm5long_s5001_d051.pdb, arm5long_s5001_d051
hide everything
show cartoon, arm5long_s5001_d051 and chain A
spectrum b, grey70 orange red, arm5long_s5001_d051 and chain A, 0, 2
show cartoon, arm5long_s5001_d051 and not chain A
color palecyan, arm5long_s5001_d051 and not chain A
set cartoon_transparency, 0.55, arm5long_s5001_d051 and not chain A
show sticks, arm5long_s5001_d051 and not polymer
color yellow, arm5long_s5001_d051 and not polymer
bg_color white
orient arm5long_s5001_d051 and chain A
print "GREY = core   ORANGE = kept but not core   RED = truncated"
print "pale cyan = HAB1 context (four chains)   yellow = ligand if present"
