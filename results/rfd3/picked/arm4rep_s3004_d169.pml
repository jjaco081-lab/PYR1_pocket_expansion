# arm4rep_s3004/d169
# cavity 294.8 A^3 (1.79x PYR1), SS 0.89, orange 11 total / 5 longest run
load arm4rep_s3004_d169.pdb, arm4rep_s3004_d169
hide everything
show cartoon, arm4rep_s3004_d169 and chain A
spectrum b, grey70 orange red, arm4rep_s3004_d169 and chain A, 0, 2
show cartoon, arm4rep_s3004_d169 and not chain A
color palecyan, arm4rep_s3004_d169 and not chain A
set cartoon_transparency, 0.55, arm4rep_s3004_d169 and not chain A
show sticks, arm4rep_s3004_d169 and not polymer
color yellow, arm4rep_s3004_d169 and not polymer
bg_color white
orient arm4rep_s3004_d169 and chain A
print "GREY = core   ORANGE = kept but not core   RED = truncated"
print "pale cyan = HAB1 context (four chains)   yellow = ligand if present"
