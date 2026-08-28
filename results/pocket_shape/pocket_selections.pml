# ---- PYR1 / 3QN1 (chain A) ----
fetch 3qn1, async=0
hide everything, 3qn1
show cartoon, 3qn1 and chain A
select pyr1_lig, 3qn1 and chain A and resn A8S
select pyr1_wall, 3qn1 and chain A and resi 59+60+61+62+81+83+87+88+89+91+92+94+108+110+115+117+120+122+141+159+160+163+164+167
select pyr1_wall_contact, 3qn1 and chain A and resi 59+61+83+87+88+89+91+92+94+108+110+115+117+120+141+159+163+164+167
show sticks, pyr1_wall or pyr1_lig
color grey80, pyr1_wall
color salmon, pyr1_wall_contact
color yellow, pyr1_lig

# ---- 3OQU (chain B) ----
fetch 3oqu, async=0
hide everything, 3oqu
show cartoon, 3oqu and chain B
select oqu_lig, 3oqu and chain B and resn A8S
select oqu_wall, 3oqu and chain B and resi 52+53+55+56+59+62+63+65+66+81+83+85+89+91+92+93+94+96+108+109+110+112+117+119+122+124+126+143+162+165+166+169
select oqu_wall_contact, 3oqu and chain B and resi 63+65+85+89+91+94+96+112+117+119+122+143+162+165+166+169
show sticks, oqu_wall or oqu_lig
color grey80, oqu_wall
color skyblue, oqu_wall_contact
color yellow, oqu_lig

# ---- superpose and compare ----
align 3oqu and chain B, 3qn1 and chain A
# cavity surfaces:
set surface_cavity_mode, 1
set surface_cavity_radius, 5
show surface, 3qn1 and chain A
show surface, 3oqu and chain B
