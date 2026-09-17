# 3oquB00
# overlap 0.85, EXPERIMENTAL: fpocket 1489.0 vs ours 644.5, ratio 2.31 -- BEST AGREEMENT, closest to PYR1's 2.18. Positive control.
#
# TWO cavity representations, loaded on top of each other so they
# can be compared directly:
#   BLUE TRANSPARENT SURFACE = PyMOL's own cavity detection
#        (surface_cavity_mode 2, radius -3, cutoff -5)
#   SOLID SPHERES = the chamber THIS PROJECT measures
# Where they disagree is the interesting part.
# F1 toggles the surface, F2 toggles the spheres.
load 3oquB00_source.cif, prot
load 3oquB00_cavity.pdb, cav
hide everything
select domain, prot and polymer and (all)
select other, prot and polymer and not domain
select ligand, prot and not polymer and not resn HOH
show cartoon, domain
color grey70, domain
show lines, other
color palecyan, other
show sticks, ligand
color yellow, ligand
util.cnc ligand
set stick_radius, 0.28, ligand

# --- PyMOL's own cavity detection, Jannis's settings ---
set surface_cavity_mode, 2
set surface_cavity_radius, -3
set surface_cavity_cutoff, -5
set surface_quality, 1
create pmolcav, prot and polymer
hide everything, pmolcav
show mesh, pmolcav
color skyblue, pmolcav
set mesh_width, 0.4, pmolcav
set transparency_mode, 3
set transparency, 0.65, pmolcav

# --- the chambers measured here: FULLY OPAQUE ---
show spheres, cav
set sphere_scale, 0.16, cav
set sphere_transparency, 0.0, cav
color firebrick, cav and chain A
color orange, cav and chain B
color yellow, cav and chain C
color green, cav and chain D
color purple, cav and chain E
color magenta, cav and chain F
color salmon, cav and chain G
color wheat, cav and chain H
set cartoon_transparency, 0.7
bg_color white
deselect
orient cav and chain A
zoom cav and chain A, 7

# --- views -------------------------------------------
python
from pymol import cmd
def v_mesh():
    cmd.hide('surface','pmolcav'); cmd.show('mesh','pmolcav')
def v_surf():
    cmd.hide('mesh','pmolcav'); cmd.show('surface','pmolcav')
def v_onlypmol():
    cmd.disable('cav'); cmd.enable('pmolcav')
def v_onlyours():
    cmd.enable('cav'); cmd.disable('pmolcav')
def v_both():
    cmd.enable('cav'); cmd.enable('pmolcav')
cmd.extend('v_mesh',v_mesh); cmd.extend('v_surf',v_surf)
cmd.extend('v_onlypmol',v_onlypmol)
cmd.extend('v_onlyours',v_onlyours)
cmd.extend('v_both',v_both)
python end
set_key F1, v_mesh
set_key F2, v_surf
set_key F3, v_onlyours
set_key F4, v_onlypmol
set_key F5, v_both
print "CHAMBERS (chain = rank, red = largest):"
print "  A    644.5 A^3   domain-lined fraction 1.00"
print "grey cartoon = CATH domain; pale lines = rest of the file; YELLOW STICKS = ligand; red spheres = reported chamber"
print "the check: does the RED blob sit on the yellow ligand, inside the grey cartoon?"
