# 61_pymol_loop_view.pml -- load a PYR1 MD trajectory and show the gate/latch
# with a flat, silhouetted look (the PyMOL equivalent of ChimeraX "simple"
# lighting + silhouettes).
#
# USAGE (local machine):
#   pymol -r 61_pymol_loop_view.pml
# or from inside PyMOL:
#   @/path/to/61_pymol_loop_view.pml
#
# EDIT THESE TWO PATHS FIRST.
# ---------------------------------------------------------------------------
# NUMBERING WARNING -- READ THIS
#   tleap renumbers residues sequentially 1..178 in the prmtop, and the crystal
#   has gaps, so the gate (native PYR1 85-89) is `resi 82-86` HERE, and the latch
#   (native 115-117) is `resi 112-114`. `resi 85-89` selects the WRONG residues
#   and looks perfectly reasonable on screen.
#
#   So the selections below are made by SEQUENCE, not by number: SGLPA (gate) and
#   HRL (latch) are each unique in this construct, which makes the selection
#   immune to the numbering entirely. The script prints what it picked -- check
#   that it says SER GLY LEU PRO ALA before believing any picture.
# ---------------------------------------------------------------------------

# ============ 1. paths ============
# Use forward slashes even on Windows.
set_wd ~/Downloads

# ============ 2. load topology + trajectory ============
# THESE ARE THE `view` FILES from 41_md_view.sh, i.e. data/md/view/*.prmtop|.nc
# They are ALREADY processed and that changes what you must not do again:
#     * waters/Na+/Cl- stripped   -> 2,839 atoms, not 46,820
#     * autoimage applied         -> no PBC teleporting; do NOT re-image
#     * CA-fit to frame 1         -> no global tumbling; do NOT intra_fit
#     * STRIDED BY 10             -> 3,000 frames at 100 ps spacing = 300 ns
#
# So do NOT pass `interval` here. Striding an already-strided file by 100 leaves
# 30 frames, which is what happens if you treat these as the raw prod.nc.
load S1_apo_open_rep0.prmtop, S1
load_traj S1_apo_open_rep0.nc, S1

# 3,000 frames x 2,839 atoms is comfortable. If playback is sluggish, stride the
# LOAD (not the file) with e.g.  interval=5  -> 600 frames at 500 ps.

# ============ 3. nothing to strip ============
# Deliberately empty. The view files carry no solvent, so `remove not polymer`
# is at best a no-op and at worst destructive: in the S2/S4 view files it would
# delete ABA and the Mn2+, which 41_md_view.sh kept on purpose.
remove hydrogens                       # optional, purely for display speed

# ============ 4. selections ============
# Verified directly against THIS topology with cpptraj `resinfo`:
#     resi 82-86  = SER GLY LEU PRO ALA  (gate,  native 85-89)
#     resi 112-114 = HIE ARG LEU         (latch, native 115-117)
# Native numbering is NOT usable here -- `resi 85-89` silently selects the wrong
# residues.
#
# The gate is also selectable by sequence, which is numbering-proof (SGLPA is
# unique in this construct). The LATCH is not: its histidine is named HIE by
# tleap, and `pepseq` may not map HIE to H, so the latch uses the verified
# numbers. Trust the printout below, not either method.
select gate,  polymer and pepseq SGLPA
select latch, polymer and resi 112-114
deselect

# Print what was actually picked. Gate MUST read SER GLY LEU PRO ALA,
# latch MUST read HIE/HIS ARG LEU. If either is wrong, stop and fix it.
/print("---- gate residues selected ----")
iterate gate and name CA, print(" gate  resi %4s  %s" % (resi, resn))
/print("---- latch residues selected ----")
iterate latch and name CA, print(" latch resi %4s  %s" % (resi, resn))

# ============ 5. representation ============
hide everything
show cartoon, polymer
color grey80, polymer

color tv_orange, gate
color slate,     latch
show sticks, (gate or latch) and sidechain
set cartoon_side_chain_helper, 1

# make the two loops read as loops, not tubes lost in the fold
set cartoon_loop_radius, 0.25
set cartoon_transparency, 0.0

# ============ 6. flat + silhouette look ============
# ChimeraX "simple"/"flat" lighting is ambient-dominated with no specular
# highlights and no shadows; silhouettes are a black outline pass.
bg_color white
set ray_opaque_background, 1

# --- flat lighting ---
set ambient, 0.45
set direct, 0.55
set reflect, 0.05
set specular, 0
set spec_reflect, 0
set shininess, 10
set light_count, 2
set ray_shadows, 0
set depth_cue, 0
set two_sided_lighting, 1

# --- silhouettes ---
# ray_trace_mode 1 = colour + black outline (closest to ChimeraX silhouettes)
#                 2 = outline only, white fill
#                 3 = cel/poster shaded, the flattest of the three
set ray_trace_mode, 1
set ray_trace_color, black
set ray_trace_gain, 0.15
set antialias, 2

# nicer cartoon geometry
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_flat_sheets, 1

# ============ 7. view + render ============
orient polymer
turn x, -10

/print("")
/print("NOTE: silhouettes only appear in RAY-TRACED output, not the live window.")
/print("      Render with:   ray 1600, 1200      then:   png gate_view.png, dpi=300")
/print("      Play the trajectory with:   mplay      (stop with mstop)")
/print("")
/print("For a fully flat, cel-shaded look instead:  set ray_trace_mode, 3")
/print("Thicker outlines:                           set ray_trace_gain, 0.3")
