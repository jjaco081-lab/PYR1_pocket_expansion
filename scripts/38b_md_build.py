#!/usr/bin/env python
"""
38b_md_build.py -- parameterise ABA and build the four solvated MD systems.

Supersedes the placeholder tleap inputs written by 38_md_prep.py: it generates
them here so the NaCl count can be set from the ACTUAL solvated water count
(two-pass), rather than guessed.

ABA parameters
--------------
data/aba_deprot.sdf is used directly as the antechamber input. It is verified to
be in the crystal frame (0.000 A per-atom deviation from aba_xtal.pdb, same atom
order), carries explicit hydrogens, and has formal charge -1 on the carboxylate,
matching the protonation decision made in scripts/00_setup_inputs.py. Because the
coordinates are already correct, the ligand is loaded with loadmol2 and needs no
separate placement step.

    antechamber -c bcc -nc -1 -at gaff2      AM1-BCC charges, GAFF2 atom types
    parmchk2    -s gaff2                     missing parameters -> frcmod

Force field: ff19SB protein / OPC water (the pairing ff19SB was fitted for) /
GAFF2 ligand / Li-Merz 12-6 ions for OPC. Mn2+ in HAB1 (3QN1) is retained.

Two passes per system:
    pass 1  solvate, count waters, compute ions for 0.15 M
    pass 2  rebuild with addionsrand, verify net charge == 0

Usage:  python 38b_md_build.py [--systems S1_apo_open,...] [--pad 12.0]
        Run on a compute node or the login node (antechamber on 38 atoms is small).
"""
import argparse, os, re, subprocess, sys, glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
MD = os.path.join(DATA, "md")
AMBER_MOD = "amber/22"

ap = argparse.ArgumentParser()
ap.add_argument("--systems", default="")
ap.add_argument("--pad", type=float, default=12.0)
ap.add_argument("--conc", type=float, default=0.15, help="NaCl molarity")
args = ap.parse_args()

LIGAND_SYSTEMS = {"S2_holo_closed", "S4_ternary"}
ALL = sorted(d for d in os.listdir(MD) if os.path.isdir(os.path.join(MD, d)))
SYSTEMS = [s.strip() for s in args.systems.split(",") if s.strip()] or ALL


def run(cmd, cwd, log):
    """run under the amber module; return (rc, output)."""
    full = f"module load {AMBER_MOD} >/dev/null 2>&1; {cmd}"
    r = subprocess.run(["bash", "-lc", full], cwd=cwd, capture_output=True, text=True)
    with open(os.path.join(cwd, log), "w") as fh:
        fh.write(f"$ {cmd}\n\n{r.stdout}\n{r.stderr}")
    return r.returncode, r.stdout + r.stderr


# ---------------- 1. ABA parameters (once) ----------------
mol2 = os.path.join(MD, "A8S.mol2")
frcmod = os.path.join(MD, "A8S.frcmod")
lib = os.path.join(MD, "A8S.lib")

if not (os.path.exists(mol2) and os.path.exists(frcmod) and os.path.exists(lib)):
    print("=" * 74)
    print("PARAMETERISING ABA (antechamber AM1-BCC / GAFF2, net charge -1)")
    print("=" * 74)
    sdf = os.path.join(DATA, "aba_deprot.sdf")
    rc, out = run(f"antechamber -i {sdf} -fi sdf -o A8S.mol2 -fo mol2 "
                  f"-c bcc -nc -1 -at gaff2 -rn A8S -s 2 -dr no",
                  MD, "antechamber.log")
    if rc != 0 or not os.path.exists(mol2):
        sys.exit(f"antechamber failed (rc={rc}); see {MD}/antechamber.log\n"
                 + out[-1500:])
    print("  antechamber OK ->", mol2)

    rc, out = run("parmchk2 -i A8S.mol2 -f mol2 -o A8S.frcmod -s gaff2",
                  MD, "parmchk2.log")
    if rc != 0:
        sys.exit(f"parmchk2 failed; see {MD}/parmchk2.log")
    # flag any parameter guessed with a high penalty
    bad = [l for l in open(frcmod) if "ATTN" in l]
    print(f"  parmchk2 OK -> {frcmod}"
          + (f"   !! {len(bad)} parameters flagged ATTN (check these)" if bad else
             "   (no ATTN flags)"))
    for l in bad[:6]:
        print("      " + l.rstrip())

    open(os.path.join(MD, "mklib.in"), "w").write(
        "source leaprc.gaff2\n"
        "loadamberparams A8S.frcmod\n"
        "A8S = loadmol2 A8S.mol2\n"
        "check A8S\n"
        "saveoff A8S A8S.lib\n"
        "quit\n")
    rc, out = run("tleap -f mklib.in", MD, "mklib.log")
    if rc != 0 or not os.path.exists(lib):
        sys.exit(f"tleap lib build failed; see {MD}/mklib.log")
    q = re.search(r"Total unperturbed charge:\s*(-?[\d.]+)", out)
    print(f"  A8S.lib OK   net charge = {q.group(1) if q else '?'} (expect -1.00)")
else:
    print("ABA parameters already present, reusing A8S.mol2 / .frcmod / .lib")


# NOTE: leaprc.water.opc already loads frcmod.ionslm_126_opc, which is the
# OPC-specific Li/Merz 12-6 set and covers mono-, di- and trivalent ions
# (Mn2+ included, Li et al. JCTC 2020, 16, 4429). Do NOT add
# frcmod.ions234lm_126_opc -- that is the TIP3P-era naming and no such file
# exists for OPC; loading it makes tleap exit with "Could not open file".
TLEAP = """source leaprc.protein.ff19SB
source leaprc.water.opc
source leaprc.gaff2
{ligload}
prot = loadpdb protein.pdb
{combine}
solvateoct sys OPCBOX {pad}
{ions}
charge sys
check sys
saveamberparm sys system.prmtop system.inpcrd
savepdb sys system_solvated.pdb
quit
"""


def build(sysname, nna=0, ncl=0, tag="pass1"):
    d = os.path.join(MD, sysname)
    lig = sysname in LIGAND_SYSTEMS
    ligload = ("loadamberparams ../A8S.frcmod\nloadoff ../A8S.lib"
               "\nlig = loadmol2 ../A8S.mol2") if lig else ""
    combine = "sys = combine { prot lig }" if lig else "sys = prot"
    ions = "addions sys Na+ 0\naddions sys Cl- 0"
    if nna or ncl:
        ions += f"\naddionsrand sys Na+ {nna} Cl- {ncl}"
    open(os.path.join(d, "tleap.in"), "w").write(TLEAP.format(
        ligload=ligload, combine=combine, pad=args.pad, ions=ions))
    rc, out = run("tleap -f tleap.in", d, f"tleap_{tag}.log")
    return rc, out


print("\n" + "=" * 74)
print("BUILDING SOLVATED SYSTEMS")
print("=" * 74)
print(f"{'system':<17}{'waters':>8}{'Na+':>6}{'Cl-':>6}{'atoms':>9}"
      f"{'charge':>9}  status")

for s in SYSTEMS:
    d = os.path.join(MD, s)
    if not os.path.exists(os.path.join(d, "protein.pdb")):
        print(f"{s:<17}   protein.pdb missing -- run 38_md_prep.py first")
        continue

    rc, out = build(s, tag="pass1")
    nwat = len(re.findall(r"\bWAT\b", out)) or None
    m = re.search(r"Added\s+(\d+)\s+residues", out)
    if m:
        nwat = int(m.group(1))
    if nwat is None:
        print(f"{s:<17}   pass1 failed -- see {d}/tleap_pass1.log")
        continue

    # 0.15 M NaCl: ions per water = conc / 55.5
    nion = int(round(nwat * args.conc / 55.5))
    rc, out = build(s, nna=nion, ncl=nion, tag="pass2")

    q = re.search(r"Total unperturbed charge:\s*(-?[\d.]+)", out)
    qq = re.findall(r"Total perturbed charge:\s*(-?[\d.]+)", out)
    natom = None
    pdb = os.path.join(d, "system_solvated.pdb")
    if os.path.exists(pdb):
        natom = sum(1 for l in open(pdb) if l.startswith(("ATOM", "HETATM")))
    ok = os.path.exists(os.path.join(d, "system.prmtop"))
    charge = q.group(1) if q else "?"
    flag = "OK" if ok and abs(float(charge if charge != '?' else 9)) < 0.01 else \
           ("built, CHECK CHARGE" if ok else "FAILED")
    print(f"{s:<17}{nwat:>8}{nion:>6}{nion:>6}{natom or 0:>9}{charge:>9}  {flag}")

print(f"""
Inputs are in {MD}/<system>/system.prmtop + system.inpcrd

Verify before running:
  * net charge must be 0.00 for every system
  * check tleap_pass2.log for 'Could not find bond parameter' or missing atoms
  * S4_ternary keeps Mn2+; confirm tleap typed it (search the log for MN)
NEXT: scripts/39_md_run.sh  (partition gpu, --gres=gpu:a100:1, pmemd.cuda)""")
