#!/bin/bash
# lib_mdinputs.sh -- THE canonical Amber input decks for this project.
#
# WHY THIS EXISTS
#   Four runners (39, 42, 57, 65) each carried their own copy of the same six input
#   decks. They happened to agree on every physics setting, but only because they
#   were copied carefully -- there was nothing stopping them drifting, and a drifted
#   `cut` or `gamma_ln` would not be visible in any output. The decks live here now.
#
# CANONICAL PROTOCOL (README 19c)
#   min1  5000 steps, solute heavy atoms restrained 10 kcal/mol/A^2
#   min2  5000 steps, unrestrained
#   heat  NVT 0->300 K over 200 ps, restraints 5.0
#   eq1   NPT 300 K, 500 ps, restraints 1.0
#   eq2   NPT 300 K, 1 ns, unrestrained
#   prod  NPT 300 K, 4 fs (HMR), SHAKE on H, PME, 10 A cutoff, Langevin gamma_ln=2.0,
#         Monte Carlo barostat, frames every 10 ps
#
#   The clock includes equilibration: 200 + 500 + 1000 = 1700 ps, so a finished
#   300 ns replicate reads 301700 ps, not 300000.
#
# ntwr: RESTART WRITE INTERVAL, and the one place this deliberately differs from
#   the older runners. 39/42/57 used ntwr=500000 (1 ns) while frames are written
#   every 10 ps, so an interrupted run re-simulates from the last restart and the
#   trajectory keeps the frames written past it -- S4 accumulated 11.6 ns of
#   duplicated frames that way (README 24e). Setting ntwr = ntwx makes restart and
#   frame writes coincide, so an interruption costs at most one frame. This changes
#   I/O only; it does not touch the trajectory or the physics, so systems built
#   before and after remain directly comparable.
#
# OVERRIDING
#   Constants read from the environment and any deviation is printed loudly.
#
# USAGE
#   source "$(dirname "$0")/lib_mdinputs.sh"
#   write_md_inputs                 # min1/min2/heat/eq1/eq2 into $PWD
#   write_prod_input <n_steps>      # prod.in for a given remaining step count

# HYDROGEN-MASS REPARTITIONING -- 4 fs production
#   lib_solvate.sh writes system_hmr.prmtop alongside system.prmtop. With mass moved
#   into the hydrogens, the fastest bond vibrations slow enough that 4 fs is stable,
#   roughly DOUBLING throughput. Water is excluded from the repartitioning (already
#   rigid under SHAKE).
#
#   WHAT IT COSTS, stated honestly: HMR is well validated for EQUILIBRIUM and
#   conformational properties -- RMSD, RMSF, populations, free energies are
#   unaffected. It is NOT kinetically faithful: the masses are deliberately wrong,
#   so diffusion and relaxation times are altered. Every observable pre-registered
#   in README 19d is structural/equilibrium, so this is a fair trade. The one number
#   it touches is the integrated autocorrelation time tau, which is used only to
#   compute effective sample size; applied consistently to every system, tau stays
#   comparable BETWEEN systems even though it is not a physical relaxation rate.
#   That was already true under Langevin gamma=2.0 (README 24e).
#
#   HEATING AND EQ1 STAY AT 2 fs. Heating from 0 K under positional restraints is
#   the least stable phase of the protocol and it costs only 700 ps to be careful.
DT_EQUIL=${DT_EQUIL:-0.002}
DT=${DT:-0.004}
STEPS_PER_PS=$(python3 -c "print(int(round(1.0/${DT})))")
CUT=${CUT:-10.0}
GAMMA_LN=${GAMMA_LN:-2.0}
TEMP0=${TEMP0:-300.0}
BAROSTAT=${BAROSTAT:-2}
TAUP=${TAUP:-2.0}
# ntwx is in STEPS, so it must be derived from dt or the frame spacing silently
# changes: 5000 steps is 10 ps at 2 fs but 20 ps at 4 fs, and script 59 assumes
# 10 ps/frame. Keep the SPACING fixed, not the step count.
FRAME_PS=${FRAME_PS:-10}
NTWX=${NTWX:-$((FRAME_PS * STEPS_PER_PS))}
NTPR_EQ=${NTPR_EQ:-5000}
NTPR_PROD=${NTPR_PROD:-25000}
NTWR=${NTWR:-$NTWX}                # see header: restart cadence == frame cadence
# The mask must exclude WATER AND WHICHEVER IONS WERE ADDED. It previously named
# Na+ only; after the switch to KCl that would have left every K+ positionally
# restrained through minimisation, heating and eq1. Both are listed so the mask is
# correct whichever salt a system was built with.
SOLUTE_MASK=${SOLUTE_MASK:-'!:WAT,Na+,K+,Cl- & !@H='}
RW_MIN=${RW_MIN:-10.0}
RW_HEAT=${RW_HEAT:-5.0}
RW_EQ1=${RW_EQ1:-1.0}
NSTLIM_HEAT=${NSTLIM_HEAT:-100000}   # 200 ps
NSTLIM_EQ1=${NSTLIM_EQ1:-250000}     # 500 ps
NSTLIM_EQ2=${NSTLIM_EQ2:-$((1000 * STEPS_PER_PS))}   # 1 ns at the production dt
EQUIL_PS=${EQUIL_PS:-1700}

md_settings_check () {
    local n=0
    [[ "$DT" == "0.004"   ]] || { echo "  !! DT overridden: $DT"; n=1; }
    [[ "$DT_EQUIL" == "0.002" ]] || { echo "  !! DT_EQUIL overridden: $DT_EQUIL"; n=1; }
    [[ "$CUT" == "10.0"   ]] || { echo "  !! CUT overridden: $CUT"; n=1; }
    [[ "$GAMMA_LN" == "2.0" ]] || { echo "  !! GAMMA_LN overridden: $GAMMA_LN"; n=1; }
    [[ "$TEMP0" == "300.0" ]] || { echo "  !! TEMP0 overridden: $TEMP0"; n=1; }
    [[ "$FRAME_PS" == "10" ]] || { echo "  !! FRAME_PS overridden: $FRAME_PS"; n=1; }
    [[ "$EQUIL_PS" == "1700" ]] || { echo "  !! EQUIL_PS overridden: $EQUIL_PS"; n=1; }
    [[ $n -eq 1 ]] && echo "  !! NON-CANONICAL MD SETTINGS -- not comparable with the rest"
    return 0
}

write_md_inputs () {
cat > min1.in <<EOF
minimisation 1, solute restrained
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=1, cut=${CUT},
  ntr=1, restraintmask='${SOLUTE_MASK}', restraint_wt=${RW_MIN},
 /
EOF
cat > min2.in <<EOF
minimisation 2, unrestrained
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=1, cut=${CUT}, ntr=0,
 /
EOF
cat > heat.in <<EOF
NVT heating 0 -> ${TEMP0} K over 200 ps
 &cntrl
  imin=0, irest=0, ntx=1,
  nstlim=${NSTLIM_HEAT}, dt=${DT_EQUIL},
  ntc=2, ntf=2, cut=${CUT},
  ntb=1, ntp=0,
  ntt=3, gamma_ln=${GAMMA_LN}, ig=-1,
  tempi=0.0, temp0=${TEMP0},
  nmropt=1,
  ntr=1, restraintmask='${SOLUTE_MASK}', restraint_wt=${RW_HEAT},
  ntpr=${NTPR_EQ}, ntwx=${NTWX}, ntwr=50000,
 /
 &wt type='TEMP0', istep1=0, istep2=80000, value1=0.0, value2=${TEMP0} /
 &wt type='TEMP0', istep1=80001, istep2=${NSTLIM_HEAT}, value1=${TEMP0}, value2=${TEMP0} /
 &wt type='END' /
EOF
cat > eq1.in <<EOF
NPT equilibration 1, 500 ps, light restraints
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=${NSTLIM_EQ1}, dt=${DT_EQUIL},
  ntc=2, ntf=2, cut=${CUT},
  ntb=2, ntp=1, barostat=${BAROSTAT}, pres0=1.0, taup=${TAUP},
  ntt=3, gamma_ln=${GAMMA_LN}, ig=-1, temp0=${TEMP0},
  ntr=1, restraintmask='${SOLUTE_MASK}', restraint_wt=${RW_EQ1},
  ntpr=${NTPR_EQ}, ntwx=${NTWX}, ntwr=50000,
 /
EOF
cat > eq2.in <<EOF
NPT equilibration 2, 1 ns, unrestrained
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=${NSTLIM_EQ2}, dt=${DT},
  ntc=2, ntf=2, cut=${CUT},
  ntb=2, ntp=1, barostat=${BAROSTAT}, pres0=1.0, taup=${TAUP},
  ntt=3, gamma_ln=${GAMMA_LN}, ig=-1, temp0=${TEMP0},
  ntr=0,
  ntpr=${NTPR_EQ}, ntwx=${NTWX}, ntwr=50000,
 /
EOF
}

write_prod_input () {   # $1 = number of steps remaining
    local NSTEPS=$1
cat > prod.in <<EOF
NPT production, ${NSTEPS} steps remaining, frames every 10 ps
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=${NSTEPS}, dt=${DT},
  ntc=2, ntf=2, cut=${CUT},
  ntb=2, ntp=1, barostat=${BAROSTAT}, pres0=1.0, taup=${TAUP},
  ntt=3, gamma_ln=${GAMMA_LN}, ig=-1, temp0=${TEMP0},
  ntr=0,
  ntpr=${NTPR_PROD}, ntwx=${NTWX}, ntwr=${NTWR},
  iwrap=1,
 /
EOF
}

# ig=-1 means pmemd draws a seed and ANNOUNCES it ("Setting random seed to N").
# The input echo still shows ig = -1, so match the announcement, and record it as
# it is produced -- a requeue overwrites .out files and the seed would be lost.
record_seed () {   # $1 = stage basename (expects $1.out)
    local s
    s=$(grep -m1 "Setting random seed to" "${1}.out" 2>/dev/null | awk '{print $NF}')
    [[ -n "$s" ]] && echo "$(date -Is) $1 ig=$s" >> seeds.txt
    return 0
}

# topology_for_run <system dir> -- the HMR prmtop when one exists, else the plain one
topology_for_run () {
    if [[ "$DT" == "0.004" ]]; then
        if [[ -s "$1/system_hmr.prmtop" ]]; then echo "$1/system_hmr.prmtop"; return 0; fi
        echo "  !! dt=4 fs requested but $1/system_hmr.prmtop is missing -- refusing" >&2
        return 1
    fi
    echo "$1/system.prmtop"
}
