module cycle_mod
use real_kind_mod
use global
use mc

implicit none

contains

!--------------------------------------------------------------------------
! Average phase durations have been computed to sum to the average cycle time.
! The average phase times are stored in ccp:
! ccp%T_G1, ccp%T_S, ccp%T_G2, ccp%T_M
! Each individual cell is randomly assigned a cycle time, initially and after
! cell division.  The cycle time is a random variate drawn from a lognormal
! distribution. Cell cycle time and phase durations are calculated in subroutine set_phase_times().
! At the same time, multiplying factors fg(:) are computed.  These are the basic
! factors multiplying the average phase times for G1, S, G2, and M. 
! assuming unconstrained growth (no checkpoint delays).
! Tinter_ave = average interphase time = ccp%T_G1 + ccp%T_S + ccp%T_G2
! Tdiv = generated cycle time for the cell (lognormal RV)
! then cp%fg = (Tdiv-cp%T_M)/Tinter_ave
!--------------------------------------------------------------------------
subroutine log_timestep(cp, ccp, dt)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
real(REAL_KIND) :: tIR, Cdrug, totDSB(2)
integer :: Nwrite, jpp
logical :: switch

tIR = (tnow - t_irradiation)/3600
10 continue
if (cp%phase == G1_phase) then
    cp%progress = cp%progress + cp%fp*dt/ccp%T_G1
    if (cp%progress >= 1) then
        cp%phase = S_phase
        cp%progress = 0
        cp%t_S_phase = tnow
    endif
elseif (cp%phase == S_phase) then
    cp%progress = cp%progress + cp%fp*dt/ccp%T_S
    if (cp%progress >= 1) then
        cp%phase = G2_phase
        cp%progress = 0
        cp%t_start_G2 = istep*DELTA_T
        if (single_cell) then
            do jpp = 1,2
                totDSB(jpp) = sum(cp%DSB(:,jpp))
            enddo
        endif
    endif
elseif (cp%phase == G2_phase) then
    if (use_Jaiswal) then
        if (is_radiation) then      ! post-IR
            tIR = (t_simulation - t_irradiation)/3600   ! time since IR, in hours
            switch = (cp%CC_act >= CC_threshold)
        else
            switch = (cp%CC_act >= CC_threshold)
        endif

        if (switch) then
            cp%phase = M_phase
            cp%progress = 0
            cp%t_mitosis = tIR
        endif
    else
        cp%progress = cp%progress + cp%fp*dt/ccp%T_G2
        if (cp%progress >= 1) then
            cp%phase = M_phase
            cp%progress = 0
        endif
    endif
elseif (cp%phase == M_phase) then
    ! We never get here - in grower() %phase is changed to dividing
    ! do nothing - in growcells the phase is immediately changed to cp%dividing, and the mitosis timer starts
endif
20  continue
! Note: if cp%phase = dividing, no repair
if (cp%phase < M_phase) then
    call updateRepair(cp, dt)
endif
end subroutine

end module

