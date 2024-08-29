! Cancer cell state development

module cellstate
use global
use cycle_mod
implicit none

integer :: kcell_dividing = 0
logical :: first

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine GrowCells(dt,t_simulation,ok)
real(REAL_KIND) :: dt, t_simulation
logical :: ok
logical :: changed		! not used
integer :: kcell, idrug, ichemo
type(cell_type),pointer :: cp

tnow = t_simulation		! now = time at the start of the timestep
ok = .true.
call grower(dt,changed,ok)
if (.not.ok) return
call CellDeath(dt,ok)
if (.not.ok) return
end subroutine

!-----------------------------------------------------------------------------------------
! Irradiate cells with dose (and duration tmin to be added)
!-----------------------------------------------------------------------------------------
subroutine Irradiation(dose,ok)
real(REAL_KIND) :: dose
logical :: ok
integer :: kcell, site(3), iv, ityp, n, kpar=0
real(REAL_KIND) :: C_O2, SER, p_death, p_recovery, R, kill_prob, tmin, Cdrug, total
real(REAL_KIND) :: SER_OER(2)
integer :: phase_count(0:4)
real(REAL_KIND) :: ph_dist(0:4), totDSB(2)
integer :: counts(8), jpp
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp

ok = .true.
Nirradiated = Ncells
t_irradiation = t_simulation
write(nflog,'(a,f8.3,i6)') 'Irradiation: t, Nirradiated: ',t_irradiation/3600,Nirradiated
call get_phase_distribution(phase_count)
total = sum(phase_count)
ph_dist = 100*phase_count/total
write(nflog,'(a,5f8.2)') 'phase distribution: ',ph_dist

nslow_sum = 0
pHR_sum = 0
pNHEJslow_sum = 0
fdecay_sum = 0

if (use_G1_CP_factor) then
    G1_CP_time = G1_CP_factor*dose*3600
endif
counts = 0
NPsurvive = 0
Napop = 0
Nmitotic = 0
tmin = 1.0      ! for now...
n = 0
do kcell = 1,nlist
    kcell_now = kcell
    cp => cell_list(kcell)
    counts(cp%phase) = counts(cp%phase) + 1
	if (cp%state == DEAD .or. cp%state == DYING) cycle
    ityp = cp%celltype
	ccp => cc_parameters(ityp)
	if (use_Iliakis) then
		if (cp%phase >= S_phase) then
			fIliakis = kIliakis**nIliakis/(kIliakis**nIliakis + (dose-dose_threshold)**nIliakis)
		else
			fIliakis = 1
		endif
	else
		fIliakis = 1.0
	endif
    call cellIrradiation(cp,dose)
enddo   
total = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    total = total + cp%totDSB0
enddo

write(*,*) 'Irradiation: phase counts: ',counts
write(nflog,*) 'At irradiation, total DSB: ',total
end subroutine

!-----------------------------------------------------------------------------------------
! Cells can be tagged to die, or finally die of anoxia or aglucosia, or they can be tagged 
! for death at division time if the drug is effective.
! Note: if simulating colony, no tagging, no death from anoxia, aglucosia
!-----------------------------------------------------------------------------------------
subroutine CellDeath(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, nlist0, site(3), i, ichemo, idrug, im, ityp, killmodel, kpar=0 
real(REAL_KIND) :: C_O2, C_glucose, Cdrug, n_O2, kmet, Kd, dMdt, kill_prob, dkill_prob, death_prob,survival_prob
real(REAL_KIND) :: t_dying
!logical :: anoxia_death, aglucosia_death
real(REAL_KIND) :: R, delayed_death_prob, factor
!type(drug_type), pointer :: dp
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp

ok = .true.
nlist0 = nlist
first = .true.
do kcell = 1,nlist
    cp => cell_list(kcell)
	ityp = cp%celltype
	ccp => cc_parameters(ityp)
	if (cp%state == DEAD) cycle
	if (cp%state == DYING) then
		cycle
	endif	
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Cell growth, death and division are handled here.  Division occurs when cell volume 
! exceeds the divide volume. 
! As the cell grows we need to adjust both Cin and Cex to maintain mass conservation.
! GROWTH DELAY
! When a cell has received a dose of radiation (or possibly drug - not yet considered)
! the cycle time is increased by an amount that depends on the dose.  The delay may be
! transmitted to progeny cells.
! 
! NOTE: now the medium concentrations are not affected by cell growth
!-----------------------------------------------------------------------------------------
subroutine grower(dt, changed, ok)
real(REAL_KIND) :: dt
logical :: changed, ok
integer :: k, kcell, nlist0, ityp, idrug, prev_phase, kpar=0
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: rr(3), c(3), rad, d_desired, R, rrsum, pdeath, mitosis_duration, f_CP
real(REAL_KIND) :: fslow(3), Cdrug
integer :: nslow(3)
integer, parameter :: MAX_DIVIDE_LIST = 100000
integer :: ndivide, divide_list(MAX_DIVIDE_LIST)
logical :: drugkilled, radkilled
logical :: divide, tagged

ok = .true.
changed = .false.
nlist0 = nlist
ndivide = 0

fslow = 0
nslow = 0
do kcell = 1,nlist0
	kcell_now = kcell
   	cp => cell_list(kcell)
    if (cp%state == DIVIDED) cycle
	if (cp%state == DEAD) cycle
	if (cp%state == DYING) then
		cycle
    endif
    
    ! start cell simulation----------------------------------------------------------
	ityp = cp%celltype
	ccp => cc_parameters(ityp)
	divide = .false.
    mitosis_duration = cp%mitosis_duration
    prev_phase = cp%phase
	if (cp%phase < M_phase) then
	    call growcell(cp,dt,f_CP)
		if (kcell == 1) then
	        nslow(cp%phase) = nslow(cp%phase) + 1
	        fslow(cp%phase) = fslow(cp%phase) + f_CP
	    endif
	endif
    call log_timestep(cp, ccp, dt)
    if (cp%phase == M_phase) then
		cp%mitosis = 0
		cp%t_start_mitosis = tnow
		ncells_mphase = ncells_mphase + 1
        cp%phase = dividing
    endif
	
    if (cp%phase == dividing) then
		cp%mitosis = (tnow - cp%t_start_mitosis)/mitosis_duration
        if (test_run) write(nflog,'(a,i6,3f10.1,f8.3)') 'grower: kcell, mitosis: ',kcell, tnow, cp%t_start_mitosis,mitosis_duration,cp%mitosis	
        if (cp%mitosis >= 1) then
			cp%G2_time = tnow - cp%t_start_G2
            if (use_SF) then
			    if (is_radiation .and. cp%Psurvive < 0) then
					if (test_run) write(nflog,'(a,i6,f6.3)') 'Exit M_phase, get P_survive: kcell, time: ',kcell,t_simulation/3600
			        call survivalProbability(cp)
			    !else
			    !    divide = .true.
			    endif
				divide = .false.	! no division when SFave is to be computed
			else
			    divide = .true.
			endif
		endif
    endif
    ! end cell simulation---------------------------------------------------------------------
    
	if (divide) then
!        if (cp%generation == 2) cycle		! we simulate cell division for PDJ, not for SFALL
		ndivide = ndivide + 1
		if (ndivide > MAX_DIVIDE_LIST) then
		    write(nflog,*) 'Error: growcells: MAX_DIVIDE_LIST exceeded: ',MAX_DIVIDE_LIST
		    ok = .false.
		    return
		endif
		divide_list(ndivide) = kcell
	endif
enddo
do k = 1,3
    if (nslow(k) > 0) then
        fslow(k) = fslow(k)/nslow(k)
    endif
enddo
do k = 1,ndivide
	changed = .true.
	kcell = divide_list(k)
   	cp => cell_list(kcell)
	call divider(kcell, ok)
	if (.not.ok) return
enddo
end subroutine


!-----------------------------------------------------------------------------------------
! Need to check phase-dependence of growth
! colony_simulation not fixed
!-----------------------------------------------------------------------------------------
subroutine growcell(cp, dt, f_CP)
type(cell_type), pointer :: cp
real(REAL_KIND) :: dt, f_CP
integer :: iph

if (cp%phase >= M_phase) then
	return
endif
! Here compute the progress factor %fp, used in cycle.f90
iph = cp%phase
f_CP = slowdown(cp)
if (.not.is_radiation .and. f_CP < 1.0) then
    write(*,*) 'no IR, growcell: f_CP < 1'
    stop
endif
cp%fp = f_CP/cp%fg(cp%phase)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine divider(kcell1, ok)
integer :: kcell1
logical :: ok
integer :: kcell2, ityp		!, nbrs0, im, imax, ipdd
type(cell_type), pointer :: cp1, cp2
type(cycle_parameters_type), pointer :: ccp
integer :: kpar = 0

ok = .true.
ityp = 1
ccp => cc_parameters(ityp)
ndivided = ndivided + 1
cp1 => cell_list(kcell1)

cp1%state = ALIVE
cp1%generation = cp1%generation + 1
cp1%birthtime = tnow
kcell_now = kcell1
cp1%mitosis_duration = get_mitosis_duration()
call set_phase_times(cp1)
cp1%mitosis = 0
cp1%phase = G1_phase
cp1%progress = 0
cp1%totMis = 0
cp1%t_divide_last = tnow

! Jaiswal
cp1%CC_act = 0	! CC_act0
if (use_cell_kcc2a_dependence) cp1%Kcc2a = get_Kcc2a(kmccp,CC_tot,CC_threshold_factor,cp1%fg(G2_phase)*ccp%T_G2/3600)
cp1%irradiated = (tnow > t_irradiation)
if (is_radiation .and. use_SF) return   ! in this case there is no need to actually have the cell divide
cp1%DSB(NHEJslow,:) = cp1%DSB(NHEJslow,:)/2		! ask Bill
cp1%DSB(NHEJfast,:) = 0
cp1%DSB(HR,:) = 0
cp1%DSB(TMEJ,:) = 0

! Second cell
if (ngaps > 0) then
    kcell2 = gaplist(ngaps)
    ngaps = ngaps - 1
else
	nlist = nlist + 1
	if (nlist > MAX_NLIST) then
		ok = .false.
		write(nflog,*) 'Error: Maximum number of cells MAX_NLIST has been exceeded.  Increase and rebuild.'
		return
	endif
	kcell2 = nlist
endif

ncells = ncells + 1
ityp = cp1%celltype
ccp => cc_parameters(ityp)
ncells_type(ityp) = ncells_type(ityp) + 1
ncells_mphase = ncells_mphase - 1
cp2 => cell_list(kcell2)
cp2 = cp1
cp2%mitosis_duration = get_mitosis_duration()
kcell_now = kcell2
call set_phase_times(cp2)
if (use_cell_kcc2a_dependence) cp2%Kcc2a = get_Kcc2a(kmccp,CC_tot,CC_threshold_factor,cp2%fg(G2_phase)*ccp%T_G2/3600)

end subroutine

end module
