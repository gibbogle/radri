! Cancer cell state development

module cellstate
use global
use cycle_mod
implicit none

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine GrowCells(dt,t_simulation,ok)
real(REAL_KIND) :: dt, t_simulation
logical :: ok
logical :: changed		! not used
type(cell_type),pointer :: cp

tnow = t_simulation		! tnow = time at the start of the timestep
ok = .true.
call grower(dt,changed,ok)
if (.not.ok) return
!call CellDeath(dt,ok)
!if (.not.ok) return
end subroutine

!-----------------------------------------------------------------------------------------
! Irradiate cells with dose
!-----------------------------------------------------------------------------------------
subroutine Irradiation(dose,ok)
real(REAL_KIND) :: dose
logical :: ok
integer :: kcell, ityp
real(REAL_KIND) :: R, total
integer :: phase_count(0:4)
real(REAL_KIND) :: ph_dist(0:4), totDSB(2)
integer :: counts(8)
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

counts = 0
NPsurvive = 0
Napop = 0
Nmitotic = 0
do kcell = 1,nlist
    kcell_now = kcell
    cp => cell_list(kcell)
    counts(cp%phase) = counts(cp%phase) + 1
	if (cp%state == DEAD .or. cp%state == DYING) cycle
    ityp = cp%celltype
	ccp => cc_parameters(ityp)
	if (use_suppression) then
		if (cp%phase >= S_phase) then
			fsup = ksup**nsup/(ksup**nsup + (dose-dose_threshold)**nsup)
		else
			fsup = 1
		endif
	else
		fsup = 1.0
	endif
    call cellIrradiation(cp,dose)
enddo   
total = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    total = total + cp%totDSB0
enddo
write(nflog,*) 'At irradiation, total DSB: ',total
end subroutine

#if 0
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine CellDeath(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, ityp
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp

ok = .true.
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
#endif

!-----------------------------------------------------------------------------------------
! Cell growth, death and division are handled here.  
!-----------------------------------------------------------------------------------------
subroutine grower(dt, changed, ok)
real(REAL_KIND) :: dt
logical :: changed, ok
integer :: k, kcell, nlist0, ityp!, kpar=0
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: mitosis_duration
integer, parameter :: MAX_DIVIDE_LIST = 100000
integer :: ndivide, divide_list(MAX_DIVIDE_LIST)
logical :: divide

ok = .true.
changed = .false.
nlist0 = nlist
ndivide = 0

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
	if (cp%phase < M_phase) then
	    call growcell(cp,dt)
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
        if (cp%mitosis >= 1) then
			cp%G2_time = tnow - cp%t_start_G2
            if (use_SF) then
			    if (is_radiation .and. cp%Psurvive < 0) then
			        call survivalProbability(cp)
			    endif
				divide = .false.	! no division when SFave is to be computed
			else
			    divide = .true.
			endif
		endif
    endif
    ! end cell simulation---------------------------------------------------------------------
    
	if (divide) then
		ndivide = ndivide + 1
		if (ndivide > MAX_DIVIDE_LIST) then
		    write(nflog,*) 'Error: growcells: MAX_DIVIDE_LIST exceeded: ',MAX_DIVIDE_LIST
		    ok = .false.
		    return
		endif
		divide_list(ndivide) = kcell
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
!-----------------------------------------------------------------------------------------
subroutine growcell(cp, dt)
type(cell_type), pointer :: cp
real(REAL_KIND) :: dt, f_CP

if (cp%phase >= M_phase) then
	return
endif
! Here compute the progress factor %fp, used in cycle.f90
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
integer :: kcell2, ityp
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
cp1%t_divide_last = tnow

! Jaiswal
cp1%CC_act = 0
if (use_cell_kcc2a_dependence) cp1%Kcc2a = get_Kcc2a(kmccp,CC_tot,CC_threshold_factor,cp1%fg(G2_phase)*ccp%T_G2/3600)
cp1%irradiated = (tnow > t_irradiation)
if (is_radiation .and. use_SF) return   ! in this case there is no need to actually have the cell divide
cp1%DSB(NHEJslow,:) = cp1%DSB(NHEJslow,:)/2
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
