!-----------------------------------------------------------------------------------------
! Units:
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!-----------------------------------------------------------------------------------------
module radri_mod
use global
use cellstate

IMPLICIT NONE

contains 

!-----------------------------------------------------------------------------------------
! This subroutine is called to initialize a simulation run. 
! infile = file with the input data
! outfile = file to hold the output 
!-----------------------------------------------------------------------------------------
subroutine Setup(infile,outfile,ok)
character*(*) :: infile, outfile
logical :: ok

ok = .true.
par_zig_init = .false.

inputfile = infile
outputfile = outfile
call ReadCellParams(ok)
if (.not.ok) return

Mnodes = 1

call ArrayInitialisation(ok)
if (.not.ok) return
write(nflog,*) 'did ArrayInitialisation'

call PlaceCells(ok)
if (.not.ok) return

istep = 0

ncells_mphase = 0

t_simulation = 0
SFdone = .false.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ArrayInitialisation(ok)
logical :: ok
integer :: x,y,z,k, ichemo
integer :: MAXX, z1, z2, nc0, inflow
integer :: cog_size
real(REAL_KIND) :: d, rr(3)

call RngInitialisation

! These are deallocated here instead of in subroutine wrapup so that when a simulation run ends
! it will still be possible to view the cell distributions and chemokine concentration fields.
if (allocated(cell_list)) deallocate(cell_list)
if (allocated(Psurvive)) deallocate(Psurvive)

nlist = 0
write(nflog,*) 'Initial count, max_nlist: ',initial_count, max_nlist
allocate(cell_list(max_nlist))
ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
! Set up Mnodes seed values for the random number generator par_zig 
!-----------------------------------------------------------------------------------------
subroutine RngInitialisation
integer, allocatable :: zig_seed(:)
integer :: i
integer :: npar, grainsize = 32

npar = Mnodes
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
par_zig_init = .true.

end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine ReadCellParams(ok)
logical :: ok
integer :: i, idrug, imetab, nmetab, im, itestcase, Nmm3, ichemo, itreatment, iuse_extra, iuse_relax, iuse_par_relax, iuse_FD
integer :: iuse_oxygen, iuse_glucose, iuse_lactate, iuse_glutamine, iuse_othernutrient, iuse_drug, iuse_metab, iV_depend
integer :: iV_random, iuse_gd_all, iuse_divide_dist, iuse_lognormal, ityp
integer :: ictype, idisplay, isconstant, ioxygengrowth, iglucosegrowth, ilactategrowth, ioxygendeath, iglucosedeath
integer :: iuse_drop, iconstant, isaveprofiledata, isaveslicedata, iusecellcycle, iusemetabolism, ifullymixed, isynchronise
logical :: use_metabolites
real(REAL_KIND) :: bdry_conc, percent, d_n_limit
real(REAL_KIND) :: sigma(2)
character*(12) :: drug_name
type(cycle_parameters_type),pointer :: ccp
logical :: write_hourly_results

ok = .true.

open(nfcell,file=inputfile,status='old')
write(*,*) 'Opened: ',trim(inputfile)
read(nfcell,*) initial_count				! initial number of tumour cells
read(nfcell,*) divide_time_median(1)
read(nfcell,*) divide_time_shape(1)
read(nfcell,*) ndays                        ! max number of days to simulate
read(nfcell,*) DELTA_T						! time step size (sec)
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
Ncelltypes = 1

call ReadCellCycleParameters(nfcell)

call ReadMcParameters(nfcell)

call ReadProtocol(nfcell)

is_radiation = .false.
close(nfcell)

! Try setting this for each cell unless use_cell_kcc_dependence
Kcc = get_Kcc(kmccp,CC_tot,CC_threshold_factor,cc_parameters(1)%T_G2/3600)
single_cell = (initial_count==1)
write(nflog,*) 'single_cell: ',single_cell

open(nfres,file='radri_ts.out',status='replace')
write(nflog,*) 'Opened radri_ts.out'

Nsteps = ndays*24*60*60/DELTA_T		! max # of steps (DELTA_T in seconds)

end subroutine

!-----------------------------------------------------------------------------------------
! The cell cycle parameters include the parameters for radiation damage and repair. 
! Time unit = hour
!-----------------------------------------------------------------------------------------
subroutine ReadCellCycleParameters(nf)
integer :: nf
type(cycle_parameters_type),pointer :: ccp
integer :: ityp=1
real(REAL_KIND) :: sgma, total

write(nflog,*) 'ReadCellCycleParameters:'
ccp => cc_parameters(1)

read(nf,*) ccp%f_G1
read(nf,*) ccp%f_S
read(nf,*) ccp%f_G2
read(nf,*) ccp%f_M

divide_dist(ityp)%class = LOGNORMAL_DIST
divide_time_median(ityp) = 60*60*divide_time_median(ityp)		! hours -> seconds
sgma = log(divide_time_shape(ityp))
divide_dist(ityp)%p1 = log(divide_time_median(ityp))	
divide_dist(ityp)%p2 = sgma
divide_time_mean(ityp) = exp(divide_dist(ityp)%p1 + 0.5*divide_dist(ityp)%p2**2)	! mean = median.exp(sigma^2/2)

call SteelMethod(ityp)

ccp%T_G1 = 3600*ccp%T_G1    ! hours -> seconds
ccp%T_S = 3600*ccp%T_S
ccp%T_G2 = 3600*ccp%T_G2
ccp%T_M = 3600*ccp%T_M

end subroutine

!-----------------------------------------------------------------------------------------
! Compute mean phase durations from phase fractions.  All corresponding to the average cycle time.
! Note that the value of ccp%T_M is overridden at when a cell is simulated 
! (unless it is a! single-cell simulation) by the value given by cp%mitosis_duration, 
! generated by get_mitosis_duration()
!-----------------------------------------------------------------------------------------
subroutine SteelMethod(ityp)
integer :: ityp
real(REAL_KIND) :: Tc, b
type(cycle_parameters_type),pointer :: ccp

ccp => cc_parameters(ityp)
Tc = divide_time_mean(ityp)/3600    ! hours
b = log(2.0)/Tc
ccp%T_G1 = -(log(1-ccp%f_G1/2))/b
ccp%T_S = -(log(1-(ccp%f_G1+ccp%f_S)/2))/b - ccp%T_G1
ccp%T_G2 = -(log(1-(ccp%f_G1+ccp%f_S+ccp%f_G2)/2))/b - ccp%T_G1 - ccp%T_S
ccp%T_M = Tc - ccp%T_G1 - ccp%T_S - ccp%T_G2
write(nflog,'(a,5f7.3)') 'modified SteelMethod mean T G1,S,G2,M, total: ',ccp%T_G1,ccp%T_S,ccp%T_G2,ccp%T_M,ccp%T_G1+ccp%T_S+ccp%T_G2+ccp%T_M
end subroutine

!-----------------------------------------------------------------------------------------
! It might be advisable to read both washout_time_h and CA_time_h ????
!-----------------------------------------------------------------------------------------
subroutine ReadProtocol(nf)
integer :: nf
real(REAL_KIND) :: halflife
character*(16) :: drugname
integer :: ndrug

read(nf,*) ndrug
if (ndrug > 0) then
    read(nf,'(a)') drugname
    read(nf,*) halflife
    write(nflog,*) 'halflife: ',halflife
    if (halflife == 0) write(nflog,*) 'No drug decay'
    read(nf,*) drug_conc0
    write(nflog,*) 'drug_conc0: ',drug_conc0
    read(nf,*) washout_time_h
    write(nflog,*) 'washout_time_h: ',washout_time_h
    if (halflife == 0) then
        Khalflife = 0
    else
        Khalflife = 0.693/halflife
    endif
    use_drug_halflife = (Khalflife > 0)
endif
read(nf,*) CA_time_h
write(nflog,*) 'CA_time_h: ',CA_time_h
read(nf,*) radiation_dose
write(nflog,*) 'radiation_dose: ',radiation_dose
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine PlaceCells(ok)
logical :: ok
integer :: kcell
type(cell_type), pointer :: cp
type(cycle_parameters_type),pointer :: ccp
	
ccp => cc_parameters(1)
t_irradiation = -1
do kcell = 1,initial_count
	call AddCell(kcell)
enddo
nlist = kcell-1
Ncells = nlist
Ncells0 = Ncells

ok = .true.
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine AddCell(kcell)
integer :: kcell
integer :: ityp, kpar = 0
real(REAL_KIND) :: R, kfactor
type(cell_type), pointer :: cp
type(cycle_parameters_type),pointer :: ccp
	
cp => cell_list(kcell)
cp%ID = kcell
cp%state = ALIVE
cp%generation = 1
cp%birthtime = 0
cp%celltype = 1
ityp = cp%celltype
ccp => cc_parameters(ityp)
cp%mitosis_duration = get_mitosis_duration()
kcell_now = kcell

! Jaiswal
R = par_uni(kpar)
kfactor = 1 + (R - 0.5)*jaiswal_std
if (single_cell) kfactor = 1
cp%kccmd = kccmd*kfactor
R = par_uni(kpar)
kfactor = 1 + (R - 0.5)*jaiswal_std
if (single_cell) kfactor = 1
cp%kccrd = kccrd*kfactor

cp%CC_act = 0
cp%ATR_act = 0
cp%ATM_act = 0
cp%G2_time = 0
CP%phase = G1_phase
cp%progress = 0
cp%Psurvive = -1    ! flags Psurvive not yet computed

end subroutine

!--------------------------------------------------------------------------------------
! Steel: 
! Probability density function of progress through cell cycle: f(t) = 2b exp(-bt) 
! => cumulative distribution function F(t) = 2(1 - exp(-bt)) = fraction of cells less than t since creation
! To generate a variate from this CDF, first generate R from U(0,1)
! R = 2(1 - exp(-bt)), exp(-bt) = 1 - R/2, -bt = ln(1 - R/2)
! t = -(1/b)ln(1 - R/2)
! If synchronisation of cell initialisation is specified in radri_main, then
! all cells start the simulation at the same point in the cell cycle,
! i.e. same phase and progress.
!--------------------------------------------------------------------------------------
subroutine SetInitialCellCycleStatus(kcell,cp)
integer :: kcell
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
integer :: ityp, kpar = 0
real(REAL_KIND) :: Tc, Tmean, scale, b, t, R, tswitch(3), fg(4), metab, f_CP, fp(4)
real(REAL_KIND) :: T_G1, T_S, T_G2, T_M, tleft, Vleft, dth

ityp = cp%celltype
ccp => cc_parameters(ityp)
Tmean = divide_time_mean(ityp)
if (single_cell .or. test_run) then
    Tc = divide_time_mean(1)
else
    Tc = cp%divide_time         ! log-normal, implies %fg
endif
scale = Tc/Tmean
fg = cp%fg
f_CP = 1.0
fp(:) = f_CP/fg(:)
T_G1 = ccp%T_G1/fp(1)
T_S = ccp%T_S/fp(2)
T_G2 = ccp%T_G2/fp(3)
if (test_run) then
    T_M = ccp%T_M/fp(4)
else
    T_M = cp%mitosis_duration
endif
if (use_cell_kcc_dependence) then
    cp%Kcc = get_Kcc(kmccp,CC_tot,CC_threshold_factor,T_G2/3600)
    cp%Kcc = min(cp%kcc, 0.9*CC_threshold)
endif

if (use_synchronise) then
    if (synch_phase == G1_phase) then
        t = synch_fraction*T_G1
    elseif (synch_phase == S_phase) then
        t = T_G1 + synch_fraction*T_S
    elseif (synch_phase == G2_phase) then
        t = T_G1 + T_S + synch_fraction*T_G2
    else
        write(*,*) 'Error: SetInitialCellCycleStatus: bad synch_phase: ',synch_phase
        write(nflog,*) 'Error: SetInitialCellCycleStatus: bad synch_phase: ',synch_phase
        stop
    endif
    cp%progress = synch_fraction
else
    b = log(2.0)/Tc
    R = par_uni(kpar)
    t = -(1/b)*log(1 - R/2)     ! cycle progression, log-normal r.v. (t/Tc = fractional progression)
endif
tswitch(1) = T_G1 
tswitch(2) = tswitch(1) + T_S
tswitch(3) = tswitch(2) + T_G2

if (t < tswitch(1)) then
    cp%phase = G1_phase
    cp%fp = fp(1)
    cp%progress = t/T_G1
elseif (t <= tswitch(2)) then
    cp%phase = S_phase
    cp%fp = fp(2)
    cp%progress = (t - tswitch(1))/T_S
elseif (t <= tswitch(3)) then
    cp%phase = G2_phase
    cp%fp = fp(3)
    cp%progress = (t - tswitch(2))/T_G2
    tleft = tswitch(3) - t
    if (use_Jaiswal) then
        cp%DSB = 0
        dth = (t - tswitch(2))/3600
        call Jaiswal_update(cp,dth)
        cp%CC_act = min(cp%CC_act,0.95*CC_threshold)     ! to prevent premature mitosis
    endif
else    ! cell in mitosis
    R = par_uni(kpar)
    cp%t_start_mitosis = -(t - tswitch(3))
    cp%progress = (t - tswitch(3))/T_M
	ncells_mphase = ncells_mphase + 1
    cp%phase = dividing
endif
if (single_cell) then
    write(*,*)
    write(*,*) 'Initial phase, progress: ',cp%phase,cp%progress
    write(nflog,*) 'Initial phase, progress: ',cp%phase,cp%progress
    write(*,*)
endif
cp%t_divide_last = -t
end subroutine

!-----------------------------------------------------------------------------------------
! Advance simulation through one time step (DELTA_T)
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step  
use, intrinsic :: iso_c_binding
integer(c_int) :: res
integer :: kcell, hour, nthour, kpar=0
integer :: nphaseh(8), iph
type(cell_type), pointer :: cp
integer :: phase_count(0:4)
real(REAL_KIND) :: total, tIR
real(REAL_KIND) :: SFtot, Pp, Pd
real(REAL_KIND) :: Cdrug
integer :: Ntot, Ndying, Ncont(5)
logical :: ok = .true. 

t_simulation = istep*DELTA_T	! seconds
nthour = 3600/DELTA_T
Cdrug = drug_conc0
if (Cdrug == 0) then
    fDNAPK = 1
else
    if (use_drug_halflife) then
        Cdrug = Cdrug*exp(-Khalflife*(t_simulation - drug_time)/3600)
    endif
    fDNAPK = logistic(Cdrug)    ! this is DNAPKact
endif

cp => cell_list(1)
tIR = istep*DELTA_T/3600.0
if (single_cell .and. (cp%phase < G2_phase)) write(nflog,'(a,f6.2,i4,5f8.3)') 'tIR,phase,progress,CC, ATR, ATM_act,fp: ',tIR,cp%phase,cp%progress,cp%CC_act,cp%ATR_act,cp%ATM_act,cp%fp

if (.not.is_radiation) then
	write(nflog,'(a,f6.1)') 'Radiation dose: ',radiation_dose
    do kcell = 1,Ncells
        cp => cell_list(kcell)
        call set_phase_times(cp)
        call SetInitialCellCycleStatus(kcell,cp)
    enddo
	call Irradiation(radiation_dose, ok)
	if (.not.ok) then
		res = 3
		return
    endif
    is_radiation = .true.
    IR_time_h = 0
endif

if (washout_time_h > 0 .and. drug_conc0 > 0) then     ! check for washout time
    if (t_simulation >= washout_time_h*3600) then
        write(nflog,'(a,i6,f8.1)') 'Drug washout: istep,time: ',istep,t_simulation/3600
        write(*,'(a,f8.1)') 'Drug washout: time: ',t_simulation/3600
        write(nflog,'(a,f8.3)') 'drug exposure time: ',(t_simulation - t_irradiation)/3600
        drug_conc0 = 0
    endif
endif
res = 0

if (t_irradiation >= 0) call GrowCells(DELTA_T,t_simulation,ok)

kcell = 1
cp => cell_list(kcell)

if (compute_cycle) then
    call get_phase_distribution(phase_count)
    total = sum(phase_count)
    phase_dist = 100*phase_count/total
endif
if (compute_cycle) then
    if (next_phase_hour > 0) then  ! check if this is a phase_hour
        if (real(istep)/nthour >= phase_hour(next_phase_hour)) then   ! record phase_dist
            write(*,*) 'Reached phase hour: ',next_phase_hour,phase_hour(next_phase_hour)
            if (next_phase_hour <= 9) then
                write(*,'(a,4i8)') 'count: ',phase_count(1:4)
                write(nflog,*) 'Reached phase hour: ',next_phase_hour,phase_hour(next_phase_hour)
                write(nflog,'(a,4i8)') 'count: ',phase_count(1:4)
                write(nflog,'(a,4f8.3)') 'dist: ',phase_dist(1:4)
            endif
	        if (compute_cycle) then
                recorded_phase_dist(next_phase_hour,1:4) = 100*phase_dist(1:4)/sum(phase_dist(1:4))
	        endif
            next_phase_hour = next_phase_hour + 1
            if (next_phase_hour > nphase_hours) next_phase_hour = 0
        endif
    endif
endif

if (mod(istep,nthour) == 0) then
    hour = istep/nthour
    nphaseh = 0
    do kcell = 1,nlist
        cp => cell_list(kcell)
        iph = cp%phase
        nphaseh(iph) = nphaseh(iph) + 1
    enddo
	if (.not. single_cell) write(*,'(a,i6,i4,4(a,i8))') 'istep, hour: ',istep,hour,' Ncells: ',Ncells   
    call get_phase_distribution(phase_count)
    total = sum(phase_count(1:4))
    phase_dist = 100*phase_count/total
    
    total = 0
    do kcell = 1,nlist
        cp => cell_list(kcell)
        total = total + cp%totDSB0
    enddo    
endif

istep = istep + 1
overstepped = (istep == maxhours*nthour)
if (overstepped) then
    write(*,*) 'overstepped the mark'
    call nondivided()
    call completed
    res = 1
    return
endif

! Check for completion of the run
if (compute_cycle .or. output_DNA_rate) then
    if (next_phase_hour == 0) then
        call completed
        res = 1
    endif
    return
endif

if (SFdone) then
    Ntot = 0
    SFtot = 0
    Ndying = 0
    Ncont = 0
    do kcell = 1,nlist
        cp => cell_list(kcell)
        Pp = cp%Psurvive
        if (Pp > 0) then
            if (cp%mitosis_time < CA_time_h*3600) then  ! adjust for 2 daughters
                Pd = 1 - sqrt(1.0 - Pp)
                Ntot = Ntot + 2
                SFtot = SFtot + 2*Pd
                if (cp%phase0 == M_phase) then
                    Ncont(3) = Ncont(3) + 2
                else
                    Ncont(1) = Ncont(1) + 2
                endif
            else                                        ! no daughter adjustment
                Ntot = Ntot + 1
                SFtot = SFtot + Pp
                if (cp%phase0 == M_phase) then
                    Ncont(4) = Ncont(4) + 1
                else
                    Ncont(2) = Ncont(2) + 1
                endif
            endif
        else                                            ! state = DYING, Psurvive = 0
            Ntot = Ntot + 1
            Ndying = Ndying + 1
            Ncont(5) = Ncont(5) + 1
        endif
    enddo
! Note that: 2 cells are contributed by 1 pre-CA interphase cell (Ncont(1)), 1 by 1 post-CA interphase cell (Ncont(2))
! 4 cells are contributed by 1 pre-CA surviving mitotic cell (Ncont(3)), 2 by a post-CA surving mitotic cells (Ncont(4))
! 2 cells are contributed by one dying mitotic cell (Ncont(5))
! Ncells0 = Ncont(1)/2 + Ncont(2) + Ncont(3)/4 + Ncont(4)/2 + Ncont(5)/2

    SFave = SFtot/Ntot
    write(*,*)
!    write(nflog,'(a,7i6)') 'Ncont, Ndying, Ncells0: ',Ncont,Ndying,Ncont(1)/2 + Ncont(2) + (Ncont(3)/2 + Ncont(4))/2 + Ncont(5)/2
!    write(nflog,'(a,i6,2x,f8.3)') 'Ntot, SFtot: ',Ntot,SFtot
    write(nflog,'(a,e12.4,f8.3)') 'SFave,log10(SFave): ',SFave,log10(SFave)
    write(*,'(a,e12.4,f8.3)') 'SFave,log10(SFave): ',SFave,log10(SFave)
    call completed
    res = 1
endif

end subroutine

!-----------------------------------------------------------------------------------------
! When the mark is overstepped, locate cells that have not reached mitosis
!-----------------------------------------------------------------------------------------
subroutine nondivided
integer :: kcell, n
type(cell_type), pointer :: cp

n = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%psurvive < 0) then
        n = n+1
        write(*,'(a,i6,i3,3f8.4)') 'nondivided: kcell, phase: ',kcell, cp%phase
        write(*,*) 'fg: ',cp%fg
        write(nflog,'(a,i6,i3,3f8.4)') 'nondivided: kcell, phase: ',kcell, cp%phase
        write(nflog,*) 'fg: ',cp%fg
    endif
enddo
write(*,*) 'Total nondivided: ',n
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine completed
integer :: kcell, ph, nir(4), nmitosis,nsum, kcellmax, i, j, k, ityp
real(REAL_KIND) :: sftot_phase(4), sfsum, sfmax
type(cycle_parameters_type), pointer :: ccp
type(cell_type), pointer :: cp
logical :: only_M_phase = .false.
logical :: PDS4 = .false.
real(REAL_KIND) :: dt, phi, PDS4_M(3) = [0.191, 0.414286, 0.732812]
real(REAL_KIND) :: normalised_phase_dist(60,0:4)   
REAL(REAL_KIND) :: ave(15), SFMave

if (overstepped) then
    SFave = 1.0E-6
    goto 99
endif
if (compute_cycle) then
    write(nflog,*) 'Completed compute cycle'
    write(*,*) 'Completed compute cycle'
    do i = 1,nphase_hours
        write(nflog,'(f6.1,4x,4f8.3)') phase_hour(i),recorded_phase_dist(i,1:4)
        write(*,'(f6.1,4x,4f8.3)') phase_hour(i),recorded_phase_dist(i,1:4)
        if (PDS4) then
            phi = phi + (recorded_phase_dist(i,4) - PDS4_M(i))**2
        endif
    enddo
    write(*,*) 'wrote recorded_phase_dist'
    write(nflog,*) 'wrote recorded_phase_dist'
    if (PDS4) then 
        write(nflog,'(a)') '    hour    expt   model   error'
        write(*,'(a)') '    hour    expt   model   error'
        do i = 1,nphase_hours
            write(nflog,'(f8.1,3f8.4)') phase_hour(i),PDS4_M(i),recorded_phase_dist(i,4),recorded_phase_dist(i,4) - PDS4_M(i)
            write(*,'(f8.1,3f8.4)') phase_hour(i),PDS4_M(i),recorded_phase_dist(i,4),recorded_phase_dist(i,4) - PDS4_M(i)
        enddo
        write(nflog,'(a,f6.3)') '    phi: ',phi
        write(*,'(a,f6.3)') '    phi: ',phi
    endif
    
    if (nphase_hours > 0) then
        if (only_M_phase) then
            write(nfres,'(20e15.6)') (recorded_phase_dist(i,4),i=1,nphase_hours)
        else
            if (normalise) then
                ityp = 1
                ccp => cc_parameters(ityp)
                control_ave(1) = 100*ccp%f_G1
                control_ave(2) = 100*ccp%f_S
                control_ave(3) = 100*ccp%f_G2
                control_ave(4) = 100*ccp%f_M

                write(*,*) 'Normalising PDs'
                write(nflog,*) 'Normalising PDs'
                write(nflog,'(a,4f8.3)') 'control: ',control_ave(1:4)
                dt = 0.5
                do i = 1,nphase_hours
                    do j = 1,4
                        normalised_phase_dist(i,j) = recorded_phase_dist(i,j)/control_ave(j)
                    enddo
                    write(nflog,'(f6.1,4f8.3)') phase_hour(i),normalised_phase_dist(i,1:4)
                enddo
                write(*,*) 'write PEST output'
                write(*,'(a,a,i6)') 'expt_tag,nphase_hours: ',expt_tag,nphase_hours
                write(nfres,'(20f8.5)') (normalised_phase_dist(i,1:4),i=1,nphase_hours)
                write(nflog,'(20f8.5)') (normalised_phase_dist(i,1:4),i=1,nphase_hours)
            else
                write(*,*) 'Not normalising PDs'
                write(nfres,'(20e15.6)') (recorded_phase_dist(i,1:4),i=1,nphase_hours)
            endif
        endif
    endif
endif
if (output_DNA_rate) then
    write(nflog,*) 'Completed'
    write(*,*) 'Completed'
    if (nphase_hours > 0) then
        write(*,*) 'write DNA_rate'
        write(nfres,'(20e15.6)') (recorded_DNA_rate(i),i=1,nphase_hours)
        do i = 1,nphase_hours
            write(nflog,'(f6.2,4x,4f6.3)') phase_hour(i),recorded_DNA_rate(i)
            write(*,'(f6.2,4x,4f6.3)') phase_hour(i),recorded_DNA_rate(i)
        enddo
    endif
    return
endif

! Look at average survival by IR phase
nir = 0
sftot_phase = 0
sfmax = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) cycle
    if (cp%totDSB0 <= 0) cycle
    if (cp%phase0 == G1_phase) then
        ph = 1
        !if (cp%Psurvive > sfmax) then
        !    sfmax = cp%Psurvive
        !    kcellmax = kcell
        !endif
    elseif (cp%phase0 == S_phase) then
        ph = 2
    elseif (cp%phase0 == G2_phase) then
        ph = 3
    else
        ph = 4
    endif
    nir(ph) = nir(ph) + 1
    sftot_phase(ph) = sftot_phase(ph) + cp%Psurvive
enddo
nmitosis = sum(nir)
!write(*,'(a,4i6)') 'nir: ',nir
!write(nflog,'(a,6f12.3)') 'totPmit, totPaber, tottotDSB: ',totPmit, totPaber, tottotDSB
!write(*,'(a,i6,5f11.1)') 'Nmitosis, totPmit, totPaber, tottotDSB: ',int(Nmitosis),totPmit, totPaber, tottotDSB
!write(*,'(a,6e12.3)') 'totPaber: ',totPaber
!write(nflog,'(a,4e12.3)') 'SFtot_phase: ',SFtot_phase
!write(nflog,'(a,i6,f10.5)') 'Nmitosis,SFtot: ',Nmitosis,sum(SFtot_phase)
! adjust for pre-rep doubling of misjoins
totNmisjoins(1) = 2*totNmisjoins(1)
write(nflog,'(a,7f9.3)') 'Ave (pre, post) NDSB, Nmisjoins: ', &
    totNDSB/nmitosis,totNmisjoins/nmitosis,sum(totNmisjoins)/nmitosis
!write(*,'(a,7f9.3)') 'Ave (pre, post) NDSB, Nmisjoins: ', &
!    totNDSB/nmitosis,totNmisjoins/nmitosis,sum(totNmisjoins)/nmitosis

99 continue
if (use_synchronise) call G2_time_distribution()

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine G2_time_distribution
integer :: kcell, n, k
integer :: cnt(40)
real(REAL_KIND) :: t
type(cell_type), pointer :: cp

n = 0
cnt = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) cycle
    t = cp%G2_time/3600
    if (t == 0) cycle
    n = n+1
    k = t + 1
    cnt(k) = cnt(k) + 1
enddo
write(*,*) 'G2_time distribution (h): n: ',n
write(nflog,*) 'G2_time distribution (h): n: ',n
do k = 1,40
    if (cnt(k) > 0) then
        write(*,'(i2,a,i2,i6,f8.3)') k-1, '-', k, cnt(k), cnt(k)/real(n)
        write(nflog,'(i2,a,i2,i6,f8.3)') k-1, '-', k, cnt(k), cnt(k)/real(n)
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Execute(infile_array,inbuflen,outfile_array,outbuflen,res) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char), intent(in) :: infile_array(*), outfile_array(*)
integer(c_int) :: inbuflen, outbuflen, res
character*(2048) :: infile, outfile, logfile
character*(13) :: fname
logical :: ok, isopen
integer :: i

res = 0
infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

inquire(unit=nflog,OPENED=isopen)
if (isopen) then
	close(nflog)
endif
i = index(infile,'.')
logfile = infile(1:i)//'log'
write(*,*) 'infile: ',trim(infile)
write(*,*) 'logfile: ',trim(logfile)
open(nflog,file=logfile,status='replace')

write(nflog,*) 'inputfile:  ', trim(infile)
write(nflog,*) 'outputfile: ', trim(outfile)
call Setup(infile,outfile,ok)
if (ok) then
	res = 0
else
    res = 1
	write(nflog,*) '=== Setup failed ==='
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run 
use, intrinsic :: iso_c_binding
integer(c_int) :: res

call Wrapup

if (res == 0) then
	write(nflog,*) ' Execution successful!'
elseif (res == -1) then
	write(nflog,*) ' Execution stopped'
elseif (res == 2) then
	write(nflog,*) ' No more live cells'
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Wrapup
integer :: ierr, ichemo, idrug
logical :: isopen

ierr = 0
!if (allocated(gaplist)) deallocate(gaplist,stat=ierr)

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)

if (par_zig_init) then
	call par_zigfree
endif
end subroutine

end module
