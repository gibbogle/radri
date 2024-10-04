!-----------------------------------------------------------------------------------------
! Units:
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM
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
character*(64) :: msg
integer :: ichemo, error, kcell, idrug, ityp
type(cycle_parameters_type),pointer :: ccp
logical :: isopen

ok = .true.
par_zig_init = .false.

inputfile = infile
outputfile = outfile
call ReadCellParams(ok)
if (.not.ok) return

start_wtime = wtime()

Mnodes = 1

call ArrayInitialisation(ok)
if (.not.ok) return
write(nflog,*) 'did ArrayInitialisation'

call PlaceCells(ok)
if (.not.ok) return

istep = 0
ndivided = 0
Ndying = 0
Ndead = 0

ncells_mphase = 0

t_simulation = 0
allocate(nphase(0:ndays*24,8))
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ArrayInitialisation(ok)
logical :: ok
integer :: x,y,z,k, ichemo
integer :: MAXX, z1, z2, nc0, inflow
integer :: cog_size
real(REAL_KIND) :: d, rr(3)

ok = .false.
call RngInitialisation

! These are deallocated here instead of in subroutine wrapup so that when a simulation run ends
! it will still be possible to view the cell distributions and chemokine concentration fields.
if (allocated(cell_list)) deallocate(cell_list)
if (allocated(gaplist)) deallocate(gaplist)
if (allocated(nphase)) deallocate(nphase)
if (allocated(Psurvive)) deallocate(Psurvive)

ngaps = 0
nlist = 0

write(nflog,*) 'Initial count, max_nlist: ',initial_count, max_nlist

allocate(cell_list(max_nlist))
allocate(gaplist(max_ngaps))

ok = .true.

end subroutine

!----------------------------------------------------------------------------------------- 
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
character*(1) :: numstr
type(cycle_parameters_type),pointer :: ccp
logical :: write_hourly_results

ok = .true.

open(nfcell,file=inputfile,status='old')
read(nfcell,*) initial_count				! initial number of tumour cells
read(nfcell,*) divide_time_median(1)
read(nfcell,*) divide_time_shape(1)
read(nfcell,*) ndays							! number of days to simulate
read(nfcell,*) DELTA_T						! time step size (sec)
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
Ncelltypes = 1

call ReadCellCycleParameters(nfcell)

call ReadMcParameters(nfcell)
call ReadProtocol(nfcell)
is_radiation = .false.
close(nfcell)

! Try setting this for each cell unless use_cell_kcc2a_dependence
Kcc2a = get_Kcc2a(kmccp,CC_tot,CC_threshold_factor,cc_parameters(1)%T_G2/3600)
single_cell = (initial_count==1)
write(nflog,*) 'single_cell: ',single_cell

if (use_PEST) then
else
    open(nfout,file=outputfile,status='replace')
endif
if (.not.(use_PEST .and. .false.)) then
    if (.not.use_PEST) then
	    open(nfres,file='radri_ts.out',status='replace')
	    write(nflog,*) 'Opened radri_ts.out'
    else
	    open(nfres,file=outputfile,status='replace')
	    write(nflog,*) 'Opened ',trim(outputfile)
    endif
endif

Nsteps = ndays*24*60*60/DELTA_T		! DELTA_T in seconds

end subroutine

!-----------------------------------------------------------------------------------------
! The cell cycle parameters include the parameters for radiation damage and repair. 
! Time unit = hour
!-----------------------------------------------------------------------------------------
subroutine ReadCellCycleParameters(nf)
integer :: nf
type(cycle_parameters_type),pointer :: ccp
integer :: ityp, i
real(REAL_KIND) :: sigma, total

write(nflog,*) 'ReadCellCycleParameters:'
do ityp = 1,1
    ccp => cc_parameters(ityp)

    read(nf,*) ccp%f_G1
    read(nf,*) ccp%f_S
    read(nf,*) ccp%f_G2
    read(nf,*) ccp%f_M

    divide_dist(ityp)%class = LOGNORMAL_DIST
    divide_time_median(ityp) = 60*60*divide_time_median(ityp)		! hours -> seconds
    sigma = log(divide_time_shape(ityp))
    divide_dist(ityp)%p1 = log(divide_time_median(ityp))	
    divide_dist(ityp)%p2 = sigma
    divide_time_mean(ityp) = exp(divide_dist(ityp)%p1 + 0.5*divide_dist(ityp)%p2**2)	! mean = median.exp(sigma^2/2)

    call SteelMethod(ityp)

    total = ccp%T_G1 + ccp%T_S + ccp%T_G2 + ccp%T_M
    write(nflog,'(a,8f8.3)') 'T_G1,T_S,T_G2,T_M, total: ',ccp%T_G1,ccp%T_S,ccp%T_G2,ccp%T_M, total
    ccp%T_G1 = 3600*ccp%T_G1    ! hours -> seconds
    ccp%T_S = 3600*ccp%T_S
    ccp%T_G2 = 3600*ccp%T_G2
    ccp%T_M = 3600*ccp%T_M
enddo

end subroutine

!-----------------------------------------------------------------------------------------
! Compute phase durations from phase fractions.  All corresponding to the average cycle time.
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
ccp%T_M = log(1 + ccp%f_M)/b
ccp%T_G2 = Tc - ccp%T_G1 - ccp%T_S - ccp%T_M
write(nflog,'(a,2f8.3)') 'SteelMethod: alt T_G2,T_M: ', -(log(1-(ccp%f_G1+ccp%f_S+ccp%f_G2)/2))/b - ccp%T_G1 - ccp%T_S,Tc - ccp%T_G1 - ccp%T_S - ccp%T_G2
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadProtocol(nf)
integer :: nf
real(REAL_KIND) :: halflife
character*(16) :: drugname

read(nf,'(a)') drugname
read(nf,*) halflife
write(*,*) 'halflife: ',halflife
read(nf,*) drug_conc
write(*,*) 'drug_conc: ',drug_conc
read(nf,*) washout_time_h
write(*,*) 'washout_time_h: ',washout_time_h
if (halflife == 0) then  ! 0 flags no decay of the drug
    Khalflife = 0
else
    Khalflife = 0.693/halflife
endif
use_drug_halflife = (Khalflife > 0)
if (washout_time_h < 0) then     ! this signals a CDTD expt for which CA_time = washout time, i.e. Cho1 only.  Otherwise CA_time takes the input parameter value.
    washout_time_h = -washout_time_h
    CA_time_h = washout_time_h
elseif (drug_conc == 0) then
    washout_time_h = 0           ! this signals that there is no washout
endif
read(nf,*) radiation_dose
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
Nviable = Ncells_type

ok = .true.
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine AddCell(kcell)
integer :: kcell
integer :: ityp, kpar = 0
real(REAL_KIND) :: R, kfactor
real(8) :: kcc2a_std = 0.7
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
Ncells_type(ityp) = Ncells_type(ityp) + 1
cp%mitosis_duration = get_mitosis_duration()
kcell_now = kcell

! Jaiswal
R = par_uni(kpar)
kfactor = 1 + (R - 0.5)*jaiswal_std
cp%kt2cc = kt2cc*kfactor
R = par_uni(kpar)
kfactor = 1 + (R - 0.5)*jaiswal_std
cp%ke2cc = ke2cc*kfactor

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
if (use_cell_kcc2a_dependence) then
    cp%Kcc2a = get_Kcc2a(kmccp,CC_tot,CC_threshold_factor,T_G2/3600)
    cp%Kcc2a = min(cp%kcc2a, 0.9*CC_threshold)
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
real(REAL_KIND) :: r(3), rmax, tstart, dt, dts, diam_um, framp, area, diam
integer :: i, ic, ichemo, ndt, iz, idrug, ityp, idiv, ndiv, NpreCA, Nd, Nnew
integer :: nvars, ns, nphaseh(8), ph
real(REAL_KIND) :: SFlive
type(cell_type), pointer :: cp
integer :: phase_count(0:4), nG2
real(REAL_KIND) :: total
real(REAL_KIND) :: fATM, fATR, fCP, ATM_DSB, DNA_rate
real(REAL_KIND) :: pATM_sum, pATR_sum, DSB_sum
real(REAL_KIND) :: SFtot, Pp, Pd, newSFtot, total_mitosis_time
logical :: PEST_OK
logical :: ok = .true.
logical :: dbug

integer :: kcell1, kcell2, iph1, iph2
real(8) :: prog1, prog2

t_simulation = istep*DELTA_T	! seconds
dbug = (istep < 0)
nthour = 3600/DELTA_T

if (ngaps > 200) then
	call squeezer
endif

cp => cell_list(39)
!write(nfres,'(a,2i6,8f8.3)') 'cell 39: DSB: ',istep,cp%phase,cp%progress,t_simulation/3600,cp%DSB(1:2,:),cp%Nmis

if (.not.is_radiation) then
	write(nflog,'(a,f6.1)') 'Radiation dose: ',radiation_dose
    write(nflog,*) 'before npar_uni, npar_rnor = ',npar_uni,npar_rnor
    do kcell = 1,Ncells
        cp => cell_list(kcell)
        call set_phase_times(cp)
        call SetInitialCellCycleStatus(kcell,cp)
    enddo
    write(nflog,*) 'after npar_uni, npar_rnor = ',npar_uni,npar_rnor
	call Irradiation(radiation_dose, ok)
	if (.not.ok) then
		res = 3
		return
    endif
    is_radiation = .true.
    IR_time_h = 0
endif

if (washout_time_h > 0 .and. drug_conc > 0) then     ! check for washout time
    if (t_simulation >= washout_time_h*3600) then
        write(nflog,'(a,i6,f8.1)') 'Drug washout: istep,time: ',istep,t_simulation/3600
        write(*,'(a,f8.1)') 'Drug washout: time: ',t_simulation/3600
        write(nflog,'(a,f8.3)') 'drug exposure time: ',(t_simulation - t_irradiation)/3600
        write(nflog,*) 'npar_uni, npar_rnor = ',npar_uni,npar_rnor
        drug_conc = 0
    endif
endif
res = 0

call GrowCells(DELTA_T,t_simulation,ok)

call getNviable

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

if (dbug .or. mod(istep,nthour) == 0) then
    hour = istep/nthour
    nphaseh = 0
    do kcell = 1,nlist
        cp => cell_list(kcell)
        i = cp%phase
        nphaseh(i) = nphaseh(i) + 1
    enddo
    nphase(hour,:) = nphaseh
	write(*,'(a,i6,i4,4(a,i8))') 'istep, hour: ',istep,hour,' Nlive: ',Ncells
	write(nflog,'(a,i6,i4,4(a,i8))') 'istep, hour: ',istep,hour,' Nlive: ',Ncells
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
PEST_OK = .true.
if (use_PEST) then  
    PEST_OK = (next_phase_hour == 0)
endif
    
if (is_radiation .and. (NPsurvive >= (Nirradiated - Napop - Nmitotic)) .and. PEST_OK) then  !!! needs to change
    ! getSFlive computes the average psurvive for all cells that reach mitosis,
    ! which number NPsurvive = Nirradiated - Napop.
    call getSFlive(SFlive)
    SFtot = SFlive*(Nirradiated - Napop)
    write(*,*)
    write(nflog,'(a,4i6,e12.3)') 'NPsurvive,Nirradiated,Napop,Nmitotic,SFtot: ',NPsurvive,Nirradiated,Napop,Nmitotic,SFtot
    write(nflog,'(a,e12.3)') 'Unadjusted SFave = SFtot/NPsurvive: ',SFtot/NPsurvive
    write(nflog,'(a,e12.3)') 'SFave including apoptosis killing: ',SFtot/Nirradiated
    if (include_daughters) then
    ! To adjust SFlive to replace a cell that reached mitosis before CA with its daughter cells.
    ! In order to compare simulated SFave with SF determined by experimental CA (clonogenic analysis),
    ! SFave needs to be calculated as the average of cells that make it to CA.
        write(nflog,*) 'Accounting for daughters: from NPsurvive: ',NPsurvive
        newSFtot = 0
        Nnew = 0
        NpreCA = 0
        Nd = 0
        total_mitosis_time = 0
        do kcell = 1,nlist
            cp => cell_list(kcell)
            if (cp%state == DEAD) cycle
            total_mitosis_time = total_mitosis_time + cp%mitosis_time
            if (cp%mitosis_time < CA_time_h*3600) then
                NpreCA = NpreCA + 1
                Pp = cp%psurvive
                Pd = 1 - sqrt(1.0 - Pp)   ! Pd = psurvive for the 2 daughters: Pp = 2Pd - Pd^2
                NPsurvive = NPsurvive + 1
                SFtot = SFtot - Pp + 2*Pd
                newSFtot = newSFtot + 2*Pd
                Nnew = Nnew + 2
            endif
        enddo
        write(nflog,*) 'to NPsurvive: ',NPsurvive
        write(*,'(a,3i6,2f10.3)') 'NpreCA, Nd, Nnew, newSFtot, newSFave: ',NpreCA,Nd,Nnew,newSFtot,newSFtot/Nnew
        write(*,'(a,3i6,f8.4)') 'nlist,Nirradiated,NPsurvive, new SFave: ',nlist,Nirradiated,NPsurvive,newSFtot/nlist
        write(*,'(a,i8,f8.2)') 'Average mitosis_time: ',nlist,total_mitosis_time/(3600*nlist)
        write(*,'(a,f8.2,i6)') 'CA_time_h, NpreCA: ',CA_time_h,NpreCA
    endif
    if (NPsurvive > 0) then
        SFave = SFtot/(NPsurvive + Napop)
    else
        SFave = 0
    endif
    write(*,'(a,i6,2e12.3)') 'NPsurvive,SFlive,SFtot: ',NPsurvive,SFlive,SFtot
    write(*,*)
    write(nflog,'(a,f8.2)') 'CA_time_h: ',CA_time_h
    write(nflog,'(a,2i6)') 'Npsurvive, Napop: ',Npsurvive,Napop
    write(nflog,'(a,e12.4,f8.3)') 'SFave,log10(SFave): ',SFave,log10(SFave)
    write(*,'(a,e12.4,f8.3)') 'SFave,log10(SFave): ',SFave,log10(SFave)
    write(*,*)

    call completed
    res = 1
endif

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getNviable
integer :: Nlive(MAX_CELLTYPES)
integer :: kcell, ityp, idrug, nd
logical :: tag
type(cell_type), pointer :: cp

Nviable = 0
Nlive = 0
nd = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) then
	    nd = nd+1
	    cycle
	endif
    ityp = cp%celltype
    Nlive(ityp) = Nlive(ityp) + 1
	if (cp%state == DYING) cycle
	Nviable(ityp) = Nviable(ityp) + 1
enddo
if (Nlive(1) /= Ncells_type(1)) then
	write(*,'(a,5i8)') 'Error: getNviable: Nlive /= Ncells_type(1), nd, Napop, Nmitotic: ',Nlive(1),Ncells_type(1),nd,Napop, Nmitotic
	write(nflog,'(a,5i8)') 'Error: getNviable: Nlive /= Ncells_type(1), nd, Napop, Nmitotic: ',Nlive(1),Ncells_type(1),nd,Napop, Nmitotic
	stop
endif
if (Nviable(1) /= Ncells_type(1) - Ndying(1)) then
	write(*,'(a,4i8)') 'Error: getNviable: Nviable /= Ncells_type(1) - Ndying, Nmitotic: ',Nviable(1),Ncells_type(1),Ndying(1),Nmitotic
	write(nflog,'(a,4i8)') 'Error: getNviable: Nviable /= Ncells_type(1) - Ndying, Nmitotic: ',Nviable(1),Ncells_type(1),Ndying(1),Nmitotic
	stop
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
    if (cp%state /= DEAD .and. cp%psurvive < 0) then
        n = n+1
        write(*,'(a,i6,i3,3f8.4)') 'nondivided: kcell, phase, kt2cc,ke2cc,kcc2a: ',kcell, cp%phase,cp%kt2cc,cp%ke2cc,cp%kcc2a
        write(*,*) 'fg: ',cp%fg
        write(nflog,'(a,i6,i3,3f8.4)') 'nondivided: kcell, phase, kt2cc,ke2cc,kcc2a: ',kcell, cp%phase,cp%kt2cc,cp%ke2cc,cp%kcc2a
        write(nflog,*) 'fg: ',cp%fg
    endif
enddo
write(*,*) 'Total nondivided: ',n
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getSFlive(SF)
real(REAL_KIND) :: SF
real(REAL_KIND) :: sfsum, totDSB, total, total0
integer :: kcell, n
type(cell_type), pointer :: cp

n = 0
sfsum = 0
total = 0
total0 = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) cycle
    if (cp%totDSB0 <= 0) then
        cycle  ! this cell was not irradiated - must be a daughter cell (can't happen)
    endif
    n = n+1
    sfsum = sfsum + cp%Psurvive
    totDSB = sum(cp%DSB)
    total = total + totDSB
    total0 = total0 + cp%totDSB0
    if (cp%Psurvive < 0) then
        write(*,*) 'ERROR: getSFlive: Psurvive: ',kcell,cp%Psurvive
        stop
    endif
enddo
SF = sfsum/n
end subroutine

!-----------------------------------------------------------------------------------------
! The selection of outputs needs to be made an input parameter.
! Choices:
! SF + distribution
! SF only
! distribution only
!
! The first two can be lumped together - if distribution is not wanted it will not be read.
! I.e. need only use_SF
! For now assume always use_SF = true, because fitting without SF is no good.
! For PEST runs, use_SF corresponds to 'M' fitting
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
        if (cp%Psurvive > sfmax) then
            sfmax = cp%Psurvive
            kcellmax = kcell
        endif
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
write(nflog,'(a,6f12.3)') 'totPmit, totPaber, tottotDSB: ',totPmit, totPaber, tottotDSB
write(*,'(a,i6,5f11.1)') 'Nmitosis, totPmit, totPaber, tottotDSB: ',int(Nmitosis),totPmit, totPaber, tottotDSB
write(*,'(a,6e12.3)') 'totPaber: ',totPaber
write(nflog,'(a,4f10.5)') 'SFtot_phase: ',SFtot_phase
write(nflog,'(a,i6,f10.5)') 'Nmitosis,SFtot: ',Nmitosis,sum(SFtot_phase)
! adjust for pre-rep doubling of misjoins
totNmisjoins(1) = 2*totNmisjoins(1)
write(nflog,'(a,7f9.3)') 'Ave (pre, post) NDSB, Nmisjoins: ', &
    totNDSB/nmitosis,totNmisjoins/nmitosis,sum(totNmisjoins)/nmitosis
write(*,'(a,7f9.3)') 'Ave (pre, post) NDSB, Nmisjoins: ', &
    totNDSB/nmitosis,totNmisjoins/nmitosis,sum(totNmisjoins)/nmitosis

#if 0
if (.false.) then   ! make this true to write BBB lines
    ! Averages
    SFMave = 0
    ave = 0
    do kcell = 1,nlist
        cp => cell_list(kcell)
        SFMave = SFMave + cp%Psurvive
        do i = 1,3
            do j = 1,2
                k = (i-1)*2 + j
                ave(k) = ave(k) + cp%DSB0(i,j)
            enddo
        enddo
        ave(7) = ave(7) + cp%t_mitosis
    enddo
    SFMave = SFMave/nlist
    ave = ave/nlist
endif
#endif

99 continue
if (use_synchronise) call G2_time_distribution()
if (use_PEST) then
    if (use_SF) then
        write(nfres,'(e15.6)') log10(SFave)
    endif
endif

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

#if 0
!-----------------------------------------------------------------------------------------
! This is for drm_monolayer_deriv.exe
!-----------------------------------------------------------------------------------------
subroutine getResults(SF, dist)
!DEC$ ATTRIBUTES DLLEXPORT :: getResults
real(REAL_KIND) :: SF
integer :: dist(:)

SF = SFave
dist = phase_dist
end subroutine
#endif

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Execute(infile_array,inbuflen,outfile_array,outbuflen,res) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char), intent(in) :: infile_array(*), outfile_array(*)
integer(c_int) :: inbuflen, outbuflen, res
character*(2048) :: infile, outfile, logfile
character*(13) :: fname
character*(1) :: numstr
logical :: ok, success, isopen
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

res = 0
write(nflog,*) 'inputfile:  ', trim(infile)
write(nflog,*) 'outputfile: ', trim(outfile)

DELTA_T = 600
nsteps = 100
res=0

call Setup(infile,outfile,ok)
if (.not. ok) then
	write(nflog,*) '=== Setup failed ==='
endif
if (ok) then
	res = 0
else
    write(nflog,*) 'Setup error'
	res = 1
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
if (allocated(gaplist)) deallocate(gaplist,stat=ierr)

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)
inquire(nfphase,OPENED=isopen)
if (isopen) close(nfphase)

if (par_zig_init) then
	call par_zigfree
endif
end subroutine

end module
