module mc
use real_kind_mod
use global
use eta_module
use km10_mod

implicit none

! There are 4 pathways:
integer, parameter :: NHEJfast = 1
integer, parameter :: NHEJslow = 2
integer, parameter :: HR = 3
integer, parameter :: TMEJ = 4  ! this is alt-EJ
character*(2) :: phaseName(4) = ['G1','S ','G2','M ']

real(8) :: Pcomplex = 0.4337    ! fraction of created DSBs that are complex (McMahon: complexFrac)
real(8) :: pJeggo = 0.9         ! fraction of post-replication DSBs that are fast. (NOT USED)
real(8) :: Kcoh                 ! cohesin effect (read in input file)
real(8) :: pHRs_max, pHRc_max, pHRs_G2, pHRc_G2
real(8) :: Preass               ! rate of reassignment to pathway 4, TMEJ (prob of reass/hour)
logical :: use_sigmoid = .true.
real(8) :: rmin = 0.1, kdecay = 0.1, ksig = 1, csig = 8.56    ! decay function parameters  TO BE RENAMED
real(8) :: sigma_NHEJ = 0.04187
real(8) :: sigma_TMEJ = 0.08
real(8) :: repRate(NP)          ! repair rate parameters
real(8) :: baseRate
real(8) :: mitRate(2)           ! (McMahon: mitoticRate) TO BE RENAMED
real(8) :: Klethal = 0.4

real(8) :: KATM1G1, KATM2G1     ! KATM parameters for G1 CP slowdown
real(8) :: KATM1G1M, KATM2G1M   ! KATM parameters for post-mitosis G1 CP slowdown
real(8) :: KATM1S, KATM2S       ! KATM parameters for S CP slowdown
real(8) :: KATR1S, KATR2S       ! KATR parameters for S ATR_act activation (when ATR_in_S = 2)

! DNA-PK inhibition parameter
real(8) :: Chalf    ! inhibitor concentration that halves repair rate 

! Jaiswal formulation (26/09/22)
real(8) :: Kcc2a, Kcc2e, Kd2e, Kd2t, Ke2cc, Kt2cc, Kti2t
real(8) :: Km1, Km10, Km10t
real(8) :: kmmp, kmmd               ! ATM production, decay
real(8) :: kmrp, kmrd               ! ATR production, decay
real(8) :: kmccp, kmccrd, kmccmd    ! CC production, ATR-driven decay, ATM-driven decay
real(8) :: CC_tot, ATR_tot, ATM_tot, CC_act0, CC_threshold, CC_threshold_factor, norm_factor
real(8) :: G1_tdelay                ! delay before ATM_act updated (hours)
logical :: use_Jaiswal = .true.
logical :: vary_km10 = .true.
real(8) :: jaiswal_std = 0.6

real(8) :: repRateFactor(NP)

real(8) :: totNmisjoins(2), totPaber(2), totPmit(2)  ! 1 = pre-rep, 2 = post-rep
real(8) :: tottotDSB, totNDSB(2)

logical :: use_addATMATRtimes = .false.
logical :: use_G2_pATM_Nindependent = .false.
logical :: output_DNA_rate = .false.
logical :: negligent_G2_CP = .false.
logical :: use_DSB_CP = .false.
logical :: use_D_model = .false.

logical :: use_exp_slowdown = .false.
logical :: use_G1_stop = .false.    ! These flags control use of either CP delay (true) or slowdown (false)
logical :: use_S_stop = .false.
logical :: use_G1_pATM = .false.
logical :: use_S_pATM = .false.

logical :: use_G2_stop = .false.                        ! because use_Jaiswal is true
logical :: use_phase_dependent_CP_parameters = .true.   ! now always true

! Normalisation
real(8) :: control_ave(4)   ! now set equal to ccp%f_G1, ...
logical :: normalise, M_only
character(6) :: expt_tag

! Iliakis
real(8) :: nIliakis
real(8) :: kIliakis     ! if kIliakis = 0, fIliakis = 1
real(8) :: fIliakis
logical :: use_Iliakis 

real(8) :: next_write_time

logical, parameter :: drop_mitotic_cells = .true.
logical, parameter :: write_nfres = .false.

contains

!--------------------------------------------------------------------------
! Note: decision about PEST output is based on expt_ID. 
!--------------------------------------------------------------------------
subroutine ReadMcParameters(nfin)
integer :: nfin
integer :: iuse_baserate, iuse_exp, expt_ID, nCPparams, iph, j
real(8) :: TMEJrep, TMEJfid, SSArep, SSAfid
real(8) :: pHR_S, pfc_S, pHR_G2, pfc_G2, k3, k4

write(*,*) 'ReadMcParameters:'
read(nfin,*) expt_ID
read(nfin,*) CA_time_h
read(nfin,*) baseRate
read(nfin,*) mitRate(1)
read(nfin,*) mitRate(2)
read(nfin,*) Klethal

read(nfin,*) KATM1G1
read(nfin,*) KATM2G1
read(nfin,*) KATM1G1M
read(nfin,*) KATM2G1M
read(nfin,*) KATM1S
read(nfin,*) KATM2S
read(nfin,*) KATR1S
read(nfin,*) KATR2S

read(nfin,*) repRate(NHEJfast)
read(nfin,*) repRate(NHEJslow)
read(nfin,*) repRate(HR)
read(nfin,*) repRate(TMEJ)
read(nfin,*) pComplex
read(nfin,*) Kcoh
read(nfin,*) pHRs_max
read(nfin,*) pHRc_max
read(nfin,*) rmin
read(nfin,*) ksig
read(nfin,*) csig
read(nfin,*) nIliakis
read(nfin,*) kIliakis
read(nfin,*) G1_tdelay
read(nfin,*) Chalf
read(nfin,*) Preass
read(nfin,*) dsigma_dt
read(nfin,*) sigma_NHEJ

read(nfin,*) Kcc2a
read(nfin,*) Kcc2e
read(nfin,*) Kd2e
read(nfin,*) Kd2t
read(nfin,*) Ke2cc
read(nfin,*) Kt2cc
read(nfin,*) Kti2t
read(nfin,*) Kmccp
read(nfin,*) Kmccmd
read(nfin,*) Kmccrd
read(nfin,*) Kmrp
read(nfin,*) Kmrd
read(nfin,*) Kmmp
read(nfin,*) Kmmd
read(nfin,*) CC_tot
read(nfin,*) ATR_tot
read(nfin,*) ATM_tot
read(nfin,*) CC_act0
read(nfin,*) CC_threshold_factor
read(nfin,*) norm_factor

CC_threshold = CC_threshold_factor*CC_tot

use_Iliakis = (kIliakis > 0)

use_SF = .false.
nphase_hours = 0
next_phase_hour = 0
phase_hour(:) = 0
output_DNA_rate = .false.
normalise = .false.
M_only = .false.
!if (expt_ID == -1) then    ! PDSN dose = 0
!    expt_tag = "PDSN0G"
!    compute_cycle = .true.
!    normalise = .true.
!    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
!    nphase_hours = 2
!    next_phase_hour = 1
!    phase_hour(1:2) = [5.0, 11.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
!    ! Note: output G1, S, G2, M
!elseif (expt_ID == -2) then    ! PDSN dose = 2
!    expt_tag = "PDSN2G"
!    compute_cycle = .true.
!    normalise = .true.
!    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
!    next_phase_hour = 1
!!    nphase_hours = 6
!!    phase_hour(1:6) = [1.0, 2.0, 3.0, 5.0, 8.5, 11.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
!    nphase_hours = 25
!    do j = 1,25
!        phase_hour(j) = j
!    enddo
!!    phase_hour(1:5) = [5.0, 8.5, 11.5, 18.0, 25.0]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
!elseif (expt_ID == -3) then    ! PDSN dose = 6
!    expt_tag = "PDSN6G"
!    compute_cycle = .true.
!    normalise = .true.
!    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
!    nphase_hours = 3
!    next_phase_hour = 1
!    phase_hour(1:3) = [5.0, 8.5, 11.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
if (expt_ID == 1) then
    use_SF = .true.     ! in this case SFave only is recorded
    compute_cycle = .false.

elseif (mod(expt_ID,10) == 2) then    ! this is the compute_cycle case for CA-135
    expt_tag = "CA-135"
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 5    !5
    next_phase_hour = 1
    phase_hour(1:5) = [5.0, 8.5, 11.5, 18.5, 24.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
    ! Note: output G1, S, G2, M
elseif (mod(expt_ID,10) == 9) then    ! this is the compute_cycle case for CC-13
    expt_tag = "CC-13"
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 8    !8
    next_phase_hour = 1
    phase_hour(1:8) = [1, 2, 3, 5, 8, 12, 18, 24]   
    ! Note: output M
elseif (mod(expt_ID,10) == 5) then    ! this is the compute_cycle case for CC-11
    expt_tag = "CC-11"
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 5
    next_phase_hour = 1
    phase_hour(1:5) = [1.0, 1.5, 2.5, 3.5, 4.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
    ! Note: output G1, S, G2

elseif (expt_ID == 6) then    ! this is the output_DNA_rate case (EDU)
    compute_cycle = .false.
    output_DNA_rate = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 18
    do j = 1,nphase_hours
        phase_hour(j) = j*0.25
    enddo
    next_phase_hour = 1
elseif (expt_ID == 3) then
    use_SF = .true.     ! in this case SFave is recorded and there are multiple phase distribution recording times
    nphase_hours = 4
    next_phase_hour = 1
    phase_hour(1:4) = [8, 12, 18, 24]
elseif (expt_ID == 4) then    ! this is the synchronised IR case
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 1
    next_phase_hour = 1
    phase_hour(1:5) = [40, 0, 0, 0, 0]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
elseif (mod(expt_ID,10) == 7) then    ! this is the compute_cycle case for multiple times, no PEST
    compute_cycle = .true.
    normalise = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 49
    next_phase_hour = 1
    do j = 1,49
        phase_hour(j) = (j-1)*0.5
    enddo
elseif (mod(expt_ID,10) == 8) then    ! this is the compute_cycle case for M%-only experiments
    compute_cycle = .true.
    use_SF = .false.    ! in this case no SFave is recorded, there are multiple phase distribution recording times
    nphase_hours = 5
    next_phase_hour = 1
    phase_hour(1:5) = [1.0, 1.5, 2.5, 3.5, 4.5]   ! these are hours post irradiation, incremented when irradiation time is known (in ReadProtocol)
else
    if (use_PEST) then
        write(*,*) 'Error: ReadMcParameters: with PEST expt_ID must be 1 - 8'
        write(nflog,*) 'Error: ReadMcParameters: with PEST expt_ID must be 1 - 8'
        stop
    endif
endif
write(*,*) 'nphase_hours: ',nphase_hours
if (nphase_hours > 0) then
    write(*,*) 'phase_hour:'
    write(*,'(10f6.2)') phase_hour(1:nphase_hours)
endif
write(*,*) 'did ReadMcParameters:'

end subroutine

!--------------------------------------------------------------------------
! f_S is the fractional progress through S-phase
! DSB0[:,:] is the initial number of DSBs on each repair pathway
! Distinction between pre- and post-DNA replication DSBs.
! pre = DSB(:,1), post = DSB(:,2)
! totDSB0 is the total initial DSB count
!--------------------------------------------------------------------------
subroutine cellIrradiation(cp, dose)
type(cell_type), pointer :: cp
real(8) :: dose
integer :: ityp, kpar = 0
integer :: phase
real(8) :: DSB0(NP,2)
real(8) :: totDSB0, baseDSB, fin, T_S, T_G2,f_S, NG1, NNHEJ, pHRs, pHRc, pHR, pNHEJ, NHRs, NHRc
real(8) :: Pbase, Pdie, R
real(8) :: DSB_Gy, L    ! = 35
real(8) :: th, Npre, Npre_s, Npre_c, Npost, Npost_s, Npost_c, Pc, x, fstart 
type(cycle_parameters_type),pointer :: ccp
logical :: use_Poisson_DSB = .true.
integer :: kcell

if (use_Poisson_DSB .AND. use_no_random) use_Poisson_DSB = .false.
cp%irradiated = .true.
next_write_time = 0
ccp => cc_parameters(1)
phase = cp%phase
cp%phase0 = min(phase, M_phase)
L = exp(-35.d0)
if (use_Poisson_DSB .and. Ncells > 1) then
    DSB_Gy = poisson_gen(L)
else
    DSB_Gy = 35
endif
NG1 = DSB_Gy*dose
DSB0 = 0

T_S = ccp%T_S*cp%fg(2)
T_G2 = ccp%T_G2*cp%fg(3)
th = 0
fstart = 0
if (constant_S_pHR) fstart = 0.8
If (phase == G1_phase) Then
    f_S = 0
ElseIf (phase == S_phase) Then
    f_S = cp%progress
    th = max(0.0d0,(cp%progress - fstart)*T_S/3600)
ElseIf (phase >= G2_phase) Then
    f_S = 1
    th = ((1.0 - fstart)*T_S + cp%progress*T_G2)/3600
End If
totDSB0 = (1 + f_S) * NG1
DSB0(TMEJ,:) = 0
If (phase == G1_phase) Then
    Npre = totDSB0
Else
    Npre = NG1 * (1 - f_S)
End If
Npost = totDSB0 - Npre
If (f_S > 0) Then
    pHRs = fIliakis*pHRs_max * ((1 - rmin) * fdecay(th) + rmin)
    pHRc = fIliakis*pHRc_max * ((1 - rmin) * fdecay(th) + rmin)
    if (single_cell) &
        write(*,'(a,i2,5f10.3)') 'phase, cp%progress, th, fdecay(th), pHRs, pHRc: ', &
                                    phase, cp%progress, th, fdecay(th), pHRs, pHRc
else
    pHRs = 0
    pHRc = 0
End If
pHR = (1 - pComplex)*pHRs + pComplex*pHRc
    DSB0(NHEJfast,1) = (1-pComplex)*Npre
    DSB0(NHEJfast,2) = (1-pComplex)*(1-pHRs)*Npost     ! fast
    DSB0(NHEJslow,1) = pComplex*Npre
    DSB0(NHEJslow,2) = pComplex*(1-pHRc)*Npost     ! slow
    DSB0(HR,1) = 0
    DSB0(HR,2) = pHR*Npost

cp%pATM = 0
cp%pATR = 0
if (phase == G1_phase) then
    ! Apoptosis in G1
    Pbase = exp(-baseRate*totDSB0)   ! this is 1 - Pdie
    Pdie = 1 - Pbase 
    if (single_cell) Pdie = 0
    R = par_uni(kpar)
    if (R < Pdie) then  ! cell dies of apoptosis
        cp%state = DEAD
        Napop = Napop + 1
        Ncells = Ncells - 1
        ityp = cp%celltype
        Ncells_type(ityp) = Ncells_type(ityp) - 1
        Nviable(ityp) = Ncells_type(ityp)
        Ndead(ityp) = Ndead(ityp) + 1
        write(nflog,*) 'apoptotic death: ',kcell_now, phase
    endif
endif

cp%DSB = DSB0
cp%DSB0 = DSB0
cp%totDSB0 = totDSB0
cp%Nmis = 0

totPmit = 0
totPaber = 0
tottotDSB = 0

totNDSB = 0
totNmisjoins = 0

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function fdecay(t) result(f)
real(8) :: t, f

if (use_sigmoid) then
    f = fsigmoid(t)
else
    f = exp(-kdecay*t)
endif
end function

!--------------------------------------------------------------------------
! From Richards's curve
! t is hours since start of S-phase, csig approx = T_S
!--------------------------------------------------------------------------
function fsigmoid(t) result(f)
real(8) :: t, f

f = 1.0/(1.0 + exp(ksig*(t - csig)))
!write(*,'(a,4e12.3)') 'fsigmoid: ksig,csig,t,f: ',ksig,csig,t,f
end function


!------------------------------------------------------------------------
! Only G2 is affected by pATR
! Best to turn off pATR in S-phase by setting katr1s = katr3s = 0
! Now slowdown is acting in G1-phase and S-phase, and only affected by pATM
! fslow in computed for G2, but not used because G2 is controlled by Jaiswal.
! Slowdown must occur only for cells that were irradiated.  
! This means for cells that were present at IR.  In other words:
! Initially, every cell has cp%irradiated = false
! At IR, every cell has cp%irradiated = true
! At cell division, each daughter cell should have cp%irradiated = false
!------------------------------------------------------------------------
subroutine get_slowdown_factors(cp,fATM,fATR)
type(cell_type), pointer :: cp
integer :: iph
real(REAL_KIND) :: fATM, fATR
real(REAL_KIND) :: pATM, pATR, k1, k2, N_DSB, atm, atr
logical :: use_ATR
logical :: OK

use_ATR = (ATR_in_S == 2)
if (.not.cp%irradiated) then
    fATM = 1
    fATR = 1
    return
endif
iph = cp%phase
if (iph > G2_phase) then
    write(*,*) 'Error: get_slowdown_factors called with iph = ',iph
    stop
endif
OK = .true.
if ((iph == G1_phase) .and. use_G1_stop) OK = .false.
if ((iph == S_phase) .and. use_S_stop) OK = .false.
if (.not.OK) then
    write(*,*) 'Error: get_slowdown_factors: iph, use_G1_stop, use_S_stop: ',iph, use_G1_stop, use_S_stop
    stop
endif

atr = 0
fATR = 1
if (iph == G1_phase) then
    if (use_G1_pATM ) then
        atm = cp%pATM
    else
        atm = cp%ATM_act
    endif
elseif (iph == S_phase) then
    if (use_S_pATM ) then
        atm = cp%pATM
    else
        atm = cp%ATM_act
        if (use_ATR) atr = cp%ATR_act
    endif
endif
if (iph == 1) then
    k1 = KATM1G1
    k2 = KATM2G1
elseif (iph == 2) then
    k1 = KATM1S
    k2 = KATM2S
endif
if (iph == G1_phase) then
    if (cp%birthtime > t_irradiation) then      ! post-mitosis
        if (use_SF) then
            write(*,*) 'get_slowdown_factors: should not get here, stopping'
            write(nflog,'(a,2i6,2f8.3,i6)') 'in get_slowdown_factors: kcell,iph,birthtime,t_irrad: ',kcell_now,iph,cp%birthtime/3600,t_irradiation/3600   !,cp%rad_state
            close(nflog)
            stop
        endif
        k1 = KATM1G1M
        k2 = KATM2G1M
    endif
endif
    
if ((k2 + atm) > 0) then
    fATM = max(0.01,1 - k1*atm/(k2 + atm))
else
    fATM = 1.0
endif
if (iph == S_phase .and. use_ATR) then
    k1 = KATR1S
    k2 = KATR2S
    fATR = max(0.01,1 - k1*atr/(k2 + atr))
endif
end subroutine

!------------------------------------------------------------------------
! Phase slowdown is handled here. 
!------------------------------------------------------------------------
function slowdown(cp) result(fslow)
type(cell_type), pointer :: cp
real(REAL_KIND) :: fslow
type(cycle_parameters_type),pointer :: ccp
integer :: iph
real(REAL_KIND) :: fATM, fATR, dt

ccp => cc_parameters(1)
if (use_phase_dependent_CP_parameters) then
    iph = cp%phase
else
    iph = 1
endif
if ((iph == 1 .and. use_G1_stop) .or. &
    (iph == 2 .and. use_S_stop) .or. &
    (iph == 3 .and. (use_G2_stop .or. use_jaiswal))) then
    fslow = 1.0
else
    dt = DELTA_T
    call get_slowdown_factors(cp,fATM, fATR)
    if (use_addATMATRtimes) then
        dt = DELTA_T
        fslow = max(0.0,fATM + fATR - 1)
    else
        fslow = fATM*fATR
    endif
endif
end function

!------------------------------------------------------------------------
! No DSB repair
!------------------------------------------------------------------------
subroutine test_Jaiswal
real(8) :: t, dth, tsim, dose, T_G2h, R, kfactor, DSB0(NP,2), kmccp_temp, Kcc2a_temp
integer :: it, Nt, i, kpar = 0
type(cell_type), pointer :: cp

tsim = 60
dth = 1/6.0
Nt = int(tsim/dth)
tnow = 0
T_G2h = 3.951   ! use mean cycle time
cp => cell_list(1)
cp%phase = G2_phase
cp%progress = 0.9
cp%t_start_G2 = 0
do i = 1,8
    kmccp_temp = 4.0 + (i-1)*0.5
    Kcc2a_temp = get_Kcc2a(kmccp_temp,CC_tot,CC_threshold_factor,T_G2h)    
    write(*,'(a,2f8.3)') 'kmccp, kcc2a: ',kmccp_temp,kcc2a_temp
enddo
cp%Kcc2a = get_Kcc2a(kmccp,CC_tot,CC_threshold_factor,T_G2h)  
write(*,*) 'Kcc2a: ',cp%Kcc2a
cp%kt2cc = kt2cc    ! use mean values
cp%ke2cc = ke2cc
cp%CC_act = 8.05    ! initialised to G2(0.9)
cp%ATM_act = 0
cp%ATR_act = 0
dose = 2
fIliakis = kIliakis**nIliakis/(kIliakis**nIliakis + (dose-dose_threshold)**nIliakis)    ! normally set in Irradiation()
cp%fg = 1
call cellIrradiation(cp, dose)
DSB0 = cp%DSB
write(*,'(a,3f8.1)') 'DSB0(1:3,1): ',DSB0(1:3,1)
write(*,'(a,3f8.1)') 'DSB0(1:3,2): ',DSB0(1:3,2)
t = 0
write(nflog,*) 'test_Jaiswal: Nt: ',Nt
write(*,*) 'G2: ATR_act, CC_act, t_G2, D_ATR, D_ATM '
do it = 1,Nt
    call Jaiswal_update(cp, dth)
    t = t + dth
    tnow = t*3600
    cp%DSB = DSB0*exp(-t*0.256) ! 0.256 OK for dose = 2 (matches model with 45a parameters)
    if (cp%CC_act > CC_threshold) exit
enddo
stop
end subroutine

!------------------------------------------------------------------------
! Note: DNA repair is happening in updateRepair (pathwayRepair), and since 
! misrepair is computed at the same time, best to leave that as it is.  This
! means that DSB(:) is being updated with the full time step, not using the
! reduced time step employed for the Jaiswal equations.  Therefore damage D
! is assumed to be fixed for this updating.
!------------------------------------------------------------------------
subroutine Jaiswal_update(cp, dth)
type(cell_type), pointer :: cp
real(8) :: dth
real(8) :: dt = 0.001
real(8) :: D_ATR, D_ATM, CC_act, ATR_act, ATM_act, CC_inact, ATR_inact, ATM_inact, tIR
real(8) :: dCC_act_dt, dATR_act_dt, dATM_act_dt, t, t_G2, Kkcc2a, DSB(NP), CC_act0, d(3)
integer :: iph, it, Nt
type(cycle_parameters_type),pointer :: ccp
logical :: use_ATR  ! ATR is used in G2, and computed in S if ATR_in_S >= 1
logical :: dbug

iph = cp%phase
if (iph >= 7) iph = iph - 6     ! checkpoint phase numbers --> phase number, to continue ATM and ATR processes through checkpoints
if (iph > G2_phase) then
    return
endif
use_ATR = (iph == 3) .or. ((iph == 2) .and. (ATR_in_S >= 1))
Nt = int(dth/dt + 0.5)

do it = 1,NP
    DSB(it) = sum(cp%DSB(it,:))     ! add pre and post
enddo
ATR_act = cp%ATR_act
CC_act = cp%CC_act
dbug = (iph == 2 .and. (kcell_now <= 0))
if (iph == G1_phase) then
    D_ATM = DSB(NHEJslow)*norm_factor
    ATM_act = cp%ATM_act
elseif (iph == S_phase) then
    D_ATM = (DSB(HR) + DSB(NHEJslow))*norm_factor
    ATM_act = cp%ATM_act
    if (use_ATR) then
        D_ATR = DSB(HR)*norm_factor
        ATR_act = cp%ATR_act
    endif
elseif (iph == G2_phase) then
    D_ATR = DSB(HR)*norm_factor
    D_ATM = (DSB(HR) + DSB(NHEJslow))*norm_factor
    CC_act = cp%CC_act
    CC_act0 = CC_act
    ATR_act = cp%ATR_act
    ATM_act = cp%ATM_act
    t_G2 = (tnow - cp%t_start_G2)/3600
    if (t_G2 > 40) then     ! force ATR_act to taper to 0 after 30h in G2
        ATR_act = 0
    elseif (t_G2 > 30) then
        ATR_act = ATR_act*(40 - t_G2)/(40 - 30) 
    endif
else
    return
endif
if (use_cell_kcc2a_dependence) then
    Kkcc2a = cp%Kcc2a
else
    Kkcc2a = Kcc2a
endif
do it = 1,Nt
    ATM_inact = ATM_tot - ATM_act
    ATR_inact = ATR_tot - ATR_act
    if (iph == G2_phase) then
        CC_inact = CC_tot - CC_act
        ATR_inact = ATR_tot - ATR_act
        dCC_act_dt = (Kkcc2a + CC_act) * CC_inact / (Kmccp + CC_inact) - cp%Kt2cc * ATM_act * CC_act / (Kmccmd + CC_act) - cp%Ke2cc * ATR_act * CC_act / (Kmccrd + CC_act)
        d(1) = (Kkcc2a + CC_act) * CC_inact / (Kmccp + CC_inact)    ! CC_act effect
        d(2) = - cp%Kt2cc * ATM_act * CC_act / (Kmccmd + CC_act)    ! ATM_act effect
        d(3) = - cp%Ke2cc * ATR_act * CC_act / (Kmccrd + CC_act)    ! ATR_act effect
        dATR_act_dt = Kd2e * D_ATR * ATR_inact / (Kmrp + ATR_inact) - Kcc2e * ATR_act * CC_act / (Kmrd + CC_act)
        
        CC_act = CC_act + dt * dCC_act_dt
        CC_act = max(CC_act, 0.0)
        ATR_act = ATR_act + dt * dATR_act_dt
        ATR_act = min(ATR_act, ATR_tot)
    elseif (use_ATR .and. D_ATR > 0) then
        dATR_act_dt = Kd2e * D_ATR * ATR_inact / (Kmrp + ATR_inact)  - Kcc2e * ATR_act * CC_act / (Kmrd + CC_act)
        ATR_act = ATR_act + dt * dATR_act_dt
        ATR_act = min(ATR_act, ATR_tot)
    endif
        dATM_act_dt = Kd2t * D_ATM * ATM_inact / (Kmmp + ATM_inact) - Kti2t * ATM_act / (Kmmd + ATM_act)   
    ATM_act = ATM_act + dt*dATM_act_dt
    ATM_act = min(ATM_act, ATM_tot)
    t = it*dt
enddo
cp%ATM_act = ATM_act
if (iph == G2_phase) then
    cp%CC_act = CC_act
    cp%ATR_act = ATR_act
    cp%dCC_act_dt = dCC_act_dt
elseif (iph == S_phase .and. use_ATR) then
    cp%ATR_act = ATR_act
endif
t = t_simulation/3600.
tIR = (tnow - t_irradiation)/3600
end subroutine

!--------------------------------------------------------------------------
! To determine repair on a given pathway
! Using parameters from mcradio, dth in hours
!--------------------------------------------------------------------------
subroutine pathwayRepair(path, dth, N0, N)
integer :: path
real(8) :: dth, N0, N, dt
integer :: Nt, it

N = N0
if (dth >= 0) then
    Nt = 1
    dt = dth/Nt
    do it = 1,Nt
        N = N*exp(-repRate(path)*dt*repRateFactor(path))
    enddo
else
    N = 0
endif
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function misrepairRate(initialBreaks, finalBreaks, eta) result(Pmis)
real(8) :: initialBreaks, finalBreaks, eta, Pmis
real(8) :: repairedBreaks, atanNum, atanDen
real(8) :: etamod, etafac = 1.0

etamod = etafac*eta
repairedBreaks = initialBreaks-finalBreaks
if (repairedBreaks < 1E-10) then
    Pmis = 0
    return
endif
atanNum = sqrt(3.0)*etamod*repairedBreaks
atanDen = 2 + etamod*(2*initialBreaks*finalBreaks*etamod + initialBreaks + finalBreaks) 
Pmis = 1 - 2 * atan(atanNum/atanDen) / (repairedBreaks*sqrt(3.0)*etamod)
end function

!------------------------------------------------------------------------
! Moved here from updateRepair
!------------------------------------------------------------------------
subroutine getRepRateFactor(cp)
type(cell_type), pointer :: cp
real(REAL_KIND) :: Cdrug

repRateFactor = 1.0
! Now the effect of inhibition is to reduce the repair rate.
! Chalf is the drug concentration that reduces repair rate by 1/2.  fRR = exp(-k*C/Chalf)
if (Chalf == 0) then
    write(nflog,*) 'ERROR: Chalf = 0'
    write(*,*) 'ERROR: Chalf = 0'
    stop
endif
! Since there is no decay in vivo, can set Cdrug = initial conc, unless washout time has been reached.
if (use_drug_halflife) then
    Cdrug = drug_conc*exp(-Khalflife*(t_simulation - drug_time)/3600)
else
    Cdrug = drug_conc
endif

repRateFactor(1:2) = exp(-0.693*Cdrug/Chalf)
end subroutine

!------------------------------------------------------------------------
! To use parameters from mcradio, need to convert time from secs to hours
! To keep track of pre- and post-rep contributions to misjoins (to cp%Nlethal))
! need: 
! Nmis(2), cp%Nmis(2)
!------------------------------------------------------------------------
subroutine updateRepair(cp,dt)
type(cell_type), pointer :: cp
real(8) :: dt
integer :: phase
real(8) :: DSB(NP,2)
real(8) :: DSB0(NP,2)
real(8) :: totDSB0, totDSB, Pmis, Nmis(2), dmis, dNmis(NP), totDSBinfid0, totDSBinfid
real(8) :: ATR_DSB, ATM_DSB, dth, binMisProb
real(8) :: Cdrug, inhibrate, Nreassign
real(8) :: f_S, eta_NHEJ, eta_TMEJ, tIR, eta, Nrep, sigma
real(8) :: th_since_IR
logical :: pathUsed(NP)
integer :: k, iph, jpp  ! jpp = 1 for pre, = 2 for post
logical :: use_DSBinfid = .true.
real(8) :: DSB_min = 1.0e-3
logical :: use_totMis = .true.      ! was false!
logical :: use_ATM != .false.
logical :: dbug
logical :: use_G1_tdelay = .false.
logical :: do_G1_Jaiswal
logical :: use_constant_V = .false.

if (cp%state == DIVIDED) return
dbug = (kcell_now == 39 .and. istep == 108)
dth = dt/3600   ! hours
use_ATM = .true.
phase = cp%phase
iph = phase
DSB = cp%DSB
call getRepRateFactor(cp)
! Preass is an input parameter = prob of reassignment per hour
if (Preass > 0 .and. phase >= S_phase) then
    do jpp = 1,2
    do k = 1,2  ! only NHEJ pathways
        Nreassign = DSB(k,jpp)*Preass*dth
        DSB(k,jpp) = DSB(k,jpp) - Nreassign
    enddo
    enddo
endif
DSB0 = DSB     ! initial DSBs for this time step

totDSB0 = sum(DSB0)
ATM_DSB = sum(DSB(NHEJslow,:)) + sum(DSB(HR,:))   ! complex DSB
ATR_DSB = sum(DSB(HR,:))

if (iph == G1_phase) then
    if (use_G1_tdelay) then
        ! We want to update Jaiswal for a cell in G1 only if G1_tdelay has elapsed since IR
        ! Current time is tnow (sec).  IR time is t_irradiation
        ! Time since IR is th_since_IR(hours)
        th_since_IR = (tnow - t_irradiation)/3600
        do_G1_Jaiswal = (th_since_IR > G1_tdelay).and.(G1_tdelay > 0)
    else
        ! We want to update Jaiswal for a cell in G1 only if the cell has undergone mitosis post-IR
        do_G1_Jaiswal = (cp%birthtime > t_irradiation)      ! (cp%generation > 1)
    endif
endif
do_G1_Jaiswal = .true.      ! Now do Jaiswal for G1 whether pre- or post-mitosis.  Use special katm parameters for G1 post-mitosis
if (((iph == G1_phase).and.do_G1_Jaiswal).or.(iph >= S_phase)) then
    if (cp%irradiated) then
        call Jaiswal_update(cp,dth)  ! try this
    endif
endif

DSB = 0
do k = 1,3
    do jpp = 1,2
        call pathwayRepair(k, dth, DSB0(k,jpp), DSB(k,jpp))
        if (DSB(k,jpp) < DSB_min) DSB(k,jpp) = 0
    enddo
enddo

if (phase == G1_phase) then
    f_S = 0.0
elseif (phase == S_phase) then
    f_S = cp%progress
elseif (phase >= G2_phase) then
    f_S = 1.0
endif
tIR = (t_simulation - t_irradiation)/3600   ! time since IR, in hours
tIR = max(tIR,0.0)
if (use_constant_V) then
    sigma = sigma_NHEJ + tIR*dsigma_dt
    eta_NHEJ = etafun(1.d0,sigma)
else
    if (use_Arnould) then
        eta_NHEJ = eta_Arnould(phase, f_S, tIR, R_Arnould, sigma_NHEJ, Kcoh)
    else
        eta_NHEJ = eta_lookup(phase, NHEJfast, f_S, tIR) 
    endif
endif

Nmis = 0
! For NHEJ pathways
! pre-rep fraction = (1 - f_S)
! post-rep fraction = f_S
! total with doubling = 2(1 - f_S) + f_S = 2 - f_S
totDSB0 = sum(DSB0(NHEJfast,:)) + sum(DSB0(NHEJslow,:))
totDSB = sum(DSB(NHEJfast,:)) + sum(DSB(NHEJslow,:))
Pmis = misrepairRate(totDSB0, totDSB, eta_NHEJ)
dmis = Pmis*(totDSB0 - totDSB)
if (isnan(dmis)) then
    write(nflog,*) 'dmis is NaN'
    write(nflog,'(a,2f8.2,e12.3)') 'totDSB0, totDSB, eta_NHEJ: ',totDSB0, totDSB, eta_NHEJ
    stop
endif
Nrep = totDSB0 - totDSB

Nmis(1) = Nmis(1) + dmis*(1 - f_S)  ! doubled at mitosis
Nmis(2) = Nmis(2) + dmis*f_S

if (sum(DSB0(TMEJ,:)) > 0) then ! not used currently
    ! For TMEJ pathway
    Pmis = misrepairRate(sum(DSB0(TMEJ,:)), sum(DSB(TMEJ,:)), eta_TMEJ)
    dmis = Pmis*(sum(DSB0(TMEJ,:)) - sum(DSB(TMEJ,:)))
    Nmis(1) = Nmis(1) + dmis*(1 - f_S)
    Nmis(2) = Nmis(2) + dmis*f_S
endif
cp%DSB = DSB
cp%Nmis = cp%Nmis + Nmis

end subroutine

!------------------------------------------------------------------------
! Clonogenic survival probability at first mitosis (McMahon, mcradio)
! cp%phase0 is the cell's phase at IR
! Note that this assumes that cells died of apoptosis in G1 at baseRate
! (see cellIrradiation())
! Now pre-rep and post-rep Nlethal, Nmisjoins, Paber
!------------------------------------------------------------------------
subroutine survivalProbability(cp)
type(cell_type), pointer :: cp
real(8) :: DSB(NP,2), totDSB(2), Nmis(2), Nlethal(2), Paber(2), Pbase, Papop, Pmit(2), Psurv
real(8) :: Nlethal_sum, Paber1_nodouble, Nmistot, tIR,totNmis
integer :: k, jpp, ityp
real(8), parameter :: kmit = 1.0    !0.033  ! Baide et al., D10 = 0.1

DSB = cp%DSB
do jpp = 1,2
    totDSB(jpp) = sum(DSB(:,jpp))
    totNDSB(jpp) = totNDSB(jpp) + totDSB(jpp)
enddo
Nmis = cp%Nmis
totNmisjoins = totNmisjoins + Nmis

do k = 1,2
    Pmit(k) = exp(-mitRate(k)*totDSB(k))
enddo
if (cp%phase0 < M_phase) then   ! G1, S, G2
    Paber(1) = exp(-2*Klethal*Nmis(1))
    Paber1_nodouble = exp(-Klethal*Nmis(1))
    Paber(2) = exp(-Klethal*Nmis(2))
    cp%Psurvive = Pmit(1)*Pmit(2)*Paber(1)*Paber(2)  
    cp%Psurvive_nodouble = Pmit(1)*Pmit(2)*Paber1_nodouble*Paber(2)
    if (single_cell) then
        write(nfres,'(a,6e12.3)') 'totNmisjoins,totNDSB: ', &
        2*totNmisjoins(1),totNmisjoins(2),2*totNmisjoins(1)+totNmisjoins(2),totNDSB,sum(totNDSB)
        write(nflog,'(a,3f8.1,4x,3f8.1)') 'mitosis: DSB, Nmis: ',totDSB,sum(totDSB),2*Nmis(1),Nmis(2),2*Nmis(1)+Nmis(2)
     endif
else    ! M_phase or dividing
    if (drop_mitotic_cells) then
        cp%state = DEAD
        cp%Psurvive = 0
        Nmitotic = Nmitotic + 1
        ityp = cp%celltype
        Ncells_type(ityp) = Ncells_type(ityp) - 1
        Nviable(ityp) = Ncells_type(ityp)
        Ndead(ityp) = Ndead(ityp) + 1
    else
        Paber = 1
        cp%Psurvive = exp(-kmit*sum(totDSB))
    endif
endif
tIR = (tnow - t_irradiation)/3600
totNmis = 2*Nmis(1)+Nmis(2)
if (cp%state /= DEAD) then
    NPsurvive = NPsurvive + 1   ! this is the count of cells for which Psurvive has been computed
    cp%mitosis_time = tnow      ! needed to determine if mitosis occurs before or after CA
    cp%state = DIVIDED      ! this actually means that the cell reached division
endif
totPmit = totPmit + Pmit
totPaber = totPaber + Paber
tottotDSB = tottotDSB + sum(totDSB)
end subroutine

!------------------------------------------------------------------------
! Evaluate SFave at the time of drug washout
!------------------------------------------------------------------------
subroutine washoutSF
type(cell_type), pointer :: cp
real(8) :: DSB(NP,2), Nmis(2), Nlethal(2), Paber(2), Pmit(2), Psurvive
real(8) :: totDSB(2), SFwave, sumDSB(2),sumNmis
integer :: icell, k, jpp

sumDSB = 0
sumNmis = 0
SFwave = 0
do icell = 1,Ncells
    cp => cell_list(icell)
    DSB = cp%DSB
    do jpp = 1,2
        totDSB(jpp) = sum(DSB(:,jpp))
    enddo
    sumDSB = sumDSB + totDSB
    Nmis = cp%Nmis
    sumNmis = sumNmis + 2*Nmis(1) + Nmis(2)
    do k = 1,2
        Pmit(k) = exp(-mitRate(k)*totDSB(k))
    enddo
    Paber(1) = exp(-2*Klethal*Nmis(1))
    Paber(2) = exp(-Klethal*Nmis(2))
    Psurvive = Pmit(1)*Pmit(2)*Paber(1)*Paber(2)
    SFwave = SFwave + Psurvive
enddo
SFwave = SFwave/Ncells
if (single_cell) then
    write(nflog,'(a,3f8.1,4x,3f8.1)') 'washout: DSB, Nmis: ',totDSB,sum(totDSB),2*Nmis(1),Nmis(2),2*Nmis(1)+Nmis(2)
else
    write(nflog,'(a,2e12.3)') 'SFwave, log10(SFwave): ',SFwave,log10(SFwave)
    write(nflog,'(a,3f11.2)') 'aveDSB, aveNmis: ',sumDSB/Ncells, sumNmis/Ncells
endif
end subroutine

end module