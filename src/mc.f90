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

! G1 checkpoint
logical :: use_G1_CP_factor = .false.
real(8) :: G1_CP_factor = 0.0
real(8) :: G1_CP_time

! Checking TMEJ misjoining
real(8) :: misjoins(2)      ! 1 = NHEJ, 2 = TMEJ

! Checking prob of slow repair in G2
logical :: check_G2_slow = .true.
integer :: nslow_sum
real(8) :: pHR_sum, pNHEJslow_sum, fdecay_sum

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
if (use_equal_mitrates) mitrate(2) = mitrate(1)
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

! Initialise misjoining counts
misjoins = 0

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
integer, parameter :: option = 2
type(cycle_parameters_type),pointer :: ccp
logical :: use_Jeggo = .true.
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
if (kcell_now <= 10) write(nflog,'(a,i6,f8.3)') 'cellIrradiation: kcell, DSB_Gy: ',kcell_now,DSB_Gy
NG1 = DSB_Gy*dose
DSB0 = 0

T_S = ccp%T_S*cp%fg(2)
T_G2 = ccp%T_G2*cp%fg(3)
th = 0
if (use_Jeggo) then     ! using this
    fstart = 0
    if (constant_S_pHR) fstart = 0.8
    If (phase == G1_phase) Then
        f_S = 0
    ElseIf (phase == S_phase) Then
        f_S = cp%progress
        th = max(0.0d0,(cp%progress - fstart)*T_S/3600)
    ElseIf (phase >= G2_phase) Then
        f_S = 1
        th = ((1.0 - fstart)*T_S + cp%progress*T_G2)/3600     ! This was an error, had T_S (undefined), changed again to replace ccp%T_S by T_S
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
    if (option == 1) then
        pNHEJ = 1 - pHR
        DSB0(NHEJfast,1) = (1 - pComplex) * Npre
        DSB0(NHEJfast,2) = pJeggo * pNHEJ * Npost     ! fast
        DSB0(NHEJslow,1) = pComplex * Npre
        DSB0(NHEJslow,2) = (1 - pJeggo) * pNHEJ * Npost     ! slow
        DSB0(HR,1) = 0
        DSB0(HR,2) = pHR * Npost
        if (phase == G2_phase) then
            nslow_sum = nslow_sum + 1
            pHR_sum = pHR_sum + pHR
            pNHEJslow_sum = pNHEJslow_sum + (1 - pJeggo) * pNHEJ
            fdecay_sum = fdecay_sum + fdecay(th)
        endif
    elseif (option == 2) then   ! USING THIS
        DSB0(NHEJfast,1) = (1-pComplex)*Npre
        DSB0(NHEJfast,2) = (1-pComplex)*(1-pHRs)*Npost     ! fast
        DSB0(NHEJslow,1) = pComplex*Npre
        DSB0(NHEJslow,2) = pComplex*(1-pHRc)*Npost     ! slow
        DSB0(HR,1) = 0
        DSB0(HR,2) = pHR*Npost
!if (kcell_now <= 10) then
!    write(nflog,'(a,2i4,4f8.3)') 'kcell,phase, Npre,Npost,totals: ',kcell_now,phase, Npre,Npost,sum(DSB0(:,:)),sum(DSB0(1:2,1:2))+pHRs*Npre+pHRc*Npost
!    write(nflog,'(a,5f8.3)') 'pHRs,pHRc,pHRs*Npre,pHRc*Npost,pHR*Npost: ',pHRs,pHRc,pHRs*Npre,pHRc*Npost,pHR*Npost
!endif
    endif
else    ! not using this
    if (phase == G1_phase) then
        f_S = 0.0
        totDSB0 = NG1
        DSB0(NHEJslow,1) = Pcomplex*totDSB0
        DSB0(NHEJfast,1) = (1 - Pcomplex)*totDSB0
    elseif (phase == S_phase) then
        th = (istep*DELTA_T - cp%t_S_phase)/3600  ! time since start of S_phase (h)
        f_S = cp%progress
        totDSB0 = (1+f_S)*NG1
        pHRs = pHRs_max*((1-rmin)*fdecay(th) +rmin)
        pHRc = pHRc_max*((1-rmin)*fdecay(th) +rmin)
        Npre = NG1*(1 - f_S)
        Npost = NG1*2*f_S
        Npre_s = Npre*(1 - pComplex)
        Npre_c = Npre*pComplex
        Npost_s = Npost*(1 - pComplex)
        Npost_c = Npost*pComplex
        NHRs = Npost_s*pHRs
        NHRc = Npost_c*pHRc
        DSB0(HR,1) = NHRs + NHRc
        DSB0(NHEJfast,1) = Npre_s + Npost_s - NHRs  ! not correct
        DSB0(NHEJslow,1) = Npre_c + Npost_c - NHRc
    elseif (phase >= G2_phase) then
        f_S = 1
        th = (istep*DELTA_T - cp%t_S_phase)/3600  ! time since start of S_phase (h)
        totDSB0 = 2*NG1
        Npost = totDSB0
        pHRs = pHRs_max*((1-rmin)*fdecay(th)+rmin)
        DSB0(HR,2) = Npost*pHRs
        DSB0(NHEJfast,2) = Npost*(1 - pHRs)
        DSB0(NHEJslow,2) = 0
    endif
endif
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
!elseif (phase == G2_phase) then
!    if (use_Jaiswal) then
!        ! nothing needed to be done here
!    elseif (use_D_model) then   ! currently false
!        cp%pATM = K_ATM(3,1)*totDSB0/(K_ATM(3,2) + totDSB0)
!        cp%pATR = K_ATR(3,1)*totDSB0/(K_ATR(3,2) + totDSB0)
!        if (cp%pATR > 0) then
!            write(*,*) 'cp%pATR: ', kcell_now,K_ATR(3,1),cp%pATR
!            stop
!        endif
!    endif
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

#if 0
!--------------------------------------------------------------------------
! NOT USED
! Jaiswal is used to compute ATM_act
!--------------------------------------------------------------------------
subroutine updateATM(iph,pATM,ATM_DSB,dth)
integer :: iph
real(8) :: pATM, ATM_DSB, dth, pATM0
real(8) :: r, t, fz, kdecay
logical :: OK

OK = .true.
if (iph == G1_phase .and. .not.use_G1_pATM) OK = .false.
if (iph == S_phase .and. .not.use_S_pATM) OK = .false.
if (iph >= G2_phase) OK = .false.
if (.not.OK) then
    write(*,*) 'ERROR: updateATM: should not get here with iph, use_G1_pATM, use_S_pATM: ',iph,use_G1_pATM,use_S_pATM
    stop
endif

pATM0 = pATM
if (use_G2_pATM_Nindependent .and. iph == G2_phase) then
    r = K_ATM(iph,1)   ! rate of production of pATM
else
    r = K_ATM(iph,1)*ATM_DSB   ! rate of production of pATM
endif
t = (tnow - t_irradiation)/3600
kdecay = max(1.0e-8,K_ATM(iph,2))  ! decay rate constant
pATM = r/kdecay + (pATM - r/kdecay)*exp(-kdecay*dth)
if (isnan(pATM)) then
    write(*,*) 'NAN in updateATM: ',iph, r, kdecay, r/kdecay
    stop
endif
end subroutine

!------------------------------------------------------------------------
! NOT USED
! since use_G1_stop = false.  Using slowdown factors
!------------------------------------------------------------------------
function G1_checkpoint_time(cp) result(t)
type(cell_type), pointer :: cp
real(REAL_KIND) :: t, th, atm
integer :: iph = 1

if (.not.use_G1_stop) then
    write(*,*) 'Error: G1_checkpoint_time: should not get here with use_G1_stop = ',use_G1_stop
    stop
endif
if (use_DSB_CP) then
    t = 0
    return
endif
if (use_G1_CP_factor) then
    t = G1_CP_time
    return
endif
if (use_G1_pATM) then
    atm = cp%pATM
else
    atm = cp%ATM_act
endif
th = K_ATM(iph,3)*(1 - exp(-K_ATM(iph,4)*atm))
if (isnan(th)) then
    write(*,*) 'NAN in G1_checkpoint_time: ',atm
    stop
endif
t = 3600*th
end function

!------------------------------------------------------------------------
! Combined effect of ATM_act (or pATM) and ATR_act (was pATR) in S
! Time computed in hours, returned in secs
! CP delay is the sum of delays created by ATM_act and by ATR_act
! NOT USED
! since use_S_stop = false.  Using slowdown factors
!------------------------------------------------------------------------
function S_checkpoint_time(cp) result(t)
type(cell_type), pointer :: cp
real(REAL_KIND) :: t, th, th_ATM, th_ATR, atm, atr
integer :: iph = 2
logical :: use_ATR

use_ATR = (ATR_in_S == 2)
if (.not.use_S_stop) then
    write(*,*) 'Error: S_checkpoint_time: should not get here with use_S_stop = ',use_S_stop
    stop
endif
if (use_S_pATM) then
    atm = cp%pATM
else
    atm = cp%ATM_act
endif
atr = cp%ATR_act
th_ATM = K_ATM(iph,3)*(1 - exp(-K_ATM(iph,4)*atm))  
if (use_ATR) then
    th_ATR = K_ATR(iph,3)*(1 - exp(-K_ATR(iph,4)*atr))
else
    th_ATR = 0
endif
th = th_ATM + th_ATR   
t = 3600*th
end function

!------------------------------------------------------------------------
! Combined effect of pATM and pATR in G2
! Time computed in hours, returned in secs
! CP delay is the sum of delays created by pATM and by pATR
! NOT USED
! G2 progression is handled by Jaiswal.
!------------------------------------------------------------------------
function G2_checkpoint_time(cp) result(t)
type(cell_type), pointer :: cp
real(REAL_KIND) :: t, th, th_ATM, th_ATR, N_DSB, k3, k4
integer :: iph = 3
real(REAL_KIND) :: G2_DSB_threshold = 15

if (use_Jaiswal) then
    write(*,*) 'ERROR: G2_checkpoint_time not used with Jaiswal'
    stop
endif
if (use_D_model) then   ! currently false
    th_ATM = cp%pATM
!    th_ATR = cp%pATR
    th_ATR = 0  ! if effect of pATR is disabled
else
    if (use_DSB_CP) then    ! currently false
        N_DSB = sum(cp%DSB)
        k3 = K_ATR(iph,3)
        k4 = K_ATR(iph,4)
        th = k3*N_DSB/(k4 + N_DSB)
        if (kcell_now < 10) write(*,'(a,4f8.3)') 'N_DSB,k3,k4,th: ',N_DSB,k3,k4,th
        t = 3600*th
        return
    endif
    if (negligent_G2_CP .and. (sum(cp%DSB) < G2_DSB_threshold)) then
        th_ATM = 0
    else
        th_ATM = K_ATM(iph,3)*(1 - exp(-K_ATM(iph,4)*cp%pATM))
    endif
    th_ATR = K_ATR(iph,3)*(1 - exp(-K_ATR(iph,4)*cp%pATR))
endif
th = th_ATM + th_ATR
totG2delay = th + totG2delay
nG2delay = nG2delay + 1
t = 3600*th
end function
#endif

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
!if (use_DSB_CP) then    ! currently false
!    N_DSB = sum(cp%DSB)
!    k3 = K_ATM(iph,3)
!    k4 = K_ATM(iph,4)
!    fATM = max(0.01,1 - k3*N_DSB/(k4 + N_DSB))
!    return
!endif
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
!k3 = K_ATM(iph,3)
!k4 = K_ATM(iph,4)
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
        !k3 = KATM3G1M
        !k4 = KATM4G1M
        k1 = KATM1G1M
        k2 = KATM2G1M
    endif
endif
    
!if ((k4 + atm) > 0) then
!    fATM = max(0.01,1 - k3*atm/(k4 + atm))
!else
!    fATM = 1.0
!endif
if ((k2 + atm) > 0) then
    fATM = max(0.01,1 - k1*atm/(k2 + atm))
else
    fATM = 1.0
endif
if (iph == S_phase .and. use_ATR) then
    !k3 = K_ATR(iph,3)
    !k4 = K_ATR(iph,4)
    !fATR = max(0.01,1 - k3*atr/(k4 + atr))
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

#if 0
!------------------------------------------------------------------------
!------------------------------------------------------------------------
subroutine get_CP_delay(cp)
type(cell_type), pointer :: cp

write(*,*) 'get_CP_delay: phase: ',cp%phase
if (cp%phase == G1_phase) then
    cp%CP_delay = G1_checkpoint_time(cp)
elseif (cp%phase == S_phase) then
    cp%CP_delay = S_checkpoint_time(cp)
elseif (cp%phase == G2_phase .and. .not.use_Jaiswal) then
    cp%CP_delay = G2_checkpoint_time(cp)
endif
end subroutine
#endif

!------------------------------------------------------------------------
! No DSB repair
!------------------------------------------------------------------------
subroutine test_Jaiswal
real(8) :: t, dth, tsim, dose, T_G2h, R, kfactor, DSB0(NP,2), kmccp_temp, Kcc2a_temp
integer :: it, Nt, i, kpar = 0
type(cell_type), pointer :: cp

tsim = 60
dth = 1/6.0     !was 0.1, 1/6 is usual model time step
Nt = int(tsim/dth)
!Nt = 1
tnow = 0
T_G2h = 3.951   ! use mean divide time
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
!R = par_uni(kpar)
!kfactor = 1 + (R - 0.5)*jaiswal_std
!cp%kt2cc = kt2cc*kfactor
!R = par_uni(kpar)
!kfactor = 1 + (R - 0.5)*jaiswal_std
!cp%ke2cc = ke2cc*kfactor
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
!    write(*,'(2f8.4)') t,cp%CC_act
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
!    write(*,'(a,i6,2e12.3)') 'Jaiswal_update: D_ATM,ATM_act: ',kcell_now,D_ATM,ATM_act
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
        ATR_act = ATR_act*(40 - t_G2)/(40 - 30) ! fixed (was (t_G2 - 30)/(40 - 30))
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
        
! Try this
!        dCC_act_dt = max(dCC_act_dt,0.0)
        CC_act = CC_act + dt * dCC_act_dt
!        write(nflog,'(a,i6,4f8.3)') 'it,CC_act,dCC_act_dt,d(1:2): ',it,CC_act,dCC_act_dt, d(1:2)
        CC_act = max(CC_act, 0.0)
        ATR_act = ATR_act + dt * dATR_act_dt
        ATR_act = min(ATR_act, ATR_tot)
!        if (single_cell .and. it == Nt) then
!            write(*,'(a,4e12.3)')  'terms:  ',(Kkcc2a + CC_act) * CC_inact / (Kmccp + CC_inact), &
!                     - cp%Kt2cc * ATM_act * CC_act / (Kmccmd + CC_act), &
!                     - cp%Ke2cc * ATR_act * CC_act / (Kmccrd + CC_act), &
!                    dCC_act_dt!
!
!            write(*,'(3i6,e12.3)') istep,Nt,it,dCC_act_dt
!        endif
    elseif (use_ATR .and. D_ATR > 0) then
        dATR_act_dt = Kd2e * D_ATR * ATR_inact / (Kmrp + ATR_inact)  - Kcc2e * ATR_act * CC_act / (Kmrd + CC_act)
        ATR_act = ATR_act + dt * dATR_act_dt
        ATR_act = min(ATR_act, ATR_tot)
!        if (single_cell) write(nflog,'(a,i4,4f8.3)') 'iph,D_ATR, Kd2e, ATR_inact, dATR_act_dt: ', &
!                    iph, D_ATR, Kd2e, ATR_inact, dATR_act_dt
    endif
!    if (iph == G2_phase) then   ! just testing to see why G2 is so extended with kiliakis = 1
        dATM_act_dt = Kd2t * D_ATM * ATM_inact / (Kmmp + ATM_inact) - Kti2t * ATM_act / (Kmmd + ATM_act)   
!    else
!        dATM_act_dt = 0
!    endif
    ATM_act = ATM_act + dt*dATM_act_dt
    ATM_act = min(ATM_act, ATM_tot)
!    ATM_act = min(ATM_act,0.5)
    !if (single_cell .and. iph == G2_phase .and. t_simulation < 3.2*3600) then
    !    write(*,'(a,4f8.4,e12.3)') 'Jaiswal: dt,D_ATR,ATR_inact,ATR_act,dATRdt: ',dt,D_ATR,ATR_inact,ATR_act,dATR_act_dt
    !endif
    t = it*dt
!    if (t_G2 <= 0.1 .and. it <= 10) write(nflog,'(a,3f8.4)') 'D_ATM,ATM_act,ATM_inact: ',D_ATM,ATM_act,ATM_inact
enddo
!if (kcell_now == 8) write(nflog,'(a,5f9.4,e12.3)') 'Jaiswal: ',D_ATM,Kd2t,Kmmp,Kti2t,Kmmd,dATM_act_dt
!if (single_cell) then
!    write(*,'(a,2i4,4f8.4)') 'kcell,iph,ATR,ATM: ',kcell_now,iph,ATR_act,ATM_act
!    write(nflog,'(a,2i4,4f8.4)') 'kcell,iph,ATR,ATM: ',kcell_now,iph,ATR_act,ATM_act
!endif
cp%ATM_act = ATM_act
if (iph == G2_phase) then
    cp%CC_act = CC_act
    cp%ATR_act = ATR_act
    cp%dCC_act_dt = dCC_act_dt
!    write(nflog,'(a,i8,2e12.3)') 'Jaiswal_update: ',kcell_now,CC_act0,CC_act
!    write(*,'(a,6f8.3)') 't_simulation, ATM_act, ATR_act, CC_act, D_ATM, D_ATR: ',t_simulation/3600,ATM_act, ATR_act, CC_act,D_ATM,D_ATR
elseif (iph == S_phase .and. use_ATR) then
    cp%ATR_act = ATR_act
endif
t = t_simulation/3600.
tIR = (tnow - t_irradiation)/3600
!if (test_run .and. (kcell_now==3 .or. kcell_now==6)) then
!    write(nflog,'(a,2i4,f8.2,2e12.3,f8.3)') 'Jaiswal_update: kcell,iph,tIR,ATM_act,ATR_act,CC_act: ',kcell_now,iph,tIR,ATM_act,ATR_act,CC_act
!endif
!if (kcell_now == 18) write(nflog,'(a,2i4,f8.2,2e12.3,f8.3)') 'Jaiswal_update: kcell,iph,tIR,ATM_act,ATR_act,CC_act: ',kcell_now,iph,tIR,ATM_act,ATR_act,CC_act
!write(nflog,'(a,9f8.3)') 'tIR,CC_act,dCC_act_dt,d(1:3), ATM_act, ATR_act: ',tIR,CC_act,dCC_act_dt, d(1:3), ATM_act, ATR_act
!write(nflog,'(a,4f8.3,e12.3)') 'tIR, D_ATM, ATR_act, ATM_act,dCC_act_dt: ',tIR,D_ATM, ATR_act, ATM_act,dCC_act_dt
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
!    Nt = 10
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

if (eta > 1.0E-3 .and. test_run) write(nflog,'(a,i6,2e12.3)') 'kcell,eta,Pmis: ',kcell_now,eta,Pmis

!if (kcell_now == 1) write(nflog,'(a,2f8.1,2e12.3)') 'initialBreaks, finalBreaks, eta, Pmis: ',initialBreaks, finalBreaks,eta,Pmis

!write(*,*) 'misrepairRate: '
!write(*,'(a,3f6.1)') 'initialBreaks, finalBreaks, repairedBreaks: ',initialBreaks, finalBreaks, repairedBreaks
!write(*,'(a,2e12.3)') 'eta, Pmis: ',eta,Pmis
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
!if (single_cell .and. phase > G2_phase) then
!    write(nflog,*) 'updateRepair: phase: ',phase
!    stop
!endif
if (cp%DSB0(TMEJ,1) /= 0) then
    write(*,*) 'DSB0(TMEJ,1): ',cp%DSB0(TMEJ,1)
    write(nflog,*) 'a DSB0(TMEJ,1): ',cp%DSB0(TMEJ,1)
    stop
endif
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
!    if (iph == G2_phase) write(*,'(a,2i6,e12.3)') 'did Jaiswal: ',kcell_now, cp%generation, cp%ATM_act
endif

!if ((iph == 1 .and. use_G1_pATM) .or. (iph == 2 .and. use_S_pATM)) then 
!    call updateATM(iph,cp%pATM,ATM_DSB,dth)     ! updates the pATM mass through time step = dth
!endif

DSB = 0
do k = 1,3
do jpp = 1,2
    call pathwayRepair(k, dth, DSB0(k,jpp), DSB(k,jpp))
    if (DSB(k,jpp) < DSB_min) DSB(k,jpp) = 0
!    if (dbug .and.k == 1 .and. jpp == 1) write(nfres,'(a,3f8.3)') 'DSB0,DSB,repratefactor: ',DSB0(k,jpp),DSB(k,jpp),repratefactor(1)
enddo
enddo
! DSB0(k) is the count before repair, DSB(k) is the count after repair

#if 0
! The following commented out code follows MEDRAS
totDSB = sum(DSB)
totDSBinfid0 = 0
totDSBinfid = 0
do k = 1,NP
    if (DSB0(k) > 0 .and. fidRate(k) < 1.0) then
        totDSBinfid0 = totDSBinfid0 + DSB0(k)
        totDSBinfid = totDSBinfid + DSB(k)
    endif
enddo

! Testing eta dependence
!eta_G1 = eta_G2 ! -> Pmis = 0.01 - 0.07
!eta_G2 = eta_G1 ! -> Pmis = 0.1 - 0.127
! Therefore small eta -> small Pmis -> smaller Nmis (Pmis is ~0.05 for eta_G2, ~0.11 for eta_G1) -> bigger SF

if (phase == G1_phase) then
    eta = eta_G1
elseif (phase == S_phase) then
    eta = eta_G1 + cp%progress*(eta_G2 - eta_G1)
elseif (phase == G2_phase) then
    eta = eta_G2
endif
if (use_DSBinfid) then
    Pmis = misrepairRate(totDSBinfid0, totDSBinfid, eta)
else
    Pmis = misrepairRate(totDSB0, totDSB, eta)
endif
cp%totMis = cp%totMis + Pmis*(totDSB0 - totDSB)
binMisProb = cp%totMis/(cp%totDSB0 - totDSB)
if (single_cell) write(nflog,'(a,2e12.3)') 'Pmis, binMisProb: ',Pmis,binMisProb
Nmis = 0
dNmis = 0
do k = 1,NP
!    if (pathUsed(k) .and. fidRate(k) < 1.0) then
    if (fidRate(k) < 1.0) then
        if (use_totMis) then
            dNmis(k) = (DSB0(k) - DSB(k))*(1 - fidRate(k)*(1-binMisProb))
            Nmis = Nmis + dNmis(k)
        else
            dNmis(k) = (DSB0(k) - DSB(k))*(1 - fidRate(k)*(1-Pmis))
            ! Approx = (number of repairs)*(1 - fidRate) = NDSB*repRate*(1 - fidRate)
!            if (k == 2 .or. k == 4) write(nfphase,'(a,i4,4f8.4)') 'k: ',k,(DSB0(k) - DSB(k)),fidRate(k),Pmis,fidRate(k)*(1-Pmis)
            Nmis = Nmis + dNmis(k)
!            if (DSB0(k) > DSB(k)) then
!                write(nflog,'(a,i3,3f8.3,e12.3)') 'dNmis: ',k,(DSB0(k) - DSB(k)),Pmis,(1 - fidRate(k)*(1-Pmis)),dNmis
!            endif
        endif
    endif
enddo
write(nflog,'(a,6f8.3)') 'dNmis, Nmis: ',dNmis,Nmis
if (isnan(Nmis)) then
    write(*,*) 'updateRepair: Nmis isnan'
    stop
endif
#endif

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
    !    eta_NHEJ = etafun(1.d0,sigma_NHEJ)
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

! Here we could increment 
! Nreptot by Nrep = totDSB0 - totDSB
! Pmistot by Nrep*Pmis

!if (single_cell .and. (tnow < CA_time_h*3600)) then
!    write(nflog,'(a,i4,f6.2,3f6.1,e12.3,3f8.4)') &
!                 'iph, pr, totDSB0, totDSB, DSB_rep, eta, Pmis, dmis, tIR: ', &
!                     phase,cp%progress,totDSB0, totDSB, totDSB0-totDSB, eta_NHEJ, Pmis, dmis, tIR
!endif
Nmis(1) = Nmis(1) + dmis*(1 - f_S)  ! doubled at mitosis
Nmis(2) = Nmis(2) + dmis*f_S
misjoins(1) = misjoins(1) + Nmis(1) + Nmis(2)


!if (single_cell) &
!write(nfout,'(a,5f8.4)') 'f_S, tIR, totDSB0, eta_NHEJ, Nmis: ',f_S, tIR, totDSB0, eta_NHEJ, Nmis
!write(nflog,'(a,4f10.4)') 'totDSB0, totDSB, eta_NHEJ, Pmis: ',totDSB0, totDSB, eta_NHEJ, Pmis
!if (tIR == 2.0) then
!    write(nflog,*) 'tIR: ',tIR
!    write(nflog,*) 'Varying eta to find how Pmis depends on eta'
!    do k = 1,11
!        eta = eta_NHEJ*(0.5 + (k-1)*0.1)
!        Pmis = misrepairRate(totDSB0, totDSB, eta)
!        write(nflog,'(a,i4,2f8.4)') 'k, eta, Pmis: ',k,eta,Pmis
!    enddo
!endif
if (sum(DSB0(TMEJ,:)) > 0) then ! not used currently
    ! For TMEJ pathway
    Pmis = misrepairRate(sum(DSB0(TMEJ,:)), sum(DSB(TMEJ,:)), eta_TMEJ)
    dmis = Pmis*(sum(DSB0(TMEJ,:)) - sum(DSB(TMEJ,:)))
    Nmis(1) = Nmis(1) + dmis*(1 - f_S)
    Nmis(2) = Nmis(2) + dmis*f_S
    misjoins(2) = misjoins(2) + Nmis(1) + Nmis(2)
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
!    write(nfres,'(a,i6,4f8.3,e12.3)') 'kcell,totDSB,Nmis,Psurvive: ', kcell_now,totDSB,Nmis,cp%Psurvive
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
        if (kcell_now == 1) write(*,'(a,2f8.3,6e12.3)') 'totDSB,Pmit,Nmis,Paber: ',&
                            totDSB,Pmit,Nmis,Paber  
        if (single_cell) write(nfres,'(a,2f8.3,6e12.3)') '(2) totDSB,Pmit,Nmis,Paber: ',totDSB,Pmit,Nmis,Paber  
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

#if 0
Nlethal_sum = Klethal*(2*Nmis(1) + Nmis(2))
if (Nlethal_sum > NMDIST*ddist_Nlethal) then
    count_Nlethal(NMDIST) = count_Nlethal(NMDIST) + 1
else
    k = Nlethal_sum/ddist_Nlethal + 1
    count_Nlethal(k) = count_Nlethal(k) + 1
endif 
if (sum(totDSB) > NMDIST*ddist_totDSB) then
    count_totDSB(NMDIST) = count_totDSB(NMDIST) + 1
else
    k = sum(totDSB)/ddist_totDSB + 1
    count_totDSB(k) = count_totDSB(k) + 1
endif
#endif

if (single_cell) then
    Nmistot = 2*Nmis(1) + Nmis(2)
    write(*,'(a,i1,6f6.1,8f8.3,3f8.4)') 'AAA ',synch_phase,synch_fraction,cp%DSB0(1:2,1),cp%DSB0(1:3,2),t_mitosis, &
            totDSB,Pmit,Nmis,Nmistot,Paber,cp%psurvive  !,cp%psurvive_nodouble
    write(*,'(a,f8.3)') 'mitosis_time: ',cp%mitosis_time/3600
! write out signalling results
!open(nfpar, file='signal.out', status = 'replace')
!do k = 1,istep_signal
!    write(nfpar,'(4f10.6)') signalling(:,k)
!enddo
!close(nfpar)
endif
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

#if 0
!------------------------------------------------------------------------
! Get average DNA growth factor for cells in S-phase
!------------------------------------------------------------------------
subroutine get_DNA_synthesis_rate(DNA_rate)
real(8) :: DNA_rate
type(cell_type), pointer :: cp
integer :: kcell, iph, cnt
real(8) :: atm, k3m, k4m, fATM, atr, k3r, k4r, fATR, rate, rate_sum, atm_ave

!write(*,'(a)') 'get_DNA_synthesis_rate'
k3m = K_ATM(S_phase,3)
k4m = K_ATM(S_phase,4)
k3r = K_ATR(S_phase,3)
k4r = K_ATR(S_phase,4)
fATR = 1.0
cnt = 0
rate_sum = 0
atm_ave = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    iph = cp%phase
    if (iph == S_phase) then
        cnt = cnt + 1
        if (use_S_pATM) then
            atm = cp%pATM
        else
            atm = cp%ATM_act
        endif
        atm_ave = atm_ave + atm
!        fATM = max(0.01,1 - k3*atm/(k4 + atm))  ! 0.01 ??
        if (use_exp_slowdown) then
            fATM = exp(-k4m*atm)
        else
            fATM = max(0.01,1 - k3m*atm/(k4m + atm))
        endif
        if (use_ATR_S) then
            atr = cp%ATR_act
            fATR = max(0.01,1 - k3r*atr/(k4r + atr))
        endif
        rate = fATM*fATR
        rate_sum = rate_sum + rate
    endif
enddo
DNA_rate = rate_sum/cnt
atm_ave = atm_ave/cnt
!write(*,'(a,3f8.3,e12.3)') 'DNA growth rate factor: ',DNA_rate,k3,k4,atm_ave
end subroutine
#endif

!------------------------------------------------------------------------
!------------------------------------------------------------------------
subroutine test_Pmis
type(cell_type), pointer :: cp
real(8) :: dose = 6, dt = 0.1, Reffmin = 0.9, sigma = 0.0413, dsigma_dt = 0.0239, Kcoh = 1.0
real(8) :: r1 = 2.081, r2 = 0.2604, p = 0.43, T_S = 9.04
real(8) :: t, N(2), dN(2), Ntot0, Ntot, R, S, eta, eta_A
real(8) :: Pmis(100), Nrep(100), f_S, fsigma, Ptot, Nreptot, Pmis_ave
integer :: it
integer :: nhours = 2

write(*,*) 'test_Pmis'
! Evaluate Pmis at intervals over a period after IR, compute weighted average 
! G1, IR at 0.0
! Initial DSBs
kcell_now = 1
cp => cell_list(1)
cp%phase = G2_phase
cp%progress = 0.3
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
write(nflog,*)
write(nflog,'(a,i4,f6.2)') 'phase, progress: ',cp%phase,cp%progress
write(nflog,'(a,5f8.3)') 'dose, DSB0: ',dose,cp%DSB0(1:2,:)
write(nflog,'(a,2f8.4)') 'Reffmin, dsigma_dt: ',Reffmin, dsigma_dt
! Repair
N(1) = cp%DSB0(1,1) + cp%DSB0(1,2)  ! NHEJfast
N(2) = cp%DSB0(2,1) + cp%DSB0(2,2)  ! NHEJslow

do it = 1,nhours*10
    t = it*dt
    if (cp%phase == 1) then
        f_S = 0
        R = (1 - Reffmin)*exp(-Kclus*t) + Reffmin
    elseif (cp%phase == 2) then
        f_S = t/T_S + cp%progress
        R = (1 - f_S)*Reffmin + f_S*1.26
    elseif (cp%phase == 3) then
        f_S = 1
        R = 1.26
    endif
    fsigma = 1 - (1 - Kcoh)*f_S
    S = fsigma*(sigma + t*dsigma_dt)
    eta = etafun(R,S)
!    eta_A = eta_Arnould(f_S, t, Rmin, sigma, Kcoh)
!    write(nflog,*) 'eta_Arnould: ',eta_A   ! checking that eta_A = eta   OK!
!    write(nflog,'(a,4f8.4,e12.3)') 'R, S, sigma,dsigma_dt,eta: ',R,S, sigma,dsigma_dt,eta
    Ntot0 = N(1) + N(2)
    dN(1) = r1*dt*N(1)
    dN(2) = r2*dt*N(2)
    N(1) = N(1) - dN(1)
    N(2) = N(2) - dN(2)
    Ntot = N(1) + N(2)
    Pmis(it) = misrepairRate(Ntot0, Ntot, eta)
    Nrep(it) = dN(1) + dN(2)
    write(nflog,'(i4,5f8.2,e12.3,f8.4)') it, f_S, R, Ntot0, Ntot, Nrep(it),eta, Pmis(it)
enddo   
Ptot = 0
Nreptot = 0
do it = 1,nhours*10
    Ptot = Ptot + Nrep(it)*Pmis(it)
    Nreptot = Nreptot + Nrep(it)
enddo
Pmis_ave = Ptot/Nreptot
write(nflog,*) 'Pmis_ave: ',Pmis_ave
return
end subroutine

end module