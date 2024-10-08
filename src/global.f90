! Global definitions

module global

use real_kind_mod
use par_zig_mod
use, intrinsic :: ISO_C_BINDING

implicit none

integer, parameter :: NP = 4

integer, parameter :: NORMAL_DIST      = 1
integer, parameter :: LOGNORMAL_DIST   = 2
integer, parameter :: EXPONENTIAL_DIST = 3
integer, parameter :: CONSTANT_DIST    = 4

integer, parameter :: ALIVE = 1
integer, parameter :: DYING = 2
integer, parameter :: DEAD = 3
integer, parameter :: DIVIDED = 4

integer, parameter :: G1_phase      = 1
integer, parameter :: S_phase       = 2
integer, parameter :: G2_phase      = 3
integer, parameter :: M_phase       = 4
integer, parameter :: dividing      = 5
! checkpoint phases are no longer used, instead phase progress is slowed
integer, parameter :: G1_checkpoint = 6
integer, parameter :: S_checkpoint  = 7
integer, parameter :: G2_checkpoint = 8

integer, parameter :: nfin=10, nfout=11, nflog=12, nfres=13, nfrun=14, nfcell=15, nftreatment=16, nfphase=17, &
					  nfpar=18, nftcp=20, nfpest = 21

integer, parameter :: MAX_CELLTYPES = 2
integer, parameter :: max_nlist = 300000
real(REAL_KIND), parameter :: PI = 4.0*atan(1.0)

type cell_type
	integer :: ID
	integer :: celltype
	logical :: active
	integer :: state
	logical :: irradiated
	integer :: generation
	real(REAL_KIND) :: divide_time      ! cycle time
	real(REAL_KIND) :: fg(4)			! to make sum(T_G1, T_S, T_G2, T_M) consistent with Tdivide
	real(REAL_KIND) :: t_divide_last	! these two values are used for colony simulation
	real(REAL_KIND) :: t_divide_next
	real(REAL_KIND) :: t_S_phase
	real(REAL_KIND) :: birthtime
    real(REAL_KIND) :: t_mitosis		! time from IR to end of mitosis
	real(REAL_KIND) :: t_start_mitosis
	real(REAL_KIND) :: t_start_G2
	real(REAL_KIND) :: G2_time
	real(REAL_KIND) :: mitosis

	! Cell cycle 
    integer :: phase
    real(REAL_KIND) :: progress, fp, mitosis_duration
    
	! DRM section
	integer :: phase0
	real(8) :: pATM, pATR, DSB(NP,2), DSB0(NP,2),totDSB0
	real(8) :: Psurvive, mitosis_time
	real(8) :: Nmis(2)
	
	! Jaiswal section (26/09/22)
	real(REAL_KIND) :: CC_act, ATR_act, ATM_act, dCC_act_dt, kt2cc, ke2cc, kcc2a
    
end type

type cycle_parameters_type
    real(REAL_KIND) :: f_G1, f_S, f_G2, f_M
    real(REAL_KIND) :: T_G1, T_S, T_G2, T_M
end type

type dist_type
	integer :: class
	real(REAL_KIND) :: p1, p2, p3
end type

type(dist_type) :: divide_dist(MAX_CELLTYPES)
type(cell_type), allocatable, target :: cell_list(:)

integer :: initial_count

integer :: nlist, Ncells, Ncells0, ncells_mphase, lastID, Ncelltypes
integer :: Ncells_type(MAX_CELLTYPES), Ndying(MAX_CELLTYPES), Nviable(MAX_CELLTYPES), Ndead(MAX_CELLTYPES)

type(cycle_parameters_type), target :: cc_parameters(MAX_CELLTYPES)

integer :: istep, ndays, nsteps 
integer :: Mnodes
real(REAL_KIND) :: DELTA_T, tnow   
real(REAL_KIND) :: divide_time_median(MAX_CELLTYPES), divide_time_shape(MAX_CELLTYPES), divide_time_mean(MAX_CELLTYPES)
real(REAL_KIND) :: t_simulation
real(REAL_KIND) :: start_wtime
real(REAL_KIND) :: IR_time_h, CA_time_h, washout_time_h

integer, allocatable :: gaplist(:)
integer :: ngaps, ndivided
integer, parameter :: max_ngaps = 200000

character*(2048) :: inputfile
character*(2048) :: outputfile

logical :: simulation_start, par_zig_init
logical :: is_radiation
logical :: use_gaplist = .true.
logical :: dbug = .false.

integer :: seed(2)
integer :: kcell_now

! PEST variables
logical :: use_PEST = .false.
character*(128) :: PEST_outputfile

logical :: use_synchronise   ! now set in main
integer :: synch_phase
real(REAL_KIND) :: synch_fraction
logical :: single_cell

! DRM section
logical :: DRM = .true.
logical, parameter :: use_Napop = .true.   ! use count of apoptosed cells in SFave calculation - true for consistency with CA
integer :: NPsurvive, Nirradiated, Napop, Nmitotic
real(REAL_KIND), allocatable :: Psurvive(:)
logical :: include_daughters = .true.
real(REAL_KIND) :: radiation_dose, t_irradiation, SFave, t_mitosis
logical :: use_SF = .true.
real(REAL_KIND) :: phase_hour(60)
integer :: nphase_hours, next_phase_hour
real(REAL_KIND) :: phase_dist(0:4)    ! % of cells in each phase
real(REAL_KIND) :: recorded_phase_dist(60,0:4)   ! % of cells in each phase phase_hour after IR
real(REAL_KIND) :: recorded_DNA_rate(60)         ! average S-phase DNA rate in each phase_hour after IR
integer, allocatable :: nphase(:,:)
real(REAL_KIND) :: totNmis = 0
integer :: maxhours = 199
logical :: overstepped

logical, parameter :: constant_S_pHR = .true.
real(REAL_KIND), parameter :: dose_threshold = 1
integer :: ATR_in_S = 1		! 0 = no ATR signalling in S, 1 = signalling, no CP effect, 2 = signalling and CP effect
logical, parameter :: use_Arnould = .true.
real(REAL_KIND) :: R_Arnould = 0.7, Kclus = 0.693	! for DSB clustering
logical :: use_cell_kcc2a_dependence = .true.

logical :: compute_cycle

real(REAL_KIND) :: mitosis_std = 0.1336*3600     ! Chao 2019, Erlang k=14, L = 28, hours -> seconds

! Drug half-life simulation
logical :: use_drug_halflife
real(REAL_KIND) :: Khalflife, drug_time, drug_conc

logical :: test_run = .false.	! to check Psurvive etc
LOGICAL :: use_no_random = .false.	! to turn off variation in cycle time, DSB_Gy


! DEBUGGING
integer :: tracked1(500)
integer :: ntrack1 = 0
integer :: tracked2(500)
integer :: ntrack2 = 0

!DEC$ ATTRIBUTES DLLEXPORT :: nsteps, use_PEST, PEST_outputfile
!DEC$ ATTRIBUTES DLLEXPORT :: use_synchronise, synch_phase, synch_fraction
contains

!-----------------------------------------------------------------------------------------
! WTIME returns a reading of the wall clock time.
!-----------------------------------------------------------------------------------------
real(DP) function wtime()
!DEC$ ATTRIBUTES DLLEXPORT :: wtime
  integer :: clock_max, clock_rate, clock_reading

  call system_clock ( clock_reading, clock_rate, clock_max )
  wtime = real(clock_reading,kind=DP)/clock_rate
end function

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
integer :: n1,n2,kpar
integer :: k,R

if (n1 == n2) then
    random_int = n1
elseif (n1 > n2) then
    write(nflog,*) 'ERROR: random_int: n1 > n2: ',n1,n2
    stop
endif
R = par_shr3(kpar)
k = abs(R)
random_int = n1 + mod(k,(n2-n1+1))
end function

!--------------------------------------------------------------------------------
! Make a random choice of an integer from 1 - N on the basis of probabilities in
! the array p(:) (assumed to be normalized).
!--------------------------------------------------------------------------------
integer function random_choice(p,N,kpar)
integer :: N,kpar
real(REAL_KIND) :: p(:)
integer :: k
real(REAL_KIND) :: R, psum

R = par_uni(kpar)
psum = 0
do k = 1,N
    psum = psum + p(k)
    if (R <= psum) then
        random_choice = k
        return
    endif
enddo
write(nflog,*) 'ERROR: random_choice: ',N,p
stop
end function

!--------------------------------------------------------------------------------
! Returns a permutation of the elements of a()
!--------------------------------------------------------------------------------
subroutine permute(a,n,kpar)
integer :: a(*),n,kpar
integer :: i,k,tmp

do i = 1,n
    k = random_int(1,n,kpar)
	tmp = a(i)
	a(i) = a(k)
	a(k) = tmp
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function norm(r)
real(REAL_KIND) :: r(3)

norm = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function norm2(r)
real(REAL_KIND) :: r(3)

norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine normalize(r)
real(REAL_KIND) :: r(3)

r = r/norm(r)
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine get_vnorm(v,vnorm)
real(REAL_KIND) :: v(3), vnorm(3)
real(REAL_KIND) :: d

d = dot_product(v,v)
vnorm = v/sqrt(d)
end subroutine

!-----------------------------------------------------------------------------------------
! Assumes unconstrained growth rate - no checkpoint delays
! M-phase duration has already been assigned as cp%mitosis_duration
! Tinter is interphase time.
!-----------------------------------------------------------------------------------------
subroutine set_phase_times(cp)
type(cell_type), pointer :: cp
real(REAL_KIND) :: Tdiv, fg(4), Tinter, Tinter_ave
real(REAL_KIND) :: T_G1, T_S, T_G2, T_M
integer :: ityp
type(cycle_parameters_type), pointer :: ccp

ityp = cp%celltype
ccp => cc_parameters(ityp)

T_M = cp%mitosis_duration
Tdiv = DivideTime(ityp)     ! log-normally distributed r.v.
Tdiv = min(Tdiv,1.2*divide_time_median(ityp))	! limit max divide time
Tinter = Tdiv - T_M
Tinter_ave = ccp%T_G1 + ccp%T_S + ccp%T_G2
T_G1 = Tinter*ccp%T_G1/Tinter_ave	! G1, S and G2 times are scaled mean phase times
T_S = Tinter*ccp%T_S/Tinter_ave
T_G2 = Tinter*ccp%T_G2/Tinter_ave
! Note: fg(phase) is actually independent of phase for G1, S, G2 = Tinter/Tinter_ave
fg(G1_phase) = T_G1/ccp%T_G1
fg(S_phase) = T_S/ccp%T_S
fg(G2_phase) = T_G2/ccp%T_G2
fg(M_phase) = 1.0
cp%divide_time = Tdiv   ! cycle time, varies with cell
cp%fg = fg
end subroutine	

!--------------------------------------------------------------------------------------
! This is actually cycle time
!--------------------------------------------------------------------------------------
real(REAL_KIND) function DivideTime(ityp)
integer :: ityp
real(REAL_KIND) :: p1, p2
integer :: kpar = 0

dividetime = 0
p1 = divide_dist(ityp)%p1
p2 = divide_dist(ityp)%p2
select case (divide_dist(ityp)%class)
case (NORMAL_DIST)
	DivideTime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
	DivideTime = rv_lognormal(p1,p2,kpar)
case (CONSTANT_DIST)
	DivideTime = p1
end select
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function DivisionTime(ityp)
integer :: ityp
integer :: kpar = 0
real(REAL_KIND), parameter :: rndfraction = 0.2

DivisionTime = rv_lognormal(divide_dist(ityp)%p1,divide_dist(ityp)%p2,kpar)
end function

!--------------------------------------------------------------------------------------
! par_rnor is N(0,1)
! p1 = mean
! p2 = std deviation
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_normal(p1,p2,kpar)
integer :: kpar
real(REAL_KIND) :: p1,p2
real(REAL_KIND) :: R

R = par_rnor(kpar)
rv_normal = p1+R*p2
end function

!--------------------------------------------------------------------------------------
! When Y is normal N(p1,p2) then X = exp(Y) is lognormal with
!   median = m = exp(p1)
!   shape  = s = exp(p2)
! Also median = m = mean/(s^2/2)
! kpar = parallel process number
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_lognormal(p1,p2,kpar)
integer :: kpar
real(REAL_KIND) :: p1,p2
real(REAL_KIND) :: R,z

R = par_rnor(kpar)
z = p1 + R*p2
rv_lognormal = exp(z)
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_exponential(lambda,kpar)
real(REAL_KIND) :: lambda
integer :: kpar
real(REAL_KIND) :: R

R = par_uni(kpar)
rv_exponential = -log(1-R)/lambda
end function

!--------------------------------------------------------------------------------------
! Cumulative probability distribution for a lognormal variate with median m, shape s
! Computes Pr(X < a) where X = exp(Y) and Y is N(p1,p2), p1 = log(m), p2 = log(s)
! Pr(X < a) = Pr(Y < log(a)) = Pr(p1 + p2*R < log(a)) = Pr(R < (log(a)-p1)/p2)
! where R is N(0,1)
!--------------------------------------------------------------------------------------
real(REAL_KIND) function cum_prob_lognormal(a,p1,p2)
real(REAL_KIND) :: a, p1, p2
real(REAL_KIND) :: b, prob

b = (log(a) - p1)/p2
prob = 0.5 + 0.5*erf(b/sqrt(2.0))
cum_prob_lognormal = prob
end function

!-----------------------------------------------------------------------------------------
! Seconds
!-----------------------------------------------------------------------------------------
function get_mitosis_duration() result(t)
real(REAL_KIND) :: t
type(cycle_parameters_type), pointer :: ccp
integer :: ityp = 1
ccp => cc_parameters(ityp)
if (single_cell .or. test_run .OR. use_no_random) then
    t = ccp%T_M
else
    t = rv_normal(ccp%T_M, mitosis_std, 0)
endif
t = max(t,0.0)
end function



!--------------------------------------------------------------------------------------
! Determine real roots r(:) of the cubic equation:
! x^3 + a.x^2 + b.x + c = 0
! If there is one real root, n=1 and the root is r(1)
! If there are three distinct real roots, n=3 and the roots are r(1), r(2), r(3)
! If there is a repeated root, n=2 and the single root is r(1), the repeated root is r(2)
!--------------------------------------------------------------------------------------
subroutine cubic_roots(a, b, c, r, n)
real(REAL_KIND) :: a, b, c, r(3)
integer :: n
real(REAL_KIND) :: QQ, RR, theta, R2, Q3, AA, BB

QQ = (a*a - 3*b)/9
RR = (2*a*a*a - 9*a*b + 27*c)/54
Q3 = QQ*QQ*QQ
R2 = RR*RR
if (R2 < Q3) then
	n = 3
	theta = acos(RR/sqrt(Q3))
	r(1) = -2*sqrt(QQ)*cos(theta/3) - a/3
	r(2) = -2*sqrt(QQ)*cos((theta+2*PI)/3) - a/3
	r(3) = -2*sqrt(QQ)*cos((theta-2*PI)/3) - a/3
else
	n = 1
	AA = -sign(1.d0,RR)*(abs(RR) + sqrt(R2 - Q3))**(1.d0/3.d0)
	if (AA == 0) then
		BB = 0
	else
		BB = QQ/AA
	endif
	r(1) = AA + BB - a/3
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Squeeze gaps out of cellist array, adjusting occupancy array.
!-----------------------------------------------------------------------------------------
subroutine squeezer()
integer :: last, kcell, site(3), indx(2), i, j, idc, n, region

if (ngaps == 0) return
last = nlist
kcell = 0
n = 0
do
    kcell = kcell+1
    if (cell_list(kcell)%state == DEAD) then    ! a gap
        if (kcell == last) exit
        do
            if (last == 0) then
                write(nflog,*) 'last = 0: kcell: ',kcell
                stop
            endif
            if (cell_list(last)%state == DEAD) then
                last = last-1
                n = n+1
                if (n == ngaps) exit
            else
                exit
            endif
        enddo
        if (n == ngaps) exit
        cell_list(kcell) = cell_list(last)
        last = last-1
        n = n+1
    endif
    if (n == ngaps) exit
enddo
nlist = nlist - ngaps
ngaps = 0

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_phase_distribution(phase_count)
integer :: phase_count(0:4)
integer :: kcell, ph
type(cell_type), pointer :: cp

phase_count = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) then
        ph = 0      ! 
    elseif (cp%phase == G1_phase .or. cp%phase == G1_checkpoint) then
        ph = 1
    elseif (cp%phase == S_phase .or. cp%phase == S_checkpoint) then
        ph = 2
    elseif (cp%phase == G2_phase .or. cp%phase == G2_checkpoint) then
        ph = 3
    else
        ph = 4
    endif
    phase_count(ph) = phase_count(ph) + 1
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! https://hpaulkeeler.com/simulating-poisson-random-variables-direct-method/
!-----------------------------------------------------------------------------------------
function poisson_gen(L) result(res)
real(8) :: L
integer :: k, res, kpar=0
real(8) :: p, u
k = 0
p = 1
do
	k = k+1
	u = par_uni(kpar)
	p = p*u
	if (p < L) exit
enddo
res = k-1
end function


end module
