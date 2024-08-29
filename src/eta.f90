module eta_module

use global

implicit none

integer, parameter :: Neta = 11
integer, parameter :: NtIR = 15
real(8) :: dsigma_dt
real(8) :: eta_table(Neta,NtIR,NP)

contains

!--------------------------------------------------------------------------
! From MEDRAS
! S is sigma
! R is Radius
! d is r = 2*Radius
!--------------------------------------------------------------------------
function thetaM(R, S) result(thetaVal)
real(8) :: R, S
real(8) :: d, d2, R2, R3, S2
real(8) :: termOne, termTwo, termThree, scaling, thetaVal

d = 2*R
d2 = d*d
R2 = R*R
R3 = R*R2
S2 = S*S
termOne = 8*sqrt(2*pi)*R3*S*erf(d/(S*sqrt(2.)))
termTwo = -exp(-d2/(2*S2))*(d2*d2 + 4*d2*(S2 - 3*R2) + 16*d*R3 + 8*S2*(S2 - 3*R2))
termThree = 8*S2*S2 - 24*R2*S2
scaling = pi*S2/(4*R3)
thetaVal = (termOne + termTwo + termThree)*scaling
end function

!--------------------------------------------------------------------------
! From Eq 15 in McMahon2021
! S is sigma
!--------------------------------------------------------------------------
function theta15(R, S) result(thetaVal)
real(8) :: R, S
real(8) :: R2, S2
real(8) :: termOne, termTwo, termThree, scaling, thetaVal

R2 = R*R
S2 = S*S
termOne   = sqrt(2*pi)*R2*S*erf(R*sqrt(2.0)/S)
termTwo   = -exp(-4*R2/(2*S2))*(S2 -R2*S2)
termThree = S2 - 3*R2*S2
scaling   = 2*pi*S2/(R*R2)
thetaVal = (termOne + termTwo + termThree)*scaling
end function	

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine check_eta(S)
real(8) :: S
real(8) :: R, eta

write(*,*) 'check_eta:'
R = 1
eta = (6/(4*pi*R**3))*theta15(R,S)
write(*,*) 'eta from theta15: ',eta
eta = (6/(4*pi*R**3))*thetaM(R,S)
write(*,*) 'eta from thetaM: ',eta
end subroutine

!--------------------------------------------------------------------------
! From Eq 14 in McMahon2021
! S is sigma
!--------------------------------------------------------------------------
function etafun(R,S) result(eta)
real(8) :: R, S, eta

!eta = (6/(4*pi*R**3))*theta15(R,S)
eta = (6/(4*pi*R**3))*thetaM(R,S)
end function	

!--------------------------------------------------------------------------
! fsmin is Kcoh
! Note that the new expression for Reff covers the whole cycle - phase not used.
!--------------------------------------------------------------------------
function eta_Arnould(phase,f_S, tIR, Reffmin, S_NHEJ, fsmin) result(eta)
integer :: phase
real(8) :: f_S, tIR, Reffmin, S_NHEJ, fsmin, eta
real(8) :: Reff, sigma, fsigma
logical, parameter :: use_old_method = .false.   ! results the same as new method

Reff = (1 - f_S)*((1 - Reffmin)*exp(-Kclus*tIR) + Reffmin) + f_S*1.26
fsigma = 1 - (1 - fsmin)*f_S
sigma = S_NHEJ + tIR*dsigma_dt
sigma = fsigma*sigma

eta = etafun(Reff,sigma)
end function

!--------------------------------------------------------------------------
! S_NHEJ, S_TMEJ are sigma_NHEJ, sigma_TMEJ
! fsmin is the minimum of the multiplying factor fsigma that reduces sigma 
! as f_S increases.  This is the cohesin effect: fsmin = f_coh
!--------------------------------------------------------------------------
subroutine make_eta_table(S_NHEJ, S_TMEJ, fsmin)
real(8) :: S_NHEJ, S_TMEJ, fsmin
real(8) :: V, R, S, eta, fsigma
integer :: k, it

write(nflog,*) 'make_eta_table: ',S_NHEJ, S_TMEJ
write(nflog,*) 'Neta, NtIR: ',Neta,NtIR
fsigma = 1
do k = 1,Neta
fsigma = 1 -(1 - fsmin)*(k-1)/(Neta - 1.0)
do it = 1,NtIR
	V = 1 + (k-1.0)/(Neta-1.0)
	R = V**(1./3.)
	S = S_NHEJ + (it-1)*dsigma_dt
	S = fsigma*S
	eta = etafun(R,S)
	eta_table(k,it,1:2) = eta
	eta = etafun(R,S_TMEJ)
	eta_table(k,it,3:NP) = eta
enddo
enddo
end subroutine

!--------------------------------------------------------------------------
! f_S = fractional progress through S-phase. =0 in G1, =1 in G2
! tIR = time since IR, in hours.  Since NtIR = 15, it is assumed 
! Two sigma values have been defined: sigma_NHEJ and sigma_TMEJ
! For each pathway, lookup tables for eta as a function of f_S, tIR have been 
! precomputed, requiring only interpolation.
! The number of table entries is Neta*NtIR (= 11*15)
! In eta_table(k, it,..) k is the f_S index, it is the tIR index
! k = f_S*Neta
! it = min(tIR,NtIR)
! NOT USED NOW
!--------------------------------------------------------------------------
function eta_lookup(phase,path,f_S,tIR) result(eta)
integer :: phase, path
real(8) :: f_S, tIR, eta
real(8) :: fk, fi
integer :: k1, k2, it1, it2

if (tIR >= NtIR - 1) then
    it1 = NtIR-1
    it2 = NtIR
    fi = 1
else
    it1 = tIR + 1   ! rounded down
    it2 = it1 + 1
    fi = tIR - it1 + 1
endif
if (phase == G1_phase) then
    k1 = 1
    k2 = 2
    fk = 0
elseif (phase > S_phase) then
    k1 = Neta - 1
    k2 = Neta
    fk = 1.0
else
    k1 = f_S*(Neta-1) + 1
    if (k1 < Neta) then
        k2 = k1+1
        fk = f_S*(Neta - 1) + 1 - k1
    else
        k1 = Neta - 1
        k2 = Neta
        fk = 1.0 
    endif
endif
! fk is the fractional weight of k2, 1-fk is the fractional weight of k1
! fi is the fractional weight of it2, 1-fi is the fractional weight of it1
eta = (1-fk)*(1-fi)*eta_table(k1,it1,path) + fk*(1-fi)*eta_table(k2,it1,path) &
    + (1-fk)*fi*eta_table(k1,it2,path) + fk*fi*eta_table(k2,it2,path)
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine test_eta(S_NHEJ, fsmin)
real(8) :: S_NHEJ, fsmin
real(8) :: f_S, tIR, dt, t, fsigma, V, eta, S, R
real(8) :: St(2,100), Rt(2,100), etat(2,100), eta_array(11)
real(8), parameter :: T_G1 = 6.35, T_S = 9.0, T_G2 = 3.2
integer :: phase, nt, it, k, NHEJfast, iR, iS

! Evaluate eta_Arnould
write(*,*) 'sigma_nhej: ',S_NHEJ
write(*,*) '           S factor'
write(*,'(a,11f8.2)') '    R     ',(1 + (iS-1)/10., iS=1,11)
do iR = 11,1,-1
    R = 0.5 + (iR-1)*0.05
    do iS = 1,11
        S = S_NHEJ*(1 +(iS-1)/10.)
        eta_array(iS) = etafun(R,S) 
    enddo
    write(*,'(f6.2,4x,11f8.5)') R,eta_array
enddo
return

! Compare R, S, eta vs tIR for two cases: IR at start of G1, IR at start of S
dt = 0.2
etat = 0
! G1 case
phase = 1
do it = 1,100
    t = (it-1)*dt
    if (phase == 1 .and. t > T_G1) phase = 2
    if (phase == 2 .and. t > (T_G1 + T_S)) phase = 3
    if (phase == 3 .and. t > (T_G1 + T_S + T_G2)) exit
    if (phase == 1) then
        fsigma = 1
        V = 1
    elseif (phase == 2) then
        fsigma = 1 - (1 - fsmin)*(t - T_G1)/T_S
        V = 1 + (t - T_G1)/T_S
    else
        fsigma = fsmin
        V = 2
    endif
	S = S_NHEJ + dsigma_dt*t
	S = fsigma*S
	R = V**(1./3.)
    eta = etafun(R,S)
    St(1,it) = S
    Rt(1,it) = R
    etat(1,it) = eta
enddo

! S case
phase = 2
do it = 1,100
    t = (it-1)*dt
!    if (phase == 1 .and. t > T_G1) phase = 2
    if (phase == 2 .and. t > (T_S)) phase = 3
    if (phase == 3 .and. t > (T_S + T_G2)) exit
    if (phase == 1) then
        fsigma = 1
        V = 1
    elseif (phase == 2) then
        fsigma = 1 - (1 - fsmin)*t/T_S
        V = 1 + t/T_S
    else
        fsigma = fsmin
        V = 2
    endif
	S = S_NHEJ + dsigma_dt*t
	S = fsigma*S
	R = V**(1./3.)
    eta = etafun(R,S)
    St(2,it) = S
    Rt(2,it) = R
    etat(2,it) = eta
!    write(nflog,'(i4,3f8.3,e12.3)') phase,t,S,R,eta
enddo
do it = 1,100
    if (etat(2,it) == 0) exit
    write(nflog,'(i4,2(2f8.3,e12.3))') it, (St(k,it),Rt(k,it),etat(k,it),k=1,2)
enddo
return

phase = 2
NHEJfast = 1
dt = 0.1
nt = T_S/dt
do it = 1,nt
    tIR = 9.5 + (it-1)*dt
    f_S = (it-1)*dt/T_S
    eta = eta_lookup(phase, NHEJfast, f_S, tIR) 
    write(nflog,'(a,3f8.3)') 'f_S, tIR, eta: ',f_S, tIR, eta
enddo

end subroutine

end module