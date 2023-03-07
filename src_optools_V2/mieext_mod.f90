!!! E. K.H. Lee: Rewrite into fortran 90 of the Grainger AOPP EDOG IDL/Fortran 77 Mie code
!!! - Orginal history: 
!     General purpose Mie scattering routine for single particles
!
!     Author: R Grainger 1990
!
!     History:
!       G Thomas, March 2005: Added calculation of Phase function and
!         code to ensure correct calculation of backscatter coeficient
!       G Thomas, November 2006: Calculation of backscatter efficiency
!         now done using the A and B values rather than the intensity at
!         180 degrees.
!       G Thomas, July 2012: Removed the maximum size parameter
!         restriction (Imaxx & Itermax set to 12e6). Also all integers
!         changed to 32-bit (* 4) to avoid overflow.
!       G Thomas/D Grainger, 2 Aug 2012: Bug fix in backscatter
!       calculation
!       G McGarragh, 29, Jul 2015: Add support to output the 11, 33, 12,
!         and 34 elements of the 4x4 single scattering phase matrix F,
!         where F21 = F12 and F43 = -F34, which is all that is required
!         for randomly oriented spherical particles.  F11 is identical
!         to the old phase function output DPh.

module mieext_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64
  !integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: Imaxx = 12000
  real(dp), parameter :: RImax = 2.5_dp ! Largest real part of refractive index
  integer, parameter :: Itermax = Imaxx * RImax ! must be large enough to cope with the largest possible nmx = x *real(scm) + 15
                                                ! or nmx =  Dx + 4.05*Dx**(1./3.) + 2.0

  public :: MieExt

contains

  subroutine MieExt(Dx, SCm, Dqxt, Dqsc, Dg, Error)

    implicit none

    real(dp), intent(in) :: Dx
    complex(dp), intent(in) :: SCm

    integer, intent(out) :: Error
    real(dp), intent(out) :: Dqxt, Dqsc, Dg

    integer :: Nstop, NmX, n
    real(dp) :: Chi, Chi0, Chi1, APsi, APsi0, APsi1
    real(dp) :: Psi, Psi0, Psi1
    complex(dp) :: Ir, Cm, A, ANM1, B, BNM1, Xi, Xi0, Xi1, Y
    complex(dp), allocatable, dimension(:) :: D

    real(dp) :: Tnp1, Tnm1
    real(dp) :: Dn
    real(dp) :: Rnx, Turbo, A2
    complex(dp) :: A1

    if (Dx > Imaxx) then
      Error = -1
      return
    end if

    Cm = SCm
    Ir = 1.0_dp/Cm
    Y = Dx * Cm

    if (Dx < 0.02_dp) then
      Nstop = 2
    else
      if (Dx <= 8.0_dp) then
        Nstop = int(Dx + 4.0_dp*Dx**(1.0_dp/3.0_dp) + 2.0_dp)
      else
        if (Dx < 4200.0_dp) then
          NStop = int(Dx + 4.05_dp*Dx**(1.0_dp/3.0_dp) + 2.0_dp)
        else
          Nstop = int(Dx + 4.0_dp*Dx**(1.0_dp/3.0_dp) + 2.0_dp)
        end if
      end if
    end if

    NmX = int(max(real(NStop,dp),real(abs(Y),dp)) + 15.0_dp)
    if (NmX > Itermax) then
      Error = -2
      return
    end if


    allocate(D(NmX))
    D(NmX) = cmplx(0.0_dp,0.0_dp,dp)
    do n = NmX-1,1,-1
      A1 = (real(n,dp)+1.0_dp)/Y
      D(n) = A1 - 1.0_dp/(A1 + D(n+1))
    end do

    Psi0 = cos(Dx)
    Psi1 = sin(Dx)
    Chi0 = -sin(Dx)
    Chi1 = cos(Dx)
    APsi0 = Psi0
    APsi1 = Psi1
    Xi0 = cmplx(APsi0,Chi0,dp)
    Xi1 = cmplx(Apsi1,Chi1,dp)
    Dqsc = 0.0_dp
    Dqxt = 0.0_dp
    Dg = 0.0_dp
    Tnp1 = 1.0_dp

    do n = 1, Nstop
      DN = real(n,dp)
      Tnp1 = Tnp1 + 2.0_dp
      Tnm1 = Tnp1 - 2.0_dp
      A2 = Tnp1 / (DN*(DN + 1.0_dp))
      Turbo = (DN + 1.0_dp)/DN
      Rnx = DN/Dx
      Psi = Tnm1*Psi1/Dx - Psi0
      APsi = Psi
      Chi = Tnm1*Chi1/Dx - Chi0
      Xi = cmplx(Apsi, Chi,dp)
      A = ((D(n)*Ir+Rnx)*APsi-APsi1) / ((D(n)*Ir+Rnx)*Xi-Xi1)
      B = ((D(n)*Cm+Rnx)*APsi-APsi1) / ((D(n)*Cm+Rnx)*Xi-Xi1)  
      Dqxt = Dqxt + Tnp1 * real(A + B,dp)
      Dqsc = Dqsc + Tnp1 * (A*conjg(A) + B*conjg(B))
      if (n > 1) then
        Dg = Dg + (DN*DN - 1.0_dp) * real(ANM1*conjg(A) + BNM1 * conjg(B),dp) / &
          & DN + Tnm1 * real(ANM1*conjg(BNM1),dp) / (DN*DN - DN)
      end if
      Anm1 = A
      Bnm1 = B
      Psi0 = Psi1
      Psi1 = Psi
      Apsi1 = Psi1
      Chi0 = Chi1
      Chi1 = Chi
      Xi1 = cmplx(APsi1,Chi1,dp)
    end do
   
    if (Dg > 0.0_dp) then
      Dg = 2.0_dp * Dg / Dqsc
    end if 
    Dqsc = 2.0_dp * Dqsc / Dx**2
    Dqxt = 2.0_dp * Dqxt / Dx**2 

    deallocate(D)

    Error = 0

  end subroutine MieExt

end module mieext_mod
