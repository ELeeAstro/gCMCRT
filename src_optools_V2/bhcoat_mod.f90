!!! E. K.H. Lee: Rewrite of the bhcoat (coated spheres) classic routine (with Draine edits).
!!! - Orginal history: 
!***********************************************************************
!
! Subroutine BHCOAT calculates Q_ext, Q_sca, Q_back, g=<cos> 
! for coated sphere.
! All bessel functions computed by upward recurrence.
! Input:
!        X = 2*PI*RCORE*REFMED/WAVEL
!        Y = 2*PI*RMANT*REFMED/WAVEL
!        RFREL1 = REFCOR/REFMED
!        RFREL2 = REFMAN/REFMED 
! where  REFCOR = complex refr.index of core)
!        REFMAN = complex refr.index of mantle)
!        REFMED = real refr.index of medium)
!        RCORE = radius of core
!        RMANT = radius of mantle
!        WAVEL = wavelength of light in ambient medium

! returns:
!        QQEXT = C_ext/pi*rmant^2
!        QQSCA = C_sca/pi*rmant^2
!        QBACK = 4*pi*(dQ_sca/dOmega)
!              = "radar backscattering efficiency factor"
!        GSCA  = <cos(theta)> for scattered power
!
! Routine BHCOAT is taken from Bohren & Huffman (1983)
! extended by Prof. Francis S. Binkowski of The University of North
! Carolina at Chapel Hill to evaluate GSCA=<cos(theta)>
! History:
! 92.11.24 (BTD) Explicit declaration of all variables
! 00.05.16 (BTD) Added IMPLICIT NONE
! 12.04.10 (FSB) Modified by Prof. Francis S. Binkowski of
!                The University of North Carolina at Chapel Hill
!                to evaluate GSCA=<cos(theta)>
! 12.06.15 (BTD) Cosmetic changes
!***********************************************************************
!!! -
module bhcoat_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: del = 1.0e-8_dp
  complex(dp), parameter :: ii = cmplx(0.0_dp,1.0_dp,dp)

  public :: BHCOAT

contains

  subroutine BHCOAT(xx, yy, rrfrl1, rrfrl2, qqext, qqsca, qback, gsca)
    implicit none

    real(dp), intent(in) :: xx, yy
    complex(dp), intent(in) :: rrfrl1, rrfrl2

    real(dp), intent(out) :: qqext, qqsca, qback, gsca

    integer :: iflag, n, nstop
    real(dp) :: chi0y, chi1y, chiy, en, psi0y, psi1y, psiy
    real(dp) :: qext, qsca, rn, x, y, ystop
    complex(dp) :: amess1, amess2, amess3, amess4, an, an1, ancap
    complex(dp) :: bn, bn1 ,bncap, brack
    complex(dp) :: chi0x2, chi0y2, chi1x2, chi1y2, chix2, chipx2, chipy2, chiy2, crack
    complex(dp) :: d0x1, d0x2 ,d0y2, d1x1, d1x2, d1y2, dnbar, gnbar
    complex(dp) :: refrel, rfrel1, rfrel2
    complex(dp) :: xback, xi0y, xi1y, xiy
    complex(dp) :: x1, x2, y2

    x = xx
    y = yy
    rfrel1 = rrfrl1
    rfrel2 = rrfrl2    

    ! -----------------------------------------------------------
    ! del is the inner sphere convergence criterion
    ! -----------------------------------------------------------
    x1 = rfrel1*x
    x2 = rfrel2*x
    y2 = rfrel2*y
    ystop = y + 4.0_dp*y**(1.0_dp/3.0_dp) + 2.0_dp
    refrel = rfrel2/rfrel1
    nstop = int(ystop) + 1

    ! -----------------------------------------------------------
    ! series terminated after nstop terms
    ! -----------------------------------------------------------
    d0x1 = cos(x1)/sin(x1)
    d0x2 = cos(x2)/sin(x2)
    d0y2 = cos(y2)/sin(y2)
    psi0y = cos(y)
    psi1y = sin(y)
    chi0y = -sin(y)
    chi1y = cos(y)
    xi0y = psi0y - ii*chi0y
    xi1y = psi1y-ii*chi1y
    chi0y2 = -sin(y2)
    chi1y2 = cos(y2)
    chi0x2 = -sin(x2)
    chi1x2 = cos(x2)
    qsca = 0.0_dp
    qext = 0.0_dp
    gsca = 0.0_dp
    xback = cmplx(0.0_dp,0.0_dp,dp)
    iflag = -1

    do n = 1, nstop
      rn = real(n,dp)
      en = rn
      psiy = (2.0_dp*rn-1.0_dp)*psi1y/y - psi0y
      chiy = (2.0_dp*rn-1.0_dp)*chi1y/y - chi0y
      xiy = psiy - ii*chiy
      d1y2 = 1.0_dp/(rn/y2-d0y2) - rn/y2

      if (iflag == -1) then

        ! calculate inner sphere ancap, bncap
        ! and brack and crack

        d1x1 = 1.0_dp/(rn/x1-d0x1) - rn/x1
        d1x2 = 1.0_dp/(rn/x2-d0x2) - rn/x2
        chix2 = (2.0_dp*rn-1.0_dp)*chi1x2/x2 - chi0x2
        chiy2 = (2.0_dp*rn-1.0_dp)*chi1y2/y2 - chi0y2
        chipx2 = chi1x2 - rn*chix2/x2
        chipy2 = chi1y2 - rn*chiy2/y2
        ancap = refrel*d1x1 - d1x2
        ancap = ancap/(refrel*d1x1*chix2-chipx2)
        ancap = ancap/(chix2*d1x2-chipx2)
        brack = ancap*(chiy2*d1y2-chipy2)
        bncap = refrel*d1x2 - d1x1
        bncap = bncap/(refrel*chipx2-d1x1*chix2)
        bncap = bncap/(chix2*d1x2-chipx2)
        crack = bncap*(chiy2*d1y2-chipy2)

        ! calculate convergence test expressions for inner sphere
        ! see pp 483-485 of bohren & huffman for definitions

        amess1 = brack*chipy2
        amess2 = brack*chiy2
        amess3 = crack*chipy2
        amess4 = crack*chiy2

      endif ! test on iflag.eq.0

      ! now test for convergence for inner sphere
      ! all four criteria must be satisfied
      ! see p 484 of Bohren & Huffman

      if (abs(amess1) < del*abs(d1y2) .and. &
        & abs(amess2) < del .and. &
        & abs(amess3) < del*abs(d1y2).and. &
        & abs(amess4) < del) then
        ! convergence for inner sphere
        brack = cmplx(0.0_dp,0.0_dp,dp)
        crack = cmplx(0.0_dp,0.0_dp,dp)
        iflag = 0
      else
        ! no convergence yet
        iflag = -1
      endif

      dnbar = d1y2 - brack*chipy2
      dnbar = dnbar/(1.0_dp-brack*chiy2)
      gnbar = d1y2 - crack*chipy2
      gnbar = gnbar/(1.0_dp-crack*chiy2)

      ! store previous values of an and bn for use in computation of 
      ! g=<cos(theta)>

      if (n > 1) then
        an1 = an
        bn1 = bn
      endif

      ! update an and bn
      an = (dnbar/rfrel2+rn/y)*psiy - psi1y
      an = an/((dnbar/rfrel2+rn/y)*xiy-xi1y)
      bn = (rfrel2*gnbar+rn/y)*psiy - psi1y
      bn = bn/((rfrel2*gnbar+rn/y)*xiy-xi1y)

      ! calculate sums for qsca,qext,xback
      qsca = qsca+(2.0_dp*rn+1.0_dp)*(abs(an)*abs(an)+abs(bn)*abs(bn))
      xback = xback+(2.0_dp*rn+1.0_dp)*(-1.0_dp)**n*(an-bn)
      qext = qext+(2.0_dp*rn+1.0_dp)*(real(an,dp)+real(bn,dp))

      ! (FSB) calculate the sum for the asymmetry factor

      gsca = gsca+((2.0_dp*en+1.0_dp)/(en*(en+1.0_dp))) * &
        &  (real(an,dp)*real(bn,dp)+aimag(an)*aimag(bn))
      if (n > 1) then
        gsca = gsca+((en-1.0_dp)*(en+1.0_dp)/en) * &
          & (real(an1,dp)*real(an,dp)+aimag(an1)*aimag(an) + &
          & real(bn1,dp)*real(bn,dp)+aimag(bn1)*aimag(bn))
      endif

      ! continue update for next iteration

      psi0y = psi1y
      psi1y = psiy
      chi0y = chi1y
      chi1y = chiy
      xi1y = psi1y - ii*chi1y
      chi0x2 = chi1x2
      chi1x2 = chix2
      chi0y2 = chi1y2
      chi1y2 = chiy2
      d0x1 = d1x1
      d0x2 = d1x2
      d0y2 = d1y2

    end do

    ! have summed sufficient terms
    ! now compute QQSCA,QQEXT,QBACK, and GSCA

    qqsca = (2.0_dp/(y*y))*qsca
    qqext = (2.0_dp/(y*y))*qext
    qback = (abs(xback))**2
    qback = (1.0_dp/(y*y))*qback
    gsca = 2.0_dp*gsca/qsca  
    
  end subroutine BHCOAT

end module bhcoat_mod
