!!! E. K.H. Lee: Rewrite into fortran 90 of the classic BHMIE f77 code (with Draine edits)
!!! - Orginal history: 
! C Subroutine BHMIE is derived from the Bohren-Huffman Mie scattering
! C     subroutine to calculate scattering and absorption by a homogenous
! C     isotropic sphere.
! C Given:
! C    X = 2*pi*a/lambda
! C    REFREL = (complex refr. index of sphere)/(real index of medium)
! C    NANG = number of angles between 0 and 90 degrees
! C           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
! C           if called with NANG<2, will set NANG=2 and will compute
! C           scattering for theta=0,90,180.
! C Returns:
! C    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
! C                                scatt. E perp. to scatt. plane)
! C    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
! C                                scatt. E parr. to scatt. plane)
! C    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
! C    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
! C    QBACK = 4.*pi*(dC_sca/domega)/pi*a**2
! C          = backscattering efficiency
! C    GSCA = <cos(theta)> for scattering
! C
! C S1 and S2 are the diagonal elements of the "amplitude scattering matrix"
! C (see eq. 3.12 of Bohren & Huffman 1983) -- the off-diagonal elements
! C vanish for a spherical target.
! C For unpolarized incident light, the intensity of scattered light a
! C distance r from the sphere is just
! C          1
! C  I_s = ------ * I_in * S_11
! C        (kr)^2
! C
! C where k=2*pi/lambda 
! C and the "Muller matrix element" S_11 = 0.5*( |S_1|^2 + |S_2|^2 )
! C
! C for incident light polarized perp to the scattering plane,
! C the scattered light is polarized perp to the scattering plane
! C with intensity I_s = I_in * |S_1|^2 / (kr)^2
! C
! C for incident light polarized parallel to the scattering plane,
! C the scattered light is polarized parallel to the scattering plane
! C with intensity I_s = I_in * |S_2|^2 / (kr)^2
! C
! C History:
! C Original program taken from Bohren and Huffman (1983), Appendix A
! C Modified by B.T.Draine, Princeton Univ. Obs., 90.10.26
! C in order to compute <cos(theta)>
! C 91.05.07 (BTD): Modified to allow NANG=1
! C 91.08.15 (BTD): Corrected error (failure to initialize P)
! C 91.08.15 (BTD): Modified to enhance vectorizability.
! C 91.08.15 (BTD): Modified to make NANG=2 if called with NANG=1
! C 91.08.15 (BTD): Changed definition of QBACK.
! C 92.01.08 (BTD): Converted to full double precision and double complex
! C                 eliminated 2 unneed lines of code
! C                 eliminated redundant variables (e.g. APSI,APSI0)
! C                 renamed RN -> EN = double precision N
! C                 Note that DOUBLE COMPLEX and DCMPLX are not part
! C                 of f77 standard, so this version may not be fully
! C                 portable.  In event that portable version is
! C                 needed, use src/bhmie_f77.f
! C 93.06.01 (BTD): Changed AMAX1 to generic function MAX
! C 98.09.17 (BTD): Added variable "SINGLE" and warning in event that
! C                 code is used with single-precision arithmetic (i.e.,
! C                 compiler does not support DOUBLE COMPLEX)
! C 99.02.17 (BTD): Replaced calls to REAL() and IMAG() by
! C                 REALPART() and IMAGPART() for compatibility with g77
! C                 Note that when code is used with standard f77 
! C                 compilers, it is now necessary to enable two lines
! C                 defining functions REALPART(X) and IMAGPART(X)
! C 99.02.19 (BTD): added lines to be enabled to properly define
! C                 REALPART() and IMAGPART() if NOT using g77
! C                 ***see below!!***
! C 01.02.16 (BTD): added IMPLICIT NONE
! C 01.02.27 (BTD): changed definition of QBACK back to convention of
! C                 Bohren & Huffman and others:
! C                 Q_back = 4.*pi*(dC_sca/dOmega)/(pi*a^2) in backward
! C                          direction
! c 02.03.09 (BTD): defined statement function REALPART_SP to
! c                 avoid warning regarding type conversion when taking
! c                 real part of S1(1) to evaluate QEXT
! c                 some cleanup regarding type conversion
! c 02.05.30 (BTD): introduced internal double complex arrays DCXS1,DCXS2
! c                 to possibly increase accuracy during summations.
! c                 After summations, output scattering amplitudes
! c                 via single complex arrays S1,S2 as before.
! c                 Usage of this routine is unaffected by change.
! c                 Note: no longer need statement function REALPART_SP
! c 02.09.18 (BTD): Error in evaluation of QBACK reported by Okada Yasuhiko
! c                 Was calculating QBACK using S1 rather than DCXS1
! c                 Corrected.
! c 02.10.16 (BTD): Added comments explaining definition of S_1 and S_2 .
! C end history
! C
!!! - 
module bhmie_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64
  !integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: pii = 4.0_dp*atan(1.0_dp)

  integer, parameter :: mxnang = 1000
  integer, parameter :: nmxx = 150000

  public :: BHMIE

contains

  subroutine BHMIE(x,refrel,nang,S1,S2,qext,qsca,qback,gsca,ierr)
    implicit none

    integer, intent(in) :: nang
    real(dp), intent(in) :: x
    complex(dp), intent(in) :: refrel

    integer :: ierr
    real(dp), intent(out) :: qext, qsca, qback, gsca
    complex(dp), dimension(2*nang-1), intent(out) :: S1, S2

    integer :: j, jj, n, nstop, nmx, nn
    real(dp) :: chi, chi0, chi1, dang, dx, en, fn, p
    real(dp) :: psi, psi0, psi1, theta, xstop, ymod
    real(dp), dimension(nang) :: amu, pi, pi0, pi1, tau
    complex(dp) :: an, an1, bn, bn1, drefrl, xi, xi1, y
    complex(dp), allocatable, dimension(:) :: d
    complex(dp), dimension(2*nang-1) :: dcxS1, dcxS2


    if (nang > mxnang) then
      ierr = -1
      return
    else if (nang < 2) then
      ierr = -1 ! For ease, nang >= 2 
      return
    end if

    dx = x
    drefrl = refrel
    y = x*drefrl
    ymod = abs(y)

    xstop = x + 4.0_dp*x**(1.0_dp/3.0_dp) + 2.0_dp
    nmx = int(max(xstop,ymod)) + 15

    nstop = int(xstop)

    if (nmx > nmxx) then
      ierr = -2
      return
    end if

    dang = 0.0_dp
    if (nang > 1) then
      dang = 0.5_dp*pii/real(nang-1,dp)
    end if

    do j = 1, nang
      theta = real(j-1,dp)*dang
      amu(j) = cos(theta)
    end do

    do j = 1, nang
      pi0(j) = 0.0_dp
      pi1(j) = 1.0_dp
    end do

    nn = 2*nang-1
    do j = 1, nn
      dcxS1(j) = cmplx(0.0_dp,0.0_dp,dp)
      dcxS1(j) = cmplx(0.0_dp,0.0_dp,dp)
    end do

    allocate(d(nmx))
    d(nmx) = cmplx(0.0_dp,0.0_dp,dp)
    nn = nmx - 1
    do n = 1, nn
      en = real(nmx-n+1,dp)
      d(nmx-n) = (en/y)-(1.0_dp/(d(nmx-n+1) + en/y))
    end do

    psi0 = cos(dx)
    psi1 = sin(dx)
    chi0 = -sin(dx)
    chi1 = cos(dx)
    xi1 = cmplx(psi1,-chi1,dp)
    qsca = 0.0_dp
    gsca = 0.0_dp
    p = -1.0_dp
    do n = 1, nstop
      en = real(n,dp)
      fn = (2.0_dp*en+1.0_dp)/(en*(en+1.0_dp))
      psi = (2.0_dp*en-1.0_dp)*psi1/dx - psi0
      chi = (2.0_dp*en-1.0_dp)*chi1/dx - chi0
      xi = cmplx(psi,-chi,dp)
      if (n > 1) then
        an1 = an
        bn1 = bn
      end if
      an = (d(n)/drefrl+en/dx)*psi - psi1
      an = an/((d(n)/drefrl+en/dx)*xi - xi1)
      bn = (drefrl*d(n)+en/dx)*psi - psi1
      bn = bn/((drefrl*d(n) + en/dx)*xi - xi1)

      qsca = qsca + real((2.0_dp*en-1.0_dp)*(abs(an)**2+abs(bn)**2),dp)
      gsca = gsca + real(((2.0_dp*en-1.0_dp)/(en*(en+1.0_dp))) * &
        & (real(an,dp)*real(bn,dp)+aimag(an)*aimag(bn)),dp)
      if (n > 1) then
        gsca = gsca + real(((en-1.0_dp)*(en+1.0_dp)/en) * & 
          & (real(an1,dp)*real(an,dp) + aimag(an1)*aimag(an) + &
          & real(bn1,dp)*real(bn,dp) + aimag(bn1)*aimag(bn)),dp)
      end if

      do j = 1, nang
        pi(j) = pi1(j)
        tau(j) = en*amu(j)*pi(j)-(en+1.0_dp)*pi0(j)
        dcxS1(j) = dcxS1(j) + fn*(an*pi(j) + bn*tau(j))
        dcxS2(j) = dcxS2(j) + fn*(an*tau(j) + bn*pi(j))
      end do

      p = -p
      do j = 1, nang-1
        jj = 2*nang - j
        dcxS1(jj) = dcxS1(jj) + fn*p*(an*pi(j) - bn*tau(j))
        dcxS2(jj) = dcxS2(jj) + fn*p*(bn*pi(j) - an*tau(j))
      end do

      psi0 = psi1
      psi1 = psi
      chi0 = chi1
      chi1 = chi
      xi1 = cmplx(psi1,-chi1,dp)

      do j = 1, nang
        pi1(j) = ((2.0_dp*en+1.0_dp)*amu(j)*pi(j) - (en+1.0_dp)*pi0(j))/en
        pi0(j) = pi(j)
      end do
    end do

    gsca = real(2.0_dp*gsca/qsca,dp)
    qsca = real((2.0_dp/dx**2)*qsca,dp)
    qext = real((4.0_dp/dx**2)*real(dcxS1(1),dp),dp)
    qback = real(4.0_dp*(abs(dcxS1(2*nang-1))/dx)**2,dp)

    do j = 1, 2*nang-1
      S1(j) = dcxS1(j)
      S2(j) = dcxS2(j)
    enddo

    deallocate(d)

    ierr = 0

  end subroutine BHMIE

end module bhmie_mod
