!!! E. K.H. Lee: Rewrite of the distribution of hollow spheres (DHS) code of M. Min.
!!! - Orginal history: 
! c This subroutine computes the extinction and scattering cross sections
! c of DHS grains with f_max=maxf.
! c e1	real part of the refractive index
! c e2	imaginary part of the refractive index
! c lam	wavelength (micron)
! c rad	volume equivalent radius (micron)
! c cext	extinction cross section (micron^2)
! c csca	scattering cross section (micron^2)
! c maxf	f_max of the DHS (between 0 and 1).
!!! - 
module dhs_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  real(dp), parameter :: eps = 3.0e-14_dp

  ! integer, parameter :: mxang = 200
  ! integer, parameter :: ll = 200000
  ! real(dp) :: toler = 1.0e-6_dp
  ! complex(dp), parameter :: ci = cmplx(0.0_dp,1.0_dp,dp)
  ! complex(dp), parameter :: czero = cmplx(0.0_dp,0.0_dp,dp)

  public :: q_dhs
  private :: gauleg

contains

  subroutine q_dhs(e1, e2, lam, rad, cext, csca, g, maxf, ierr)
    implicit none

    real(dp), intent(in) :: e1, e2, lam, rad, maxf

    integer, intent(out) :: ierr
    real(dp), intent(out) :: cext, csca, g

    integer :: n, i
    real(dp), dimension(50) :: f, wf
    complex(dp) :: m_in, m
    real(dp) :: tot, cabss
    real(dp) :: rad0, qe1, qa1, qs1, wvno, qbs, gqsc, rcore
    real(dp), dimension(10) :: mu
    real(dp), dimension(10,2) :: m1, m2, d21, s21

    mu(1) = 0.1_dp
    mu(2) = 0.4_dp
    cext = 0.0_dp
    csca = 0.0_dp
    cabss = 0.0_dp
    g = 0.0_dp
    m = cmplx(e1,e2,dp)
    tot = 0.0_dp

    n = 20

    call gauleg(0.0_dp,maxf,f,wf,n)

    do i = 1, n
      tot = tot + wf(i)
    end do

    do i = 1, n
      m = cmplx(e1,-e2,dp)
      m_in = cmplx(1.0_dp,0.0_dp,dp)
      wvno = 2.0_dp*pi/lam
      rad0 = rad/((1.0_dp-f(i))**(1.0_dp/3.0_dp))
      rcore = rad0*f(i)**(1.0_dp/3.0_dp)
      ierr = 0
      call DMiLay(rcore, rad0, wvno, m, m_in, mu, &
        & 2, qe1, qs1, qbs, gqsc, &
        & m1, m2, s21, d21, 10 , ierr)
      if (ierr < 0) then
        return
      end if
      qa1 = qe1 - qs1
      if (qe1 <= 0.0_dp) then
        qe1 = 0.0_dp
      end if
      if (qs1 <= 0.0_dp) then
        qs1 = 0.0_dp
      end if
      if (qa1 <= 0.0_dp) then
        qa1 = 0.0_dp
      end if
      cext = cext + pi*rad0**2*wf(i)*qe1/tot
      csca = csca + pi*rad0**2*wf(i)*qs1/tot
      cabss = cabss + pi*rad0**2*wf(i)*qa1/tot
      g = g + pi*rad0**2*wf(i)*gqsc/tot
    end do

    g = g/csca

  end subroutine q_dhs

  subroutine gauleg(x1,x2,x,w,n)
    implicit none
 
    integer, intent(in) :: n
    real(dp), intent(in) :: x1, x2

    real(dp), dimension(n), intent(out) :: x, w 

    integer :: i, j, m, ic
    real(dp) :: p1, p2, p3, pp, xl, xm, z, z1
    logical :: conv

    m = (n+1)/2
    xm = 0.5_dp*(x2+x1)
    xl = 0.5_dp*(x2-x1)

    do i = 1, m
      z = cos(pi*(real(i,dp)-0.25_dp)/(real(n,dp)+0.5_dp))
      ic = 0
      do while (ic == 0)
        p1 = 1.0_dp
        p2 = 0.0_dp
        do j = 1, n
          p3 = p2
          p2 = p1
          p1 = ((2.0_dp*real(j,dp)-1.0_dp)*z*p2 - (real(j,dp)-1.0_dp)*p3)/real(j,dp)
        end do
        pp = real(n,dp)*(z*p1-p2)/(z*z - 1.0_dp)
        z1 = z
        z = z1-p1/pp
        if (abs(z-z1) < eps) then
          ic = 1
        end if
      end do
      x(i) = xm-xl*z
      x(n+1-i) = xm+xl*z
      w(i) = 2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
      w(n+1-i) = w(i)
    end do

  end subroutine gauleg
  
end module dhs_mod
