module mie_approx_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

  public :: adt, rayleigh, rayleigh_gans, geo_optics

contains

  subroutine adt(x, ri, q_abs, q_sca, q_ext)
    implicit none

    real(dp), intent(in) :: x
    complex(dp), intent(in) :: ri

    real(dp), intent(out) :: q_abs, q_sca, q_ext

    real(dp) :: rho1, rho2, rho0, beta0, beta, fac, fac2
    complex(dp) ::  rho

    rho = 2.0_dp * x * (ri - 1.0_dp)
    rho0 = abs(rho)
    rho1 = real(rho, dp)
    rho2 = aimag(rho)

    if (abs(rho1) > 0.0_dp) then
      beta0 = atan(abs(rho2)/abs(rho1))
      if (rho1 < 0.0_dp .and. rho2 > 0.0_dp) then
        beta = pi - beta0
      else if (rho1 < 0.0_dp .and. rho2 < 0.0_dp) then
        beta = pi + beta0
      else if (rho1 > 0.0_dp .and. rho2 < 0.0_dp) then 
        beta = 2.0_dp*pi - beta0
      else 
        beta = beta0 
      endif     
    else
      if (rho2 > 0.0_dp) then
        beta = 0.5_dp*pi 
      else
        beta = 1.5_dp*pi
      endif
    end if

    if (rho0 < 1.0e-3_dp) then
      q_ext = (4.0_dp/3.0_dp)*rho2 + 0.5_dp*(rho1**2 - rho2**2)
      q_abs = (4.0_dp/3.0_dp)*rho2 - rho2**2
      q_sca = 0.5_dp*rho0**2
    else
      fac = exp(-rho2)
      fac2 = fac**2
      q_ext = 2.0_dp + (4.0_dp/rho0**2)*(cos(2.0_dp*beta) - fac*(cos(rho1 - 2*beta) + rho0*sin(rho1 - beta)))
      q_abs = 1.0_dp + fac2/rho2 + (fac2 - 1.0_dp)/(2.0_dp*rho2**2)
      q_sca = q_ext - q_abs
    end if

    if (x >= 100.0_dp) then
      q_ext = q_ext + 2.0_dp * x**(-2.0_dp/3.0_dp)
    end if

  end subroutine adt

  ! Small dielectric sphere approximation, small particle limit x << 1
  subroutine rayleigh(x, ri, q_abs, q_sca, q_ext)
    implicit none

    real(dp), intent(in) :: x
    complex(dp), intent(in) :: ri

    real(dp), intent(out) :: q_abs, q_sca, q_ext

    complex(dp) :: alp

    alp = (ri**2 - 1.0_dp)/(ri**2 + 2.0_dp)

    q_sca = 8.0_dp/3.0_dp * x**4 * abs(alp)**2
    q_ext = 4.0_dp * x * aimag(alp * (1.0_dp + x**2/15.0_dp*alp * ((ri**4+27.0_dp*ri**2+38.0_dp)/(2.0_dp*ri**2+3.0_dp)))) + &
     & 8.0_dp/3.0_dp * x**4 * real(alp**2,dp)

    q_abs = q_ext - q_sca

  end subroutine rayleigh

  subroutine rayleigh_gans(x, ri, q_abs, q_sca, q_ext)
    implicit none

    real(dp), intent(in) :: x
    complex(dp), intent(in) :: ri

    real(dp), intent(out) :: q_abs, q_sca, q_ext

    real(dp) :: phix

    if (x <= 1.0_dp) then
      phix = 1.185_dp*x**4*(1.0_dp - 0.4_dp*x**2 + 0.096_dp*x**4) 
    else if (x > 1.0_dp .and. x <= 2.0_dp) then
      phix = 1.92_dp*x**2 - 1.084*x ! Note, think typo in textbooks, x**2 fits better than x**4
    else if (x > 2.0_dp .and. x <= 12.5_dp) then
      phix = 2.122_dp*x**2 - 1.456_dp*x
    else
      phix = 2.0_dp*x**2
    end if

    q_abs = 8.0_dp/3.0_dp * x * aimag(ri - 1.0_dp)
    q_sca = phix *  abs(ri - 1.0_dp)**2

    q_ext = q_sca + q_abs

  end subroutine rayleigh_gans

  subroutine geo_optics(x, ri, q_abs, q_sca, q_ext)
    implicit none

    real(dp), intent(in) :: x
    complex(dp), intent(in) :: ri

    real(dp), intent(out) :: q_abs, q_sca, q_ext

    q_ext = 2.0_dp*(1.0_dp + x**(-2.0_dp/3.0_dp))
    q_abs = 1.0_dp
    q_sca = 1.0_dp

  end subroutine geo_optics

end module mie_approx_mod
