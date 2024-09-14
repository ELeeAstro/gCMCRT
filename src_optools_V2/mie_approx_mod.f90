module mie_approx_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  real(dp), parameter :: gam = 0.577215665_dp

  public :: adt, rayleigh, rayleigh_gans, geo_optics

contains

  !! anomalous diffraction theory (ADT) valid for x >> 1, |m - 1| << 1
  !! We follow the Draine et al. formulation
  !! (Also known as the van de Hulst (1957) approximation)
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

    if (rho0 < 1.0e-1_dp) then
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

    if (x >= 10.0_dp) then
      q_ext = q_ext + 2.0_dp * x**(-2.0_dp/3.0_dp)
    end if

    ! First Mie bump only:
    !q_ext = q_ext * (1.1_dp + (real(ri,dp) - 1.2_dp)/3.0_dp)

  end subroutine adt

  !! Rayleigh scattering regime x << 1, |mx| << 1
  !! i.e. small particles with minimal field changes
  subroutine rayleigh(x, ri, q_abs, q_sca, q_ext)
    implicit none

    real(dp), intent(in) :: x
    complex(dp), intent(in) :: ri

    real(dp), intent(out) :: q_abs, q_sca, q_ext

    complex(dp) :: alp

    alp = (ri**2 - 1.0_dp)/(ri**2 + 2.0_dp)

    !! Follow Bohren and Huffman (1983) approximations
    q_sca = 8.0_dp/3.0_dp * x**4 * abs(alp)**2
    q_abs = 4.0_dp * x * &
      & aimag(alp * (1.0_dp + x**2/15.0_dp*alp * ((ri**4+27.0_dp*ri**2+38.0_dp)/(2.0_dp*ri**2+3.0_dp))))

    q_ext = q_abs + q_sca

  end subroutine rayleigh

  !! Rayleigh-Gans approximation - uses fitting function for phi(x)
  !! valid: |m - 1| << 1, 2x(m - 1) << 1 
  !! i.e. `soft particles with not too big x and phase shifts'
  subroutine rayleigh_gans(x, ri, q_abs, q_sca, q_ext, g)
    implicit none

    real(dp), intent(in) :: x
    complex(dp), intent(in) :: ri

    real(dp), intent(out) :: q_abs, q_sca, q_ext, g

    real(dp) :: phix, Ci, hx, x4

    !! Shifrin and Ston (1976) function for phi(x)
    ! if (x <= 1.0_dp) then
    !   phix = 1.185_dp*x**4*(1.0_dp - 0.4_dp*x**2 + 0.096_dp*x**4) 
    ! else if (x > 1.0_dp .and. x <= 2.0_dp) then
    !   phix = 1.92_dp*x**2 - 1.084*x ! Note, think typo in textbooks, x**2 fits better than x**4
    ! else if (x > 2.0_dp .and. x <= 12.5_dp) then
    !   phix = 2.122_dp*x**2 - 1.456_dp*x
    ! else
    !   phix = 2.0_dp*x**2
    ! end if

    x4 = 4.0_dp * x

    !! Attempt using the cosine integral
    call cisib(x4, Ci)
    phix = 2.5_dp + 2.0_dp*x**2 - (sin(x4)/(x4)) - &
      & ((7.0_dp*(1.0_dp - cos(x4)))/(16.0_dp*x**2)) + &
      & (1.0_dp/(2.0_dp*x**2) - 2.0_dp) * (gam + log(x4) - Ci) 

    q_ext = phix *  abs(ri - 1.0_dp)**2

    q_abs = 8.0_dp/3.0_dp * x * aimag(ri)

    q_sca = q_ext - q_abs

    !! Calculate (estimate) the asymmetry factor
    if (x < 1e-2_dp) then
      !! Avoid numerical error
      g = 0.0_dp
    else
      hx = 4.0_dp/x*4 * (((9.0_dp/128.0_dp) - (5.0_dp*x**2/64.0_dp)) * (x4*sin(x4)+cos(x4)) + &
        & x**6/2.0_dp + 11.0_dp*x**4/8.0_dp - 31.0_dp*x**2/64.0_dp - 9.0_dp/128_dp + &
        & x**2*(x**2 - 3.0_dp/8.0_dp)*(Ci - gam - log(x4)))

      g = min(hx/phix,1.0_dp) ! Limit maximum to 1

    end if

  end subroutine rayleigh_gans

  !! Geometric optics regime x >> 1, 2x(m - 1) >> 1
  !! i.e. large particles with large phase shifts
  subroutine geo_optics(x, ri, q_abs, q_sca, q_ext)
    implicit none

    real(dp), intent(in) :: x
    complex(dp), intent(in) :: ri

    real(dp), intent(out) :: q_abs, q_sca, q_ext

    real(dp) :: n, k, z, t, rho, c, b , TT, Sn, S1, S2, S3, phin, q_e

    n = real(ri,dp)
    k = aimag(ri)

    ! if (n > 3.0_dp) then
    !   !! More complex theory 
    !   !! `Hard' particles only (n = 1.2-1.55)
    !   z = log(k)
    !   c = 4.0_dp * k * x
    !   b = sqrt(1.0_dp - n**(-2))
    !   rho = 2.0_dp*x*(n - 1.0_dp)
    !   t = (21.2_dp + 20.1_dp*z + 11.1*z**2 + z**3)**(-1)
    !   TT = 1.0_dp + (n - 1.0_dp)*(1.0_dp - exp(-1.0_dp/(t*rho)))
    !   S1 = (8.0_dp * n**4 * (n**4 + 1.0_dp))/((n**4 - 1.0_dp)**2*(n**2 + 1.0_dp))
    !   S2 = (n**2*(n**2 - 1.0_dp)**2)/((n**2 + 1.0_dp)**3)
    !   S3 = (3.0_dp*n**7 - 7.0_dp*n*6 - 13.0_dp*n**5 - 9.0_dp*n**4 - 7.0_dp*n**3 - 3.0*n**2 - n - 1.0_dp) / &
    !     & (3.0_dp *(n**4 - 1.0_dp)*(n**2 + 1.0_dp)*(n + 1.0_dp))
    !   Sn = S1 * log(n) - S2 * log((n + 1.0_dp)/(n - 1.0_dp)) + S3
    !   q_abs = TT * (1.0_dp - 2.0*n**2/c**2 * (exp(-c*b)*(1.0_dp + c*b) - exp(-c)*(1.0_dp + c)) - &
    !     & Sn * (1.0_dp - exp(-c*b))**2)
    ! else
      !! complex angular momentum theory approximation for q_abs
      !! `Soft' particles only (n < 1.2)
      b = sqrt(1.0_dp - n**(-2))
      phin = n**2 * (1.0_dp - b**3)
      c = 4.0_dp * k * x
      !q_e = n*c*(acos(1.0_dp/n) - 1.0_dp/n*sqrt(1.0_dp - (1.0_dp/n)**2))
      !q_abs = 2.0_dp/3.0_dp * c * phin + q_e

      q_abs = 1.0_dp - exp(-2.0_dp/3.0_dp * phin * c) ! Shifrin and Tonna (1992) approximation
    !end if

    !! Fringing (edge effects) approximation
    q_ext = 2.0_dp*(1.0_dp + x**(-2.0_dp/3.0_dp))

    q_sca = q_ext - q_abs

  end subroutine geo_optics

  !! Cosine integral calculation adpated from special_functions.f90 (cisib)
  !! https://people.math.sc.edu/Burkardt/f_src/special_functions/special_functions.f90
  !! Shanjie Zhang, Jianming Jin, Computation of Special Functions, Wiley, 1996
  subroutine cisib(x, ci)
    implicit none

    real(dp), intent(in) :: x

    real(dp), intent(out) :: ci

    real(dp) :: fx, gx, x2

    x2 = x * x

    if (x == 0.0_dp) then

      ci = -1.0e300_dp

    else if (x <= 1.0_dp) then

      ci = (((( -3.0e-8_dp         * x2 &
        &       + 3.10e-6_dp     ) * x2 &
        &       - 2.3148e-4_dp   ) * x2 &
        &       + 1.041667e-2_dp ) * x2 &
        &       - 0.25_dp        ) * x2 &
        & + 0.577215665_dp + log(x)

    else

      fx = (((( x2               &
        & + 38.027264_dp  ) * x2 &
        & + 265.187033_dp ) * x2 &
        & + 335.67732_dp  ) * x2 &
        & + 38.102495_dp  ) /    &
        & (((( x2                &
        & + 40.021433_dp  ) * x2 &
        & + 322.624911_dp ) * x2 &
        & + 570.23628_dp  ) * x2 &
        & + 157.105423_dp )
  
      gx = (((( x2                &
        & + 42.242855_dp  ) * x2  &
        & + 302.757865_dp ) * x2  &
        & + 352.018498_dp ) * x2  &
        & + 21.821899_dp ) /      &
        & (((( x2                 &
        & + 48.196927_dp   ) * x2 &
        & + 482.485984_dp  ) * x2 &
        & + 1114.978885_dp ) * x2 &
        & + 449.690326_dp  ) / x

      ci = fx * sin(x) / x - gx * cos(x) / x

    end if

  end subroutine cisib

end module mie_approx_mod
