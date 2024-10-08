module cloud_tables_dist
  use optools_data_mod
  use optools_aux, only : trapz
  use cloud_tables_mie, only : cl_mie
  use ieee_arithmetic
  implicit none

  private
  public :: dist_cl

contains

  subroutine dist_cl(z, l, eps_comb, cl_out_k, cl_out_a, cl_out_g)
    implicit none

    integer, intent(in) :: z, l
    complex(kind=dp), intent(in) :: eps_comb
    real(kind=dp), intent(out) :: cl_out_k, cl_out_a, cl_out_g
    real(kind=dp), dimension(ndist) :: nd_dist, ifunc_k, ifunc_a, ifunc_g

    integer :: m, idist_m
    real(kind=dp), dimension(3) :: cl_out_k2, cl_out_a2, cl_out_g2
    real(kind=dp) :: beta, alpha,  aeff, const, Ev, Var, muu, sigg, lam


    Ev = a_cl_lay(z)
    Var = var_cl_lay(z)

    idist_m = idist

    select case(idist_m)
    case(1)

      ! Single particle size (idist = 1, nmode = 1)
      call cl_mie(l,nd_cl_lay(z),a_cl_lay(z),eps_comb,cl_out_k,cl_out_a,cl_out_g)

      !print*, l, nd_cl_lay(z),a_cl_lay(z),eps_comb,cl_out_k,cl_out_a,cl_out_g

    case(2)

      ! Small 3 size peaked size distribution @ 1% delta a to smooth Mie ressonance bumps
      call cl_mie(l,nd_cl_lay(z),a_cl_lay(z)*0.99_dp,eps_comb,cl_out_k2(1),cl_out_a2(1),cl_out_g2(1))
      call cl_mie(l,nd_cl_lay(z),a_cl_lay(z),eps_comb,cl_out_k2(2),cl_out_a2(2),cl_out_g2(2))
      call cl_mie(l,nd_cl_lay(z),a_cl_lay(z)*1.01_dp,eps_comb,cl_out_k2(3),cl_out_a2(3),cl_out_g2(3))

      cl_out_k = sum(cl_out_k2(:))/3.0_dp
      cl_out_a = sum(cl_out_a2(:))/3.0_dp
      cl_out_g = sum(cl_out_g2(:))/3.0_dp

    case(3)

      !! log normal distribution - particle size in prf is the median/geometric mean size
      !! sig = std. deviation - lsig = ln(sig), typically 1 < lsig < 2
      !! Here, use Ev as the median grain size and Var as the sigma

      do m = 1, ndist
        ! Distribution in cm-3 cm-1
        nd_dist(m) = (nd_cl_lay(z)  / (a_dist(m) * sqrt(twopi) * log(Var))) * &
          & exp(-(log(a_dist(m)/Ev)**2/(2.0_dp * log(Var)**2)))

        if ((ieee_is_nan(nd_dist(m)) .eqv. .True.) .or. (ieee_is_finite(nd_dist(m)) .eqv. .False.)) then
          nd_dist(m) = 1.0e-99_dp
        end if

        ! Limiter for very low numbers
        nd_dist(m) = max(nd_dist(m),1.0e-99_dp)

        ! Call mie theory routine for this distribution point
        call cl_mie(l,nd_dist(m),a_dist(m),eps_comb,ifunc_k(m),ifunc_a(m),ifunc_g(m))

        !print*, l, m, nd_dist(m), ,a_dist(m), 

      end do

    case(4)

      !! Gamma distribution - particle size in prf sets the parameters of the distribution
      !! Use an eff_fac to give varience as width of mean particle size, typically 0 < eff_fac << 1 (~0.1)

      alpha = Ev**2/Var
      beta = Ev/Var

      const = nd_cl_lay(z) * (beta**(alpha)/gamma(alpha))

      do m = 1, ndist
        ! Distribution in cm-3 cm-1
        nd_dist(m) = const * a_dist(m)**(alpha-1.0_dp) * exp(-beta*a_dist(m))

        if ((ieee_is_nan(nd_dist(m)) .eqv. .True.) .or. (ieee_is_finite(nd_dist(m)) .eqv. .False.)) then
          nd_dist(m) = 1.0e-99_dp
        end if

        ! Limiter for very low numbers
        nd_dist(m) = max(nd_dist(m),1.0e-99_dp)

        ! Call mie theory routine for this distribution point
        call cl_mie(l,nd_dist(m),a_dist(m),eps_comb,ifunc_k(m),ifunc_a(m),ifunc_g(m))

      end do

    case(5)

      !! Inverse-Gamma distribution - particle size in prf sets the paramaters of the distribution
      !! Use an eff_fac to give varience as width of mean particle size, typically 0 < eff_fac << 1 (~0.1)

      alpha = Ev**2/Var + 2.0_dp
      beta = Ev*(alpha - 1.0_dp)

      const = nd_cl_lay(z) * (beta**(alpha)/gamma(alpha))

      do m = 1, ndist
        ! Distribution in cm-3 cm-1
        nd_dist(m) = const * (1.0_dp/a_dist(m))**(alpha+1.0_dp) * exp(-beta/a_dist(m))

        if ((ieee_is_nan(nd_dist(m)) .eqv. .True.) .or. (ieee_is_finite(nd_dist(m)) .eqv. .False.)) then
          nd_dist(m) = 1.0e-99_dp
        end if

        ! Limiter for very low numbers
        nd_dist(m) = max(nd_dist(m),1.0e-99_dp)

        ! Call mie theory routine for this distribution point
        call cl_mie(l,nd_dist(m),a_dist(m),eps_comb,ifunc_k(m),ifunc_a(m),ifunc_g(m))

      end do

    case(6)

      !! Rayleigh distribution - particle size in prf sets the sigma of the distribution

      !! sig is directly related to the distribution mean or varience
      sig = Ev/sqrt(pi/2.0_dp)
      !sig = sqrt(Var/(2.0_dp - pi/2.0_dp))

      do m = 1, ndist
        ! Distribution in cm-3 cm-1
        nd_dist(m) = nd_cl_lay(z) * (a_dist(m)/sig**2) * exp(-a_dist(m)**2/(2.0_dp * sig**2))

        if ((ieee_is_nan(nd_dist(m)) .eqv. .True.) .or. (ieee_is_finite(nd_dist(m)) .eqv. .False.)) then
          nd_dist(m) = 1.0e-99_dp
        end if

        ! Limiter for very low numbers
        nd_dist(m) = max(nd_dist(m),1.0e-99_dp)

        ! Call mie theory routine for this distribution point
        call cl_mie(l,nd_dist(m),a_dist(m),eps_comb,ifunc_k(m),ifunc_a(m),ifunc_g(m))

      end do

    case(7)

      !! Hansen distribution - particle size in prf sets the effective size of the distribution
      !! According to Hansen - 0 < veff < 0.5 typically

      !! variables related to the effective size and effective varience
      !! veff is a namelist variable!

      aeff = a_cl_lay(z)
      veff = var_cl_lay(z)

      const = nd_cl_lay(z) / (gamma((1.0_dp - 2.0_dp*veff)/veff) * &
        & (aeff*veff)**((2.0_dp*veff - 1.0_dp)/veff))

      do m = 1, ndist
        ! Distribution in cm-3 cm-1
        nd_dist(m) =  const * a_dist(m)**((1.0_dp - 3.0_dp*veff)/veff) * &
          & exp(-(a_dist(m)/(aeff*veff)))

        if ((ieee_is_nan(nd_dist(m)) .eqv. .True.) .or. (ieee_is_finite(nd_dist(m)) .eqv. .False.)) then
          nd_dist(m) = 1.0e-99_dp
        end if

        ! Limiter for very low numbers
        nd_dist(m) = max(nd_dist(m),1.0e-99_dp)

        ! Call mie theory routine for this distribution point
        call cl_mie(l,nd_dist(m),a_dist(m),eps_comb,ifunc_k(m),ifunc_a(m),ifunc_g(m))

      end do

    case(8)

      !! Exponential distribution

      lam = 1.0_dp/Ev

      do m = 1, ndist
        ! Distribution in cm-3 cm-1
        nd_dist(m) = nd_cl_lay(z) * (lam * exp(-lam*a_dist(m)))

        if ((ieee_is_nan(nd_dist(m)) .eqv. .True.) .or. (ieee_is_finite(nd_dist(m)) .eqv. .False.)) then
          nd_dist(m) = 1.0e-99_dp
        end if

        ! Limiter for very low numbers
        nd_dist(m) = max(nd_dist(m),1.0e-99_dp)

        ! Call mie theory routine for this distribution point
        call cl_mie(l,nd_dist(m),a_dist(m),eps_comb,ifunc_k(m),ifunc_a(m),ifunc_g(m))

      end do

    case default
      print*, 'ERROR - idist size distribution selection integer not valid - STOPPING'
      print*, 'idist: ', idist_m
      stop
    end select

    ! Check if integration of the function can == nd_cl_lay
    ! cl_out_k = trapz(a_dist(:),nd_dist(:))
    ! print*, cl_out_k, nd_cl_lay(1,z)
    ! stop


    ! Lastly, if we have a distribution, integrate the distribution to find:
    ! 1. Total kappa value, 2. effective SSA, 3. effective g

    if (idist_m == 0) then

    else if (idist_m > 2) then

      select case(idist_int)

      case(1)
        ! Use trapezoid rule - function in [cm-3 cm-1]
        ! Total extinction = integral over all sizes
        cl_out_k = trapz(a_dist(:),ifunc_k(:))
        ! effective SSA = integral for k_sca / k_ext =  integral(k_ext * a) / k_ext
        cl_out_a = trapz(a_dist(:),ifunc_k(:) * ifunc_a(:)) ! Store intermediate result for cl_out_g
        ! effective g = itergral for k_sca * g / k_sca = integeral(k_sca * g) / k_sca
        cl_out_g = trapz(a_dist(:),ifunc_k(:) * ifunc_a(:) * ifunc_g(:))/cl_out_a

        cl_out_a = cl_out_a / cl_out_k ! Albedo is scattering/extinction

      case(2)
        ! Use simpson 1/3 rule - function in [cm-3 um-1]

      case(3)
        ! Use simpson 3/8 rule - function in [cm-3 um-1]

      case(4)

        ! piecewise cubic Hermitian interpolation and integration
        ! cl_out_k = PCHI(a_dist(:), ifunc_k(:), amin, amax)
        ! cl_out_a = PCHI(a_dist(:), ifunc_k(:) * ifunc_a(:), amin, amax) ! Store intermediate result for cl_out_g
        ! cl_out_g = PCHI(a_dist(:), ifunc_k(:) * ifunc_a(:) * ifunc_g(:), amin, amax) /cl_out_a
        !
        ! cl_out_a = cl_out_a / cl_out_k

      case default
        print*, 'ERROR - idist > 2, but idist_int method not valid - STOPPING'
        print*, 'idist, idist_int: ', idist, idist_int
        stop
      end select


    end if


  end subroutine dist_cl

end module cloud_tables_dist
