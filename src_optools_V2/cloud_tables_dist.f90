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
    real(kind=dp), dimension(ndist) :: z_tr, r_dist, wext, wsca

    integer :: m, idist_m
    real(kind=dp), dimension(3) :: cl_out_k2, cl_out_a2, cl_out_g2
    real(kind=dp) :: beta, alpha, aeff, const, Ev, Var, muu, sigg, lam, veff_l
    real(kind=dp) :: z_val, ext_int, sca_int, g_int


    Ev = a_cl_lay(z)
    Var = var_cl_lay(z)

    idist_m = idist

    select case(idist_m)
    case(1)

      ! Single particle size (idist = 1, nmode = 1)
      call cl_mie(z,l,nd_cl_lay(z),a_cl_lay(z),eps_comb,cl_out_k,cl_out_a,cl_out_g)

      !print*, l, nd_cl_lay(z),a_cl_lay(z),eps_comb,cl_out_k,cl_out_a,cl_out_g

    case(2)

      ! Narrow three-size distribution with 1% radius spacing to smooth Mie-resonance features.
      call cl_mie(z,l,nd_cl_lay(z),a_cl_lay(z)*0.99_dp,eps_comb,cl_out_k2(1),cl_out_a2(1),cl_out_g2(1))
      call cl_mie(z,l,nd_cl_lay(z),a_cl_lay(z),eps_comb,cl_out_k2(2),cl_out_a2(2),cl_out_g2(2))
      call cl_mie(z,l,nd_cl_lay(z),a_cl_lay(z)*1.01_dp,eps_comb,cl_out_k2(3),cl_out_a2(3),cl_out_g2(3))

      ext_int = sum(cl_out_k2(:))
      sca_int = sum(cl_out_k2(:) * cl_out_a2(:))
      g_int = sum(cl_out_k2(:) * cl_out_a2(:) * cl_out_g2(:))
      call set_cloud_moments(z,l,eps_comb,ext_int/3.0_dp,ext_int,sca_int,g_int, &
        & cl_out_k,cl_out_a,cl_out_g)

    case(3)

      !! log-normal distribution - transformed trapezoid rule
      !! Substitution: r(z) = exp(ln(Ev) + sqrt(2)*ln(sig_g)*z), z in [-5, 5]
      !! The log-normal integral transforms to:
      !!   integral = (N/sqrt(pi)) * trapz( exp(-z^2) * C_ext(r(z)), z )
      !! This removes the 1/r singularity and sig_g prefactor, giving a
      !! well-conditioned Gaussian-weighted integral over uniform z-space.
      !! cl_mie is called with nd=1 to return pure cross-sections [cm2];
      !! nd_dist holds the Gaussian shape weight exp(-z^2), and N/sqrt(pi)
      !! is applied at the integration step. SSA and g are pure ratios so
      !! the N/sqrt(pi) prefactor cancels there.
      !! Ev = median/geometric mean radius; Var = geometric standard deviation (sig_g)

      do m = 1, ndist
        ! Uniformly spaced z from -5 to 5
        z_tr(m) = -5.0_dp + 10.0_dp * real(m-1, kind=dp) / real(ndist-1, kind=dp)
        ! Transformed radius
        r_dist(m) = exp(log(Ev) + sqrt(2.0_dp) * log(Var) * z_tr(m))
        ! Gaussian shape weight (N/sqrt(pi) applied at integration step below)
        nd_dist(m) = exp(-z_tr(m)**2)

        ! Call Mie routine with unit number density to obtain pure cross-sections
        call cl_mie(z,l,1.0_dp,r_dist(m),eps_comb,ifunc_k(m),ifunc_a(m),ifunc_g(m))

      end do

      ! Pre-compute weighted extinction and scattering integrands once
      wext(:) = ifunc_k(:) * nd_dist(:)
      wsca(:) = wext(:) * ifunc_a(:)
      ext_int  = trapz(z_tr(:), wext(:))
      sca_int  = trapz(z_tr(:), wsca(:))
      g_int = trapz(z_tr(:), wsca(:) * ifunc_g(:))

      ! Integrate in z-space; N/sqrt(pi) only enters cl_out_k, cancels in SSA and g
      call set_cloud_moments(z,l,eps_comb,(nd_cl_lay(z)/sqrt(pi))*ext_int, &
        & ext_int,sca_int,g_int,cl_out_k,cl_out_a,cl_out_g)

    case(4)

      !! Gamma distribution: the particle size in the PRF sets the distribution parameters.
      !! Use eff_fac to set the variance relative to the mean particle size, typically 0 < eff_fac << 1 (~0.1).

      alpha = Ev**2/Var
      beta = Ev/Var

      const = nd_cl_lay(z) * exp(alpha*log(beta) - log_gamma(alpha))

      do m = 1, ndist
        ! Distribution in cm-3 cm-1
        nd_dist(m) = const * a_dist(m)**(alpha-1.0_dp) * exp(-beta*a_dist(m))

        if ((ieee_is_nan(nd_dist(m)) .eqv. .True.) .or. (ieee_is_finite(nd_dist(m)) .eqv. .False.)) then
          nd_dist(m) = 1.0e-99_dp
        end if

        ! Limiter for very low numbers
        nd_dist(m) = max(nd_dist(m),1.0e-99_dp)

        ! Call mie theory routine for this distribution point
        call cl_mie(z,l,nd_dist(m),a_dist(m),eps_comb,ifunc_k(m),ifunc_a(m),ifunc_g(m))

      end do

    case(5)

      !! Inverse-gamma distribution: the particle size in the PRF sets the distribution parameters.
      !! Use eff_fac to set the variance relative to the mean particle size, typically 0 < eff_fac << 1 (~0.1).

      alpha = Ev**2/Var + 2.0_dp
      beta = Ev*(alpha - 1.0_dp)

      const = nd_cl_lay(z) * exp(alpha*log(beta) - log_gamma(alpha))

      do m = 1, ndist
        ! Distribution in cm-3 cm-1
        nd_dist(m) = const * (1.0_dp/a_dist(m))**(alpha+1.0_dp) * exp(-beta/a_dist(m))

        if ((ieee_is_nan(nd_dist(m)) .eqv. .True.) .or. (ieee_is_finite(nd_dist(m)) .eqv. .False.)) then
          nd_dist(m) = 1.0e-99_dp
        end if

        ! Limiter for very low numbers
        nd_dist(m) = max(nd_dist(m),1.0e-99_dp)

        ! Call mie theory routine for this distribution point
        call cl_mie(z,l,nd_dist(m),a_dist(m),eps_comb,ifunc_k(m),ifunc_a(m),ifunc_g(m))

      end do

    case(6)

      !! Rayleigh distribution: the particle size in the PRF sets the distribution sigma.

      !! sigg is directly related to the distribution mean or variance.
      sigg = Ev/sqrt(pi/2.0_dp)
      !sigg = sqrt(Var/(2.0_dp - pi/2.0_dp))

      do m = 1, ndist
        ! Distribution in cm-3 cm-1
        nd_dist(m) = nd_cl_lay(z) * (a_dist(m)/sigg**2) * exp(-a_dist(m)**2/(2.0_dp * sigg**2))

        if ((ieee_is_nan(nd_dist(m)) .eqv. .True.) .or. (ieee_is_finite(nd_dist(m)) .eqv. .False.)) then
          nd_dist(m) = 1.0e-99_dp
        end if

        ! Limiter for very low numbers
        nd_dist(m) = max(nd_dist(m),1.0e-99_dp)

        ! Call mie theory routine for this distribution point
        call cl_mie(z,l,nd_dist(m),a_dist(m),eps_comb,ifunc_k(m),ifunc_a(m),ifunc_g(m))

      end do

    case(7)

      !! Hansen distribution: the particle size in the PRF sets the effective distribution size.
      !! According to Hansen - 0 < veff < 0.5 typically

      !! Variables related to the effective size and effective variance.
      aeff = a_cl_lay(z)
      veff_l = var_cl_lay(z)

      lam = (1.0_dp - 2.0_dp*veff_l)/veff_l
      const = nd_cl_lay(z) * exp(-log_gamma(lam) - lam*log(aeff*veff_l))

      do m = 1, ndist
        ! Distribution in cm-3 cm-1
        nd_dist(m) =  const * a_dist(m)**((1.0_dp - 3.0_dp*veff_l)/veff_l) * &
          & exp(-(a_dist(m)/(aeff*veff_l)))

        if ((ieee_is_nan(nd_dist(m)) .eqv. .True.) .or. (ieee_is_finite(nd_dist(m)) .eqv. .False.)) then
          nd_dist(m) = 1.0e-99_dp
        end if

        ! Limiter for very low numbers
        nd_dist(m) = max(nd_dist(m),1.0e-99_dp)

        ! Call mie theory routine for this distribution point
        call cl_mie(z,l,nd_dist(m),a_dist(m),eps_comb,ifunc_k(m),ifunc_a(m),ifunc_g(m))

      end do

    case(8)

      !! Exponential distribution - log-z transformed trapezoid rule
      !! Substitution: u = log(z/z_center), z = r/Ev, z_center = 2.0
      !! weight = exp(-z)*z (Jacobian from log-z sub)
      !! cl_mie called with nd=1 to return pure cross-sections [cm2];
      !! nd_dist holds the shape weight exp(-z)*z, and N is applied at
      !! the integration step. SSA and g are pure ratios so N cancels there.
      !! Matches Python: int_exp_ana_trapz_logz_psd.py

      do m = 1, ndist
        ! Uniform grid in u = log(z/z_center), z_min=0.05, z_max=15.0, z_center=2.0
        z_tr(m) = log(0.05_dp/2.0_dp) + (log(15.0_dp/2.0_dp) - log(0.05_dp/2.0_dp)) &
                  * real(m-1, kind=dp) / real(ndist-1, kind=dp)
        ! z = z_center * exp(u),  r = Ev * z
        z_val = 2.0_dp * exp(z_tr(m))
        r_dist(m) = Ev * z_val
        ! Shape weight exp(-z)*z (N applied at integration step below)
        nd_dist(m) = exp(-z_val) * z_val

        ! Call Mie routine with unit number density to obtain pure cross-sections
        call cl_mie(z,l,1.0_dp,r_dist(m),eps_comb,ifunc_k(m),ifunc_a(m),ifunc_g(m))

      end do

      ! Pre-compute weighted extinction and scattering integrands once
      wext(:) = ifunc_k(:) * nd_dist(:)
      wsca(:) = wext(:) * ifunc_a(:)
      ext_int  = trapz(z_tr(:), wext(:))
      sca_int  = trapz(z_tr(:), wsca(:))
      g_int = trapz(z_tr(:), wsca(:) * ifunc_g(:))

      ! Integrate over u; N only enters cl_out_k, cancels in SSA and g
      call set_cloud_moments(z,l,eps_comb,nd_cl_lay(z)*ext_int, &
        & ext_int,sca_int,g_int,cl_out_k,cl_out_a,cl_out_g)

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

    else if (idist_m > 2 .and. idist_m /= 3 .and. idist_m /= 8) then

      select case(idist_int)

      case(1)
        ! Use trapezoid rule - function in [cm-3 cm-1]
        ! Pre-compute scattering integrand once to avoid duplicate array products
        wsca(:) = ifunc_k(:) * ifunc_a(:)
        ext_int = trapz(a_dist(:), ifunc_k(:))
        sca_int  = trapz(a_dist(:), wsca(:))
        g_int = trapz(a_dist(:), wsca(:) * ifunc_g(:))
        call set_cloud_moments(z,l,eps_comb,ext_int,ext_int,sca_int,g_int, &
          & cl_out_k,cl_out_a,cl_out_g)

      case(2:4)
        print*, 'ERROR - Cloud integration method is not implemented - STOPPING'
        print*, 'idist, idist_int: ', idist, idist_int
        stop

      case default
        print*, 'ERROR - idist > 2, but idist_int method not valid - STOPPING'
        print*, 'idist, idist_int: ', idist, idist_int
        stop
      end select


    end if

    call validate_cloud_output(z,l,eps_comb,cl_out_k,cl_out_a,cl_out_g)

  end subroutine dist_cl

  subroutine set_cloud_moments(z,l,eps_comb,opacity,ext_moment,sca_moment, &
    & g_moment,cl_out_k,cl_out_a,cl_out_g)
    implicit none

    integer, intent(in) :: z, l
    complex(kind=dp), intent(in) :: eps_comb
    real(kind=dp), intent(in) :: opacity, ext_moment, sca_moment, g_moment
    real(kind=dp), intent(out) :: cl_out_k, cl_out_a, cl_out_g

    if (.not. ieee_is_finite(opacity) .or. .not. ieee_is_finite(ext_moment) .or. &
      & .not. ieee_is_finite(sca_moment) .or. .not. ieee_is_finite(g_moment) .or. &
      & opacity < 0.0_dp .or. ext_moment < 0.0_dp .or. sca_moment < 0.0_dp) then
      print*, 'ERROR - Corrupted cloud distribution moments - STOPPING'
      print*, 'Layer, wavelength index, wavelength [um]: ', z, l, wl(l)
      print*, 'Layer number density, radius [cm], variance: ', &
        & nd_cl_lay(z), a_cl_lay(z), var_cl_lay(z)
      print*, 'Refractive index, idist, idist_int, imie: ', &
        & eps_comb, idist, idist_int, imie
      print*, 'Opacity, extinction, scattering, scattering-g moments: ', &
        & opacity, ext_moment, sca_moment, g_moment
      stop
    end if

    cl_out_k = max(opacity,1.0e-199_dp)
    if (ext_moment <= 0.0_dp .or. sca_moment <= 0.0_dp) then
      cl_out_a = 0.0_dp
      cl_out_g = 0.0_dp
    else
      cl_out_a = min(1.0_dp,max(0.0_dp,sca_moment/ext_moment))
      cl_out_g = min(1.0_dp,max(-1.0_dp,g_moment/sca_moment))
    end if

  end subroutine set_cloud_moments

  subroutine validate_cloud_output(z,l,eps_comb,cl_out_k,cl_out_a,cl_out_g)
    implicit none

    integer, intent(in) :: z, l
    complex(kind=dp), intent(in) :: eps_comb
    real(kind=dp), intent(inout) :: cl_out_k, cl_out_a, cl_out_g

    if (.not. ieee_is_finite(cl_out_k) .or. .not. ieee_is_finite(cl_out_a) .or. &
      & .not. ieee_is_finite(cl_out_g) .or. cl_out_k < 0.0_dp) then
      print*, 'ERROR - Corrupted cloud distribution output - STOPPING'
      print*, 'Layer, wavelength index, wavelength [um]: ', z, l, wl(l)
      print*, 'Layer number density, radius [cm], variance: ', &
        & nd_cl_lay(z), a_cl_lay(z), var_cl_lay(z)
      print*, 'Refractive index, idist, idist_int, imie: ', &
        & eps_comb, idist, idist_int, imie
      print*, 'Opacity, albedo, asymmetry: ', cl_out_k, cl_out_a, cl_out_g
      stop
    end if

    cl_out_k = max(cl_out_k,1.0e-199_dp)
    cl_out_a = min(1.0_dp,max(0.0_dp,cl_out_a))
    if (cl_out_a <= 0.0_dp) then
      cl_out_g = 0.0_dp
    else
      cl_out_g = min(1.0_dp,max(-1.0_dp,cl_out_g))
    end if

  end subroutine validate_cloud_output

end module cloud_tables_dist
