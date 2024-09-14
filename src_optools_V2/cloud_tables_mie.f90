module cloud_tables_mie
  use optools_data_mod
  use mie_routines, only : shexqnn2
  use mieext_mod, only : MieExt
  use bhmie_mod, only : BHMIE
  use bhcoat_mod, only : BHCOAT
  use dhs_mod, only : q_dhs
  use lxmie_mod, only : lxmie
  use mie_approx_mod, only : madt, rayleigh, rayleigh_gans, geo_optics
  use ieee_arithmetic
  implicit none

  public :: cl_mie

contains


  subroutine cl_mie(l,nd,a,eps,cl_out_k,cl_out_a,cl_out_g)
    implicit none

    integer, intent(in) :: l
    real(kind=dp), intent(in) :: nd, a
    complex(kind=dp), intent(in) :: eps
    real(kind=dp), intent(out) :: cl_out_k, cl_out_a, cl_out_g

    real(kind=dp) :: rQext, rQsca, rQabs, rQbk, rQpr, ralbedo, rg
    real(kind=dp) :: x, xsec, xsec_ext, xsec_sca

    !! Parameters for MIEX - set as paramaters to avoid messing up for now
    integer, parameter :: rnang = 1
    integer :: rier
    complex(kind=dp), dimension(rnang) :: rSA1, rSA2
    logical, parameter :: rdoSA = .False.

    !! Parameters for MieExt
    complex(dp) :: eps_in

    !! Parameters for BHMIE
    integer, parameter :: nang = 2
    complex(kind=dp), dimension(2*nang-1) :: S1, S2

    !! Parameters for BHCOAT
    real(dp) :: x_core, x_mant
    complex(dp) :: eps_core, eps_mant

    !! Variables for q_DHS
    real(kind=dp) :: a_mant, e1, e2

    xsec = pi * a**2

    x = (twopi * a)/wl_cm(l)

    if (x < 1.00001e-6_dp) then
      x = 1.00001e-6_dp
    end if

    rier = 0

    select case(imie)

    case(0)
      !! Special limiting cases for efficency
      x = (twopi * a)/wl_cm(l) ! Unlimit size parameter

      if (x < 0.01_dp) then
        !! Use rayleigh scattering approximation
        call rayleigh(x, eps, rQabs, rQsca, rQext, rg)
      else if (x > 100.0_dp) then
        call madt(x, eps, rQabs, rQsca, rQext, rg)
      else
        !! Call LX-MIE with negative k value
        eps_in = cmplx(real(eps,dp),-aimag(eps))
        call lxmie(eps_in, x, rQext, rQsca, rQabs, rg)
      end if

      cl_out_k = rQext * xsec * nd
      cl_out_a = rQsca/rQext
      cl_out_g = rg

      if ((ieee_is_nan(cl_out_a) .eqv. .True.) .or. (rier /= 0)) then
        print*, 'Special: NaN: ', l, real(wl(l)), rQext, a, x, eps
      end if

    case(1)

      !! Use the MIEX Mie theory routine
      ! - careful with memory in parallel for large x
      call shexqnn2(eps, x, rQext, rQsca, rQabs, rQbk, rQpr, ralbedo, rg, &
        & rier, rSA1, rSA2, rdoSA, rnang)

      cl_out_k = rQext * xsec * nd
      cl_out_a = ralbedo
      cl_out_g = rg
     
      if ((ieee_is_nan(cl_out_a) .eqv. .True.) .or. (rier /= 0)) then
        if (x > 100.0_dp) then
           print*, 'MieX: lx approx: ', l, real(wl(l)), rQext, a, x, eps
           cl_out_k = 2.0_dp * xsec * nd ; cl_out_a = 0.9_dp ; cl_out_g = 0.9_dp
        else
          print*, 'MieX: zero approx: ', l, real(wl(l)), rQext, a, x, eps
          cl_out_k = 0.0_dp ; cl_out_a = 0.0_dp ; cl_out_g = 0.0_dp
        end if
      end if

    case(2)

      !! Use the MieExt routine
      eps_in = cmplx(real(eps,dp),-aimag(eps))

      call MieExt(x, eps_in, rQext, rQsca, rg, rier)

      cl_out_k = rQext * xsec * nd
      cl_out_a = rQsca/rQext
      cl_out_g = rg

      if ((ieee_is_nan(cl_out_k) .eqv. .True.) .or. (rier /= 0)) then
        if (x > 100.0_dp) then
           print*, 'MieExt: lx approx: ', l, real(wl(l)), rQext, a, x, eps_in
           cl_out_k = 2.0_dp * xsec * nd ; cl_out_a = 0.9_dp ; cl_out_g = 0.9_dp
        else
          print*, 'MieExt: zero approx: ', l, real(wl(l)), rQext, a, x, eps_in
          cl_out_k = 0.0_dp ; cl_out_a = 0.0_dp ; cl_out_g = 0.0_dp
        end if
      end if

    case(3)
      
      !! Use the classic BHMIE routine (w. Draine edits)
      call BHMIE(x, eps, nang, S1, S2, rQext, rQsca, rQbk, rg, rier)

      cl_out_k = rQext * xsec * nd
      cl_out_a = rQsca/rQext
      cl_out_g = rg

      if ((ieee_is_nan(cl_out_a) .eqv. .True.) .or. (rier /= 0)) then
        if (x > 100.0_dp) then
           print*, 'BHMIE: lx approx: ', l, real(wl(l)), rQext, a, x, eps
           cl_out_k = 2.0_dp * xsec * nd ; cl_out_a = 0.9_dp ; cl_out_g = 0.9_dp
        else
          print*, 'BHMIE: zero approx: ', l, real(wl(l)), rQext, a, x, eps
          cl_out_k = 0.0_dp ; cl_out_a = 0.0_dp ; cl_out_g = 0.0_dp
        end if
      end if

    case(4)

      a_mant = (a*1e4_dp) / (1.0_dp - fmax)**(third) ! Convert to um for q_dhs

      e1 = real(eps,dp)
      e2 = aimag(eps)

      !! Use DHS (Distribution of hollow spheres) routine - note: output is cross section in um^2
      !! So Q is misleading, should be C
      call q_dhs(e1, e2, wl(l), a_mant, rQext, rQsca, rg, fmax, rier)

      cl_out_k = rQext*1e-8_dp * nd
      cl_out_a = rQsca/rQext
      cl_out_g = rg
      !print*, cl_out_k, cl_out_a, cl_out_g

      if ((ieee_is_nan(cl_out_a) .eqv. .True.) .or. (rier /= 0)) then
        if (x > 100.0_dp) then
           print*, 'DHS: lx approx: ', l, real(wl(l)), rQext, a, x, eps
           cl_out_k = 2.0_dp * xsec * nd ; cl_out_a = 0.9_dp ; cl_out_g = 0.9_dp
        else
          print*, 'DHS: zero approx: ', l, real(wl(l)), rQext, a, x, eps
          cl_out_k = 0.0_dp ; cl_out_a = 0.0_dp ; cl_out_g = 0.0_dp
        end if
      end if

    case(5)

      !! Set core and mantle values
      x_core = 0.1_dp * x
      x_mant = x

      eps_core = eps
      eps_mant = eps

      !! Use the classic BHCOAT (coated spheres) routine (w. Draine edits)
      call BHCOAT(x_core, x_mant, eps_core, eps_mant, rQext, rQsca, rQbk, rg)

      cl_out_k = rQext * xsec * nd
      cl_out_a = rQsca/rQext
      cl_out_g = rg

      if ((ieee_is_nan(cl_out_a) .eqv. .True.) .or. (rier /= 0)) then
        if (x > 100.0_dp) then
           print*, 'BHCOAT: lx approx: ', l, real(wl(l)), rQext, a, x, eps
           cl_out_k = 2.0_dp * xsec * nd ; cl_out_a = 0.9_dp ; cl_out_g = 0.9_dp
        else
          print*, 'BHCOAT: zero approx: ', l, real(wl(l)), rQext, a, x, eps
          cl_out_k = 0.0_dp ; cl_out_a = 0.0_dp ; cl_out_g = 0.0_dp
        end if
      end if

    case(6)

      !! use the LX-MIE routine translated from Kitzmann & Heng (2018)
      eps_in = cmplx(real(eps,dp),-aimag(eps))
      call lxmie(eps_in, x, rQext, rQsca, rQabs, rg)

      cl_out_k = rQext * xsec * nd
      cl_out_a = rQsca/rQext
      cl_out_g = rg

      if ((ieee_is_nan(cl_out_a) .eqv. .True.) .or. (rier /= 0)) then
        if (x > 100.0_dp) then
           print*, 'LX-MIE: lx approx: ', l, real(wl(l)), rQext, rQsca, a, x, eps
           cl_out_k = 2.0_dp * xsec * nd ; cl_out_a = 0.9_dp ; cl_out_g = 0.9_dp
        else
          print*, 'LX-MIE: zero approx: ', l, real(wl(l)), rQext, rQsca, a, x, eps
          cl_out_k = 0.0_dp ; cl_out_a = 0.0_dp ; cl_out_g = 0.0_dp
        end if
      end if     
     
    case default
      print*, 'Invalid Mie theory method selected: ',imie
      stop
    end select

  end subroutine cl_mie

end module cloud_tables_mie
