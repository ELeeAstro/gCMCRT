module cloud_tables_mie
  use optools_data_mod
  use mie_routines, only : shexqnn2
  use ieee_arithmetic
  implicit none


  private
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

    !! Variables for q_DHS
    real(kind=dp) :: a_mant

    !! Parameters for DHS

    xsec = pi * (a*1.0e-4_dp)**2

    x = (twopi * a)/wl(l)
    if (x <= 1.00001e-6_dp) then
      x = 1.00001e-6_dp
    else if (x >= 100000.0_dp) then
      x = 100000.0_dp
    end if

    select case(imie)

    case(1)

      !! Use the MIEX Mie theory routine
      ! - careful with memory in parallel for large x
      call shexqnn2(eps, x, rQext, rQsca, rQabs, rQbk, rQpr, ralbedo, rg, &
        & rier, rSA1, rSA2, rdoSA, rnang)

        cl_out_k = rQext * xsec * nd
        cl_out_a = ralbedo
        cl_out_g = rg
     
      if ((ieee_is_nan(cl_out_k) .eqv. .True.) .or. (rier /= 0)) then
        print*, 'cl: NaN in cl mie theory: ', imie, l, wl(l), rQext, xsec, nd, a, x, eps
      end if

    case(2)

       a_mant = a / (1.0_dp - fmax)**(third)

      !! Use DHS (Distribution of hollow spheres) routine
      call q_dhs(real(real(eps),kind=dp), real(aimag(eps),kind=dp), wl(l), a_mant, xsec_ext, xsec_sca, rg, fmax)

      cl_out_k = xsec_ext * nd
      cl_out_a = xsec_sca/xsec_ext
      cl_out_g = rg

    case default

    end select

  end subroutine cl_mie


end module cloud_tables_mie
