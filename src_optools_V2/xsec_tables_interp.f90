module xsec_tables_interp
  use optools_data_mod
  use optools_aux, only : locate_triplet, locate_inside, Bezier_interp
  use ieee_arithmetic
  implicit none

  private
  public :: interp_xsec_tables_Bezier


contains

  subroutine interp_xsec_tables_Bezier(l,z,xsec_work)
    implicit none

    integer, intent(in) :: l, z
    real(dp), dimension(nxsec), intent(out) :: xsec_work

    integer :: iwl_idx(3), wl_region, iwl_exact, s
    real(dp) :: wll
    real(dp), dimension(3) :: lxa, wla

   ! Interpolate each species xsec table to the wavelength
    xsec_work(:) = 0.0_dp

    wll = wl(l)

    do s = 1, nxsec

      call locate_triplet(xsec_tab(s)%wl(:),wll,iwl_idx,wl_region,iwl_exact)

      ! Cross-sections are zero outside their tabulated wavelength range.
      if (wl_region /= locate_inside) cycle

      if (iwl_exact > 0) then
        xsec_work(s) = 10.0_dp**xsec_tab(s)%lx_abs(iwl_exact)
        cycle
      end if

      ! Wavelength is within table range, perform Bezier interpolation
      wla(:) = xsec_tab(s)%wl(iwl_idx(:))
      lxa(:) = xsec_tab(s)%lx_abs(iwl_idx(:))
      call Bezier_interp(wla,lxa,wll,xsec_work(s))
      xsec_work(s) = 10.0_dp**xsec_work(s)

      !print*, s, wll, xsec_tab(s)%wl(iwl2),xsec_work(s)
      !print*,  lxa(:), wla(:)

    end do

  end subroutine interp_xsec_tables_Bezier  


end module xsec_tables_interp
