module xsec_tables_interp
  use optools_data_mod
  use optools_aux, only : locate, Bezier_interp
  use ieee_arithmetic
  implicit none

  private
  public :: interp_xsec_tables_Bezier


contains

  subroutine interp_xsec_tables_Bezier(l,z,xsec_work)
    implicit none

    integer, intent(in) :: l, z
    real(dp), dimension(nxsec), intent(out) :: xsec_work

    integer :: iwl1, iwl2, iwl3, s
    real(dp) :: wll
    real(dp), dimension(3) :: lxa, wla

   ! Interpolate each species xsec table to the wavelength
    xsec_work(:) = 0.0_dp

    wll = wl(l)

    do s = 1, nxsec

      ! Check in wavelength within bounds
      if ((wll > xsec_tab(s)%wl(xsec_tab(s)%nwl)) .or. (wll < xsec_tab(s)%wl(1))) then
        cycle
      end if
      
      ! Find wavelength grid index
      call locate(xsec_tab(s)%wl(:),wll,iwl2)
      iwl1 = iwl2 - 1
      iwl3 = iwl2 + 1

      if (iwl1 <= 0) then
        iwl1 = 1
        iwl2 = 2
        iwl3 = 3
      else if (iwl3 > xsec_tab(s)%nwl) then
        iwl1 = xsec_tab(s)%nwl - 2
        iwl2 = xsec_tab(s)%nwl - 1
        iwl3 = xsec_tab(s)%nwl
      end if

      ! Wavelength is within table range, perform Bezier interpolation
      wla(1) = xsec_tab(s)%wl(iwl1)
      wla(2) = xsec_tab(s)%wl(iwl2)
      wla(3) = xsec_tab(s)%wl(iwl3)

      lxa(1) = xsec_tab(s)%lx_abs(iwl1)
      lxa(2) = xsec_tab(s)%lx_abs(iwl2)
      lxa(3) = xsec_tab(s)%lx_abs(iwl3)
      call Bezier_interp(wla(:), lxa(:), 3, wll, xsec_work(s))
      xsec_work(s) = 10.0_dp**xsec_work(s)

      !print*, s, wll, xsec_tab(s)%wl(iwl2),xsec_work(s)
      !print*,  lxa(:), wla(:)

    end do

  end subroutine interp_xsec_tables_Bezier  


end module xsec_tables_interp
