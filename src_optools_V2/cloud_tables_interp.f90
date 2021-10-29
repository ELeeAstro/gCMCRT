module cloud_tables_interp
  use optools_data_mod
  use optools_aux, only : locate, linear_log_interp
  use ieee_arithmetic
  implicit none


  private
  public :: interp_cl_tables

contains

  subroutine interp_cl_tables(l,n_work,k_work)
    implicit none

    integer, intent(in) :: l
    integer :: s
    integer :: iwl, iwl1, iwl_ex
    real(kind=dp) :: wl_ex, fac
    real(kind=dp) :: xval, yvaln, yvalk
    real(kind=dp) :: x0, x1, y0n, y1n, y0k, y1k

    real(kind=dp), dimension(ncl), intent(out) :: n_work, k_work
    logical :: extrap_flag

    xval = wl(l)

    n_work(:) = 0.0_dp
    k_work(:) = 0.0_dp

    extrap_flag = .False.

    do s = 1, ncl

      ! Find wavelength grid index
      call locate(cl_tab(s)%wl(:),wl(l),iwl)
      iwl1 = iwl + 1
      ! Check in wavelength within bounds.
      if (iwl1 > cl_tab(s)%nwl) then
        extrap_flag = .True.
        ! wavelength is above bounds of table, extrapolate to value
        ! depending on conducting flag
        select case(cl_tab(s)%conducting)
        case(.False.)
          ! Linear decrease in k proportional to wl
          ! n assumed constant at highest wl value
          n_work(s) = cl_tab(s)%n(nwl)
          k_work(s) = cl_tab(s)%k(nwl) * (cl_tab(s)%wl(nwl)/wl(l))
        case(.True.)
          ! ln-ln extrapolate n and k to higher wavelengths
          wl_ex = 0.7_dp * cl_tab(s)%wl(nwl)
          call locate(cl_tab(s)%wl(:),wl_ex,iwl_ex)
          ! Check that at least 3 points back is used for extrapolation
          if ((iwl_ex > cl_tab(s)%nwl - 3) .or. (iwl_ex == 0)) then
            iwl_ex = cl_tab(s)%nwl - 3
          end if
          fac = log(wl(l)/cl_tab(s)%wl(nwl))/log(cl_tab(s)%wl(iwl_ex)/cl_tab(s)%wl(nwl))
          n_work(s) = exp(log(cl_tab(s)%n(nwl)) &
            & + fac * log(cl_tab(s)%n(iwl_ex)/cl_tab(s)%n(nwl)))
          k_work(s) = exp(log(cl_tab(s)%k(nwl)) &
            & + fac * log(cl_tab(s)%k(iwl_ex)/cl_tab(s)%k(nwl)))
        end select

      else if (iwl < 1) then
        extrap_flag = .True.
        ! wavelength is below bounds of table, extrapolate to value
        ! n and k are equal to the lowest wl value
        n_work(s) = cl_tab(s)%n(1)
        k_work(s) = cl_tab(s)%k(1)

      else
        extrap_flag = .False.
        ! wavelength is within bounds of table, interpolate to values
        x0 = cl_tab(s)%wl(iwl) ; x1 = cl_tab(s)%wl(iwl1)

        y0n = cl_tab(s)%n(iwl) ; y1n = cl_tab(s)%n(iwl1)
        call linear_log_interp(xval, x0, x1, y0n, y1n, yvaln)
        n_work(s) = yvaln

        y0k = cl_tab(s)%k(iwl) ; y1k = cl_tab(s)%k(iwl1)
        call linear_log_interp(xval, x0, x1, y0k, y1k, yvalk)
        k_work(s) = yvalk

        ! Check for NaN's from interpolation
        if ((ieee_is_nan(n_work(s)) .eqv. .True.) .or. (ieee_is_nan(k_work(s)) .eqv. .True.)) then
          print*, 'cl: NaN in cl table: ',  l, s, cl_tab(s)%sp, cl_tab(s)%conducting
          print*, '---', xval, yvaln, yvalk, x0, x1, y0n, y1n, y0k, y1k, extrap_flag, '---'
        end if

      end if

      n_work(s) = max(n_work(s),1e-99_dp)
      k_work(s) = max(k_work(s),1e-99_dp)

      !print*, l, wl(l), s, cl_tab(s)%sp, n_work(s), k_work(s), extrap_flag

    end do

  end subroutine interp_cl_tables

end module cloud_tables_interp
