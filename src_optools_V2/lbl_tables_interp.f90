module lbl_tables_interp
  use optools_data_mod
  use optools_aux, only : locate, locate_triplet, locate_below, locate_above, &
    & linear_log_interp, bilinear_log_interp, Bezier_interp
  use ieee_arithmetic
  implicit none


  private
  public :: interp_lbl_tables, interp_lbl_tables_Bezier

contains

  subroutine interp_lbl_tables(l,z,lbl_work)
    implicit none

    integer, intent(in) :: l, z
    integer :: s
    integer :: iwl, iwl1, iT, iT1, iP, iP1
    real(kind=dp) :: T, P, wll
    real(kind=dp) :: aval
    real(kind=dp) :: xval, yval, zval, x0, x1, y0, y1, z0, z1
    real(kind=dp) :: a000, a010, a100, a110, a001, a011, a101, a111, aval0, aval1
    real(kind=dp) :: a00, a01, a10, a11
    real(kind=dp) :: a0, a1
    logical :: edge_case

    real(kind=dp), dimension(nlbl), intent(out) :: lbl_work

    wll = wl(l)
    T = TG_lay(z)
    P = PG_lay(z)

    xval = wll
    yval = T
    zval = P

    ! Interpolate each species lbl table to the layer properties
    lbl_work(:) = 0.0_dp

    do s = 1, nlbl

      ! Change edge_case flag
      edge_case = .False.

      ! Find wavelength grid index
      call locate(lbl_tab(s)%wl(:),wll,iwl)
      iwl1 = iwl + 1
      ! Check in wavelength within bounds
      if ((iwl1 > lbl_tab(s)%nwl) .or. (iwl < 1)) then
        cycle
      end if
      ! Find temperature grid index
      call locate(lbl_tab(s)%T(:),T,iT)
      iT1 = iT + 1
      ! Find pressure grid index
      call locate(lbl_tab(s)%P(:),P,iP)
      iP1 = iP + 1

      !! Handle the pressure and temperature boundary cases explicitly.
      !! Check temperature corner cases first (most likely)
      if (iT == lbl_tab(s)%nT) then
        edge_case = .True.
        ! Temperature is too high in layer, calculate at highest available T
        x0 = lbl_tab(s)%wl(iwl) ; x1 = lbl_tab(s)%wl(iwl1)
        if (iP == lbl_tab(s)%nP) then
          ! Pressure is also too high in layer, calculate at highest P
          a0 = lbl_tab(s)%k_abs(iwl,iP,iT) ; a1 = lbl_tab(s)%k_abs(iwl1,iP,iT)
          call linear_log_interp(xval, x0, x1, a0, a1, aval)
          lbl_work(s) = aval
        else if (iP == 0) then
          ! Pressure is too low in layer, calculate at lowest P
          a0 = lbl_tab(s)%k_abs(iwl,1,iT) ; a1 = lbl_tab(s)%k_abs(iwl1,1,iT)
          call linear_log_interp(xval, x0, x1, a0, a1, aval)
          lbl_work(s) = aval
        else
          ! Pressure is within the table; perform bilinear interpolation.
          z0 = lbl_tab(s)%P(iP) ; z1 = lbl_tab(s)%P(iP1)
          a00 = lbl_tab(s)%k_abs(iwl,iP,iT) ; a01 = lbl_tab(s)%k_abs(iwl,iP1,iT)
          a10 = lbl_tab(s)%k_abs(iwl1,iP,iT) ; a11 = lbl_tab(s)%k_abs(iwl1,iP1,iT)
          call bilinear_log_interp(xval, zval, x0, x1, z0, z1, a00, a10, a01, a11, aval)
          lbl_work(s) = aval
        end if
      else if (iT == 0) then
        edge_case = .True.
        ! Temperature is too low in layer, calculate at lowest available T
        x0 = lbl_tab(s)%wl(iwl) ; x1 = lbl_tab(s)%wl(iwl1)
        if (iP == lbl_tab(s)%nP) then
          ! Pressure is also too high in layer, calculate at highest P
          a0 = lbl_tab(s)%k_abs(iwl,iP,1) ; a1 = lbl_tab(s)%k_abs(iwl1,iP,1)
          call linear_log_interp(xval, x0, x1, a0, a1, aval)
          lbl_work(s) = aval
        else if (iP == 0) then
          ! Pressure is too low in layer, calculate at lowest P
          a0 = lbl_tab(s)%k_abs(iwl,1,1) ; a1 = lbl_tab(s)%k_abs(iwl1,1,1)
          call linear_log_interp(xval, x0, x1, a0, a1, aval)
          lbl_work(s) = aval
        else
          ! Pressure is within the table; perform bilinear interpolation.
          z0 = lbl_tab(s)%P(iP) ; z1 = lbl_tab(s)%P(iP1)
          a00 = lbl_tab(s)%k_abs(iwl,iP,1) ; a01 = lbl_tab(s)%k_abs(iwl,iP1,1)
          a10 = lbl_tab(s)%k_abs(iwl1,iP,1) ; a11 = lbl_tab(s)%k_abs(iwl1,iP1,1)
          call bilinear_log_interp(xval, zval, x0, x1, z0, z1, a00, a10, a01, a11, aval)
          lbl_work(s) = aval
        end if
      end if

      !! Check after T edge case
      if (edge_case .eqv. .True.) then
        ! Check for NaNs from interpolation.
        if (ieee_is_nan(lbl_work(s)) .eqv. .True.) then
          print*, 'lbl: NaN in lbl table temperature edge case: ', l, z, s, lbl_tab(s)%sp
          print*, '---', iwl, iwl1, iT, iT1, iP, iP1, '---'
        end if
        cycle
      end if

      !! Check pressure corner case
      if (iP == lbl_tab(s)%nP) then
        edge_case = .True.
        ! Pressure is too high for table
        ! Temperature is within table (thanks to previous test), perform bilinear interpolation
        x0 = lbl_tab(s)%wl(iwl) ; x1 = lbl_tab(s)%wl(iwl1)
        y0 = lbl_tab(s)%T(iT) ; y1 = lbl_tab(s)%T(iT1)
        a00 = lbl_tab(s)%k_abs(iwl,iP,iT) ; a01 = lbl_tab(s)%k_abs(iwl,iP,iT1)
        a10 = lbl_tab(s)%k_abs(iwl1,iP,iT) ; a11 = lbl_tab(s)%k_abs(iwl1,iP,iT1)
        call bilinear_log_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)
        lbl_work(s) = aval
      else if (iP == 0) then
        edge_case = .True.
        ! Pressure is too low for table
        ! Temperature is within table (thanks to previous test), perform bilinear interpolation
        x0 = lbl_tab(s)%wl(iwl) ; x1 = lbl_tab(s)%wl(iwl1)
        y0 = lbl_tab(s)%T(iT) ; y1 = lbl_tab(s)%T(iT1)
        a00 = lbl_tab(s)%k_abs(iwl,1,iT) ; a01 = lbl_tab(s)%k_abs(iwl,1,iT1)
        a10 = lbl_tab(s)%k_abs(iwl1,1,iT) ; a11 = lbl_tab(s)%k_abs(iwl1,1,iT1)
        call bilinear_log_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)
        lbl_work(s) = aval
      end if

      !! Check after P edge case
      if (edge_case .eqv. .True.) then
        ! Check for NaNs from interpolation.
        if (ieee_is_nan(lbl_work(s)) .eqv. .True.) then
          print*, 'lbl: NaN in lbl table pressure edge case: ', l, z, s, lbl_tab(s)%sp
          print*, '---', iwl, iwl1, iT, iT1, iP, iP1, '---'
        end if
        cycle
      end if

      !! Finally, after all the checks the point is within the table range
      !! Perform trilinear interpolation.

      x0 = lbl_tab(s)%wl(iwl) ; x1 = lbl_tab(s)%wl(iwl1)
      y0 = lbl_tab(s)%T(iT) ; y1 = lbl_tab(s)%T(iT1)
      z0 = lbl_tab(s)%P(iP) ; z1 = lbl_tab(s)%P(iP1)

      ! The point is within the table; perform trilinear interpolation.
      a000 = lbl_tab(s)%k_abs(iwl,iP,iT)
      a100 = lbl_tab(s)%k_abs(iwl1,iP,iT)
      a010 = lbl_tab(s)%k_abs(iwl,iP,iT1)
      a110 = lbl_tab(s)%k_abs(iwl1,iP,iT1)
      a001 = lbl_tab(s)%k_abs(iwl,iP1,iT)
      a101 = lbl_tab(s)%k_abs(iwl1,iP1,iT)
      a011 = lbl_tab(s)%k_abs(iwl,iP1,iT1)
      a111 = lbl_tab(s)%k_abs(iwl1,iP1,iT1)

      call bilinear_log_interp(xval, yval, x0, x1, y0, y1, a000, a100, a010, a110, aval0)
      call bilinear_log_interp(xval, yval, x0, x1, y0, y1, a001, a101, a011, a111, aval1)
      call linear_log_interp(zval, z0, z1, aval0, aval1, aval)

      lbl_work(s) = aval

      ! Check for NaNs from interpolation.
      if (ieee_is_nan(lbl_work(s)) .eqv. .True.) then
        print*, 'lbl: NaN in lbl table tri-linear_log_interp: ', l, z, s, lbl_tab(s)%sp
        print*, '---', xval, yval, zval, x0, x1, y0, y1, z0, z1, &
          & a000, a100, a010, a110, a001, a101, a011, a111, aval0, aval1, aval, '---'
      end if

    end do

  end subroutine interp_lbl_tables

  subroutine interp_lbl_tables_Bezier(l,z,lbl_work)
    implicit none

    integer, intent(in) :: l, z

    real(kind=dp), dimension(nlbl), intent(out) :: lbl_work

    integer :: s
    integer :: iT1, iT2, iT3, iP1, iP2, iP3
    integer :: iT_idx(3), iP_idx(3), T_region, P_region
    integer :: iT_exact, iP_exact, iT_fixed, iP_fixed
    real(kind=dp) :: T, P, lT, lP
    real(kind=dp), dimension(3) :: lTa, lPa, lka, lka_lbl

    ! Interpolate each species lbl table to the layer properties
    lbl_work(:) = 0.0_dp

    T = TG_lay(z)
    P = PG_lay(z)
    lT = log10(TG_lay(z))
    lP = log10(PG_lay(z))

    do s = 1, nlbl

      ! Find temperature grid index triplet
      call locate_triplet(lbl_tab(s)%T(:),T,iT_idx,T_region,iT_exact)
      iT1 = iT_idx(1)
      iT2 = iT_idx(2)
      iT3 = iT_idx(3)

      lTa(1) = lbl_tab(s)%lT(iT1)
      lTa(2) = lbl_tab(s)%lT(iT2)
      lTa(3) = lbl_tab(s)%lT(iT3)

      ! Find the upper and lower T and P triplet indices.
      call locate_triplet(lbl_tab(s)%P(:),P,iP_idx,P_region,iP_exact)
      iP1 = iP_idx(1)
      iP2 = iP_idx(2)
      iP3 = iP_idx(3)

      lPa(1) = lbl_tab(s)%lP(iP1)
      lPa(2) = lbl_tab(s)%lP(iP2)
      lPa(3) = lbl_tab(s)%lP(iP3)

      iT_fixed = iT_exact
      if (T_region == locate_below) iT_fixed = 1
      if (T_region == locate_above) iT_fixed = lbl_tab(s)%nT

      iP_fixed = iP_exact
      if (P_region == locate_below) iP_fixed = 1
      if (P_region == locate_above) iP_fixed = lbl_tab(s)%nP

      if (iT_fixed > 0 .and. iP_fixed > 0) then
        lbl_work(s) = 10.0_dp**lbl_tab(s)%lk_abs(l,iP_fixed,iT_fixed)
      else if (iT_fixed > 0) then
        lka(:) = lbl_tab(s)%lk_abs(l,iP_idx(:),iT_fixed)
        call Bezier_interp(lPa,lka,lP,lbl_work(s))
        lbl_work(s) = 10.0_dp**lbl_work(s)
      else if (iP_fixed > 0) then
        lka(:) = lbl_tab(s)%lk_abs(l,iP_fixed,iT_idx(:))
        call Bezier_interp(lTa,lka,lT,lbl_work(s))
        lbl_work(s) = 10.0_dp**lbl_work(s)
      else
        lka(:) = lbl_tab(s)%lk_abs(l,iP_idx(:),iT1)
        call Bezier_interp(lPa,lka,lP,lka_lbl(1))
        lka(:) = lbl_tab(s)%lk_abs(l,iP_idx(:),iT2)
        call Bezier_interp(lPa,lka,lP,lka_lbl(2))
        lka(:) = lbl_tab(s)%lk_abs(l,iP_idx(:),iT3)
        call Bezier_interp(lPa,lka,lP,lka_lbl(3))
        call Bezier_interp(lTa,lka_lbl,lT,lbl_work(s))
        lbl_work(s) = 10.0_dp**lbl_work(s)
      end if

    end do


  end subroutine interp_lbl_tables_Bezier

end module lbl_tables_interp
