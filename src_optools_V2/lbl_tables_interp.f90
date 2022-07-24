module lbl_tables_interp
  use optools_data_mod
  use optools_aux, only : locate, linear_log_interp, bilinear_log_interp, Bezier_interp
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

      !! Do large corner cases situation checks - this sucks but has to be done !!
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
          ! Prssure is within table, perform bilinear interpolation
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
          ! Prssure is within table, perform bilinear interpolation
          z0 = lbl_tab(s)%P(iP) ; z1 = lbl_tab(s)%P(iP1)
          a00 = lbl_tab(s)%k_abs(iwl,iP,1) ; a01 = lbl_tab(s)%k_abs(iwl,iP1,1)
          a10 = lbl_tab(s)%k_abs(iwl1,iP,1) ; a11 = lbl_tab(s)%k_abs(iwl1,iP1,1)
          call bilinear_log_interp(xval, zval, x0, x1, z0, z1, a00, a10, a01, a11, aval)
          lbl_work(s) = aval
        end if
      end if

      !! Check after T edge case
      if (edge_case .eqv. .True.) then
        ! Check for NaN's from interpolation
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
        ! Check for NaN's from interpolation
        if (ieee_is_nan(lbl_work(s)) .eqv. .True.) then
          print*, 'lbl: NaN in lbl table pressure edge case: ', l, z, s, lbl_tab(s)%sp
          print*, '---', iwl, iwl1, iT, iT1, iP, iP1, '---'
        end if
        cycle
      end if

      !! Finally, after all the checks the point is within the table range
      !! Perform a tri-linear interpolation

      x0 = lbl_tab(s)%wl(iwl) ; x1 = lbl_tab(s)%wl(iwl1)
      y0 = lbl_tab(s)%T(iT) ; y1 = lbl_tab(s)%T(iT1)
      z0 = lbl_tab(s)%P(iP) ; z1 = lbl_tab(s)%P(iP1)

      ! Point is within table, perform tri-linear interpolation
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

      ! Check for NaN's from interpolation
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
      call locate(lbl_tab(s)%T(:),T,iT2)
      iT1 = iT2 - 1
      iT3 = iT2 + 1

      if (iT1 <= 0) then
        iT1 = 1
        iT2 = 2
        iT3 = 3
      else if (iT3 > lbl_tab(s)%nT) then
        iT1 = lbl_tab(s)%nT - 2
        iT2 = lbl_tab(s)%nT - 1
        iT3 = lbl_tab(s)%nT
      end if

      lTa(1) = lbl_tab(s)%lT(iT1)
      lTa(2) = lbl_tab(s)%lT(iT2)
      lTa(3) = lbl_tab(s)%lT(iT3)

      ! Find upper and lower T and P triplet indexes
      call locate(lbl_tab(s)%P(:),P,iP2)
      iP1 = iP2 - 1
      iP3 = iP2 + 1 

      if (iP1 <= 0) then
        iP1 = 1
        iP2 = 2
        iP3 = 3
      else if (iP3 > lbl_tab(s)%nP) then
        iP1 = lbl_tab(s)%nP - 2
        iP2 = lbl_tab(s)%nP - 1
        iP3 = lbl_tab(s)%nP
      end if

      lPa(1) = lbl_tab(s)%lP(iP1)
      lPa(2) = lbl_tab(s)%lP(iP2)
      lPa(3) = lbl_tab(s)%lP(iP3)

      if (T >= lbl_tab(s)%T(lbl_tab(s)%nT)) then
        ! Temperature is too high, outside table range - use highest available T data
        if (P >= lbl_tab(s)%P(lbl_tab(s)%nP)) then
           ! Pressure is too high, outside table range - use highest available P data
           lbl_work(s) = 10.0_dp**lbl_tab(s)%lk_abs(l,lbl_tab(s)%nP,lbl_tab(s)%nT)
        else if (P <= lbl_tab(s)%P(1)) then
           ! Pressure is too low, outside table range - use lowest available P data
           lbl_work(s)  = 10.0_dp**lbl_tab(s)%lk_abs(l,1,lbl_tab(s)%nT)
        else
          ! Pressure is within table range, perform Bezier interpolation at highest T
          lka(1) = lbl_tab(s)%lk_abs(l,iP1,lbl_tab(s)%nT)
          lka(2) = lbl_tab(s)%lk_abs(l,iP2,lbl_tab(s)%nT)
          lka(3) = lbl_tab(s)%lk_abs(l,iP3,lbl_tab(s)%nT)
          call Bezier_interp(lPa(:), lka(:), 3, lP, lbl_work(s))
          lbl_work(s) = 10.0_dp**lbl_work(s)
        end if
      else if (T <= lbl_tab(s)%T(1)) then
        ! Temperature is too low, outside table range - use lowest available T data
        if (P >= lbl_tab(s)%P(lbl_tab(s)%nP)) then
           ! Pressure is too high, outside table range - use highest available P data
            lbl_work(s) = 10.0_dp**lbl_tab(s)%lk_abs(l,lbl_tab(s)%nP,1)
        else if (P <= lbl_tab(s)%P(1)) then
           ! Pressure is too low, outside table range - use lowest available P data
           lbl_work(s) = 10.0_dp**lbl_tab(s)%lk_abs(l,1,1)
        else
          ! Pressure is within table range, perform linear interpolation at lowest T
          lka(1) = lbl_tab(s)%lk_abs(l,iP1,1)
          lka(2) = lbl_tab(s)%lk_abs(l,iP2,1)
          lka(3) = lbl_tab(s)%lk_abs(l,iP3,1)
          call Bezier_interp(lPa(:), lka(:), 3, lP, lbl_work(s))
          lbl_work(s) = 10.0_dp**lbl_work(s)
        end if
      else
        ! Temperature is within the normal range
        if (P >= lbl_tab(s)%P(lbl_tab(s)%nP)) then
          ! Pressure is too high, outside table range - use highest availible P
          lka(1) = lbl_tab(s)%lk_abs(l,lbl_tab(s)%nP,iT1)
          lka(2) = lbl_tab(s)%lk_abs(l,lbl_tab(s)%nP,iT2)
          lka(3) = lbl_tab(s)%lk_abs(l,lbl_tab(s)%nP,iT3)
          call Bezier_interp(lTa(:), lka(:), 3, lT, lbl_work(s))
          lbl_work(s) = 10.0_dp**lbl_work(s)
        else if (P <= lbl_tab(s)%P(1)) then
          ! Pressure is too low, outside table range - use lowest availible P
          lka(1) = lbl_tab(s)%lk_abs(l,1,iT1)
          lka(2) = lbl_tab(s)%lk_abs(l,1,iT2)
          lka(3) = lbl_tab(s)%lk_abs(l,1,IT3)
          call Bezier_interp(lTa(:), lka(:), 3, lT, lbl_work(s))
          lbl_work(s) = 10.0_dp**lbl_work(s)
        else
          ! Both pressure and temperature are within table bounds, perform Bezier interpolation 4 times
          lka(1) = lbl_tab(s)%lk_abs(l,iP1,iT1)
          lka(2) = lbl_tab(s)%lk_abs(l,iP2,iT1)
          lka(3) = lbl_tab(s)%lk_abs(l,iP3,iT1)
          call Bezier_interp(lPa(:), lka(:), 3, lP, lka_lbl(1)) ! Result at T1, P_in
          lka(1) = lbl_tab(s)%lk_abs(l,iP1,iT2)
          lka(2) = lbl_tab(s)%lk_abs(l,iP2,iT2)
          lka(3) = lbl_tab(s)%lk_abs(l,iP3,iT3)
          call Bezier_interp(lPa(:), lka(:), 3, lP, lka_lbl(2)) ! Result at T2, P_in
          lka(1) = lbl_tab(s)%lk_abs(l,iP1,iT3)
          lka(2) = lbl_tab(s)%lk_abs(l,iP2,iT3)
          lka(3) = lbl_tab(s)%lk_abs(l,iP3,iT3)
          call Bezier_interp(lPa(:), lka(:), 3, lP, lka_lbl(3)) ! Result at T3, P_in
          call Bezier_interp(lTa(:), lka_lbl(:), 3, lT, lbl_work(s)) ! Result at T_in, P_in
          lbl_work(s) = 10.0_dp**lbl_work(s)
        end if
      end if

    end do


  end subroutine interp_lbl_tables_Bezier

end module lbl_tables_interp
