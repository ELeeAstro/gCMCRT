module CK_tables_interp
  use optools_data_mod
  use optools_aux, only : locate, bilinear_log_interp, linear_log_interp
  use ieee_arithmetic
  implicit none


contains

  subroutine interp_CK_tables(l,z,nG,CK_work)
    implicit none

    integer, intent(in) :: l, z, nG
    integer :: s, g
    integer :: iwl, iwl1, iT, iT1, iP, iP1
    real(kind=dp), dimension(nCK,nG), intent(out) :: CK_work
    real(kind=dp) :: T, P, wll
    real(kind=dp) :: aval
    real(kind=dp) :: xval, yval, zval, x0, x1, y0, y1, z0, z1
    real(kind=dp) :: a000, a010, a100, a110, a001, a011, a101, a111, aval0, aval1
    real(kind=dp) :: a00, a01, a10, a11
    real(kind=dp) :: a0, a1
    logical :: edge_case

    T = TG_lay(z)
    P = PG_lay(z)

    yval = T
    zval = P

    ! Interpolate from each CK table to the layer
    CK_work(:,:) = 0.0_dp

    do s = 1, nCK

      ! Change edge_case flag
      edge_case = .False.

      ! Wavelength is known, index is just l
      iwl = l
      ! Find temperature grid index
      call locate(CK_tab(s)%T(:),T,iT)
      iT1 = iT + 1
      ! Find pressure grid index
      call locate(CK_tab(s)%P(:),P,iP)
      iP1 = iP + 1

      !! Do large corner cases situation checks - this sucks but has to be done !!
      !! Check temperature corner cases first (most likely)
      if (iT == CK_tab(s)%nT) then
        edge_case = .True.
        ! Temperature is too high in layer, calculate at highest available T
        if (iP == CK_tab(s)%nP) then
          ! Pressure is also too high in layer, calculate at highest P
          CK_work(s,:) = CK_tab(s)%k_abs(iwl,iP,iT,:)
        else if (iP == 0) then
          ! Pressure is too low in layer, calculate at lowest P
          CK_work(s,:) = CK_tab(s)%k_abs(iwl,1,iT,:)
        else
          ! Prssure is within table, perform linear interpolation
          z0 = CK_tab(s)%P(iP) ; z1 = CK_tab(s)%P(iP1)
          do g = 1, CK_tab(s)%nG
            a0 = CK_tab(s)%k_abs(iwl,iP,iT,g) ; a1 = CK_tab(s)%k_abs(iwl,iP1,iT,g)
            call linear_log_interp(zval, z0, z1, a0, a1, aval)
            CK_work(s,g) = aval
          end do
        end if
      else if (iT == 0) then
        edge_case = .True.
        ! Temperature is too low in layer, calculate at lowest available T
        if (iP == CK_tab(s)%nP) then
          ! Pressure is also too high in layer, calculate at highest P
          CK_work(s,:) = CK_tab(s)%k_abs(iwl,iP,1,:)
        else if (iP == 0) then
          ! Pressure is too low in layer, calculate at lowest P
          CK_work(s,:) = CK_tab(s)%k_abs(iwl,1,1,:)
        else
          ! Prssure is within table, perform linear interpolation
          z0 = CK_tab(s)%P(iP) ; z1 = CK_tab(s)%P(iP1)
          do g = 1, CK_tab(s)%nG
            a0 = CK_tab(s)%k_abs(iwl,iP,1,g) ; a1 = CK_tab(s)%k_abs(iwl,iP1,1,g)
            call linear_log_interp(zval, z0, z1, a0, a1, aval)
            CK_work(s,g) = aval
          end do
        end if
      end if

      !! Check after T edge case
      if (edge_case .eqv. .True.) then
        ! Check for NaN's from interpolation
        do g = 1, CK_tab(s)%nG
          if (ieee_is_nan(CK_work(s,g)) .eqv. .True.) then
            print*, 'CK: NaN in CK table temperature edge case: ', l, z, g, s, CK_tab(s)%sp
            print*, '---', iwl, iwl1, iT, iT1, iP, iP1, '---'
          end if
        end do
        cycle
      end if

      !! Check pressure corner case
      if (iP == CK_tab(s)%nP) then
        edge_case = .True.
        ! Pressure is too high for table
        ! Temperature is within table (thanks to previous test), perform linear interpolation
        y0 = CK_tab(s)%T(iT) ; y1 = CK_tab(s)%T(iT1)
        do g = 1, CK_tab(s)%nG
          a0 = CK_tab(s)%k_abs(iwl1,iP,iT,g) ; a1 = CK_tab(s)%k_abs(iwl1,iP,iT1,g)
          call linear_log_interp(yval, y0, y1, a0, a1, aval)
          CK_work(s,g) = aval
        end do
      else if (iP == 0) then
       edge_case = .True.
        ! Pressure is too low for table
        ! Temperature is within table (thanks to previous test), perform linear interpolation
        y0 = CK_tab(s)%T(iT) ; y1 = CK_tab(s)%T(iT1)
        do g = 1, CK_tab(s)%nG
          a0 = CK_tab(s)%k_abs(iwl,1,iT,g) ; a1 = CK_tab(s)%k_abs(iwl,1,iT1,g)
          call linear_log_interp(yval, y0, y1, a0, a1, aval)
          CK_work(s,g) = aval
        end do
      end if

      !! Check after P edge case
      if (edge_case .eqv. .True.) then
        ! Check for NaN's from interpolation
        do g = 1, CK_tab(s)%nG
          if (ieee_is_nan(CK_work(s,g)) .eqv. .True.) then
            print*, 'CK: NaN in CK table pressure edge case: ', l, z, g, s, CK_tab(s)%sp
            print*, '---', iwl, iwl1, iT, iT1, iP, iP1, '---'
          end if
        end do
        cycle
      end if

      !! Finally, after all the checks the point is within the table range
      !! Perform a bi-linear interpolation

      y0 = CK_tab(s)%T(iT) ; y1 = CK_tab(s)%T(iT1)
      z0 = CK_tab(s)%P(iP) ; z1 = CK_tab(s)%P(iP1)
      ! Interpolate each g-ordinate
      do g = 1, CK_tab(s)%nG
        ! Point is within table, perform bi-linear interpolation
        a00 = CK_tab(s)%k_abs(iwl,iP,iT,g) ; a10 = CK_tab(s)%k_abs(iwl,iP1,iT,g)
        a01 = CK_tab(s)%k_abs(iwl,iP,iT1,g) ; a11 = CK_tab(s)%k_abs(iwl,iP1,iT1,g)

        call bilinear_log_interp(yval, zval, y0, y1, z0, z1, a00, a10, a01, a11, aval)
        CK_work(s,g) = aval

        ! Check for NaN's from interpolation
        if (ieee_is_nan(CK_work(s,g)) .eqv. .True.) then
          print*, 'CK: NaN in CK table bi-linear_log_interp: ', l, z, g, s, CK_tab(s)%sp
          print*, '---', xval, yval, zval,y0, y1, z0, z1, &
          & a00, a10, a01, a11,  aval, '---'
        end if

      end do

    end do

  end subroutine interp_CK_tables

  subroutine interp_CK_tables_wl(l,z,nG,CK_work)
    implicit none

    integer, intent(in) :: l, z, nG
    integer :: s, g
    integer :: iwl, iwl1, iT, iT1, iP, iP1
    real(kind=dp), dimension(nCK,nG), intent(out) :: CK_work
    real(kind=dp) :: T, P, wll
    real(kind=dp) :: aval
    real(kind=dp) :: xval, yval, zval, x0, x1, y0, y1, z0, z1
    real(kind=dp) :: a000, a010, a100, a110, a001, a011, a101, a111, aval0, aval1
    real(kind=dp) :: a00, a01, a10, a11
    real(kind=dp) :: a0, a1
    logical :: edge_case

    wll = wl(l)
    T = TG_lay(z)
    P = PG_lay(z)

    xval = wll
    yval = T
    zval = P

    ! Interpolate from each CK table to the layer
    CK_work(:,:) = 0.0_dp

    do s = 1, nCK

      ! Change edge_case flag
      edge_case = .False.

      ! Find wavelength grid index
      call locate(CK_tab(s)%wl(:),wll,iwl)
      iwl1 = iwl + 1
      ! Check in wavelength within bounds
      if ((iwl1 > CK_tab(s)%nwl) .or. (iwl < 1)) then
        cycle
      end if
      ! Find temperature grid index
      call locate(CK_tab(s)%T(:),T,iT)
      iT1 = iT + 1
      ! Find pressure grid index
      call locate(CK_tab(s)%P(:),P,iP)
      iP1 = iP + 1

      !! Do large corner cases situation checks - this sucks but has to be done !!
      !! Check temperature corner cases first (most likely)
      if (iT == CK_tab(s)%nT) then
        edge_case = .True.
        ! Temperature is too high in layer, calculate at highest available T
        x0 = CK_tab(s)%wl(iwl) ; x1 = CK_tab(s)%wl(iwl1)
        if (iP == CK_tab(s)%nP) then
          ! Pressure is also too high in layer, calculate at highest P
          do g = 1, CK_tab(s)%nG
            a0 = CK_tab(s)%k_abs(iwl,iP,iT,g) ; a1 = CK_tab(s)%k_abs(iwl1,iP,iT,g)
            call linear_log_interp(xval, x0, x1, a0, a1, aval)
            CK_work(s,g) = aval
          end do
        else if (iP == 0) then
          ! Pressure is too low in layer, calculate at lowest P
          do g = 1, CK_tab(s)%nG
            a0 = CK_tab(s)%k_abs(iwl,1,iT,g) ; a1 = CK_tab(s)%k_abs(iwl1,1,iT,g)
            call linear_log_interp(xval, x0, x1, a0, a1, aval)
            CK_work(s,g) = aval
          end do
        else
          ! Prssure is within table, perform bilinear interpolation
          z0 = CK_tab(s)%P(iP) ; z1 = CK_tab(s)%P(iP1)
          do g = 1, CK_tab(s)%nG
            a00 = CK_tab(s)%k_abs(iwl,iP,iT,g) ; a01 = CK_tab(s)%k_abs(iwl,iP1,iT,g)
            a10 = CK_tab(s)%k_abs(iwl1,iP,iT,g) ; a11 = CK_tab(s)%k_abs(iwl1,iP1,iT,g)
            call bilinear_log_interp(xval, zval, x0, x1, z0, z1, a00, a10, a01, a11, aval)
            CK_work(s,g) = aval
          end do
        end if
      else if (iT == 0) then
        edge_case = .True.
        ! Temperature is too low in layer, calculate at lowest available T
        x0 = CK_tab(s)%wl(iwl) ; x1 = CK_tab(s)%wl(iwl1)
        if (iP == CK_tab(s)%nP) then
          ! Pressure is also too high in layer, calculate at highest P
          do g = 1, CK_tab(s)%nG
            a0 = CK_tab(s)%k_abs(iwl,iP,1,g) ; a1 = CK_tab(s)%k_abs(iwl1,iP,1,g)
            call linear_log_interp(xval, x0, x1, a0, a1, aval)
            CK_work(s,g) = aval
          end do
        else if (iP == 0) then
          ! Pressure is too low in layer, calculate at lowest P
          do g = 1, CK_tab(s)%nG
            a0 = CK_tab(s)%k_abs(iwl,1,1,g) ; a1 = CK_tab(s)%k_abs(iwl1,1,1,g)
            call linear_log_interp(xval, x0, x1, a0, a1, aval)
            CK_work(s,g) = aval
          end do
        else
          ! Prssure is within table, perform bilinear interpolation
          z0 = CK_tab(s)%P(iP) ; z1 = CK_tab(s)%P(iP1)
          do g = 1, CK_tab(s)%nG
            a00 = CK_tab(s)%k_abs(iwl,iP,1,g) ; a01 = CK_tab(s)%k_abs(iwl,iP1,1,g)
            a10 = CK_tab(s)%k_abs(iwl1,iP,1,g) ; a11 = CK_tab(s)%k_abs(iwl1,iP1,1,g)
            call bilinear_log_interp(xval, zval, x0, x1, z0, z1, a00, a10, a01, a11, aval)
            CK_work(s,g) = aval
          end do
        end if
      end if

      !! Check after T edge case
      if (edge_case .eqv. .True.) then
        ! Check for NaN's from interpolation
        do g = 1, CK_tab(s)%nG
          if (ieee_is_nan(CK_work(s,g)) .eqv. .True.) then
            print*, 'CK: NaN in CK table temperature edge case: ', l, z, g, s, CK_tab(s)%sp
            print*, '---', iwl, iwl1, iT, iT1, iP, iP1, '---'
          end if
        end do
        cycle
      end if

      !! Check pressure corner case
      if (iP == CK_tab(s)%nP) then
        edge_case = .True.
        ! Pressure is too high for table
        ! Temperature is within table (thanks to previous test), perform bilinear interpolation
        x0 = CK_tab(s)%wl(iwl) ; x1 = CK_tab(s)%wl(iwl1)
        y0 = CK_tab(s)%T(iT) ; y1 = CK_tab(s)%T(iT1)
        do g = 1, CK_tab(s)%nG
          a00 = CK_tab(s)%k_abs(iwl,iP,iT,g) ; a01 = CK_tab(s)%k_abs(iwl,iP,iT1,g)
          a10 = CK_tab(s)%k_abs(iwl1,iP,iT,g) ; a11 = CK_tab(s)%k_abs(iwl1,iP,iT1,g)
          call bilinear_log_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)
          CK_work(s,g) = aval
        end do
      else if (iP == 0) then
        edge_case = .True.
        ! Pressure is too low for table
        ! Temperature is within table (thanks to previous test), perform bilinear interpolation
        x0 = CK_tab(s)%wl(iwl) ; x1 = CK_tab(s)%wl(iwl1)
        y0 = CK_tab(s)%T(iT) ; y1 = CK_tab(s)%T(iT1)
        do g = 1, CK_tab(s)%nG
          a00 = CK_tab(s)%k_abs(iwl,1,iT,g) ; a01 = CK_tab(s)%k_abs(iwl,1,iT1,g)
          a10 = CK_tab(s)%k_abs(iwl1,1,iT,g) ; a11 = CK_tab(s)%k_abs(iwl1,1,iT1,g)
          call bilinear_log_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)
          CK_work(s,g) = aval
        end do
      end if

      !! Check after P edge case
      if (edge_case .eqv. .True.) then
        ! Check for NaN's from interpolation
        do g = 1, CK_tab(s)%nG
          if (ieee_is_nan(CK_work(s,g)) .eqv. .True.) then
            print*, 'CK: NaN in CK table pressure edge case: ', l, z, g, s, CK_tab(s)%sp
            print*, '---', iwl, iwl1, iT, iT1, iP, iP1, '---'
          end if
        end do
        cycle
      end if

      !! Finally, after all the checks the point is within the table range
      !! Perform a tri-linear interpolation

      x0 = CK_tab(s)%wl(iwl) ; x1 = CK_tab(s)%wl(iwl1)
      y0 = CK_tab(s)%T(iT) ; y1 = CK_tab(s)%T(iT1)
      z0 = CK_tab(s)%P(iP) ; z1 = CK_tab(s)%P(iP1)

      ! Interpolate each g-ordinate
      do g = 1, CK_tab(s)%nG
        ! Point is within table, perform tri-linear interpolation
        a000 = CK_tab(s)%k_abs(iwl,iP,iT,g)
        a100 = CK_tab(s)%k_abs(iwl1,iP,iT,g)
        a010 = CK_tab(s)%k_abs(iwl,iP,iT1,g)
        a110 = CK_tab(s)%k_abs(iwl1,iP,iT1,g)
        a001 = CK_tab(s)%k_abs(iwl,iP1,iT,g)
        a101 = CK_tab(s)%k_abs(iwl1,iP1,iT,g)
        a011 = CK_tab(s)%k_abs(iwl,iP1,iT1,g)
        a111 = CK_tab(s)%k_abs(iwl1,iP1,iT1,g)

        call bilinear_log_interp(xval, yval, x0, x1, y0, y1, a000, a100, a010, a110, aval0)
        call bilinear_log_interp(xval, yval, x0, x1, y0, y1, a001, a101, a011, a111, aval1)
        call linear_log_interp(zval, z0, z1, aval0, aval1, aval)

        CK_work(s,g) = aval

        ! Check for NaN's from interpolation
        if (ieee_is_nan(CK_work(s,g)) .eqv. .True.) then
          print*, 'CK: NaN in CK table tri-linear_log_interp: ', l, z, g, s, CK_tab(s)%sp
          print*, '---', xval, yval, zval, x0, x1, y0, y1, z0, z1, &
            & a000, a100, a010, a110, a001, a101, a011, a111, aval0, aval1, aval, '---'
        end if

      end do

    end do


  end subroutine interp_CK_tables_wl

end module CK_tables_interp
