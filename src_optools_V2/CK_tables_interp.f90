module CK_tables_interp
  use optools_data_mod
  use optools_aux, only : locate, locate_triplet, locate_below, locate_above, &
    & bilinear_log_interp, linear_log_interp, Bezier_interp
  use ieee_arithmetic
  implicit none


contains

  subroutine interp_CK_tables_Bezier(l,z,nG,CK_work)
    implicit none

    integer, intent(in) :: l, z, nG

    real(kind=dp), dimension(nCK,nG), intent(out) :: CK_work

    integer :: s, g
    integer :: iT1, iT2, iT3, iP1, iP2, iP3
    integer :: iT_idx(3), iP_idx(3), T_region, P_region
    integer :: iT_exact, iP_exact, iT_fixed, iP_fixed
    real(kind=dp) :: T, P, lT, lP
    real(kind=dp), dimension(3) :: lTa, lPa, lka, lka_ck

    ! Interpolate from each CK table to the layer
    CK_work(:,:) = 0.0_dp

    !! A large if block statement is used to determine if the T-p point
    !! is within the ck table, then interpolates.
    !! If T > Tmax, then T = Tmax, if P > Pmax, P = Pmax and if P < Pmin, P = Pmin
    !! This avoids failures when k tables do not cover the intended T-P range.
    !! K TABLES SHOULD BE CREATED OVER AN APPROPRIATE T-P RANGE.

    T = TG_lay(z)
    P = PG_lay(z)
    lT = log10(T)
    lP = log10(P)

    do s = 1, nCK

      ! Find temperature grid index triplet
      call locate_triplet(CK_tab(s)%T(:),T,iT_idx,T_region,iT_exact)
      iT1 = iT_idx(1)
      iT2 = iT_idx(2)
      iT3 = iT_idx(3)

      lTa(1) = CK_tab(s)%lT(iT1)
      lTa(2) = CK_tab(s)%lT(iT2)
      lTa(3) = CK_tab(s)%lT(iT3)

      ! Find pressure grid index triplet
      call locate_triplet(CK_tab(s)%P(:),P,iP_idx,P_region,iP_exact)
      iP1 = iP_idx(1)
      iP2 = iP_idx(2)
      iP3 = iP_idx(3)

      lPa(1) = CK_tab(s)%lP(iP1)
      lPa(2) = CK_tab(s)%lP(iP2)
      lPa(3) = CK_tab(s)%lP(iP3)

      iT_fixed = iT_exact
      if (T_region == locate_below) iT_fixed = 1
      if (T_region == locate_above) iT_fixed = CK_tab(s)%nT

      iP_fixed = iP_exact
      if (P_region == locate_below) iP_fixed = 1
      if (P_region == locate_above) iP_fixed = CK_tab(s)%nP

      if (iT_fixed > 0 .and. iP_fixed > 0) then
        CK_work(s,:) = 10.0_dp**CK_tab(s)%lk_abs(l,iP_fixed,iT_fixed,:)
      else if (iT_fixed > 0) then
        do g = 1, CK_tab(s)%nG
          lka(:) = CK_tab(s)%lk_abs(l,iP_idx(:),iT_fixed,g)
          call Bezier_interp(lPa,lka,lP,CK_work(s,g))
          CK_work(s,g) = 10.0_dp**CK_work(s,g)
        end do
      else if (iP_fixed > 0) then
        do g = 1, CK_tab(s)%nG
          lka(:) = CK_tab(s)%lk_abs(l,iP_fixed,iT_idx(:),g)
          call Bezier_interp(lTa,lka,lT,CK_work(s,g))
          CK_work(s,g) = 10.0_dp**CK_work(s,g)
        end do
      else
        do g = 1, CK_tab(s)%nG
          lka(:) = CK_tab(s)%lk_abs(l,iP_idx(:),iT1,g)
          call Bezier_interp(lPa,lka,lP,lka_ck(1))
          lka(:) = CK_tab(s)%lk_abs(l,iP_idx(:),iT2,g)
          call Bezier_interp(lPa,lka,lP,lka_ck(2))
          lka(:) = CK_tab(s)%lk_abs(l,iP_idx(:),iT3,g)
          call Bezier_interp(lPa,lka,lP,lka_ck(3))
          call Bezier_interp(lTa,lka_ck,lT,CK_work(s,g))
          CK_work(s,g) = 10.0_dp**CK_work(s,g)
        end do
      end if

    end do

  end subroutine interp_CK_tables_Bezier

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

      !! Handle the pressure and temperature boundary cases explicitly.
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
          ! Pressure is within the table; perform linear interpolation.
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
          ! Pressure is within the table; perform linear interpolation.
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
        ! Check for NaNs from interpolation.
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
          a0 = CK_tab(s)%k_abs(iwl,iP,iT,g) ; a1 = CK_tab(s)%k_abs(iwl,iP,iT1,g)
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
        ! Check for NaNs from interpolation.
        do g = 1, CK_tab(s)%nG
          if (ieee_is_nan(CK_work(s,g)) .eqv. .True.) then
            print*, 'CK: NaN in CK table pressure edge case: ', l, z, g, s, CK_tab(s)%sp
            print*, '---', iwl, iwl1, iT, iT1, iP, iP1, '---'
          end if
        end do
        cycle
      end if

      !! Finally, after all the checks the point is within the table range
      !! Perform bilinear interpolation.

      y0 = CK_tab(s)%T(iT) ; y1 = CK_tab(s)%T(iT1)
      z0 = CK_tab(s)%P(iP) ; z1 = CK_tab(s)%P(iP1)
      ! Interpolate each g-ordinate
      do g = 1, CK_tab(s)%nG
        ! The point is within the table; perform bilinear interpolation.
        a00 = CK_tab(s)%k_abs(iwl,iP,iT,g) ; a10 = CK_tab(s)%k_abs(iwl,iP1,iT,g)
        a01 = CK_tab(s)%k_abs(iwl,iP,iT1,g) ; a11 = CK_tab(s)%k_abs(iwl,iP1,iT1,g)

        call bilinear_log_interp(yval, zval, y0, y1, z0, z1, a00, a10, a01, a11, aval)
        CK_work(s,g) = aval

        ! Check for NaNs from interpolation.
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

      !! Handle the pressure and temperature boundary cases explicitly.
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
          ! Pressure is within the table; perform bilinear interpolation.
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
          ! Pressure is within the table; perform bilinear interpolation.
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
        ! Check for NaNs from interpolation.
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
        ! Check for NaNs from interpolation.
        do g = 1, CK_tab(s)%nG
          if (ieee_is_nan(CK_work(s,g)) .eqv. .True.) then
            print*, 'CK: NaN in CK table pressure edge case: ', l, z, g, s, CK_tab(s)%sp
            print*, '---', iwl, iwl1, iT, iT1, iP, iP1, '---'
          end if
        end do
        cycle
      end if

      !! Finally, after all the checks the point is within the table range
      !! Perform trilinear interpolation.

      x0 = CK_tab(s)%wl(iwl) ; x1 = CK_tab(s)%wl(iwl1)
      y0 = CK_tab(s)%T(iT) ; y1 = CK_tab(s)%T(iT1)
      z0 = CK_tab(s)%P(iP) ; z1 = CK_tab(s)%P(iP1)

      ! Interpolate each g-ordinate
      do g = 1, CK_tab(s)%nG
        ! The point is within the table; perform trilinear interpolation.
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

        ! Check for NaNs from interpolation.
        if (ieee_is_nan(CK_work(s,g)) .eqv. .True.) then
          print*, 'CK: NaN in CK table tri-linear_log_interp: ', l, z, g, s, CK_tab(s)%sp
          print*, '---', xval, yval, zval, x0, x1, y0, y1, z0, z1, &
            & a000, a100, a010, a110, a001, a101, a011, a111, aval0, aval1, aval, '---'
        end if

      end do

    end do


  end subroutine interp_CK_tables_wl

end module CK_tables_interp
