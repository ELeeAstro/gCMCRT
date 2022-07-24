module CIA_tables_interp
  use optools_data_mod
  use optools_aux, only : locate, linear_log_interp, bilinear_log_interp, Bezier_interp
  use ieee_arithmetic
  implicit none

  private :: CIA_special
  public :: interp_CIA_tables, interp_CIA_tables_Bezier


contains

  subroutine interp_CIA_tables_Bezier(l,z,CIA_work)
    implicit none

    integer, intent(in) :: l, z
    real(kind=dp), intent(out) :: CIA_work

    integer :: s, sn, j
    integer :: iwn1, iwn2, iwn3, iT1, iT2, iT3
    real(kind=dp) :: T, iCIA
    real(kind=dp), dimension(3) :: Ta, wna, ka, ka_cia


    CIA_work = 0.0_dp

    T = TG_lay(z)

    do s = 1, nCIA

      ! Special species check
      if (CIA_tab(s)%form /= 4) then
        call CIA_special(s,l,z,CIA_work)
        cycle
      end if

      sn = 0
      if (CIA_tab(s)%nset > 1) then
        do j = 1, CIA_tab(s)%nset
          if (wn(l) > CIA_tab(s)%wn_s(j) .and. wn(l) < CIA_tab(s)%wn_e(j)) then
            if (T > CIA_tab(s)%Tmin(j) .and. T < CIA_tab(s)%Tmax(j)) then
              sn = j
              exit
            end if
          end if
        end do

        if (sn == 0) then

          do j = 1, CIA_tab(s)%nset
            if (wn(l) > CIA_tab(s)%wn_s(j) .and. wn(l) < CIA_tab(s)%wn_e(j)) then
              if (T < CIA_tab(s)%Tmin(j)) then
                sn = j
                exit
              end if
              if (T > CIA_tab(s)%Tmax(j)) then
                sn = j
                exit
              end if
            end if
          end do

          if (sn == 0) then
            cycle
          end if

        end if

      else
       sn = 1
      end if

      ! Locate required wn indexes in CIA wn array
      call locate(CIA_tab(s)%wn(sn,1:CIA_tab(s)%irec(sn)),wn(l),iwn2)

      ! Check in wavenumber within bounds
      if ((iwn2+1 > CIA_tab(s)%irec(sn)) .or. (iwn2 < 1)) then
        cycle
      else
        iwn1 = iwn2 - 1
        iwn3 = iwn2 + 1
        if (iwn1 <= 0) then
          iwn1 = 1
          iwn2 = 2
          iwn3 = 3
        else if (iwn3 > CIA_tab(s)%irec(sn)) then
          iwn1 = CIA_tab(s)%irec(sn) - 2
          iwn2 = CIA_tab(s)%irec(sn) - 1
          iwn3 = CIA_tab(s)%irec(sn)
        end if
      end if

      wna(1) = CIA_tab(s)%wn(sn,iwn1)
      wna(2) = CIA_tab(s)%wn(sn,iwn2)
      wna(3) = CIA_tab(s)%wn(sn,iwn3)

      ! Locate required T indexes in CIA wn array for layer temperature
      call locate(CIA_tab(s)%T(sn,1:CIA_tab(s)%nT(sn)),T,iT2)
      iT1 = iT2 - 1
      iT3 = iT2 + 1

      if (iT1 <= 0) then
        iT1 = 1
        iT2 = 2
        iT3 = 3
      else if (iT3 > CIA_tab(s)%nT(sn)) then
        iT1 = CIA_tab(s)%nT(sn) - 2
        iT2 = CIA_tab(s)%nT(sn) - 1
        iT3 = CIA_tab(s)%nT(sn)
      end if

      Ta(1) = CIA_tab(s)%T(sn,iT1)
      Ta(2) = CIA_tab(s)%T(sn,iT2)
      Ta(3) = CIA_tab(s)%T(sn,iT3)

      if (T <= CIA_tab(s)%Tmin(sn)) then
        ! Interpolate at minimum T for set
        ka(1) = CIA_tab(s)%tab(sn,iwn1,1)
        ka(2) = CIA_tab(s)%tab(sn,iwn2,1)
        ka(3) = CIA_tab(s)%tab(sn,iwn3,1)
        call Bezier_interp(wna(:), ka(:), 3, wn(l), iCIA)
        CIA_work = CIA_work +  iCIA &
          & * VMR_lay(CIA_tab(s)%iVMR(1),z) * N_lay(z) &
          & * VMR_lay(CIA_tab(s)%iVMR(2),z) * N_lay(z)
      else if (T >= CIA_tab(s)%Tmax(sn)) then
        ! Interpolate at maximum T for set
        ka(1) = CIA_tab(s)%tab(sn,iwn1,CIA_tab(s)%nT(1))
        ka(2) = CIA_tab(s)%tab(sn,iwn2,CIA_tab(s)%nT(1))
        ka(3) = CIA_tab(s)%tab(sn,iwn3,CIA_tab(s)%nT(1))
        call Bezier_interp(wna(:), ka(:), 3, wn(l), iCIA)
        CIA_work = CIA_work +  iCIA &
          & * VMR_lay(CIA_tab(s)%iVMR(1),z) * N_lay(z) &
          & * VMR_lay(CIA_tab(s)%iVMR(2),z) * N_lay(z)
      else
        ! Both wavenumber and temperature are within table bounds, perform Bezier interpolation 4 times
        ka(1) = CIA_tab(s)%tab(sn,iwn1,iT1)
        ka(2) = CIA_tab(s)%tab(sn,iwn2,iT1)
        ka(3) = CIA_tab(s)%tab(sn,iwn3,iT1)
        call Bezier_interp(wna(:), ka(:), 3, wn(l), ka_cia(1)) ! Result at T1, wn_in
        ka(1) = CIA_tab(s)%tab(sn,iwn1,iT2)
        ka(2) = CIA_tab(s)%tab(sn,iwn2,iT2)
        ka(3) = CIA_tab(s)%tab(sn,iwn3,iT2)
        call Bezier_interp(wna(:), ka(:), 3, wn(l), ka_cia(2)) ! Result at T2, wn_in
        ka(1) = CIA_tab(s)%tab(sn,iwn1,iT3)
        ka(2) = CIA_tab(s)%tab(sn,iwn2,iT3)
        ka(3) = CIA_tab(s)%tab(sn,iwn3,iT3)
        call Bezier_interp(wna(:), ka(:), 3, wn(l), ka_cia(3)) ! Result at T3, wn_in
        call Bezier_interp(Ta(:), ka_cia(:), 3, T, iCIA) ! Result at T_in, wn_in
        CIA_work = CIA_work +  iCIA &
          & * VMR_lay(CIA_tab(s)%iVMR(1),z) * N_lay(z) &
          & * VMR_lay(CIA_tab(s)%iVMR(2),z) * N_lay(z)
      end if
      
    end do

  end subroutine interp_CIA_tables_Bezier

  subroutine interp_CIA_tables(l,z,CIA_work)
    implicit none

    integer, intent(in) :: l, z
    real(kind=dp), intent(out) :: CIA_work
    integer :: s, sn, j
    integer :: iwn, iwn1, iT, iT1
    real(kind=dp) :: xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval


    CIA_work = 0.0_dp
    do s = 1, nCIA

      ! Special species check
      if (CIA_tab(s)%form /= 4) then
        call CIA_special(s,l,z,CIA_work)
        cycle
      end if

      sn = 0
      if (CIA_tab(s)%nset > 1) then
        do j = 1, CIA_tab(s)%nset
          if (wn(l) > CIA_tab(s)%wn_s(j) .and. wn(l) < CIA_tab(s)%wn_e(j)) then
            if (TG_lay(z) > CIA_tab(s)%Tmin(j) .and. TG_lay(z) < CIA_tab(s)%Tmax(j)) then
              sn = j
              exit
            end if
          end if
        end do

        if (sn == 0) then

          do j = 1, CIA_tab(s)%nset
            if (wn(l) > CIA_tab(s)%wn_s(j) .and. wn(l) < CIA_tab(s)%wn_e(j)) then
              if (TG_lay(z) < CIA_tab(s)%Tmin(j)) then
                sn = j
                exit
              end if
              if (TG_lay(z) > CIA_tab(s)%Tmax(j)) then
                sn = j
                exit
              end if
            end if
          end do

          if (sn == 0) then
            cycle
          end if

        end if

      else
       sn = 1
      end if

      ! Locate required wn indexes in CIA wn array
      call locate(CIA_tab(s)%wn(sn,1:CIA_tab(s)%irec(sn)),wn(l),iwn)
      iwn1 = iwn + 1
      ! Check in wavenumber within bounds
      if ((iwn1 > CIA_tab(s)%irec(sn)) .or. (iwn < 1)) then
        cycle
      end if
      ! Locate required T indexes in CIA wn array for layer temperature
      call locate(CIA_tab(s)%T(sn,1:CIA_tab(s)%nT(sn)),TG_lay(z),iT)
      iT1 = iT + 1

      !! Perform temperature edge case check
      if (iT < 1) then

        ! Temperature of layer is outside lower bounds of table
        !print*, 'CIA: TG_lay < minval(T) @: ', CIA_tab(s)%sp, z, TG_lay(z), minval(CIA_tab(s)%T(:)), 'Assuming = minval(T)'

        ! Perform wn linear interp to minval(T)
        xval = wn(l) ; x0 = CIA_tab(s)%wn(sn,iwn) ; x1 = CIA_tab(s)%wn(sn,iwn1)
        y0 = CIA_tab(s)%tab(sn,iwn,1) ; y1 = CIA_tab(s)%tab(sn,iwn1,1)

        ! Perform log linear interpolation
        call linear_log_interp(xval, x0, x1, y0, y1, yval)

        ! Check for NaN's from interpolation
        if (ieee_is_nan(yval) .eqv. .True.) then
          print*, 'CIA: NaN in CIA table linear_log_interp: ', l, z, CIA_tab(s)%sp
          print*, '--', xval, yval, x0, x1, y0, y1
        end if

        ! Add to result to work variable in units of [cm-1]
        CIA_work = CIA_work + yval &
          & * VMR_lay(CIA_tab(s)%iVMR(1),z) * N_lay(z) &
          & * VMR_lay(CIA_tab(s)%iVMR(2),z) * N_lay(z)

      else if (iT1 > CIA_tab(s)%nT(sn)) then

        ! Temperature of layer is outside upper bounds of table
        !print*, 'CIA: TG_lay > maxval(T) @: ', CIA_tab(s)%sp, z, TG_lay(z), maxval(CIA_tab(s)%T(:)), 'Assuming = maxval(T)'

        ! Perform wn linear interp to maxval(T)
        xval = wn(l) ; x0 = CIA_tab(s)%wn(sn,iwn) ; x1 = CIA_tab(s)%wn(sn,iwn1)
        y0 = CIA_tab(s)%tab(sn,iwn,CIA_tab(s)%nT(1)) ; y1 = CIA_tab(s)%tab(sn,iwn1,CIA_tab(s)%nT(1))

        ! Perform log linear interpolation
        call linear_log_interp(xval, x0, x1, y0, y1, yval)

        ! Check for NaN's from interpolation
        if (ieee_is_nan(yval) .eqv. .True.) then
          print*, 'CIA: NaN in CIA table linear_log_interp: ', l, z, CIA_tab(s)%sp
          print*, '--', xval, yval, x0, x1, y0, y1
        end if

        ! Add to result to work variable in units of [cm-1]
        CIA_work = CIA_work + yval &
          & * VMR_lay(CIA_tab(s)%iVMR(1),z) * N_lay(z) &
          & * VMR_lay(CIA_tab(s)%iVMR(2),z) * N_lay(z)

      else

        !! wn and T are within the table bounds
        xval = wn(l) ; x0 = CIA_tab(s)%wn(sn,iwn) ; x1 = CIA_tab(s)%wn(sn,iwn1)
        yval = TG_lay(z) ; y0 = CIA_tab(s)%T(sn,iT) ; y1 = CIA_tab(s)%T(sn,iT1)
        a00 = CIA_tab(s)%tab(sn,iwn,iT) ; a10 = CIA_tab(s)%tab(sn,iwn1,iT)
        a01 = CIA_tab(s)%tab(sn,iwn,iT1) ; a11 = CIA_tab(s)%tab(sn,iwn1,iT1)

        ! Perform bi-linear interpolation
        call bilinear_log_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)

        ! Check for NaN's from bi-linear interpolation
        if (ieee_is_nan(aval) .eqv. .True.) then
          print*, 'CIA: NaN in CIA table bilinear_log_interp: ', l, z, CIA_tab(s)%sp
          print*, '--', xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval
        end if

        ! Add to result to work variable in units of [cm-1]
        CIA_work = CIA_work + aval &
          & * VMR_lay(CIA_tab(s)%iVMR(1),z) * N_lay(z) &
          & * VMR_lay(CIA_tab(s)%iVMR(2),z) * N_lay(z)

      end if

    end do

  end subroutine interp_CIA_tables

  subroutine CIA_special(s,l,z,CIA_work)
    use CIA_tables_Hminus
    use CIA_tables_Heminus
    use CIA_tables_H2minus
    use CIA_tables_fake_H2O_special
    implicit none

    integer, intent(in) :: s, l, z
    real(kind=dp), intent(inout) :: CIA_work
    real(kind=dp) :: CIA_spec


    select case (CIA_tab(s)%sp)

    case ('H-')

      call CIA_Hminus(s,l,z,CIA_spec)

    case ('He-')

      if (CIA_tab(s)%form == 2) then
        call CIA_Heminus_Bell(s,l,z,CIA_spec)
      else
        call CIA_Heminus(s,l,z,CIA_spec)
      end if

    case ('H2-')

      call CIA_H2minus_Bell(s,l,z,CIA_spec)
 
    case ('H2O')

      call Fake_H2O_special(s,l,z,CIA_spec)

    case default

      print*, 'ERROR - CIA species not found in CIA_special - STOPPING'
      print*, 'Species: ', CIA_tab(s)%sp, CIA_tab(s)%form
      stop

    end select

    ! Check for NaN's from interpolation
    if (ieee_is_nan(CIA_spec) .eqv. .True.) then
      print*, 'CIA: NaN in CIA table special: ', s, l, z, CIA_tab(s)%sp
      print*, '--', CIA_spec
    end if

    CIA_work = CIA_work + CIA_spec

  end subroutine CIA_special


end module CIA_tables_interp
