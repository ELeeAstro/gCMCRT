module CIA_tables_interp
  use optools_data_mod
  use optools_aux, only : locate, linear_log_interp, bilinear_log_interp
  use ieee_arithmetic
  implicit none

  private :: CIA_special
  public :: interp_CIA_tables


contains

  subroutine interp_CIA_tables(l,z,CIA_work)
    implicit none

    integer, intent(in) :: l, z
    real(kind=dp), intent(out) :: CIA_work
    integer :: s
    integer :: iwn, iwn1, iT, iT1
    real(kind=dp) :: xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval


    CIA_work = 0.0_dp
    do s = 1, nCIA

      ! Special species check
      if (CIA_tab(s)%form == 0) then
        call CIA_special(s,l,z,CIA_work)
        cycle
      end if

      ! Locate required wn indexes in CIA wn array
      call locate(CIA_tab(s)%wn(:),wn(l),iwn)
      iwn1 = iwn + 1
      ! Check in wavenumber within bounds
      if ((iwn1 > CIA_tab(s)%nwl) .or. (iwn < 1)) then
        cycle
      end if
      ! Locate required T indexes in CIA wn array for layer temperature
      call locate(CIA_tab(s)%T(:),TG_lay(z),iT)
      iT1 = iT + 1

      !! Perform temperature edge case check
      if (iT < 1) then

        ! Temperature of layer is outside lower bounds of table
        !print*, 'CIA: TG_lay < minval(T) @: ', CIA_tab(s)%sp, z, TG_lay(z), minval(CIA_tab(s)%T(:)), 'Assuming = minval(T)'

        ! Perform wn linear interp to minval(T)
        xval = wn(l) ; x0 = CIA_tab(s)%wn(iwn) ; x1 = CIA_tab(s)%wn(iwn1)
        y0 = CIA_tab(s)%tab(iwn,1) ; y1 = CIA_tab(s)%tab(iwn1,1)

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

      else if (iT1 > CIA_tab(s)%nT) then

        ! Temperature of layer is outside upper bounds of table
        !print*, 'CIA: TG_lay > maxval(T) @: ', CIA_tab(s)%sp, z, TG_lay(z), maxval(CIA_tab(s)%T(:)), 'Assuming = maxval(T)'

        ! Perform wn linear interp to maxval(T)
        xval = wn(l) ; x0 = CIA_tab(s)%wn(iwn) ; x1 = CIA_tab(s)%wn(iwn1)
        y0 = CIA_tab(s)%tab(iwn,CIA_tab(s)%nT) ; y1 = CIA_tab(s)%tab(iwn1,CIA_tab(s)%nT)

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
        xval = wn(l) ; x0 = CIA_tab(s)%wn(iwn) ; x1 = CIA_tab(s)%wn(iwn1)
        yval = TG_lay(z) ; y0 = CIA_tab(s)%T(iT) ; y1 = CIA_tab(s)%T(iT1)
        a00 = CIA_tab(s)%tab(iwn,iT) ; a10 = CIA_tab(s)%tab(iwn1,iT)
        a01 = CIA_tab(s)%tab(iwn,iT1) ; a11 = CIA_tab(s)%tab(iwn1,iT1)

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
    implicit none

    integer, intent(in) :: s, l, z
    real(kind=dp), intent(inout) :: CIA_work
    real(kind=dp) :: CIA_spec


    select case (CIA_tab(s)%sp)

    case ('H-')

      call CIA_Hminus(s,l,z,CIA_spec)

    case ('He-')

      call CIA_Heminus(s,l,z,CIA_spec)

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
