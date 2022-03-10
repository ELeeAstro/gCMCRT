module CIA_tables_Heminus
  use optools_data_mod
  use optools_aux, only : bilinear_interp, locate
  implicit none


contains

  ! Follows Kurucz (1970) parameters - Carbon, Gingerich & Latham (1969) polynomial fit to John (1966)
  subroutine CIA_Heminus(s,l,z,CIA_spec)
    implicit none

    integer, intent(in) :: s, l, z
    integer :: n
    real(kind=dp), intent(out) :: CIA_spec
    real(kind=dp) :: kff
    real(kind=dp) :: T, a, b, c

    ! Temperature [K]
    T = TG_lay(z)

    ! Polynomial coefficents, frequency f in [Hz]
    a = 3.397e-46_dp + (-5.216e-31_dp + 7.039e-15_dp/freq(l))/freq(l)
    b = -4.116e-42_dp + (1.067e-26_dp + 8.135e-11_dp/freq(l))/freq(l)
    c = 5.081e-37_dp + (-8.724e-23_dp - 5.659e-8_dp/freq(l))/freq(l)
    kff = a*T + b + c/T

    kff = kff * (VMR_lay(CIA_tab(s)%iVMR(1),z) * N_lay(z) * VMR_lay(CIA_tab(s)%iVMR(2),z) * N_lay(z)) !  He [cm-3] * e- [cm-3]

    ! End units [cm-1]

    CIA_spec = kff

  end subroutine CIA_Heminus

  ! Interpolates from Bell table
  subroutine CIA_Heminus_Bell(s,l,z,CIA_spec)
    implicit none

    integer, intent(in) :: s, l, z
    real(kind=dp), intent(out) :: CIA_spec

    integer :: iwl, iwl1, iT, iT1
    real(kind=dp) :: T5040
    real(kind=dp) :: xval, x0, x1, yval, y0, y1
    real(kind=dp) :: a00, a10, a01, a11, aval

    T5040 = 5040.0_dp / TG_lay(z)

    call locate(CIA_tab(s)%wl(1,:),wl_A(l),iwl)
    iwl1 = iwl + 1

    if (iwl == 0 .or. iwl1 > CIA_tab(s)%nwl) then
       CIA_spec = 0.0_dp
       return
    end if

    call locate(CIA_tab(s)%T(1,:),T5040,iT)
    iT1 = iT + 1

    if (iT == 0 .or. iT1 > CIA_tab(s)%nT(1)) then
       CIA_spec = 0.0_dp
       return
    end if

    !! wn and T are within the table bounds
    xval = wl_A(l) ; x0 = CIA_tab(s)%wl(1,iwl) ; x1 = CIA_tab(s)%wl(1,iwl1)
    yval = T5040 ; y0 = CIA_tab(s)%T(1,iT) ; y1 = CIA_tab(s)%T(1,iT1)
    a00 = CIA_tab(s)%tab(1,iwl,iT) ; a10 = CIA_tab(s)%tab(1,iwl1,iT)
    a01 = CIA_tab(s)%tab(1,iwl,iT1) ; a11 = CIA_tab(s)%tab(1,iwl1,iT1)

    ! Perform bi-linear interpolation
    call bilinear_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)

    CIA_spec = 1.0e-26_dp * aval * (VMR_lay(CIA_tab(s)%iVMR(1),z) * N_lay(z)) * (VMR_lay(CIA_tab(s)%iVMR(2),z) * N_lay(z) &
    & * kb * TG_lay(z)) !! * P(e-) [dyne cm-2] * He [cm-3]

  end subroutine CIA_Heminus_Bell

end module CIA_tables_Heminus
