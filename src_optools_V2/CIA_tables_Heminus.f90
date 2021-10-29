module CIA_tables_Heminus
  use optools_data_mod
  implicit none

  ! Follows Kurucz (1970) parameters - Carbon, Gingerich & Latham (1969) polynomial fit to John (1966)

contains

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

end module CIA_tables_Heminus
