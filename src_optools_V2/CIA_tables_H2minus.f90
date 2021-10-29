module CIA_tables_H2minus
  use optools_data_mod
  implicit none

  ! Bell (1980) paramaters
  real(kind=dp), dimension(8), parameter :: &
    & theta = (/0.5_dp, 0.8_dp, 1.0_dp, 1.2_dp, 1.6_dp, 2.0_dp, 2.8_dp, 3.6_dp /) &

  real(kind=dp), dimension(18), parameter :: &
    & wlA = (/3505_dp, 4142_dp, 5063_dp, 5696_dp, 6509_dp, 7594_dp, 9113_dp, 11391_dp, 15188_dp, 18226_dp, 22783_dp, 30377_dp, &
    &  36452_dp, 45565_dp, 60753_dp, 91130_dp, 113913_dp, 151883_dp/)

  real(kind=dp), dimension(8,18), parameter :: &
    & kff = (//)

  real(kind=dp), dimension(8,18), parameter :: &
    & kff12 = (//)

contains

  subroutine CIA_H2minus(s,l,z,CIA_spec)
    implicit none

    integer, intent(in) :: s, l, z
    integer :: n
    real(kind=dp), intent(out) :: CIA_spec
    real(kind=dp) :: kff, kbf, fbf, xbf, sff
    real(kind=dp) :: T, T5040

    T = TG_lay(z)
    T5040 = 5040.0_dp / T

    ! Do bound free calculation
    if (wl(l) > lam_0) then
      xbf = 0.0_dp
      !kbf = 0.0_dp
    else
      fbf = 0.0_dp
      do n = 1, 6
        fbf = fbf + Cn_bf(n) * (1.0_dp/wl(l) - 1.0_dp/lam_0)**((real(n,kind=dp)-1.0_dp)/2.0_dp)
      end do
      xbf = 1.0e-18_dp * wl(l)**3 * (1.0_dp/wl(l) - 1.0_dp/lam_0)**(3.0_dp/2.0_dp) * fbf
      !kbf = 0.750_dp * T**(-5.0_dp/2.0_dp) * exp(alf/(lam_0*T)) * (1.0_dp - exp(-alf/(wl(l)*T))) * xbf
    end if

    ! Do free free calculation
    sff = 0.0_dp
    if (wl(l) >= 0.3645_dp) then
      do n = 1, 6
        sff = sff + T5040**((real(n,kind=dp)+1.0_dp)/2.0_dp) &
                & * (wl(l)**2*An_ff2(n) + Bn_ff2(n)  + Cn_ff2(n)/wl(l) + Dn_ff2(n)/wl(l)**2 + En_ff2(n)/wl(l)**3 &
                & + Fn_ff2(n)/wl(l)**4)
      end do
      kff = 1.0e-29_dp * sff
    else if ((wl(l) < 0.3645_dp) .and. (wl(l) > 0.1823_dp)) then
      do n = 1, 6
        sff = sff + T5040**((real(n,kind=dp)+1.0_dp)/2.0_dp) &
                & * (wl(l)**2*An_ff1(n) + Bn_ff1(n)  + Cn_ff1(n)/wl(l) + Dn_ff1(n)/wl(l)**2 + En_ff1(n)/wl(l)**3 &
                & + Fn_ff1(n)/wl(l)**4)
      end do
      kff = 1.0e-29_dp * sff
    else
      kff = 0.0_dp
    end if

    ! xbf is in [cm2 molecule-1] and kff is in [cm4 dyne-1] - convert to [cm-1] before CMCRT output
    kbf = xbf * VMR_lay(CIA_tab(s)%iVMR_3(1),z) * N_lay(z) ! * H- [molecule cm-3]
    !kbf = kbf * (VMR_lay(CIA_tab(s)%iVMR_3(2),z) * N_lay(z) * VMR_lay(CIA_tab(s)%iVMR_3(3),z) * N_lay(z)) * kb * T ! * P(e-) [dyne cm-2] * H [cm-3]
    kff = kff * (VMR_lay(CIA_tab(s)%iVMR_3(2),z) * N_lay(z) * VMR_lay(CIA_tab(s)%iVMR_3(3),z) * N_lay(z)) * kb * T ! * P(e-) [dyne cm-2] * H [cm-3]

    CIA_spec = kbf + kff

  end subroutine CIA_H2minus

end module CIA_tables_H2minus
