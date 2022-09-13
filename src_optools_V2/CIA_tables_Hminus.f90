module CIA_tables_Hminus
  use optools_data_mod
  use ieee_arithmetic
  implicit none

  ! John (1988) paramaters
  real(kind=dp), dimension(6), parameter :: &
    & An_ff1 = (/518.1021_dp, 472.2636_dp, -482.2089_dp, 115.5291_dp, 0.0_dp, 0.0_dp/), &
    & Bn_ff1 = (/-734.8666_dp, 1443.4137_dp, -737.1616_dp, 169.6374_dp, 0.0_dp, 0.0_dp/), &
    & Cn_ff1 = (/1021.1775_dp, -1977.3395_dp, 1096.8827_dp, -245.6490_dp, 0.0_dp, 0.0_dp/), &
    & Dn_ff1 = (/-479.0721_dp, 922.3575_dp, -521.1341_dp, 114.2430_dp, 0.0_dp, 0.0_dp/), &
    & En_ff1 = (/93.1373_dp, -178.9275_dp, 101.7963_dp, -21.9972_dp, 0.0_dp, 0.0_dp/), &
    & Fn_ff1 = (/-6.4285_dp, 12.3600_dp, -7.0571_dp, 1.5097_dp, 0.0_dp, 0.0_dp/)

  real(kind=dp), dimension(6), parameter :: &
    & An_ff2 = (/0.0_dp, 2483.3460_dp, -3449.8890_dp, 2200.0400_dp, -696.2710_dp, 88.2830_dp/), &
    & Bn_ff2 = (/0.0_dp, 285.8270_dp, -1158.3820_dp, 2427.7190_dp, -1841.4000_dp, 444.5170_dp/), &
    & Cn_ff2 = (/0.0_dp, -2054.2910_dp, 8746.5230_dp, -13651.1050_dp, 8642.9700_dp, -1863.8640_dp/), &
    & Dn_ff2 = (/0.0_dp, 2827.7760_dp, -11485.6320_dp, 16755.5240_dp, -10051.5300_dp, 2095.2880_dp/), &
    & En_ff2 = (/0.0_dp, -1341.5370_dp, 5303.6090_dp, -7510.4940_dp, 4400.0670_dp, -901.7880_dp/), &
    & Fn_ff2 = (/0.0_dp, 208.9520_dp, -812.9390_dp, 1132.7380_dp, -655.0200_dp, 132.9850_dp/)

  real(kind=dp), dimension(6), parameter :: &
    & Cn_bf = (/152.519_dp, 49.534_dp, -118.858_dp, 92.536_dp, -34.194_dp, 4.982_dp /)

  real(kind=dp), parameter :: alf = 1.439e8_dp, lam_0 = 1.6419_dp, lam_min = 0.125_dp

contains

  subroutine CIA_Hminus(s,l,z,CIA_spec)
    implicit none

    integer, intent(in) :: s, l, z
    integer :: n
    real(kind=dp), intent(out) :: CIA_spec
    real(kind=dp) :: kff, kbf, fbf, xbf, sff
    real(kind=dp) :: T, T5040

    T = TG_lay(z)
    T5040 = 5040.0_dp / T

    ! Do bound free calculation
    if ((wl(l) > lam_0) .or. (wl(l) < lam_min)) then
      xbf = 0.0_dp
      kbf = 0.0_dp
    else
      fbf = 0.0_dp
      do n = 1, 6
        fbf = fbf + Cn_bf(n) * (1.0_dp/wl(l) - 1.0_dp/lam_0)**((real(n,kind=dp)-1.0_dp)/2.0_dp)
      end do
      ! xbf is in [cm2 molecule-1] and kff is in [cm4 dyne-1] - convert to [cm-1] before CMCRT output
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



    kbf = xbf * VMR_lay(CIA_tab(s)%iVMR_3(1),z) * N_lay(z) ! * H- [molecule cm-3]
    !kbf = kbf * (VMR_lay(CIA_tab(s)%iVMR_3(2),z) * N_lay(z) * VMR_lay(CIA_tab(s)%iVMR_3(3),z) * N_lay(z)) * kb * T ! * P(e-) [dyne cm-2] * H [cm-3]
    kff = kff * (VMR_lay(CIA_tab(s)%iVMR_3(2),z) * N_lay(z) * VMR_lay(CIA_tab(s)%iVMR_3(3),z) * N_lay(z)) * kb * T ! * P(e-) [dyne cm-2] * H [cm-3]

    CIA_spec = kbf + kff

  end subroutine CIA_Hminus

end module CIA_tables_Hminus
