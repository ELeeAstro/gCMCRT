module CIA_tables_fake_H2O_special
  use optools_data_mod
  implicit none

  real(dp), parameter :: sig0 = 1e-30_dp
  real(dp), parameter :: lam1 = 0.3_dp , lam2 = 10.0_dp
  real(dp), parameter :: alph = 4000.0_dp

contains

  subroutine Fake_H2O_special(s,l,z,CIA_spec)
    implicit none

    integer, intent(in) :: s, l, z
    real(dp), intent(out) :: CIA_spec 
   

    CIA_spec = sig0 * (1.0_dp + alph * (wl(l) - lam1)/(lam2 - lam1)) * 1e4_dp

    CIA_spec = CIA_spec * (VMR_lay(CIA_tab(s)%iVMR(1),z) * N_lay(z))

  end subroutine Fake_H2O_special

end module CIA_tables_fake_H2O_special
