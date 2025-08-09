!! Reads the main namnelist and adds variables to global data module
subroutine read_namelist()
  use mc_data_mod
  implicit none

  open(newunit=u_nml, file='CMCRT.nml', status='old', action='read')

  read(u_nml, nml=main)

  ! Send some flags to the device
  do_images_d = do_images ! images switch
  do_moments_d = do_moments ! calculate J, H, K moments
  do_trans_d = do_trans ! For peeloff - include transmision spectra counters
  do_cf_d = do_cf ! For contribution functions

  do_scat_loop_d = do_scat_loop

  do_g_bias_d = do_g_bias

  wght_deg_d = wght_deg ! Include weight degradiation

  do_LD_d = do_LD
  ilimb_d = ilimb
  LD_c_d(:) = LD_c(:)
  Rs_d = Rs
  inc_d = inc
  phase_d = phase
  sm_ax_d = sm_ax

  do_surf_d = do_surf
  alb_surf_d = alb_surf

  LHS_d = LHS

end subroutine read_namelist


subroutine read_CMCRT_par()
  use mc_data_mod
  implicit none

  integer :: u_par

  open(newunit=u_par, file='CMCRT.par', status='old', action='read')



end subroutine read_CMCRT_par
