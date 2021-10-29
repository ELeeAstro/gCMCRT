module mc_k_gord_samp
  use mc_precision
  use mc_data_mod
  use mc_class_pac
  use mc_class_grid
  use cudafor
  use curand_device
  use mc_k_aux
  implicit none



contains

  attributes(device) subroutine gord_samp(ph)
    implicit none

    type(pac), intent(inout) :: ph

    real(dp) :: g_test

    if (ng_d == 1) then
      ph%ig = 1
    else  
      ! Non-biased sampling (reflection + transmission spectra)
      g_test = curand_uniform(ph%iseed)
      call locate(gord_cdf_d(:),g_test,ph%ig)
      if (ph%ig < 1) then
        ph%ig = 1
      else if (ph%ig > ng_d) then
        ph%ig = ng_d
      end if
    end if
  
  end subroutine gord_samp


  attributes(device) subroutine gord_samp_bias(ph)
    implicit none

    type(pac), intent(inout) :: ph

    real(dp) :: g_test

    if (ng_d == 1) then
      ph%ig = 1
    else
      g_test = curand_uniform(ph%iseed)
      ph%ig = int(g_test * real(ng_d,dp)) + 1
      if (ph%ig < 1) then
        ph%ig = 1
      else if (ph%ig > ng_d) then
        ph%ig = ng_d
      end if
      ph%wght = ph%wght * (real(ng_d,dp) * gord_w_d(ph%ig))
    end if

  end subroutine gord_samp_bias

  attributes(device) subroutine gord_samp_cell(ph)
    implicit none
    type(pac), intent(inout) :: ph

    real(dp) :: g_test

    ! Non-biased sampling (emission spectra), has to be weighted by cell luminosity cdf per g-ordinance
    if (ng_d == 1) then
      ph%ig = 1
    else
      g_test = curand_uniform(ph%iseed)
      call locate(cell_gord_cdf_d(:,ph%c(1),ph%c(2),ph%c(3)),g_test,ph%ig)
      if (ph%ig < 1) then
        ph%ig = 1
      else if (ph%ig > ng_d) then
        ph%ig = ng_d
      end if
    end if

  end subroutine gord_samp_cell

  attributes(device) subroutine gord_samp_cell_bias(ph)
    implicit none

    type(pac), intent(inout) :: ph

    real(dp) :: g_test
 
    ! Biased sampling (emission spectra)   
    if (ng_d == 1) then
      ph%ig = 1
    else if (curand_uniform(ph%iseed) < (1.0_dp - xi_g)) then
      g_test = curand_uniform(ph%iseed)
      call locate(cell_gord_cdf_d(:,ph%c(1),ph%c(2),ph%c(3)),g_test,ph%ig)
    else
      g_test = curand_uniform(ph%iseed)
      ph%ig = int(g_test * real(ng_d,dp)) + 1
    end if

    if (ph%ig < 1) then
      ph%ig = 1
    else if (ph%ig > ng_d) then
      ph%ig = ng_d
    end if

    ph%wght = ph%wght * cell_gord_wght_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))

  end subroutine gord_samp_cell_bias


end module mc_k_gord_samp
