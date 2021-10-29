module mc_k_RR
  use mc_precision
  use mc_data_mod
  use mc_class_pac
  use cudafor
  use curand_device
  implicit none

  real(dp), parameter :: RR_thresh = 0.1_dp
  real(dp), parameter :: RR_boost = 10.0_dp
  real(dp), parameter :: iRR_cut = 1.0e-3_dp

contains

  attributes(device) subroutine RR_test(ph)
    implicit none

    type(pac), intent(inout) :: ph

    if (ph%wght > iRR_cut) then
      ph%p_flag = 0
    else if ((ph%wght < iRR_cut) .and. (curand_uniform(ph%iseed) < RR_thresh)) then
      ph%wght = ph%wght * RR_boost
      ph%p_flag = 0
    else
      ph%p_flag = -44
    end if


  end subroutine RR_test

end module mc_k_RR
