module mc_k_limb_dark
  use mc_precision
  use mc_data_mod
  use mc_class_pac
  use mc_k_aux
  use curand_device
  use cudafor
  implicit none

  integer, device :: limb_dark_d
  real(dp), device :: delmu_d
  real(dp), dimension(:), allocatable, device :: xmu_d, limb_cdf_d

contains

  attributes(device) subroutine darkening(ph)
    implicit none

    type(pac), intent(inout) :: ph
    integer :: n
    real(dp) :: ran


    ran = curand_uniform(ph%iseed)

    call locate(limb_cdf_d,ran,n)

    ph%cost = xmu_d(n+1)-(limb_cdf_d(n+1)-ran)/ &
      & (limb_cdf_d(n+1)-limb_cdf_d(n))*delmu_d

    if (ph%cost > 1.0_dp .or. ph%cost < 0.0_dp) then
      print*,'dammit'
    endif

    ph%sint = sqrt(1.0_dp - ph%cost**2)

  end subroutine darkening

end module mc_k_limb_dark
