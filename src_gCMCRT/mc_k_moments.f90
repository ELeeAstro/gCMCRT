module mc_k_moments
  use mc_precision
  use mc_data_mod
  use mc_class_pac
  use mc_class_grid
  use cudafor
  implicit none

  real(dp), dimension(:), allocatable :: jp, hp, kp
  real(dp), dimension(:), allocatable :: jm, hm, km

  real(dp), dimension(:), allocatable, device :: jp_d, hp_d, kp_d
  real(dp), dimension(:), allocatable , device :: jm_d, hm_d, km_d

contains

  attributes(device) subroutine moments_1D(ph)
    implicit none

    type(pac), intent(in) :: ph

    real(dp) :: cost_mom, rp2
    integer :: istat

    ! If in spherical grid must scale to the radius of the grid size
    if (ph%geo == 2) then
      rp2 = ph%xp**2 + ph%yp**2 + ph%zp**2
      cost_mom = (ph%xp*ph%nxp + ph%yp*ph%nyp + ph%zp*ph%nzp)/sqrt(rp2)
    else
      cost_mom = ph%cost
    end if

    if (cost_mom > 0.0_dp) then
      ! Packet travelling upward ('+ve'), add to plus moments at this level

      istat = atomicadd(jp_d(ph%c(3)), 1.0_dp/cost_mom)
      istat = atomicadd(hp_d(ph%c(3)), 1.0_dp)
      istat = atomicadd(kp_d(ph%c(3)), cost_mom)
    else if (cost_mom < 0.0_dp) then
      ! Packet travelling downward ('-ve'), add to minus moments at this level

      istat = atomicadd(jm_d(ph%c(3)+1), -1.0_dp/cost_mom)
      istat = atomicadd(hm_d(ph%c(3)+1), -1.0_dp)
      istat = atomicadd(km_d(ph%c(3)+1), -cost_mom)
    end if

  end subroutine moments_1D

end module mc_k_moments
