module mc_k_vol_samp
  use mc_precision
  use mc_data_mod
  use mc_class_pac
  use mc_class_grid
  use curand_device
  use cudafor
  use libm
  implicit none


contains

  attributes(device) subroutine sph_samp_3D(i,j,k,ph)
    implicit none

    type(pac), intent(inout) :: ph
    integer, intent(in) :: i, j, k

    real(dp) :: r_samp, phi_samp, theta_samp

    ! set photon cell location
    ph%c(1) = i
    ph%c(2) = j
    ph%c(3) = k

    ! Sample a random position inside the cell
    r_samp = r_d(i) + cbrt(curand_uniform(ph%iseed)) * (r_d(i+1) - r_d(i))
    phi_samp = phi_d(j) + curand_uniform(ph%iseed) * (phi_d(j+1) - phi_d(j))
    theta_samp = acos(cos(theta_d(k)) + curand_uniform(ph%iseed) * (cos(theta_d(k+1)) - cos(theta_d(k))))

    ! Need to find xp, yp, zp from cell spherical coordinates
    ph%xp = r_samp * cos(phi_samp) * sin(theta_samp)
    ph%yp = r_samp * sin(phi_samp) * sin(theta_samp)
    ph%zp = r_samp * cos(theta_samp)


  end subroutine sph_samp_3D

  attributes(device) subroutine cart_samp_3D(i,j,k,ph)
    implicit none

    type(pac), intent(inout) :: ph
    integer, intent(in) :: i, j, k

    ph%xp = x_d(i) + &
      & curand_uniform(ph%iseed)*(x_d(i+1)-x_d(i)) - grid_d%x_max
    ph%yp = y_d(j) + &
      & curand_uniform(ph%iseed)*(y_d(j+1)-y_d(j)) - grid_d%y_max
    ph%zp = z_d(k) + &
      & curand_uniform(ph%iseed)*(z_d(k+1)-z_d(k)) - grid_d%z_max


    ph%c(1) = int(grid_d%n_x*(ph%xp+grid_d%x_max)/(2.0_dp*grid_d%x_max)) + 1
    ph%c(2) = int(grid_d%n_y*(ph%yp+grid_d%y_max)/(2.0_dp*grid_d%y_max)) + 1
    ph%c(3) = int(grid_d%n_z*(ph%zp+grid_d%z_max)/(2.0_dp*grid_d%z_max)) + 1


  end subroutine cart_samp_3D

end module mc_k_vol_samp
