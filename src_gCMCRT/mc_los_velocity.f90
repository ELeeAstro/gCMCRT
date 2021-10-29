module mc_los_velocity
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_class_imag
  implicit none


contains

  subroutine compute_vlos()
    implicit none

    integer :: i,j,k
    real(dp) :: u, v, w
    real(dp) :: phi_cell, theta_cell, Rz_cell, phi_v, theta_v
    real(dp) :: angfreq, vorb, vsys, vlos, orb_phase

    ! Compute line of sight velocity for each grid cell
    allocate(v_los(grid%n_lay,grid%n_phi-1,grid%n_theta-1))

    do k = 1, grid%n_theta-1
      do j = 1, grid%n_phi-1
        do i = 1, grid%n_lay

          ! Wind components for cell (i,j,k), in [cm/s]
          u = u_wind(i,j,k)
          v = v_wind(i,j,k)
          w = w_wind(i,j,k)

          phi_cell = phiarr(j) + 0.5_dp * dphi           ! Central phi coordinate of cell
          theta_cell = pi - (thetarr(k) + 0.5_dp*dtheta) ! Converted to -/+ pi radian
          Rz_cell = 0.5_dp*(H(i) + H(i+1))               ! Center of cell in radius assumed

          ! Viewing vector components in radians
          phi_v = im%vphi*pi/180.0_dp
          theta_v = im%vtheta*pi/180.0_dp

          ! Compute angular frequency, orbital velocity and systemic velocity
          angfreq = twopi / (orbital_period * daysec)  ! units of [rad/s]
          vorb = (twopi * Au * sm_ax) / (orbital_period * daysec) ! units of [cm/s]
          vsys = systemic_velocity * 1e5_dp ! units of [cm/s]

          vlos = 0.0_dp

          ! Add velocity component due to atmospheric winds (equation 7)
          if (winds_on .eqv. .True.) then
              vlos = vlos + u*sin(theta_v)*sin(phi_cell - phi_v) &
                  & + ( v*cos(theta_cell) - w*sin(theta_cell) )*sin(theta_v)*cos(phi_cell - phi_v) &
                  & - (v*sin(theta_cell) + w*cos(theta_cell) )*cos(theta_v)
          end if

          ! Add velocity component due to planetary rotation (equation 10)
          if (rotation_on .eqv. .True.) then
              vlos = vlos + angfreq*Rz_cell*sin(theta_cell)*sin(theta_v)*sin(phi_cell - phi_v)
          end if

          ! Add velocity component due to orbital motion (equation 12)
          if (orbit_on .eqv. .True.) then
              orb_phase = pi - phi_v  ! phi_v already expressed in radians
              vlos = vlos + vsys + vorb*sin(theta_v)*sin(orb_phase)
          end if

          v_los(i,j,k) = vlos

        end do
      end do
    end do


  end subroutine compute_vlos

end module mc_los_velocity
