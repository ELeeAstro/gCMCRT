module mc_k_findcell
  use mc_data_mod
  use mc_class_pac
  use mc_class_grid
  use mc_k_aux
  use cudafor
  implicit none

contains

  attributes(device) subroutine findcell(ph)

    implicit none

    type(pac), intent(inout) :: ph

    real(dp) :: rp, rphi, rtheta
    real(dp) :: mu

    ! Find the cell number of current packet position.
    rp = sqrt(ph%xp**2 + ph%yp**2 + ph%zp**2)

    ! Guard against pathological zero-radius packets.
    if (rp <= 0.0_dp) then
       ph%c(1) = -1
       ph%c(2) = -1
       ph%c(3) = -1
       return
    end if

    ! Radial cell.
    call locate(r_d(:), rp, ph%c(1))

    ! Longitude cell.
    rphi = atan2(ph%yp, ph%xp)

    if (rphi < 0.0_dp) then
       rphi = rphi + twopi
    end if

    call locate(phi_d(:), rphi, ph%c(2))

    ! Theta cell.
    ! Clip to avoid acos(NaN) from tiny round-off outside [-1, 1].
    mu = ph%zp / rp
    mu = max(-1.0_dp, min(1.0_dp, mu))

    rtheta = acos(mu)

    call locate(theta_d(:), rtheta, ph%c(3))

  end subroutine findcell  

end module  mc_k_findcell
