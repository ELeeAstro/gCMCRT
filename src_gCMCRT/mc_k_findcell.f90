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

    ! Find the cell number of current packet position

    rp = sqrt(ph%xp**2 + ph%yp**2 + ph%zp**2) ! Radial distance
    call locate(r_d(:),rp,ph%c(1)) ! Find radial cell

    rphi = atan2(ph%yp,ph%xp) ! phi (longitude) value
    if (rphi < 0.0_dp) then  !if negative, add 2*pi
       rphi = rphi + twopi
    endif
    call locate(phi_d(:),rphi,ph%c(2)) !Find phi cell

    rtheta = acos(ph%zp/rp)  ! theta (laittude) value
    call locate(theta_d(:),rtheta,ph%c(3)) !Find theta cell


  end subroutine findcell

end module  mc_k_findcell
