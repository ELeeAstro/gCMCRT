module mc_k_findwall_cart
  use mc_precision
  use mc_data_mod
  use mc_class_pac
  use mc_class_grid
  implicit none

contains

  attributes(device) subroutine findwall_cart(ph, ioffset, dcell)
    implicit none

    type(pac), intent(in) :: ph

    integer, dimension(3), intent(out) :: ioffset
    real(dp), intent(out) :: dcell

    real(dp) :: dx, dy, dz

    !! Find the distance to the x,y,z edges in the direction of travel
    if (ph%nxp > 0.0_dp) then
      dx = (x_d(ph%c(1)+1)-ph%xp)/ph%nxp
      ioffset(1) = 1
    else if (ph%nxp < 0.0_dp) then
      dx = (x_d(ph%c(1))-ph%xp)/ph%nxp
      ioffset(1) = -1
    else if (ph%nxp == 0.0_dp) then
      dx = 1.0e2_dp*grid_d%x_max
      ioffset(1) = 0
    end if

    if (ph%nyp > 0.0_dp) then
      dy = (y_d(ph%c(2)+1)-ph%yp)/ph%nyp
      ioffset(2) = 1
    else if(ph%nyp < 0.0_dp) then
      dy = (y_d(ph%c(2))-ph%yp)/ph%nyp
      ioffset(2) = -1
    else if(ph%nyp == 0.0_dp) then
      dy = 1.0e2_dp*grid_d%y_max
      ioffset(2) = 0
    end if

    if (ph%nzp > 0.0_dp) then
      dz = (z_d(ph%c(3)+1)-ph%zp)/ph%nzp
      ioffset(3) = 1
    else if(ph%nzp < 0.0_dp) then
      dz = (z_d(ph%c(3))-ph%zp)/ph%nzp
      ioffset(3) = -1
    else if(ph%nzp == 0.0_dp) then
      dz = 1.0e2_dp*grid_d%z_max
      ioffset(3) = 0
    end if

    !! Check for zero distances
    if (dx == 0.0_dp) then
      dx = 1.0e2_dp*grid_d%x_max
      ioffset(1) = 0
    end if
    if (dy == 0.0_dp) then
      dy = 1.0e2_dp*grid_d%y_max
      ioffset(2) = 0
    end if
    if (dz == 0.0_dp) then
      dz = 1.0e2_dp*grid_d%z_max
      ioffset(3) = 0
    end if

    ! Distance to cell wall is minium of the 3 distances
    dcell = min(dx,dy,dz)

    if(dx < 0.0_dp) dcell=min(dy,dz)
    if(dy < 0.0_dp) dcell=min(dx,dz)
    if(dz < 0.0_dp) dcell=min(dx,dy)


  end subroutine findwall_cart


end module mc_k_findwall_cart
