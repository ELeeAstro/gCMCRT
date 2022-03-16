module mc_k_raytrace
  use mc_precision
  use mc_data_mod
  use mc_class_pac
  use mc_class_grid
  use mc_class_imag
  use mc_k_findwall_cart
  use mc_k_findwall_sph
  use ieee_arithmetic
  use cudafor
  implicit none

contains

  attributes(device) subroutine raytrace_sph_3D(ray)
    implicit none

    type(pac), intent(inout) :: ray

    integer :: n
    integer, dimension(3) :: ioffset
    real(dp) :: taurun,taucell,d,dcell
    real(dp) :: r2,rcosa,r2min,smax,d1, deps

    !! Begin the tau integration
    ray%tau = 0.0_dp
    ray%p_flag = 0
    d = 0.0_dp

    !Ray direction is now observation direction
    ray%nxp = im_d%obsx; ray%nyp = im_d%obsy ; ray%nzp = im_d%obsz

    ! calculate smax -- maximum distance photon can travel.  smax is
    !  either to the stellar surface or to the edge of the envelope
    !  depending on the photon direction
    r2 = ray%xp**2 + ray%yp**2 + ray%zp**2
    ! rcoa towards observation direction
    rcosa = ray%xp * ray%nxp + ray%yp * ray%nyp + ray%zp * ray%nzp
    r2min = max(r2-rcosa**2, 0.0_dp)
    if (rcosa > 0.0_dp) then  ! is photon traveling outward?
      smax = sqrt(grid_d%r2_max-r2min)-rcosa !  yes, find dist to rmax
    else !no, photon traveling inward
      if (r2min > 1.0_dp) then  !does sight-line hit star?
        smax = sqrt(grid_d%r2_max-r2min)-rcosa !  no, find dist to rmax
      else
        smax = -rcosa - sqrt(1.0_dp-r2min)  !  yes, find dist to star
         ray%p_flag = -1  ! i.e., photon occulted
        return
      endif
    endif

    ! tau = 0 if max distance is small
    if (smax < 1.0e-12_dp) then
      ray%p_flag = -2
      return
    endif

    ! integrate through grid ********
    do while (ray%p_flag == 0 .and. d < (0.999990_dp * smax))

      ! find distance to next cell, dcell  *
      call findwall_sph(ray, ioffset, dcell)

      if (dcell <= 0.0_dp .or. ieee_is_nan(dcell)) then
        ray%p_flag = -22
        print*, 'tauint peeloff: dcell < 0', ray%id, dcell
        exit
      end if

      ! update total optical depth.  optical depth to next cell wall is
      ! taucell= (distance to cell)*(opacity of current cell)
      taucell = dcell * rhokap_d(ray%ig,ray%c(1),ray%c(2),ray%c(3))

      ! Small offset for dcell
      deps = (r_d(ray%c(1)+1) - r_d(ray%c(1)))*1.0e-12_dp
      deps = max(deps, 1.0e-12_dp)
      d1 = dcell + deps

        ray%xp = ray%xp + d1 * ray%nxp
        ray%yp = ray%yp + d1 * ray%nyp
        ray%zp = ray%zp + d1 * ray%nzp
        ! update packet cell
        ray%c(1) = ray%c(1) + ioffset(1)
        ray%c(2) = ray%c(2) + ioffset(2)
        ray%c(3) = ray%c(3) + ioffset(3)

        ray%tau = ray%tau + taucell
        d = d + d1

      ! Detect if entering surface or escaping atmosphere
      if (ray%c(1) < 1) then
        ray%p_flag = -11
        exit
      else if (ray%c(1) >= grid_d%n_lev) then
        ray%p_flag = 0
        exit
      end if

      ! Detect if longitude passing from 1->NY or NY->1
      if (ray%c(2) >= grid_d%n_phi) then
         ray%c(2) = 1
      else if (ray%c(2) < 1) then
         ray%c(2) = grid_d%n_phi-1
      end if

      ! Detect if out of latitude bounds
      if ((ray%c(3) >= grid_d%n_theta) .or. (ray%c(3) < 1)) then
        ! Terminate integration of packet
        ray%p_flag = -3
        exit
      end if

    end do

  end subroutine raytrace_sph_3D

  attributes(device) subroutine raytrace_cart_3D(ray)
    implicit none

    type(pac), intent(inout) :: ray

    integer, dimension(3) :: ioffset
    real(dp) :: d, d1, dcell, taucell, dsx, dsy, dsz, smax

    !Ray direction is now observation direction
    ray%nxp = im_d%obsx; ray%nyp = im_d%obsy ; ray%nzp = im_d%obsz

    ray%xp = ray%xp + grid_d%x_max
    ray%yp = ray%yp + grid_d%y_max
    ray%zp = ray%zp + grid_d%z_max

    !! Begin the tau integration
    ray%tau = 0.0_dp
    ray%p_flag = 0

    if (ray%nxp > 0.0_dp) then
      dsx = (2.0_dp*grid_d%x_max-ray%xp)/ray%nxp
    else if (ray%nxp < 0.0_dp) then
      dsx = -ray%xp/ray%nxp
    else if(ray%nxp == 0.0_dp) then
      dsx = 1.0e2_dp*grid_d%x_max
    endif

    if (ray%nyp > 0.0_dp) then
      dsy = (2.0_dp*grid_d%y_max-ray%yp)/ray%nyp
    else if (ray%nyp < 0.0_dp) then
      dsy = -ray%yp/ray%nyp
    else if(ray%nyp == 0.0_dp) then
      dsy = 1.0e2_dp*grid_d%y_max
    endif

    if (ray%nzp > 0.0_dp) then
      dsz = (2.0_dp*grid_d%z_max-ray%zp)/ray%nzp
    else if (ray%nzp < 0.0_dp) then
      dsz = -ray%zp/ray%nzp
    else if(ray%nzp == 0.0_dp) then
      dsz = 1.0e2_dp*grid_d%z_max
    endif

    smax = min(dsx,dsy,dsz)

    if (smax < 1.0e-12_dp) then
      ray%tau = 0.0_dp
      ray%p_flag = -1
      return
    endif

    d = 0.0_dp
    ! integrate through grid ********
    do while (d < 0.999_dp*smax)

      ! find distance to next cell, dcell  *
      call findwall_cart(ray, ioffset, dcell)

      ! update total optical depth.  optical depth to next cell wall is
      ! taucell= (distance to cell)*(opacity of current cell)
      taucell = dcell * rhokap_d(ray%ig,ray%c(1),ray%c(2),ray%c(3))

      d1 = dcell
      ray%xp = ray%xp + d1 * ray%nxp
      ray%yp = ray%yp + d1 * ray%nyp
      ray%zp = ray%zp + d1 * ray%nzp

      ray%c(1) = int(grid_d%n_x*(ray%xp)/(2.0_dp*grid_d%x_max)) + 1
      ray%c(2) = int(grid_d%n_y*(ray%yp)/(2.0_dp*grid_d%y_max)) + 1
      ray%c(3) = int(grid_d%n_z*(ray%zp)/(2.0_dp*grid_d%z_max)) + 1

      ! update ray cell
      ! ray%c(1) = ray%c(1) + ioffset(1)
      ! ray%c(2) = ray%c(2) + ioffset(2)
      ! ray%c(3) = ray%c(3) + ioffset(3)

      ray%tau = ray%tau + taucell
      d = d + d1

      if ((ray%c(1) >= grid_d%n_x) .or. (ray%c(1) < 1)) then
        ray%p_flag = 0
        exit
      end if
      if ((ray%c(2) >= grid_d%n_y) .or. (ray%c(2) < 1)) then
        ray%p_flag = 0
        exit
      end if
      if ((ray%c(3) >= grid_d%n_z) .or. (ray%c(3) < 1)) then
        ray%p_flag = 0
        exit
      end if

    end do

  end subroutine raytrace_cart_3D

end module mc_k_raytrace
