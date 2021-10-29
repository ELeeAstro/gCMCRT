module mc_k_tauint
  use mc_precision
  use mc_data_mod
  use mc_class_pac
  use mc_class_grid
  use mc_k_moments
  use mc_k_findcell
  use mc_k_findwall_cart
  use mc_k_findwall_sph
  use ieee_arithmetic
  use cudafor
  implicit none

contains

  attributes(device) subroutine tauint_1D_sph(ph)
    implicit none

    type(pac), intent(inout) :: ph

    real(dp) :: dr, taucell, d1
    integer :: zoffset

    !! Begin the tau integration
    ph%tau = 0.0_dp
    ph%p_flag = 0

    do while (ph%tau < ph%tau_p)

      call findwall_1D_sph(ph,zoffset,dr)

      !! Calculate optical depth to level
      taucell = dr * rhokap_1D_d(ph%c(3))

      !! Check if packet ends path in this layer
      if ((ph%tau + taucell) >= ph%tau_p) then
        ! Packet stops in this cell - move distance then exit loop
        d1 = (ph%tau_p-ph%tau)/rhokap_1D_d(ph%c(3))
        ph%xp = ph%xp + d1 * ph%nxp
        ph%yp = ph%yp + d1 * ph%nyp
        ph%zp = ph%zp + d1 * ph%nzp
        exit
      else
        ! Packet continues to level edge - update position, cell index and tau
        ! slightly increase the distance to avoid roundoff error
        ph%xp = ph%xp + (dr + 1.0e-12_dp) * ph%nxp
        ph%yp = ph%yp + (dr + 1.0e-12_dp) * ph%nyp
        ph%zp = ph%zp + (dr + 1.0e-12_dp) * ph%nzp

        ph%c(3) = ph%c(3) + zoffset

        if (do_moments_d .eqv. .True.) then
          call moments_1D(ph)
        end if

        ph%tau = ph%tau + taucell

      end if

      ! Check is packet has exited the domain
      if (ph%c(3) >= grid_d%n_lev) then
        ph%p_flag = 1
        exit
      else if (ph%c(3) < 1) then
        ph%p_flag = 2
        exit
      end if

    end do

  end subroutine tauint_1D_sph

  attributes(device) subroutine tauint_1D_pp(ph)
    implicit none

    type(pac), intent(inout) :: ph

    real(dp) :: dsz, taucell, d1
    integer :: zoffset

    !! Begin the tau integration
    ph%tau = 0.0_dp
    ph%p_flag = 0

    do while (ph%tau < ph%tau_p)

      !! Calculate dsz, the distance to the next vertical level
      if (ph%nzp > 0.0_dp) then
        ! Packet travelling upward, find distance to upper level
        dsz = (z_d(ph%c(3)+1)-ph%zp)/ph%nzp
        zoffset = 1
      else if (ph%nzp < 0.0_dp) then
        ! Packet travelling downward, find distance to lower level
        dsz = (z_d(ph%c(3))-ph%zp)/ph%nzp
        zoffset = -1
      else if(ph%nzp == 0.0_dp) then
        ! Packet travelling directly in z plane, use a large distance to move packet
        dsz = 1.0e2_dp*grid_d%z_max
        zoffset = 0
      endif

      !! Calculate optical depth to level
      taucell = dsz * rhokap_1D_d(ph%c(3))

      !! Check if packet ends path in this layer
      if ((ph%tau + taucell) >= ph%tau_p) then
        ! Packet stops in this cell - move distance then exit loop
        d1 = (ph%tau_p-ph%tau)/rhokap_1D_d(ph%c(3))
        ph%xp = ph%xp + d1 * ph%nxp
        ph%yp = ph%yp + d1 * ph%nyp
        ph%zp = ph%zp + d1 * ph%nzp
        exit
      else
        ! Packet continues to level edge - update position, cell index and tau
        ph%xp = ph%xp + dsz * ph%nxp
        ph%yp = ph%yp + dsz * ph%nyp
        ph%zp = ph%zp + dsz * ph%nzp

        ph%c(3) = ph%c(3) + zoffset

        if (do_moments_d .eqv. .True.) then
          call moments_1D(ph)
        end if

        ph%tau = ph%tau + taucell

      end if

      ! Check is packet has exited the domain
      if (ph%c(3) >= grid_d%n_lev) then
        ph%p_flag = 1
        exit
      else if (ph%c(3) < 1) then
        ph%p_flag = 2
        exit
      end if

    end do

  end subroutine tauint_1D_pp

  attributes(device) subroutine tauint_cart_3D(ph)
    implicit none

    type(pac), intent(inout) :: ph

    integer, dimension(3) :: ioffset
    real(dp) :: taurun,taucell,d,dcell
    real(dp) :: r2,rcosa,r2min,smax
    real(dp) :: deps, dsx, dsy, dsz, xcur, ycur, zcur, d1

    !! Begin the tau integration
    ph%tau = 0.0_dp
    ph%p_flag = 0

    ph%xp = ph%xp + grid_d%x_max
    ph%yp = ph%yp + grid_d%y_max
    ph%zp = ph%zp + grid_d%z_max


    ! integrate through grid ********
    do while (ph%tau < ph%tau_p)

      ! find distance to next cell, dcell  *
      call findwall_cart(ph, ioffset, dcell)

      ! update total optical depth.  optical depth to next cell wall is
      ! taucell= (distance to cell)*(opacity of current cell)
      taucell = dcell * rhokap_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))

      if ((ph%tau + taucell) >= ph%tau_p) then
        d1 = (ph%tau_p-ph%tau)/rhokap_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))
        ph%xp = ph%xp + d1 * ph%nxp
        ph%yp = ph%yp + d1 * ph%nyp
        ph%zp = ph%zp + d1 * ph%nzp
        ph%c(1) = int(grid_d%n_x*(ph%xp)/(2.0_dp*grid_d%x_max)) + 1
        ph%c(2) = int(grid_d%n_y*(ph%yp)/(2.0_dp*grid_d%y_max)) + 1
        ph%c(3) = int(grid_d%n_z*(ph%zp)/(2.0_dp*grid_d%z_max)) + 1
        exit
      else
        d1 = dcell + (dcell * 1e-12_dp)
        ph%xp = ph%xp + d1 * ph%nxp
        ph%yp = ph%yp + d1 * ph%nyp
        ph%zp = ph%zp + d1 * ph%nzp
        ! update packet cell
        ! ph%c(1) = ph%c(1) + ioffset(1)
        ! ph%c(2) = ph%c(2) + ioffset(2)
        ! ph%c(3) = ph%c(3) + ioffset(3)

        ph%c(1) = int(grid_d%n_x*(ph%xp)/(2.0_dp*grid_d%x_max)) + 1
        ph%c(2) = int(grid_d%n_y*(ph%yp)/(2.0_dp*grid_d%y_max)) + 1
        ph%c(3) = int(grid_d%n_z*(ph%zp)/(2.0_dp*grid_d%z_max)) + 1

        ph%tau = ph%tau + taucell
      end if

      if ((ph%c(1) >= grid_d%n_x) .or. (ph%c(1) < 1)) then
        ph%p_flag = 1
        exit
      end if
      if ((ph%c(2) >= grid_d%n_y) .or. (ph%c(2) < 1)) then
        ph%p_flag = 2
        exit
      end if
      if ((ph%c(3) >= grid_d%n_z) .or. (ph%c(3) < 1)) then
        ph%p_flag = 3
        exit
      end if

    end do

    ph%xp = ph%xp - grid_d%x_max
    ph%yp = ph%yp - grid_d%y_max
    ph%zp = ph%zp - grid_d%z_max

    ph%c(1) = int(grid_d%n_x*(ph%xp+grid_d%x_max)/(2.0_dp*grid_d%x_max)) + 1
    ph%c(2) = int(grid_d%n_y*(ph%yp+grid_d%y_max)/(2.0_dp*grid_d%y_max)) + 1
    ph%c(3) = int(grid_d%n_z*(ph%zp+grid_d%z_max)/(2.0_dp*grid_d%z_max)) + 1



  end subroutine tauint_cart_3D


  attributes(device) subroutine tauint_sph_3D(ph)
    implicit none

    type(pac), intent(inout) :: ph

    integer, dimension(3) :: ioffset
    real(dp) :: taurun,taucell,dcell
    real(dp) :: r2,rcosa,r2min,smax
    real(dp) :: xcur, ycur, zcur, d1, d, deps

    !! Begin the tau integration
    ph%tau = 0.0_dp
    ph%p_flag = 0
    d = 0.0_dp

    ! Calc smax - max distance photon can travel
    ! To core of planet or to outer envelope - direction dependent
    r2 = ph%xp**2 + ph%yp**2 + ph%zp**2
    rcosa = ph%xp*ph%nxp + ph%yp*ph%nyp + ph%zp*ph%nzp
    r2min = max(r2 - rcosa**2, 0.0_dp)

    if (rcosa > 0.0_dp) then  !Photon traveling outward
      smax = sqrt(grid_d%r2_max-r2min) - rcosa !yes, find dist to rmax
    else ! Photon traveling inward
      if (r2min > 1.0_dp) then !does sight-line hit star?
        smax = sqrt(grid_d%r2_max-r2min) - rcosa ! no, find dist to rmax
      else
        smax = -rcosa-sqrt(1.0_dp-r2min) !Yes, find dist to core
      end if
    end if

    ! if smax is small, tauflag = 1
    ! Assume absoption/ exiting atmosphere
    if (smax < 1.0e-12_dp) then
      ph%p_flag = 555
      return
    end if


    ! integrate through grid ********
    do while (ph%tau < ph%tau_p .and. (d < (0.999990_dp*smax)))

      ! find distance to next cell, dcell  *
      call findwall_sph(ph, ioffset, dcell)

      ! Distance is negative
      if (dcell <= 0.0_dp .or. ieee_is_nan(dcell)) then
        ph%p_flag = -2
        print*, 'tauint : dcell < 0', dcell, ph%c(1), ioffset(1)
        return
      end if

      ! update total optical depth.  optical depth to next cell wall is
      ! taucell= (distance to cell)*(opacity of current cell)
      taucell = dcell * rhokap_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))

      if ((ph%tau + taucell) >= ph%tau_p) then
        d1 = (ph%tau_p-ph%tau)/rhokap_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))
        ph%xp = ph%xp + d1 * ph%nxp
        ph%yp = ph%yp + d1 * ph%nyp
        ph%zp = ph%zp + d1 * ph%nzp

        if (wght_deg_d .eqv. .True.) then
          ph%wght = ph%wght * &
          & exp(-(rhokap_d(ph%ig,ph%c(1),ph%c(2),ph%c(3)) * (1.0_dp - ssa_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))) * d1))
        end if

        exit
      else
        ! Small offset for dcell
        deps = (r_d(ph%c(1)+1) - r_d(ph%c(1)))*1.0e-12_dp
        deps = max(deps, 1.0e-12_dp)
        d1 = dcell + deps

        if (wght_deg_d .eqv. .True.) then
          ph%wght = ph%wght * &
          & exp(-(rhokap_d(ph%ig,ph%c(1),ph%c(2),ph%c(3)) * (1.0_dp - ssa_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))) * d1))
        end if

        ph%xp = ph%xp + d1 * ph%nxp
        ph%yp = ph%yp + d1 * ph%nyp
        ph%zp = ph%zp + d1 * ph%nzp

        ! update packet cell
        ph%c(1) = ph%c(1) + ioffset(1)
        ph%c(2) = ph%c(2) + ioffset(2)
        ph%c(3) = ph%c(3) + ioffset(3)

        ph%tau = ph%tau + taucell
        d = d + d1
      end if

      ! Detect if enetering surface or escaping atmosphere
      if (ph%c(1) < 1) then
        ph%p_flag = 1
        exit
      else if (ph%c(1) >= grid_d%n_lev) then
        ph%p_flag = 2
        exit
      end if

      ! Detect if longitude passing from 1->NY or NY->1
      if (ph%c(2) >= grid_d%n_phi) then
         ph%c(2) = 1
      else if (ph%c(2) < 1) then
         ph%c(2) = grid_d%n_phi-1
      end if

      ! Detect if out of latitude bounds
      if ((ph%c(3) >= grid_d%n_theta) .or. (ph%c(3) < 1)) then
        ! Terminate integration of packet
        ph%p_flag = -4
        exit
      end if

    end do

    ! Detect if out of latitude bounds
    if ((ph%c(3) >= grid_d%n_theta) .or. (ph%c(3) < 1)) then
      ! Terminate integration of packet
      ph%p_flag = -4
    end if


  end subroutine tauint_sph_3D

end module mc_k_tauint
