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

    logical :: hit_surface
    integer, dimension(3) :: ioffset
    real(dp) :: taurun,taucell,d,dcell
    real(dp) :: r2,rcosa,r2min,smax
    real(dp) :: deps, dsx, dsy, dsz, xcur, ycur, zcur, d1
    real(dp) :: r_surf2

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

    logical :: hit_surface
    integer, dimension(3) :: ioffset

    real(dp) :: taucell
    real(dp) :: dcell, dstep, dleft, dmove
    real(dp) :: r2, rcosa, r2min, smax
    real(dp) :: d, deps, d_absorb
    real(dp) :: r_surf2, eps_path
    real(dp) :: kappa, kappa_abs

    ! Begin the tau integration.
    ph%tau    = 0.0_dp
    ph%p_flag = 0
    d         = 0.0_dp

    ! Calculate maximum geometric distance.
    r2     = ph%xp**2 + ph%yp**2 + ph%zp**2
    rcosa  = ph%xp*ph%nxp + ph%yp*ph%nyp + ph%zp*ph%nzp
    r2min  = max(r2 - rcosa**2, 0.0_dp)

    r_surf2 = grid_d%r_min * grid_d%r_min

    if (rcosa > 0.0_dp) then
       hit_surface = .False.
       smax = sqrt(max(grid_d%r2_max - r2min, 0.0_dp)) - rcosa
    else
       if (r2min > r_surf2) then
          hit_surface = .False.
          smax = sqrt(max(grid_d%r2_max - r2min, 0.0_dp)) - rcosa
       else
          hit_surface = .True.
          smax = -rcosa - sqrt(max(r_surf2 - r2min, 0.0_dp))
       end if
    end if

    eps_path = 1.0e-12_dp * max(1.0_dp, abs(smax))

    if (smax <= eps_path) then
       ph%p_flag = 555
       return
    end if

    do while (ph%tau < ph%tau_p .and. d < smax - eps_path)

       ! Defensive guard before accessing opacity and scattering arrays.
       if ((ph%c(1) < 1) .or. (ph%c(1) >= grid_d%n_lev) .or. &
           (ph%c(2) < 1) .or. (ph%c(2) >= grid_d%n_phi) .or. &
           (ph%c(3) < 1) .or. (ph%c(3) >= grid_d%n_theta)) then
          ph%p_flag = -777
          return
       end if

       call findwall_sph(ph, ioffset, dcell)

       if (dcell <= 0.0_dp .or. ieee_is_nan(dcell)) then
          ph%p_flag = -2
          print*, 'tauint_sph_3D: invalid dcell', dcell, &
                  ph%c(1), ph%c(2), ph%c(3), &
                  ioffset(1), ioffset(2), ioffset(3)
          return
       end if

       dleft = smax - d
       dstep = min(dcell, dleft)

       kappa = rhokap_d(ph%ig, ph%c(1), ph%c(2), ph%c(3))
       taucell = dstep * kappa

       ! Sampled interaction occurs inside the current physical segment.
       if ((ph%tau + taucell) >= ph%tau_p) then

          if (kappa <= 0.0_dp) then
             ph%p_flag = -23
             print*, 'tauint_sph_3D: non-positive kappa at sampled interaction', &
                     ph%id, kappa, ph%c(1), ph%c(2), ph%c(3)
             return
          end if

          d_absorb = (ph%tau_p - ph%tau) / kappa

          ph%xp = ph%xp + d_absorb * ph%nxp
          ph%yp = ph%yp + d_absorb * ph%nyp
          ph%zp = ph%zp + d_absorb * ph%nzp

          if (wght_deg_d .eqv. .True.) then
             kappa_abs = kappa * (1.0_dp - ssa_d(ph%ig, ph%c(1), ph%c(2), ph%c(3)))
             ph%wght = ph%wght * exp(-kappa_abs * d_absorb)
          end if

          ph%tau = ph%tau_p
          d = d + d_absorb

          exit

       end if

       ! No sampled interaction in this physical segment.
       ph%tau = ph%tau + taucell

       if (wght_deg_d .eqv. .True.) then
          kappa_abs = kappa * (1.0_dp - ssa_d(ph%ig, ph%c(1), ph%c(2), ph%c(3)))
          ph%wght = ph%wght * exp(-kappa_abs * dstep)
       end if

       ! If this segment ended at the geometric boundary rather than a cell
       ! wall, classify the packet and stop.
       if (dcell >= dleft) then

          ph%xp = ph%xp + dstep * ph%nxp
          ph%yp = ph%yp + dstep * ph%nyp
          ph%zp = ph%zp + dstep * ph%nzp

          d = smax

          if (hit_surface .eqv. .True.) then
             ph%p_flag = 1      ! inner boundary / surface
          else
             ph%p_flag = 2      ! escaped outer boundary
          end if

          return

       end if

       ! Normal cell-wall crossing: use a small numerical nudge for position
       ! and cell-index transition only. Do not count it as physical optical
       ! depth or physical path length.
       deps = (r_d(ph%c(1)+1) - r_d(ph%c(1))) * 1.0e-12_dp
       deps = max(deps, 1.0e-12_dp)

       dmove = dstep + deps

       ph%xp = ph%xp + dmove * ph%nxp
       ph%yp = ph%yp + dmove * ph%nyp
       ph%zp = ph%zp + dmove * ph%nzp

       ph%c(1) = ph%c(1) + ioffset(1)
       ph%c(2) = ph%c(2) + ioffset(2)
       ph%c(3) = ph%c(3) + ioffset(3)

       d = d + dstep

       ! Radial boundary tests.
       if (ph%c(1) < 1) then
          ph%p_flag = 1
          exit
       else if (ph%c(1) >= grid_d%n_lev) then
          ph%p_flag = 2
          exit
       end if

       ! Periodic longitude wrapping.
       if (ph%c(2) >= grid_d%n_phi) then
          ph%c(2) = 1
       else if (ph%c(2) < 1) then
          ph%c(2) = grid_d%n_phi - 1
       end if

       ! Theta boundary test.
       if ((ph%c(3) >= grid_d%n_theta) .or. (ph%c(3) < 1)) then
          ph%p_flag = -4
          exit
       end if

    end do

    ! If the loop ended geometrically before the sampled tau was reached,
    ! classify by the precomputed geometric endpoint.
    if ((ph%p_flag == 0) .and. (ph%tau < ph%tau_p) .and. (d >= smax - eps_path)) then
       if (hit_surface .eqv. .True.) then
          ph%p_flag = 1
       else
          ph%p_flag = 2
       end if
    end if

  end subroutine tauint_sph_3D

end module mc_k_tauint
