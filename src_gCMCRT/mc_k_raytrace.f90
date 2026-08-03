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

  ! Spherical ray-trace status contract.  Zero is the only successful status.
  ! Positive values are expected physical terminations; negative values are
  ! numerical/geometry errors.
  integer, parameter :: RAYTRACE_SUCCESS                    = 0
  integer, parameter :: RAYTRACE_OCCULTED_CORE              = 10
  integer, parameter :: RAYTRACE_ERR_NO_FORWARD_PATH        = -101
  integer, parameter :: RAYTRACE_ERR_INVALID_START_CELL     = -102
  integer, parameter :: RAYTRACE_ERR_INVALID_WALL_DISTANCE  = -103
  integer, parameter :: RAYTRACE_ERR_CROSSED_INNER_BOUNDARY = -104
  integer, parameter :: RAYTRACE_ERR_THETA_BOUNDARY         = -105
  integer, parameter :: RAYTRACE_ERR_INVALID_OPACITY        = -106
  integer, parameter :: RAYTRACE_ERR_INVALID_TAU            = -107
  integer, parameter :: RAYTRACE_ERR_INVALID_CONTRIBUTION   = -108
  integer, parameter :: RAYTRACE_ERR_INVALID_GEOMETRY       = -109

  integer, parameter :: N_RAYTRACE_DIAGNOSTICS = 11

  ! Shared diagnostics for direct transmission rays.  The first failure of
  ! each type retains a representative packet snapshot for host reporting.
  integer, device :: raytrace_diag_count_d(N_RAYTRACE_DIAGNOSTICS)
  integer(8), device :: raytrace_diag_packet_d(N_RAYTRACE_DIAGNOSTICS)
  integer, device :: raytrace_diag_cell_d(3,N_RAYTRACE_DIAGNOSTICS)
  real(dp), device :: raytrace_diag_tau_d(N_RAYTRACE_DIAGNOSTICS)
  real(dp), device :: raytrace_diag_value_d(N_RAYTRACE_DIAGNOSTICS)

contains

  attributes(device) subroutine raytrace_sph_3D(ray)

    implicit none

    type(pac), intent(inout) :: ray

    integer, dimension(3) :: ioffset
    real(dp) :: taucell
    real(dp) :: d, dcell, dstep, dleft, dmove
    real(dp) :: r2, rcosa, r2min, smax, n2
    real(dp) :: deps, kappa
    real(dp) :: r_surf2, eps_path

    ! Begin the tau integration
    ray%tau    = 0.0_dp
    ray%d      = 0.0_dp
    ray%p_flag = RAYTRACE_SUCCESS
    d          = 0.0_dp

    ! Ray direction is now the observation direction
    ray%nxp = im_d%obsx
    ray%nyp = im_d%obsy
    ray%nzp = im_d%obsz

    n2 = ray%nxp**2 + ray%nyp**2 + ray%nzp**2
    if (ieee_is_nan(n2) .or. abs(n2) > huge(n2) .or. &
        abs(n2 - 1.0_dp) > 1.0e-10_dp) then
       ray%d = n2
       ray%p_flag = RAYTRACE_ERR_INVALID_GEOMETRY
       return
    end if

    ! Maximum geometric distance to the outer boundary, unless the line of
    ! sight intersects the inner spherical boundary first.
    r2     = ray%xp**2 + ray%yp**2 + ray%zp**2
    rcosa  = ray%xp*ray%nxp + ray%yp*ray%nyp + ray%zp*ray%nzp
    r2min  = max(r2 - rcosa**2, 0.0_dp)

    r_surf2 = grid_d%r_min * grid_d%r_min

    if (rcosa > 0.0_dp) then
       ! Ray points outward: only the outer boundary matters.
       smax = sqrt(max(grid_d%r2_max - r2min, 0.0_dp)) - rcosa
    else
       ! Ray initially points inward.
       if (r2min > r_surf2) then
          ! The chord misses the inner boundary and exits the outer boundary.
          smax = sqrt(max(grid_d%r2_max - r2min, 0.0_dp)) - rcosa
       else
          ! The chord hits the inner boundary: occulted.
          smax = -rcosa - sqrt(max(r_surf2 - r2min, 0.0_dp))
          ray%d = smax
          ray%p_flag = RAYTRACE_OCCULTED_CORE
          return
       end if
    end if

    if (ieee_is_nan(smax) .or. abs(smax) > huge(smax)) then
       ray%d = smax
       ray%p_flag = RAYTRACE_ERR_INVALID_GEOMETRY
       return
    end if

    eps_path = 1.0e-12_dp * max(1.0_dp, abs(smax))

    if (smax <= eps_path) then
       ray%d = smax
       ray%p_flag = RAYTRACE_ERR_NO_FORWARD_PATH
       return
    end if

    do while (ray%p_flag == RAYTRACE_SUCCESS .and. d < smax - eps_path)

       ! Defensive guard before accessing opacity arrays.
       if ((ray%c(1) < 1) .or. (ray%c(1) >= grid_d%n_lev) .or. &
           (ray%c(2) < 1) .or. (ray%c(2) >= grid_d%n_phi) .or. &
           (ray%c(3) < 1) .or. (ray%c(3) >= grid_d%n_theta)) then
          ray%d = d
          ray%p_flag = RAYTRACE_ERR_INVALID_START_CELL
          return
       end if
       if ((ray%ig < 1) .or. (ray%ig > ng_d)) then
          ray%d = real(ray%ig,dp)
          ray%p_flag = RAYTRACE_ERR_INVALID_OPACITY
          return
       end if

       call findwall_sph(ray, ioffset, dcell)

       if (dcell <= 0.0_dp .or. ieee_is_nan(dcell) .or. &
           abs(dcell) > huge(dcell)) then
          ray%d = dcell
          ray%p_flag = RAYTRACE_ERR_INVALID_WALL_DISTANCE
          return
       end if

       dleft = smax - d

       ! If the next cell wall is beyond the final geometric boundary,
       ! integrate only the terminal physical segment and escape.
       if (dcell >= dleft) then
          dstep = dleft

          kappa = rhokap_d(ray%ig, ray%c(1), ray%c(2), ray%c(3))
          if (kappa < 0.0_dp .or. ieee_is_nan(kappa) .or. &
              abs(kappa) > huge(kappa)) then
             ray%d = kappa
             ray%p_flag = RAYTRACE_ERR_INVALID_OPACITY
             return
          end if
          taucell = dstep * kappa
          if (taucell < 0.0_dp .or. ieee_is_nan(taucell) .or. &
              abs(taucell) > huge(taucell)) then
             ray%d = taucell
             ray%p_flag = RAYTRACE_ERR_INVALID_TAU
             return
          end if

          ray%xp = ray%xp + dstep * ray%nxp
          ray%yp = ray%yp + dstep * ray%nyp
          ray%zp = ray%zp + dstep * ray%nzp

          ray%tau = ray%tau + taucell
          if (ray%tau < 0.0_dp .or. ieee_is_nan(ray%tau) .or. &
              abs(ray%tau) > huge(ray%tau)) then
             ray%d = ray%tau
             ray%p_flag = RAYTRACE_ERR_INVALID_TAU
             return
          end if
          d = smax
          ray%d = d

          ! Successful, visible peel-off ray. Leave p_flag = 0.
          exit
       end if

       ! Normal cell-wall crossing.
       dstep = dcell

       kappa = rhokap_d(ray%ig, ray%c(1), ray%c(2), ray%c(3))
       if (kappa < 0.0_dp .or. ieee_is_nan(kappa) .or. &
           abs(kappa) > huge(kappa)) then
          ray%d = kappa
          ray%p_flag = RAYTRACE_ERR_INVALID_OPACITY
          return
       end if
       taucell = dstep * kappa
       if (taucell < 0.0_dp .or. ieee_is_nan(taucell) .or. &
           abs(taucell) > huge(taucell)) then
          ray%d = taucell
          ray%p_flag = RAYTRACE_ERR_INVALID_TAU
          return
       end if

       ! Numerical nudge only for position/cell crossing, not for optical-depth
       ! or geometric path-length accounting.
       deps = (r_d(ray%c(1)+1) - r_d(ray%c(1))) * 1.0e-12_dp
       deps = max(deps, 1.0e-12_dp)
       dmove = dstep + deps

       ray%xp = ray%xp + dmove * ray%nxp
       ray%yp = ray%yp + dmove * ray%nyp
       ray%zp = ray%zp + dmove * ray%nzp

       ray%c(1) = ray%c(1) + ioffset(1)
       ray%c(2) = ray%c(2) + ioffset(2)
       ray%c(3) = ray%c(3) + ioffset(3)

       ray%tau = ray%tau + taucell
       if (ray%tau < 0.0_dp .or. ieee_is_nan(ray%tau) .or. &
           abs(ray%tau) > huge(ray%tau)) then
          ray%d = ray%tau
          ray%p_flag = RAYTRACE_ERR_INVALID_TAU
          return
       end if
       d = d + dstep
       ray%d = d

       ! Entered inner radial boundary.
       if (ray%c(1) < 1) then
          ray%p_flag = RAYTRACE_ERR_CROSSED_INNER_BOUNDARY
          exit
       end if

       ! Escaped outer radial boundary.
       if (ray%c(1) >= grid_d%n_lev) then
          ray%p_flag = RAYTRACE_SUCCESS
          exit
       end if

       ! Periodic longitude wrapping.
       if (ray%c(2) >= grid_d%n_phi) then
          ray%c(2) = 1
       else if (ray%c(2) < 1) then
          ray%c(2) = grid_d%n_phi - 1
       end if

       ! Out of theta bounds.
       if ((ray%c(3) >= grid_d%n_theta) .or. (ray%c(3) < 1)) then
          ray%p_flag = RAYTRACE_ERR_THETA_BOUNDARY
          exit
       end if

    end do

  end subroutine raytrace_sph_3D

  attributes(device) subroutine checked_transmission_contribution(ray, contri, valid)
    implicit none

    type(pac), intent(in) :: ray
    real(dp), intent(out) :: contri
    logical, intent(out) :: valid

    contri = 0.0_dp
    valid = .False.

    if (ray%p_flag /= RAYTRACE_SUCCESS) then
       call record_raytrace_diagnostic(ray, ray%p_flag, ray%d)
       return
    end if

    if (ray%tau < 0.0_dp .or. ieee_is_nan(ray%tau) .or. &
        abs(ray%tau) > huge(ray%tau)) then
       call record_raytrace_diagnostic(ray, RAYTRACE_ERR_INVALID_TAU, ray%tau)
       return
    end if

    contri = ray%wght * (1.0_dp - exp(-ray%tau))
    if (contri < 0.0_dp .or. ieee_is_nan(contri) .or. &
        abs(contri) > huge(contri)) then
       call record_raytrace_diagnostic(ray, RAYTRACE_ERR_INVALID_CONTRIBUTION, contri)
       contri = 0.0_dp
       return
    end if

    valid = .True.

  end subroutine checked_transmission_contribution

  attributes(device) subroutine record_raytrace_diagnostic(ray, status, value)
    implicit none

    type(pac), intent(in) :: ray
    integer, intent(in) :: status
    real(dp), intent(in) :: value

    integer :: slot, previous_count

    select case(status)
    case(RAYTRACE_OCCULTED_CORE)
       slot = 1
    case(RAYTRACE_ERR_NO_FORWARD_PATH)
       slot = 2
    case(RAYTRACE_ERR_INVALID_START_CELL)
       slot = 3
    case(RAYTRACE_ERR_INVALID_WALL_DISTANCE)
       slot = 4
    case(RAYTRACE_ERR_CROSSED_INNER_BOUNDARY)
       slot = 5
    case(RAYTRACE_ERR_THETA_BOUNDARY)
       slot = 6
    case(RAYTRACE_ERR_INVALID_OPACITY)
       slot = 7
    case(RAYTRACE_ERR_INVALID_TAU)
       slot = 8
    case(RAYTRACE_ERR_INVALID_CONTRIBUTION)
       slot = 9
    case(RAYTRACE_ERR_INVALID_GEOMETRY)
       slot = 10
    case default
       slot = 11
    end select

    previous_count = atomicadd(raytrace_diag_count_d(slot), 1)
    if (previous_count == 0) then
       raytrace_diag_packet_d(slot) = ray%id
       raytrace_diag_cell_d(1,slot) = ray%c(1)
       raytrace_diag_cell_d(2,slot) = ray%c(2)
       raytrace_diag_cell_d(3,slot) = ray%c(3)
       raytrace_diag_tau_d(slot) = ray%tau
       raytrace_diag_value_d(slot) = value
    end if

  end subroutine record_raytrace_diagnostic

  subroutine reset_raytrace_diagnostics()
    implicit none

    integer :: count_zero(N_RAYTRACE_DIAGNOSTICS)
    integer(8) :: packet_zero(N_RAYTRACE_DIAGNOSTICS)
    integer :: cell_zero(3,N_RAYTRACE_DIAGNOSTICS)
    real(dp) :: real_zero(N_RAYTRACE_DIAGNOSTICS)

    count_zero(:) = 0
    packet_zero(:) = 0_8
    cell_zero(:,:) = 0
    real_zero(:) = 0.0_dp

    raytrace_diag_count_d(:) = count_zero(:)
    raytrace_diag_packet_d(:) = packet_zero(:)
    raytrace_diag_cell_d(:,:) = cell_zero(:,:)
    raytrace_diag_tau_d(:) = real_zero(:)
    raytrace_diag_value_d(:) = real_zero(:)

  end subroutine reset_raytrace_diagnostics

  subroutine report_raytrace_diagnostics(context, phase_index, wavelength_index, &
      wavelength, packet_count)
    implicit none

    character(len=*), intent(in) :: context
    integer, intent(in) :: phase_index, wavelength_index, packet_count
    real(dp), intent(in) :: wavelength

    integer :: counts(N_RAYTRACE_DIAGNOSTICS)
    integer(8) :: packets(N_RAYTRACE_DIAGNOSTICS)
    integer :: cells(3,N_RAYTRACE_DIAGNOSTICS)
    real(dp) :: taus(N_RAYTRACE_DIAGNOSTICS)
    real(dp) :: values(N_RAYTRACE_DIAGNOSTICS)
    integer :: i, total

    counts(:) = raytrace_diag_count_d(:)
    total = sum(counts)
    if (total == 0) return

    packets(:) = raytrace_diag_packet_d(:)
    cells(:,:) = raytrace_diag_cell_d(:,:)
    taus(:) = raytrace_diag_tau_d(:)
    values(:) = raytrace_diag_value_d(:)

    print*, 'WARNING: invalid direct transmission rays were skipped.'
    print*, '  estimator: ', trim(context)
    print*, '  phase, wavelength index, wavelength: ', phase_index, wavelength_index, wavelength
    print*, '  requested packets, skipped invalid rays: ', packet_count, total
    if (packet_count > 0) then
       print*, '  skipped fraction: ', real(total,dp)/real(packet_count,dp)
    end if
    print*, '  The estimator retains the requested-packet normalisation.'

    do i = 1, N_RAYTRACE_DIAGNOSTICS
       if (counts(i) == 0) cycle
       print*, '  status code, count: ', raytrace_diagnostic_code(i), counts(i)
       print*, '    meaning: ', trim(raytrace_diagnostic_description(i))
       print*, '    first packet, cell: ', packets(i), cells(:,i)
       print*, '    first tau: ', taus(i)
       print*, '    ', trim(raytrace_diagnostic_value_label(i)), ': ', values(i)
    end do

  end subroutine report_raytrace_diagnostics

  integer function raytrace_diagnostic_code(slot)
    implicit none
    integer, intent(in) :: slot

    select case(slot)
    case(1)
       raytrace_diagnostic_code = RAYTRACE_OCCULTED_CORE
    case(2)
       raytrace_diagnostic_code = RAYTRACE_ERR_NO_FORWARD_PATH
    case(3)
       raytrace_diagnostic_code = RAYTRACE_ERR_INVALID_START_CELL
    case(4)
       raytrace_diagnostic_code = RAYTRACE_ERR_INVALID_WALL_DISTANCE
    case(5)
       raytrace_diagnostic_code = RAYTRACE_ERR_CROSSED_INNER_BOUNDARY
    case(6)
       raytrace_diagnostic_code = RAYTRACE_ERR_THETA_BOUNDARY
    case(7)
       raytrace_diagnostic_code = RAYTRACE_ERR_INVALID_OPACITY
    case(8)
       raytrace_diagnostic_code = RAYTRACE_ERR_INVALID_TAU
    case(9)
       raytrace_diagnostic_code = RAYTRACE_ERR_INVALID_CONTRIBUTION
    case(10)
       raytrace_diagnostic_code = RAYTRACE_ERR_INVALID_GEOMETRY
    case default
       raytrace_diagnostic_code = -999
    end select

  end function raytrace_diagnostic_code

  character(len=80) function raytrace_diagnostic_description(slot)
    implicit none
    integer, intent(in) :: slot

    select case(slot)
    case(1)
       raytrace_diagnostic_description = 'ray intersected the opaque inner boundary'
    case(2)
       raytrace_diagnostic_description = 'no positive forward path to the outer boundary'
    case(3)
       raytrace_diagnostic_description = 'ray started from an invalid spherical cell'
    case(4)
       raytrace_diagnostic_description = 'wall finder returned a non-positive or non-finite distance'
    case(5)
       raytrace_diagnostic_description = 'annulus ray unexpectedly crossed the inner boundary'
    case(6)
       raytrace_diagnostic_description = 'ray unexpectedly left the theta grid'
    case(7)
       raytrace_diagnostic_description = 'g index or sampled opacity was invalid'
    case(8)
       raytrace_diagnostic_description = 'optical-depth increment or accumulated tau was invalid'
    case(9)
       raytrace_diagnostic_description = 'blocked-light contribution was negative or non-finite'
    case(10)
       raytrace_diagnostic_description = 'ray position/direction produced invalid spherical geometry'
    case default
       raytrace_diagnostic_description = 'unknown non-zero ray-trace status'
    end select

  end function raytrace_diagnostic_description

  character(len=48) function raytrace_diagnostic_value_label(slot)
    implicit none
    integer, intent(in) :: slot

    select case(slot)
    case(1)
       raytrace_diagnostic_value_label = 'path length to opaque core'
    case(2)
       raytrace_diagnostic_value_label = 'computed forward path length'
    case(3)
       raytrace_diagnostic_value_label = 'path length before invalid cell'
    case(4)
       raytrace_diagnostic_value_label = 'invalid wall distance'
    case(5:6)
       raytrace_diagnostic_value_label = 'path length at invalid crossing'
    case(7)
       raytrace_diagnostic_value_label = 'invalid g index or opacity'
    case(8)
       raytrace_diagnostic_value_label = 'invalid optical-depth value'
    case(9)
       raytrace_diagnostic_value_label = 'invalid blocked-light contribution'
    case(10)
       raytrace_diagnostic_value_label = 'invalid geometry value'
    case default
       raytrace_diagnostic_value_label = 'status-dependent diagnostic value'
    end select

  end function raytrace_diagnostic_value_label

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
