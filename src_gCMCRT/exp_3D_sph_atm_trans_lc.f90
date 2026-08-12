module exp_3D_sph_atm_trans_lc_kernel
  use mc_precision
  use mc_class_pac
  use mc_class_grid
  use mc_k_source_pac_inc
  use mc_k_findcell
  use mc_k_gord_samp
  use mc_k_raytrace, only : raytrace_sph_3D, checked_transmission_contribution, &
    & record_raytrace_diagnostic, RAYTRACE_ERR_INVALID_START_CELL, &
    & RAYTRACE_ERR_INVALID_GEOMETRY
  use cudafor
  use curand_device
  implicit none

  type(curandStateMRG32k3a), allocatable, dimension(:), device :: iseed
  integer, device :: transit_tangent_chord_d

contains

  attributes(global) subroutine set_iseed(Nph)
    implicit none

    integer, intent(in) :: Nph
    integer(8) :: id, seed, seq, offset

    id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (id > Nph) return

    seed = random_seed_d
    seq = id - 1
    offset = 0
    call curand_init(seed, seq, offset, iseed(id))

  end subroutine set_iseed


  attributes(global) subroutine exp_3D_sph_atm_trans_lc_k(l, Nph)
    implicit none

    integer, intent(in) :: l, Nph
    type(pac) :: ph, ray
    integer :: count_stat
    real(dp) :: contri, rstat, lon
    real(dp) :: t_ca, xca, yca
    logical :: valid_direct_ray

    ph%id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (ph%id > Nph) return

    ph%iseed = iseed(ph%id)
    ph%wght = 1.0_dp
    ph%geo = 2
    ph%wl = wl_d(l)

    ! Sample the phase-appropriate projected overlap: the full planet disk in
    ! total transit or the star/planet lens in ingress/egress.  The source also
    ! applies the limb-darkening weight.
    call source_pac_inc_3D_transit(ph)
    if (ph%p_flag == TRANSIT_SOURCE_TANGENT_CHORD) then
      count_stat = atomicadd(transit_tangent_chord_d, 1)
      iseed(ph%id) = ph%iseed
      return
    else if (ph%p_flag == TRANSIT_SOURCE_INVALID) then
      ray = ph
      ray%c(:) = 0
      ray%tau = 0.0_dp
      ray%d = ph%bp
      ray%p_flag = RAYTRACE_ERR_INVALID_GEOMETRY
      call record_raytrace_diagnostic(ray, ray%p_flag, ray%d)
      iseed(ph%id) = ph%iseed
      return
    end if
    if (ph%wght <= 0.0_dp) then
      iseed(ph%id) = ph%iseed
      return
    end if

    ! --- Opaque planet body (impact parameter below the solid radius) -------
    ! The solid-body occultation is handled deterministically on the host with
    ! a batman-style limb-darkened overlap integral.  Do not add opaque-core
    ! packets to the Monte Carlo accumulators.
    if (ph%bp <= grid_d%r_min) then
      iseed(ph%id) = ph%iseed
      return
    end if

    ! --- Atmospheric annulus: trace the slant optical depth -----------------
    call locate(r_d, ph%bp, ph%b_idx)
    if ((ph%b_idx < 1) .or. (ph%b_idx > grid_d%n_lay)) then
      ray = ph
      ray%c(:) = 0
      ray%tau = 0.0_dp
      ray%d = ph%bp
      ray%p_flag = RAYTRACE_ERR_INVALID_GEOMETRY
      call record_raytrace_diagnostic(ray, ray%p_flag, ray%d)
      iseed(ph%id) = ph%iseed
      return
    end if

    call findcell(ph)
    if ((ph%c(1) < 1) .or. (ph%c(1) >= grid_d%n_lev) .or. &
        (ph%c(2) < 1) .or. (ph%c(2) >= grid_d%n_phi) .or. &
        (ph%c(3) < 1) .or. (ph%c(3) >= grid_d%n_theta)) then
      ray = ph
      ray%tau = 0.0_dp
      ray%d = 0.0_dp
      ray%p_flag = RAYTRACE_ERR_INVALID_START_CELL
      call record_raytrace_diagnostic(ray, ray%p_flag, ray%d)
      iseed(ph%id) = ph%iseed
      return
    end if

    if (do_g_bias_d .eqv. .True.) then
      call gord_samp_bias(ph)
    else
      call gord_samp(ph)
    end if

    ray = ph
    call raytrace_sph_3D(ray)
    call checked_transmission_contribution(ray, contri, valid_direct_ray)
    if (valid_direct_ray .eqv. .False.) then
      iseed(ph%id) = ph%iseed
      return
    end if

    ! Atmosphere-only blocked flux.  The opaque body is added on the host.
    if (use_block_accum_d .eqv. .True.) then
      rstat = atomicadd(block_accum_d(blockIdx%x,BLOCK_ACC_TOTAL),contri)
    else
      rstat = atomicadd(A_atm_d,contri)
    end if

    ! East/west limb split from the tangent (closest-approach) point of the
    ! chord in the planet frame, NOT the upstream launch position.  The ray
    ! starts at ph and travels along ph%n; the closest approach to the planet
    ! centre is at parameter t_ca = -(ph . n).
    t_ca = -(ph%xp*ph%nxp + ph%yp*ph%nyp + ph%zp*ph%nzp)
    xca = ph%xp + t_ca*ph%nxp
    yca = ph%yp + t_ca*ph%nyp
    lon = atan2(yca, xca)
    if (sin(lon) >= 0.0_dp) then
      if (use_block_accum_d .eqv. .True.) then
        rstat = atomicadd(block_accum_d(blockIdx%x,BLOCK_ACC_EAST),contri)
      else
        rstat = atomicadd(A_atm_east_d,contri)
      end if
    else
      if (use_block_accum_d .eqv. .True.) then
        rstat = atomicadd(block_accum_d(blockIdx%x,BLOCK_ACC_WEST),contri)
      else
        rstat = atomicadd(A_atm_west_d,contri)
      end if
    end if

    iseed(ph%id) = ph%iseed

  end subroutine exp_3D_sph_atm_trans_lc_k

end module exp_3D_sph_atm_trans_lc_kernel


subroutine exp_3D_sph_atm_trans_lightcurve()
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_class_imag
  use exp_3D_sph_atm_trans_lc_kernel
  use mc_k_raytrace, only : reset_raytrace_diagnostics, report_raytrace_diagnostics
  use mc_opacset
  use mc_read_prf
  use mc_views
  use LHS_sampling_mod, only : LHS_sample
  use random_cpu
  use cudafor
  implicit none

  character(len=8) :: fmt
  character(len=3) :: n_str
  integer :: Nph, l, iscat, istat, n, s_wl, n_trans, n_ingress, n_loop
  integer :: n_wl_out
  integer :: nml_iostat
  integer :: transit_tangent_chord
  integer, allocatable, dimension(:) :: uT
  integer, device :: l_d, Nph_d
  integer :: n_theta, n_phi, n_lay
  real(dp) :: pl, pc, sc
  real(dp) :: LD_norm, norm0, norm_phase, area_phase, R_top, R_s_cm
  real(dp) :: geom_tol
  real(dp) :: rp_body, d_star
  real(dp) :: A_lens, lens_psi, lens_xlo, lens_xhi, lens_ymax
  real(dp) :: A_atm, A_atm_east, A_atm_west
  real(dp) :: depth, depth_atm, depth_atm_east, depth_atm_west, depth_opaque
  real(dp), allocatable, dimension(:,:) :: trans_block_accum
  logical :: use_block_accum
  character(len=512) :: nml_iomsg
  type(dim3) :: blocks, threads

  namelist /sph_3D_trans_lc/ Nph, s_wl, n_wl, n_trans, n_ingress, pl, pc, sc, &
    & n_theta, n_phi, n_lay, iscat, use_block_accum

  use_block_accum = .True.
  read(u_nml, nml=sph_3D_trans_lc, iostat=nml_iostat, iomsg=nml_iomsg)
  if (nml_iostat /= 0) then
    print*, 'ERROR reading &sph_3D_trans_lc from CMCRT.nml: ', trim(nml_iomsg)
    stop
  end if

  if (Nph < 1) then
    print*, 'ERROR: Nph must be positive in 3D_sph_trans_lc.'
    stop
  end if

  fmt = '(I3.3)'

  grid%n_lay = n_lay
  grid%n_lev = n_lay + 1
  grid%n_theta = n_theta
  grid%n_phi = n_phi

  pl_d = pl
  pc_d = pc
  sc_d = sc
  iscat_d = iscat

  threads = dim3(128, 1, 1)
  blocks = dim3(Nph / threads%x, 1, 1)
  if (mod(Nph,threads%x) /= 0) blocks%x = blocks%x + 1
  allocate(iseed(Nph))
  use_block_accum_d = use_block_accum
  if (use_block_accum .eqv. .True.) then
    allocate(trans_block_accum(blocks%x,N_BLOCK_ACC))
    allocate(block_accum_d(blocks%x,N_BLOCK_ACC))
  end if
  Nph_d = Nph
  call set_iseed<<<blocks, threads>>>(Nph_d)

  print*, 'Transit light-curve per-block accumulators: ', use_block_accum

  call read_1D_prf()
  call read_wl(s_wl)
  n_wl_out = n_wl - s_wl + 1
  if (ck .eqv. .True.) then
    call read_g_ord()
  else
    ng = 1
    ng_d = 1
  end if

  call set_grid()
  call set_transit_views(n_trans, n_ingress, n_loop)

  grid_d = grid

  allocate(uT(n_loop))

  ! Stellar normalisation: F_star = 2*pi*R_s^2 * LD_norm, with
  ! LD_norm = int_0^1 I(mu) mu dmu for the active limb-darkening law.
  ! A sampled region of area A_phase, sampled by Nph uniform packets giving the
  ! atmosphere-only MC sum S = sum[ I(mu)*blocking ], contributes
  !   depth_atm = (A_phase * S/Nph) / (2*pi*R_s^2*LD_norm)
  !             = norm0 * A_phase * S,
  ! with norm0 = 1/(2*pi*R_s^2*LD_norm*Nph).  A_phase = pi*R_top^2 at full
  ! transit (the whole planet disk) or the star/planet lens area at
  ! ingress/egress (see transit_lens_geom / source_pac_inc_3D_transit).
  LD_norm = stellar_LD_norm()
  R_s_cm = Rs * Rsun
  R_top = H(grid%n_lev)
  geom_tol = transit_geometry_tolerance(R_top, R_s_cm, sm_ax*Au)
  rp_body = H(1) / R_s_cm
  norm0 = 1.0_dp / (2.0_dp * pi * R_s_cm**2 * LD_norm * real(Nph, dp))

  ! Resolve any roundoff-level contact phases before writing the phase state to
  ! the output headers.  The same tolerance is used by set_transit_views.
  do n = 1, n_loop
    if (trans_state(n) == 1) then
      call transit_lens_geom(R_top, R_s_cm, xstar(n), ystar(n), geom_tol, &
        & A_lens, lens_psi, lens_xlo, lens_xhi, lens_ymax)
      if (A_lens <= 0.0_dp) trans_state(n) = 0
    end if
  end do

  do n = 1, n_loop
    write(n_str, fmt) n
    open(newunit=uT(n), file='Transit_'//trim(n_str)//'.txt', status='replace', action='write')
    ! header: output wavelength count, radii, view angles, phase, state, star pos
    write(uT(n),*) n_wl_out, H(1), H(grid%n_lev), viewphi_n(n), viewthet_n(n), &
      & phi_key(n), trans_state(n), xstar(n), ystar(n), zstar(n)
    ! column key: wl, depth, depth_atm, depth_atm_east, depth_atm_west, depth_opaque
    call flush(uT(n))
  end do

  if (LHS .eqv. .True.) then
    allocate(x_ran(Nph), y_ran(Nph), z_ran(Nph), x_ran_d(Nph), y_ran_d(Nph), z_ran_d(Nph))
    call rng_seed(random_seed)
  end if

  call read_next_opac(s_wl)

  do l = s_wl, n_wl
    call set_grid_opac(iscat)

    do n = 1, n_loop
      im%vphi = viewphi_n(n)
      im%vtheta = viewthet_n(n)
      call set_image()

      im%fsum = 0.0_dp
      im%qsum = 0.0_dp
      im%usum = 0.0_dp
      im%fail_pscat = 0
      im%fail_pemit = 0

      ! Per-phase blocked-region area, and (partial phases only) the lens
      ! sampler geometry consumed by source_pac_inc_3D_transit.
      if (trans_state(n) == 1) then
        call transit_lens_geom(R_top, R_s_cm, xstar(n), ystar(n), geom_tol, &
          & A_lens, lens_psi, lens_xlo, lens_xhi, lens_ymax)
        trans_delta_d = sqrt(xstar(n)**2 + ystar(n)**2)
        trans_psi_d  = lens_psi
        trans_xlo_d  = lens_xlo
        trans_xhi_d  = lens_xhi
        trans_ymax_d = lens_ymax
        area_phase = A_lens
      else if (trans_state(n) == 2) then
        area_phase = pi * R_top**2
      else
        area_phase = 0.0_dp
      end if
      norm_phase = norm0 * area_phase
      trans_state_d = trans_state(n)
      xstar_d = xstar(n)
      ystar_d = ystar(n)

      A_atm_d = 0.0_dp
      A_atm_east_d = 0.0_dp
      A_atm_west_d = 0.0_dp
      if (use_block_accum .eqv. .True.) then
        block_accum_d(:,:) = 0.0_dp
      end if
      transit_tangent_chord_d = 0

      l_d = l
      im_d = im
      call reset_raytrace_diagnostics()

      if (trans_state(n) /= 0) then
        if (LHS .eqv. .True.) then
          call LHS_sample(Nph, 2, x_ran, y_ran, z_ran, .False.)
          x_ran_d(:) = x_ran(:)
          y_ran_d(:) = y_ran(:)
          z_ran_d(:) = z_ran(:)
        end if

        call exp_3D_sph_atm_trans_lc_k<<<blocks, threads>>>(l_d, Nph_d)

        istat = cudaDeviceSynchronize()
        if (istat /= 0) then
          print*, 'ERROR after exp_3D_sph_atm_trans_lc:', istat
          stop
        end if
        call report_raytrace_diagnostics('3D_sph_trans_lc', n, l, wl(l), Nph)
      end if

      if (use_block_accum .eqv. .True.) then
        trans_block_accum(:,:) = block_accum_d(:,:)
        A_atm = sum(trans_block_accum(:,BLOCK_ACC_TOTAL))
        A_atm_east = sum(trans_block_accum(:,BLOCK_ACC_EAST))
        A_atm_west = sum(trans_block_accum(:,BLOCK_ACC_WEST))
      else
        A_atm = A_atm_d
        A_atm_east = A_atm_east_d
        A_atm_west = A_atm_west_d
      end if
      transit_tangent_chord = transit_tangent_chord_d
      if (transit_tangent_chord > 0) then
        print*, '3D_sph_trans_lc: unresolved outer-tangent chords omitted'
        print*, '  phase, wavelength index, count: ', n, l, transit_tangent_chord
        print*, '  These lie inside the launch guard shell and are assigned zero contribution.'
      end if

      if (trans_state(n) == 0) then
        depth_opaque = 0.0_dp
      else
        d_star = sqrt(xstar(n)**2 + ystar(n)**2) / R_s_cm
        depth_opaque = opaque_ld_depth_batman(rp_body, d_star, LD_norm)
      end if

      depth_atm = norm_phase * A_atm
      depth_atm_east = norm_phase * A_atm_east
      depth_atm_west = norm_phase * A_atm_west
      depth = depth_opaque + depth_atm

      write(uT(n),*) wl(l), depth, depth_atm, depth_atm_east, depth_atm_west, depth_opaque

      print*, n, l, wl(l), phi_key(n), trans_state(n), depth, depth_atm, depth_opaque
    end do

    call read_next_opac(l+1)
  end do

contains

  real(dp) function stellar_intensity_mu(mu) result(Imu)
    implicit none
    real(dp), intent(in) :: mu
    real(dp) :: mus

    mus = max(0.0_dp, min(1.0_dp, mu))
    if (do_LD .eqv. .False.) then
      Imu = 1.0_dp
      return
    end if

    select case(ilimb)
    case(1)
      Imu = 1.0_dp - LD_c(1)*(1.0_dp - mus)
    case(2)
      Imu = 1.0_dp - LD_c(1)*(1.0_dp - mus) - LD_c(2)*(1.0_dp - mus)**2
    case(3)
      Imu = 1.0_dp - LD_c(1)*(1.0_dp - mus) - LD_c(2)*(1.0_dp - sqrt(mus))
    case(4)
      if (mus > 0.0_dp) then
        Imu = 1.0_dp - LD_c(1)*(1.0_dp - mus) - LD_c(2)*mus*log(mus)
      else
        Imu = 1.0_dp - LD_c(1)
      end if
    case(5)
      ! Nonzero LD_c(2) is rejected at startup because this law is singular
      ! at the stellar limb.  The remaining expression is the linear limit.
      Imu = 1.0_dp - LD_c(1)*(1.0_dp - mus)
    case(6)
      Imu = 1.0_dp - LD_c(1)*(1.0_dp - mus) - LD_c(2)*(1.0_dp - mus**(1.5_dp)) &
        & - LD_c(3)*(1.0_dp - mus**2)
    case(7)
      Imu = 1.0_dp - LD_c(1)*(1.0_dp - sqrt(mus)) - LD_c(2)*(1.0_dp - mus) &
        & - LD_c(3)*(1.0_dp - mus**(1.5_dp)) - LD_c(4)*(1.0_dp - mus**2)
    case(8)
      Imu = 1.0_dp - LD_c(1)*(1.0_dp - mus**LD_c(2))
    case default
      Imu = 1.0_dp
    end select

  end function stellar_intensity_mu


  real(dp) function circle_overlap_area_norm(x, rp, d) result(area)
    implicit none
    real(dp), intent(in) :: x, rp, d
    real(dp) :: u, v, w, eps_geom

    eps_geom = 100.0_dp * epsilon(1.0_dp)

    if ((x <= 0.0_dp) .or. (rp <= 0.0_dp)) then
      area = 0.0_dp
    else if (d >= x + rp - eps_geom) then
      area = 0.0_dp
    else if (d <= abs(x - rp) + eps_geom) then
      area = pi * min(x, rp)**2
    else
      u = (d**2 + x**2 - rp**2) / (2.0_dp*d*x)
      v = (d**2 + rp**2 - x**2) / (2.0_dp*d*rp)
      w = (-d + x + rp)*(d + x - rp)*(d - x + rp)*(d + x + rp)
      area = x**2 * acos(min(1.0_dp, max(-1.0_dp, u))) &
        & + rp**2 * acos(min(1.0_dp, max(-1.0_dp, v))) &
        & - 0.5_dp * sqrt(max(0.0_dp, w))
    end if

  end function circle_overlap_area_norm


  real(dp) function opaque_ld_depth_batman(rp, d, LDfac) result(depth_opaque)
    implicit none
    real(dp), intent(in) :: rp, d, LDfac
    integer, parameter :: n_opaque_ld = 2000
    integer :: i
    real(dp) :: x0, x1, dx, xa, xb, xm, mu, area_a, area_b

    if ((rp <= 0.0_dp) .or. (LDfac <= 0.0_dp)) then
      depth_opaque = 0.0_dp
      return
    end if

    x0 = max(d - rp, 0.0_dp)
    x1 = min(d + rp, 1.0_dp)
    if (x1 <= x0) then
      depth_opaque = 0.0_dp
      return
    end if

    dx = (x1 - x0) / real(n_opaque_ld, dp)
    depth_opaque = 0.0_dp
    area_a = 0.0_dp
    do i = 1, n_opaque_ld
      xa = x0 + dx * real(i - 1, dp)
      xb = x0 + dx * real(i, dp)
      xm = 0.5_dp * (xa + xb)
      mu = sqrt(max(1.0_dp - xm**2, 0.0_dp))
      area_b = circle_overlap_area_norm(xb, rp, d)
      depth_opaque = depth_opaque + stellar_intensity_mu(mu) * (area_b - area_a)
      area_a = area_b
    end do

    depth_opaque = depth_opaque / (2.0_dp * pi * LDfac)

  end function opaque_ld_depth_batman


  !! Disk integral  LDfac = int_0^1 I(mu) mu dmu  for the active limb-darkening
  !! law (ilimb / LD_c from the &main namelist).  Stellar flux is
  !! F_star = 2*pi*R_s^2*LDfac.  Must match limb_darkening_mu in
  !! mc_k_source_pac_inc.f90.  For a uniform disk (do_LD = .False.) this is 1/2.
  real(dp) function stellar_LD_norm() result(LDfac)
    implicit none
    integer :: im_mu, n_mu
    real(dp) :: mu, dmu, Imu

    if (do_LD .eqv. .False.) then
      LDfac = 0.5_dp
      return
    end if

    n_mu = 2000
    dmu = 1.0_dp/real(n_mu, dp)
    LDfac = 0.0_dp
    do im_mu = 1, n_mu
      mu = (real(im_mu, dp) - 0.5_dp)*dmu     ! midpoint, avoids mu = 0
      Imu = stellar_intensity_mu(mu)
      LDfac = LDfac + Imu*mu*dmu
    end do

  end function stellar_LD_norm

end subroutine exp_3D_sph_atm_trans_lightcurve
