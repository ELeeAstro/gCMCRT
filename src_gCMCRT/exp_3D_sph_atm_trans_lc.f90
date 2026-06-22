module exp_3D_sph_atm_trans_lc_kernel
  use mc_precision
  use mc_class_pac
  use mc_class_grid
  use mc_k_source_pac_inc
  use mc_k_findcell
  use mc_k_gord_samp
  use mc_k_raytrace
  use cudafor
  use curand_device
  implicit none

  type(curandStateMRG32k3a), allocatable, dimension(:), device :: iseed

contains

  attributes(global) subroutine set_iseed(Nph)
    implicit none

    integer, intent(in) :: Nph
    integer(8) :: id, seed
    integer :: seq, offset

    id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (id > Nph) return

    seed = id + id**2 + id/2
    seq = 0
    offset = 0
    call curand_init(seed, seq, offset, iseed(id))

  end subroutine set_iseed


  attributes(global) subroutine exp_3D_sph_atm_trans_lc_k(l, Nph)
    implicit none

    integer, intent(in) :: l, Nph
    type(pac) :: ph, ray
    real(dp) :: contri, rstat, lon, blocking
    real(dp) :: t_ca, xca, yca

    ph%id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (ph%id > Nph) return

    ph%iseed = iseed(ph%id)
    ph%wght = 1.0_dp
    ph%geo = 2
    ph%wl = wl_d(l)

    ! Sample a point on the full projected planet disk.  source applies the
    ! finite-stellar-disk mask and limb-darkening weight (ph%wght = 0 if the
    ! point falls off the stellar disk or the planet is out of transit).
    call source_pac_inc_3D_transit(ph)
    if (ph%wght <= 0.0_dp) then
      iseed(ph%id) = ph%iseed
      return
    end if

    ! --- Opaque planet body (impact parameter below the solid radius) -------
    ! Blocks the stellar beam completely at every wavelength / g-ordinate, so
    ! no ray trace and no corr-k weight are needed.  ph%wght already carries the
    ! limb-darkening weight and on-star mask.
    if (ph%bp < grid_d%r_min) then
      rstat = atomicadd(A_block_d, ph%wght)
      iseed(ph%id) = ph%iseed
      return
    end if

    ! --- Atmospheric annulus: trace the slant optical depth -----------------
    call locate(r_d, ph%bp, ph%b_idx)
    if ((ph%b_idx < 1) .or. (ph%b_idx > grid_d%n_lay)) then
      iseed(ph%id) = ph%iseed
      return
    end if

    call findcell(ph)
    if ((ph%c(1) < 1) .or. (ph%c(1) >= grid_d%n_lev) .or. &
        (ph%c(2) < 1) .or. (ph%c(2) >= grid_d%n_phi) .or. &
        (ph%c(3) < 1) .or. (ph%c(3) >= grid_d%n_theta)) then
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
    if (ray%p_flag /= 0) then
      iseed(ph%id) = ph%iseed
      return
    end if

    blocking = 1.0_dp - exp(-ray%tau)
    contri = ph%wght * blocking

    ! Total blocked flux (body + atmosphere) and the atmosphere-only part.
    rstat = atomicadd(A_block_d, contri)
    rstat = atomicadd(A_atm_d, contri)

    ! East/west limb split from the tangent (closest-approach) point of the
    ! chord in the planet frame, NOT the upstream launch position.  The ray
    ! starts at ph and travels along ph%n; the closest approach to the planet
    ! centre is at parameter t_ca = -(ph . n).
    t_ca = -(ph%xp*ph%nxp + ph%yp*ph%nyp + ph%zp*ph%nzp)
    xca = ph%xp + t_ca*ph%nxp
    yca = ph%yp + t_ca*ph%nyp
    lon = atan2(yca, xca)
    if (sin(lon) >= 0.0_dp) then
      rstat = atomicadd(A_atm_east_d, contri)
    else
      rstat = atomicadd(A_atm_west_d, contri)
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
  integer, allocatable, dimension(:) :: uT
  integer, device :: l_d, Nph_d
  integer :: n_theta, n_phi, n_lay
  real(dp) :: pl, pc, sc
  real(dp) :: LD_norm, norm0, norm_phase, area_phase, R_top
  real(dp) :: A_lens, lens_psi, lens_xlo, lens_xhi, lens_ymax
  real(dp) :: A_block, A_atm, A_atm_east, A_atm_west
  real(dp) :: depth, depth_atm, depth_atm_east, depth_atm_west
  type(dim3) :: blocks, threads

  namelist /sph_3D_trans_lc/ Nph, s_wl, n_wl, n_trans, n_ingress, pl, pc, sc, &
    & n_theta, n_phi, n_lay, iscat

  read(u_nml, nml=sph_3D_trans_lc)

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
  blocks = dim3(ceiling(real(Nph,dp)/threads%x), 1, 1)
  allocate(iseed(Nph))
  Nph_d = Nph
  call set_iseed<<<blocks, threads>>>(Nph_d)

  call read_1D_prf()
  call read_wl()
  call read_g_ord()

  call set_grid()
  call set_transit_views(n_trans, n_ingress, n_loop)

  grid_d = grid

  allocate(uT(n_loop))

  ! Stellar normalisation: F_star = 2*pi*R_s^2 * LD_norm, with
  ! LD_norm = int_0^1 I(mu) mu dmu for the active limb-darkening law.
  ! A blocked region of area A_phase, sampled by Nph uniform packets giving the
  ! MC sum S = sum[ I(mu)*blocking ], contributes transit depth
  !   depth = (A_phase * S/Nph) / (2*pi*R_s^2*LD_norm) = norm0 * A_phase * S,
  ! with norm0 = 1/(2*pi*R_s^2*LD_norm*Nph).  A_phase = pi*R_top^2 at full
  ! transit (the whole planet disk) or the star/planet lens area at
  ! ingress/egress (see transit_lens_geom / source_pac_inc_3D_transit).
  LD_norm = stellar_LD_norm()
  R_top = H(grid%n_lev)
  norm0 = 1.0_dp / (2.0_dp * pi * (Rs*Rsun)**2 * LD_norm * real(Nph, dp))

  do n = 1, n_loop
    write(n_str, fmt) n
    open(newunit=uT(n), file='Transit_'//trim(n_str)//'.txt', action='readwrite')
    ! header: n_wl, inner/outer radius, view angles, phase, transit state, star pos
    write(uT(n),*) n_wl, H(1), H(grid%n_lev), viewphi_n(n), viewthet_n(n), &
      & phi_key(n), trans_state(n), xstar(n), ystar(n), zstar(n)
    ! column key: wl, depth, depth_atm, depth_atm_east, depth_atm_west
    call flush(uT(n))
  end do

  if (LHS .eqv. .True.) then
    allocate(x_ran(Nph), y_ran(Nph), z_ran(Nph), x_ran_d(Nph), y_ran_d(Nph), z_ran_d(Nph))
    call rng_seed(123)
  end if

  call read_next_opac(s_wl)

  do l = s_wl, n_wl
    call set_grid_opac()

    do n = 1, n_loop
      im%vphi = viewphi_n(n)
      im%vtheta = viewthet_n(n)
      call set_image()

      im%fsum = 0.0_dp
      im%qsum = 0.0_dp
      im%usum = 0.0_dp
      im%fail_pscat = 0
      im%fail_pemit = 0

      trans_state_d = trans_state(n)
      xstar_d = xstar(n)
      ystar_d = ystar(n)

      ! Per-phase blocked-region area, and (partial phases only) the lens
      ! sampler geometry consumed by source_pac_inc_3D_transit.
      if (trans_state(n) == 1) then
        call transit_lens_geom(R_top, Rs*Rsun, xstar(n), ystar(n), &
          & A_lens, lens_psi, lens_xlo, lens_xhi, lens_ymax)
        trans_delta_d = sqrt(xstar(n)**2 + ystar(n)**2)
        trans_psi_d  = lens_psi
        trans_xlo_d  = lens_xlo
        trans_xhi_d  = lens_xhi
        trans_ymax_d = lens_ymax
        area_phase = A_lens
      else
        area_phase = pi * R_top**2
      end if
      norm_phase = norm0 * area_phase

      A_block_d = 0.0_dp
      A_atm_d = 0.0_dp
      A_atm_east_d = 0.0_dp
      A_atm_west_d = 0.0_dp

      l_d = l
      im_d = im

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

      A_block = A_block_d
      A_atm = A_atm_d
      A_atm_east = A_atm_east_d
      A_atm_west = A_atm_west_d

      depth = norm_phase * A_block
      depth_atm = norm_phase * A_atm
      depth_atm_east = norm_phase * A_atm_east
      depth_atm_west = norm_phase * A_atm_west

      write(uT(n),*) wl(l), depth, depth_atm, depth_atm_east, depth_atm_west

      print*, n, l, wl(l), phi_key(n), trans_state(n), depth, depth_atm
    end do

    call read_next_opac(l+1)
  end do

contains

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
      select case(ilimb)
      case(1)
        Imu = 1.0_dp - LD_c(1)*(1.0_dp - mu)
      case(2)
        Imu = 1.0_dp - LD_c(1)*(1.0_dp - mu) - LD_c(2)*(1.0_dp - mu)**2
      case(3)
        Imu = 1.0_dp - LD_c(1)*(1.0_dp - mu) - LD_c(2)*(1.0_dp - sqrt(mu))
      case(4)
        Imu = 1.0_dp - LD_c(1)*(1.0_dp - mu) - LD_c(2)*mu*log(mu)
      case(5)
        Imu = 1.0_dp - LD_c(1)*(1.0_dp - mu) - LD_c(2)/(1.0_dp - exp(mu))
      case(6)
        Imu = 1.0_dp - LD_c(1)*(1.0_dp - mu) - LD_c(2)*(1.0_dp - mu**(1.5_dp)) &
          & - LD_c(3)*(1.0_dp - mu**2)
      case(7)
        Imu = 1.0_dp - LD_c(1)*(1.0_dp - sqrt(mu)) - LD_c(2)*(1.0_dp - mu) &
          & - LD_c(3)*(1.0_dp - mu**(1.5_dp)) - LD_c(4)*(1.0_dp - mu**2)
      case(8)
        Imu = 1.0_dp - LD_c(1)*(1.0_dp - mu**LD_c(2))
      case default
        Imu = 1.0_dp
      end select
      LDfac = LDfac + Imu*mu*dmu
    end do

  end function stellar_LD_norm

end subroutine exp_3D_sph_atm_trans_lightcurve
