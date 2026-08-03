module exp_3D_sph_atm_transmission_kernel
  use mc_precision
  use mc_class_pac
  use mc_class_grid
  use mc_k_source_pac_inc
  use mc_k_tauint
  use mc_k_scatt
  use mc_k_peeloff_scatt
  use mc_k_findcell
  use mc_k_tau_samp
  use mc_k_RR
  use mc_k_gord_samp
  use mc_k_raytrace, only : raytrace_sph_3D, checked_transmission_contribution
  use cudafor
  use curand_device
  implicit none

  integer :: nscat_tot
  integer, device :: nscat_tot_d

  type(curandStateMRG32k3a), allocatable, dimension(:), device :: iseed


contains

  attributes(global) subroutine set_iseed(Nph)
    implicit none

    integer, intent(in) :: Nph
    integer(8) :: id, seed, seq, offset

    id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (id > Nph) then
      return
    end if

    seed = random_seed_d
    seq = id - 1
    offset = 0
    call curand_init(seed, seq, offset, iseed(id))

  end subroutine set_iseed


  attributes(global) subroutine exp_3D_sph_atm_transmission_k(l,Nph)
    implicit none

    integer, intent(in) :: l, Nph
    type(pac) :: ph, ray
    integer :: seq, offset, i, n, istat, nscat
    integer :: b_cf_idx
    real(dp) :: contri, rstat, lon, b_phys
    logical :: valid_direct_ray

    ! Set a random seed for this packet
    ph%id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (ph%id > Nph) then
      return
    end if
    ph%iseed = iseed(ph%id)


    ph%wght = 1.0_dp
    ph%geo = 2
    ph%wl = wl_d(l)

    do
      call source_pac_inc_3D(ph)
      call locate(r_d,ph%bp,ph%b_idx)
      if ((ph%b_idx < 1) .or. (ph%b_idx > grid_d%n_lay)) then
         !print*, ph%id, ph%bp, ph%b_idx
         ph%wght = 1.0_dp ! have to return weight due to possible limb darkening
        cycle
      else
        exit
      end if
    enddo

    call findcell(ph)


    ! Sample a g-ordinance value (for corr-k)
    if (do_g_bias_d .eqv. .True.) then
      call gord_samp_bias(ph)
    else
      call gord_samp(ph)
    end if


    ! Perform peeloff at starting location
    ray = ph
    call raytrace_sph_3D(ray)

    call checked_transmission_contribution(ray, contri, valid_direct_ray)

    ! Exclude an invalid direct observer ray from every direct estimator.  The
    ! incident packet can still have a valid scattering history below.
    if (valid_direct_ray .eqv. .True.) then
      if (use_block_accum_d .eqv. .True.) then
        rstat = atomicadd(block_accum_d(blockIdx%x,BLOCK_ACC_TOTAL),contri)
      else
        rstat = atomicadd(T_trans_d,contri)
      end if

      ! Longitude of the initial packet position
      lon = atan2(ph%yp, ph%xp)

      ! Add contributions to the east and west limbs
      if (sin(lon) >= 0.0_dp) then
        if (use_block_accum_d .eqv. .True.) then
          rstat = atomicadd(block_accum_d(blockIdx%x,BLOCK_ACC_EAST),contri)
        else
          rstat = atomicadd(T_trans_east_d,contri)
        end if
      else
        if (use_block_accum_d .eqv. .True.) then
          rstat = atomicadd(block_accum_d(blockIdx%x,BLOCK_ACC_WEST),contri)
        else
          rstat = atomicadd(T_trans_west_d,contri)
        end if
      end if

      ! Do the contibution function for the binned b
      if (do_cf_d .eqv. .True.) then
        b_phys = ph%bp * H_d(1)
        if (b_phys <= b_cf_grid_d(1)) then
          b_cf_idx = 1
        else if (b_phys >= b_cf_grid_d(size(b_cf_grid_d))) then
          b_cf_idx = size(b_cf_d)
        else
          call locate(b_cf_grid_d,b_phys,b_cf_idx)
        end if
        rstat = atomicadd(b_cf_d(b_cf_idx),contri)
        istat = atomicadd(b_n_cf_d(b_cf_idx),1)
      end if
    end if


    if (do_scat_loop_d .eqv. .True.) then
      ph%p_flag = 0
    else
      ph%p_flag = -1
    end if

    ! Begin scattering loop

    nscat = 0

    !! Enter scattering loop
    do while (ph%p_flag == 0)
      !! Sample a tau for the packet
      ph%tau_p = -log(curand_uniform(ph%iseed))
      !call tau_force_scatt(ph)
      if (ph%p_flag /= 0) then
        exit
      end if
      !! Move packet for sampled tau distance
      call tauint_sph_3D(ph)

      if (ph%p_flag /= 0) then
        !print*, ph%id, ph%wght, ph%p_flag, 'died'
        exit
      end if

      if (curand_uniform(ph%iseed) <= dorg_d(ph%c(1),ph%c(2),ph%c(3))) then
        ! Gas scattering - Rayleigh scattering
        ph%wght = ph%wght * ssa_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))
        ph%iscatt = 3
      else
        ! Cloud scattering - do given scattering phase function
        ph%wght = ph%wght * ssa_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))
        ph%iscatt = iscat_d
      end if

      call peeloff_scatt(ph)
      call scatt_pac_2(ph)
      call RR_test(ph)

      nscat =  nscat + 1

    end do

    ! Add number of scatterings to total
    if (use_block_accum_d .eqv. .True.) then
      istat = atomicadd(block_nscat_accum_d(blockIdx%x),nscat)
    else
      istat = atomicadd(nscat_tot_d,nscat)
    end if

    ! Give back iseed to saved device array for next iteration with this ph%id
    iseed(ph%id) = ph%iseed

  end subroutine exp_3D_sph_atm_transmission_k


end module exp_3D_sph_atm_transmission_kernel

subroutine exp_3D_sph_atm_transmission()
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_class_imag
  use exp_3D_sph_atm_transmission_kernel
  use mc_k_raytrace, only : reset_raytrace_diagnostics, report_raytrace_diagnostics
  use mc_opacset
  use mc_read_prf
  use LHS_sampling_mod, only : LHS_sample
  use random_cpu
  use mc_k_scatt, only : validate_volume_scattering
  use cudafor
  implicit none


  integer :: Nph, l, uT, iscat, ucf, nb_cf, i, s_wl, n_wl_out
  integer, device :: l_d, Nph_d
  integer :: n_theta, n_phi, n_lay
  real(dp) :: dH, norm
  real(dp) :: viewthet, viewphi
  real(dp) :: pl, pc, sc
  real(dp), allocatable, dimension(:,:) :: trans_block_accum
  integer, allocatable, dimension(:) :: block_nscat_accum
  logical :: use_block_accum
  character(len=512) :: nml_iomsg

  integer :: istat, nml_iostat
  type(dim3) :: blocks, threads

  namelist /sph_3D_trans/ Nph, s_wl, n_wl, pl, pc, sc, n_theta, n_phi, n_lay, &
    & viewthet, viewphi, iscat, nb_cf, use_block_accum

  use_block_accum = .True.
  read(u_nml, nml=sph_3D_trans, iostat=nml_iostat, iomsg=nml_iomsg)
  if (nml_iostat /= 0) then
    print*, 'ERROR reading &sph_3D_trans from CMCRT.nml: ', trim(nml_iomsg)
    stop
  end if

  if (Nph < 1) then
    print*, 'ERROR: Nph must be positive in 3D_sph_trans.'
    stop
  end if
  call validate_volume_scattering(iscat, '3D_sph_trans')

  ! Give namelist paramaters to equilvanet values inside gCMCRT
  grid%n_lay = n_lay
  grid%n_lev = n_lay + 1
  grid%n_theta = n_theta
  grid%n_phi = n_phi

  im%vtheta = viewthet
  im%vphi = viewphi

  pl_d = pl
  pc_d = pc
  sc_d = sc
  iscat_d = iscat

  threads = dim3(128, 1, 1)
  blocks = dim3(ceiling(real(Nph,dp)/threads%x),1,1)
  allocate(iseed(Nph))
  Nph_d = Nph
  call set_iseed<<<blocks, threads>>>(Nph_d)

  call read_1D_prf()
  print*, 'read wl'
  call read_wl(s_wl)
  n_wl_out = n_wl - s_wl + 1
  print*, 'read gord'
  if (ck .eqv. .True.) then
    call read_g_ord()
  else
    ng = 1
    ng_d = 1
  end if


  call set_grid()
  call set_image()


  ! Send data to GPU data containers
  im_d = im
  grid_d = grid

  allocate(T_trans(n_wl), T_trans_east(n_wl), T_trans_west(n_wl))

  ! Grid for GPU threads/blocks
  threads = dim3(128,1,1)
  blocks = dim3(ceiling(real(Nph)/threads%x),1,1)

  use_block_accum_d = use_block_accum
  if (use_block_accum .eqv. .True.) then
    allocate(trans_block_accum(blocks%x,N_BLOCK_ACC))
    allocate(block_accum_d(blocks%x,N_BLOCK_ACC))
    allocate(block_nscat_accum(blocks%x))
    allocate(block_nscat_accum_d(blocks%x))
  end if

  print*, Nph, threads, blocks
  print*, 'Transmission per-block accumulators: ', use_block_accum


  open(newunit=uT,file='Transmission.txt',status='replace',action='write')
  write(uT,*) n_wl_out, H(1), H(grid%n_lev)

  if (do_cf .eqv. .True.) then
    if (nb_cf < 1) then
      print*, 'ERROR: nb_cf must be at least 1 when do_cf is true.'
      stop
    end if

    allocate(b_cf_grid(nb_cf+1),b_cf_grid_d(nb_cf+1))
    allocate(b_cf(nb_cf),b_cf_d(nb_cf))
    allocate(b_n_cf(nb_cf),b_n_cf_d(nb_cf))
    dH = (H(grid%n_lev) - H(1))/real(nb_cf,dp)
    do i = 1, nb_cf+1
      b_cf_grid(i) = H(1) + real(i-1,dp) * dH
    end do
    b_cf_grid(nb_cf+1) = H(grid%n_lev)
    b_cf_grid_d(:) = b_cf_grid(:)
    b_cf(:) = 0.0_dp
    b_n_cf(:) = 0
    open(newunit=ucf,file='cf_trans.txt',status='replace',action='write')
    write(ucf,*) n_wl_out, nb_cf, H(1), H(grid%n_lev), dH
  end if

  call read_next_opac(s_wl)

  do l = s_wl, n_wl

    call set_grid_opac(iscat)

    im%fsum = 0.0_dp
    im%qsum = 0.0_dp
    im%usum = 0.0_dp
    im%fsum_occ = 0.0_dp
    im%fail_pscat = 0
    im%fail_pemit = 0

    nscat_tot = 0
    nscat_tot_d = nscat_tot

    if (do_cf .eqv. .True.) then
      b_cf_d(:) = b_cf(:)
      b_n_cf_d(:) = b_n_cf(:)
    end if

    T_trans(l) = 0.0_dp
    T_trans_d = T_trans(l)
    T_trans_east(l) = 0.0_dp
    T_trans_east_d = T_trans_east(l)
    T_trans_west(l) = 0.0_dp
    T_trans_west_d = T_trans_west(l)
    if (use_block_accum .eqv. .True.) then
      block_accum_d(:,:) = 0.0_dp
      block_nscat_accum_d(:) = 0
    end if

    l_d = l
    im_d = im
    call reset_raytrace_diagnostics()

    if (LHS .eqv. .True.) then
      if (l == s_wl) then
        ! Allocate CPU and GPU arrays if first call
        allocate(x_ran(Nph),y_ran(Nph),z_ran(Nph),x_ran_d(Nph),y_ran_d(Nph),z_ran_d(Nph))
        call rng_seed(random_seed)
      end if
      ! Generate Nph samples using Latin Hypercube Sampling 
      call LHS_sample(Nph, 2, x_ran, y_ran, z_ran, .False.)
      ! Send samples to GPU memory
      x_ran_d(:) = x_ran(:)
      y_ran_d(:) = y_ran(:)
      z_ran_d(:) = z_ran(:)
    end if

    call exp_3D_sph_atm_transmission_k<<<blocks, threads>>>(l_d, Nph_d)

    call read_next_opac(l+1)

    istat = cudaDeviceSynchronize()
    if (istat /= 0) then
      print*, 'ERROR after exp_3D_sph_atm_transmission:', istat
      stop
    end if

    im = im_d
    if (use_block_accum .eqv. .True.) then
      trans_block_accum(:,:) = block_accum_d(:,:)
      im%fsum = sum(trans_block_accum(:,BLOCK_ACC_F))
      im%qsum = sum(trans_block_accum(:,BLOCK_ACC_Q))
      im%usum = sum(trans_block_accum(:,BLOCK_ACC_U))
      im%fsum_occ = sum(trans_block_accum(:,BLOCK_ACC_F_OCC))
      block_nscat_accum(:) = block_nscat_accum_d(:)
      nscat_tot = sum(block_nscat_accum(:))
    else
      nscat_tot = nscat_tot_d
    end if
    call report_raytrace_diagnostics('3D_sph_trans', 1, l, wl(l), Nph)

    norm = (H(grid%n_lev)**2 - H(1)**2) / (2.0_dp * real(Nph, dp))

    if (do_cf .eqv. .True.) then
      b_cf(:) = b_cf_d(:)
      b_n_cf(:) = b_n_cf_d(:)
      b_cf(:) = norm * b_cf(:)
      write(ucf,*) wl(l),  sum(b_cf(:)),  b_cf(:)
      !call flush(ucf)
      b_cf(:) = 0.0_dp
      b_n_cf(:) = 0
    end if

    if (use_block_accum .eqv. .True.) then
      T_trans(l) = sum(trans_block_accum(:,BLOCK_ACC_TOTAL))
      T_trans_east(l) = sum(trans_block_accum(:,BLOCK_ACC_EAST))
      T_trans_west(l) = sum(trans_block_accum(:,BLOCK_ACC_WEST))
    else
      T_trans(l) = T_trans_d
      T_trans_east(l) = T_trans_east_d
      T_trans_west(l) = T_trans_west_d
    end if

    T_trans(l)      = norm * T_trans(l)
    T_trans_east(l) = norm * T_trans_east(l)
    T_trans_west(l) = norm * T_trans_west(l)

    write(uT,*) wl(l), T_trans(l), T_trans_east(l), T_trans_west(l)
    call flush(uT)

    print*, l, wl(l), T_trans(l), T_trans_east(l), T_trans_west(l)
    print*, 'pscat failures and nscat_tot: ', im%fail_pscat, nscat_tot

  end do

end subroutine exp_3D_sph_atm_transmission
