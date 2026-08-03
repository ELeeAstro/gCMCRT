module exp_3D_sph_atm_trans_hires_kernel
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
  use mc_k_raytrace, only : raytrace_sph_3D, checked_transmission_contribution
  use cudafor
  use curand_device
  implicit none

  integer :: nscat_tot
  integer, device :: nscat_tot_d

  type(curandStateMRG32k3a), allocatable, dimension(:), device :: iseed


contains

  attributes(global) subroutine set_iseed(Nph_pad)
    implicit none

    integer, intent(in) :: Nph_pad
    integer(8) :: id, seed, seq, offset

    id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (id > Nph_pad) then
      return
    end if

    seed = random_seed_d
    seq = id - 1
    offset = 0
    call curand_init(seed, seq, offset, iseed(id))

  end subroutine set_iseed


  attributes(global) subroutine exp_3D_sph_atm_trans_hires_k(l, Nph)

    integer, intent(in) :: l, Nph
    type(pac) :: ph, ray
    integer :: seq, offset, i, n, istat, nscat
    real(dp) :: contri, rstat
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
    ph%ig = 1

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

    ! Perform peeloff at starting location
    ray = ph
    call raytrace_sph_3D(ray)

    call checked_transmission_contribution(ray, contri, valid_direct_ray)
    if (valid_direct_ray .eqv. .True.) then
      if (use_block_accum_d .eqv. .True.) then
        rstat = atomicadd(block_accum_d(blockIdx%x,BLOCK_ACC_TOTAL),contri)
      else
        rstat = atomicadd(T_trans_d,contri)
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

      if(ph%p_flag /= 0) then
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

  end subroutine exp_3D_sph_atm_trans_hires_k


end module exp_3D_sph_atm_trans_hires_kernel

subroutine exp_3D_sph_atm_trans_hires()
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_class_imag
  use exp_3D_sph_atm_trans_hires_kernel
  use mc_k_raytrace, only : reset_raytrace_diagnostics, report_raytrace_diagnostics
  use mc_opacset
  use mc_read_prf
  use LHS_sampling_mod, only : LHS_sample
  use random_cpu
  use mc_k_scatt, only : validate_volume_scattering
  use cudafor
  use mc_los_velocity
  implicit none


  character (len=8) :: fmt
  character (len=3) :: n_str
  character (len=512) :: nml_iomsg
  integer :: Nph, l, iscat, istat, n, nml_iostat
  integer, allocatable, dimension(:) :: uT
  integer, device :: l_d, Nph_d
  integer :: n_theta, n_phi, n_lay
  real(dp) :: viewthet
  real(dp), allocatable, dimension(:) :: viewphi
  real(dp) :: pl, pc, sc
  real(dp), allocatable, dimension(:,:) :: trans_block_accum
  integer, allocatable, dimension(:) :: block_nscat_accum
  logical :: use_block_accum

  type(dim3) :: blocks, threads

  namelist /sph_3D_trans_hires/ Nph, n_wl, pl, pc, sc, n_theta, n_phi, &
    & viewthet, viewphi, n_lay, iscat, use_block_accum

  if (n_phase < 1) then
    print*, 'ERROR: n_phase must be positive in 3D_sph_trans_hi.'
    stop
  end if
  allocate(uT(n_phase),viewphi(n_phase))

  use_block_accum = .True.
  read(u_nml, nml=sph_3D_trans_hires, iostat=nml_iostat, iomsg=nml_iomsg)
  if (nml_iostat /= 0) then
    print*, 'ERROR reading &sph_3D_trans_hires from CMCRT.nml: ', trim(nml_iomsg)
    stop
  end if

  if (Nph < 1) then
    print*, 'ERROR: Nph must be positive in 3D_sph_trans_hi.'
    stop
  end if
  call validate_volume_scattering(iscat, '3D_sph_trans_hi')

  if ((lbl .eqv. .False.) .or. (ck .eqv. .True.)) then
    print*, 'ERROR: 3D_sph_trans_hi currently supports line-by-line opacity only.'
    print*, 'Set lbl = .True. and ck = .False.'
    stop
  end if

  if (doppler_on .eqv. .True.) then
    if (inc_lbl .eqv. .False.) then
      print*, 'ERROR: Doppler high-resolution transmission currently requires inc_lbl = .True.'
      stop
    end if
    if (inc_xsec .eqv. .True.) then
      print*, 'ERROR: inc_xsec is not supported by the Doppler opacity reader.'
      stop
    end if
  end if

  ng = 1
  ng_d = 1

  fmt = '(I3.3)'

  ! Give namelist paramaters to equilvanet values inside gCMCRT
  grid%n_lay = n_lay
  grid%n_lev = n_lay + 1
  grid%n_theta = n_theta
  grid%n_phi = n_phi

  if (cmd_vphi .eqv. .False.) then
    print*, 'Using namelist vphi'
    im%vphi = viewphi(1)
    !write(vphi_arg , *) viewphi
  else
    print*, 'Using cmdline vphi'
    viewphi(:) = im%vphi
  end if
  print*, im%vphi, viewphi(:)

  im%vtheta = viewthet

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
  call read_wl()

  call set_grid()

  if (LHS .eqv. .True.) then
    ! Allocate the LHS storage once.  Samples are regenerated below for each
    ! phase/wavelength, but the allocatables retain one well-defined lifetime.
    allocate(x_ran(Nph), y_ran(Nph), z_ran(Nph))
    allocate(x_ran_d(Nph), y_ran_d(Nph), z_ran_d(Nph))
    call rng_seed(random_seed)
  end if

  if (doppler_on .eqv. .True.) then
    call read_1D_wprf()
    call compute_vlos(viewphi(:))
  end if

  ! Send data to GPU data containers
  im_d = im
  grid_d = grid

  allocate(T_trans(n_wl))

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
  print*, 'High-resolution transmission per-block accumulators: ', use_block_accum


  do n = 1, n_phase
    write(n_str,fmt) n
    open(newunit=uT(n),file='Transmission_'//trim(n_str)//'.txt',status='replace',action='write')
    write(uT(n),*) n_wl, H(1), H(grid%n_lev),  viewphi(n)
    call flush(uT(n))
  end do


  if (doppler_on .eqv. .True.) then
    call read_next_opac_doppler(1)
  else
    call read_next_opac(1)
  end if

  do l = 1, n_wl

    do n = 1, n_phase


      if (cmd_vphi .eqv. .False.) then
        im%vphi = viewphi(n)
      else
        viewphi(n) = im%vphi
      end if
      im%vtheta = viewthet

      call set_image()

      if (doppler_on .eqv. .True.) then
        call shift_opac(n,l)
      end if

      call set_grid_opac(iscat)

      im%fsum = 0.0_dp
      im%qsum = 0.0_dp
      im%usum = 0.0_dp
      im%fsum_occ = 0.0_dp
      im%fail_pscat = 0
      im%fail_pemit = 0

      nscat_tot = 0
      nscat_tot_d = nscat_tot

      T_trans(l) = 0.0_dp
      T_trans_d = T_trans(l)
      if (use_block_accum .eqv. .True.) then
        block_accum_d(:,:) = 0.0_dp
        block_nscat_accum_d(:) = 0
      end if

      l_d = l
      im_d = im
      call reset_raytrace_diagnostics()

      if (LHS .eqv. .True.) then
        ! Generate Nph samples using Latin Hypercube Sampling
        call LHS_sample(Nph, 2, x_ran, y_ran, z_ran, .False.)
        ! Send samples to GPU memory
        x_ran_d(:) = x_ran(:)
        y_ran_d(:) = y_ran(:)
        z_ran_d(:) = z_ran(:)
      end if

      call exp_3D_sph_atm_trans_hires_k<<<blocks, threads>>>(l_d,Nph_d)

      if (n == n_phase) then
        if (doppler_on .eqv. .True.) then
          call read_next_opac_doppler(l+1)
        else
          call read_next_opac(l+1)
        end if
      end if

      istat = cudaDeviceSynchronize()
      if (istat /= 0) then
        print*, 'ERROR after exp_3D_sph_atm_trans_hires:', istat
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
      call report_raytrace_diagnostics('3D_sph_trans_hi', n, l, wl(l), Nph)

      ! Give T_trans_d back to CPU
      if (use_block_accum .eqv. .True.) then
        T_trans(l) = sum(trans_block_accum(:,BLOCK_ACC_TOTAL))
      else
        T_trans(l) = T_trans_d
      end if

      T_trans(l) = (H(grid%n_lev)**2 - H(1)**2) / (2.0_dp * real(Nph,dp)) * T_trans(l)
      write(uT(n),*) wl(l), T_trans(l)
      !call flush(uT(n))

      print*, n, l, wl(l), T_trans(l)
      print*, 'pscat failures and nscat_tot: ', im%fail_pscat, nscat_tot

    end do

  end do

end subroutine exp_3D_sph_atm_trans_hires
