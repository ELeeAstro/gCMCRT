module exp_3D_sph_atm_albedo_kernel
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

    ! Get packet id
    id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (id > Nph) then
      return
    end if

    ! Set seed for packet
    seed = random_seed_d
    seq = id - 1
    offset = 0
    call curand_init(seed, seq, offset, iseed(id))

  end subroutine set_iseed


  attributes(global) subroutine exp_3D_sph_atm_albedo_k(l,Nph)
    implicit none

    integer, intent(in) :: l, Nph
    type(pac) :: ph, ray
    integer :: seq, offset, istat, nscat
    real(dp) :: contri

    ! Get the iseed for this packet
    ph%id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (ph%id > Nph) then
      return
    end if
    ph%iseed = iseed(ph%id)

    ! Initial weight, geometry and wavelength
    ph%wght = 1.0_dp
    ph%geo = 2
    ph%wl = wl_d(l)

    ! Source packet from external source
    call source_pac_inc_3D(ph)

    ! Find the initial cell number
    call findcell(ph)

    ! Sample a g-ordinance value (for corr-k)
    if (do_g_bias_d .eqv. .True.) then
      call gord_samp_bias(ph)
    else
      call gord_samp(ph)
    end if

    if (do_scat_loop_d .eqv. .True.) then
      ph%p_flag = 0
    else
      ph%p_flag = -1
    end if

    ! Begin scattering loop

    ! Number of scattering events counter
    nscat = 0

    !! Enter scattering loop
    do while (ph%p_flag == 0)
      !! Sample a tau for the packet
      !if (nscat > 0) then
       ph%tau_p = -log(curand_uniform(ph%iseed))
      !else
        !call tau_force_scatt(ph)
        !call tau_force_stretch(ph)
      !end if

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

  end subroutine exp_3D_sph_atm_albedo_k


end module exp_3D_sph_atm_albedo_kernel

subroutine exp_3D_sph_atm_albedo()
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_class_imag
  use exp_3D_sph_atm_albedo_kernel
  use mc_opacset
  use mc_read_prf
  use LHS_sampling_mod, only : LHS_sample
  use random_cpu
  use mc_k_scatt, only : validate_volume_scattering
  use cudafor
  implicit none


  integer :: Nph, l, iscat, istat, s_wl, n, nml_iostat, n_wl_out
  integer, allocatable, dimension(:) :: uT
  integer, device :: l_d, Nph_d
  integer :: n_theta, n_phi, n_lay
  character (len=8) :: fmt
  character (len=3) :: n_str
  character (len=512) :: nml_iomsg
  real(dp) :: viewthet
  real(dp) :: pl, pc, sc
  real(dp), allocatable, dimension(:) :: viewphi
  real(dp), allocatable, dimension(:,:) :: block_accum
  integer, allocatable, dimension(:) :: block_nscat_accum
  logical :: use_block_accum

  type(dim3) :: blocks, threads


  namelist /sph_3D_alb/ Nph, s_wl, n_wl, pl, pc, sc, n_theta, n_phi, n_lay, &
    & viewthet, viewphi, iscat, use_block_accum

  if (n_phase < 1) then
    print*, 'ERROR: n_phase must be positive in 3D_sph_alb.'
    stop
  end if
  allocate(uT(n_phase),viewphi(n_phase))

  use_block_accum = .True.
  read(u_nml, nml=sph_3D_alb, iostat=nml_iostat, iomsg=nml_iomsg)
  if (nml_iostat /= 0) then
    print*, 'ERROR reading &sph_3D_alb from CMCRT.nml: ', trim(nml_iomsg)
    stop
  end if

  if (Nph < 1) then
    print*, 'ERROR: Nph must be positive in 3D_sph_alb.'
    stop
  end if
  call validate_volume_scattering(iscat, '3D_sph_alb')

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

  im%vtheta = viewthet

  pl_d = pl
  pc_d = pc
  sc_d = sc
  iscat_d = iscat

  threads = dim3(128, 1, 1)
  blocks = dim3(Nph / threads%x,1,1)
  if (mod(Nph,threads%x) /= 0) blocks%x = blocks%x + 1
  allocate(iseed(Nph))
  use_block_accum_d = use_block_accum
  if (use_block_accum .eqv. .True.) then
    allocate(block_accum(blocks%x,N_BLOCK_ACC))
    allocate(block_accum_d(blocks%x,N_BLOCK_ACC))
    allocate(block_nscat_accum(blocks%x))
    allocate(block_nscat_accum_d(blocks%x))
  end if
  Nph_d = Nph
  call set_iseed<<<blocks, threads>>>(Nph_d)

  print*, 'Albedo per-block accumulators: ', use_block_accum

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

  ! Send data to GPU data containers
  grid_d = grid

  allocate(alb_out(n_wl),alb_out_d(n_wl))

  do n = 1, n_phase
    write(n_str,fmt) n
    open(newunit=uT(n),file='Alb_'//trim(n_str)//'.txt',status='replace',action='write')
    write(uT(n),*) n_wl_out, H(1), H(grid%n_lev), viewphi(n)
    call flush(uT(n))
  end do

  if (LHS .eqv. .True.) then
    allocate(x_ran(Nph),y_ran(Nph),z_ran(Nph),x_ran_d(Nph),y_ran_d(Nph),z_ran_d(Nph))
    call rng_seed(random_seed)
  end if

  call read_next_opac(s_wl)

  print*, 'starting loop'

  do l = s_wl, n_wl
 
    call set_grid_opac(iscat)

    do n = 1, n_phase

      if (cmd_vphi .eqv. .False.) then
        im%vphi = viewphi(n)
      else
        viewphi(n) = im%vphi
      end if
      im%vtheta = viewthet

      call set_image()

      im%fsum = 0.0_dp
      im%qsum = 0.0_dp
      im%usum = 0.0_dp
      im%fsum_occ = 0.0_dp
      im%fail_pscat = 0
      im%fail_pemit = 0

      im_d = im
       
      nscat_tot = 0
      nscat_tot_d = nscat_tot
      if (use_block_accum .eqv. .True.) then
        block_accum_d(:,:) = 0.0_dp
        block_nscat_accum_d(:) = 0
      end if

      alb_out(l) = 0.0_dp
      alb_out_d(l) = alb_out(l)

      if (do_images .eqv. .True.) then
        f(:,:) = 0.0_dp ; q(:,:) = 0.0_dp ; u(:,:) = 0.0_dp ; im_err(:,:) = 0.0_dp
        f_d(:,:) = f(:,:) ; q_d(:,:) = q(:,:) ; u_d(:,:) = u(:,:) ; im_err_d(:,:) = im_err(:,:)
      end if

      l_d = l

      if (LHS .eqv. .True.) then
        ! Generate Nph samples using Latin Hypercube Sampling
        call LHS_sample(Nph, 2, x_ran, y_ran, z_ran, .False.)
        ! Send samples to GPU memory
        x_ran_d(:) = x_ran(:)
        y_ran_d(:) = y_ran(:)
        z_ran_d(:) = z_ran(:)
      end if

      call exp_3D_sph_atm_albedo_k<<<blocks, threads>>>(l_d, Nph_d)

      if (n == n_phase) then
        call read_next_opac(l+1)
      end if

      istat = cudaDeviceSynchronize()
      if (istat /= 0) then
        print*, 'ERROR after exp_3D_sph_atm_albedo_k:', istat
        stop
      end if

      im = im_d
      if (use_block_accum .eqv. .True.) then
        block_accum(:,:) = block_accum_d(:,:)
        im%fsum = sum(block_accum(:,BLOCK_ACC_F))
        im%qsum = sum(block_accum(:,BLOCK_ACC_Q))
        im%usum = sum(block_accum(:,BLOCK_ACC_U))
        im%fsum_occ = sum(block_accum(:,BLOCK_ACC_F_OCC))
        block_nscat_accum(:) = block_nscat_accum_d(:)
        nscat_tot = sum(block_nscat_accum(:))
      else
        nscat_tot = nscat_tot_d
      end if

      ! Give fsum back to CPU
      alb_out(l) = im%fsum / real(Nph,dp) * pi

      write(uT(n),*) wl(l), alb_out(l)
      !call flush(uT)

      if (do_images .eqv. .True.) then
        f(:,:) = f_d(:,:)/real(Nph,dp) ; q(:,:) = q_d(:,:)/real(Nph,dp)
        u(:,:) = u_d(:,:)/real(Nph,dp) ; im_err(:,:) = im_err_d(:,:)
        call output_im(n,l,n_phase)
      end if

      print*, n, l, real(wl(l)), alb_out(l)
      print*, n, 'pscat failures and nscat_tot: ', im%fail_pscat, nscat_tot

    end do
  end do

end subroutine exp_3D_sph_atm_albedo
