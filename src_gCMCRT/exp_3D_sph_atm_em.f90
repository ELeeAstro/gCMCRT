module exp_3D_sph_atm_em_kernel
  use mc_precision
  use mc_class_pac
  use mc_class_grid
  use mc_k_tauint
  use mc_k_scatt
  use mc_k_tau_samp
  use mc_k_RR
  use mc_k_gord_samp
  use mc_k_emit_iso
  use mc_k_peeloff_emit
  use mc_k_peeloff_scatt
  use mc_k_vol_samp
  use cudafor
  use curand_device
  implicit none

  integer :: nscat_tot
  integer, device :: nscat_tot_d

  integer,allocatable,dimension(:) :: Nph_i, Nph_j, Nph_k
  integer,allocatable,dimension(:),device :: Nph_i_d, Nph_j_d, Nph_k_d

  real(dp), allocatable, dimension(:,:,:) :: wght_start
  real(dp), allocatable, dimension(:,:,:), device :: wght_start_d



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

  attributes(global) subroutine exp_3D_sph_atm_em_k(l,Nph_sum)
    implicit none
    integer, intent(in) :: l, Nph_sum
    type(pac) :: ph
    integer :: seq, offset, nscat, istat
    integer :: i, j, k

    ! Set a random seed for this packet
    ph%id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (ph%id > Nph_sum) then
      return
    end if
    ph%iseed = iseed(ph%id)

    i = Nph_i_d(ph%id)
    j = Nph_j_d(ph%id)
    k = Nph_k_d(ph%id)

    ph%wght = wght_start_d(i,j,k)
    ph%geo = 2
    ph%wl = wl_d(l)

    call sph_samp_3D(i,j,k,ph)

    ! Sample a g-ordinance value (for corr-k)
    if (do_g_bias_d .eqv. .True.) then
      call gord_samp_cell_bias(ph)
    else
      call gord_samp_cell(ph)
    end if

    call emit_iso(ph)

    ! Perform emission peeloff at starting location
    call peeloff_emit(ph)

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

      ! Defensive guard before indexing device arrays.
      if ((ph%c(1) < 1) .or. (ph%c(1) >= grid_d%n_lev) .or. &
          (ph%c(2) < 1) .or. (ph%c(2) >= grid_d%n_phi) .or. &
          (ph%c(3) < 1) .or. (ph%c(3) >= grid_d%n_theta)) then
        ph%p_flag = -777
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

      nscat = nscat + 1

      if (nscat > 1e4) then
        exit
      end if 

    end do

    ! Add number of scatterings to total
    if (use_block_accum_d .eqv. .True.) then
      istat = atomicadd(block_nscat_accum_d(blockIdx%x),nscat)
    else
      istat = atomicadd(nscat_tot_d,nscat)
    end if

    ! Give back iseed to saved device array for next iteration with this ph%id
    iseed(ph%id) = ph%iseed

  end subroutine exp_3D_sph_atm_em_k


end module exp_3D_sph_atm_em_kernel

subroutine exp_3D_sph_atm_em()
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_class_imag
  use exp_3D_sph_atm_em_kernel
  use mc_opacset
  use mc_read_prf
  use mc_set_em
  use mc_k_scatt, only : validate_volume_scattering
  use mc_views
  use LHS_sampling_mod, only : LHS_sample
  use random_cpu
  use cudafor
  implicit none

  integer :: Nph_tot, Nph_sum, l, Nph_pad, n_lay, iscat, n_loop
  character (len=8) :: fmt
  character (len=3) :: n_str
  character (len=512) :: nml_iomsg
  integer, allocatable, dimension(:) :: uT
  integer, device :: Nph_pad_d
  integer :: i, j, k, n, nn, s_wl, n_wl_out
  integer, device :: Nph_sum_d, l_d
  integer :: n_theta, n_phi, istat, nml_iostat
  real(dp) :: viewphi, viewthet
  real(dp) :: pl, pc, sc
  real(dp) :: xi_emb
  real(dp),allocatable,dimension(:) :: em_out, em_out_occ
  real(dp),allocatable,dimension(:,:) :: block_accum
  integer,allocatable,dimension(:) :: block_nscat_accum
  logical :: use_block_accum

  integer :: id
  integer,allocatable,dimension(:,:,:) :: Nph_cell

  type(dim3) :: blocks, threads


  namelist /sph_3D_em/ Nph_tot, s_wl, n_wl, pl, pc, sc, n_theta, n_phi, &
    & viewthet, viewphi, n_lay, xi_emb, iscat, use_block_accum


  use_block_accum = .True.
  read(u_nml, nml=sph_3D_em, iostat=nml_iostat, iomsg=nml_iomsg)
  if (nml_iostat /= 0) then
    print*, 'ERROR reading &sph_3D_em from CMCRT.nml: ', trim(nml_iomsg)
    stop
  end if

  if (Nph_tot < 1) then
    print*, 'ERROR: Nph_tot must be positive in 3D_sph_em.'
    stop
  end if
  if ((xi_emb < 0.0_dp) .or. (xi_emb > 1.0_dp)) then
    print*, 'ERROR: xi_emb must lie in [0,1] in 3D_sph_em.'
    stop
  end if
  call validate_volume_scattering(iscat, '3D_sph_em')

  fmt = '(I3.3)'

  ! Give namelist paramaters to equilvanet values inside gCMCRT
  grid%n_lay = n_lay
  grid%n_lev = n_lay + 1
  grid%n_theta = n_theta
  grid%n_phi = n_phi

  pl_d = pl
  pc_d = pc
  sc_d = sc
  iscat_d = iscat

  Nph_pad = Nph_tot
  threads = dim3(128, 1, 1)
  blocks = dim3(ceiling(real(Nph_pad,dp)/threads%x),1,1)
  allocate(iseed(Nph_pad))
  allocate(Nph_i(Nph_pad),Nph_j(Nph_pad),Nph_k(Nph_pad))
  allocate(Nph_i_d(Nph_pad),Nph_j_d(Nph_pad),Nph_k_d(Nph_pad))
  use_block_accum_d = use_block_accum
  if (use_block_accum .eqv. .True.) then
    allocate(block_accum(blocks%x,N_BLOCK_ACC))
    allocate(block_accum_d(blocks%x,N_BLOCK_ACC))
    allocate(block_nscat_accum(blocks%x))
    allocate(block_nscat_accum_d(blocks%x))
  end if
  Nph_pad_d = Nph_pad
  call set_iseed<<<blocks, threads>>>(Nph_pad_d)

  print*, 'Emission per-block accumulators: ', use_block_accum

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

  allocate(Nph_cell(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  allocate(wght_start(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  allocate(wght_start_d(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  allocate(em_out(n_wl), em_out_occ(n_wl))

  ! Calculate viewphi from n_phase
  if (do_phase .eqv. .False. .or. n_phase <= 1) then
    ! Single phase, use namelist values - need to check for command line passage
    im%vphi = viewphi
    im%vtheta = viewthet
    allocate(viewphi_n(1), viewthet_n(1), phi_key(1), occ_state(1))
    viewphi_n(1) = viewphi
    viewthet_n(1) = viewthet
    phi_key(1) = -1.0
    occ_state(1) = 0
    occ_state_d = 0
    n_loop = 1
  else
    ! Do full phase curve operation - splitting to regular section and eclipse section
    call set_phase_views_em(n_phase, n_ecl, n_loop)
  end if

  ! Write output files ready for input
  allocate(uT(n_loop))
  do n = 1, n_loop
    write(n_str,fmt) n
    open(newunit=uT(n),file='Em_'//trim(n_str)//'.txt',status='replace',action='write')
    write(uT(n),*) n_wl_out, H(1), H(grid%n_lev), viewphi_n(n), viewthet_n(n), phi_key(n)
    call flush(uT(n))
  end do


  call read_next_opac(s_wl)

  call rng_seed(random_seed)

  print*, 'starting loop'

  do l = s_wl, n_wl

    call set_grid_opac(iscat)
    call set_grid_em(l)

    ! The intrinsic emission proposal depends on wavelength and atmospheric
    ! state, not viewing phase.  Build and transfer it once per wavelength so
    ! every phase reuses the same packet-to-cell allocation.
    call allocate_emission_packets(Nph_tot, xi_emb, Nph_cell, wght_start, Nph_sum)

    if (Nph_sum < 128) then
      threads = dim3(Nph_sum, 1, 1)
      blocks = dim3(1,1,1)
    else
      threads = dim3(128, 1, 1)
      blocks = dim3(ceiling(real(Nph_sum,dp)/threads%x),1,1)
    end if

    !! Flatten the packet allocation into the reusable packet-to-cell buffers.
    id = 1
    do k = 1, grid%n_theta-1
      do j = 1, grid%n_phi-1
        do i = 1, grid%n_lay
          do nn = 1, Nph_cell(i,j,k)
            Nph_i(id) = i
            Nph_j(id) = j
            Nph_k(id) = k
            id = id + 1
          end do
        end do
      end do
    end do
    if (id /= Nph_sum + 1) then
      print*, 'ERROR: flattened emission packet count is inconsistent.'
      print*, 'Expected final id, actual id: ', Nph_sum + 1, id
      stop
    end if

    l_d = l
    Nph_sum_d = Nph_sum
    Nph_i_d(1:Nph_sum) = Nph_i(1:Nph_sum)
    Nph_j_d(1:Nph_sum) = Nph_j(1:Nph_sum)
    Nph_k_d(1:Nph_sum) = Nph_k(1:Nph_sum)
    wght_start_d(:,:,:) = wght_start(:,:,:)

    do n = 1, n_loop


      im%vphi = viewphi_n(n)
      im%vtheta = viewthet_n(n)

      call set_image()


      im%fsum = 0.0_dp
      im%qsum = 0.0_dp
      im%usum = 0.0_dp
      im%fail_pscat = 0
      im%fail_pemit = 0
      im%fsum_occ = 0.0_dp

      im_d = im

      if (do_phase .eqv. .True.) then
        occ_state_d = occ_state(n)
      end if

      if (occ_state(n) == 1) then
        xstar_d = xstar(n)
        ystar_d = ystar(n)
      end if

      nscat_tot = 0
      nscat_tot_d = nscat_tot
      if (use_block_accum .eqv. .True.) then
        block_accum_d(:,:) = 0.0_dp
        block_nscat_accum_d(:) = 0
      end if

      if (do_images .eqv. .True.) then
        f(:,:) = 0.0_dp ; q(:,:) = 0.0_dp ; u(:,:) = 0.0_dp ; im_err(:,:) = 0.0_dp
        f_d(:,:) = f(:,:) ; q_d(:,:) = q(:,:) ; u_d(:,:) = u(:,:) ; im_err_d(:,:) = im_err(:,:)
      end if

      call exp_3D_sph_atm_em_k<<<blocks, threads>>>(l_d,Nph_sum_d)

      if (n == n_loop) then
        call read_next_opac(l+1)
      end if

      istat = cudaDeviceSynchronize()
      if (istat /= 0) then
        print*, 'ERROR after exp_3D_sph_atm_em:', istat
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

      em_out(l) = im%fsum / real(Nph_tot,dp)
      em_out_occ(l) = im%fsum_occ / real(Nph_tot,dp)

      write(uT(n),*) wl(l), em_out(l), em_out_occ(l), grid%lumtot
      !call flush(uT)

      if (do_cf .eqv. .True.) then
        cf(:,:,:) = cf_d(:,:,:)
        call output_cf(n,l,n_loop)
        cf_d(:,:,:) = 0.0_dp
      end if

      if (do_images .eqv. .True.) then
        f(:,:) = f_d(:,:)/real(Nph_tot,dp) ; q(:,:) = q_d(:,:)/real(Nph_tot,dp)
        u(:,:) = u_d(:,:)/real(Nph_tot,dp) ; im_err(:,:) = im_err_d(:,:)
        call output_im(n,l,n_loop)
      end if


      print*, l, real(wl(l)), Nph_tot, Nph_sum, real(em_out(l)), real(em_out_occ(l)), grid%lumtot
      print*, n, viewphi_n(n), viewthet_n(n), phi_key(n)
      print*, n, 'pemit, pscat failures and nscat_tot: ',  im%fail_pemit, im%fail_pscat, nscat_tot

    end do
  
  end do

  deallocate(Nph_i,Nph_j,Nph_k)
  deallocate(Nph_i_d,Nph_j_d,Nph_k_d)

end subroutine exp_3D_sph_atm_em
