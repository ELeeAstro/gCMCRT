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
  use mc_k_raytrace
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
    integer(8) :: id, seed
    integer :: seq, offset

    id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (id > Nph) then
      return
    end if

    seed = id + id**2 + id/2
    seq = 0
    offset = 0
    call curand_init(seed, seq, offset, iseed(id))

  end subroutine set_iseed


  attributes(global) subroutine exp_3D_sph_atm_transmission_k(l,Nph)
    implicit none

    integer, intent(in) :: l, Nph
    type(pac) :: ph, ray
    integer :: seq, offset, i, n, istat, nscat
    integer :: b_cf_idx
    real(dp) :: contri

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

    contri = ray%wght * (1.0_dp - exp(-ray%tau)) * ray%bp * H_d(1)
    istat = atomicadd(T_trans_d(l),contri)

    ! Do the contibution function for the binned b
    if (do_cf_d .eqv. .True.) then
      call locate(b_cf_grid_d,ray%bp*H_d(1),b_cf_idx)
      istat = atomicadd(b_cf_d(b_cf_idx),ray%wght * (1.0_dp - exp(-ray%tau)) * ray%bp * H_d(1))
      istat = atomicadd(b_n_cf_d(b_cf_idx),1)
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

      if (curand_uniform(ph%iseed) < dorg_d(ph%c(1),ph%c(2),ph%c(3))) then
        ! Gas scattering - Rayleigh scattering
        ph%wght = ph%wght * ssa_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))
        ph%iscatt = 3
      else
        ! Cloud scattering - do given scattering phase function
        ph%wght = ph%wght * ssa_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))
        ph%iscatt = iscat_d
      end if

      call peeloff_scatt(ph)
      call scatt_pac(ph)
      call RR_test(ph)

      nscat =  nscat + 1

    end do

    ! Add number of scatterings to total
    istat = atomicadd(nscat_tot_d, nscat)

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
  use mc_opacset
  use mc_read_prf
  use LHS_sampling_mod, only : LHS_sample_2D
  use cudafor
  implicit none


  integer :: Nph, l, uT, iscat, ucf, nb_cf, i, s_wl
  integer, device :: l_d, Nph_d, nb_cf_d
  integer :: n_theta, n_phi, n_lay
  real(dp) :: dH
  real(dp) :: viewthet, viewphi
  real(dp) :: pl, pc, sc

  integer :: istat
  type(dim3) :: blocks, threads

  namelist /sph_3D_trans/ Nph, s_wl, n_wl, pl, pc, sc, n_theta, n_phi, n_lay, viewthet, viewphi, iscat, nb_cf

  read(u_nml, nml=sph_3D_trans)

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
  call read_wl()
  print*, 'read gord'
  call read_g_ord()


  call set_grid()
  call set_image()


  ! Send data to GPU data containers
  im_d = im
  grid_d = grid

  allocate(T_trans(n_wl),T_trans_d(n_wl))

  ! Grid for GPU threads/blocks
  threads = dim3(128,1,1)
  blocks = dim3(ceiling(real(Nph)/threads%x),1,1)

  print*, Nph, threads, blocks


  open(newunit=uT,file='Transmission.txt',action='readwrite')
  write(uT,*) n_wl, H(1), H(grid%n_lev)

  if (do_cf .eqv. .True.) then
    allocate(b_cf_grid(nb_cf),b_cf_grid_d(nb_cf))
    allocate(b_cf(nb_cf),b_cf_d(nb_cf))
    allocate(b_n_cf(nb_cf),b_n_cf_d(nb_cf))
    dH = (H(grid%n_lev) - H(1))/real(nb_cf-1,dp)
    b_cf_grid(1) = H(1)
    do i = 2, nb_cf
        b_cf_grid(i) = b_cf_grid(i-1) + dH
       !print*, i, b_cf_grid(i), H(1), H(grid%n_lev)
    end do
    b_cf_grid_d(:) = b_cf_grid(:)
    b_cf(:) = 0.0_dp
    b_n_cf(:) = 0
    open(newunit=ucf,file='cf_trans.txt',action='readwrite')
    write(ucf,*) n_wl, nb_cf, H(1), H(grid%n_lev), dH
  end if

  call read_next_opac(s_wl)

  do l = s_wl, n_wl

    call set_grid_opac()

    im%fsum = 0.0_dp
    im%qsum = 0.0_dp
    im%usum = 0.0_dp
    im%fail_pscat = 0
    im%fail_pemit = 0

    nscat_tot = 0
    nscat_tot_d = nscat_tot

    if (do_cf .eqv. .True.) then
      b_cf_d(:) = b_cf(:)
      b_n_cf_d(:) = b_n_cf(:)
    end if

    T_trans(l) = 0.0_dp
    T_trans_d(l) = T_trans(l)

    l_d = l
    im_d = im

    if (LHS .eqv. .True.) then
      if (l == s_wl) then
        ! Allocate CPU and GPU arrays if first call
        allocate(x_ran(Nph),y_ran(Nph),x_ran_d(Nph),y_ran_d(Nph))
        call random_seed()
      end if
      ! Generate Nph samples using Latin Hypercube Sampling 
      call LHS_sample_2D(Nph, x_ran, y_ran, 'lhs', .False., 1000, 'random_cd')
      ! Send samples to GPU memory
      x_ran_d(:) = x_ran(:)
      y_ran_d(:) = y_ran(:)
    end if

    call exp_3D_sph_atm_transmission_k<<<blocks, threads>>>(l_d, Nph_d)

    call read_next_opac(l+1)

    istat = cudaDeviceSynchronize()

    im = im_d
    nscat_tot = nscat_tot_d

    ! Give T_trans_d back to CPU
    T_trans(l) = T_trans_d(l)

    if (do_cf .eqv. .True.) then
      b_cf(:) = b_cf_d(:)
      b_n_cf(:) = b_n_cf_d(:)
      do i = 1, nb_cf-1
        b_cf(i) = dH/ real(b_n_cf(i),dp) * b_cf(i)
      end do
      write(ucf,*) wl(l),  sum(b_cf(:)),  b_cf(:)
      !call flush(ucf)
      b_cf(:) = 0.0_dp
      b_n_cf(:) = 0
    end if

    T_trans(l) = (H(grid%n_lev) - H(1)) / real(Nph,dp) * T_trans(l)
    write(uT,*) wl(l), T_trans(l)
    call flush(uT)

    print*, l, wl(l), T_trans(l)
    print*, 'pscat failures and nscat_tot: ', im%fail_pscat, nscat_tot

  end do

end subroutine exp_3D_sph_atm_transmission
