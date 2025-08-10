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
    integer(8) :: id, seed
    integer :: seq, offset

    ! Get packet id
    id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (id > Nph) then
      return
    end if

    ! Set seed for packet
    seed = id + id**2 + id/2
    seq = 0
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
  use cudafor
  implicit none


  integer :: Nph, l, iscat, istat, s_wl, n
  integer, allocatable, dimension(:) :: uT
  integer, device :: l_d, Nph_d
  integer :: n_theta, n_phi, n_lay
  character (len=8) :: fmt
  character (len=3) :: n_str
  real(dp) :: viewthet
  real(dp) :: pl, pc, sc
  real(dp), allocatable, dimension(:) :: viewphi

  type(dim3) :: blocks, threads


  namelist /sph_3D_alb/ Nph, s_wl, n_wl, pl, pc, sc, n_theta, n_phi, n_lay, viewthet, viewphi, iscat

  allocate(uT(n_phase),viewphi(n_phase))

  read(u_nml, nml=sph_3D_alb)

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
  blocks = dim3(ceiling(real(Nph,dp)/threads%x),1,1)
  allocate(iseed(Nph))
  Nph_d = Nph
  call set_iseed<<<blocks, threads>>>(Nph_d)

  call read_1D_prf()
  call read_wl()
  call read_g_ord()

  call set_grid()

  ! Send data to GPU data containers
  grid_d = grid

  allocate(alb_out(n_wl),alb_out_d(n_wl))

  do n = 1, n_phase
    write(n_str,fmt) n
    open(newunit=uT(n),file='Alb_'//trim(n_str)//'.txt',action='readwrite')
    write(uT(n),*) n_wl, H(1), H(grid%n_lev), viewphi(n)
    call flush(uT(n))
  end do

  call read_next_opac(s_wl)

  print*, 'starting loop'

  do l = s_wl, n_wl
 
    call set_grid_opac()

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
      im%fail_pscat = 0
      im%fail_pemit = 0

      im_d = im
       
      nscat_tot = 0
      nscat_tot_d = nscat_tot

      alb_out(l) = 0.0_dp
      alb_out_d(l) = alb_out(l)

      f(:,:) = 0.0_dp ; q(:,:) = 0.0_dp ; u(:,:) = 0.0_dp ; im_err(:,:) = 0.0_dp
      f_d(:,:) = f(:,:) ; q_d(:,:) = q(:,:) ; u_d(:,:) = u(:,:) ; im_err_d(:,:) = im_err(:,:)

      l_d = l

      if (LHS .eqv. .True.) then
        if (l == s_wl) then
          ! Allocate CPU and GPU arrays if first call
          allocate(x_ran(Nph),y_ran(Nph),z_ran(Nph),x_ran_d(Nph),y_ran_d(Nph),z_ran_d(Nph))
          !call random_seed()
          call rng_seed(213)
        end if
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

      im = im_d
      nscat_tot = nscat_tot_d

      ! Give fsum back to CPU
      alb_out(l) = im%fsum / real(Nph,dp) * pi

      write(uT(n),*) wl(l), alb_out(l)
      !call flush(uT)

      if (do_images .eqv. .True.) then
        f(:,:) = f_d(:,:)/real(Nph,dp) ; q(:,:) = q_d(:,:)/real(Nph,dp)
        u(:,:) = u_d(:,:)/real(Nph,dp) ; im_err(:,:) = im_err_d(:,:)
        call output_im(n,l)
      end if

      print*, n, l, real(wl(l)), alb_out(l)
      print*, n, 'pscat failures and nscat_tot: ', im%fail_pscat, nscat_tot

    end do
  end do

end subroutine exp_3D_sph_atm_albedo
