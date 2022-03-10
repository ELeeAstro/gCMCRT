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
  use mc_k_raytrace
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
    integer(8) :: id, seed
    integer :: seq, offset

    id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (id > Nph_pad) then
      return
    end if

    seed = id + id**2 + id/2
    seq = 0
    offset = 0
    call curand_init(seed, seq, offset, iseed(id))

  end subroutine set_iseed


  attributes(global) subroutine exp_3D_sph_atm_trans_hires_k(l, Nph)

    integer, intent(in) :: l, Nph
    type(pac) :: ph, ray
    integer :: seq, offset, i, n, istat, nscat
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

    contri = ray%wght * (1.0_dp - exp(-ray%tau)) * ray%bp * H_d(1)
    istat = atomicadd(T_trans_d(l),contri)

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

  end subroutine exp_3D_sph_atm_trans_hires_k


end module exp_3D_sph_atm_trans_hires_kernel

subroutine exp_3D_sph_atm_trans_hires()
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_class_imag
  use exp_3D_sph_atm_trans_hires_kernel
  use mc_opacset
  use mc_read_prf
  use cudafor
  use mc_los_velocity
  implicit none


  character (len=8) :: fmt
  character (len=3) :: n_str
  integer :: Nph, l, iscat, istat, n
  integer, allocatable, dimension(:) :: uT
  integer, device :: l_d, Nph_d
  integer :: n_theta, n_phi, n_lay
  real(dp) :: viewthet
  real(dp), allocatable, dimension(:) :: viewphi
  real(dp) :: pl, pc, sc

  type(dim3) :: blocks, threads

  namelist /sph_3D_trans_hires/ Nph, n_wl, pl, pc, sc, n_theta, n_phi, viewthet, viewphi, n_lay, iscat

  allocate(uT(n_phase),viewphi(n_phase))

  read(u_nml, nml=sph_3D_trans_hires)

  fmt = '(I3.3)'

  ! Give namelist paramaters to equilvanet values inside gCMCRT
  grid%n_lay = n_lay
  grid%n_lev = n_lay + 1
  grid%n_theta = n_theta
  grid%n_phi = n_phi

  im%vtheta = viewthet
  im%vphi = viewphi(1)

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
  call read_1D_wprf()
  call read_wl()

  call set_grid()
  call set_image()

  if (doppler_on .eqv. .True.) then
    call compute_vlos(viewphi(:))
  end if

  ! Send data to GPU data containers
  im_d = im
  grid_d = grid

  allocate(T_trans(n_wl),T_trans_d(n_wl))

  ! Grid for GPU threads/blocks
  threads = dim3(128,1,1)
  blocks = dim3(ceiling(real(Nph)/threads%x),1,1)

  print*, Nph, threads, blocks


  do n = 1, n_phase
    write(n_str,fmt) n
    open(newunit=uT(n),file='Transmission_'//trim(n_str)//'.txt',action='readwrite')
    write(uT(n),*) n_wl, H(1), H(grid%n_lev), viewphi(n)
  end do


  if (doppler_on .eqv. .True.) then
    call read_next_opac_doppler(1)
  else
    call read_next_opac(1)
  end if

  do l = 1, n_wl

    do n = 1, n_phase

      im%vphi = viewphi(n)
      im%vtheta = viewthet

      if (doppler_on .eqv. .True.) then
        call shift_opac(n,l)
      end if

      call set_grid_opac()

      im%fsum = 0.0_dp
      im%qsum = 0.0_dp
      im%usum = 0.0_dp
      im%fail_pscat = 0
      im%fail_pemit = 0

      nscat_tot = 0
      nscat_tot_d = nscat_tot

      T_trans(l) = 0.0_dp
      T_trans_d(l) = T_trans(l)

      l_d = l
      im_d = im
      call exp_3D_sph_atm_trans_hires_k<<<blocks, threads>>>(l_d,Nph_d)

      !istat = cudaDeviceSynchronize()

      im = im_d
      nscat_tot = nscat_tot_d

      ! Give T_trans_d back to CPU
      T_trans(l) = T_trans_d(l)

      T_trans(l) = (H(grid%n_lev) - H(1)) / real(Nph,dp) * T_trans(l)
      write(uT(n),*) wl(l), T_trans(l)
      call flush(uT(n))

      print*, n, l, wl(l), T_trans(l)
      print*, 'pscat failures and nscat_tot: ', im%fail_pscat, nscat_tot

    end do

    if (doppler_on .eqv. .True.) then
      call read_next_opac_doppler(l+1)
    else
      call read_next_opac(l+1)
    end if

  end do

end subroutine exp_3D_sph_atm_trans_hires
