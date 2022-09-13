module exp_3D_sph_atm_pol_kernel
  use mc_precision
  use mc_class_pac
  use mc_class_grid
  use mc_k_source_pac_inc
  use mc_k_emit_iso
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

    id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (id > Nph) then
      return
    end if

    seed = id + id**2 + id/2
    seq = 0
    offset = 0
    call curand_init(seed, seq, offset, iseed(id))

  end subroutine set_iseed


  attributes(global) subroutine exp_3D_sph_atm_pol_k(l,Nph)
    implicit none

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

    call source_pac_inc_3D(ph)

    call findcell(ph)

    ! Begin scattering loop
    ph%p_flag = 0

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

      !print*, ph%c(1),ph%c(2),ph%c(3)

      if (ph%p_flag == 1) then
        ph%iscatt = 2
        call peeloff_scatt(ph)
      !  call emit_iso_surf_sph(ph)
        call scatt_pac(ph)
        ph%p_flag = 0
        nscat = nscat + 1
        cycle
      else
        !exit
        !print*, ph%c(1),ph%c(2),ph%c(3)
      end if

      if (ph%p_flag /= 0) then
        !print*, ph%id, ph%wght, ph%p_flag, 'died'
        exit
      end if

      ! Cloud scattering - do given scattering phase function
      ph%wght = ph%wght * ssa_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))
      ph%iscatt = 3

      call peeloff_scatt(ph)
      call scatt_pac(ph)
      call RR_test(ph)

      nscat =  nscat + 1

    end do

    istat = atomicadd(nscat_tot_d, nscat)

    ! Give back iseed to saved device array for next iteration with this ph%id
    iseed(ph%id) = ph%iseed

  end subroutine exp_3D_sph_atm_pol_k


end module exp_3D_sph_atm_pol_kernel

subroutine exp_3D_sph_atm_pol()
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_class_imag
  use mc_Draine_G
  use exp_3D_sph_atm_pol_kernel
  use mc_opacset
  use mc_read_prf
  use cudafor
  implicit none


  integer :: Nph, l, uT, i, istat
  integer, device :: l_d, Nph_d
  integer :: n_theta, n_phi, n_lay
  real(dp) :: viewthet, viewphi
  real(dp) :: pl, pc, sc, tau0, dtau, dr

  type(dim3) :: blocks, threads


  namelist /sph_3D_alb/ Nph, n_wl, pl, pc, sc, n_theta, n_phi, n_lay, viewthet, viewphi

  read(u_nml, nml=sph_3D_alb)

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

  threads = dim3(128, 1, 1)
  blocks = dim3(ceiling(real(Nph,dp)/threads%x),1,1)
  allocate(iseed(Nph))
  Nph_d = Nph
  call set_iseed<<<blocks, threads>>>(Nph_d)

  ! Calculate height replacing - call read_1D_prf()
  allocate(H(grid%n_lev),H_d(grid%n_lev))
  do i = 1, grid%n_lev
    H(i) = i
  end do
  H_d(:) = real(H(:),dp)

  call read_wl()

  call set_grid()
  call set_image()


  ! Send data to GPU data containers
  im_d = im
  grid_d = grid

  allocate(alb_out(n_wl),alb_out_d(n_wl))

  ! Grid for GPU threads/blocks
  threads = dim3(128,1,1)
  blocks = dim3(ceiling(real(Nph)/threads%x),1,1)

  print*, Nph, threads, blocks


  open(newunit=uT,file='Albedo.txt',action='readwrite')
  write(uT,*) n_wl, H(1), H(grid%n_lev)


  ng = 1
  allocate(rhokap_d(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  allocate(ssa_d(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  allocate(g_d(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  allocate(dorg_d(grid%n_lay,grid%n_phi-1,grid%n_theta-1))


  ! Homogenous opacity
  tau0 = 0.001_dp
  dtau = tau0/real(grid%n_lay,dp) !(grid%r_max - grid%r_min)/tau0
  dr = (grid%r_max - grid%r_min)!/real(grid%n_lay,dp)


  do i = 1, grid%n_lay
    rhokap_d(1,i,:,:) = dtau/(r(i+1) - r(i)) * grid%r_del
  end do

  ssa_d(:,:,:,:) = 1.0_dp !0.9999_dp
  g_d(:,:,:) = 0.0_dp
  dorg_d(:,:,:) = 0.0_dp


  do l = 1, n_wl

    !call set_grid_opac()

    im%fsum = 0.0_dp
    im%qsum = 0.0_dp
    im%usum = 0.0_dp
    im%fail_pscat = 0
    im%fail_pemit = 0

    nscat_tot = 0
    nscat_tot_d = nscat_tot

    alb_out(l) = 0.0_dp
    alb_out_d(l) = alb_out(l)

    f(:,:) = 0.0_dp ; q(:,:) = 0.0_dp ; u(:,:) = 0.0_dp ; im_err(:,:) = 0.0_dp
    f_d(:,:) = f(:,:) ; q_d(:,:) = q(:,:) ; u_d(:,:) = u(:,:) ; im_err_d(:,:) = im_err(:,:)

    l_d = l
    im_d = im
    call exp_3D_sph_atm_pol_k<<<blocks, threads>>>(l_d, Nph_d)

    istat = cudaDeviceSynchronize()

    im = im_d
    nscat_tot = nscat_tot_d

    ! Give fsum back to CPU
    im%fsum = im%fsum/real(Nph,dp)
    im%qsum = im%qsum/real(Nph,dp)
    im%usum = im%usum/real(Nph,dp)

    alb_out(l) = im%fsum * pi

    write(uT,*) wl(l), alb_out(l), im%fsum, im%qsum, im%usum
    call flush(uT)

    if (do_images .eqv. .True.) then
      f(:,:) = f_d(:,:)/real(Nph,dp)
      q(:,:) = q_d(:,:)/real(Nph,dp)
      u(:,:) = u_d(:,:)/real(Nph,dp)
      im_err(:,:) = im_err_d(:,:)/real(Nph,dp)
      call output_im(1,l)
    end if

    print*, l, wl(l), alb_out(l), im%qsum*pi, im%usum*pi
    print*, im%fsum, im%qsum,im%usum,sqrt(im%qsum**2 + im%usum**2)/im%fsum, -im%qsum/im%fsum


    print*, 'pscat failures and nscat_tot: ', im%fail_pscat, nscat_tot

  end do

end subroutine exp_3D_sph_atm_pol
