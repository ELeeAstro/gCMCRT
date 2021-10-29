module exp_3D_sph_atm_kernel
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
    if (id > Nph) then
      return
    end if

    seed = id + id**2 + id/2
    seq = 0
    offset = 0
    call curand_init(seed, seq, offset, iseed(id))

  end subroutine set_iseed

  attributes(global) subroutine exp_3D_sph_atm_k(Nph)
    implicit none

    integer, intent(in) :: Nph
    type(pac) :: ph, ray
    integer :: seq, offset, i, n, istat
    real(dp) :: contri

    ! Set a random seed for this packet
    ph%id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (ph%id > Nph) then
      return
    end if
    ph%iseed = iseed(ph%id)


    ph%wght = 1.0_dp
    ph%geo = 2
    ph%ig = 1

    call source_pac_inc_3D(ph)
    call findcell(ph)

    !call gord_samp(ph)

    ! Begin scattering loop
    ph%p_flag = 0

    n = 0
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

      ph%wght = ph%wght * ssa_d(ph%ig,ph%c(1),ph%c(2),ph%c(3))

      !if (n == 0) then
        ph%iscatt = 3
      !else
        !ph%iscatt = 1
      !end if

      call peeloff_scatt(ph)
      call scatt_pac(ph)
      call RR_test(ph)

      n = n + 1

    end do

    ! Give back iseed to saved device array for next iteration with this ph%id
    iseed(ph%id) = ph%iseed

  end subroutine exp_3D_sph_atm_k


end module exp_3D_sph_atm_kernel


subroutine exp_3D_sph_atm()
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_class_imag
  use cudafor
  use exp_3D_sph_atm_kernel
  implicit none

  integer :: Nph, l, uT, i, j, k
  integer, device :: l_d, Nph_d
  integer :: n_theta, n_phi, n_lay
  real(dp) :: viewthet, viewphi, dr, dr3, v_tot, tau0, dtau
  real(dp) :: pl, pc, sc

  type(dim3) :: blocks, threads

  pl = 0.51
  pc = 0.39
  sc = 1.0

  pl_d = pl
  pc_d = pc
  sc_d = sc

  !!Nph = 1e7
  !Nph = 1024000 * 1
  Nph = 1024

  threads = dim3(128, 1, 1)
  blocks = dim3(ceiling(real(Nph,dp)/threads%x),1,1)
  allocate(iseed(Nph))
  Nph_d = Nph
  call set_iseed<<<blocks, threads>>>(Nph_d)

  !! Calculate how many threads and blocks we need
  ! Nph = 10
  ! threads = dim3(Nph,1,1)
  ! blocks = dim3(1,1,1)

  print*, 'Npackets: '
  print*, Nph
  print*, 'GPU info: '
  print*, 'threads: ', threads
  print*, 'blocks: ', blocks

  print*, '================='

  grid%n_phi = 145
  grid%n_theta = 45
  grid%n_lev = 11
  grid%n_lay = grid%n_lev - 1

  allocate(thetarr(grid%n_theta),phiarr(grid%n_phi),r(grid%n_lev))

  grid%r_min = 1.0_dp
  grid%r_max = 11.0_dp
  grid%r2_max = grid%r_max**2

  dr = (grid%r_max - grid%r_min)/real(grid%n_lay,dp)
  r(1) = grid%r_min
  do i = 2, grid%n_lay
    r(i) = r(i-1) + dr
    !print*, i, r(i), dr
  end do
  r(grid%n_lev) =  grid%r_max

  !! phi (longitude) grid set up
  ! 1st grid = 0 degrees
  allocate(aarr(grid%n_phi),barr(grid%n_phi))
  phiarr(1) = 0.0_dp
  aarr(1) = sin(phiarr(1))
  barr(1) = -cos(phiarr(1))

  ! Spacing in longitude
  dphi = twopi / real(grid%n_phi-1,kind=dp)
  do j = 2, grid%n_phi - 1
    phiarr(j) = phiarr(j-1) + dphi
    aarr(j) = sin(phiarr(j))
    barr(j) = -cos(phiarr(j))
    !print*, j, phiarr(j), aarr(j), barr(j), dphi
  end do

  ! last grid = 360 degrees
  phiarr(grid%n_phi) = twopi
  aarr(grid%n_phi) = sin(phiarr(grid%n_phi))
  barr(grid%n_phi) = -cos(phiarr(grid%n_phi))

  !! theta (latitude) grid set up

  ! Set up faces of theta-grid.
  ! Also set up tan(theta)**2 arrays for finding constant theta surfaces.
  ! If theta is close to 90degrees, set arrays to -1 -- used in findwall.

  ! Spacing in latitude
  dtheta = pi / real(grid%n_theta-1,kind=dp)

  allocate(tan2thet(grid%n_theta))
  ! 1st gird = 0 degrees
  thetarr(1) = 0.0_dp
  tan2thet(1) = 0.0_dp

  do k = 2, grid%n_theta-1
    thetarr(k) = thetarr(k-1) + dtheta
    if ((thetarr(k) > (1.57060_dp)).and.(thetarr(k) < (1.57090_dp))) then
      tan2thet(k) = -1.0_dp
      !tan2thet(k) = tan(thetarr(k))**2
    else
      tan2thet(k) = tan(thetarr(k))**2
    endif
    !print*, k, thetarr(k), tan2thet(k), dtheta
  end do

  ! last grid = pi
  thetarr(grid%n_theta) = pi
  tan2thet(grid%n_theta) = 0.0_dp

  ! Find volume of each cell - spherical coordinates
  allocate(v_cell(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  do i = 1, grid%n_lay
    dr3 = r(i+1)**3 - r(i)**3
    do j = 1, grid%n_phi-1
      dphi = phiarr(j+1) - phiarr(j)
      do k = 1, grid%n_theta-1
        dtheta = cos(thetarr(k)) - cos(thetarr(k+1))
        v_cell(i,j,k) = (dr3 * dphi * dtheta) / 3.0_dp
        !print*, i, j, k, v_cell(i,j,k), dr3, dphi, dtheta
      end do
    end do
  end do

  v_tot = (4.0_dp/3.0_dp) * pi * (r(grid%n_lev)**3 - r(1)**3)
  print*, 'Volumes Check:', v_tot - sum(v_cell), v_tot, sum(v_cell), sum(v_cell)/v_tot


  ! Homogenous opacity
  tau0 = 5.0_dp
  dtau = (grid%r_max - grid%r_min)/tau0


  allocate(rhokap(1,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  do i = 1, grid%n_lay
    rhokap(1,i,:,:) = dtau/dr
  end do


  allocate(ssa(1,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  ssa(1,:,:,:) = 0.99_dp!0.9999_dp!0.999_dp!1.0_dp

  allocate(gg(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  gg(:,:,:) =  0.0_dp


  !! Observation direction and image set up
  im%vtheta = 90.0_dp
  im%vphi = 90.0_dp
  xpix = 200
  ypix = 200
  rimage = 1.1_dp

  call set_image()

  im_d = im

  allocate(r_d(grid%n_lev),theta_d(grid%n_theta),phi_d(grid%n_phi))
  allocate(tan2thet_d(grid%n_theta))
  allocate(aarr_d(grid%n_phi),barr_d(grid%n_phi))
  allocate(rhokap_d(1,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  allocate(g_d(grid%n_lay,grid%n_phi-1,grid%n_theta-1),ssa_d(1,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  r_d(:) = r(:)
  phi_d(:) = phiarr(:)
  theta_d(:) = thetarr(:)
  tan2thet_d(:) = tan2thet(:)
  aarr_d(:) = aarr(:)
  barr_d(:) = barr(:)
  rhokap_d(1,:,:,:) = rhokap(1,:,:,:)
  ssa_d(1,:,:,:) = ssa(1,:,:,:)
  g_d(:,:,:) = gg(:,:,:)
  grid_d = grid

  call exp_3D_sph_atm_k<<<blocks, threads>>>(Nph_d)

  im = im_d
  f(:,:) = f_d(:,:)/real(Nph,dp)
  q(:,:) = q_d(:,:)/real(Nph,dp)
  u(:,:) = u_d(:,:)/real(Nph,dp)

  im%fsum = im%fsum/real(Nph,dp)
  im%qsum = im%qsum/real(Nph,dp)
  im%usum = im%usum/real(Nph,dp)

  print*, 'Emitted: ', Nph
  print*, 'Asked for: ', Nph
  print*, 'Total f, q, u :', im%fsum, im%qsum, im%usum

  print*, 'Check: ', ssa(1,1,1,1), im%fsum * pi, im%qsum * pi

  print*, im%vphi, im%fsum, im%qsum, im%usum

  print*, 'degree %: ', im%qsum/im%fsum * 100.0_dp, sqrt(im%qsum**2 + im%usum**2)/im%fsum * 100.0_dp

  print*, 'fimage'

  ! Output images
  open(newunit=uT,file='fimage.dat',status='unknown',form='unformatted')
  write(uT) real(f)
  close(uT)

  print*, 'qimage'


  open(newunit=uT,file='qimage.dat',status='unknown',form='unformatted')
  write(uT) real(q)
  close(uT)

  print*, 'uimage'


  open(newunit=uT,file='uimage.dat',status='unknown',form='unformatted')
  write(uT) real(u)
  close(uT)

  print*, 'fimage_ascii'

  ! Output images
  open(newunit=uT,file='fimage_ascii.dat',action='readwrite',form='formatted')
  write(uT,*)  im%fsum
  do i = 1, im%x_pix
    write(uT,*) (f(i,j), j = 1, im%y_pix)
  end do
  close(uT)

  print*, 'qimage_ascii'


  open(newunit=uT,file='qimage_ascii.dat',action='readwrite',form='formatted')
  write(uT,*)  im%qsum
  do i = 1, im%x_pix
    write(uT,*) (q(i,j), j = 1, im%y_pix)
  end do
  close(uT)

  print*, 'uimage_ascii'


  open(newunit=uT,file='uimage_ascii.dat',action='readwrite',form='formatted')
  write(uT,*)  im%usum
  do i = 1, im%x_pix
    write(uT,*) (u(i,j), j = 1, im%y_pix)
  end do
  close(uT)


end subroutine exp_3D_sph_atm
