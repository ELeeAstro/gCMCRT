module exp_1D_pp_atm_kernel
  use mc_precision
  use mc_class_pac
  use mc_class_grid
  use mc_k_source_pac_inc
  use mc_k_emit_iso
  use mc_k_tauint
  use mc_k_scatt
  use cudafor
  use curand_device
  implicit none


  real(dp), dimension(:), allocatable, device :: energy_d, erri_d

contains

  attributes(global) subroutine exp_1D_pp_atm_k(mu_s, n_mu)
    implicit none

    integer, intent(in) :: n_mu
    real(dp), intent(in) :: mu_s

    type(pac) :: ph
    integer :: ii, ij, seq, offset
    integer :: l, istat

    ! Set a random seed for this packet
    ph%id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    ph%seed = ph%id
    seq = 0
    offset = 0
    call curand_init(ph%seed, seq, offset, ph%iseed)

    ph%wght = 1.0_dp
    ph%geo = 1

    !call emit_iso_surf(ph)
    call source_pac_inc_1D(ph, mu_s)

    ph%p_flag = 0
    !! Enter scattering loop
    do while (ph%p_flag == 0)
      !! Sample a tau for the packet
      ph%tau_p = -log(curand_uniform(ph%iseed))

      !! Move packet for sampled tau distance
      call tauint_1D_pp(ph)

      !! Packet exited lower boundary - reemit at surface
      ! if (ph%p_flag == 2) then
      !   call emit_iso_surf(ph)
      !   ph%p_flag = 0
      !   cycle
      ! end if

      if (ph%p_flag == 2) then
        ! Scatter packet with Lambertian
        ph%iscatt = 2
        call scatt_pac(ph)
        !call emit_iso_surf(ph)
        ph%p_flag = 0
        cycle
      else
        exit
      end if



      !! Test for scattering probability
      ! if (curand_uniform(ph%iseed) < ssa_1D_d(ph%c(3))) then
      !   !! Scatter packet isotropically
      !   ph%iscatt = 1
      !   call scatt_pac(ph)
      ! else
      !   !! Remove packet from simulation
      !   exit
      ! end if

    end do

    !! If packet exited top of atmosphere, collect it's angular distribution
    if (ph%p_flag == 1) then
      l = int(real(n_mu,dp)*ph%cost) + 1
      istat = atomicadd(erri_d(l), 1.0_dp)
      istat = atomicadd(energy_d(l), 1.0_dp)
    end if

  end subroutine exp_1D_pp_atm_k

end module exp_1D_pp_atm_kernel


subroutine exp_1D_pp_atm()
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_k_moments
  use exp_1D_pp_atm_kernel
  use cudafor
  implicit none

  integer :: Nph, i, u1, u2
  real(dp) :: mu_s
  real(dp) :: tau0, dz, dtau
  real(dp) :: dthet, halfw, rj

  integer :: n_mu
  real(dp), dimension(:), allocatable :: energy, sigmai, erri, intensity, theta

  type(dim3) :: blocks, threads

  integer, device :: n_mu_d
  real(dp), device :: mu_s_d

  !Nph = 1e7
  Nph = 1024000 * 1
  !Nph = 10

  !! Calculate how many threads and blocks we need
  !threads = dim3(Nph,1,1)
  !blocks = dim3(1,1,1)

  threads= dim3(256,1,1)
  blocks = dim3(ceiling(real(Nph)/threads%x),1,1)

  print*, 'Npackets: '
  print*, Nph
  print*, 'GPU info: '
  print*, 'threads: ', threads
  print*, 'blocks: ', blocks

  print*, '================='

  n_mu = 90 ; n_mu_d = n_mu
  allocate(energy(n_mu), sigmai(n_mu), erri(n_mu), intensity(n_mu))
  energy(:) = 0.0_dp ; sigmai(:) = 0.0_dp ; erri(:) = 0.0_dp ; intensity(:) = 0.0_dp

  dthet = 1.0_dp/real(n_mu,dp)
  halfw = 0.5_dp*dthet
  allocate(theta(n_mu))
  do i = 1, n_mu
    rj = real(i-1,dp)
    theta(i) = acos(rj*dthet+halfw) * 180.0_dp/pi
  end do

  allocate(energy_d(n_mu), erri_d(n_mu))
  energy_d(:) = energy(:)
  erri_d(:) = erri(:)

  mu_s = 1.0_dp
  mu_s_d = mu_s

  grid%z_max = 10.0_dp
  grid%z_min = 0.0_dp
  grid%r_min = 0.0_dp
  grid%n_lay = 10
  grid%n_lev = grid%n_lay + 1

  allocate(rhokap_1D(grid%n_lay),ssa_1D(grid%n_lay),g_1D(grid%n_lay))
  allocate(z(grid%n_lev))

  !! Grid spacing
  dz = grid%z_max/real(grid%n_lay,dp)
  z(1) = 0.0_dp
  do i = 2, grid%n_lev
    z(i) = z(i-1) + dz
  end do

  !! Grid opacity
  tau0 = 0.01_dp
  dtau = tau0/real(grid%n_lay,dp)
  do i = 1, grid%n_lay
    rhokap_1D(i) = dtau/dz
  end do

  ssa_1D(:) = 1.0_dp
  g_1D(:) = 0.0_dp

  !! Allocate arrys and send all data to device
  allocate(rhokap_1D_d(grid%n_lay),ssa_1D_d(grid%n_lay),g_1D_d(grid%n_lay))
  allocate(z_d(grid%n_lev))
  rhokap_1D_d(:) = rhokap_1D(:)
  ssa_1D_d(:) = ssa_1D(:)
  g_1D_d(:) = g_1D(:)
  z_d(:) = z(:)
  grid_d = grid

  !! Allocate moments
  allocate(jp(grid%n_lev),hp(grid%n_lev),kp(grid%n_lev))
  jp(:) = 0.0_dp ; hp(:) = 0.0_dp ; kp(:) = 0.0_dp
  allocate(jp_d(grid%n_lev),hp_d(grid%n_lev),kp_d(grid%n_lev))
  jp_d(:) = jp(:) ; hp_d(:) = hp(:) ; kp_d(:) = kp(:)
  allocate(jm(grid%n_lev),hm(grid%n_lev),km(grid%n_lev))
  jm(:) = 0.0_dp ; hm(:) = 0.0_dp ; km(:) = 0.0_dp
  allocate(jm_d(grid%n_lev),hm_d(grid%n_lev),km_d(grid%n_lev))
  jm_d(:) = jm(:) ; hm_d(:) = hm(:) ; km_d(:) = km(:)


  call exp_1D_pp_atm_k<<<blocks, threads>>>(mu_s_d,n_mu_d)

  ! Give back device data to host
  jp(:) = jp_d(:)/real(Nph,dp);  hp(:) = hp_d(:)/real(Nph,dp); kp(:) = kp_d(:)/real(Nph,dp)
  jm(:) = jm_d(:)/real(Nph,dp);  hm(:) = hm_d(:)/real(Nph,dp); km(:) = km_d(:)/real(Nph,dp)

  energy(:) = energy_d(:)
  erri(:) = erri_d(:)

  intensity(:) = energy(:)/(2.0_dp*real(Nph,dp)*cos(theta(:)*pi/180.0_dp))*real(n_mu,dp)
  sigmai(:) = sqrt(erri(:))/real(Nph,dp)
  energy(:) = energy(:)/real(Nph,dp)

  do i = 1, grid%n_lev
    print*, i, jp(i), hp(i), kp(i), jm(i), hm(i), km(i)
  end do
  print*, '=========='
  do i = 1, n_mu
    print*, i, theta(i), intensity(i), sigmai(i), energy(i)
  end do

  open(newunit=u1,file='moment_gpu_pp.txt',action='readwrite')
  do i = 1, grid%n_lev
    write(u1,*) i, jp(i), hp(i), kp(i), jm(i), hm(i), km(i)
  end do

  open(newunit=u2,file='inten_gpu_pp.txt',action='readwrite')
  do i = 1, n_mu
    write(u2,*) i, theta(i), intensity(i), sigmai(i), energy(i)
  end do



end subroutine exp_1D_pp_atm
