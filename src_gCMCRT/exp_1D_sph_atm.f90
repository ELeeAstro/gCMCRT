module exp_1D_sph_atm_kernel
  use mc_precision
  use mc_class_pac
  use mc_class_grid
  use mc_k_source_pac_inc
  use mc_k_emit_iso
  use mc_k_tauint
  use mc_k_scatt
  use mc_k_limb_dark
  use mc_k_moments
  use cudafor
  use curand_device
  implicit none


  real(dp), device :: energy_tot_d, erri_tot_d
  real(dp), dimension(:), allocatable, device :: energy_d, erri_d
  real(dp), dimension(:,:), allocatable, device :: image_sph_d

contains

  attributes(global) subroutine exp_1D_sph_atm_k(mu_s, n_mu, fractx, nx)
    implicit none

    integer, intent(in) :: n_mu, nx
    real(dp), intent(in) :: mu_s, fractx

    type(pac) :: ph
    integer :: seq, offset
    integer :: l, istat, ix, iy
    real(dp) :: xim, yim

    ! Set a random seed for this packet
    ph%id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    ph%seed = ph%id !+ ph%id**2
    seq = 0
    offset = 0
    !call curand_init(ph%seed, seq, offset, ph%iseed)

    ph%wght = 1.0_dp
    ph%geo = 2

    if (limb_dark_d == 0) then
      call emit_iso_surf(ph)
    else
      call emit_iso_surf_sph(ph)
    end if

    ph%p_flag = 0
    !! Enter scattering loop
    do while (ph%p_flag == 0)
      !! Sample a tau for the packet
      ph%tau_p = -log(curand_uniform(ph%iseed))

      !! Move packet for sampled tau distance
      call tauint_1D_sph(ph)

      !! Packet exited lower boundary - reemit at surface
      ! if (ph%p_flag == 2) then
      !   call emit_iso_surf(ph)
      !   ph%zp = grid%r_min
      !   ! Add the positive moment values from the surface if required
      !   if (do_moments_d .eqv. .True.) then
      !     call moments_1D(ph)
      !   end if
      !   cycle
      ! end if

      !! Test for scattering probability
      if (curand_uniform(ph%iseed) < ssa_1D_d(ph%c(3))) then
        !! Scatter packet isotropically
        ph%iscatt = 1
        call scatt_pac(ph)
      else
        !! Remove packet from simulation
        exit
      end if

    end do

    !! If packet exited top of atmosphere, collect it's angular distribution
    !! Add packet to image data
    if (ph%p_flag == 1) then
      istat = atomicadd(erri_tot_d, 1.0_dp)
      istat = atomicadd(energy_tot_d, 1.0_dp)

      l = int(real(n_mu,dp)*ph%cost) + 1
      istat = atomicadd(erri_d(l), 1.0_dp)
      istat = atomicadd(energy_d(l), 1.0_dp)

      xim = ph%yp*ph%cosp-ph%xp*ph%sinp
      yim = ph%zp*ph%sint-ph%yp*ph%cost*ph%sinp-ph%xp*ph%cost*ph%cosp
      ix = int((xim+grid_d%r_max)/fractx) + 1
      iy = int((yim+grid_d%r_max)/fractx) + 1
      if (ix<=nx.and.iy<=nx.and.ix>0.and.iy>0) then
         istat = atomicadd(image_sph_d(ix,iy), 1.0_dp)
      else
         print*,'error in making image ', ix, iy
      endif

    end if

  end subroutine exp_1D_sph_atm_k

end module exp_1D_sph_atm_kernel


subroutine exp_1D_sph_atm()
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_k_moments
  use mc_k_limb_dark
  use exp_1D_sph_atm_kernel
  use cudafor
  implicit none

  integer :: Nph, i, j, u1, u2, u3
  real(dp) :: mu_s
  real(dp) :: tau0, dz, dtau
  real(dp) :: dthet, halfw, rj

  integer :: rhoexp, limb_dark

  integer :: opb, nx, ntab
  real(dp) :: rtau, iopb, b, c, dri, dr, fractx, delmu
  real(dp), dimension(:), allocatable :: xmu, prob

  integer :: n_mu
  real(dp), dimension(:), allocatable :: energy, sigmai, erri, intensity, theta, tau
  real(dp), dimension(:,:), allocatable :: image_sph
  real(dp) :: energy_tot, erri_tot

  type(dim3) :: blocks, threads

  integer, device :: n_mu_d, nx_d
  real(dp), device :: mu_s_d, fractx_d

  !Nph = 1e7
  Nph = 1024000 * 1
  !Nph = 10

  !! Calculate how many threads and blocks we need
  !threads = dim3(Nph,1,1)
  !blocks = dim3(1,1,1)

  !threads= dim3(1024,1,1)
  threads= dim3(256,1,1)
  !threads= dim3(128,1,1)
  blocks = dim3(ceiling(real(Nph)/threads%x),1,1)

  print*, 'Npackets: '
  print*, Nph
  print*, 'GPU info: '
  print*, 'threads: ', threads
  print*, 'blocks: ', blocks

  print*, '================='

  n_mu = 30 ; n_mu_d = n_mu
  allocate(energy(n_mu), sigmai(n_mu), erri(n_mu), intensity(n_mu))
  energy(:) = 0.0_dp ; sigmai(:) = 0.0_dp ; erri(:) = 0.0_dp ; intensity(:) = 0.0_dp
  energy_tot = 0.0_dp
  erri_tot = 0.0_dp

  allocate(energy_d(n_mu), erri_d(n_mu))
  energy_d(:) = energy(:)
  erri_d(:) = erri(:)
  energy_tot_d = energy_tot
  erri_tot_d = erri_tot


  dthet = 1.0_dp/real(n_mu,dp)
  halfw = 0.5_dp*dthet
  allocate(theta(n_mu))
  do i = 1, n_mu
    rj = real(i-1,dp)
    theta(i) = acos(rj*dthet+halfw) * 180.0_dp/pi
  end do


  !! Grid values here
  ! rstar=1.0_dp
  ! rmine=1.0_dp
  ! rmaxe=40.0_dp
  ! taumax=10.0_dp
  rhoexp=-1
  ! albedo=1.0_dp
  limb_dark = 0
  nx = 99
  ntab = 110
  opb = 1
  nx_d = 99

  allocate(image_sph(nx,nx))
  image_sph(:,:) = 0.0_dp
  allocate(image_sph_d(nx,nx))
  image_sph_d(:,:) = image_sph(:,:)

  mu_s = 1.0_dp
  mu_s_d = mu_s

  !! --- spherical stellar atmosphere basic set-up --- !!
  !opb determines grid spacing; r(i)=c*i**iopb
  iopb = 1.0_dp/opb
  print*, 'iopb: ', iopb

  grid%n_lev = 20
  grid%n_lay = grid%n_lev - 1
  grid%r_min = 1.0_dp
  grid%r_max = 40.0_dp
  tau0 = 10.0_dp


  dri = (grid%r_max-grid%r_min)/(real(grid%n_lev,dp)-1.0_dp)**opb
  print*, 'dri',dri
  allocate(r(grid%n_lev))
  do i = 1, grid%n_lev
    r(i) = grid%r_min + dri*(real(i,dp)-1.0_dp)**opb
  enddo


  !set one of the grid points to rmine, for higher accuracy
  r(grid%n_lev) = grid%r_max
  ! if (rmine.gt.rmin) then
  !   call locate(ri,nri,rmine,i)
  !   ri(i+1)=rmine
  ! endif

  ! density goes as rho**(rhoexp)
  ! calculate constant of proportionality so envelope tau=taumax
  if (rhoexp /= -1.0_dp) then
    c = 1.0_dp+rhoexp
    b = tau0*c/(grid%r_max**c-grid%r_min**c)
  else
    b = tau0/(log(grid%r_max)-log(grid%r_min))
  endif

  allocate(rhokap_1D(grid%n_lay),ssa_1D(grid%n_lay),g_1D(grid%n_lay))
  allocate(tau(grid%n_lev))

  rtau = 0.0_dp
  tau(1) = rtau
  do i = 1, grid%n_lay
    dr = r(i+1) - r(i)
    ! since density is described by an analytic equation easy
    ! to integrate, set optical depth of cell equal to
    ! integral(rho*dr)
    if (r(i) >= grid%r_min) then
      if (rhoexp /= -1.0_dp) then
        rhokap_1D(i) = b/c*(r(i+1)**c-r(i)**c)/dr
      else
        rhokap_1D(i) = b*(log(r(i+1))-log(r(i)))/dr
      endif
      rtau = rtau + rhokap_1D(i)*dr
      tau(i+1) = rtau
    else
      rhokap_1D(i) = 0.0_dp
    endif
    !print*, i, z(i), rhokap_1D(i), rtau
  end do

  print*, 'optical depth through spherical envelope = ',rtau
  print*, 'should equal taumax ',tau0

  print*,'done with opacity grid '
  print*, ' '

  fractx = 2.0_dp*grid%r_max/real(nx,dp)
  fractx_d = fractx

  !! Do the limb darkening probability ditribution
  delmu = 1.0_dp/(real(ntab,dp))
  allocate(xmu(ntab),prob(ntab))
  xmu(1) = 0.0_dp
  prob(1) = 0.0_dp
  do i = 2, ntab
    xmu(i) = delmu*real(i,dp)
    prob(i) = 0.5_dp*(xmu(i)**3+xmu(i)**2)
  end do

  print*, 'in table, prob(n), should equal 1, ',prob(ntab)
  allocate(limb_cdf_d(ntab),xmu_d(ntab))
  limb_dark_d = limb_dark
  limb_cdf_d(:) = prob(:)
  xmu_d(:) = xmu(:)
  delmu_d = delmu


  !! ---------------------- !!

  !! --- simple homogenous opacity distribution (same as default pp test) --- !!
  ! grid%r_max = 10.0_dp
  ! grid%r_min = 0.0_dp
  ! grid%n_lay = 10
  ! grid%n_lev = grid%n_lay + 1

  !allocate(rhokap_1D(grid%n_lay),ssa_1D(grid%n_lay),g_1D(grid%n_lay))
  ! allocate(r(grid%n_lev))

  !! Grid spacing
  ! dr = grid%r_max/real(grid%n_lay,dp)
  ! r(1) = 0.0_dp
  ! do i = 2, grid%n_lev
  !   r(i) = r(i-1) + dr
  ! end do

  !! Grid opacity
  ! tau0 = 10.0_dp
  ! dtau = tau0/real(grid%n_lay,dp)
  ! do i = 1, grid%n_lay
  !   rhokap_1D(i) = dtau/dr
  ! end do
  !! ---------------------- !!

  ssa_1D(:) = 1.0_dp
  g_1D(:) = 0.0_dp

  !! Allocate arrys and send all data to device
  allocate(rhokap_1D_d(grid%n_lay),ssa_1D_d(grid%n_lay),g_1D_d(grid%n_lay))
  allocate(r_d(grid%n_lev))
  rhokap_1D_d(:) = rhokap_1D(:)
  ssa_1D_d(:) = ssa_1D(:)
  g_1D_d(:) = g_1D(:)
  r_d(:) = r(:)
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


  call exp_1D_sph_atm_k<<<blocks, threads>>>(mu_s_d,n_mu_d,fractx_d,nx_d)

  ! Give back device data to host
  jp(:) = jp_d(:)/real(Nph,dp);  hp(:) = hp_d(:)/real(Nph,dp); kp(:) = kp_d(:)/real(Nph,dp)
  jm(:) = jm_d(:)/real(Nph,dp);  hm(:) = hm_d(:)/real(Nph,dp); km(:) = km_d(:)/real(Nph,dp)

  energy(:) = energy_d(:)
  erri(:) = erri_d(:)
  image_sph(:,:) = image_sph_d(:,:)

  energy_tot = energy_tot_d/real(Nph,dp)
  erri_tot = erri_tot_d

  intensity(:) = energy(:)/(2.0_dp*real(Nph,dp)*cos(theta(:)*pi/180.0_dp))*real(n_mu,dp)
  sigmai(:) = sqrt(erri(:))/real(Nph,dp)
  energy(:) = energy(:)/real(Nph,dp)
  image_sph(:,:) = image_sph(:,:)/real(Nph,dp)

  do i = 1, grid%n_lev
    print*, i, r(i),(tau0-tau(i)), jp(i), hp(i), kp(i), jm(i), hm(i), km(i)
  end do
  print*, '=========='
  do i = 1, n_mu
    print*, i, theta(i), intensity(i), sigmai(i), energy(i)
  end do
  print*, '=========='
  print*, '(energy per steradian)*4*pi ', ', error'
  print*, energy_tot, energy_tot * 1.0_dp/sqrt(erri_tot)

  open(newunit=u1,file='moment_gpu_sph.txt',action='readwrite')
  do i = 1, grid%n_lev
    write(u1,*) i, r(i), (tau0-tau(i)), jp(i), hp(i), kp(i), jm(i), hm(i), km(i)
  end do

  open(newunit=u2,file='inten_gpu_sph.txt',action='readwrite')
  do i = 1, n_mu
    write(u2,*) i, theta(i), intensity(i), sigmai(i), energy(i)
  end do

  open(unit=u3,file='image_sph.txt',action='readwrite')
  do j = 1, nx
     write(u3,*) (image_sph(i,j),i=1,nx)
  enddo
  close(u3)




end subroutine exp_1D_sph_atm
