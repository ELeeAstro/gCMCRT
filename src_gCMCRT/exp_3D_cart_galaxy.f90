!!! Kernel routine for the 3D diffuse galaxy test - this is run on the device (GPU)
module diffuse_kernel
  use mc_precision
  use mc_class_pac
  use mc_class_grid
  use mc_class_imag
  use mc_data_mod
  use cudafor
  use curand_device
  use mc_k_vol_samp
  use mc_k_emit_iso
  use mc_k_scatt
  use mc_k_tau_samp
  use mc_k_tauint
  use mc_k_peeloff_emit
  use mc_k_peeloff_scatt
  !use mc_k_RR
  implicit none


contains


  attributes(global) subroutine diffuse_k(Nph, i, j, k)
    implicit none

    integer, intent(in) :: Nph, i, j, k

    integer :: l, ii, ij, seq, offset, aflag
    type(pac) :: ph

    ! Set a random seed for this packet
    ph%id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    ph%seed = ph%id + ph%id**2 + ph%id/2
    seq = 0
    offset = 0
    call curand_init(ph%seed, seq, offset, ph%iseed)

    ! Give some information to the packet class
    ph%wght = 1.0_dp
    ph%geo = 1

    ! First find random position in this cell
    call cart_samp_3D(i,j,k,ph)

    ! Emit the packet isotropically
    call emit_iso(ph)

    ! Do the emission peeloff
    !call peeloff_emit(ph)

    ph%p_flag = 0
    !! Scattering loop - packet scatters until exits simulation
    do while (ph%p_flag == 0)
      !ph%tau = -log(curand_uniform(ph%iseed))
      call tau_force_scatt(ph)
      if (ph%p_flag /= 0) then
        exit
      end if

      ! Integrate tau
      call tauint_cart_3D(ph)
      !print*, ph%seed, ph%cell(1), ph%cell(2), ph%cell(3), aflag
      if (ph%p_flag /= 0) then
        exit
      end if

      if (curand_uniform(ph%iseed) < ssa_d(ph%c(1),ph%c(2),ph%c(3))) then
        !ph%wght = ph%wght * albedo
        ! Scatter peeloff
        ph%iscatt = 4
        call peeloff_scatt(ph)
        ! New direction for packet
        ph%iscatt = 4
        call scatt_pac(ph)
        !call scatt_pac(ph)
      else
        exit
      end if

      ! call RR_test(ph)
      ! if (ph%p_flag /= 0) then
      !   exit
      ! end if
    end do



  end subroutine diffuse_k


end module diffuse_kernel

subroutine exp_3D_cart_galaxy()
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_class_imag
  use cudafor
  use diffuse_kernel
  implicit none

  integer :: i,j,k

  integer :: Nph_tot, nxim, nyim, nxg, nyg, nzg
  real(dp) :: kappa, albedo, hgg, pl, pc, sc
  real(dp) :: xmax, ymax, zmax, rimage, viewthet, viewphi

  real(dp) :: dx, dy, dz

  integer :: jcount, totscatt
  real(dp) :: rhosum, xx, yy, zz, rho, emis, dV
  real(dp) :: taueq, taupole, delta, g2
  real(dp) :: temp, diff, rand

  integer :: Nph, ncell, Nph_sum
  type(dim3) :: blocks, threads

  integer :: seq, offset, seed = 3898
  type(curandStateXORWOW) :: iseed

  integer, device :: i_d, j_d, k_d
  integer, device :: Nph_d

  namelist /diffuse_nml/ Nph_tot,kappa,albedo,hgg,pl,pc,sc, &
    & xmax,ymax,zmax,rimage,viewthet,viewphi, &
    & nxim, nyim, nxg, nyg, nzg

  print*, 'In diffuse'

  !! Read paramaters from the namelist for this experiment
  read(u_nml, nml=diffuse_nml)

  pl_d = pl
  pc_d = pc
  sc_d = sc

  print*, Nph, kappa, albedo, hgg, pl, pc, sc
  print*, xmax,ymax,zmax,rimage,viewthet,viewphi
  print*, nxim, nyim, nxg, nyg, nzg

  !allocate(im(1))

  !! Observation direction and image set up
  im%x_pix = nxim
  im%y_pix = nyim
  im%sinto = sin(viewthet*pi/180.0_dp)
  im%costo = cos(viewthet*pi/180.0_dp)
  im%sinpo = sin(viewphi*pi/180.0_dp)
  im%cospo = cos(viewphi*pi/180.0_dp)
  im%phio = atan2(im%sinpo,im%cospo)
  im%obsx = im%sinto*im%cospo
  im%obsy = im%sinto*im%sinpo
  im%obsz = im%costo

  allocate(f(nxim,nyim),q(nxim,nyim),u(nxim,nyim))
  f(:,:) = 0.0_dp ; q(:,:) = 0.0_dp ; u(:,:) = 0.0_dp
  im%fsum = sum(f(:,:)) ; im%qsum = sum(q(:,:)) ;  im%usum = sum(u(:,:))
  im%rimage = rimage

  !! Set up grid
  grid%n_x = nxg; grid%n_y = nyg; grid%n_z = nzg
  grid%x_max = xmax; grid%y_max = ymax; grid%z_max = zmax
  allocate(x(grid%n_x+1), y(grid%n_y+1),z(grid%n_z+1))
  x(:) = 0.0_dp ; y(:) = 0.0_dp ; z(:) = 0.0_dp

  dx = grid%x_max/real(grid%n_x,dp)
  do i = 1, grid%n_x+1
    x(i) = real((i-1),dp)*2.0_dp*dx
  end do
  dy = grid%y_max/real(grid%n_y,dp)
  do i = 1, grid%n_y+1
    y(i) = real((i-1),dp)*2.0_dp*dy
  end do
  dz = grid%z_max/real(grid%n_z,dp)
  do i = 1, grid%n_z+1
    z(i) = real((i-1),dp)*2.0_dp*dz
  end do


  allocate(rhokap(grid%n_x,grid%n_y,grid%n_z),emit(grid%n_x,grid%n_y,grid%n_z))
  rhokap(:,:,:) = 0.0_dp ; emit(:,:,:) = 0.0_dp

  rhosum = 0.0_dp
  grid%lumtot = 0.0_dp
  do i = 1, grid%n_x
    xx = x(i) - grid%x_max + dx
    do j = 1, grid%n_y
      yy = y(j) - grid%y_max + dy
      do k = 1, grid%n_z
        zz = z(k) - grid%z_max + dz
        call density(xx,yy,zz,rho)
        rhokap(i,j,k) = rho*kappa*3.09e21_dp ! rho*kappa*R, R=1kpc=3.09e21cm
        rhosum = rhosum + rho !rhokap(i,j,k)/(kappa*3.09e21_dp) /2.0_dp
        call setemit(xx,yy,zz,emis)
        emit(i,j,k) = emis
        grid%lumtot = grid%lumtot + emit(i,j,k)
        !print*, i, j, k, rhokap(i,j,k), emit(i,j,k)

      end do
    end do
  end do


  dV = 2.0_dp*dx
  dV = dV*2.0_dp*dy
  dV = dV*2.0_dp*dz
  dV = dV*3.09e21_dp
  dV = dV/2.0e33_dp  ! dV/M_sun
  dV = dV*3.09e21_dp
  dV = dV*3.09e21_dp
  print*, 'Mass (M_sun) = ',rhosum*dV
  print*, 'Luminosity = ',grid%lumtot

  taueq = 0.0_dp
  taupole = 0.0_dp
  do i = 1, grid%n_x
     taueq = taueq + rhokap(i,(grid%n_y)/2,(grid%n_z)/2)
  enddo
  do i = 1, grid%n_z
     taupole = taupole + rhokap((grid%n_x)/2,(grid%n_y)/2,i)
  enddo
  taueq = taueq * 2.0_dp * dx
  taupole = taupole * 2.0_dp * dz
  print*, 'taueq = ',taueq,'  taupole = ',taupole


  print*, '==================='

  seq = 0
  offset = 0

  !call curandinit(seed, seq, offset, iseed)


  !! Pass grid and image to device global memory
  im_d%x_pix = im%x_pix
  im_d%y_pix = im%y_pix
  allocate(f_d(im%x_pix,im%y_pix),q_d(im%x_pix,im%y_pix),u_d(im%x_pix,im%y_pix))
  f_d(:,:) = f(:,:)
  q_d(:,:) = q(:,:)
  u_d(:,:) = u(:,:)
  im_d = im

  allocate(x_d(grid%n_x+1),y_d(grid%n_y+1),z_d(grid%n_z+1))
  allocate(rhokap_d(grid%n_x,grid%n_y,grid%n_z))
  allocate(g_d(grid%n_x,grid%n_y,grid%n_z),ssa_d(grid%n_x,grid%n_y,grid%n_z))
  x_d(:) = x(:)
  y_d(:) = y(:)
  z_d(:) = z(:)
  rhokap_d(:,:,:) = rhokap(:,:,:)
  ssa_d(:,:,:) = albedo
  g_d(:,:,:) = hgg
  grid_d = grid


  !! Call the experiment main kernel routine
  !! This replaces the big photon loop in the CPU code
  Nph_sum = 0
  do i = 1, grid%n_x
    do j = 1, grid%n_y
      do k = 1, grid%n_z

        temp = Nph_tot*emit(i,j,k)/grid%lumtot
        diff = temp - int(temp)
        call random_number(rand)
        if (rand < diff .and. diff > 0.0_dp) then
           Nph = int(temp)+1
        else
          Nph = int(temp)
        end if

        !Nph = Nph_tot*emit(i,j,k)/grid%lumtot
        if (Nph == 0) then
          cycle
        end if

        !print*, i, j, k, Nph, temp, rand, diff

        if (Nph < 256) then
          threads = dim3(Nph, 1, 1)
          blocks = dim3(1,1,1)
        else
          threads = dim3(256, 1, 1)
          blocks = dim3(ceiling(real(Nph)/threads%x),1,1)
        end if

        Nph_d = Nph
        Nph_sum = Nph_sum + Nph
        i_d = i ; j_d = j ;  k_d = k
        call diffuse_k<<<blocks, threads>>>(Nph_d, i_d, j_d, k_d)

      end do
    end do
  end do

  im = im_d
  f(:,:) = f_d(:,:)/real(Nph_sum,dp)
  q(:,:) = q_d(:,:)/real(Nph_sum,dp)
  u(:,:) = u_d(:,:)/real(Nph_sum,dp)

  im%fsum = im%fsum/real(Nph_sum,dp)
  im%qsum = im%qsum/real(Nph_sum,dp)
  im%usum = im%usum/real(Nph_sum,dp)

  print*, 'Emitted: ', Nph_sum
  print*, 'Asked for: ', Nph_tot
  print*, 'Total f, q, u :', im%fsum, im%qsum, im%usum


  ! Output images
  open(unit=10,file='fimage.dat',status='unknown',form='unformatted')
  write(10) real(f)
  close(10)

  open(unit=10,file='qimage.dat',status='unknown',form='unformatted')
  write(10) real(q)
  close(10)

  open(unit=10,file='uimage.dat',status='unknown',form='unformatted')
  write(10) real(u)
  close(10)

  ! Output images
  open(unit=10,file='fimage_ascii.dat',action='readwrite',form='formatted')
  write(10,*)  im%fsum
  do i = 1, im%x_pix
    write(10,*) (f(j,i), j = 1, im%y_pix)
  end do
  close(10)

  open(unit=10,file='qimage_ascii.dat',action='readwrite',form='formatted')
  write(10,*)  im%qsum
  do i = 1, im%x_pix
    write(10,*) (q(j,i), j = 1, im%y_pix)
  end do
  close(10)

  open(unit=10,file='uimage_ascii.dat',action='readwrite',form='formatted')
  write(10,*)  im%fsum
  do i = 1, im%x_pix
    write(10,*) (u(j,i), j = 1, im%y_pix)
  end do
  close(10)


contains

  subroutine density(x,y,z,rho)
     implicit none

     real(dp), intent(in) :: x, y, z
     real(dp), intent(out) :: rho

     real(dp) ::  w ,w2, r, r2, h, rr

      w2=x*x+y*y
      w=sqrt(w2)
      r2=w2+z*z
      r=sqrt(r2)

      h = 0.25_dp
      rr = 5.0_dp
      if (r < 20.0_dp) then
        rho = exp(-r/rr)*exp(-abs(z)/h)*0.66e-23_dp  ! rho in g/cm^3
      else
        rho = 0.0_dp
      endif

      !print*, x, y, z, w2, w, r2, r, rho

  end subroutine density

  subroutine setemit(x,y,z,emis)
    implicit none

    real(dp), intent(in) :: x,y,z
    real(dp), intent(out) :: emis

    real(dp) ::  w,w2,r,r2,h,rr

      w2=x*x+y*y
      w=sqrt(w2)
      r2=w2+z*z
      r=sqrt(r2)

      h = 0.5_dp
      rr = 5.0_dp
      if (r < 20.0_dp) then
        emis = exp(-r/rr)*exp(-abs(z)/h)
      else
        emis = 0.0_dp
      endif

  end subroutine setemit

end subroutine exp_3D_cart_galaxy
