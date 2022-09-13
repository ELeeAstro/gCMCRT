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
  use cudafor
  implicit none

  integer :: Nph_tot, Nph_sum,  Nph, l, Nph_pad, n_lay, iscat
  character (len=8) :: fmt
  character (len=3) :: n_str
  integer, allocatable, dimension(:) :: uT
  integer, device :: Nph_pad_d
  integer :: i, j, k, n, nn, s_wl
  integer, device :: Nph_sum_d, l_d
  integer :: n_theta, n_phi, istat
  real(dp) :: viewthet
  real(dp), allocatable, dimension(:) :: viewphi
  real(dp) :: pl, pc, sc
  real(dp) :: diff, temp, rand, temp2, xi_emb
  real(dp),allocatable,dimension(:) :: em_out

  integer :: id
  integer,allocatable,dimension(:,:,:) :: Nph_cell

  type(dim3) :: blocks, threads


  namelist /sph_3D_em/ Nph_tot, s_wl, n_wl, pl, pc, sc, n_theta, n_phi, viewthet, viewphi, n_lay, xi_emb, iscat

  allocate(uT(n_phase),viewphi(n_phase))

  read(u_nml, nml=sph_3D_em)

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

  pl_d = pl
  pc_d = pc
  sc_d = sc
  iscat_d = iscat

  Nph_pad = int(real(Nph_tot*1.10_dp,dp))
  threads = dim3(128, 1, 1)
  blocks = dim3(ceiling(real(Nph_pad,dp)/threads%x),1,1)
  allocate(iseed(Nph_pad))
  Nph_pad_d = Nph_pad
  call set_iseed<<<blocks, threads>>>(Nph_pad_d)

  call read_1D_prf()
  call read_wl()
  call read_g_ord()

  call set_grid()

  ! Send data to GPU data containers
  grid_d = grid

  allocate(Nph_cell(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  allocate(wght_start(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  allocate(wght_start_d(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
  allocate(em_out(n_wl))

  do n = 1, n_phase
    write(n_str,fmt) n
    open(newunit=uT(n),file='Em_'//trim(n_str)//'.txt',action='readwrite')
    write(uT(n),*) n_wl, H(1), H(grid%n_lev), viewphi(n)
    call flush(uT(n))
  end do

  call read_next_opac(s_wl)

  print*, 'starting loop'

  do l = s_wl, n_wl

    call set_grid_opac()
    call set_grid_em(l)

    do n = 1, n_phase

      if (cmd_vphi .eqv. .False.) then
        im%vphi = viewphi(n)
      else
        viewphi(n) = im%vphi
      end if
      im%vtheta = viewthet

      call set_image()

      do k = 1, grid%n_theta-1
        do j = 1, grid%n_phi-1
          do i = 1, grid%n_lay

            if (l_cell(i,j,k) == 0.0_dp) then
              Nph_cell(i,j,k) = 0
              wght_start(i,j,k) = 0.0_dp
              cycle
            end if

            temp = real(Nph_tot,dp) * (1.0_dp - xi_emb) * l_cell(i,j,k)/grid%lumtot
            temp2 = real(Nph_tot,dp) * (xi_emb / n_cells)
            diff = (temp + temp2) - int(temp + temp2)

            call random_number(rand)
            if (rand < diff .and. diff > 0.0_dp) then
              Nph = int(temp+temp2)+1
            else
              Nph = int(temp+temp2)
            end if

            !Nph = Nph_tot*int(l_cell(i,j,k)/grid%lumtot)
            Nph_cell(i,j,k) = Nph

            wght_start(i,j,k) = 1.0_dp / ((1.0_dp - xi_emb) + &
            & xi_emb*(grid%lumtot/n_cells/l_cell(i,j,k)))

          end do
        end do
      end do

      Nph_sum = sum(Nph_cell(:,:,:))

      if (Nph_sum < 128) then
        threads = dim3(Nph_sum, 1, 1)
        blocks = dim3(1,1,1)
      else
        threads = dim3(128, 1, 1)
        blocks = dim3(ceiling(real(Nph_sum,dp)/threads%x),1,1)
      end if

      !! Now for optimisation we have to 'flatten' the array, by assosiating each packet with a cell
      allocate(Nph_i(Nph_sum),Nph_j(Nph_sum),Nph_k(Nph_sum))
      id = 1
      do k = 1, grid%n_theta-1
        do j = 1, grid%n_phi-1
          do i = 1, grid%n_lay
            do nn = 1, Nph_cell(i,j,k)
              Nph_i(id) = i
              Nph_j(id) = j
              Nph_k(id) = k
              !print*, id, Nph_sum , n, Nph_cell(i,j,k), Nph_i(id),Nph_j(id),Nph_k(id)
              id = id + 1
            end do
          end do
        end do
      end do


      im%fsum = 0.0_dp
      im%qsum = 0.0_dp
      im%usum = 0.0_dp
      im%fail_pscat = 0
      im%fail_pemit = 0

      im_d = im

      nscat_tot = 0
      nscat_tot_d = nscat_tot

      f(:,:) = 0.0_dp ; q(:,:) = 0.0_dp ; u(:,:) = 0.0_dp ; im_err(:,:) = 0.0_dp
      f_d(:,:) = f(:,:) ; q_d(:,:) = q(:,:) ; u_d(:,:) = u(:,:) ; im_err_d(:,:) = im_err(:,:)

      l_d = l
      Nph_sum_d = Nph_sum
      if (Nph_sum > 0) then
        allocate(Nph_i_d(Nph_sum),Nph_j_d(Nph_sum),Nph_k_d(Nph_sum))
        Nph_i_d(:) = Nph_i(:)
        Nph_j_d(:) = Nph_j(:)
        Nph_k_d(:) = Nph_k(:)
        wght_start_d(:,:,:) = wght_start(:,:,:)
        call exp_3D_sph_atm_em_k<<<blocks, threads>>>(l_d,Nph_sum_d)
      end if

      if (n == n_phase) then
        call read_next_opac(l+1)
      end if

      !istat = cudaDeviceSynchronize()

      im = im_d
      nscat_tot = nscat_tot_d

      em_out(l) = im%fsum / real(Nph_sum,dp)

      write(uT(n),*) wl(l), em_out(l), grid%lumtot
      !call flush(uT)

      if (do_cf .eqv. .True.) then
        cf(:,:,:) = cf_d(:,:,:)
        call output_cf(n,l)
        cf_d(:,:,:) = 0.0_dp
      end if

      if (do_images .eqv. .True.) then
        f(:,:) = f_d(:,:) ; q(:,:) = q_d(:,:) ; u(:,:) = u_d(:,:) ; im_err(:,:) = im_err_d(:,:)
        call output_im(n,l)
      end if


      print*, n, l, real(wl(l)), Nph_tot, Nph_sum, real(em_out(l)), grid%lumtot
      print*, n, 'pemit, pscat failures and nscat_tot: ',  im%fail_pemit, im%fail_pscat, nscat_tot


      deallocate(Nph_i,Nph_j,Nph_k)
      deallocate(Nph_i_d,Nph_j_d,Nph_k_d)

    end do
  
  end do


end subroutine exp_3D_sph_atm_em
