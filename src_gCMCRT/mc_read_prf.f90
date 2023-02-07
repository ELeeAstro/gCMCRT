module mc_read_prf
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  implicit none


contains

  subroutine read_1D_prf()
    implicit none

    integer :: u, i, j, k, jj, icount, nlay, nlev, ngas, n, l
    real(kind=dp) :: mu
    !real(kind=dp), allocatable, dimension(:) :: VMR_dum

    print*, ' - Reading in '//trim(exp_name)//'.prf file - '

    open(newunit=u, file=trim(exp_name)//'.prf', status='old', action='read', form='formatted')

    read(u,*); read(u,*)
    read(u,*) nlay

    nlev = nlay + 1
    nlay_1D = nlay

    allocate(lay(nlay), PG_1D(nlay), MOL_W_1D(nlay), TG_1D(nlay), &
    & RH_1D(nlay))

    read(u,*)
    read(u,*) ngas

    !allocate(VMR_1D(ngas,nlay),  VMR_dum(ngas))

    do i = 1, ngas
      read(u,*)
      !print*, g_name(i)
    end do

    read(u,*); read(u,*)
    do i = 1, nlay
      read(u,*) lay(i), PG_1D(i), TG_1D(i), MOL_W_1D(i)!, (VMR_dum(j), j=1,ngas)
      !VMR_1D(:,i) = VMR_dum(:)
      PG_1D(i) = PG_1D(i) * bar
      RH_1D(i) = (PG_1D(i) * MOL_W_1D(i) * amu) / (kb * TG_1D(i))
    end do

    close(u)

    print*, ' - Complete - '

    allocate(RH(grid%n_lay,grid%n_phi-1,grid%n_theta-1), TG(grid%n_lay,grid%n_phi-1,grid%n_theta-1), &
      & PG(grid%n_lay,grid%n_phi-1,grid%n_theta-1), MOL_W(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
    !allocate(VMR(ngas,grid%n_lay,grid%n_phi-1,grid%n_theta-1))

    if (oneD .eqv. .True. ) then
      do i = 1, nlay
        RH(i,:,:) = RH_1D(i)
        TG(i,:,:) = TG_1D(i)
        PG(i,:,:) = PG_1D(i)
        MOL_W(i,:,:) = MOL_W_1D(i)
        !do j = 1, ngas
        !  VMR(j,i,:,:) = VMR_1D(j,i)
        !end do
      end do
    else if (threeD .eqv. .True.) then
      ! Unpack the 1D profile to the 3D grid
      n = 1
      do k = 1, grid%n_theta-1
        do j = 1, grid%n_phi-1
          do i = 1, grid%n_lay
            RH(i,j,k) = RH_1D(n)
            TG(i,j,k) = TG_1D(n)
            PG(i,j,k) = PG_1D(n)
            MOL_W(i,j,k) = MOL_W_1D(n)
            !do l = 1, ngas
            !  VMR(l,i,j,k) = VMR_1D(l,n)
            !end do
            n = n + 1
          end do
        end do
      end do

    end if

    open(newunit=u,file=trim(exp_name)//'.hprf',action='read',status='old',form='formatted')
    allocate(H(grid%n_lev))
    do i = 1, grid%n_lev
      read(u,*) lay(i), H(i)
    end do

    allocate(H_d(grid%n_lev))
    H_d(:) = H(:)

    ! 1D arrays have served their purpose, deallocate to save memory
    deallocate(lay,PG_1D,MOL_W_1D,TG_1D,RH_1D)
    !deallocate(VMR_dum,VMR_1D)

    close(u)

  end subroutine read_1D_prf

  subroutine read_1D_wprf()
    implicit none

    integer :: u, i, j, k, dum, n

    print*, ' - Reading in '//trim(exp_name)//'.wprf file - '

    open(newunit=u, file=trim(exp_name)//'.wprf', status='old', action='read', form='formatted')

    allocate(u_wind_1D(nlay_1D),v_wind_1D(nlay_1D),w_wind_1D(nlay_1D))

    do i = 1, nlay_1D
      read(u,*) dum, u_wind_1D(i), v_wind_1D(i), w_wind_1D(i)
    end do

    allocate(u_wind(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
    allocate(v_wind(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
    allocate(w_wind(grid%n_lay,grid%n_phi-1,grid%n_theta-1))

    if (oneD .eqv. .True. ) then
      do i = 1, grid%n_lay
        u_wind(i,:,:) = u_wind_1D(i)
        v_wind(i,:,:) = v_wind_1D(i)
        w_wind(i,:,:) = w_wind_1D(i)
      end do
    else if (threeD .eqv. .True.) then
      ! Unpack the 1D profile to the 3D grid
      n = 1
      do k = 1, grid%n_theta-1
        do j = 1, grid%n_phi-1
          do i = 1, grid%n_lay
            u_wind(i,j,k) = u_wind_1D(n)
            v_wind(i,j,k) = v_wind_1D(n)
            w_wind(i,j,k) = w_wind_1D(n)
            n = n + 1
          end do
        end do
      end do
    end if

    close(u)

    deallocate(u_wind_1D,v_wind_1D,w_wind_1D)

    print*, ' - Complete - '


  end subroutine read_1D_wprf

  subroutine read_wl()
    implicit none

    integer :: iwl, num_wl, l, uwl

    open(newunit=uwl,file='wavelengths.wl',action='read',status='old',form='formatted')

    read(uwl,*) num_wl
    allocate(wl(num_wl))
    do l = 1, num_wl
      read(uwl,*) iwl, wl(l)
    end do

    ! Calculate the wavelength band edges
    allocate(wl_e(num_wl+1))
    wl_e(1) = 0.3_dp
    wl_e(num_wl+1) = 30.0_dp
    do l = 2, num_wl
      wl_e(l) = (wl(l-1) + wl(l))/2.0_dp
    end do

    !! Send wl array to gpu memory
    allocate(wl_d(num_wl))
    wl_d(:) = wl(:)

    close(uwl)

  end subroutine read_wl

  subroutine read_g_ord()
    implicit none

    integer :: u, g

    ! Read g-ordinances and delg with weights to file g.ord
    open(newunit=u, file='gord.cmcrt', status='old', action='read',form='formatted')

    read(u,*) ng

    allocate(gord_cdf(ng), gord_x(ng), gord_w(ng))
    do g = 1, ng
      read(u,*) gord_x(g), gord_w(g)
      print*, gord_x(g), gord_w(g)
    end do

    close(u)

    print*, 'Constructing g_ord_CDF'

    gord_cdf(1) = 0.0_dp
    do g = 1, ng-1
      gord_cdf(g+1) = gord_cdf(g) + gord_w(g)
    end do

    ! Send data to GPU
    allocate(gord_cdf_d(ng), gord_x_d(ng), gord_w_d(ng))
    ng_d = ng
    gord_cdf_d(:) = gord_cdf(:)
    gord_x_d(:) = gord_x(:)
    gord_w_d(:) = gord_w(:)

    print*, '- Complete -'

  end subroutine read_g_ord


end module mc_read_prf
