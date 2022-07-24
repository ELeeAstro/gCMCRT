module CK_tables_read
  use optools_data_mod
  implicit none

  private :: read_CK_NEMESIS!, read_CK_HELIOS
  public :: read_CK_tables

contains

  subroutine read_CK_tables(pre_mixed)
    implicit none

    logical, intent(in) :: pre_mixed
    integer :: s


    print*, ' ~~ Reading in CK tables ~~ '

    do s = 1, nCK

      select case (CK_tab(s)%form)
      case(0)
        print*, ' - Skipping special species: ', s,  CK_tab(s)%sp, CK_tab(s)%form, CK_tab(s)%iVMR

      case(1)
        print*, ' - Reading NEMESIS CK table: ', s,  CK_tab(s)%sp, CK_tab(s)%form, CK_tab(s)%iVMR
        call read_CK_NEMESIS(s)
      case(2)

        print*, ' - Reading CMCRT CK table: ', s, CK_tab(s)%sp, CK_tab(s)%form, CK_tab(s)%iVMR
        call read_CK_CMCRT(s,pre_mixed)

      case(3)
        print*, ' - Reading SOCRATES CK table: ', s, CK_tab(s)%sp, CK_tab(s)%form, CK_tab(s)%iVMR
        !call read_CK_SOCRATES(s)

      ! case(4)
      !   print*, 'Reading HITRAN CK table for species: ', s, CIA_tab(s)%sp, CIA_tab(s)%form
      !   call read_CK_HITRAN(s)

      case default
        print*, 'ERROR - CK table format integer not valid - STOPPING'
        print*, 'Species: ', s, CK_tab(s)%sp, CK_tab(s)%form
        stop
      end select
    end do

    print*, ' ~~ Quest completed  ~~ '

  end subroutine read_CK_tables

  !! Reads a CMCRT formatted ascii CK table
  subroutine read_CK_CMCRT(s,pre_mixed)
    implicit none

    integer, intent(in) :: s
    logical, intent(in) :: pre_mixed

    integer :: np, nt, n_bins, ng,  l, u, i, j, p, t, g
    character(len=10) sp
    real(kind=dp), allocatable, dimension(:) :: press, temp, wll, gx, gw
    real(kind=dp), allocatable, dimension(:) :: k_abs

    ! Open CMCRT formatted file
    print*, ' - CK - CMCRT format reading: ', s, CK_tab(s)%sp, trim(CK_tab(s)%path)
    open(newunit=u, file=trim(CK_tab(s)%path),form='formatted', status='old',action='read')

    ! Read header
    read(u,*)
    read(u,*) nt, np, n_bins, ng

    ! Give data to table object
    CK_tab(s)%nT = nt
    CK_tab(s)%nP = np
    CK_tab(s)%nwl = n_bins
    CK_tab(s)%nG = ng

    ! Allocate local arrays
    allocate(press(np),temp(nt),wll(n_bins),k_abs(ng),gx(ng),gw(ng))

    ! Allocate table arrays
    allocate(CK_tab(s)%T(CK_tab(s)%nT),CK_tab(s)%lT(CK_tab(s)%nT))
    allocate(CK_tab(s)%P(CK_tab(s)%nP),CK_tab(s)%lP(CK_tab(s)%nP))
    allocate(CK_tab(s)%wl(CK_tab(s)%nwl),CK_tab(s)%wn(CK_tab(s)%nwl))
    allocate(CK_tab(s)%Gx(CK_tab(s)%nG),CK_tab(s)%Gw(CK_tab(s)%nG))
    allocate(CK_tab(s)%k_abs(CK_tab(s)%nwl,CK_tab(s)%nP,CK_tab(s)%nT,CK_tab(s)%nG))
    allocate(CK_tab(s)%lk_abs(CK_tab(s)%nwl,CK_tab(s)%nP,CK_tab(s)%nT,CK_tab(s)%nG))


    ! Read pressure, temperature and wavenumbers
    do t = 1, nt
      read(u,*) temp(t)
    end do
    do p = 1, np
      read(u,*) press(p)
    end do
    read(u,*) (wll(l),l = 1, CK_tab(s)%nwl)
    read(u,*) (gx(g),g = 1, CK_tab(s)%nG)
    read(u,*) (gw(g),g = 1, CK_tab(s)%nG)

    print*,'Min, max T: ', temp(1), temp(nt)
    print*,'Min, max P: ', press(1), press(np)
    print*,'Min, max wl: ', wll(1), wll(n_bins)
    print*,'1, ng gx: ', gx(1), gx(ng)
    print*,'1, ng gw: ', gw(1), gw(ng)

    read(u,*)

    ! Give pressure, temperature and wavelengths to table object
    CK_tab(s)%T(:) = temp(:)
    CK_tab(s)%P(:) = press(:) * bar

    CK_tab(s)%lT(:) = log10(CK_tab(s)%T(:))
    CK_tab(s)%lP(:) = log10(CK_tab(s)%P(:))

    ! Give g ordinance x and weights to object
    CK_tab(s)%Gx(:) = gx(:)
    CK_tab(s)%Gw(:) = gw(:)


    ! Reverse l index as table is in wavenumber order
    do l = 1, CK_tab(s)%nwl
      CK_tab(s)%wl(l) = wll(l)
    end do

    ! Add wavenumber information to k table [cm-1]
    CK_tab(s)%wn(:) = 1.0_dp/(CK_tab(s)%wl(:) * 1.0e-4_dp)

    ! print*, CK_tab(s)%nP, CK_tab(s)%nT, CK_tab(s)%nwl
    ! print*, CK_tab(s)%P(:)/bar
    ! print*, CK_tab(s)%T(:)
    ! print*, CK_tab(s)%wl(:)

    ! Read in the table data - the CMCRT format from the HELIOS-K data is a bit weird, wavelength is in the inner loop
    if (pre_mixed .eqv. .True.) then
      ! The premixed table goes as pressure first due to GGchem constraints
      do j = 1, CK_tab(s)%nP
        do i = 1, CK_tab(s)%nT
          do l = CK_tab(s)%nwl, 1, -1
            read(u,*) (k_abs(g), g = 1, CK_tab(s)%nG)
            ! Reverse l index as table is in wavenumber order from HELIOS-K
            CK_tab(s)%k_abs(l,j,i,:) = max(k_abs(:),1.0e-99_dp)
            CK_tab(s)%lk_abs(l,j,i,:) = log10(CK_tab(s)%k_abs(l,j,i,:))
            !print*, l,j,i, CK_tab(s)%k_abs(l,j,i,1),CK_tab(s)%k_abs(l,j,i,CK_tab(s)%nG)
          end do
        end do
      end do
    else
      do i = 1, CK_tab(s)%nT
        do j = 1, CK_tab(s)%nP
          do l = CK_tab(s)%nwl, 1, -1
            read(u,*) (k_abs(g), g = 1, CK_tab(s)%nG)
            ! Reverse l index as table is in wavenumber order from HELIOS-K
            CK_tab(s)%k_abs(l,j,i,:) = max(k_abs(:),1.0e-99_dp)
            CK_tab(s)%lk_abs(l,j,i,:) = log10(CK_tab(s)%k_abs(l,j,i,:))
            !print*, l,j,i, CK_tab(s)%k_abs(l,j,i,1),CK_tab(s)%k_abs(l,j,i,CK_tab(s)%nG)
          end do
        end do
      end do
    end if

    ! Deallocate work arrays and close units
    deallocate(press,temp,wll,k_abs,gx,gw)
    close(u)


  end subroutine read_CK_CMCRT


  !! Read NEMESIS format CK table
  subroutine read_CK_NEMESIS(s)
    implicit none

    integer, intent(in) :: s
    integer :: u, i, j, g, l
    integer :: irec, nrec, jrec
    integer :: irec0, idgas, isogas
    integer :: np, nt, ng, npoint
    real(kind=sp) :: vmin, delv, fwhm
    real(kind=sp) :: g1_dum, g2_dum, p_dum, t_dum, k_dum
    real(kind=sp) :: v_dum, vmax

    ! Open NEMESIS .kta file
    print*, ' - CK - NEMESIS Reading: ', s, CK_tab(s)%sp, trim(CK_tab(s)%path)
    open(newunit=u, file=trim(CK_tab(s)%path), status='old', access='direct', action='read',recl=4)

    ! Read header information
    read(u, rec=1) irec0
    read(u, rec=2) npoint
    read(u, rec=3) vmin
    read(u, rec=4) delv
    read(u, rec=5) fwhm
    read(u, rec=6) np
    read(u, rec=7) nt
    read(u, rec=8) ng
    read(u, rec=9) idgas
    read(u, rec=10) isogas

    ! Print header
    print*, ' - CK - NEMESIS Header: ', irec0, npoint, vmin, delv, fwhm, np, nt, ng ,idgas, isogas

    ! Allocate dimensions to CK_tab object arrays
    CK_tab(s)%nwl = npoint
    allocate(CK_tab(s)%wl(CK_tab(s)%nwl),CK_tab(s)%wn(CK_tab(s)%nwl))
    CK_tab(s)%nP = np
    allocate(CK_tab(s)%P(CK_tab(s)%nP), CK_tab(s)%lP(CK_tab(s)%nP))
    CK_tab(s)%nT = nt
    allocate(CK_tab(s)%T(CK_tab(s)%nT),CK_tab(s)%lT(CK_tab(s)%nT))
    CK_tab(s)%nG = ng
    allocate(CK_tab(s)%Gx(CK_tab(s)%nG),CK_tab(s)%Gw(CK_tab(s)%nG))
    allocate(CK_tab(s)%k_abs(CK_tab(s)%nwl,CK_tab(s)%nP,CK_tab(s)%nT,CK_tab(s)%nG))
    allocate(CK_tab(s)%lk_abs(CK_tab(s)%nwl,CK_tab(s)%nP,CK_tab(s)%nT,CK_tab(s)%nG))

    !! Start reading data from NEMESIS table

    ! Wavelength grid for now (TODO: use NEMESIS code to check for wn or wl grid in k_table - Ryan + Katy use wl grid.
    vmax = vmin + (npoint - 1)*delv
    CK_tab(s)%wl(1) = real(vmin, kind=dp)
    do l = 2, CK_tab(s)%nwl
      CK_tab(s)%wl(l) = CK_tab(s)%wl(l-1) + real(delv,kind=dp)
    end do

    irec = 11

    ! Read in g-ordinances
    do g = 1, CK_tab(s)%nG
      read(u, rec=irec) g1_dum
      CK_tab(s)%Gx(g) = real(g1_dum, kind=dp)
      !print*, g, CK_tab(s)%Gx(g)
      irec = irec + 1
    end do

    ! Read in g-weights
    do g = 1, CK_tab(s)%nG
      read(u, rec=irec) g2_dum
      CK_tab(s)%Gw(g) = real(g2_dum, kind=dp)
      !print*, g, CK_tab(s)%Gw(g)
      irec = irec + 1
    end do

    irec = 11 + 2*CK_tab(s)%nG + 2

    ! Read in pressure points
    do j = 1, CK_tab(s)%nP
      read(u, rec=irec) p_dum
      CK_tab(s)%P(j) = real(p_dum,kind=dp) * atm ! Convert atm to dyne
      CK_tab(s)%lP(j) = log10(CK_tab(s)%P(j))
      !print*, j, CK_tab(s)%P(j) / bar, log10(CK_tab(s)%P(j)/bar)
      irec = irec + 1
    end do

    ! Read in temperature points
    do j = 1, CK_tab(s)%nT
      read(u, rec=irec) t_dum
      CK_tab(s)%T(j) = real(t_dum,kind=dp)
      CK_tab(s)%lT(j) = log10(CK_tab(s)%T(j))
      !print*, j, CK_tab(s)%T(j)
      irec = irec + 1
    end do

    ! Check if custion wavelength grid is to read
    if (delv < 0.0_dp) then
      do l = 1, CK_tab(s)%nwl
        read(u,rec=irec) v_dum
        CK_tab(s)%wl(l) = real(v_dum,kind=dp)
        !print*, l, CK_tab(s)%wl(l)
        irec = irec + 1
      end do
    end if

    ! Add wavenumber information to k table [cm-1]
    CK_tab(s)%wn(:) = 1.0_dp/(CK_tab(s)%wl(:) * 1.0e-4_dp)

    print*,'Min, max T: ', CK_tab(s)%T(1), CK_tab(s)%T(CK_tab(s)%nT)
    print*,'Min, max P: ', CK_tab(s)%P(1)/1.0e6_dp, CK_tab(s)%P(CK_tab(s)%nP)/1.0e6_dp
    print*,'Min, max wl: ', CK_tab(s)%wl(1), CK_tab(s)%wl(CK_tab(s)%nwl)
    print*,'1, ng gx: ', CK_tab(s)%Gx(1), CK_tab(s)%Gx(CK_tab(s)%nG)
    print*,'1, ng gw: ', CK_tab(s)%Gw(1), CK_tab(s)%Gw(CK_tab(s)%nG)

    ! Some record keeping
    if (delv < 0.0_dp) then
      nrec = 11 + 2*ng + 2 + np + nt + npoint ! number of records read so far
    else
      nrec = 11 + 2*ng + 2 + np + nt
    end if
    jrec = irec0 - nrec      ! number of records to skip to get to first k-table record
    irec = nrec
    !print*, irec, irec0, nrec, jrec

    ! Skip some records
    do j = 1, jrec
      read(u, rec=irec)
      irec = irec + 1
    end do

    ! Read in the k values
    do l = 1, CK_tab(s)%nwl
      do j = 1, CK_tab(s)%nP
        do i = 1, CK_tab(s)%nT
          do g = 1, CK_tab(s)%nG
            read(u,rec=irec) k_dum
            ! De-scale 1e20 NEMESIS scaling factor and check maximum (i.e. no zeros) units [cm2 molecule-1]
            CK_tab(s)%k_abs(l,j,i,g) = max(real(k_dum, kind=dp) * 1.0e-20_dp,1.0e-99_dp)
            CK_tab(s)%lk_abs(l,j,i,g) = log10(CK_tab(s)%k_abs(l,j,i,g))
            irec = irec + 1
          end do
        end do
      end do
    end do

    print*, s, ' - Complete - '

    close(u)

  end subroutine read_CK_NEMESIS

end module CK_tables_read
