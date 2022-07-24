module lbl_tables_read
  use optools_data_mod
  implicit none

  private :: read_lbl_Joost, read_lbl_CMCRT
  public :: read_lbl_tables

contains

  subroutine read_lbl_tables()
    implicit none

    integer :: s

    print*, ' ~~ Reading in lbl tables ~~ '

    do s = 1, nlbl

      select case (lbl_tab(s)%form)
      case(0)
        print*, ' - Read Joost lbl table: ', s,  lbl_tab(s)%sp, lbl_tab(s)%form, lbl_tab(s)%iVMR
        call read_lbl_Joost(s)

      case(1)
        print*, ' - Reading CMCRT lbl table: ', s,  lbl_tab(s)%sp, lbl_tab(s)%form, lbl_tab(s)%iVMR
        call read_lbl_CMCRT(s)

      case(2)
        print*, ' - Reading HELIOS-K lbl table: ', s, lbl_tab(s)%sp, lbl_tab(s)%form, lbl_tab(s)%iVMR
        !call read_lbl_HELIOS(s)

      case(3)
        print*, ' - Reading SOCRATES lbl table: ', s, lbl_tab(s)%sp, lbl_tab(s)%form, lbl_tab(s)%iVMR
        !call read_lbl_SOCRATES(s)

      ! case(4)
      !   print*, 'Reading HITRAN lbl table for species: ', s, CIA_tab(s)%sp, CIA_tab(s)%form
      !   call read_lbl_HITRAN(s)

      case default
        print*, 'ERROR - lbl table format integer not valid - STOPPING'
        print*, 'Species: ', s, lbl_tab(s)%sp, lbl_tab(s)%form
        stop
      end select
    end do

    print*, ' ~~ Quest completed  ~~ '


  end subroutine read_lbl_tables

  !! Read HELIOS-K CMCRT format lbl table
  subroutine read_lbl_CMCRT(s)
    implicit none

    integer, intent(in) :: s
    integer :: np, nt, n_bins, l, u, i, j, p, t
    character(len=10) sp
    real(kind=dp), allocatable, dimension(:) :: press, temp, wn
    real(kind=dp), allocatable, dimension(:) :: k_abs

    ! Open gCMCRT's formatted file
    print*, ' - lbl - CMCRT Reading: ', s, lbl_tab(s)%sp, trim(lbl_tab(s)%path)
    open(newunit=u, file=trim(lbl_tab(s)%path),form='formatted', status='old',action='read')

    ! Read header
    read(u,*)
    read(u,*) nt, np, n_bins

    ! Give data to table object
    lbl_tab(s)%nP = np
    lbl_tab(s)%nT = nt
    lbl_tab(s)%nwl = n_bins

    ! Allocate local arrays
    allocate(press(np),temp(nt),wn(n_bins),k_abs(n_bins))

    ! Allocate table arrays
    allocate(lbl_tab(s)%T(lbl_tab(s)%nT),lbl_tab(s)%lT(lbl_tab(s)%nT))
    allocate(lbl_tab(s)%P(lbl_tab(s)%nP),lbl_tab(s)%lP(lbl_tab(s)%nP))
    allocate(lbl_tab(s)%wl(lbl_tab(s)%nwl))
    allocate(lbl_tab(s)%lk_abs(lbl_tab(s)%nwl,lbl_tab(s)%nP,lbl_tab(s)%nT))

    ! Read pressure, temperature and wavenumbers
    do t = 1, nt
      read(u,*) temp(t)
    end do
    do p = 1, np
      read(u,*) press(p)
    end do
    read(u,*) (wn(l),l = 1, lbl_tab(s)%nwl)

    print*,'Min, max T: ', temp(1), temp(nt)
    print*,'Min, max P: ', press(1), press(np)
    print*,'Min, max wn: ', wn(1), wn(n_bins)

    read(u,*)

    ! Give pressure, temperature and wavelengths to table object
    lbl_tab(s)%P(:) = press(:) * bar
    lbl_tab(s)%T(:) = temp(:)

    lbl_tab(s)%lP(:) = log10(lbl_tab(s)%P(:))
    lbl_tab(s)%lT(:) = log10(lbl_tab(s)%T(:))

    ! Reverse l index as table is in wavenumber order - convert wn to um
    do l = 1, lbl_tab(s)%nwl
      lbl_tab(s)%wl(l) = 1.0_dp/wn(n_bins-l+1) * 1e4_dp
    end do

    ! print*, lbl_tab(s)%nP, lbl_tab(s)%nT, lbl_tab(s)%nwl
    ! print*, lbl_tab(s)%P(:)/bar
    ! print*, lbl_tab(s)%T(:)
    ! print*, lbl_tab(s)%wl(:)

    ! Read in the table data
    do i = 1, lbl_tab(s)%nT
      do j = 1, lbl_tab(s)%nP
        read(u,*) (k_abs(l), l = 1, n_bins)
        ! Reverse l index as table is in wavenumber order & log values
        do l = 1, n_bins
          lbl_tab(s)%lk_abs(l,j,i) = log10(max(k_abs(n_bins-l+1),1e-99_dp))
        end do
      end do
    end do

    ! Deallocate work arrays and close units
    deallocate(press,temp,wn,k_abs)
    close(u)

  end subroutine read_lbl_CMCRT

  !! Read Joost format lbl table
  subroutine read_lbl_Joost(s)
    implicit none

    integer, intent(in) :: s
    integer :: np, nt, n_bins, l, u, i, j
    real(kind=dp), allocatable, dimension(:) :: press, temp, wl

    ! Open Joost's formatted file
    print*, ' - lbl - Joost Reading: ', s, lbl_tab(s)%sp, trim(lbl_tab(s)%path)
    open(newunit=u, file=trim(lbl_tab(s)%path), action='read')

    ! Read header
    read(u,*) np, nt, n_bins

    ! Give data to table object
    lbl_tab(s)%nP = np
    lbl_tab(s)%nT = nt
    lbl_tab(s)%nwl = n_bins

    ! Allocate local arrays
    allocate(press(np),temp(nt),wl(n_bins))

    ! Allocate table arrays
    allocate(lbl_tab(s)%T(lbl_tab(s)%nT),lbl_tab(s)%lT(lbl_tab(s)%nT))
    allocate(lbl_tab(s)%P(lbl_tab(s)%nP),lbl_tab(s)%lP(lbl_tab(s)%nP))
    allocate(lbl_tab(s)%wl(lbl_tab(s)%nwl))
    allocate(lbl_tab(s)%lk_abs(lbl_tab(s)%nwl,lbl_tab(s)%nP,lbl_tab(s)%nT))

    ! Read pressure, temperature and wavelengths
    read(u,*) press(:)
    read(u,*) temp(:)
    read(u,*) wl(:)

    ! Give pressure, temperature and wavelengths to table object
    lbl_tab(s)%P(:) = press(:) * bar
    lbl_tab(s)%T(:) = temp(:)

    lbl_tab(s)%lP(:) = log10(lbl_tab(s)%P(:))
    lbl_tab(s)%lT(:) = log10(lbl_tab(s)%T(:))

    lbl_tab(s)%wl(:) = wl(:)

    ! print*, lbl_tab(s)%nP, lbl_tab(s)%nT, lbl_tab(s)%nwl
    ! print*, lbl_tab(s)%P(:)/bar
    ! print*, lbl_tab(s)%T(:)
    ! print*, lbl_tab(s)%wl(:)

    ! Read in the table data
    do i = 1, lbl_tab(s)%nP
      do j = 1, lbl_tab(s)%nT
        read(u,*) lbl_tab(s)%k_abs(:,i,j)
        lbl_tab(s)%lk_abs(:,i,j) = log10(max(lbl_tab(s)%k_abs(:,i,j),1.0e-99_dp))
        !print*, (lbl_tab(s)%k_abs(l,i,j), l = 1, lbl_tab(s)%nwl)
      end do
    end do

    ! Deallocate work arrays and close units
    deallocate(press,temp,wl)
    close(u)

  end subroutine read_lbl_Joost


end module lbl_tables_read
