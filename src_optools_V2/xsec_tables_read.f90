module xsec_tables_read
  use optools_data_mod
  implicit none

  private :: read_xsec_VULCAN, read_xsec_PSG
  public :: read_xsec_tables

contains

  subroutine read_xsec_tables()
    implicit none

    integer :: s

    print*, ' ~~ Reading in xsec tables ~~ '

    do s = 1, nxsec
      select case (xsec_tab(s)%form)
      case(1)
        print*, ' - Reading VULCAN xsec table: ', s,  xsec_tab(s)%sp, xsec_tab(s)%form, xsec_tab(s)%iVMR
        call read_xsec_VULCAN(s)
      case(2)
        print*, ' - Reading PSG xsec table: ', s,  xsec_tab(s)%sp, xsec_tab(s)%form, xsec_tab(s)%iVMR
        call read_xsec_PSG(s)
      case default
        print*, 'ERROR - xsec table format integer not valid - STOPPING'
        print*, 'Species: ', s, xsec_tab(s)%sp, xsec_tab(s)%form
        stop
      end select
    end do

  end subroutine read_xsec_tables

  subroutine read_xsec_VULCAN(s)
    implicit none

    integer, intent(in) :: s

    integer :: u, io, nwl, l
    real(dp) :: xwl, cross_sec

    ! Open VULCAN xsec formatted file
    print*, ' - xsec - VULCAN Reading: ', s, xsec_tab(s)%sp, trim(xsec_tab(s)%path)
    open(newunit=u, file=trim(xsec_tab(s)%path),form='formatted', status='old',action='read')

    ! Read header
    read(u,*)

    ! Find number of lines in file
    nwl = 0 
    do 
      read(u,*,iostat=io)
      if (is_iostat_end(io) .eqv. .True.) then
        exit
      else
        nwl = nwl + 1
      end if
    end do

    xsec_tab(s)%nwl = nwl

    allocate(xsec_tab(s)%wl(xsec_tab(s)%nwl))
    allocate(xsec_tab(s)%lx_abs(xsec_tab(s)%nwl))

    ! Read in xsec data
    rewind(u)
    read(u,*)
    do l = 1, xsec_tab(s)%nwl
      read(u,*) xwl, cross_sec
      xsec_tab(s)%wl(l) = xwl * 1e-3_dp
      xsec_tab(s)%lx_abs(l) = log10(max(cross_sec,1.0e-99_dp))
      !print*, xsec_tab(s)%sp, l, xsec_tab(s)%wl(l), xsec_tab(s)%lx_abs(l)
    end do

  end subroutine read_xsec_VULCAN

 subroutine read_xsec_PSG(s)
    implicit none

    integer, intent(in) :: s

    integer :: u, io, nwl, l
    real(dp) :: xwl, cross_sec

    ! Open PSG xsec formatted file
    print*, ' - xsec - PSG Reading: ', s, xsec_tab(s)%sp, trim(xsec_tab(s)%path)
    open(newunit=u, file=trim(xsec_tab(s)%path),form='formatted', status='old',action='read')

    ! Read header
    read(u,*); read(u,*); read(u,*); read(u,*); read(u,*);

    ! Find number of lines in file
    nwl = 0
    do
      read(u,*,iostat=io)
      if (is_iostat_end(io) .eqv. .True.) then
        exit
      else
        nwl = nwl + 1
      end if
    end do

    xsec_tab(s)%nwl = nwl

    allocate(xsec_tab(s)%wl(xsec_tab(s)%nwl))
    allocate(xsec_tab(s)%lx_abs(xsec_tab(s)%nwl))

    ! Read in xsec data
    rewind(u)
    read(u,*); read(u,*); read(u,*); read(u,*); read(u,*);
    do l = 1, xsec_tab(s)%nwl
      read(u,*) xwl, cross_sec
      xsec_tab(s)%wl(l) = xwl
      xsec_tab(s)%lx_abs(l) = log10(max(cross_sec,1.0e-99_dp))
      !print*, xsec_tab(s)%sp, l, xsec_tab(s)%wl(l), xsec_tab(s)%lx_abs(l)
    end do

  end subroutine read_xsec_PSG


end module xsec_tables_read
