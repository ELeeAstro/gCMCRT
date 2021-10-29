module cloud_tables_read
  use optools_data_mod
  implicit none



  private :: read_nk_DIHRT
  public :: read_cl_nk

contains

  subroutine read_cl_nk()
    implicit none

    integer :: s

    print*, ' ~~ Reading in nk tables ~~'

    do s = 1, ncl

      select case (cl_tab(s)%form)

      case(5)
        print*, ' - Reading DIHRT nk table: ', cl_tab(s)%sp, cl_tab(s)%form, cl_tab(s)%iVMR
        call read_nk_DIHRT(s)
      case default
        print*, 'ERROR - nk table format integer not valid - STOPPING'
        print*, 'Species: ', s, cl_tab(s)%sp, cl_tab(s)%form
        stop
      end select
    end do

    print*, ' ~~ Quest completed  ~~ '

  end subroutine read_cl_nk

  subroutine read_nk_DIHRT(s)
    implicit none

    integer, intent(in) :: s
    integer :: l, u, nlines
    logical :: conducting_flag

    ! Open DOHRT formatted nk file
    print*, ' - cl - DIHRT Reading: ', s, cl_tab(s)%sp, trim(cl_tab(s)%path)
    open(newunit=u, file=trim(cl_tab(s)%path), form='formatted', status='old', action='read')

    ! Read number of lines and conducting flag
    read(u,*) nlines, conducting_flag

    ! Put data into cl_tab container
    cl_tab(s)%nwl = nlines
    allocate(cl_tab(s)%wl(cl_tab(s)%nwl))
    allocate(cl_tab(s)%n(cl_tab(s)%nwl),cl_tab(s)%k(cl_tab(s)%nwl))

    cl_tab(s)%conducting = conducting_flag

    ! Read 4 blank lines
    read(u,*) ; read(u,*); read(u,*) ; read(u,*)

    ! Read wavelengths, n and k values
    do l = 1, cl_tab(s)%nwl
       read(u,*) cl_tab(s)%wl(l), cl_tab(s)%n(l),cl_tab(s)%k(l)
        cl_tab(s)%n(l) = max(0.0_dp, cl_tab(s)%n(l))
        cl_tab(s)%k(l) = max(0.0_dp, cl_tab(s)%k(l))
        !print*, l, cl_tab(s)%wl(l), cl_tab(s)%n(l),cl_tab(s)%k(l)
    end do

  end subroutine read_nk_DIHRT


end module cloud_tables_read
