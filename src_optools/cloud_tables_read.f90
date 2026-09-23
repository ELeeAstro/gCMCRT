module cloud_tables_read
  use optools_data_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
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

    ! Open DIHRT formatted nk file
    print*, ' - cl - DIHRT Reading: ', s, cl_tab(s)%sp, trim(cl_tab(s)%path)
    open(newunit=u, file=trim(cl_tab(s)%path), form='formatted', status='old', action='read')

    ! Read number of lines and conducting flag
    read(u,*) nlines, conducting_flag

    if (conducting_flag .and. nlines < 4) then
      print*, 'ERROR - Conducting nk table requires at least 4 wavelength points - STOPPING'
      print*, 'Species, path, points: ', cl_tab(s)%sp, trim(cl_tab(s)%path), nlines
      stop
    end if


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

      if (.not. ieee_is_finite(cl_tab(s)%wl(l)) .or. cl_tab(s)%wl(l) <= 0.0_dp) then
        print*, 'ERROR - Cloud nk wavelengths must be finite and positive - STOPPING'
        print*, 'Species, path, point, wavelength: ', cl_tab(s)%sp, trim(cl_tab(s)%path), l, cl_tab(s)%wl(l)
        stop
      end if

      if (l > 1) then
        if (cl_tab(s)%wl(l) <= cl_tab(s)%wl(l-1)) then
          print*, 'ERROR - Cloud nk wavelengths must be strictly increasing - STOPPING'
          print*, 'Species, path, points: ', cl_tab(s)%sp, trim(cl_tab(s)%path), l-1, l
          print*, 'Wavelengths: ', cl_tab(s)%wl(l-1), cl_tab(s)%wl(l)
          stop
        end if
      end if

      cl_tab(s)%n(l) = max(1.0e-12_dp, cl_tab(s)%n(l))
      cl_tab(s)%k(l) = max(1.0e-12_dp, cl_tab(s)%k(l))
      !print*, l, cl_tab(s)%wl(l), cl_tab(s)%n(l),cl_tab(s)%k(l)
    end do

    close(u)

  end subroutine read_nk_DIHRT


end module cloud_tables_read
