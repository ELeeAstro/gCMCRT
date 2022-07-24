module CIA_tables_mod
  use optools_data_mod
  use optools_table_class
  use CIA_tables_read, only : read_CIA_tables
  use CIA_tables_interp, only : interp_CIA_tables, interp_CIA_tables_Bezier
  implicit none

  logical :: first_call = .True.

  real(kind=dp), allocatable, dimension(:) :: CIA_out
  real(kind=sp), allocatable, dimension(:) :: CIA_write

  ! Namelist variables
  integer :: iopts
  integer, allocatable, dimension(:) :: form
  character(len=150), allocatable, dimension(:) :: paths

  namelist /CIA_nml/ iopts, form, paths

  private :: find_CIA_consituents, output_CIA_table
  public :: calc_CIA_table

contains

  !! Driver routine for CIA table calculation

  subroutine calc_CIA_table()
    implicit none

    integer :: s, i, j, l, z, ni
    logical :: exists
    real(kind=dp) :: CIA_work

    CIA_work = 0.0_dp

    !Note, the order the array allocations is important

    ! Allocate number of CIA tables
    allocate(CIA_tab(nCIA))

    ! Allocate required number of arrays from namelist options
    allocate(form(nCIA),paths(nCIA))

    ! Read CIA namelist parameters
    read(u_nml, nml=CIA_nml)

    ! Allocate work arrays
    allocate(CIA_out(nlay),CIA_write(nlay))

    ! Give the classes some global data from par and namelists
    CIA_tab(:)%sp = CIA_name(:)
    CIA_tab(:)%form = form(:)
    CIA_tab(:)%path = paths(:)

    ! Find the CIA constituents from lookup table
    call find_CIA_consituents()

    ! Find the prf VMR_indexes of the CIA consituent species
    do s = 1, nCIA

      ! Check for 3 species special
      if (CIA_tab(s)%i3 .eqv. .True.) then
        ni = 3
      else
        ni = 2
      end if

      do i = 1, ni
        exists = .False.

        do j = 1, ngas
          if (ni == 2) then
            if (CIA_tab(s)%sp_con(i) == g_name(j)) then
              CIA_tab(s)%iVMR(i) = j
              exists = .True.
              exit
            end if
          else if (ni == 3) then
            if (CIA_tab(s)%sp_con_3(i) == g_name(j)) then
              CIA_tab(s)%iVMR_3(i) = j
              exists = .True.
              exit
            end if
          end if
        end do

        if (exists .eqv. .False.) then
          print*, 'ERROR - Specified CIA species component not found in prf VMR list - STOPPING'
          if (ni == 2) then
            print*, 'Species 2 part: ', CIA_tab(s)%sp, CIA_tab(s)%sp_con(i)
          else if (ni == 3) then
            print*, 'Species 3 part: ', CIA_tab(s)%sp, CIA_tab(s)%sp_con_3(i)
          end if
          stop
        end if

      end do
    end do

    ! Read the CIA tables
    call read_CIA_tables()

    print*, ' ~~ Performing CIA interpolation and output ~~ '
    print*, ' ~~ Please wait... ~~ '

    !! Begin openMP loops
    !$omp parallel default (none), &
    !$omp& private (l,z), &
    !$omp& shared (nwl,nlay,CIA_out,RH_lay,wl), &
    !$omp& firstprivate(CIA_work)

    ! Perform CIA table interpolation to
    ! Species loops are inside subroutines
    do l = 1, nwl
      !$omp single
      if (mod(l,nwl/10) == 0) then
        print*, l, wl(l), nwl
      end if
      !$omp end single
      !$omp do schedule (dynamic)
      do z = 1, nlay

        ! Find the CIA opacity for this layer from tables
        !call interp_CIA_tables(l,z,CIA_work)
        call interp_CIA_tables_Bezier(l,z,CIA_work)
        ! Convert interpolated result to cm2 g-1 of atmosphere and add to output array
        CIA_out(z) = CIA_work/RH_lay(z)

      end do
      !$omp end do

      !$omp single
      ! Output CMCRT formatted CIA table for layers
      call output_CIA_table(l)
      !$omp end single

    end do
    !$omp end parallel

    !deallocate all allocated arrays
    deallocate(CIA_out,CIA_write)
    deallocate(CIA_tab)
    deallocate(form,paths)
    ! Close uCIA i/o unit
    close(uCIA)

    print*, ' ~~ Quest completed  ~~ '

  end subroutine calc_CIA_table

  subroutine find_CIA_consituents()
    implicit none

    integer :: s

    do s = 1, nCIA

      select case(CIA_tab(s)%sp)

      case('H2-H2')
        CIA_tab(s)%sp_con(1) = 'H2'
        CIA_tab(s)%sp_con(2) = 'H2'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 113

      case('H2-He','He-H2')
        CIA_tab(s)%sp_con(1) = 'H2'
        CIA_tab(s)%sp_con(2) = 'He'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 334

      case('H-')
        CIA_tab(s)%sp_con_3(1) = 'H-'
        CIA_tab(s)%sp_con_3(2) = 'H'
        CIA_tab(s)%sp_con_3(3) = 'e-'

        CIA_tab(s)%i3 = .True.

      case('He-')
        CIA_tab(s)%sp_con(1) = 'He'
        CIA_tab(s)%sp_con(2) = 'e-'

      case('H2-')
        CIA_tab(s)%sp_con(1) = 'H2'
        CIA_tab(s)%sp_con(2) = 'e-'

      case('H2-H','H-H2')
        CIA_tab(s)%sp_con(1) = 'H2'
        CIA_tab(s)%sp_con(2) = 'H'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 4

      case('H-He','He-H')
        CIA_tab(s)%sp_con(1) = 'He'
        CIA_tab(s)%sp_con(2) = 'H'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 10

      case('CO2-CO2')
        CIA_tab(s)%sp_con(1) = 'CO2'
        CIA_tab(s)%sp_con(2) = 'CO2'

        CIA_tab(s)%nset = 3
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 10
        CIA_tab(s)%nT(2) = 6
        CIA_tab(s)%nT(3) = 3
        !CIA_tab(s)%nT(4) = 1

      case('N2-H2O','H2O-N2')
        CIA_tab(s)%sp_con(1) = 'N2'
        CIA_tab(s)%sp_con(2) = 'H2O'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 21

      case('N2-N2')
        CIA_tab(s)%sp_con(1) = 'N2'
        CIA_tab(s)%sp_con(2) = 'N2'

        CIA_tab(s)%nset = 6
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 14
        CIA_tab(s)%nT(2) = 10
        CIA_tab(s)%nT(3) = 10
        CIA_tab(s)%nT(4) = 5
        CIA_tab(s)%nT(5) = 5
        CIA_tab(s)%nT(6) = 14

      case('H2O-H2O')
        CIA_tab(s)%sp_con(1) = 'H2O'
        CIA_tab(s)%sp_con(2) = 'H2O'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 21

      case('CH4-CH4')
        CIA_tab(s)%sp_con(1) = 'CH4'
        CIA_tab(s)%sp_con(2) = 'CH4'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 7

      case('H2-CH4','CH4-H2')
        CIA_tab(s)%sp_con(1) = 'H2'
        CIA_tab(s)%sp_con(2) = 'CH4'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 10

      case('CH4-He','He-CH4')
        CIA_tab(s)%sp_con(1) = 'He'
        CIA_tab(s)%sp_con(2) = 'CH4'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 10

      case('O2-CO2','CO2-O2')
        CIA_tab(s)%sp_con(1) = 'O2'
        CIA_tab(s)%sp_con(2) = 'CO2'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 1

      case('O2-N2','N2-O2')
        CIA_tab(s)%sp_con(1) = 'O2'
        CIA_tab(s)%sp_con(2) = 'N2'

        CIA_tab(s)%nset = 5
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 7
        CIA_tab(s)%nT(2) = 5
        CIA_tab(s)%nT(3) = 5
        CIA_tab(s)%nT(4) = 1
        CIA_tab(s)%nT(5) = 1

      case('O2-O2')
        CIA_tab(s)%sp_con(1) = 'O2'
        CIA_tab(s)%sp_con(2) = 'O2'

        CIA_tab(s)%nset = 8
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 15
        CIA_tab(s)%nT(2) = 1
        CIA_tab(s)%nT(3) = 1
        CIA_tab(s)%nT(4) = 1
        CIA_tab(s)%nT(5) = 1
        CIA_tab(s)%nT(6) = 1
        CIA_tab(s)%nT(7) = 4
        CIA_tab(s)%nT(8) = 5

      case('H2O')
        CIA_tab(s)%sp_con(1) = 'H2O'
        CIA_tab(s)%sp_con(2) = 'H2O'

      case default
        print*, 'ERROR - CIA species constituents could not be found - STOPPING'
        print*, 'Species: ', CIA_tab(s)%sp
        stop
      end select

    end do

  end subroutine find_CIA_consituents

  subroutine output_CIA_table(l)
    implicit none

    integer, intent(in) :: l
    integer :: z

    if (first_call .eqv. .True.) then
      !print*, 'Outputing CIA.cmcrt'
      ! Output k-table in 1D or flattened 3D CMCRT format k_CMCRT.ktb (single precision)
      open(newunit=uCIA, file='CIA.cmcrt', action='readwrite', &
      & form='unformatted',status='replace',access='stream')
      write(uCIA) nlay, nwl
      first_call = .False.
    end if

    ! Convert to single precision on output, also care for underfloat
    CIA_write(:) = real(max(CIA_out(:),1.0e-30_dp),kind=sp)
    write(uCIA) CIA_write

  end subroutine output_CIA_table

end module CIA_tables_mod
