module lbl_tables_mod
  use optools_data_mod
  use lbl_tables_read
  use lbl_tables_interp
  use lbl_tables_combine
  implicit none

  logical :: first_call = .True.

  real(kind=dp), allocatable, dimension(:) :: lbl_out
  real(kind=sp), allocatable, dimension(:) :: lbl_write

  ! Namelist options
  integer :: iopts
  integer, allocatable, dimension(:) :: form
  character(len=150), allocatable, dimension(:) :: paths
  logical :: interp_wl

  namelist /lbl_nml/ iopts, form, paths, interp_wl

  private :: output_lbl_table
  public :: calc_lbl_table


contains

  subroutine calc_lbl_table()
    implicit none

    integer :: s, j, l, z
    logical :: exists
    real(kind=dp), allocatable, dimension(:) :: lbl_work
    real(kind=dp) :: lbl_comb

    ! Note, the order of the array allocations is important (due to namelist order)

    ! Allocate number of lbl tables
    allocate(lbl_tab(nlbl))

    ! Allocate required number of arrays from namelist options
    allocate(form(nlbl),paths(nlbl))

    ! Read lbl namelist parameters
    read(u_nml, nml=lbl_nml)

    ! Allocate work arrays
    allocate(lbl_out(nlay),lbl_write(nlay))

    ! Allocate private work arrays and initialise
    allocate(lbl_work(nlbl))
    lbl_work(:) = 0.0_dp
    lbl_comb = 0.0_dp

    ! Give the classes some global data from par and namelists
    lbl_tab(:)%sp = lbl_name(:)
    lbl_tab(:)%form = form(:)
    lbl_tab(:)%path = paths(:)

    ! Find the prf VMR_indexes of the lbl species
    do s = 1, nlbl
      exists = .False.
      do j = 1, ngas
        if (lbl_tab(s)%sp == g_name(j)) then
          lbl_tab(s)%iVMR = j
          exists = .True.
          exit
        end if
      end do
      if (exists .eqv. .False.) then
        print*, 'ERROR - Specifed lbl species not found in prf VMR list - STOPPING'
        print*, 'Species: ', lbl_tab(s)%sp
        stop
      end if
    end do

    ! Read the lbl tables
    call read_lbl_tables()

    print*, ' ~~ Performing lbl interpolation, combining and output ~~ '
    print*, ' ~~ Please wait... ~~ '

    !! Begin openMP loops
    !$omp parallel default (none), &
    !$omp& private (l,z), &
    !$omp& shared (nwl,wl,nlay,lbl_out,RH_lay,interp_wl), &
    !$omp& firstprivate(lbl_work, lbl_comb)

    ! Perform lbl table interpolation to layer T,p
    ! Species loops are inside subroutines
    do l = 1, nwl
      !$omp single
      if (mod(l,nwl/10) == 0) then
        print*, l, wl(l), nwl
      end if
      !$omp end single
      !$omp do schedule (dynamic)
      do z = 1, nlay

        ! Find the lbl opacity for this layer from tables
        if (interp_wl .eqv. .True.) then
           call interp_lbl_tables(l,z,lbl_work(:))
        else
          call interp_lbl_tables_Bezier(l,z,lbl_work(:))
        end if

        ! Combine interpolated lbl opacity for each species with VMR of species
        call combine_lbl_opacity(z,lbl_work(:),lbl_comb)

        ! Convert interpolated result to cm2 g-1 of atmosphere and add to output array
        lbl_out(z) = lbl_comb/RH_lay(z)

      end do
      !$omp end do

      !$omp single
      ! Output CMCRT formatted lbl table for layers
      call output_lbl_table(l)
      !$omp end single

    end do
    !$omp end parallel

    print*, ' ~~ Quest completed ~~ '

    !deallocate all allocated arrays
    deallocate(lbl_tab)
    deallocate(form,paths)
    deallocate(lbl_out,lbl_write)
    deallocate(lbl_work)
    ! Close ulbl i/o unit
    close(ulbl)
    close(ulbl_for)

  end subroutine calc_lbl_table


  subroutine output_lbl_table(l)
    implicit none

    integer, intent(in) :: l
    integer :: z, reclen

    if (first_call .eqv. .True.) then
      inquire(iolength=reclen) lbl_write
      ! Output lbl-table in 1D flattened 3D CMCRT format lbl.cmcrt (single precision)
      open(newunit=ulbl, file='lbl.cmcrt', action='readwrite', &
      & form='unformatted',status='replace',access='direct',recl=reclen)
      first_call = .False.
    end if

    ! Convert to single precision on output, also care for underfloat
    lbl_write(:) = real(max(lbl_out(:),1.0e-30_dp),kind=sp)
    write(ulbl,rec=l) lbl_write

  end subroutine output_lbl_table

end module lbl_tables_mod
