module lbl_tables_mod
  use optools_data_mod
  use lbl_tables_read
  use lbl_tables_interp
  use lbl_tables_combine
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
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

  private :: output_lbl_table, validate_lbl_tables
  public :: calc_lbl_table


contains

  subroutine calc_lbl_table()
    implicit none

    integer :: s, j, l, z
    logical :: exists
    real(kind=dp), allocatable, dimension(:) :: lbl_work
    real(kind=dp) :: lbl_comb

    ! The array-allocation order is important because it follows the namelist order.

    ! Allocate number of lbl tables
    allocate(lbl_tab(nlbl))

    ! Allocate required number of arrays from namelist options
    allocate(form(nlbl),paths(nlbl))

    ! Read lbl namelist parameters
    read(u_nml, nml=lbl_nml)

    if (interp_wl .eqv. .True.) then
      print*, 'ERROR - LBL wavelength interpolation is not supported - STOPPING'
      print*, 'Set interp_wl = .False. and use the calculation wavelength grid in the LBL tables.'
      stop
    end if

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

    ! Find the PRF VMR indices of the LBL species.
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
        print*, 'ERROR - Specified lbl species not found in prf VMR list - STOPPING'
        print*, 'Species: ', lbl_tab(s)%sp
        stop
      end if
    end do

    ! Read the lbl tables
    call read_lbl_tables()

    call validate_lbl_tables()

    print*, ' ~~ Performing lbl interpolation, combining and output ~~ '
    print*, ' ~~ Please wait... ~~ '

    !! Begin OpenMP loops.
    !$omp parallel default (none), &
    !$omp& private (l,z), &
    !$omp& shared (nwl,wl,nlay,lbl_out,RH_lay,interp_wl), &
    !$omp& firstprivate(lbl_work, lbl_comb)

    ! Perform lbl table interpolation to layer T,p
    ! Species loops are inside subroutines
    do l = 1, nwl
      !$omp single
      if (mod(l,max(1,nwl/10)) == 0) then
        print*, l, wl(l), nwl
      end if
      !$omp end single
      !$omp do schedule (dynamic)
      do z = 1, nlay

        ! Find the lbl opacity for this layer from tables
        call interp_lbl_tables_Bezier(l,z,lbl_work(:))

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
    ! Close the LBL I/O unit.
    close(ulbl)

  end subroutine calc_lbl_table

  subroutine validate_lbl_tables()
    implicit none

    integer :: s, i
    real(kind=dp), parameter :: range_tol = 1.0e-2_dp
    real(kind=dp) :: scale

    do s = 1, nlbl
      if (.not. allocated(lbl_tab(s)%T) .or. .not. allocated(lbl_tab(s)%lT) .or. &
        & .not. allocated(lbl_tab(s)%P) .or. .not. allocated(lbl_tab(s)%lP) .or. &
        & .not. allocated(lbl_tab(s)%wl) .or. .not. allocated(lbl_tab(s)%lk_abs)) then
        print*, 'ERROR - LBL table did not provide all required interpolation data - STOPPING'
        print*, 'Species, path: ', lbl_tab(s)%sp, trim(lbl_tab(s)%path)
        stop
      end if

      if (lbl_tab(s)%nT < 3 .or. lbl_tab(s)%nP < 3) then
        print*, 'ERROR - LBL Bezier interpolation requires at least 3 T and P points - STOPPING'
        print*, 'Species, nT, nP: ', lbl_tab(s)%sp, lbl_tab(s)%nT, lbl_tab(s)%nP
        stop
      end if

      if (lbl_tab(s)%nwl /= nwl) then
        print*, 'ERROR - LBL table wavelength count does not match wavelengths.wl - STOPPING'
        print*, 'Species, table nwl, calculation nwl: ', lbl_tab(s)%sp, lbl_tab(s)%nwl, nwl
        stop
      end if

      if (any(.not. ieee_is_finite(lbl_tab(s)%T)) .or. &
        & any(.not. ieee_is_finite(lbl_tab(s)%P)) .or. &
        & any(.not. ieee_is_finite(lbl_tab(s)%wl)) .or. &
        & any(.not. ieee_is_finite(lbl_tab(s)%lk_abs))) then
        print*, 'ERROR - LBL interpolation data must be finite - STOPPING'
        print*, 'Species, path: ', lbl_tab(s)%sp, trim(lbl_tab(s)%path)
        stop
      end if

      if (any(lbl_tab(s)%T <= 0.0_dp) .or. any(lbl_tab(s)%P <= 0.0_dp) .or. &
        & any(lbl_tab(s)%wl <= 0.0_dp)) then
        print*, 'ERROR - LBL temperature, pressure, and wavelength grids must be positive - STOPPING'
        print*, 'Species, path: ', lbl_tab(s)%sp, trim(lbl_tab(s)%path)
        stop
      end if

      do i = 2, lbl_tab(s)%nT
        if (lbl_tab(s)%T(i) <= lbl_tab(s)%T(i-1)) then
          print*, 'ERROR - LBL temperature grid must be strictly increasing - STOPPING'
          print*, 'Species, index, values: ', lbl_tab(s)%sp, i, &
            & lbl_tab(s)%T(i-1), lbl_tab(s)%T(i)
          stop
        end if
      end do

      do i = 2, lbl_tab(s)%nP
        if (lbl_tab(s)%P(i) <= lbl_tab(s)%P(i-1)) then
          print*, 'ERROR - LBL pressure grid must be strictly increasing - STOPPING'
          print*, 'Species, index, values: ', lbl_tab(s)%sp, i, &
            & lbl_tab(s)%P(i-1), lbl_tab(s)%P(i)
          stop
        end if
      end do

      ! Bin centres can legitimately differ slightly between the LBL table
      ! and wavelengths.wl (e.g. different centring conventions), so only
      ! check the overall range rather than an exact per-point match.
      scale = max(1.0_dp, abs(lbl_tab(s)%wl(1)), abs(wl(1)))
      if (abs(lbl_tab(s)%wl(1) - wl(1)) > range_tol*scale) then
        print*, 'ERROR - LBL wavelength grid range does not match wavelengths.wl - STOPPING'
        print*, 'Species, table wl(1), calculation wl(1): ', lbl_tab(s)%sp, &
          & lbl_tab(s)%wl(1), wl(1)
        stop
      end if

      scale = max(1.0_dp, abs(lbl_tab(s)%wl(nwl)), abs(wl(nwl)))
      if (abs(lbl_tab(s)%wl(nwl) - wl(nwl)) > range_tol*scale) then
        print*, 'ERROR - LBL wavelength grid range does not match wavelengths.wl - STOPPING'
        print*, 'Species, table wl(nwl), calculation wl(nwl): ', lbl_tab(s)%sp, &
          & lbl_tab(s)%wl(nwl), wl(nwl)
        stop
      end if
    end do

  end subroutine validate_lbl_tables


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

    ! Convert to single precision on output and protect against underflow.
    lbl_write(:) = real(max(lbl_out(:),1.0e-30_dp),kind=sp)
    write(ulbl,rec=l) lbl_write

  end subroutine output_lbl_table

end module lbl_tables_mod
