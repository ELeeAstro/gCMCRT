module xsec_tables_mod
  use optools_data_mod
  use xsec_tables_read
  use xsec_tables_interp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none

  logical :: first_call = .True.

  real(kind=dp), allocatable, dimension(:) :: xsec_out
  real(kind=sp), allocatable, dimension(:) :: xsec_write

  ! Namelist variables
  integer :: iopts
  integer, allocatable, dimension(:) :: form
  character(len=150), allocatable, dimension(:) :: paths

  namelist /xsec_nml/ iopts, form, paths

  private :: output_xsec_table, combine_xsec_opacity, validate_xsec_tables
  public :: calc_xsec_table

contains

  subroutine calc_xsec_table()
   implicit none
   
    integer :: s, j, l, z
    logical :: exists
    real(kind=dp), allocatable, dimension(:) :: xsec_work
    real(kind=dp) :: xsec_comb

    ! Allocate number of xsec tables
    allocate(xsec_tab(nxsec))

    ! Allocate required number of arrays from namelist options
    allocate(form(nxsec),paths(nxsec))

    ! Read xsec namelist parameters
    read(u_nml, nml=xsec_nml)

    ! Allocate work arrays
    allocate(xsec_out(nlay),xsec_write(nlay))

    ! Allocate private work arrays and initialise
    allocate(xsec_work(nxsec))
    xsec_work(:) = 0.0_dp
    xsec_comb = 0.0_dp

    ! Give the classes some global data from par and namelists
    xsec_tab(:)%sp = xsec_name(:)
    xsec_tab(:)%form = form(:)
    xsec_tab(:)%path = paths(:)

    ! Find the PRF VMR indices of the cross-section species.
    do s = 1, nxsec
      exists = .False.
      do j = 1, ngas
        if (xsec_tab(s)%sp == g_name(j)) then
          xsec_tab(s)%iVMR = j
          exists = .True.
          exit
        end if
      end do
      if (exists .eqv. .False.) then
        print*, 'ERROR - Specified xsec species not found in prf VMR list - STOPPING'
        print*, 'Species: ', xsec_tab(s)%sp
        stop
      end if
    end do  

    ! Read the xsec tables
    call read_xsec_tables()

    call validate_xsec_tables()

    print*, ' ~~ Performing xsec interpolation, combining and output ~~ '
    print*, ' ~~ Please wait... ~~ '

    !! Begin OpenMP loops.
    !$omp parallel default (none), &
    !$omp& private (l,z), &
    !$omp& shared (nwl,wl,nlay,xsec_out,RH_lay), &
    !$omp& firstprivate(xsec_work, xsec_comb)

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

        ! Find the xsec opacity for this layer from tables
        call interp_xsec_tables_Bezier(l,z,xsec_work(:))

        ! Combine interpolated xsec opacity for each species with VMR of species
        call combine_xsec_opacity(z,xsec_work(:),xsec_comb)

        ! Convert interpolated result to cm2 g-1 of atmosphere and add to output array
        xsec_out(z) = xsec_comb/RH_lay(z)

      end do
      !$omp end do
  
      !$omp single
      ! Output CMCRT formatted xsec table for layers
      call output_xsec_table(l)
      !$omp end single

    end do
    !$omp end parallel

    deallocate(xsec_tab)
    deallocate(form,paths)
    deallocate(xsec_out,xsec_write)
    deallocate(xsec_work)
    close(uxsec)

  end subroutine calc_xsec_table

  subroutine validate_xsec_tables()
    implicit none

    integer :: s, i

    do s = 1, nxsec
      if (.not. allocated(xsec_tab(s)%wl) .or. .not. allocated(xsec_tab(s)%lx_abs)) then
        print*, 'ERROR - Xsec table did not provide all required interpolation data - STOPPING'
        print*, 'Species, path: ', xsec_tab(s)%sp, trim(xsec_tab(s)%path)
        stop
      end if

      if (xsec_tab(s)%nwl < 3) then
        print*, 'ERROR - Xsec Bezier interpolation requires at least 3 wavelength points - STOPPING'
        print*, 'Species, nwl: ', xsec_tab(s)%sp, xsec_tab(s)%nwl
        stop
      end if

      if (any(.not. ieee_is_finite(xsec_tab(s)%wl)) .or. &
        & any(.not. ieee_is_finite(xsec_tab(s)%lx_abs))) then
        print*, 'ERROR - Xsec interpolation data must be finite - STOPPING'
        print*, 'Species, path: ', xsec_tab(s)%sp, trim(xsec_tab(s)%path)
        stop
      end if

      if (any(xsec_tab(s)%wl <= 0.0_dp)) then
        print*, 'ERROR - Xsec wavelength grid must be positive - STOPPING'
        print*, 'Species, path: ', xsec_tab(s)%sp, trim(xsec_tab(s)%path)
        stop
      end if

      do i = 2, xsec_tab(s)%nwl
        if (xsec_tab(s)%wl(i) <= xsec_tab(s)%wl(i-1)) then
          print*, 'ERROR - Xsec wavelength grid must be strictly increasing - STOPPING'
          print*, 'Species, index, values: ', xsec_tab(s)%sp, i, &
            & xsec_tab(s)%wl(i-1), xsec_tab(s)%wl(i)
          stop
        end if
      end do
    end do

  end subroutine validate_xsec_tables

  subroutine combine_xsec_opacity(z,xsec_work,xsec_comb)
    implicit none

    integer, intent(in) :: z
    real(kind=dp), dimension(nxsec), intent(in) :: xsec_work
    real(kind=dp), intent(out) :: xsec_comb

    integer :: s

    ! Combine all the xsec species weighted by the layer VMR
    ! Convert to units of [cm-1] by * layer number density

    xsec_comb = 0.0_dp
    do s = 1, nxsec
      xsec_comb = xsec_comb + xsec_work(s) * N_lay(z) * VMR_lay(xsec_tab(s)%iVMR,z)
    end do

  end subroutine combine_xsec_opacity

  subroutine output_xsec_table(l)
    implicit none

    integer, intent(in) :: l
    integer :: z, reclen

    if (first_call .eqv. .True.) then
      inquire(iolength=reclen) xsec_write
      !print*, 'Outputting xsec.cmcrt'
      ! Output table in 1D or flattened 3D CMCRT format k_CMCRT.ktb (single precision)
      open(newunit=uxsec, file='xsec.cmcrt', action='readwrite',&
        & form='unformatted',status='replace', access='direct',recl=reclen)
      first_call = .False.
    end if

    ! Convert to single precision on output and protect against underflow.
    xsec_write(:) = real(max(xsec_out(:),1.0e-30_dp),kind=sp)
    write(uxsec,rec=l) xsec_write

  end subroutine output_xsec_table

end module xsec_tables_mod
