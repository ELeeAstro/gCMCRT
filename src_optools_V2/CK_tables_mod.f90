module CK_tables_mod
  use optools_data_mod
  use CK_tables_read
  use CK_tables_interp
  use CK_table_RO
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  !use CK_table_rebin
  implicit none

  logical :: first_call = .True.

  real(kind=dp), allocatable, dimension(:,:) :: CK_out
  real(kind=sp), allocatable, dimension(:,:) :: CK_write

  real(kind=dp), allocatable, dimension(:) :: Gw, Gx

  ! Namelist variables
  integer :: iopts, nG, gdist
  real(kind=dp) :: gmin1, gmax1, gmin2, gmax2
  integer, allocatable, dimension(:) :: form
  character(len=150), allocatable, dimension(:) :: paths
  logical :: pre_mixed, rebin
  integer :: nrebin
  logical :: interp_wl

  namelist /CK_nml/ iopts, form, paths, nG, gmin1, gmax1, gmin2, gmax2, &
    & pre_mixed, rebin, nrebin, interp_wl

  private :: output_CK_table, output_CK_gord, validate_CK_tables
  public :: calc_CK_table


contains

  subroutine calc_CK_table()
    implicit none

    integer :: s, j, l, z
    logical :: exists
    real(kind=dp), allocatable, dimension(:) :: CK_RO
    real(kind=dp), allocatable, dimension(:,:) :: CK_work

    ! The array-allocation order is important because it follows the namelist order.

    ! Allocate number of CK tables
    allocate(CK_tab(nCK))

    ! Allocate required number of arrays from namelist options
    allocate(form(nCK),paths(nCK))

    ! Read CK namelist parameters
    read(u_nml, nml=CK_nml)

    if (rebin .eqv. .True.) then
      print*, 'ERROR - CK rebinning is not implemented - STOPPING'
      print*, 'Set rebin = .False. in CK_nml.'
      stop
    end if

    ! Allocate work arrays
    allocate(CK_out(nG,nlay),CK_write(nG,nlay))

    allocate(Gx(nG),Gw(nG))

    ! Allocate private work arrays and initialise
    allocate(CK_work(nCK,nG),CK_RO(nG))
    CK_work(:,:) = 0.0_dp
    CK_RO(:) = 0.0_dp

    ! Give the classes some global data from par and namelists
    CK_tab(:)%sp = CK_name(:)
    CK_tab(:)%form = form(:)
    CK_tab(:)%path = paths(:)

    ! Find the PRF VMR indices of the CK species.
    if (pre_mixed .eqv. .False.) then
      do s = 1, nCK
        exists = .False.
        do j = 1, ngas
          if (CK_tab(s)%sp == g_name(j)) then
            CK_tab(s)%iVMR = j
            exists = .True.
            exit
          end if
        end do
        if (exists .eqv. .False.) then
          print*, 'ERROR - Specified CK species not found in prf VMR list - STOPPING'
          print*, 'Species: ', CK_tab(s)%sp
          stop
        end if
      end do
    end if

    ! Read the CK tables
    call read_CK_tables(pre_mixed)

    ! Validate table dimensions and grids before interpolation or random overlap
    call validate_CK_tables()

    ! Rebin each k table if requested.
    !if (rebin .eqv. .True.) then
      !call rebin_CK_tables(nrebin,nG)
    !end if

    ! Use the first table's g grid and weights.
    Gx(:) = CK_tab(1)%Gx(:)
    Gw(:) = CK_tab(1)%Gw(:)

    if (pre_mixed .eqv. .False.) then
      print*, ' ~~ Performing CK interpolation, RO and output ~~ '
    else
      print*, ' ~~ Performing CK premixed interpolation and output ~~ '
    end if
    print*, ' ~~ Please wait... ~~ '

    !! Begin OpenMP loops.
    !$omp parallel default (none), &
    !$omp& private (l,z), &
    !$omp& shared (nwl,wl,nlay,CK_out,RH_lay,nG,Gw,Gx,pre_mixed,N_lay,interp_wl), &
    !$omp& firstprivate(CK_work,CK_RO)


    ! Interpolate the CK tables to the layer temperature and pressure.
    ! Species loops are contained within the interpolation subroutines.
    do l = 1, nwl
      !$omp single
      if (mod(l,max(1,nwl/10)) == 0) then
        print*, l, wl(l), nwl
      end if
      !$omp end single
      !$omp do schedule (dynamic)
      do z = 1, nlay

        ! Determine the CK opacity for this layer from the input tables.
        if (interp_wl .eqv. .True.) then
          call interp_CK_tables_wl(l,z,nG,CK_work(:,:))
        else
          !call interp_CK_tables(l,z,nG,CK_work(:,:))
          call interp_CK_tables_Bezier(l,z,nG,CK_work(:,:))
        end if

        if (pre_mixed .eqv. .True.) then
          ! Return the interpolated pre-mixed CK table to the output array.
          CK_out(:,z) = (CK_work(1,:)*N_lay(z))/RH_lay(z)
          !CK_out(:,z) = (CK_work(1,:))/RH_lay(z)
        else
          ! Perform random overlap with resorting and rebinning for all species.
          call RO_CK_RORR(z,nG,Gw(:),CK_work(:,:),CK_RO(:))
          !call RO_CK_2(z,nG,Gw(:),Gx(:),CK_work(:,:),CK_RO(:))
          !call RO_CK(z,nG,Gw(:),CK_work(:,:),CK_RO(:))
          ! Convert the overlapped result to cm^2 g^-1 of atmosphere.
          CK_out(:,z) = CK_RO(:)/RH_lay(z)
        end if


      end do
      !$omp end do

      !$omp single
      ! Write the CMCRT-formatted CK table for this wavelength bin.
      call output_CK_table(l)
      !$omp end single

    end do
    !$omp end parallel

    print*, ' ~~ Quest completed ~~ '

    ! Write the g ordinates and weights.
    call output_CK_gord()

    ! Deallocate all local arrays.
    deallocate(CK_out,CK_write)
    deallocate(CK_work)
    deallocate(CK_tab)
    deallocate(form,paths)
    ! Close the CK output unit.
    close(uCK)

  end subroutine calc_CK_table

  subroutine validate_CK_tables()
    implicit none

    integer :: s, i
    real(kind=dp), parameter :: grid_tol = 1.0e-10_dp
    real(kind=dp) :: scale

    if (nCK < 1) then
      print*, 'ERROR - Correlated-k opacity is enabled but no CK tables are configured - STOPPING'
      stop
    end if

    do s = 1, nCK
      if (.not. allocated(CK_tab(s)%Gx) .or. .not. allocated(CK_tab(s)%Gw) .or. &
        & .not. allocated(CK_tab(s)%wl) .or. .not. allocated(CK_tab(s)%T) .or. &
        & .not. allocated(CK_tab(s)%P)) then
        print*, 'ERROR - CK table did not provide all required grids - STOPPING'
        print*, 'Species, path: ', CK_tab(s)%sp, trim(CK_tab(s)%path)
        stop
      end if

      if (CK_tab(s)%nG /= nG) then
        print*, 'ERROR - CK table nG does not match the CK namelist - STOPPING'
        print*, 'Species, file nG, namelist nG: ', CK_tab(s)%sp, CK_tab(s)%nG, nG
        stop
      end if

      if (CK_tab(s)%nwl /= nwl) then
        print*, 'ERROR - CK table wavelength count does not match wavelengths.wl - STOPPING'
        print*, 'Species, table nwl, calculation nwl: ', CK_tab(s)%sp, CK_tab(s)%nwl, nwl
        stop
      end if

      if (CK_tab(s)%nT < 3 .or. CK_tab(s)%nP < 3) then
        print*, 'ERROR - CK Bezier interpolation requires at least 3 T and P points - STOPPING'
        print*, 'Species, nT, nP: ', CK_tab(s)%sp, CK_tab(s)%nT, CK_tab(s)%nP
        stop
      end if

      if (any(CK_tab(s)%T <= 0.0_dp) .or. any(CK_tab(s)%P <= 0.0_dp)) then
        print*, 'ERROR - CK temperature and pressure grids must be positive - STOPPING'
        print*, 'Species, path: ', CK_tab(s)%sp, trim(CK_tab(s)%path)
        stop
      end if

      if (any(.not. ieee_is_finite(CK_tab(s)%T)) .or. &
        & any(.not. ieee_is_finite(CK_tab(s)%P)) .or. &
        & any(.not. ieee_is_finite(CK_tab(s)%lk_abs))) then
        print*, 'ERROR - CK interpolation data must be finite - STOPPING'
        print*, 'Species, path: ', CK_tab(s)%sp, trim(CK_tab(s)%path)
        stop
      end if

      do i = 2, CK_tab(s)%nT
        if (CK_tab(s)%T(i) <= CK_tab(s)%T(i-1)) then
          print*, 'ERROR - CK temperature grid must be strictly increasing - STOPPING'
          print*, 'Species, index, values: ', CK_tab(s)%sp, i, &
            & CK_tab(s)%T(i-1), CK_tab(s)%T(i)
          stop
        end if
      end do

      do i = 2, CK_tab(s)%nP
        if (CK_tab(s)%P(i) <= CK_tab(s)%P(i-1)) then
          print*, 'ERROR - CK pressure grid must be strictly increasing - STOPPING'
          print*, 'Species, index, values: ', CK_tab(s)%sp, i, &
            & CK_tab(s)%P(i-1), CK_tab(s)%P(i)
          stop
        end if
      end do

      do i = 1, nwl
        scale = max(1.0_dp, abs(CK_tab(s)%wl(i)), abs(wl(i)))
        if (abs(CK_tab(s)%wl(i) - wl(i)) > grid_tol*scale) then
          print*, 'ERROR - CK wavelength grid does not match wavelengths.wl - STOPPING'
          print*, 'Species, index, table wl, calculation wl: ', CK_tab(s)%sp, i, &
            & CK_tab(s)%wl(i), wl(i)
          stop
        end if
      end do
    end do

    do s = 2, nCK
      do i = 1, nG
        scale = max(1.0_dp, abs(CK_tab(s)%Gx(i)), abs(CK_tab(1)%Gx(i)))
        if (abs(CK_tab(s)%Gx(i) - CK_tab(1)%Gx(i)) > grid_tol*scale) then
          print*, 'ERROR - CK Gx grids differ between species - STOPPING'
          print*, 'Species, index, value, reference: ', CK_tab(s)%sp, i, &
            & CK_tab(s)%Gx(i), CK_tab(1)%Gx(i)
          stop
        end if

        scale = max(1.0_dp, abs(CK_tab(s)%Gw(i)), abs(CK_tab(1)%Gw(i)))
        if (abs(CK_tab(s)%Gw(i) - CK_tab(1)%Gw(i)) > grid_tol*scale) then
          print*, 'ERROR - CK Gw grids differ between species - STOPPING'
          print*, 'Species, index, value, reference: ', CK_tab(s)%sp, i, &
            & CK_tab(s)%Gw(i), CK_tab(1)%Gw(i)
          stop
        end if
      end do
    end do

  end subroutine validate_CK_tables

  subroutine output_CK_table(l)
    implicit none

    integer, intent(in) :: l
    integer :: z, g, reclen

    if (first_call .eqv. .True.) then
      ! Output k-table in 1D flattened 3D CMCRT format CK.cmcrt (single precision)
      inquire(iolength=reclen) CK_write
      open(newunit=uCK, file='CK.cmcrt', action='readwrite', &
              & form='unformatted', status='replace', access='direct',recl=reclen)
      first_call = .False.
    end if

    ! Convert to single precision on output and protect against underflow.
    CK_write = real(max(CK_out,1.0e-30_dp),kind=sp)
    write(uCK,rec=l) CK_write

  end subroutine output_CK_table

  subroutine output_CK_gord()
    implicit none
    integer :: g, u_g
    real(kind=dp) :: sum1

    print*, ' ~~ Outputting gord.cmcrt ~~ '

    ! Write the g ordinates and their weights to gord.cmcrt.
    open(newunit=u_g, file='gord.cmcrt', action='readwrite',form='formatted')
    write(u_g,*) nG

    sum1 = 0.0_dp
    do g = 1, nG
      write(u_g,*) Gx(g), Gw(g)
      sum1 = sum1 + Gw(g)
      !print*, g,  Gx(g), Gw(g), sum1
    end do

    print*,'G-ordinance sums:', sum(Gx(:)), sum(Gw(:))

    close(u_g)

    print*, ' ~~ Quest completed  ~~'

end subroutine output_CK_gord


end module CK_tables_mod
