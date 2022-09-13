module CK_tables_mod
  use optools_data_mod
  use CK_tables_read
  use CK_tables_interp
  use CK_table_RO
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

  private :: output_CK_table, output_CK_gord
  public :: calc_CK_table


contains

  subroutine calc_CK_table()
    implicit none

    integer :: s, j, l, z
    logical :: exists
    real(kind=dp), allocatable, dimension(:) :: CK_RO
    real(kind=dp), allocatable, dimension(:,:) :: CK_work

    ! Note, the order the array allocations is important (due to namelist order)

    ! Allocate number of CK tables
    allocate(CK_tab(nCK))

    ! Allocate required number of arrays from namelist options
    allocate(form(nCK),paths(nCK))

    ! Read CK namelist parameters
    read(u_nml, nml=CK_nml)

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

    ! Find the prf VMR_indexes of the CK species
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
          print*, 'ERROR - Specifed CK species not found in prf VMR list - STOPPING'
          print*, 'Species: ', CK_tab(s)%sp
          stop
        end if
      end do
    end if

    ! Read the CK tables
    call read_CK_tables(pre_mixed)

    ! Rebin each k-table if wanted
    !if (rebin .eqv. .True.) then
      !call rebin_CK_tables(nrebin,nG)
    !end if

    ! make weights equal to 1st table for now.
    Gx(:) = CK_tab(1)%Gx(:)
    Gw(:) = CK_tab(1)%Gw(:)

    if (pre_mixed .eqv. .False.) then
      print*, ' ~~ Performing CK interpolation, RO and output ~~ '
    else
      print*, ' ~~ Performing CK premixed interpolation and output ~~ '
    end if
    print*, ' ~~ Please wait... ~~ '

    !! Begin openMP loops
    !$omp parallel default (none), &
    !$omp& private (l,z), &
    !$omp& shared (nwl,wl,nlay,CK_out,RH_lay,nG,Gw,Gx,pre_mixed,N_lay,interp_wl), &
    !$omp& firstprivate(CK_work,CK_RO)


    ! Perform CK table interpolation to layer T,p
    ! Species loops are inside subroutines
    do l = 1, nwl
      !$omp single
      if (mod(l,nwl/10) == 0) then
        print*, l, wl(l), nwl
      end if
      !$omp end single
      !$omp do schedule (dynamic)
      do z = 1, nlay

        ! Find the CK opacity for this layer from tables
        if (interp_wl .eqv. .True.) then
          call interp_CK_tables_wl(l,z,nG,CK_work(:,:))
        else
          !call interp_CK_tables(l,z,nG,CK_work(:,:))
          call interp_CK_tables_Bezier(l,z,nG,CK_work(:,:))
        end if

        if (pre_mixed .eqv. .True.) then
          ! Pre-mixed, give back interpolated ck table to output array
          CK_out(:,z) = (CK_work(1,:)*N_lay(z))/RH_lay(z)
          !CK_out(:,z) = (CK_work(1,:))/RH_lay(z)
        else
          ! Perform the random overlap for all species
          call RO_CK_2(z,nG,Gw(:),Gx(:),CK_work(:,:),CK_RO(:))
          !call RO_CK(z,nG,Gw(:),CK_work(:,:),CK_RO(:))
          ! Convert overlapped result to cm2 g-1 of atmosphere and add to output array
          CK_out(:,z) = CK_RO(:)/RH_lay(z)
        end if


      end do
      !$omp end do

      !$omp single
      ! Output CMCRT formatted CK table for layers
      call output_CK_table(l)
      !$omp end single

    end do
    !$omp end parallel

    print*, ' ~~ Quest completed ~~ '

    ! Output g ordinances and weights
    call output_CK_gord()

    !deallocate all allocated arrays
    deallocate(CK_out,CK_write)
    deallocate(CK_work)
    deallocate(CK_tab)
    deallocate(form,paths)
    ! Close uCK i/o unit
    close(uCK)

  end subroutine calc_CK_table

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

    ! Convert to single precision on output, also care for underfloat
    CK_write = real(max(CK_out,1.0e-30_dp),kind=sp)
    write(uCK,rec=l) CK_write

  end subroutine output_CK_table

  subroutine output_CK_gord()
    implicit none
    integer :: g, u_g
    real(kind=dp) :: sum1

    print*, ' ~~ Outputing gord.cmcrt ~~ '

    ! Output g-ordinances and delg with weights to file g.ord
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
