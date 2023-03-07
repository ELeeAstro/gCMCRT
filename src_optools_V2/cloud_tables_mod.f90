module cloud_tables_mod
  use optools_data_mod
  use cloud_tables_read
  use cloud_tables_interp
  use cloud_tables_emt
  use cloud_tables_dist
  use cloud_tables_carma
  implicit none

  logical :: first_call = .True.

  real(kind=dp), allocatable, dimension(:) :: cl_out_k, cl_out_a, cl_out_g
  real(kind=sp), allocatable, dimension(:) :: cl_write

  ! Namelist options
  integer :: iopts
  integer, allocatable, dimension(:) :: form
  character(len=150), allocatable, dimension(:) :: paths

  character(len=20), dimension(0:8) :: cdist
  character(len=20), dimension(6) :: mie_meth

  logical :: cld_tab_read = .False.

  namelist /cl_nml/ iopts, form, paths, imix, sig, idist, imie, &
    & amin, amax, ndist, idist_int, eff_fac, veff, fmax, cld_tab_read

  private :: output_cl_table
  public :: calc_cloud_table


contains

  subroutine calc_cloud_table()
    implicit none

    integer :: s, j, l, z, a, m
    integer :: uext, ua, ug, n
    logical :: exists
    real(kind=dp), allocatable, dimension(:) :: n_work, k_work
    complex(kind=dp) :: eps_comb

    ! Note, the order of the array allocations is important (due to namelist order)

    ! Allocate number of cl species nk tables
    allocate(cl_tab(ncl))

    ! Allocate required number of arrays from namelist options
    allocate(form(ncl),paths(ncl))

    ! Read lbl namelist parameters
    read(u_nml, nml=cl_nml)

    ! Allocate work arrays
    allocate(cl_out_k(nlay), cl_out_a(nlay), cl_out_g(nlay), cl_write(nlay))

    if (cld_tab_read .eqv. .True.) then
      open(newunit=uext,file='cld_ext.txt',action='read')
      open(newunit=ua,file='cld_a.txt',action='read') 
      open(newunit=ug,file='cld_g.txt',action='read')

      do l = 1, nwl
       read(uext,*) n, cl_out_k(:)
       cl_out_k(:) = cl_out_k(:)!/RH_lay(:)
       read(ua,*) n, cl_out_a(:)
       read(ug,*) n, cl_out_g(:)
       call output_cl_table(l)
      end do
      return
    end if

    ! Allocate private work arrays and initialise
    allocate(n_work(ncl),k_work(ncl))
    n_work(:) = 0.0_dp
    k_work(:) = 0.0_dp
    eps_comb = (0.0_dp,0.0_dp)

    ! Give the classes some global data from par and namelists
    cl_tab(:)%sp = cl_name(:)
    cl_tab(:)%form = form(:)
    cl_tab(:)%path = paths(:)

    ! Find the clprf VMR_indexes of the lbl species
    do s = 1, ncl
      exists = .False.
      do j = 1, ndust
        if (cl_tab(s)%sp == d_name(j)) then
          cl_tab(s)%iVMR = j
          exists = .True.
          exit
        end if
      end do
      if (exists .eqv. .False.) then
        print*, 'ERROR - Specifed cloud species not found in clprf VMR list - STOPPING'
        print*, 'Species: ', cl_tab(s)%sp
        stop
      end if
    end do

    ! Some preparation for distribution calculations
    if (idist > 2) then
      allocate(a_dist(ndist))
      !! sample ndist grain sizes between amin and amax in log 10 space
      amin = amin * 1e-4_dp
      amax = amax * 1e-4_dp
      lamin = log10(amin)
      lamax = log10(amax)
      do a = 1, ndist
        a_dist(a) = 10.0_dp**((lamax-lamin) * real(a-1,kind=dp) / real(ndist-1,kind=dp) + lamin)
      end do

      ! ln the sigma (or not, might be easier not to)
      lsig = log(sig)
 
    end if

    ! Read the lbl tables
    call read_cl_nk()

    print*, ' ~~ Performing nk interpolation, EMT, Mie calculation and output ~~ '
    print*, ' ~~ Please wait... ~~ '

    cdist(0) = 'full' ; cdist(1) = 'delta peak' ; cdist(2) = 'tri- well peaked'
    cdist(3) = 'log-normal' ; cdist(4) = 'gamma' ; cdist(5) = 'inverse-gamma'
    cdist(6) = 'Rayleigh' ; cdist(7) = 'Hansen' ; cdist(8) = 'exponential'
    print*, ' -- Assuming a '//trim(cdist(idist))//' distribution: ', idist

    mie_meth(1) = 'MieX' ; mie_meth(2) = 'MieExt' ;  mie_meth(3) = 'BHMIE'
    mie_meth(4) = 'DHS'; mie_meth(5) = 'BHCOAT' ; mie_meth(6) = 'LX-MIE'
    print*, ' -- Using the '//trim(mie_meth(imie))//' method: ', imie

    !! Begin openMP loops
    !$omp parallel default (none), &
    !$omp& private (l,z), &
    !$omp& shared (nwl,nlay, wl, &
    !$omp& n_work,k_work,cl_out_k,cl_out_a,cl_out_g,RH_lay,nd_cl_lay,idist), &
    !$omp& firstprivate(eps_comb)

    ! Perform cl table interpolation to wavalength
    ! Species loops are inside subroutines
    do l = 1, nwl
      !$omp single
      if (mod(l,nwl/10) == 0) then
        print*, l, wl(l), nwl
      end if
      !$omp end single
      ! Find the nk constants for this wavelength
      !$omp single
      call interp_cl_tables(l,n_work(:),k_work(:))
      !$omp end single

      !$omp do schedule (dynamic)
      do z = 1, nlay

        if (idist > 0) then

          if (nd_cl_lay(z) < 1e-10_dp) then
            cl_out_k(z) = 1e-99_dp
            cl_out_a(z) = 0.0_dp
            cl_out_g(z) = 0.0_dp
            cycle
          end if

          ! Combine the n and k values for each species using EMT theory
          call emt_cl(z,n_work(:),k_work(:),eps_comb)
          ! perform Mie theory calculation
          call dist_cl(z, l, eps_comb, cl_out_k(z), cl_out_a(z), cl_out_g(z))
        else
          call calc_carma(z,l,n_work(:),k_work(:), cl_out_k(z), cl_out_a(z), cl_out_g(z))
        end if

        ! Convert result to cm2 g-1 of atmosphere and add to output array
        cl_out_k(z) = cl_out_k(z)/RH_lay(z)

      end do
      !$omp end do
 
      !$omp single
      ! Output CMCRT formatted cl table for layers
      call output_cl_table(l)
      !$omp end single

    end do
    !$omp end parallel

    print*, ' ~~ Quest completed ~~ '

    !deallocate all allocated arrays
    deallocate(cl_tab)
    deallocate(form,paths)
    deallocate(cl_out_k,cl_out_a,cl_out_g,cl_write)
    deallocate(n_work,k_work)
    ! Close ucl i/o units
    close(ucl_k) ; close(ucl_a) ; close(ucl_g)

  end subroutine calc_cloud_table

  subroutine output_cl_table(l)
    implicit none

    integer, intent(in) :: l
    integer :: z, reclen

    if (first_call .eqv. .True.) then
      inquire(iolength=reclen) cl_write
      ! Output lbl-table in 1D flattened 3D CMCRT format lbl.cmcrt (single precision)
      open(newunit=ucl_k, file='cl_k.cmcrt', action='readwrite', &
      & form='unformatted', status='replace', access='direct',recl=reclen)
      open(newunit=ucl_a, file='cl_a.cmcrt', action='readwrite', &
      & form='unformatted', status='replace', access='direct',recl=reclen)
      open(newunit=ucl_g, file='cl_g.cmcrt', action='readwrite', &
      & form='unformatted', status='replace', access='direct',recl=reclen)
      first_call = .False.
    end if

    ! Write the cloud/haze opacity, single scattering albedo and g to files
    ! Convert to single precision on output, also care for underfloat
    cl_write(:) = real(max(cl_out_k(:),1.0e-30_dp),kind=sp)
    write(ucl_k,rec=l) cl_write

    cl_write(:) = real(cl_out_a(:),kind=sp)
    write(ucl_a,rec=l) cl_write

    cl_write(:) = real(cl_out_g(:),kind=sp)
    write(ucl_g,rec=l) cl_write

  end subroutine output_cl_table

end module cloud_tables_mod
