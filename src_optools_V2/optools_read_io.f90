module read_io
  use optools_data_mod
  implicit none

  private
  public :: read_optools_par, read_wavelengths, read_prf, read_clprf

contains

  ! ------ Read optools.par file --------
  subroutine read_optools_par()
    implicit none

    integer :: u
    integer :: k, r, c

    print*, ' ~~ Reading in optools.par file ~~ '

    open(newunit=u, file='optools.par', form='formatted', status='old', action='read')

    read(u,*)

    read(u,*); read(u,*); read(u,*)
    read(u,'(A)') exp_name ; exp_name = trim(exp_name)

    read(u,*); read(u,*); read(u,*)
    read(u,*) corr_k
    read(u,*) lbl
    read(u,*) conti
    read(u,*) Ray_scat
    read(u,*) cloud_opc
    read(u,*) xsec_opc

    read(u,*); read(u,*); read(u,*)
    read(u,*) nCK

    allocate(CK_name(nCK))
    do k = 1, nCK
      read(u,'(A)') CK_name(k)
    end do

    read(u,*); read(u,*); read(u,*)
    read(u,*) nLbl

    allocate(Lbl_name(nLbl))
    do k = 1, nLbl
      read(u,'(A)') Lbl_name(k)
    end do

    read(u,*); read(u,*); read(u,*)
    read(u,*) nRay
    allocate(Ray_name(nRay))
    do r = 1, nRay
      read(u,*) Ray_name(r)
    end do

    read(u,*); read(u,*); read(u,*)
    read(u,*) nCIA
    allocate(CIA_name(nCIA))
    do c = 1, nCIA
      read(u,*) CIA_name(c)
    end do

    read(u,*); read(u,*); read(u,*)
    read(u,*) ncl
    allocate(cl_name(ncl))
    do c = 1, ncl
      read(u,*) cl_name(c)
    end do

    read(u,*); read(u,*); read(u,*)
    read(u,*) nxsec
    allocate(xsec_name(ncl))
    do c = 1, nxsec
      read(u,*) xsec_name(c)
    end do

    close(u)

    print*, '-----'

    print*, 'exp_name: ', exp_name

    print*, '-----'

    print*, corr_k, 'corr_k'
    print*, lbl, 'lbl'
    print*, conti, 'conti'
    print*, Ray_scat, 'Ray_scat'
    print*, cloud_opc, 'cloud_opc'
    print*, xsec_opc, 'xsec_opac'

    print*, '-----'

    print*, 'CK species: '
    do k = 1, nCK
      print*, CK_name(k)
    end do

    print*, '-----'

    print*, 'Lbl species: '
    do k = 1, nLbl
      print*, Lbl_name(k)
    end do

    print*, '-----'

    print*, 'Rayleigh scattering species: '
    do r = 1, nRay
      print*, Ray_name(r)
    end do

    print*, '-----'

    print*, 'CIA species: '
    do c = 1, nCIA
      print*, CIA_name(c)
    end do

    print*, '-----'

    print*, 'Cloud species: '
    do c = 1, ncl
      print*, cl_name(c)
    end do

    print*, 'Xsec species: '
    do c = 1, nxsec
      print*, xsec_name(c)
    end do

    print*, '-----'

    print*, ' ~~ Quest completed ~~ '

    print*, ' ~~ Opening optools.nml Namelist ~~ '
    open(newunit=u_nml, file='optools.nml', action='read')
    print*, ' ~~ Quest completed ~~ '

  end subroutine read_optools_par
  !-------------------- ENDS --------------------------------


  ! ------ Read optools wavelengths.wl file --------
  subroutine read_wavelengths()
    implicit none

    integer :: u, l


    print*, ' ~~ Reading in wavelengths.wl file ~~ '

    open(newunit=u, file='wavelengths.wl', form='formatted', status='old', action='read')

    read(u,*) nwl

    allocate(iwl(nwl),wl(nwl),wl_cm(nwl),wl_A(nwl),wn(nwl),freq(nwl))

    do l = 1, nwl
      ! Read in wavelengths [um]
      read(u,*) iwl(l), wl(l)
      ! Convert to cm
      wl_cm(l) = wl(l) * 1.0e-4_dp
      ! Convert to A
      wl_A(l) = wl(l) * 1.0e4_dp
      ! Convert to wavenumber [cm-1]
      wn(l) = 1.0_dp / wl_cm(l)
      ! Convert to frequency [Hz]
      freq(l) = c_s * wn(l)
    end do

    close(u)

    print*, ' ~~ Quest completed ~~ '

  end subroutine read_wavelengths
  !-------------------- ENDS --------------------------------


  ! ------ Read 1D profile *.prf file --------
  subroutine read_prf()
    implicit none

    integer :: u, i, j
    real(kind=dp), allocatable, dimension(:) :: VMR_dum

    print*, ' ~~ Reading in '//trim(exp_name)//'.prf file ~~ '

    open(newunit=u, file=trim(exp_name)//'.prf', form='formatted', status='old', action='read')

    read(u,*); read(u,*)
    read(u,*) nlay

    allocate(ilay(nlay), PG_lay(nlay), mu_lay(nlay), TG_lay(nlay), &
    & RH_lay(nlay), N_lay(nlay))

    read(u,*)
    read(u,*) ngas

    allocate(VMR_lay(ngas,nlay), g_name(ngas),  VMR_dum(ngas))

    do i = 1, ngas
      read(u,*) g_name(i)
      ! Switch gas name if electron named differently
      if (g_name(i) == 'el') then
        g_name(i) = 'e-'
      end if
      print*, g_name(i)
    end do

    read(u,*); read(u,*)
    do i = 1, nlay
      read(u,*) ilay(i), PG_lay(i), TG_lay(i), mu_lay(i), (VMR_dum(j), j=1,ngas)
      VMR_lay(:,i) = VMR_dum(:)
      PG_lay(i) = PG_lay(i) * bar ! Convert to dyne
      RH_lay(i) = (PG_lay(i) * mu_lay(i)) / (Rgas * TG_lay(i)) ! Find density [g cm-3]
      N_lay(i) = PG_lay(i) / (kb * TG_lay(i)) ! Find number density [cm-3]
    end do

    close(u)

    print*, ' ~~ Quest completed ~~ '

    deallocate(VMR_dum)

  end subroutine read_prf

!-------------------- ENDS --------------------------------

  subroutine read_clprf()
    implicit none

    integer :: u, i, j, k, l, sdum
    real(kind=dp) :: a_dum, var_dum, nd_dum
    real(kind=dp), allocatable, dimension(:) :: VMR_dum, a_dum_C, nd_dum_C


    print*, ' ~~ Reading in '//trim(exp_name)//'.clprf file ~~ '

    open(newunit=u, file=trim(exp_name)//'.clprf', form='formatted', status='old', action='read')

    read(u,*); read(u,*)
    read(u,*) nlay, nmode

    read(u,*)
    read(u,*) ndust

    allocate(ilay2(nlay), nd_cl_lay(nlay), a_cl_lay(nlay),var_cl_lay(nlay))
    allocate(VMR_cl_lay(ndust,nlay), d_name(ndust))
    allocate(VMR_dum(ndust))

    do i = 1, ndust
      read(u,*) d_name(i)
      print*, d_name(i)
    end do

    read(u,*); read(u,*)
    if (nmode == 1) then
      do i = 1, nlay
        read(u,*) ilay2(i), a_dum, var_dum, nd_dum , (VMR_dum(l), l=1,ndust)
        VMR_cl_lay(:,i) = VMR_dum(:)
        a_cl_lay(i) = a_dum * 1e-4_dp
        var_cl_lay(i) = var_dum
        nd_cl_lay(i) = nd_dum
      end do
    else
      allocate(a_dum_C(nmode),nd_dum_C(nmode))
      allocate(nd_C_cl_lay(ndust,nmode,nlay))
      allocate(a_C_cl_lay(ndust,nmode))
      do i = 1, ndust
         read(u,*) (a_dum_C(j), j=1,nmode)
         a_C_cl_lay(i,:) = a_dum_C(:) * 1e-4_dp
         !print*, a_C_cl_lay(i,:)
      end do
      do i = 1, nlay
        do j = 1, ndust
          read(u,*) ilay2(i), sdum,  (nd_dum_C(k), k=1,nmode)
          nd_C_cl_lay(j,:,i) = nd_dum_C(:)
          !print*, nd_C_cl_lay(j,:,i)
        end do
      end do
    end if


    close(u)

    print*, ' ~~ Quest completed ~~ '

    deallocate(VMR_dum)

  end subroutine read_clprf


end module read_io
