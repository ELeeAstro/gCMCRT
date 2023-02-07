module Ray_tables_mod
  use optools_data_mod
  implicit none

  ! Parameters for H2O
  real(kind=dp), dimension(8), parameter :: a_i = (/ 0.244257733_dp, 0.974634476e-2_dp, &
  & -0.373234996e-2_dp, 0.268678472e-3_dp, 0.158920570e-2_dp, 0.245934259e-2_dp, &
  & 0.900704920_dp, -0.166626219e-1_dp /)
  real(kind=dp), parameter :: rho_s = 1000.0_dp
  real(kind=dp), parameter :: lam_s = 0.589_dp
  real(kind=dp), parameter :: T_s = 273.15_dp
  real(kind=dp), parameter :: lam_ir = 5.432937_dp, lam_uv = 0.229202_dp

  ! Parametr for H
  real(kind=dp), parameter :: wl_ly = 121.567_dp * 1.0e-7_dp ! Lyman alpha wavelength [cm]
  real(kind=dp), parameter :: f_ly = c_s/wl_ly
  real(kind=dp), parameter :: w_l = (2.0_dp * pi * f_ly) / 0.75_dp
  real(kind=dp), dimension(10), parameter :: cp = (/1.26537_dp,3.73766_dp,8.8127_dp,19.1515_dp, &
  &  39.919_dp,81.1018_dp,161.896_dp,319.001_dp,622.229_dp,1203.82_dp/)

  ! Global work arrays
  real(kind=dp), allocatable, dimension(:) :: Ray_work
  real(kind=dp), allocatable, dimension(:) :: n_ref, King, nd_stp
  integer, allocatable, dimension(:) :: iVMR
  real(kind=dp) :: Ray_xsec
  logical :: first_call = .True.

  ! Output array
  real(kind=dp), allocatable, dimension(:) :: Ray_out
  real(kind=sp), allocatable, dimension(:) :: Ray_write

  ! Namelist variables
  integer :: iopts

  namelist /Rayleigh_nml/ iopts

  private :: Ray_xsec_calc
  public :: calc_Ray_table, refrac_index_calc

contains

  !! Driver routine for Rayleigh scattering table calculation

  subroutine calc_Ray_table()
    implicit none

    integer :: i, j, l, z, s
    logical :: exists
    real(kind=dp) :: Ray_H2O

    ! Allocate work arrays
    allocate(Ray_work(nRay), n_ref(nRay), iVMR(nRay), King(nRay), nd_stp(nRay))
    Ray_work(:) = 0.0_dp

    ! Allocate sp CMCRT output array
    allocate(Ray_out(nlay),Ray_write(nlay))

    ! Read Rayleigh namelist parameters
    read(u_nml, nml=Rayleigh_nml)

    ! Find the prf VMR_indexes of the Rayleigh species
    do i = 1, nRay
      exists = .False.
      do j = 1, ngas
        if (Ray_name(i) == g_name(j)) then
          iVMR(i) = j
          exists = .True.
          exit
        end if
      end do
      if (exists .eqv. .False.) then
        print*, 'ERROR - Specifed Rayleigh species not found in prf VMR list - STOPPING'
        print*, 'Species: ', Ray_name(i)
        stop
      end if
    end do

    print*, ' ~~ Performing Rayleigh calculation and output ~~ '
    print*, ' ~~ Please wait... ~~ '

    !! Begin openMP loops
    !$omp parallel default (none), &
    !$omp& private (l,z,s,Ray_H2O), &
    !$omp& shared (nwl,nlay,nRay,Ray_out,VMR_lay,iVMR,N_lay,RH_lay,Ray_work,wl,Ray_name)

    ! Perform Rayleigh scattering cross section calculation
    ! Species loops are inside subroutines
    do l = 1, nwl
      !$omp single
      if (mod(l,nwl/10) == 0) then
        print*, l, wl(l), nwl
      end if
      !$omp end single
      ! Find the refractive index and King factor for this wavelength
      call refrac_index_calc(l)
      ! Find the cross section for each species for this wavelength
      call Ray_xsec_calc(l)

      !$omp do schedule (static)
      do z = 1, nlay
        ! Find the Rayleigh cross section opacity for this layer
        Ray_out(z) = 0.0_dp
        do s = 1, nRay
          ! Special routine for H2O Rayleigh cross-sections
          if (Ray_name(s) == 'H2O') then
            call Ray_xsec_calc_H2O(s,l,z,Ray_H2O)
            Ray_out(z) = Ray_out(z) + VMR_lay(iVMR(s),z) * N_lay(z) * Ray_H2O
            cycle
          end if
          ! Sum species cross sections and *VMR*N_lay to convert from cm2 molecule-1 to cm-1 of atmosphere
          Ray_out(z) = Ray_out(z) + VMR_lay(iVMR(s),z) * N_lay(z) * Ray_work(s)
        end do
        ! Convert to cm2 g-1 of atmosphere
        Ray_out(z) = Ray_out(z)/RH_lay(z)
      end do
      !$omp end do

      !$omp single
      ! Output CMCRT formatted Rayleigh scattering table for layers
      call output_Ray_table(l)
      !$omp end single

    end do
    !$omp end parallel



    ! Deallocate arrays on exit
    deallocate(Ray_work, n_ref, iVMR, King)
    deallocate(Ray_out,Ray_write)
    ! Close uRay i/o unit
    close(uRay)

    print*, ' ~~ Quest completed  ~~ '

  end subroutine calc_Ray_table

  !! Find the refractive index of species

  subroutine refrac_index_calc(l)
    implicit none

    integer, intent(in) :: l
    integer :: s
    real(kind=dp) :: A, B, C, DPol

    do s = 1, nRay
      ! Find refractive index of species by name
      select case(Ray_name(s))

      case('H2')
        ! Use Dalgarno & Williams (1962) expression
        ! n_ref(s) = -1
        ! King(s) = 1.0_dp

        ! Use Irwin (2009)+ paramaters (same as Cox 2000)
        A = 13.58e-5_dp ; B = 7.52e-3_dp
        n_ref(s) = n_func2(wl(l),A,B)

        DPol = 0.0221_dp
        !King(s) = King_from_Dpol(DPol)
        King(s) = 1.0_dp !King_from_Dpol(DPol)

        nd_stp(s) = 2.65163e19_dp

      case('H')

        ! Can use Kurucz (1970), Ferland (2001) or Lee & Kim (2004) expressions
        n_ref(s) = -1
        King(s) = 1.0_dp

      case('He')
        ! Use Kurucz (1970) expression
        !n_ref(s) = -1
        !King(s) = 1.0_dp

        ! Use Irwin (2009)+ paramaters
        !A = 3.48e-5_dp ; B = 2.3e-3_dp ; DPol = 0.025
        !n_ref(s) = n_func2(wl,A,B)
        !King(s) = King_from_Dpol(DPol)

        ! Use Thalman et al. (2014) parameters
        A = 2283.0_dp ; B = 1.8102e13_dp ; C = 1.5342e10_dp
        n_ref(s) = n_func(wn(l),A,B,C)
        King(s) = 1.0_dp

        nd_stp(s) = 2.546899e19_dp

      case('e-','el')
        ! Use Thomson scattering cross section
        n_ref(s) = -1
        King(s) = 1.0_dp

      case('CO')
        ! Use Sneeps & Ubachs (2005) expression
        ! Typo error in Sneeps & Ubachs corrected (Kitzmann per.com.)
        A = 22851.0_dp; B = 0.456e14_dp ; C = 71427.0**2
        n_ref(s) = n_func(wn(l),A,B,C)
        DPol = 0.0048_dp
        !King(s) = King_from_Dpol(DPol)
        King(s) = 1.0_dp

        nd_stp(s) = 2.546899e19_dp

      case('CO2')
        ! Use Sneeps & Ubachs (2005) expression
        ! Error in Sneeps & Ubachs last term is 0.1218145e−6 not 0.1218145e−4 (See Kitzmann 2017 p 4 foootnote)
        ! Error also in the multiplication factor (See Kitzmann Rayleigh pdf)
        n_ref(s) = 1.1427e3_dp * ((5799.25_dp / (128908.9_dp**2 - wn(l)**2)) + (120.05_dp / (89223.8_dp**2 - wn(l)**2)) &
          & + (5.3334_dp / (75037.5_dp**2 - wn(l)**2)) + (4.3244_dp / (67837.7_dp**2 - wn(l)**2)) &
          & + (0.1218145e-6_dp / (2418.136_dp**2 - wn(l)**2)))
        n_ref(s) = n_ref(s) + 1.0_dp
        King(s) = 1.1364_dp + 25.3e-12_dp*wn(l)**2

        nd_stp(s) = 2.546899e19_dp

      case('CH4')
        ! Use Sneeps & Ubachs (2005) expression
        n_ref(s) = (46662.0_dp + 4.02e-6_dp*wn(l)**2)/1e8_dp + 1.0_dp
        King(s) = 1.0_dp

        nd_stp(s) = 2.546899e19_dp

      case('H2O')

        ! Is a special species - just calculate King factor for now
        Dpol = 3.0e-4_dp ! Murphy (1977)
        King(s) = King_from_Dpol_1(DPol)

      case('O2')

        ! Use Thalman et al. (2014) parameters
        A = 20564.8_dp ; B = 2.480899e13_dp ; C = 4.09e9_dp
        n_ref(s) = n_func(wn(l),A,B,C)
        King(s) = 1.09_dp + 1.385e-11_dp*wn(l)**2 + 1.448e-20_dp*wn(l)**4

        nd_stp(s) = 2.68678e19_dp

      case('N2')

        ! Use Thalman et al. (2014) parameters -
        ! note typo in Thalman et al. (2014) e13 error - used Sneep & Ubachs (2005) expression
        if (wn(l) > 21360.0_dp) then
          A = 5677.465_dp ; B = 318.81874e12_dp ; C = 14.4e9_dp
          n_ref(s) = n_func(wn(l),A,B,C)
        else ! wn(l) < 21360.0
          A = 6498.2_dp ; B = 307.4335e12_dp ; C = 14.4e9_dp
          n_ref(s) = n_func(wn(l),A,B,C)
        end if

        King(s) = 1.034_dp + 3.17e-12_dp*wn(l)

        nd_stp(s) = 2.546899e19_dp

     case('NH3')

        ! Use Irwin (2009)+ paramaters
        A = 37.0e-5_dp ; B = 12.0e-3_dp ; Dpol = 0.0922_dp
        n_ref(s) = n_func2(wl(l),A,B)
        King(s) = King_from_Dpol_1(DPol)

        nd_stp(s) = N_stp

      case('Ar')

        ! Use Thalman et al. (2014) parameters
        A = 6432.135_dp ; B = 286.06021e12_dp ; C = 14.4e9_dp
        n_ref(s) = n_func(wn(l),A,B,C)
        King(s) = 1.0_dp

        nd_stp(s) = 2.546899e19_dp

      case('N2O')
        ! Use Sneeps & Ubachs (2005) expression
        n_ref(s) = (46890.0_dp + 4.12e-6_dp*wn(l)**2)/1e8_dp + 1.0_dp
        Dpol = 0.0577_dp + 11.8e-12_dp*wn(l)**2
        King(s) = King_from_Dpol_2(DPol) 

        nd_stp(s) = 2.546899e19_dp

      case('SF6')
        ! Use Sneeps & Ubachs (2005) expression
        n_ref(s) = (71517.0_dp + 4.996e-6_dp*wn(l)**2)/1e8_dp + 1.0_dp
        King(s) = 1.0_dp

        nd_stp(s) = 2.546899e19_dp

      case default
        print*, 'ERROR - Rayleigh species not found in refrace_index_calc - STOPPING'
        print*, 'Species: ', Ray_name(s)
        stop
      end select

    end do

  contains

    real(kind=dp) function n_func(wn,A,B,C)
      implicit none

      real(kind=dp), intent(in) :: wn, A, B, C
      real(kind=dp) :: nm1

      nm1 = A + (B / (C - wn**2))
      n_func = nm1/1.0e8_dp + 1.0_dp

    end function n_func

    real(kind=dp) function n_func2(wl,A,B)
      implicit none

      real(kind=dp), intent(in) :: wl, A, B
      real(kind=dp) :: nm1

      nm1 = A * (1.0_dp + B/wl**2)
      n_func2 = nm1 + 1.0_dp

    end function n_func2

    real(kind=dp) function King_from_Dpol_1(DPol)
      implicit none

      real(kind=dp), intent(in) :: DPol

      King_from_Dpol_1 = (6.0_dp  + 3.0_dp * DPol) / (6.0_dp - 7.0_dp * DPol)

    end function King_from_Dpol_1

    real(kind=dp) function King_from_Dpol_2(DPol)
      implicit none

      real(kind=dp), intent(in) :: DPol

      King_from_Dpol_2 = (3.0_dp  + 6.0_dp * DPol) / (3.0_dp - 4.0_dp * DPol)

    end function King_from_Dpol_2

  end subroutine refrac_index_calc

  subroutine Ray_xsec_calc(l)
    implicit none

    integer, intent(in) :: l
    integer :: s
    real(kind=dp) :: Ray_spec

    ! Find the cross section for each species following Sneep & Ubachs (2005)
    do s = 1, nRay

      if (n_ref(s) < 0) then
        ! Special cases
        call Ray_xsec_special(l, s, Ray_spec)
        Ray_work(s) = Ray_spec
      else
        ! Normal cases
        Ray_work(s) = ((24.0_dp * pi**3 * wn(l)**4)/(nd_stp(s)**2)) &
          & * ((n_ref(s)**2 - 1.0_dp)/(n_ref(s)**2 + 2.0_dp))**2  * King(s)
      end if

      Ray_work(s) = max(Ray_work(s),1.0e-99_dp)

    end do

  end subroutine Ray_xsec_calc

  subroutine Ray_xsec_special(l,s,xsec)
    implicit none

    integer, intent(in) :: l, s
    real(kind=dp), intent(out) :: xsec

    integer :: p
    real(kind=dp) :: w, wwl, wb

    select case(Ray_name(s))

    case('H2')
      ! Dalgarno & Williams (1962)
      xsec = 8.14e-13_dp/wl_A(l)**4 + 1.28e-6_dp/wl_A(l)**6 + 1.61_dp/wl_A(l)**8

    case('He')
      ! Kurucz (1970)
      xsec = 5.484e-14_dp/wl_A(l)**4 * (1.0_dp + 2.44e5_dp/wl_A(l)**2 &
        & + (5.94e-10_dp/(wl_A(l)**2 * (wl_A(l)**2 - 2.950e5_dp))))**2

    case('H')

      w = 2.0_dp * pi * freq(l)
      wwl = w/w_l

      !! Choose expression, all three are very similar

      ! Kurucz (1970)
      !xsec = 5.799e-13_dp/wl_A(l)**4 + 1.422e-6_dp/wl_A(l)**6 + 2.784_dp/wl_A(l)**8

      ! Ferland (2001)
      !xsec = 8.41e-25_dp*(wwl)**4 + 3.37e-24*(wwl)**6 + 4.71e-22_dp*(wwl)**14

      ! Lee and Kim (2004)
      if (wwl <= 0.6_dp) then
        ! Low energy limit
        xsec = 0.0_dp
        do p = 0, 9
          xsec = xsec + (cp(p+1) * wwl**(2 * p))
        end do
        xsec = xsec * wwl**4
      else
        ! High energy limit (approaching Lyman alpha wavelengths)
        wb = (w - 0.75_dp*w_l)/(0.75_dp*w_l)
        xsec = (0.0433056_dp/wb**2)*(1.0_dp - 1.792_dp*wb - 23.637_dp*wb**2 - 83.1393_dp*wb**3 &
        & - 244.1453_dp*wb**4 - 699.473_dp*wb**5)
      end if
      ! Multiply by Thomson x-section
      xsec = xsec * sigT

    case('e-','el')
      ! Thomson scattering cross section expression
      xsec = sigT

    case default
      print*, 'ERROR - Rayleigh species not found in Ray_xsec_special - STOPPING'
      print*, 'Species: ', Ray_name(s)
      stop

    end select

  end subroutine Ray_xsec_special

  subroutine Ray_xsec_calc_H2O(s,l,z,Ray_H2O)
    implicit none

    integer, intent(in) :: s, l, z
    real(kind=dp), intent(out) :: Ray_H2O

    real(kind=dp) :: sig, lam, theta, A, n_H2O

    ! Validity checks to avoid singularities
    if (wl(l) <= 0.2_dp .or. wl(l) >= 2.5_dp) then
      Ray_H2O = 0.0_dp
    !else if (TG_lay(z) <= 261.15_dp .or. TG_lay(z) >= 773.15_dp) then
    !  Ray_H2O = 0.0_dp
    !else if (RH_lay(z) <= 0.0_dp .or. RH_lay(z) >= 1.060_dp) then
    !  Ray_H2O = 0.0_dp
    else

      sig = (RH_lay(z)*1e3_dp)/rho_s
      lam = wl(l)/lam_s
      theta = TG_lay(z)/T_s

      A = sig * (a_i(1) + a_i(2)*sig + a_i(3)*theta + a_i(4)*lam**2*theta + a_i(5)/lam**2 + &
      &  a_i(6)/(lam**2 - lam_uv**2) + a_i(7)/(lam**2 - lam_ir**2) + a_i(8)*sig**2)

      n_H2O = sqrt((2.0_dp * A + 1.0_dp)/(1.0_dp - A))

      Ray_H2O = ((24.0_dp * pi**3 * wn(l)**4)/(N_lay(z)**2)) &
      & * ((n_H2O**2 - 1.0_dp)/(n_H2O**2 + 2.0_dp))**2  * King(s)

    end if

    Ray_H2O = max(Ray_H2O,1.0e-99_dp)


  end subroutine Ray_xsec_calc_H2O


  subroutine output_Ray_table(l)
    implicit none

    integer, intent(in) :: l
    integer :: z, reclen

    if (first_call .eqv. .True.) then
      inquire(iolength=reclen) Ray_write
      !print*, 'Outputing Rayleigh.cmcrt'
      ! Output k-table in 1D or flattened 3D CMCRT format k_CMCRT.ktb (single precision)
      open(newunit=uRay, file='Rayleigh.cmcrt', action='readwrite',&
        & form='unformatted',status='replace', access='direct',recl=reclen)
      first_call = .False.
    end if

    ! Convert to single precision on output, also care for underfloat
    Ray_write(:) = real(max(Ray_out(:),1.0e-30_dp),kind=sp)
    write(uRay,rec=l) Ray_write

  end subroutine output_Ray_table

end module Ray_tables_mod
