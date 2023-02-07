module optools_data_mod
  use optools_table_class
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Single and Double precision kinds
  integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64

  !! Common constants
  real(kind=dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  real(kind=dp), parameter :: twopi = 2.0_dp*pi, fourpi = 4.0_dp*pi
  real(kind=dp), parameter :: half = 1.0_dp/2.0_dp
  real(kind=dp), parameter :: third = 1.0_dp/3.0_dp, twothird = 2.0_dp/3.0_dp

  !! Common physical constants - CODATA 2020
  real(kind=dp), parameter :: kb = 1.380649e-16_dp ! erg K-1 - Boltzmann's constant
  real(kind=dp), parameter :: sb_c = 5.670374419e-5_dp ! erg cm-2 s-1 K-4 - Stefan-Boltzmann's constant
  real(kind=dp), parameter :: hp = 6.62607015e-27_dp ! erg s - Planck's constant
  real(kind=dp), parameter :: hp_bar = hp/twopi
  real(kind=dp), parameter :: c_s = 2.99792458e10_dp ! cm s-1 - Vacuum speed of light
  real(kind=dp), parameter :: amu = 1.66053906660e-24_dp ! g - Atomic mass unit
  real(kind=dp), parameter :: m_el = 9.1093837015e-28_dp ! g - electron mass
  real(kind=dp), parameter :: N_A = 6.02214076e23_dp ! mol-1 - Avogadro's constant
  real(kind=dp), parameter :: a_fine = 7.2973525693e-3_dp ! fine structure constant
  real(kind=dp), parameter :: Comp_e = hp/(m_el*c_s) ! Compton wavelength for electrons
  real(kind=dp), parameter :: sigT = ((8.0_dp*pi)/3.0_dp) * ((a_fine * Comp_e)/twopi)**2 ! Thomson cross section [cm2]

  !! Common unit conversions
  real(kind=dp), parameter :: pa = 10.0_dp ! Pa to dyne
  real(kind=dp), parameter :: bar = 1.0e6_dp ! bar to dyne
  real(kind=dp), parameter :: atm = 1.01325e6_dp ! atm to dyne

  !! STP Constants
  real(kind=dp), parameter :: T_stp = 273.15_dp ! Standard temperature [K]
  real(kind=dp), parameter :: P_stp = atm ! Standard pressure [dyne]
  real(kind=dp), parameter :: N_stp = P_stp / (kb * T_stp) ! Standard number density [cm-3]

  !--------------- Global paramaters ---------------------

  ! ---- Main Switches ---- !
  logical :: corr_k
  logical :: lbl
  logical :: conti
  logical :: Ray_scat
  logical :: cloud_opc
  logical :: xsec_opc

  ! ---- Name of gas and dust species to calculate ---- !
  character(len=10), allocatable, dimension(:) :: g_name, CK_name, lbl_name, CIA_name, Ray_name, xsec_name
  character(len=10), allocatable, dimension(:) :: cl_name, d_name

  ! ---- Output file units ---- !
  integer :: uCK, ulbl, ulbl_for, uRay, ucl_k, ucl_a, ucl_g, uCIA, u_nml, uxsec

  ! ---- Number of species for each source ---- !
  integer :: nCK, nlbl, nCIA, nRay, ncl, nxsec

  ! ---- Experiment name ---- !
  character(len=50) :: exp_name

  ! ---- Cloud distribution paramaters and options ---- !
  integer :: imix, idist, imie, ndist, idist_int
  real(kind=dp) :: sig, lsig, eff_fac, veff
  real(kind=dp) :: amin, amax, lamin, lamax
  real(kind=dp), allocatable, dimension(:) :: a_dist, la_dist
  real(kind=dp) :: fmax

  !--------------- wavelength arrays---------------------
  integer :: nwl
  integer, allocatable, dimension(:) :: iwl
  real(kind=dp), allocatable, dimension(:) :: wl, wl_cm, wl_A ! Wavelength [um, cm, A]
  real(kind=dp), allocatable, dimension(:) :: wn ! Wavenumber [cm-1]
  real(kind=dp), allocatable, dimension(:) :: freq ! Frequency [Hz]

  !--------------- 1D grid arrays---------------------
  integer :: nlay, ngas, ndust, nmode
  integer, allocatable, dimension(:) :: ilay, ilay2
  real(kind=dp), allocatable, dimension(:) :: PG_lay
  real(kind=dp), allocatable, dimension(:) :: TG_lay
  real(kind=dp), allocatable, dimension(:) :: RH_lay
  real(kind=dp), allocatable, dimension(:) :: N_lay
  real(kind=dp), allocatable, dimension(:) :: mu_lay
  real(kind=dp), allocatable, dimension(:,:) :: VMR_lay
  real(kind=dp), allocatable, dimension(:) :: nd_cl_lay, a_cl_lay, la_cl_lay
  real(kind=dp), allocatable, dimension(:,:) :: VMR_cl_lay
  real(kind=dp), allocatable, dimension(:,:) :: a_C_cl_lay
  real(kind=dp), allocatable, dimension(:,:,:) :: nd_C_cl_lay

  ! --- Global CK table data --- !
  type(CK_table), allocatable, dimension(:) :: CK_tab

  ! --- Global LBL table data --- !
  type(lbl_table), allocatable, dimension(:) :: lbl_tab

  ! --- Global CIA table data --- !
  type(CIA_table), allocatable, dimension(:) :: CIA_tab

  ! --- Global cloud table data --- !
  type(nk_table), allocatable, dimension(:) :: cl_tab

  ! --- Global xsec table data --- !
  type(xsec_table), allocatable, dimension(:) :: xsec_tab

end module optools_data_mod
