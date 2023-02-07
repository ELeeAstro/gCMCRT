module mc_data_mod
  use mc_precision
  implicit none

  !! Common constants
  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp*pi, fourpi = 4.0_dp*pi
  real(dp), parameter :: half = 1.0_dp/2.0_dp
  real(dp), parameter :: third = 1.0_dp/3.0_dp, twothird = 2.0_dp/3.0_dp

  !! Common physical constants - CODATA 2014 (NOTE: Due update in 2018)
  real(dp), parameter :: kb = 1.380649e-16_dp ! erg K-1 - Boltzmann's constant
  real(dp), parameter :: sb_c = 5.670374419e-5_dp ! erg cm-2 s-1 K-4 - Stefan-Boltzmann's constant
  real(dp), parameter :: hpl = 6.62607015e-27_dp ! erg s - Planck's constant
  real(dp), parameter :: c_s = 2.99792458e10_dp ! cm s-1 - Vacuum speed of light
  real(dp), parameter :: amu = 1.660539040e-24_dp ! g - Atomic mass unit
  real(dp), parameter :: N_A = 6.022140857e23_dp ! mol-1 - Avogadro's constant

  !! Constants for BB_band function
!  real(dp), parameter :: kb_mks = 1.380649e-23_dp ! J K-1
!  real(dp), parameter :: h_mks = 6.62607015e-34_dp ! J s
!  real(dp), parameter :: cs_mks = 2.99792458e8_dp ! m s-1
!  real(dp), parameter :: c1 = (hpl * c_s) / kb
 ! real(dp), parameter :: c2 = c_s**2
 ! real(dp), parameter :: n2 = 2.0_dp * hpl * c2

  !! Common astrophysical physical constants - IAU 2015 values
  real(dp), parameter :: Rsun = 6.957e10_dp ! cm
  real(dp), parameter :: Ssun = 1361.0e3_dp ! erg s^-1 cm^-2
  real(dp), parameter :: Lsun = 3.828e33_dp ! erg s^-1
  real(dp), parameter :: Tsun = 5772.0_dp   ! K
  real(dp), parameter :: Rjup = 7.1492e9_dp ! cm - Equatorial radius !6.9911e9 !BENCHMARK VALUE
  real(dp), parameter :: Rearth = 6.3781e8  ! cm - Equatorial radius
  real(dp), parameter :: Au = 1.495978707e13_dp ! cm
  real(dp), parameter :: daysec = 86400.0_dp! s - seconds in a day

  !! Common unit conversions
  real(dp), parameter :: bar = 1e6_dp ! bar to dyne
  real(dp), parameter :: atm = 1.01325e6_dp ! atm to dyne

  !! System parameters from namelist & calculated
  real(dp) :: orbital_period, systemic_velocity, sm_ax
  real(dp), device :: sm_ax_d

  character(len=20) :: exp_name

  !! Switches for command line arguments vs namelists
  logical :: cmd_vphi = .False.
  character(6) :: vphi_arg

  character(50) :: xper


  integer :: u_nml

  logical :: lbl, ck
  logical :: do_infslab, do_diffuse, do_cart_3D

  logical :: oneD, threeD

  logical :: inc_ck, inc_lbl, inc_CIA, inc_Ray, inc_cld, inc_xsec

  logical :: do_images
  logical, device :: do_images_d

  logical :: do_moments
  logical, device :: do_moments_d

  logical :: do_trans
  logical, device :: do_trans_d

  logical :: do_cf
  logical, device :: do_cf_d

  logical :: do_g_bias
  logical, device :: do_g_bias_d

  logical :: do_scat_loop
  logical, device :: do_scat_loop_d

  logical :: wght_deg
  logical, device :: wght_deg_d

  logical :: winds_on, rotation_on, orbit_on, doppler_on

  integer :: xpix, ypix
  real(dp) :: rimage

  logical :: do_Draine
  real(dp) :: Draine_alp
  real(dp), device :: Draine_alp_d

  real(dp), parameter :: xi_g = 0.99_dp
  integer :: ng
  integer, device :: ng_d

  integer, device :: iscat_d
  real(dp), device :: pl_d, pc_d, sc_d

  integer :: n_wl
  real(dp), allocatable, dimension(:) :: wl
  real(dp), allocatable, dimension(:) :: wl_e
  real(dp), allocatable, dimension(:), device :: wl_d

  logical :: do_BB_band

  !! corr-k related arrays
  real(dp), dimension(:), allocatable :: gord_cdf, gord_x, gord_w
  real(dp), dimension(:), allocatable, device :: gord_cdf_d, gord_x_d, gord_w_d
  real(dp), dimension(:,:,:,:), allocatable :: cell_gord_cdf
  real(dp), dimension(:,:,:,:), allocatable, device :: cell_gord_cdf_d
  real(dp), dimension(:,:,:,:), allocatable :: cell_gord_wght
  real(dp), dimension(:,:,:,:), allocatable, device :: cell_gord_wght_d


  !! Transmission optical depth arrays
  real(dp), allocatable, dimension(:) :: T_trans
  real(dp), allocatable, dimension(:), device :: T_trans_d

  !! Albedo spectrum arrays
  real(dp), allocatable, dimension(:) :: alb_out
  real(dp), allocatable, dimension(:), device :: alb_out_d

  !! Limb darkening scheme
  logical :: do_LD
  logical, device :: do_LD_d
  integer :: ilimb
  logical, device :: ilimb_d
  real(dp), dimension(4) :: LD_c
  real(dp), dimension(4), device :: LD_c_d
  real(dp) :: Rs, inc, phase
  real(dp), device :: Rs_d, inc_d, phase_d

  !! Transmission contribution function variables
  integer, allocatable, dimension(:) :: b_n_cf
  integer, allocatable, dimension(:), device :: b_n_cf_d
  real(dp), allocatable, dimension(:) :: b_cf, b_cf_grid
  real(dp), allocatable, dimension(:), device :: b_cf_d, b_cf_grid_d

  real(dp), allocatable, dimension(:) :: p_cf, p_cf_grid
  real(dp), allocatable, dimension(:), device :: p_cf_d, p_cf_grid_d

  integer :: n_phase

  !! Surface properties
  logical :: do_surf
  logical, device :: do_surf_d
  real(dp) :: T_surf, emis_surf, alb_surf
  real(dp), device :: alb_surf_d

  namelist /main/ xper, exp_name, oneD, threeD, do_infslab, do_diffuse, do_cart_3D, do_images, do_moments, do_trans &
    & lbl, ck, orbital_period, sm_ax, systemic_velocity, winds_on, rotation_on, orbit_on, doppler_on, &
    & inc_ck, inc_lbl, inc_CIA, inc_Ray, inc_cld, inc_xsec, do_cf, xpix, ypix, wght_deg, Draine_alp, do_Draine, &
    & do_LD, ilimb, LD_c, Rs, inc, phase, do_g_bias, do_scat_loop, do_BB_band, n_phase, do_surf, T_surf, emis_surf, alb_surf

end module mc_data_mod
