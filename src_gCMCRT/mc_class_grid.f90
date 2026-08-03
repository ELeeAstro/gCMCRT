module mc_class_grid
  use mc_precision
  use mc_data_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none

  type geo

    !! Cartesian variables
    integer :: n_x, n_y, n_z
    real(dp) :: x_min, x_max, y_min, y_max, z_min, z_max

    !! Spherical variables
    integer :: n_lev, n_lay, n_cell
    integer :: n_theta, n_phi
    real(dp) :: r_max, r_min, r2_max
    real(dp) :: r_del

    !! Common variables
    real(dp) :: lumtot, lumtot_surf

  end type geo

  type(geo) :: grid        ! Host grid data
  type(geo), device :: grid_d  ! Device grid data

  !! 1D profile number of lines
  integer :: nlay_1D
  !!  dphi spacing, dtheta spacing
  real(dp) :: dphi, dtheta

  !! mut for Gibbs sampling of Draine phase function
  real(dp) :: mut
  real(dp), device :: mut_d

  !! Cell dimension properties
  real(dp), allocatable, dimension(:,:) :: a_surf
  real(dp), allocatable, dimension(:,:), device :: a_surf_d
  real(dp), allocatable, dimension(:,:,:) :: v_cell, a_cell
  real(dp), allocatable, dimension(:,:,:), device :: v_cell_d, a_cell_d

  !! Allocatable arrays for 1D grid
  integer, dimension(:), allocatable :: lay
  real(dp), dimension(:), allocatable :: PG_1D, MOL_W_1D, TG_1D, RH_1D
  real(dp), dimension(:), allocatable :: u_wind_1D, v_wind_1D, w_wind_1D
  !real(dp), dimension(:,:), allocatable :: VMR_1D
  real(dp), dimension(:), allocatable :: rhokap_1D, ssa_1D, g_1D

  !! Allocatabe arrays for 3D grid
  real(dp), dimension(:,:,:), allocatable :: PG, MOL_W, TG, RH
  real(dp), dimension(:,:,:), allocatable :: u_wind, v_wind, w_wind
  real(dp), dimension(:,:,:,:), allocatable :: v_los
  !real(dp), dimension(:,:,:,:), allocatable :: VMR

  !! Device arrays for 1D grid
  real(dp), dimension(:), allocatable, device :: rhokap_1D_d, ssa_1D_d, g_1D_d


  !! Alloctable Arrays for 3D grid
  real(dp), dimension(:), allocatable :: x, y, z, r, H
  real(dp), dimension(:), allocatable :: thetarr, phiarr, tan2thet, aarr, barr
  real(dp), dimension(:,:,:,:), allocatable :: rhokap, ssa
  real(dp), dimension(:,:,:), allocatable :: gg, Dgg
  real(dp), dimension(:,:,:), allocatable :: cf
  real(dp), dimension(:,:,:), allocatable :: l_cell
  real(dp), dimension(:,:,:,:), allocatable :: l_cell_g
  real(dp), dimension(:,:,:), allocatable :: dorg

  ! Surface properties
  real(dp), allocatable, dimension(:,:) :: surf_ssa

  !! Device arrays for 3D grid
  real(dp), dimension(:), allocatable, device :: x_d, y_d, z_d, r_d, H_d
  real(dp), dimension(:), allocatable, device :: theta_d, phi_d, tan2thet_d, aarr_d, barr_d
  real(dp), dimension(:,:,:,:), allocatable, device :: rhokap_d, ssa_d
  real(dp), dimension(:,:,:), allocatable, device :: g_d, Dgg_d
  real(dp), dimension(:,:,:), allocatable, device :: cf_d
  real(dp), dimension(:,:,:), allocatable, device :: dorg_d



  ! Cloud opacities and properties
  real(dp), allocatable, dimension(:,:,:) :: cld_ext ! cm^2 g^-1
  real(dp), allocatable, dimension(:,:,:) :: cld_g
  real(dp), allocatable, dimension(:,:,:) :: cld_ssa

  ! Cloud opacities and properties
  real(dp), allocatable, dimension(:,:,:), device  :: cld_ext_d ! cm^2 g^-1
  real(dp), allocatable, dimension(:,:,:), device  :: cld_g_d
  real(dp), allocatable, dimension(:,:,:), device  :: cld_ssa_d

  ! Surface properties
  real(dp), allocatable, dimension(:,:), device  :: surf_ssa_d


  integer, allocatable, dimension(:) :: u_cf

contains

  subroutine validate_spherical_grid_dimensions()
    implicit none

    if (grid%n_lay < 1) then
      print*, 'ERROR: spherical grid requires at least one radial layer.'
      print*, 'n_lay: ', grid%n_lay
      stop
    end if

    if (grid%n_lev /= grid%n_lay + 1) then
      print*, 'ERROR: spherical grid requires n_lev = n_lay + 1.'
      print*, 'n_lay, n_lev: ', grid%n_lay, grid%n_lev
      stop
    end if

    if (grid%n_theta < 3) then
      print*, 'ERROR: spherical grid requires at least three theta boundaries.'
      print*, 'n_theta: ', grid%n_theta
      stop
    end if

    if (mod(grid%n_theta,2) == 0) then
      print*, 'ERROR: spherical grid requires an odd n_theta.'
      print*, 'The equator must be a theta boundary, not the centre of a cell.'
      print*, 'n_theta: ', grid%n_theta
      stop
    end if

    if (grid%n_phi < 4) then
      print*, 'ERROR: spherical grid requires at least four phi boundaries.'
      print*, 'Smaller grids have coincident longitude-face planes.'
      print*, 'n_phi: ', grid%n_phi
      stop
    end if

  end subroutine validate_spherical_grid_dimensions

  subroutine set_grid()
    implicit none

    integer :: i, j, k, equator_idx
    real(dp) :: v_tot, dr3, a_tot
    real(dp) :: dphi_cell, dmu_cell

    print*, 'Setting grid'

    call validate_spherical_grid_dimensions()

    if (.not. allocated(H)) then
      print*, 'ERROR: height boundaries must be allocated before setting the grid.'
      stop
    end if

    if (size(H) /= grid%n_lev) then
      print*, 'ERROR: height-boundary count does not match grid%n_lev.'
      print*, 'Height count, n_lev: ', size(H), grid%n_lev
      stop
    end if

    if (any(.not. ieee_is_finite(H))) then
      print*, 'ERROR: spherical height boundaries must be finite.'
      stop
    end if

    if (H(1) <= 0.0_dp) then
      print*, 'ERROR: the inner spherical radius must be positive.'
      print*, 'H(1): ', H(1)
      stop
    end if

    if (any(H(2:grid%n_lev) <= H(1:grid%n_lev-1))) then
      print*, 'ERROR: spherical height boundaries must be strictly increasing.'
      stop
    end if

    if (oneD .eqv. .True.) then
      grid%n_cell = grid%n_lay
    else if (threeD .eqv. .True.) then
      grid%n_cell = grid%n_lay * (grid%n_theta - 1) * (grid%n_phi - 1)
    end if

    allocate(thetarr(grid%n_theta),phiarr(grid%n_phi),r(grid%n_lev))


    grid%r_min = 1.0_dp
    grid%r_max = H(grid%n_lev)/H(1)
    grid%r2_max = grid%r_max**2

    r(1) = grid%r_min
    do i = 2, grid%n_lev-1
      r(i) = H(i)/H(1)
      !print*, i, r(i)
    end do
    r(grid%n_lev) =  grid%r_max
    rimage = r(grid%n_lev) * 1.01_dp

    grid%r_del = (H(grid%n_lev) - H(1))/(r(grid%n_lev) - r(1))

    !! Phi (longitude) grid set up.
    allocate(aarr(grid%n_phi),barr(grid%n_phi))
    dphi = twopi / real(grid%n_phi-1,kind=dp)
    do j = 1, grid%n_phi
      phiarr(j) = real(j-1,kind=dp) * dphi
      aarr(j) = sin(phiarr(j))
      barr(j) = -cos(phiarr(j))
    end do
    phiarr(1) = 0.0_dp
    phiarr(grid%n_phi) = twopi
    aarr(1) = 0.0_dp
    barr(1) = -1.0_dp
    aarr(grid%n_phi) = 0.0_dp
    barr(grid%n_phi) = -1.0_dp
    if (mod(grid%n_phi-1,2) == 0) then
      j = (grid%n_phi + 1) / 2
      phiarr(j) = pi
      aarr(j) = 0.0_dp
      barr(j) = 1.0_dp
    end if

    !! Theta grid set up. Odd n_theta puts one boundary exactly at pi/2.
    dtheta = pi / real(grid%n_theta-1,kind=dp)
    allocate(tan2thet(grid%n_theta))
    equator_idx = (grid%n_theta + 1) / 2
    do k = 1, grid%n_theta
      thetarr(k) = real(k-1,kind=dp) * dtheta
      if (k == 1) then
        thetarr(k) = 0.0_dp
        tan2thet(k) = 0.0_dp
      else if (k == grid%n_theta) then
        thetarr(k) = pi
        tan2thet(k) = 0.0_dp
      else if (k == equator_idx) then
        thetarr(k) = 0.5_dp * pi
        tan2thet(k) = -1.0_dp
      else
        tan2thet(k) = tan(thetarr(k))**2
      end if
    end do

    if (count(tan2thet < 0.0_dp) /= 1) then
      print*, 'ERROR: spherical theta grid must contain exactly one equatorial plane.'
      stop
    end if

    ! Find volume of each cell - spherical coordinates
    allocate(v_cell(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
    do i = 1, grid%n_lay
      dr3 = H(i+1)**3 - H(i)**3
      do j = 1, grid%n_phi-1
        dphi_cell = phiarr(j+1) - phiarr(j)
        do k = 1, grid%n_theta-1
          dmu_cell = cos(thetarr(k)) - cos(thetarr(k+1))
          v_cell(i,j,k) = (dr3 * dphi_cell * dmu_cell) / 3.0_dp
          !print*, i, j, k, v_cell(i,j,k), dr3, dphi_cell, dmu_cell
        end do
      end do
    end do

    v_tot = (4.0_dp/3.0_dp) * pi * (H(grid%n_lev)**3 - H(1)**3)
    print*, 'Volumes Check:', v_tot - sum(v_cell), v_tot, sum(v_cell), sum(v_cell)/v_tot

    ! Find area of each cell at the surface - spherical coordinates
    allocate(a_surf(grid%n_phi-1,grid%n_theta-1))
    dr3 = H(1)**2
    do j = 1, grid%n_phi-1
      dphi_cell = phiarr(j+1) - phiarr(j)
      do k = 1, grid%n_theta-1
        dmu_cell = cos(thetarr(k)) - cos(thetarr(k+1))
        a_surf(j,k) = dr3 * dphi_cell * dmu_cell
        !print*, j, k, a_cell(j,k), dr3, dphi_cell, dmu_cell
      end do
    end do

    a_tot = 4.0_dp * pi * H(1)**2
    print*, 'Surface area Check:', a_tot - sum(a_surf), a_tot, sum(a_surf), sum(a_surf)/a_tot

    !! Allocate device arrays and send to gpu
    allocate(r_d(grid%n_lev),theta_d(grid%n_theta),phi_d(grid%n_phi))
    allocate(tan2thet_d(grid%n_theta))
    allocate(aarr_d(grid%n_phi),barr_d(grid%n_phi))
    r_d(:) = r(:)
    phi_d(:) = phiarr(:)
    theta_d(:) = thetarr(:)
    tan2thet_d(:) = tan2thet(:)
    aarr_d(:) = aarr(:)
    barr_d(:) = barr(:)

    if (do_cf .eqv. .True.) then
      allocate(cf(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(cf_d(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      cf(:,:,:) = 0.0_dp
      cf_d(:,:,:) = cf(:,:,:)
    end if

    print*, ' - Complete - '

  end subroutine set_grid


  subroutine output_cf(n,l,n_output)
    implicit none

    integer, intent(in) :: n, l, n_output
    logical, save :: first_call = .True.
    integer :: nn
    character (len=8) :: fmt
    character (len=3) :: n_str

    if (n_output < 1) then
      print*, 'ERROR: contribution-function output count must be positive.'
      print*, 'n_output: ', n_output
      stop
    end if

    if (first_call .eqv. .True.) then
      allocate(u_cf(n_output))
      fmt = '(I3.3)'      
      do nn = 1, n_output
        write(n_str,fmt) nn
        open(newunit=u_cf(nn), file='cf_'//trim(n_str)//'.txt', status='replace', action='write',form='unformatted')
      end do

      first_call = .False.
    end if

    if (size(u_cf) /= n_output) then
      print*, 'ERROR: contribution-function output count changed after initialization.'
      print*, 'Initial, requested counts: ', size(u_cf), n_output
      stop
    end if

    if ((n < 1) .or. (n > size(u_cf))) then
      print*, 'ERROR: contribution-function output index is out of range.'
      print*, 'Index, available outputs: ', n, size(u_cf)
      stop
    end if

    write(u_cf(n)) real(cf(:,:,:))

  end subroutine output_cf

end module mc_class_grid
