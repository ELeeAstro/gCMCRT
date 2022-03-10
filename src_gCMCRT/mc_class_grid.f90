module mc_class_grid
  use mc_precision
  use mc_data_mod
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
    real(dp) :: lumtot

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
  real(dp), allocatable, dimension(:,:,:) :: v_cell, a_cell
  real(dp), allocatable, dimension(:,:,:), device :: v_cell_d, a_cell_d

  !! Allocatable arrays for 1D grid
  integer, dimension(:), allocatable :: lay
  real(dp), dimension(:), allocatable :: PG_1D, MOL_W_1D, TG_1D, RH_1D
  real(dp), dimension(:), allocatable :: u_wind_1D, v_wind_1D, w_wind_1D
  real(dp), dimension(:,:), allocatable :: VMR_1D
  real(dp), dimension(:), allocatable :: rhokap_1D, ssa_1D, g_1D

  !! Allocatabe arrays for 3D grid
  real(dp), dimension(:,:,:), allocatable :: PG, MOL_W, TG, RH
  real(dp), dimension(:,:,:), allocatable :: u_wind, v_wind, w_wind
  real(dp), dimension(:,:,:,:), allocatable :: v_los
  real(dp), dimension(:,:,:,:), allocatable :: VMR

  !! Device arrays for 1D grid
  real(dp), dimension(:), allocatable, device :: rhokap_1D_d, ssa_1D_d, g_1D_d


  !! Alloctable Arrays for 3D grid
  real(dp), dimension(:), allocatable :: x, y, z, r, H
  real(dp), dimension(:), allocatable :: thetarr, phiarr, tan2thet, aarr, barr
  real(dp), dimension(:,:,:,:), allocatable :: rhokap, ssa
  real(dp), dimension(:,:,:), allocatable :: gg, Dgg, Dmut
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
  real(dp), dimension(:,:,:), allocatable, device :: g_d, Dgg_d, Dmut_d
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


contains

  subroutine set_grid()
    implicit none

    integer :: i,j,k
    real(dp) :: v_tot, dr3

    print*, 'Setting grid'

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

    !! phi (longitude) grid set up
    ! 1st grid = 0 degrees
    allocate(aarr(grid%n_phi),barr(grid%n_phi))
    phiarr(1) = 0.0_dp
    aarr(1) = sin(phiarr(1))
    barr(1) = -cos(phiarr(1))

    ! Spacing in longitude
    dphi = twopi / real(grid%n_phi-1,kind=dp)
    do j = 2, grid%n_phi - 1
      phiarr(j) = phiarr(j-1) + dphi
      aarr(j) = sin(phiarr(j))
      barr(j) = -cos(phiarr(j))
      !print*, j, phiarr(j), aarr(j), barr(j), dphi
    end do

    ! last grid = 360 degrees
    phiarr(grid%n_phi) = twopi
    aarr(grid%n_phi) = sin(phiarr(grid%n_phi))
    barr(grid%n_phi) = -cos(phiarr(grid%n_phi))

    ! Spacing in latitude
    dtheta = pi / real(grid%n_theta-1,kind=dp)

    allocate(tan2thet(grid%n_theta))
    ! 1st gird = 0 degrees
    thetarr(1) = 0.0_dp
    tan2thet(1) = 0.0_dp

    do k = 2, grid%n_theta-1
      thetarr(k) = thetarr(k-1) + dtheta
      if ((thetarr(k) > (1.57060_dp)).and.(thetarr(k) < (1.57090_dp))) then
        tan2thet(k) = -1.0_dp
        !tan2thet(k) = tan(thetarr(k))**2
      else
        tan2thet(k) = tan(thetarr(k))**2
      endif
      !print*, k, thetarr(k), tan2thet(k), dtheta
    end do

    ! last grid = pi
    thetarr(grid%n_theta) = pi
    tan2thet(grid%n_theta) = 0.0_dp

    ! Find volume of each cell - spherical coordinates
    allocate(v_cell(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
    do i = 1, grid%n_lay
      dr3 = H(i+1)**3 - H(i)**3
      do j = 1, grid%n_phi-1
        dphi = phiarr(j+1) - phiarr(j)
        do k = 1, grid%n_theta-1
          dtheta = cos(thetarr(k)) - cos(thetarr(k+1))
          v_cell(i,j,k) = (dr3 * dphi * dtheta) / 3.0_dp
          !print*, i, j, k, v_cell(i,j,k), dr3, dphi, dtheta
        end do
      end do
    end do

    v_tot = (4.0_dp/3.0_dp) * pi * (H(grid%n_lev)**3 - H(1)**3)
    print*, 'Volumes Check:', v_tot - sum(v_cell), v_tot, sum(v_cell), sum(v_cell)/v_tot

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
      call output_cf(1)
    end if

    print*, ' - Complete - '

  end subroutine set_grid


  subroutine output_cf(l)
    implicit none

    integer, intent(in) :: l
    logical, save :: first_call = .True.
    integer, save :: u_cf

    if (first_call .eqv. .True.) then
      open(newunit=u_cf, file='cf.out', action='readwrite',form='unformatted')

      first_call = .False.
      return
    end if

    write(u_cf) real(cf(:,:,:))

  end subroutine output_cf

end module mc_class_grid
