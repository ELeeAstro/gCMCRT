module mc_class_pac
  use mc_precision
  use curand_device
  implicit none

  private
  public :: pac

  type pac
    integer :: geo                         ! Geometry of the grid (1 = cartesian, 2 = spherical)
    integer(8) :: id                          ! Photon number unique id
    integer(8) :: seed                        ! Unique random seed of packet
    integer :: p_flag                      ! Integer flag giving state of packet
    !type(curandStateXORWOW) :: iseed       ! the current iseed integer of the packet (changes each random call)
    integer :: ig
    type(curandStateMRG32k3a) :: iseed
    !type(curandStateMtgp32) :: iseed
    !type(curandStatePhilox4_32_10) :: iseed
    ! integer :: iseed
    ! real(dp), dimension(97) :: u
    ! integer :: i
    ! integer :: j
    ! real(dp) :: c_ran
    real(dp) :: wght                       ! Packet weight
    real(dp) :: wl                         ! Wavelength of packet
    integer, dimension(3) :: c             ! Current X,Y,Z index of cell
    integer :: iscatt                      ! Integer flag denoting type of scattering to undergo
    real(dp) :: rp, rp2                    ! Radial position (spherical grid) and position squared
    real(dp) :: bp                         ! Impact paramater position
    integer :: b_idx
    real(dp) :: xp, yp, zp                 ! Geometric x,y,z position
    real(dp) :: nxp, nyp, nzp              ! Directional x,y,z
    real(dp) :: sint, cost, sinp, cosp     ! Directional theta,phi cosines
    real(dp) :: phi                        ! Packet frame phi angle
    real(dp) :: fi, fq, fu, fv             ! Stokes polarisation state
    real(dp) :: tau_p                      ! Optical depth sampled for packet
    real(dp) :: tau                        ! Optical depth experienced
    real(dp) :: d                          ! Distance traveled by packet
    real(dp) :: en                         ! Energy of packet

  end type pac

end module mc_class_pac
