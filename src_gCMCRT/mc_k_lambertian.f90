module mc_k_lambertian
  use, intrinsic :: iso_fortran_env
  use cudafor
  use curand_device
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp * pi
  real(dp), parameter :: tiny_norm = 1.0e-30_dp

contains

  attributes(device) subroutine normalise3(x, y, z)
    real(dp), intent(inout) :: x, y, z
    real(dp) :: norm

    norm = sqrt(x*x + y*y + z*z)

    if (norm > tiny_norm) then
      x = x / norm
      y = y / norm
      z = z / norm
    else
      x = 0.0_dp
      y = 0.0_dp
      z = 1.0_dp
    end if
  end subroutine normalise3


  attributes(device) subroutine make_basis_from_normal(nx, ny, nz, &
                                                       e1x, e1y, e1z, &
                                                       e2x, e2y, e2z)
    real(dp), intent(in)  :: nx, ny, nz
    real(dp), intent(out) :: e1x, e1y, e1z
    real(dp), intent(out) :: e2x, e2y, e2z

    real(dp) :: hx, hy, hz
    real(dp) :: nxx, nyy, nzz

    nxx = nx
    nyy = ny
    nzz = nz
    call normalise3(nxx, nyy, nzz)

    ! Choose helper vector not close to parallel with n.
    if (abs(nzz) < 0.9_dp) then
      hx = 0.0_dp
      hy = 0.0_dp
      hz = 1.0_dp
    else
      hx = 1.0_dp
      hy = 0.0_dp
      hz = 0.0_dp
    end if

    ! e1 = helper x n
    e1x = hy*nzz - hz*nyy
    e1y = hz*nxx - hx*nzz
    e1z = hx*nyy - hy*nxx

    call normalise3(e1x, e1y, e1z)

    ! e2 = n x e1
    e2x = nyy*e1z - nzz*e1y
    e2y = nzz*e1x - nxx*e1z
    e2z = nxx*e1y - nyy*e1x

    call normalise3(e2x, e2y, e2z)
  end subroutine make_basis_from_normal


  attributes(device) subroutine sample_lambertian_dir(nx, ny, nz, iseed, &
                                                      kx, ky, kz, &
                                                      mu_local)
    real(dp), intent(in) :: nx, ny, nz
    type(curandStateMRG32k3a), intent(inout) :: iseed

    real(dp), intent(out) :: kx, ky, kz
    real(dp), intent(out) :: mu_local

    real(dp) :: e1x, e1y, e1z
    real(dp) :: e2x, e2y, e2z
    real(dp) :: nxx, nyy, nzz

    real(dp) :: u1, u2
    real(dp) :: phi, cphi, sphi
    real(dp) :: sin_local

    nxx = nx
    nyy = ny
    nzz = nz
    call normalise3(nxx, nyy, nzz)

    call make_basis_from_normal(nxx, nyy, nzz, e1x, e1y, e1z, e2x, e2y, e2z)

    ! Lambertian cosine-law:
    ! p(Omega) = mu / pi
    ! p(mu) = 2 mu
    ! CDF(mu) = mu^2
    ! therefore mu = sqrt(U)
    u1 = curand_uniform(iseed)
    u2 = curand_uniform(iseed)

    ! Defensive clipping, because some RNGs can technically return 1.
    u1 = min(max(u1, 0.0_dp), 1.0_dp)
    u2 = min(max(u2, 0.0_dp), 1.0_dp)

    mu_local = sqrt(u1)
    sin_local = sqrt(max(1.0_dp - mu_local*mu_local, 0.0_dp))

    phi = twopi * u2
    cphi = cos(phi)
    sphi = sin(phi)

    kx = mu_local*nxx + sin_local*cphi*e1x + sin_local*sphi*e2x
    ky = mu_local*nyy + sin_local*cphi*e1y + sin_local*sphi*e2y
    kz = mu_local*nzz + sin_local*cphi*e1z + sin_local*sphi*e2z

    call normalise3(kx, ky, kz)
  end subroutine sample_lambertian_dir


  attributes(device) function lambertian_peeloff_pdf(nx, ny, nz, ox, oy, oz) result(wfac)
    real(dp), intent(in) :: nx, ny, nz
    real(dp), intent(in) :: ox, oy, oz
    real(dp) :: wfac
    real(dp) :: nxx, nyy, nzz
    real(dp) :: oxx, oyy, ozz
    real(dp) :: mu_obs

    nxx = nx
    nyy = ny
    nzz = nz
    call normalise3(nxx, nyy, nzz)

    oxx = ox
    oyy = oy
    ozz = oz
    call normalise3(oxx, oyy, ozz)

    mu_obs = nxx*oxx + nyy*oyy + nzz*ozz

    if (mu_obs > 0.0_dp) then
      wfac = mu_obs / pi
    else
      wfac = 0.0_dp
    end if
  end function lambertian_peeloff_pdf

end module mc_k_lambertian