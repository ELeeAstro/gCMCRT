module mc_k_refraction
  use mc_precision
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  implicit none

  private

  integer, parameter, public :: REFRACT_SUCCESS = 0
  integer, parameter, public :: REFRACT_TIR = 1
  integer, parameter, public :: REFRACT_INVALID_INPUT = -1

  public :: refract_radial_direction, radial_snell_invariant

contains

  attributes(host,device) subroutine refract_radial_direction(x, y, z, &
      ux, uy, uz, n_old, n_new, status)
    implicit none

    real(dp), intent(in) :: x, y, z, n_old, n_new
    real(dp), intent(inout) :: ux, uy, uz
    integer, intent(out) :: status

    real(dp) :: radius2, radius, direction2, direction_norm
    real(dp) :: rx, ry, rz, mx, my, mz
    real(dp) :: uix, uiy, uiz, cos_incident
    real(dp) :: ratio, sin2_transmitted, cos_transmitted, tir_tolerance
    real(dp) :: tx, ty, tz, outx, outy, outz, out2, out_norm

    status = REFRACT_INVALID_INPUT

    if (ieee_is_nan(x) .or. abs(x) > huge(x) .or. &
        ieee_is_nan(y) .or. abs(y) > huge(y) .or. &
        ieee_is_nan(z) .or. abs(z) > huge(z) .or. &
        ieee_is_nan(ux) .or. abs(ux) > huge(ux) .or. &
        ieee_is_nan(uy) .or. abs(uy) > huge(uy) .or. &
        ieee_is_nan(uz) .or. abs(uz) > huge(uz) .or. &
        ieee_is_nan(n_old) .or. abs(n_old) > huge(n_old) .or. &
        ieee_is_nan(n_new) .or. abs(n_new) > huge(n_new)) return
    if (n_old <= 0.0_dp .or. n_new <= 0.0_dp) return

    radius2 = x*x + y*y + z*z
    direction2 = ux*ux + uy*uy + uz*uz
    if (ieee_is_nan(radius2) .or. abs(radius2) > huge(radius2) .or. &
        radius2 <= 0.0_dp) return
    if (ieee_is_nan(direction2) .or. abs(direction2) > huge(direction2) .or. &
        direction2 <= 0.0_dp) return

    radius = sqrt(radius2)
    direction_norm = sqrt(direction2)
    rx = x/radius
    ry = y/radius
    rz = z/radius
    uix = ux/direction_norm
    uiy = uy/direction_norm
    uiz = uz/direction_norm

    ! Point the normal from the old medium into the new medium.  Its direction
    ! is therefore outward for an outward-moving ray and inward otherwise.
    cos_incident = uix*rx + uiy*ry + uiz*rz
    if (cos_incident >= 0.0_dp) then
      mx = rx
      my = ry
      mz = rz
    else
      mx = -rx
      my = -ry
      mz = -rz
      cos_incident = -cos_incident
    end if
    cos_incident = min(1.0_dp,max(0.0_dp,cos_incident))

    ratio = n_old/n_new
    sin2_transmitted = ratio*ratio * &
      max(0.0_dp,1.0_dp-cos_incident*cos_incident)
    if (ieee_is_nan(sin2_transmitted) .or. &
        abs(sin2_transmitted) > huge(sin2_transmitted)) return

    tir_tolerance = 256.0_dp*epsilon(1.0_dp) * &
      max(1.0_dp,abs(sin2_transmitted))
    if (sin2_transmitted > 1.0_dp+tir_tolerance) then
      ! Specular reflection at the radial interface.  The caller must retain
      ! the old cell and nudge the ray back into that medium.
      outx = uix - 2.0_dp*cos_incident*mx
      outy = uiy - 2.0_dp*cos_incident*my
      outz = uiz - 2.0_dp*cos_incident*mz
      status = REFRACT_TIR
    else
      sin2_transmitted = min(1.0_dp,max(0.0_dp,sin2_transmitted))
      cos_transmitted = sqrt(max(0.0_dp,1.0_dp-sin2_transmitted))

      tx = uix - cos_incident*mx
      ty = uiy - cos_incident*my
      tz = uiz - cos_incident*mz
      outx = ratio*tx + cos_transmitted*mx
      outy = ratio*ty + cos_transmitted*my
      outz = ratio*tz + cos_transmitted*mz
      status = REFRACT_SUCCESS
    end if

    out2 = outx*outx + outy*outy + outz*outz
    if (ieee_is_nan(out2) .or. abs(out2) > huge(out2) .or. &
        out2 <= 0.0_dp) then
      status = REFRACT_INVALID_INPUT
      return
    end if
    out_norm = sqrt(out2)
    outx = outx/out_norm
    outy = outy/out_norm
    outz = outz/out_norm
    if (ieee_is_nan(outx) .or. abs(outx) > huge(outx) .or. &
        ieee_is_nan(outy) .or. abs(outy) > huge(outy) .or. &
        ieee_is_nan(outz) .or. abs(outz) > huge(outz)) then
      status = REFRACT_INVALID_INPUT
      return
    end if

    ux = outx
    uy = outy
    uz = outz

  end subroutine refract_radial_direction


  attributes(host,device) real(dp) function radial_snell_invariant(x, y, z, &
      ux, uy, uz, refractive_index) result(invariant)
    implicit none

    real(dp), intent(in) :: x, y, z, ux, uy, uz, refractive_index

    real(dp) :: radius2, direction2, dot_position_direction
    real(dp) :: impact2

    invariant = -1.0_dp

    if (ieee_is_nan(x) .or. abs(x) > huge(x) .or. &
        ieee_is_nan(y) .or. abs(y) > huge(y) .or. &
        ieee_is_nan(z) .or. abs(z) > huge(z) .or. &
        ieee_is_nan(ux) .or. abs(ux) > huge(ux) .or. &
        ieee_is_nan(uy) .or. abs(uy) > huge(uy) .or. &
        ieee_is_nan(uz) .or. abs(uz) > huge(uz) .or. &
        ieee_is_nan(refractive_index) .or. &
        abs(refractive_index) > huge(refractive_index)) return
    if (refractive_index <= 0.0_dp) return

    radius2 = x*x + y*y + z*z
    direction2 = ux*ux + uy*uy + uz*uz
    if (ieee_is_nan(radius2) .or. abs(radius2) > huge(radius2) .or. &
        radius2 <= 0.0_dp) return
    if (ieee_is_nan(direction2) .or. abs(direction2) > huge(direction2) .or. &
        direction2 <= 0.0_dp) return

    dot_position_direction = (x*ux + y*uy + z*uz)/sqrt(direction2)
    impact2 = radius2-dot_position_direction*dot_position_direction
    if (ieee_is_nan(impact2) .or. abs(impact2) > huge(impact2)) return

    invariant = refractive_index*sqrt(max(0.0_dp,impact2))

  end function radial_snell_invariant

end module mc_k_refraction
