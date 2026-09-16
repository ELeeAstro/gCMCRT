module mc_k_refraction
  use mc_precision
  implicit none

  private

  integer, parameter, public :: REFRACT_SUCCESS = 0
  integer, parameter, public :: REFRACT_TIR = 1
  integer, parameter, public :: REFRACT_INVALID_INPUT = -1

  public :: refract_radial_direction, radial_snell_invariant

contains

  ! ieee_arithmetic exposes separate host and device specifics in NVFORTRAN.
  ! Calling its generic predicates from an attributes(host,device) routine is
  ! therefore ambiguous.  This ordered comparison is false for both NaN and
  ! infinity and is valid in host and device code.
  attributes(host,device) logical function is_finite_dp(value)
    implicit none

    real(dp), intent(in) :: value

    is_finite_dp = abs(value) <= huge(value)

  end function is_finite_dp

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

    if (.not. is_finite_dp(x) .or. .not. is_finite_dp(y) .or. &
        .not. is_finite_dp(z) .or. .not. is_finite_dp(ux) .or. &
        .not. is_finite_dp(uy) .or. .not. is_finite_dp(uz) .or. &
        .not. is_finite_dp(n_old) .or. .not. is_finite_dp(n_new)) return
    if (n_old <= 0.0_dp .or. n_new <= 0.0_dp) return

    radius2 = x*x + y*y + z*z
    direction2 = ux*ux + uy*uy + uz*uz
    if (.not. is_finite_dp(radius2) .or. radius2 <= 0.0_dp) return
    if (.not. is_finite_dp(direction2) .or. direction2 <= 0.0_dp) return

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
    if (.not. is_finite_dp(sin2_transmitted)) return

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
    if (.not. is_finite_dp(out2) .or. out2 <= 0.0_dp) then
      status = REFRACT_INVALID_INPUT
      return
    end if
    out_norm = sqrt(out2)
    outx = outx/out_norm
    outy = outy/out_norm
    outz = outz/out_norm
    if (.not. is_finite_dp(outx) .or. .not. is_finite_dp(outy) .or. &
        .not. is_finite_dp(outz)) then
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

    if (.not. is_finite_dp(x) .or. .not. is_finite_dp(y) .or. &
        .not. is_finite_dp(z) .or. .not. is_finite_dp(ux) .or. &
        .not. is_finite_dp(uy) .or. .not. is_finite_dp(uz) .or. &
        .not. is_finite_dp(refractive_index)) return
    if (refractive_index <= 0.0_dp) return

    radius2 = x*x + y*y + z*z
    direction2 = ux*ux + uy*uy + uz*uz
    if (.not. is_finite_dp(radius2) .or. radius2 <= 0.0_dp) return
    if (.not. is_finite_dp(direction2) .or. direction2 <= 0.0_dp) return

    dot_position_direction = (x*ux + y*uy + z*uz)/sqrt(direction2)
    impact2 = radius2-dot_position_direction*dot_position_direction
    if (.not. is_finite_dp(impact2)) return

    invariant = refractive_index*sqrt(max(0.0_dp,impact2))

  end function radial_snell_invariant

end module mc_k_refraction
