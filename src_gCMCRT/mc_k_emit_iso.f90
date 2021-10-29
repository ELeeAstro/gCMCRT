!! Module containing device level subroutines (GPU) that emit packets isotropically
module mc_k_emit_iso
  use mc_precision
  use mc_class_pac
  use mc_class_grid
  use mc_data_mod
  use mc_k_moments
  use mc_k_limb_dark
  use curand_device
  use cudafor
  implicit none

contains


  attributes(device) subroutine emit_iso(ph)
    implicit none

    type(pac), intent(inout) :: ph

    ph%cost = 2.0_dp*curand_uniform(ph%iseed) - 1.0_dp
    ph%sint = 1.0_dp - ph%cost**2
    if (ph%sint <= 0.0_dp) then
      ph%sint = 0.0_dp
    else
      ph%sint = sqrt(ph%sint)
    end if

    ph%phi = twopi*curand_uniform(ph%iseed)
    ph%cosp = cos(ph%phi)
    ph%sinp = sin(ph%phi)

    ! Cartesian directional vectors
    ph%nxp = ph%sint * ph%cosp
    ph%nyp = ph%sint * ph%sinp
    ph%nzp = ph%cost

    ! Depolarised packet
    ph%fi = 1.0_dp
    ph%fq = 0.0_dp
    ph%fu = 0.0_dp
    ph%fv = 0.0_dp


  end subroutine emit_iso

  !! Performs a source of a packet at a surface with a Lambertian probability distribution
  attributes(device) subroutine emit_iso_surf(ph)
    implicit none

    type(pac), intent(inout) :: ph

    ph%cost = sqrt(curand_uniform(ph%iseed))
    ph%sint = sqrt(1.0_dp - ph%cost**2)

    ph%phi = twopi*curand_uniform(ph%iseed)
    ph%cosp = cos(ph%phi)
    ph%sinp = sin(ph%phi)

    ! Cartesian directional vectors
    ph%nxp = ph%sint * ph%cosp
    ph%nyp = ph%sint * ph%sinp
    ph%nzp = ph%cost

    ! Cell vertical number is one - x and y set to 1
    ph%c(1) = 1
    ph%c(2) = 1
    ph%c(3) = 1

    ! All positions are = 0
    ph%xp = 0.0_dp
    ph%yp = 0.0_dp
    ph%zp = grid_d%r_min

    ! Non polarised packet
    ph%fi = 1.0_dp
    ph%fq = 0.0_dp
    ph%fu = 0.0_dp
    ph%fv = 0.0_dp

    ! Add the positive moment values from the surface if required
    if (do_moments_d .eqv. .True.) then
      call moments_1D(ph)
    end if


  end subroutine emit_iso_surf

  !! Performs a source of a packet at a 1D finite spherical surface with a Lambertian distribution
  attributes(device) subroutine emit_iso_surf_sph(ph)
    implicit none

    type(pac), intent(inout) :: ph

    real(dp) :: cos_lat, sin_lat, lon_p, ran, rp
    real(dp) :: costi, sinti, cospn, phinew

    ! Latitude and longitude of packet
    rp = grid_d%r_min + 1.0e-6_dp
    cos_lat = ph%zp/rp
    sin_lat = sqrt(1.0_dp - cos_lat**2)

    ! Sample a spherical longitude
    lon_p = atan2(ph%yp, ph%xp)

    ! isotropic intensity, sample from n=mu
    ph%cost = sqrt(curand_uniform(ph%iseed))
    ph%sint = sqrt(1.0_dp - ph%cost**2)

    ph%phi = twopi*curand_uniform(ph%iseed)
    ph%cosp = cos(ph%phi)
    ph%sinp = sin(ph%phi)

    !transform to coordinate system of star
    if (abs(cos_lat) < 0.9999999_dp) then
      costi = ph%cost*cos_lat - ph%sint*sin_lat*ph%cosp
      if (abs(costi) > 1.0_dp) then
        print*,  'initp: costi gt 1',costi
        if(costi.gt.1.0_dp) costi=1.0_dp
        if(costi.lt.-1.0_dp) costi=-1.0_dp
      end if
      sinti=sqrt(1.0_dp-costi**2)
      if (sinti.gt.0.0000001_dp)then
        cospn = (ph%cost-costi*cos_lat)/(sinti*sin_lat)
        if (abs(cospn) > 1.0_dp) then
          if (abs(cospn) > 1.01_dp) print*, 'cospn gt 1 in initp',cospn
          if (cospn < 0.0_dp) then
            cospn=-1.0_dp
          else
            cospn=1.0_dp
          end if
        end if
        if (ph%phi < pi) then
          phinew = acos(cospn)
          if (phinew > pi) then
            print*, 'phinew wrong in initp'
            print*,  'phinew,pi',phinew,pi
            phinew = pi
          end if
        else
          phinew = twopi-acos(cospn)
          if (phinew < pi) then
            print*, 'phinew wrong in initp'
            print*, 'phinew,pi',phinew,pi
            phinew=pi
          end if
        end if
      else
        phinew=0.0_dp
      end if
      ph%phi=phinew
    else
      if (cos_lat < 0.0_dp) costi=-ph%cost
    end if
    ph%phi = ph%phi + lon_p

    if(ph%phi > twopi) ph%phi = ph%phi - twopi
    !if(ph%phi < 0.0_dp) print*, 'phi lt 0',ph%phi
    !if(ph%phi > twopi) print*, 'phi gt 2pi',ph%phi
    ph%sinp=sin(ph%phi)
    ph%cosp=cos(ph%phi)

    ph%cost = costi
    ph%sint = sinti

    !x,y,z in units of stellar radius, add a small amount so
    ! numerical precision doesn't place photon beneath surface
    ph%xp = sin_lat*cos(lon_p)*grid_d%r_min + 1e-6_dp
    ph%yp = sin_lat*sin(lon_p)*grid_d%r_min + 1e-6_dp
    ph%zp = cos_lat*grid_d%r_min + 1e-6_dp

    !print*, ph%id, ph%xp, ph%yp, ph%zp

    ! Cartesian directional vectors
    ph%nxp = ph%sint * ph%cosp
    ph%nyp = ph%sint * ph%sinp
    ph%nzp = ph%cost

    ! Cell vertical number is one
    ph%c(1) = 1


    ! Non polarised packet
    ph%fi = 1.0_dp
    ph%fq = 0.0_dp
    ph%fu = 0.0_dp
    ph%fv = 0.0_dp

  end subroutine emit_iso_surf_sph


end module mc_k_emit_iso
