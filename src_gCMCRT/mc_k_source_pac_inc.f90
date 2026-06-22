module mc_k_source_pac_inc
  use mc_precision
  use mc_data_mod
  use mc_class_pac
  use mc_class_grid
  use mc_class_imag
  use mc_k_moments
  use cudafor
  use curand_device
  implicit none

contains

  attributes(device) subroutine source_pac_inc_1D(ph, mu_s)
    implicit none

    type(pac), intent(inout) :: ph
    real(dp), intent(in) :: mu_s

    !! For incident packets in 1D the incident angle is the zenith angle
    !! Down is negative, so minus mu_s
    ph%cost = -mu_s

    ph%sint = 1.0_dp - ph%cost**2
    if (ph%sint <= 0.0_dp)then
      ph%sint = 0.0_dp
    else
      ph%sint = sqrt(ph%sint)
    endif

    ph%phi = twopi * curand_uniform(ph%iseed)
    ph%sinp = sin(ph%phi)
    ph%cosp = cos(ph%phi)

    ph%nxp = ph%sint * ph%cosp
    ph%nyp = ph%sint * ph%sinp
    ph%nzp = ph%cost

    ! Packet z position is at zmax - xp and yp set to 0
    ph%xp = 0.0_dp
    ph%yp = 0.0_dp
    ph%zp = grid_d%z_max

    ! Packet z cell number is uppermost layer index - x and y cell = 1
    ph%c(1) = 1
    ph%c(2) = 1
    ph%c(3) = grid_d%n_lay

    ! Depolarised packet
    ph%fi = 1.0_dp
    ph%fq = 0.0_dp
    ph%fu = 0.0_dp
    ph%fv = 0.0_dp

    ! Add the negative moment values at the top of the atmosphere if required
    if (do_moments_d .eqv. .True.) then
      call moments_1D(ph)
    end if

  end subroutine source_pac_inc_1D


  attributes(device) subroutine source_pac_inc_3D(ph)
    implicit none

    type(pac), intent(inout) :: ph
    real(dp) :: rr2, r_num, ann_theta
    real(dp) :: xp_help, yp_help, rot_angle

    !! For incident packets in 3D the incident angle is from the direction of the star


    !! Use LHS samples via ph index if optioned
    if (LHS_d .eqv. .False.) then

      ann_theta = curand_uniform(ph%iseed) * twopi

      if (do_trans_d .eqv. .True.) then
       ! Sample annulus
       rr2 = sqrt(grid_d%r_min**2 + (grid_d%r_max**2 - grid_d%r_min**2)*curand_uniform(ph%iseed))
       !rr2 = grid_d%r_min/grid_d%r_max + (1.0_dp - grid_d%r_min/grid_d%r_max)*sqrt(curand_uniform(ph%iseed))
      else
        ! Sample uniform sphere
        rr2 = sqrt(grid_d%r_max**2*curand_uniform(ph%iseed))
        ! For lambertian sphere sampling
        !rr2 = sqrt(grid_d%r_min**2*curand_uniform(ph%iseed))
      end if
    else if (LHS_d .eqv. .True.) then

      ann_theta = x_ran_d(ph%id) * twopi

      if (do_trans_d .eqv. .True.) then
       ! Sample annulus
       rr2 = sqrt(grid_d%r_min**2 + (grid_d%r_max**2 - grid_d%r_min**2)*y_ran_d(ph%id))
       !rr2 = grid_d%r_min/grid_d%r_max + (1.0_dp - grid_d%r_min/grid_d%r_max)*sqrt(curand_uniform(ph%iseed))
      else
        ! Sample uniform sphere
        rr2 = sqrt(grid_d%r_max**2*y_ran_d(ph%id))
        ! For lambertian sphere sampling
        !rr2 = sqrt(grid_d%r_min**2*curand_uniform(ph%iseed))
      end if
    end if


    ph%zp = rr2 * cos(ann_theta)
    ph%yp = rr2 * sin(ann_theta)
    ph%xp = sqrt(grid_d%r_max**2 - ph%zp**2 - ph%yp**2) - 1.0e-12_dp

    ! Packet impact parameter
    ph%bp = sqrt(ph%zp**2 + ph%yp**2)

    if (do_LD_d .eqv. .True.) then
      ! We calculate the limb darkening coefficent and return a new packet weight
      call limb_darkening(rr2, ann_theta,ph)
    end if

    ! Rotate illumination for planet at a phase not at 0
    ! rotate according to given viewing angle
    if (im_d%vphi /= 180.0_dp .and. do_trans_d .eqv. .True.) then
      rot_angle = im_d%vphi * pi/180.0_dp - pi
      xp_help = ph%xp*cos(rot_angle) - ph%yp*sin(rot_angle)
      yp_help = ph%xp*sin(rot_angle) + ph%yp*cos(rot_angle)
      ph%xp = xp_help
      ph%yp = yp_help
    end if

    !! NEEDS TO BE CHANGED FOR MULTIPLE SCATTERING MODE - RAY TRACING SHOULD WORK FINE
    ! Incident radiation plane parallel from the x direction, so cost = 0, and towards phi = pi degrees
    ph%cost = 0.0_dp
    ph%sint = sqrt(1.0_dp-ph%cost**2)
    ph%phi = pi
    ph%cosp = -1.0_dp
    ph%sinp = sqrt(1.0_dp-ph%cosp**2)

    ph%nxp = ph%sint * ph%cosp
    ph%nyp = ph%sint * ph%sinp
    ph%nzp = ph%cost

    ! Depolarised packet
    ph%fi = 1.0_dp
    ph%fq = 0.0_dp
    ph%fu = 0.0_dp
    ph%fv = 0.0_dp

  end subroutine source_pac_inc_3D


  attributes(device) subroutine source_pac_inc_3D_transit(ph)
    implicit none

    type(pac), intent(inout) :: ph
    real(dp) :: rr2, ann_theta, b1, b2, s
    real(dp) :: e1x, e1y, e1z, e2x, e2y, e2z
    real(dp) :: kx, ky, kz
    real(dp) :: r2_star, mu, xb, yb, cps, sps, R_top_sq

    ! Out of transit: the planet does not cover any part of the stellar disk
    ! (delta >= R_s + R_p) or the star is on the observer side (zstar >= 0).
    ! No stellar light is intercepted, so kill the packet up front.
    if (trans_state_d == 0) then
      ph%wght = 0.0_dp
      return
    end if

    R_top_sq = grid_d%r2_max * H_d(1)**2     ! (top-of-atmosphere radius)^2 in cm^2

    if (trans_state_d == 2) then
      ! Full transit: the planet disk is entirely on the star, so sample the
      ! FULL projected disk (0 -> r_max) uniform in area (every point is on-star).
      ! The opaque-core / annulus split is done by the calling kernel.
      if (LHS_d .eqv. .False.) then
        ann_theta = curand_uniform(ph%iseed) * twopi
        rr2 = sqrt(grid_d%r_max**2 * curand_uniform(ph%iseed))
      else
        ann_theta = x_ran_d(ph%id) * twopi
        rr2 = sqrt(grid_d%r_max**2 * y_ran_d(ph%id))
      end if
      b1 = rr2 * cos(ann_theta)
      b2 = rr2 * sin(ann_theta)
      ! Distance^2 from the stellar centre (cm) for the limb-darkening weight.
      r2_star = (b1*H_d(1) - xstar_d)**2 + (b2*H_d(1) - ystar_d)**2
    else
      ! Partial transit (ingress/egress): sample uniformly inside the
      ! star/planet overlap lens by bounding-box rejection, in a frame with the
      ! star centre on +x at distance trans_delta_d.  All accepted points are
      ! on-star by construction, so no rejection of off-star packets is needed.
      ! (LHS does not apply on this branch.)
      do
        xb = trans_xlo_d + (trans_xhi_d - trans_xlo_d) * curand_uniform(ph%iseed)
        yb = trans_ymax_d * (2.0_dp*curand_uniform(ph%iseed) - 1.0_dp)
        r2_star = (xb - trans_delta_d)**2 + yb**2
        if ((xb*xb + yb*yb <= R_top_sq) .and. (r2_star <= R_s_sq_d)) exit
      end do
      ! Rotate from the star-aligned frame back to the (e1,e2) sky basis and
      ! normalise to planet radii for the ray construction below.
      cps = cos(trans_psi_d)
      sps = sin(trans_psi_d)
      b1 = (xb*cps - yb*sps) / H_d(1)
      b2 = (xb*sps + yb*cps) / H_d(1)
    end if

    ! Use the same image-plane basis as the peel-off image projection.
    e1x = -im_d%costo * im_d%cospo
    e1y = -im_d%costo * im_d%sinpo
    e1z =  im_d%sinto

    e2x = -im_d%sinpo
    e2y =  im_d%cospo
    e2z =  0.0_dp

    kx = im_d%obsx
    ky = im_d%obsy
    kz = im_d%obsz

    s = sqrt(max(grid_d%r_max**2 - (b1*b1 + b2*b2), 0.0_dp)) - 1.0e-12_dp

    ph%xp = b1*e1x + b2*e2x - s*kx
    ph%yp = b1*e1y + b2*e2y - s*ky
    ph%zp = b1*e1z + b2*e2z - s*kz

    ph%nxp = kx
    ph%nyp = ky
    ph%nzp = kz
    ph%cost = max(-1.0_dp, min(1.0_dp, ph%nzp))
    ph%sint = sqrt(max(1.0_dp - ph%cost**2, 0.0_dp))
    if (ph%sint > 1.0e-300_dp) then
      ph%cosp = ph%nxp / ph%sint
      ph%sinp = ph%nyp / ph%sint
    else
      ph%cosp = 1.0_dp
      ph%sinp = 0.0_dp
    end if
    ph%phi = atan2(ph%sinp, ph%cosp)

    ph%bp = sqrt(b1*b1 + b2*b2)

    if (do_LD_d .eqv. .True.) then
      mu = sqrt(max(1.0_dp - r2_star/R_s_sq_d, 0.0_dp))
      call limb_darkening_mu(mu, ph)
    end if

    ph%fi = 1.0_dp
    ph%fq = 0.0_dp
    ph%fu = 0.0_dp
    ph%fv = 0.0_dp

  end subroutine source_pac_inc_3D_transit


  attributes(device) subroutine limb_darkening_mu(mus, ph)
    implicit none

    type(pac), intent(inout) :: ph
    real(dp), intent(in) :: mus
    real(dp) :: Imus

    select case(ilimb_d)
    case(1)
      Imus = 1.0_dp - LD_c_d(1)*(1.0_dp - mus)
    case(2)
      Imus = 1.0_dp - LD_c_d(1)*(1.0_dp - mus) - LD_c_d(2)*(1.0_dp - mus)**2
    case(3)
      Imus = 1.0_dp  - LD_c_d(1)*(1.0_dp  - mus) - LD_c_d(2)*(1.0_dp  - sqrt(mus))
    case(4)
      Imus = 1.0_dp  - LD_c_d(1)*(1.0_dp  - mus) - LD_c_d(2)*mus*log(mus)
    case(5)
      Imus = 1.0_dp  - LD_c_d(1)*(1.0_dp  - mus) - LD_c_d(2)/(1.0 - exp(mus))
    case(6)
      Imus = 1.0_dp  - LD_c_d(1)*(1.0_dp  - mus) - LD_c_d(2)*(1.0_dp  - mus**(1.5_dp)) &
      & - LD_c_d(3)*(1.0_dp  - mus**2)
    case(7)
      Imus = 1.0_dp  - LD_c_d(1)*(1.0_dp  - sqrt(mus)) - LD_c_d(2)*(1.0_dp  - mus) &
      & - LD_c_d(3)*(1.0_dp  - mus**(1.5_dp)) - LD_c_d(4)*(1.0_dp  - mus**2)
    case(8)
      Imus = 1.0_dp - LD_c_d(1)*(1.0_dp - mus**LD_c_d(2))
    case default
      Imus = 1.0_dp
    end select

    ph%wght = ph%wght * Imus

  end subroutine limb_darkening_mu

  attributes(device) subroutine limb_darkening(rr2, ann_theta, ph)
    implicit none

    type(pac), intent(inout) :: ph
    real(dp), intent(in) :: rr2, ann_theta
    real(dp) :: phase_n, inc_n
    real(dp) :: b, p, a_sm, Rstar, rr2s
    real(dp) :: zp, yp, zs, ys, xs, thetas, phis, mus, Imus

    real(dp) :: zs_cent, ys_cent, xs_cent, theta_cent, phi_cent, mu_cent, Imu_cent

    a_sm = sm_ax_d * Au
    Rstar = Rs_d * Rsun
    rr2s = rr2 * H_d(1)

    phase_n = 90.0_dp - (im_d%vphi - 180.0_dp)

    inc_n = inc_d * pi/180.0_dp
    phase_n = phase_n * pi/180.0_dp

    b = -(a_sm * cos(inc_n))/Rstar
    if (inc_d < 0.0) then
      b = abs(b)
    end if
    p = (a_sm * cos(phase_n))/Rstar

    zp = rr2s * cos(ann_theta)
    yp = rr2s * sin(ann_theta)

    ! zs_cent = b * Rstar
    ! ys_cent = p * Rstar
    ! xs_cent = sqrt(Rstar**2 - zs_cent**2 - ys_cent**2)
    ! theta_cent = acos(zs_cent/Rstar) - pi/2.0_dp
    ! phi_cent = atan2(ys_cent, xs_cent)
    ! mu_cent = cos(theta_cent) * cos(phi_cent)
    ! Imu_cent = 1.0_dp - LD_c_d(1)*(1.0_dp - mu_cent) - LD_c_d(2)*(1.0_dp - mu_cent)**2
    ! theta_cent = theta_cent * 180.0_dp/pi
    ! phi_cent = phi_cent * 180.0_dp/pi
    !
    ! print*,' Central (phi, theta), mu: ', phi_cent, theta_cent, mu_cent

    zs = zp + b * Rstar
    ys = yp  + p * Rstar
    xs = Rstar**2 - zs**2 - ys**2
    if (xs > 0.0_dp) then
      xs = sqrt(xs)
    else
      ph%wght = 0.0_dp
      return
    end if

    ! Find longitude and latitude coordinate on star (radians)
    thetas = acos(zs/Rstar) - pi/2.0_dp
    phis = atan2(ys, xs)

    mus = cos(thetas) * cos(phis)

    select case(ilimb_d)
      ! Laws and references taken from John Southworth's (Keele) website
    case(1)
      ! Linear - Schwarzschild (1906)
      Imus = 1.0_dp - LD_c_d(1)*(1.0_dp - mus)
    case(2)
      ! Quadratic - Kopal (1950)
      Imus = 1.0_dp - LD_c_d(1)*(1.0_dp - mus) - LD_c_d(2)*(1.0_dp - mus)**2
    case(3)
      ! Square root law - Díaz-Cordovés & Giménez (1992)
      Imus = 1.0_dp  - LD_c_d(1)*(1.0_dp  - mus) - LD_c_d(2)*(1.0_dp  - sqrt(mus))
    case(4)
      ! Logarithmic - Klinglesmith & Sobieski (1970)
      Imus = 1.0_dp  - LD_c_d(1)*(1.0_dp  - mus) - LD_c_d(2)*mus*log(mus)
    case(5)
      ! Exponential law - Claret & Hauschildt (2003)
      Imus = 1.0_dp  - LD_c_d(1)*(1.0_dp  - mus) - LD_c_d(2)/(1.0 - exp(mus))
    case(6)
      ! Three paramater - Sing (2009)
      Imus = 1.0_dp  - LD_c_d(1)*(1.0_dp  - mus) - LD_c_d(2)*(1.0_dp  - mus**(1.5_dp)) &
      & - LD_c_d(3)*(1.0_dp  - mus**2)
    case(7)
      ! Four parameter (non-linear) - Claret (2000)
      Imus = 1.0_dp  - LD_c_d(1)*(1.0_dp  - sqrt(mus)) - LD_c_d(2)*(1.0_dp  - mus) &
      & - LD_c_d(3)*(1.0_dp  - mus**(1.5_dp)) - LD_c_d(4)*(1.0_dp  - mus**2)
    case(8)
      ! power2 law - Morello et al. 2017, Hestroffer 1997
      Imus = 1.0_dp - LD_c_d(1)*(1.0_dp - mus**LD_c_d(2))
    case default
      print*, 'Invalid ilimb: ', ilimb_d
      !stop
    end select

    ph%wght = ph%wght * Imus

    !print*, thetas* 180.0/pi, phis* 180.0/pi, Imus


  end subroutine limb_darkening

end module mc_k_source_pac_inc
