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

    !! For incident packets in 3D the incident angle is from the direction of the star

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

    ph%zp = rr2 * cos(ann_theta)
    ph%yp = rr2 * sin(ann_theta)
    ph%xp = sqrt(grid_d%r_max**2 - ph%zp**2 - ph%yp**2) - 1.0e-12_dp

    ! Packet impact parameter
    ph%bp = sqrt(ph%zp**2 + ph%yp**2)

    if (do_LD_d .eqv. .True.) then
      ! We calculate the limb darkening coefficent and return a new packet weight
      call limb_darkening(rr2, ann_theta,ph)
    end if


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

    if (phase_d <= 0.5_dp) then
      phase_n = 90.0_dp - phase_d*360.0_dp
    else
      phase_n = -90.0_dp - (360.0_dp - phase_d*360.0_dp)
    end if

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
      ! Four parameter - Claret (2000)
      Imus = 1.0_dp  - LD_c_d(1)*(1.0_dp  - sqrt(mus)) - LD_c_d(2)*(1.0_dp  - mus) &
      & - LD_c_d(3)*(1.0_dp  - mus**(1.5_dp)) - LD_c_d(4)*(1.0_dp  - mus**2)
    case default
      print*, 'Invalid ilimb: ', ilimb_d
      stop
    end select

    ph%wght = ph%wght * Imus

    !print*, thetas* 180.0/pi, phis* 180.0/pi, Imus


  end subroutine limb_darkening

end module mc_k_source_pac_inc
