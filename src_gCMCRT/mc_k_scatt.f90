!! Module containing device level subroutines (GPU) that scatter packets isotropically
module mc_k_scatt
  use mc_precision
  use mc_class_pac
  use mc_class_grid
  use mc_data_mod
  use mc_k_scatt_mat
  use mc_k_lambertian, only: sample_lambertian_dir
  use curand_device
  use cudafor
  use ieee_arithmetic
  use libm
  implicit none


contains

  attributes(device) subroutine scatt_pac_iso(ph)
    implicit none

    type(pac), intent(inout) :: ph


    ph%cost = 2.0_dp * curand_uniform(ph%iseed) - 1.0_dp
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

  end subroutine scatt_pac_iso

  attributes(device) subroutine depolarise_pac(ph)
    implicit none

    type(pac), intent(inout) :: ph

    ph%fi = 1.0_dp
    ph%fq = 0.0_dp
    ph%fu = 0.0_dp
    ph%fv = 0.0_dp

  end subroutine depolarise_pac

  attributes(device) subroutine set_angles_from_dir(ph)
    implicit none

    type(pac), intent(inout) :: ph
    real(dp) :: norm

    norm = sqrt(ph%nxp*ph%nxp + ph%nyp*ph%nyp + ph%nzp*ph%nzp)

    if (norm > 1.0e-300_dp) then
      ph%nxp = ph%nxp / norm
      ph%nyp = ph%nyp / norm
      ph%nzp = ph%nzp / norm
    else
      ph%nxp = 0.0_dp
      ph%nyp = 0.0_dp
      ph%nzp = 1.0_dp
    end if

    ph%cost = max(-1.0_dp, min(1.0_dp, ph%nzp))
    ph%sint = sqrt(max(1.0_dp - ph%cost*ph%cost, 0.0_dp))

    if (ph%sint > 1.0e-14_dp) then
      ph%cosp = ph%nxp / ph%sint
      ph%sinp = ph%nyp / ph%sint
      ph%phi = atan2(ph%sinp, ph%cosp)
      if (ph%phi < 0.0_dp) ph%phi = ph%phi + twopi
    else
      ph%cosp = 1.0_dp
      ph%sinp = 0.0_dp
      ph%phi = 0.0_dp
    end if

  end subroutine set_angles_from_dir

  attributes(device) subroutine rotate_direction_about_current(ph, bmu)
    implicit none

    type(pac), intent(inout) :: ph
    real(dp), intent(in) :: bmu

    real(dp) :: kx, ky, kz
    real(dp) :: e1x, e1y, e1z
    real(dp) :: e2x, e2y, e2z
    real(dp) :: hx, hy, hz
    real(dp) :: norm
    real(dp) :: mu, sin_theta, phi, cphi, sphi

    kx = ph%nxp
    ky = ph%nyp
    kz = ph%nzp

    norm = sqrt(kx*kx + ky*ky + kz*kz)

    if (norm > 1.0e-300_dp) then
      kx = kx / norm
      ky = ky / norm
      kz = kz / norm
    else
      kx = 0.0_dp
      ky = 0.0_dp
      kz = 1.0_dp
    end if

    if (abs(kz) < 0.9_dp) then
      hx = 0.0_dp
      hy = 0.0_dp
      hz = 1.0_dp
    else
      hx = 1.0_dp
      hy = 0.0_dp
      hz = 0.0_dp
    end if

    e1x = hy*kz - hz*ky
    e1y = hz*kx - hx*kz
    e1z = hx*ky - hy*kx

    norm = sqrt(e1x*e1x + e1y*e1y + e1z*e1z)

    if (norm > 1.0e-300_dp) then
      e1x = e1x / norm
      e1y = e1y / norm
      e1z = e1z / norm
    else
      e1x = 1.0_dp
      e1y = 0.0_dp
      e1z = 0.0_dp
    end if

    e2x = ky*e1z - kz*e1y
    e2y = kz*e1x - kx*e1z
    e2z = kx*e1y - ky*e1x

    norm = sqrt(e2x*e2x + e2y*e2y + e2z*e2z)

    if (norm > 1.0e-300_dp) then
      e2x = e2x / norm
      e2y = e2y / norm
      e2z = e2z / norm
    else
      e2x = 0.0_dp
      e2y = 1.0_dp
      e2z = 0.0_dp
    end if

    mu = max(-1.0_dp, min(1.0_dp, bmu))
    sin_theta = sqrt(max(1.0_dp - mu*mu, 0.0_dp))

    phi = twopi * curand_uniform(ph%iseed)
    cphi = cos(phi)
    sphi = sin(phi)

    ph%nxp = mu*kx + sin_theta*cphi*e1x + sin_theta*sphi*e2x
    ph%nyp = mu*ky + sin_theta*cphi*e1y + sin_theta*sphi*e2y
    ph%nzp = mu*kz + sin_theta*cphi*e1z + sin_theta*sphi*e2z

    call set_angles_from_dir(ph)

  end subroutine rotate_direction_about_current

  attributes(device) subroutine sample_hg_mu(ph, hgg, bmu)
    implicit none

    type(pac), intent(inout) :: ph
    real(dp), intent(in) :: hgg
    real(dp), intent(out) :: bmu

    real(dp) :: g_safe, g2, xi, denom

    g_safe = max(-0.999999999999_dp, min(0.999999999999_dp, hgg))

    if (abs(g_safe) < 1.0e-8_dp) then
      bmu = 2.0_dp * curand_uniform(ph%iseed) - 1.0_dp
    else
      g2 = g_safe*g_safe
      xi = curand_uniform(ph%iseed)
      denom = 1.0_dp - g_safe + 2.0_dp*g_safe*xi

      if (abs(denom) <= 1.0e-300_dp) then
        bmu = sign(1.0_dp, g_safe)
      else
        bmu = ((1.0_dp + g2) - ((1.0_dp - g2)/denom)**2) / (2.0_dp*g_safe)
      end if
    end if

    bmu = max(-1.0_dp, min(1.0_dp, bmu))

  end subroutine sample_hg_mu

  attributes(device) subroutine scatt_pac_2(ph)
    implicit none

    type(pac), intent(inout) :: ph

    real(dp) :: bmu
    real(dp) :: hgg, gf1, gf2, gb1, gb2, alph
    real(dp) :: q, u, xi
    real(dp) :: dg1, dg2, dg4
    real(dp) :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t9
    real(dp) :: rad1, rad2
    real(dp) :: r_here, nx_norm, ny_norm, nz_norm, mu_lam

    select case(ph%iscatt)

    case(1)

      call scatt_pac_iso(ph)
      call depolarise_pac(ph)
      return

    case(2)

      r_here = sqrt(ph%xp*ph%xp + ph%yp*ph%yp + ph%zp*ph%zp)

      if (r_here > 1.0e-300_dp) then
        nx_norm = ph%xp / r_here
        ny_norm = ph%yp / r_here
        nz_norm = ph%zp / r_here
      else
        nx_norm = 0.0_dp
        ny_norm = 0.0_dp
        nz_norm = 1.0_dp
      end if

      call sample_lambertian_dir(nx_norm, ny_norm, nz_norm, ph%iseed, &
                                ph%nxp, ph%nyp, ph%nzp, mu_lam)
      call set_angles_from_dir(ph)
      call depolarise_pac(ph)
      return

    case(3)

      q = 4.0_dp * curand_uniform(ph%iseed) - 2.0_dp
      u = cbrt(-q + sqrt(1.0_dp + q*q))

      if (abs(u) <= 1.0e-300_dp) then
        bmu = 0.0_dp
      else
        bmu = u - 1.0_dp/u
      end if

      bmu = max(-1.0_dp, min(1.0_dp, bmu))
      call rotate_direction_about_current(ph, bmu)
      call depolarise_pac(ph)
      return

    case(4)

      hgg = g_d(ph%c(1), ph%c(2), ph%c(3))
      call sample_hg_mu(ph, hgg, bmu)
      call rotate_direction_about_current(ph, bmu)
      call depolarise_pac(ph)
      return

    case(5)

      hgg = g_d(ph%c(1), ph%c(2), ph%c(3))

      if (hgg >= 0.0_dp) then
        gf1 = hgg
        gb1 = -gf1/2.0_dp
      else
        gb1 = hgg
        gf1 = -gb1/2.0_dp
      end if

      gf2 = gf1*gf1
      gb2 = gb1*gb1
      alph = min(1.0_dp, max(0.0_dp, 1.0_dp - gb2))

      if (curand_uniform(ph%iseed) < alph) then
        call sample_hg_mu(ph, gf1, bmu)
      else
        call sample_hg_mu(ph, gb1, bmu)
      end if

      call rotate_direction_about_current(ph, bmu)
      call depolarise_pac(ph)
      return

    case(6)

      dg1 = Dgg_d(ph%c(1), ph%c(2), ph%c(3))

      if (abs(dg1) < 1.0e-8_dp) then
        bmu = 2.0_dp * curand_uniform(ph%iseed) - 1.0_dp
      else
        dg2 = dg1*dg1
        dg4 = dg2*dg2
        alph = Draine_alp_d

        t0 = alph * (1.0_dp - dg2)
        t1 = alph * (dg4 - 1.0_dp)
        t2 = -3.0_dp * (4.0_dp*(dg4 - dg2) + t1*(1.0_dp + dg2))
        t3 = dg1 * (2.0_dp*curand_uniform(ph%iseed) - 1.0_dp)
        t4 = 3.0_dp*dg2*(1.0_dp + t3) + &
             alph*(2.0_dp + dg2*(1.0_dp + (1.0_dp + 2.0_dp*dg2)*t3))
        t5 = t0*(t1*t2 + t4*t4) + t1**3
        t6 = t0*4.0_dp*(dg4 - dg2)
        rad1 = t5*t5 - t6**3

        if ((abs(t0) <= 1.0e-300_dp) .or. (rad1 < 0.0_dp)) then
          call sample_hg_mu(ph, dg1, bmu)
        else
          t7 = cbrt(t5 + sqrt(rad1))

          if (abs(t7) <= 1.0e-300_dp) then
            call sample_hg_mu(ph, dg1, bmu)
          else
            t8 = 2.0_dp * (t1 + t6/t7 + t7) / t0
            t9 = sqrt(max(6.0_dp*(1.0_dp + dg2) + t8, 0.0_dp))

            if (abs(t9) <= 1.0e-300_dp) then
              call sample_hg_mu(ph, dg1, bmu)
            else
              rad2 = 6.0_dp*(1.0_dp + dg2) - t8 + 8.0_dp*t4/(t0*t9)
              rad2 = max(rad2, 0.0_dp)

              bmu = dg1/2.0_dp + &
                    (1.0_dp/(2.0_dp*dg1) - &
                    1.0_dp/(8.0_dp*dg1) * (sqrt(rad2) - t9)**2)
            end if
          end if
        end if
      end if

      bmu = max(-1.0_dp, min(1.0_dp, bmu))
      call rotate_direction_about_current(ph, bmu)
      call depolarise_pac(ph)
      return

    case default

      call scatt_pac_iso(ph)
      call depolarise_pac(ph)
      return

    end select

  end subroutine scatt_pac_2

  attributes(device) subroutine scatt_pac(ph)
    implicit none

    type(pac), intent(inout) :: ph
    real(dp) :: p1,p2,p3,p4
    real(dp) :: si,sq,su,sv
    real(dp) :: bmu, b, cosb2, costp, sintp, phip, sinbt, sinb2
    real(dp) :: a11,a12,a13,a21,a22,a23,a24,a31,a32,a33,a34,a42,a43,a44
    real(dp) :: ri1,ri3,cosi3,sini3,sini2,bott,cosdph
    real(dp) :: cosi2,sin2i3,sin2i2,cos2i3,cos2i2,sin2,cos2,sin2cos1
    real(dp) :: cos2sin1,cosi1,sini1,sin2i1,cos2i1
    real(dp) :: hgg,g2,alph, q, u
    real(dp) :: f_norm, a_norm, rprob

    real(dp) :: w1, w2, w3, g, h, l, w, g1, gf1, gf2, gb1, gb2

    real(dp) :: G4, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9
    
    real(dp) :: cost_norm, sint_norm, rp

    real(dp) :: r_here
    real(dp) :: nx_norm, ny_norm, nz_norm
    real(dp) :: mu_lam

    ! Normalise polarisation fractions to I stokes parameter
    f_norm = ph%fi
    ph%fi = 1.0_dp
    ph%fq = ph%fq/f_norm
    ph%fu = ph%fu/f_norm
    ph%fv = ph%fv/f_norm

    ! Do rejection methods for each type of scattering:
    ! Isotropic, Rayleigh, HG, TTHG or Mie scattering
    ! Then find correct scattering matrix

    select case(ph%iscatt)
    case(1)

      ! Isotropic scattering

      call scatt_pac_iso(ph)

      ! De-polarised packet
      ph%fi = 1.0_dp
      ph%fq = 0.0_dp
      ph%fu = 0.0_dp
      ph%fv = 0.0_dp

      return

    case(2)

      ! Lambertian surface scattering.
      !
      ! For a spherical lower boundary, the outward surface normal is r / |r|.
      ! This assumes xp, yp, zp are the packet position at the surface.
      r_here = sqrt(ph%xp*ph%xp + ph%yp*ph%yp + ph%zp*ph%zp)

      if (r_here > 1.0e-30_dp) then
        nx_norm = ph%xp / r_here
        ny_norm = ph%yp / r_here
        nz_norm = ph%zp / r_here
      else
        nx_norm = 0.0_dp
        ny_norm = 0.0_dp
        nz_norm = 1.0_dp
      end if

      call sample_lambertian_dir(nx_norm, ny_norm, nz_norm, ph%iseed, &
                                ph%nxp, ph%nyp, ph%nzp, mu_lam)

      ! Reconstruct angular variables from the sampled Cartesian direction.
      ph%cost = max(-1.0_dp, min(1.0_dp, ph%nzp))
      ph%sint = sqrt(max(1.0_dp - ph%cost*ph%cost, 0.0_dp))

      if (ph%sint > 1.0e-14_dp) then
        ph%cosp = ph%nxp / ph%sint
        ph%sinp = ph%nyp / ph%sint
        ph%phi = atan2(ph%sinp, ph%cosp)

        if (ph%phi < 0.0_dp) then
          ph%phi = ph%phi + twopi
        end if
      else
        ph%cosp = 1.0_dp
        ph%sinp = 0.0_dp
        ph%phi = 0.0_dp
      end if

      ! Do not reconstruct nxp/nyp/nzp here.
      ! sample_lambertian_dir already returned the normalised global direction.

      ! De-polarised packet
      ph%fi = 1.0_dp
      ph%fq = 0.0_dp
      ph%fu = 0.0_dp
      ph%fv = 0.0_dp

      return

    case(3)

      ! Rayleigh scattering via direct spherical coordinate sampling
      ! Assumes non-polarised incident packet
      q = 4.0_dp*curand_uniform(ph%iseed) - 2.0_dp
      u = cbrt((-q + sqrt(1.0_dp + q**2)))
      bmu = u - 1.0_dp/u

      call Raymat(p1,p2,p3,p4,bmu,bmu**2)

      !! Rejection method attempt below


      ! Rayleigh scattering via rejection method including polarisation
      ! do
      !   bmu = 1.0_dp - 2.0_dp*curand_uniform(ph%iseed)
      !   if (abs(bmu) > 1.0_dp) then
      !     if (bmu > 1.0_dp) then
      !       bmu = 1.0_dp
      !       cosb2 = 1.0_dp
      !       b = 0.0_dp
      !       sinbt = 0.0_dp
      !     else
      !       bmu = -1.0_dp
      !       cosb2 = 1.0_dp
      !       b = 0.0_dp
      !       sinbt = 0.0_dp
      !     end if
      !   else
      !     cosb2 = bmu**2
      !     b = cosb2 - 1.0_dp
      !     sinbt = sqrt(1.0_dp - cosb2)
      !   end if
      !   call Raymat(p1,p2,p3,p4,bmu,cosb2)
      !   ri1 = twopi * curand_uniform(ph%iseed)
      !   if (ri1 > pi) then
      !     ri3 = twopi - ri1
      !     cosi3 = cos(ri3)
      !     sini3 = sin(ri3)
      !     sin2i3 = 2.0_dp * sini3 * cosi3
      !     cos2i3 = 2.0_dp * cosi3**2 - 1.0_dp
      !     a11 = p1
      !     a12 = p2*cos2i3
      !     a13 = p2*sin2i3
      !   else !(ri1 <= pi)
      !     cosi1 = cos(ri1)
      !     sini1 = sin(ri1)
      !     sin2i1 = 2.0_dp * sini1 * cosi1
      !     cos2i1 = 2.0_dp * cosi1**2 - 1.0_dp
      !     a11 = p1
      !     a12 = p2 * cos2i1
      !     a13 = -p2 * sin2i1
      !  end if
      !
      !   rprob = (a11*ph%fi+a12*ph%fq+a13*ph%fu)
      !   if (rprob > 2.0_dp) then
      !     rprob = 2.0_dp
      !   end if
      !   !print*, ph%id, ph%fi, ph%fq, ph%fu
      !   if (ieee_is_nan(rprob) .eqv. .True.) then
      !     ph%fi = 1.0_dp
      !     ph%fq = 0.0_dp
      !     ph%fu = 0.0_dp
      !     ph%fv = 0.0_dp
      !     cycle
      !   end if
      !   if ((2.0_dp*curand_uniform(ph%iseed)) > rprob) then
      !     cycle
      !   else
      !     exit
      !   end if
      ! end do
      ! !
      ! !
      ! a_norm = rprob

      !call Raymat_depol(p1,p2,p3,p4,bmu,cosb2,sinb2)

    case(4)
        ! HG function
        hgg = g_d(ph%c(1),ph%c(2),ph%c(3))

        ! If near isotropic (g = 0), sample from isotropic phase function
        if (abs(hgg) < 1e-4_dp) then

          call scatt_pac_iso(ph)

          ! De-polarised packet
          ph%fi = 1.0_dp
          ph%fq = 0.0_dp
          ph%fu = 0.0_dp
          ph%fv = 0.0_dp

          return

        else
          ! Sample HG function directly
          g2 = hgg**2
          bmu = ((1.0_dp + g2) - &
           & ((1.0_dp - g2) / (1.0_dp - hgg + 2.0_dp * hgg * curand_uniform(ph%iseed)))**2) &
           & / (2.0_dp*hgg)

          call dustmat_HG(p1,p2,p3,p4,bmu,bmu**2,hgg,g2)
       end if

    case(5)

        ! TTHG function
        !Cahoy et al. (2010)
        !hggb = -hgg/2.0_dp
        !alpha = 1 - hggb**2
        !beta = hggb**2

        hgg = g_d(ph%c(1),ph%c(2),ph%c(3))
        if (hgg >= 0.0_dp) then
          gf1 = hgg
          gf2 = gf1**2
          gb1 = -gf1/2.0_dp
          gb2 = gb1**2
          alph = 1.0_dp - gb2
        else
          gb1 = hgg
          gb2 = gb1**2
          gf1 = -gb1/2.0_dp
          gf2 = gf1**2
          alph = 1.0_dp - gb2
        end if

        if (curand_uniform(ph%iseed) < alph) then
          bmu = ((1.0_dp + gf2) - &
            &((1.0_dp - gf2) / (1.0_dp - gf1 + 2.0_dp * gf1 * curand_uniform(ph%iseed)))**2) &
            & / (2.0_dp*gf1)
        else
          bmu = ((1.0_dp + gb2) - &
            &((1.0_dp - gb2) / (1.0_dp - gb1 + 2.0_dp * gb1 * curand_uniform(ph%iseed)))**2) &
            & / (2.0_dp*gb1)
        end if

        call dustmat_TTHG(p1,p2,p3,p4,bmu,bmu**2,gf1,gf2,gb1,gb2,alph)

      case(6)

        ! Combined HG and Rayleigh function from Draine (2003) -  = Cornette and Shanks 1992 when alpha = 1
        ! Sample from Draine (2003) phase function using the analytical method from
        ! Jendersie & d'Eon (2023) - An Approximate Mie Scattering Function for Fog and Cloud Rendering

        G1 = Dgg_d(ph%c(1),ph%c(2),ph%c(3))
        G2 = Dgg_d(ph%c(1),ph%c(2),ph%c(3))**2
        G4 = Dgg_d(ph%c(1),ph%c(2),ph%c(3))**4

        alph = Draine_alp_d

        t0 = alph - alph*G2
        t1 = alph*G4 - alph
        t2 = -3.0*(4.0*(G4 - G2) + t1*(1.0 + G2))
        t3 = G1*(2.0*curand_uniform(ph%iseed) - 1.0)
        t4 = 3.0*G2*(1.0+t3) + alph*(2.0+G2*(1.0+(1.0+2.0*G2)*t3))
        t5 = t0*(t1*t2+t4**2)+t1**3
        t6 = t0*4.0*(G4 - G2)
        t7 = cbrt(t5 + sqrt(t5**2 - t6**3))
        t8 = 2.0*(t1+t6/t7+t7)/t0
        t9 = sqrt(6.0*(1.0+G2) + t8)

        bmu = G1/2.0 + (1.0/(2.0*G1) - 1.0/(8.0*G1)*(sqrt(6.0*(1.0+G2) -  t8 + 8.0*t4/(t0*t9)) - t9)**2)

        call dustmat_Draine(p1,p2,p3,p4,bmu,bmu**2,G1,G2,alph) ! Use HG function for now

      case default

        print*, "Can't do this yet!", ph%iscatt
        !stop

      end select


    ! work variables
    costp = ph%cost
    sintp = ph%sint
    phip = ph%phi

    !! Begin scattering matrix calculations
    if (abs(bmu) > 1.0_dp) then
      if (bmu > 1.0_dp) then
        bmu = 1.0_dp
        cosb2 = 1.0_dp
        b = 0.0_dp
        sinbt = 0.0_dp
      else
        bmu = -1.0_dp
        cosb2 = 1.0_dp
        b = 0.0_dp
        sinbt = 0.0_dp
      end if
    else
      cosb2 = bmu**2
      b = cosb2 - 1.0_dp
      sinbt = sqrt(1.0_dp - cosb2)
    end if

    ri1 = twopi * curand_uniform(ph%iseed)

    if (ri1 > pi) then
      ri3 = twopi - ri1
      cosi3 = cos(ri3)
      sini3 = sin(ri3)
      sin2i3 = 2.0_dp * sini3 * cosi3
      cos2i3 = 2.0_dp * cosi3**2 - 1.0_dp
      a11 = p1
      a12 = p2*cos2i3
      a13 = p2*sin2i3

      if (bmu == 1.0_dp) then
        return
      else if (bmu == -1.0_dp) then
          sini3=1.0_dp
          sini2=1.0_dp
          cosi3=0.0_dp
          cosi2=0.0_dp
          ph%cost=-costp
          ph%sint=sintp
      else

        ph%cost = costp * bmu + sintp * sinbt * cosi3
        if (abs(ph%cost) < 1.0_dp) then
          ph%sint = abs(sqrt(1.0_dp - ph%cost**2))
          sini2 = sini3 * sintp / ph%sint
          bott = ph%sint * sinbt
          cosi2 = costp / bott - ph%cost * bmu / bott
          if (abs(cosi2) > 1.001_dp) then
            print*, 'cosi2',ph%id, cosi2
          end if
        else
          ph%sint = 0.0_dp
          sini2 = 0.0_dp
          if (ph%cost >= 1.0_dp) then
            cosi2 = -1.0_dp
          else if (ph%cost <= -1.0_dp) then
            cosi2 = 1.0_dp
          end if
        end if
      end if

      cosdph = -cosi2 * cosi3 + sini2 * sini3 * bmu
      if (abs(cosdph) > 1.0_dp) then
        if (abs(cosdph) > 1.001_dp) then
          print*, 'cosdph',ph%id, cosdph
        end if
        if (cosdph > 1.0_dp) then
          cosdph = 1.0_dp
        else
          cosdph = -1.0_dp
        end if
      end if

      ph%phi = phip + acos(cosdph)
      if (ph%phi > twopi) then
        ph%phi = ph%phi - twopi
      else if (ph%phi < 0.0_dp) then
        ph%phi = ph%phi + twopi
      end if

      sin2i2 = 2.0_dp * sini2 * cosi2
      cos2i2 = 2.0_dp * cosi2**2 - 1.0_dp
      sin2 = sin2i2 * sin2i3
      cos2 = cos2i2 * cos2i3
      sin2cos1 = sin2i2 * cos2i3
      cos2sin1 = cos2i2 * sin2i3

      a21 = p2 * cos2i2
      a22 = p1 * cos2 - p3 * sin2
      a23 = p1 * cos2sin1 + p3 * sin2cos1
      a24 = -p4 * sin2i2
      a31 = -p2 * sin2i2
      a32 = -p1 * sin2cos1 - p3 * cos2sin1
      a33 = -p1 * sin2 + p3 * cos2
      a34 = -p4 * cos2i2
      a42 = -p4 * sin2i3
      a43 = p4 * cos2i3
      a44 = p3

    else !(ri1 <= pi)

      cosi1 = cos(ri1)
      sini1 = sin(ri1)
      sin2i1 = 2.0_dp * sini1 * cosi1
      cos2i1 = 2.0_dp * cosi1**2 - 1.0_dp
      a11 = p1
      a12 = p2 * cos2i1
      a13 = -p2 * sin2i1

      if (bmu == 1.0_dp) then
        return
      else if (bmu == -1.0_dp) then
        sini3=1.0_dp
        sini2=1.0_dp
        cosi3=0.0_dp
        cosi2=0.0_dp
        ph%cost=-costp
        ph%sint=sintp
      else

        ph%cost = costp * bmu + sintp * sinbt * cosi1
        if (abs(ph%cost) < 1.0_dp) then
          ph%sint = abs(sqrt(1.0_dp - ph%cost**2))
          sini2 = sini1 * sintp / ph%sint
          bott = ph%sint * sinbt
          cosi2 = costp / bott - ph%cost * bmu / bott
        else
          ph%sint = 0.0_dp
          sini2 = 0.0_dp
          if (ph%cost >= 1.0_dp) then
            cosi2 = -1.0_dp
          else if (ph%cost <= -1.0_dp) then
            cosi2 = 1.0_dp
          end if
        end if
      end if

        cosdph = -cosi1 * cosi2 + sini1 * sini2 * bmu
        if (abs(cosdph) > 1.0_dp) then
          if (abs(cosdph) > 1.001_dp) then
            print*, 'cosdph',ph%id, cosdph
          end if
          if (cosdph > 1.0_dp) then
            cosdph = 1.0_dp
          else
            cosdph = -1.0_dp
          end if
        end if

        ph%phi = phip - acos(cosdph)
        if (ph%phi > twopi) then
          ph%phi = ph%phi - twopi
        else if (ph%phi < 0.0_dp) then
          ph%phi = ph%phi + twopi
        end if

        sin2i2 = 2.0_dp * sini2 * cosi2
        cos2i2 = 2.0_dp * cosi2 * cosi2 - 1.0_dp
        sin2 = sin2i2 * sin2i1
        cos2 = cos2i2 * cos2i1
        sin2cos1 = sin2i2 * cos2i1
        cos2sin1 = cos2i2 * sin2i1

        a21 = p2 * cos2i2
        a22 = p1 * cos2 - p3 * sin2
        a23 = -p1 * cos2sin1 - p3 * sin2cos1
        a24 = p4 * sin2i2
        a31 = p2 * sin2i2
        a32 = p1 * sin2cos1 + p3 * cos2sin1
        a33 = -p1 * sin2 + p3 * cos2
        a34 = -p4 * cos2i2
        a42 = p4 * sin2i1
        a43 = p4 * cos2i1
        a44 = p3

      end if

    if (ph%iscatt >= 4) then
      a_norm = p1
    else
      a_norm = (a11*ph%fi+a12*ph%fq+a13*ph%fu)
    end if

    si = (a11*ph%fi+a12*ph%fq+a13*ph%fu)/a_norm
    sq = (a21*ph%fi+a22*ph%fq+a23*ph%fu+a24*ph%fv)/a_norm
    su = (a31*ph%fi+a32*ph%fq+a33*ph%fu+a34*ph%fv)/a_norm
    sv = (a42*ph%fq+a43*ph%fu+a44*ph%fv)/a_norm

    ph%fi = si * f_norm
    ph%fq = sq * f_norm
    ph%fu = su * f_norm
    ph%fv = sv * f_norm

    if (ieee_is_nan(ph%fi) .eqv. .True.) then
      ph%fi = 1.0_dp
      ph%fq = 0.0_dp
      ph%fu = 0.0_dp
      ph%fv = 0.0_dp
    end if
    if (ieee_is_nan(ph%fq) .eqv. .True.) then
      ph%fi = 1.0_dp
      ph%fq = 0.0_dp
      ph%fu = 0.0_dp
      ph%fv = 0.0_dp
    end if
    if (ieee_is_nan(ph%fu) .eqv. .True.) then
      ph%fi = 1.0_dp
      ph%fq = 0.0_dp
      ph%fu = 0.0_dp
      ph%fv = 0.0_dp
    end if
    if (ieee_is_nan(ph%fv) .eqv. .True.) then
      ph%fi = 1.0_dp
      ph%fq = 0.0_dp
      ph%fu = 0.0_dp
      ph%fv = 0.0_dp
    end if

    ph%cosp = cos(ph%phi)
    ph%sinp = sin(ph%phi)

    ph%nxp = ph%sint * ph%cosp
    ph%nyp = ph%sint * ph%sinp
    ph%nzp = ph%cost

  end subroutine scatt_pac

end module mc_k_scatt
