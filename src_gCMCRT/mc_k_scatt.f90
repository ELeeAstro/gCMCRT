!! Module containing device level subroutines (GPU) that scatter packets isotropically
module mc_k_scatt
  use mc_precision
  use mc_class_pac
  use mc_class_grid
  use mc_data_mod
  use mc_k_scatt_mat
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

      ! De-polarised packet
      ph%fi = 1.0_dp
      ph%fq = 0.0_dp
      ph%fu = 0.0_dp
      ph%fv = 0.0_dp

      return

    case(2)

      ! Lambertian surface scattering

      ! Cell vertical number is 1
      ph%c(1) = 1

      ! Grid vertical distance is at minimum distance
      ph%xp = sqrt(grid_d%r_min**2 - ph%yp**2 - ph%zp**2) + 1.0e-6_dp

      ! Sampled cosine from the normal
      bmu = sqrt(curand_uniform(ph%iseed))
      rp = grid_d%r_min + 1.0e-6_dp
      cost_norm = ph%zp/rp
      sint_norm = sqrt(1.0_dp - cost_norm**2)

      ph%phi = twopi * curand_uniform(ph%iseed)
      ph%sinp = sin(ph%phi)
      ph%cosp = cos(ph%phi)

      ph%cost = bmu*cost_norm - ph%sint*sint_norm*ph%cosp
      ph%sint = sqrt(1.0_dp - ph%cost**2)

      ph%nxp = ph%sint * ph%cosp
      ph%nyp = ph%sint * ph%sinp
      ph%nzp = ph%cost

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
        stop

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
