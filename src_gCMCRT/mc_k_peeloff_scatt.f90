module mc_k_peeloff_scatt
  use mc_precision
  use mc_data_mod
  use mc_class_pac
  use mc_class_grid
  use mc_class_imag
  use mc_k_raytrace
  use mc_k_scatt_mat
  use ieee_arithmetic
  use cudafor
  implicit none

contains

  attributes(device) subroutine peeloff_scatt(ph)
    implicit none


    type(pac), intent(in) :: ph

    type(pac) :: ray
    integer :: xl, yl, istat
    real(dp) :: rstat
    real(dp) :: wfac, peel_fac, contri
    real(dp) :: phot, photq, photu, ximage, yimage

    !! We copy all data from the packet to a 'virtual' ray to trace
    !! leaving the origonal packet untouched
    ray = ph

    wfac = 0.0_dp

    !! Find the weighting for wfac for this scattering type
    call scat_peel_test(ray, wfac)

    ! Find tau from position to observation direction
    if (ray%geo == 1) then
      ! Do the 3D carestian grid ray tracing
      call raytrace_cart_3D(ray)
    else if (ray%geo == 2) then
      ! Do the 3D spherical grid ray tracing
      call raytrace_sph_3D(ray)
    end if

    ! Invalid tau path
    if (ray%p_flag /= 0) then
      !print *,'peeloff_scatt: tflag', ray%id, ray%p_flag
      istat = atomicadd(im_d%fail_pscat, 1)
      return
    endif

    ! Zero or negative tau_peel
    if ((ray%tau < 0.0_dp) .or. (ieee_is_nan(ray%tau) .eqv. .True.)) then
      !print *,'peeloff_scatt: tau_peel negative/NaN',ray%tau
      istat = atomicadd(im_d%fail_pscat, 1)
      return
    end if

    ! Peeloff weighting
    peel_fac = wfac * exp(-ray%tau) * ph%wght

    if (ieee_is_nan(peel_fac) .eqv. .True.) then
      !print *,'peeloff_scatt: peel_fac NaN',peel_fac, wfac, exp(-ray%tau), ph%wght
      istat = atomicadd(im_d%fail_pscat, 1)
      return
    end if

    if (do_trans_d .eqv. .True.) then
      ray%bp = sqrt(ray%zp**2 + ray%yp**2)
      contri = peel_fac * ray%bp * H_d(1)
      rstat = atomicadd(T_trans_d,-contri)
    end if

    phot = peel_fac * ray%fi
    photq = peel_fac * ray%fq
    photu = peel_fac * ray%fu

    if ((ieee_is_nan(phot) .or. ieee_is_nan(photq) .or. ieee_is_nan(photu)) .eqv. .True.) then
      !print*, 'peeloff_scatt: NaN phot'
      istat = atomicadd(im_d%fail_pscat, 1)
      return
    endif

    !! Add energy to total counters
    !istat = atomicadd(p_noise(ph%bin_idx,na), 1)
    rstat = atomicadd(im_d%fsum, phot)
    rstat = atomicadd(im_d%qsum, photq)
    rstat = atomicadd(im_d%usum, photu)

    !! Add weighted peeloff energy to images
    if (do_images_d .eqv. .True.) then
      !! Bin the photon into the image according to its position and
      ! direction of travel.
      ximage = im_d%rimage + ph%zp * im_d%sinto &
        & - ph%yp * im_d%costo * im_d%sinpo - ph%xp * im_d%costo * im_d%cospo
      yimage = im_d%rimage + ph%yp * im_d%cospo - ph%xp * im_d%sinpo
      xl = int(im_d%x_pix * ximage / (2.0_dp * im_d%rimage)) + 1
      yl = int(im_d%y_pix * yimage / (2.0_dp * im_d%rimage)) + 1

      if ((xl > im_d%x_pix) .or. (xl < 1)) then
        print*, 'peeloff_scatt: xl out of bounds',xl
        print*, ph%xp,ph%yp,ph%zp,ray%tau
        return
      endif
      if ((yl > im_d%y_pix) .or. (yl < 1)) then
        print*, 'peeloff_scatt: yl out of bounds',yl
        print*,  ph%xp,ph%yp,ph%zp,ph%tau
        return
      endif
      !print*, 'p', f_d(xl,yl), xl, yl, phot
      rstat = atomicadd(f_d(xl,yl), phot)
      rstat = atomicadd(q_d(xl,yl), photq)
      rstat = atomicadd(u_d(xl,yl), photu)
    end if

  end subroutine peeloff_scatt

  attributes(device) subroutine scat_peel_test(ray, wfac)
    implicit none

    type(pac), intent(inout) :: ray
    real(dp), intent(out) :: wfac
  end subroutine scat_peel_test


  attributes(device) subroutine scat_peel(ray, wfac)
    implicit none

    type(pac), intent(inout) :: ray
    real(dp), intent(out) :: wfac

    real(dp) :: p1,p2,p3,p4
    real(dp) :: si,sq,su,sv
    real(dp) :: calpha, cosmu, sinmu, cosmu2, sinmu2, nsinmu2
    real(dp) :: a11,a12,a13,a21,a22,a23,a24,a31,a32,a33,a34,a42,a43,a44
    real(dp) :: ri1,ri3,cosi3,sini3,sini2,bott
    real(dp) :: cosi2,sin2i3,sin2i2,cos2i3,cos2i2,sin2,cos2,sin2cos1
    real(dp) :: cos2sin1,cosi1,sini1,sin2i1,cos2i1
    real(dp) :: hgg,g2,hgga,hggb,g2a,g2b,alph,beta
    real(dp) :: f_norm, a_norm
    real(dp) :: rp, cost_norm, sint_norm, phi_norm, cosp_norm, sinp_norm
    real(dp) :: norm_nxp, norm_nyp, norm_nzp, cos_norm


    real(dp) :: w1, w2, w3, g, h, l, w, g1, gf1, gf2, gb1, gb2

    ! Normalise polarisation fractions to I stokes parameter
    f_norm = ray%fi
    ray%fi = 1.0_dp
    ray%fq = ray%fq/f_norm
    ray%fu = ray%fu/f_norm
    ray%fv = ray%fv/f_norm

    ! calculate calpha=cos(alpha), where alpha is angle between incident
    ! and outgoing (i.e., observed) photon direction
    calpha = ray%nxp * im_d%obsx + ray%nyp * im_d%obsy  + ray%nzp * im_d%obsz
    ! change of variable: mu = alpha

    ! work variables
    cosmu = calpha
    !! Start scattering matrix calculations
    ! Limit cosmu to 1 and -1
      if (cosmu >= 1.0_dp) then
        cosmu = 1.0_dp
        cosmu2 = 1.0_dp
        sinmu = 0.0_dp
        sinmu2 = 0.0_dp
        nsinmu2 = 0.0_dp
      else if (cosmu <= -1.0_dp) then
        cosmu = -1.0_dp
        cosmu2 = 1.0_dp
        sinmu = 0.0_dp
        sinmu2 = 0.0_dp
        nsinmu2 = 0.0_dp
      else
        cosmu2 = cosmu**2
        sinmu = sqrt(1.0_dp - cosmu2)
        sinmu2 = sinmu**2
        nsinmu2 = -sinmu2
      end if

    ! Find weight photon for:
    ! Isotropic, Rayleigh, Thompson (electron), HG, TTHG or Mie scattering
    ! Then use correct scattering matrix for each type of scattering

    select case(ray%iscatt)

    case(1)

      wfac = 1.0_dp / fourpi
      call isomat(p1,p2,p3,p4)

      ! De-polarised packet
      ray%fi = 1.0_dp
      ray%fq = 0.0_dp
      ray%fu = 0.0_dp
      ray%fv = 0.0_dp

      return

    case(2)

      ! Cell vertical number is 1
      ray%c(1) = 1
      ! Grid vertical distance is at minimum distance
      !ray%xp = grid_d%r_min
      ray%xp = sqrt(grid_d%r_min**2 - ray%yp**2 - ray%zp**2) + 1.0e-6_dp

      rp = grid_d%r_min + 1.0e-6_dp
      ! get phi coordinate in the sphere
      cost_norm = ray%zp/rp
      sint_norm = sqrt(1.0_dp - cost_norm**2)

      phi_norm = atan2(ray%yp,ray%xp)
      cosp_norm = cos(phi_norm)
      sinp_norm = sin(phi_norm)

      ! Cartesian directional vectors
      norm_nxp = sint_norm  * cosp_norm
      norm_nyp = sint_norm  * sinp_norm
      norm_nzp = cost_norm

      !! Now we need the angle between the normal direction and the detector
      cos_norm = norm_nxp * im_d%obsx + norm_nyp * im_d%obsy  + norm_nzp * im_d%obsz
      if (cos_norm >= 0.0_dp) then
        wfac = cos_norm / pi
      else
        wfac = 0.0_dp
      end if

      ! De-polarised packet
      ray%fi = 1.0_dp
      ray%fq = 0.0_dp
      ray%fu = 0.0_dp
      ray%fv = 0.0_dp

      return

    case(3)

      wfac = (0.75_dp * (1.0_dp + cosmu2)) / fourpi

      call Raymat(p1,p2,p3,p4,cosmu,cosmu2)
      !call Raymat_depol(p1,p2,p3,p4,bmu,cosb2,sinb2)

    case(4)

        hgg = g_d(ray%c(1),ray%c(2),ray%c(3))
        g2 = hgg**2
        wfac = ((1.0_dp - g2)/(1.0_dp + g2 - 2.0_dp*hgg*cosmu)**(1.5_dp)) / fourpi

        call dustmat_HG(p1,p2,p3,p4,cosmu,cosmu2,hgg,g2)

    case(5)

        !Cahoy et al. (2010)
        !hggb = -hgg/2.0_dp
        !alpha = 1 - hggb**2
        !beta = hggb**2

        hgg = g_d(ray%c(1),ray%c(2),ray%c(3))
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

        wfac = ((alph * ((1.0_dp - gf2)/(1.0_dp + gf2 - 2.0_dp*gf1*cosmu)**1.5_dp)) &
          & + ((1.0_dp - alph) * ((1.0_dp - gb2)/(1.0_dp + gb2 - 2.0_dp*gb1*cosmu)**1.5_dp))) &
          & / fourpi

        call dustmat_TTHG(p1,p2,p3,p4,cosmu,cosmu2,gf1,gf2,gb1,gb2,alph)

      case(6)

        G1 = Dgg_d(ray%c(1),ray%c(2),ray%c(3))
        G2 = Dgg_d(ray%c(1),ray%c(2),ray%c(3))**2

        wfac = (((1.0_dp - G2)/(1.0_dp + Draine_alp_d*(1.0_dp + 2.0_dp*G2)/3.0_dp)) &
        & * ((1.0_dp + Draine_alp_d*cosmu2)/(1.0_dp + G2 - 2.0_dp*G1*cosmu)**1.5_dp)) &
        & / fourpi

        alph = Draine_alp_d

        call dustmat_Draine(p1,p2,p3,p4,cosmu,cosmu2,G1,G2,alph)

      case default

        print*, "Can't do this yet!", ray%iscatt
        !stop

      end select

    ! Find i1 angle in radian
    ! Spherical triginometry to find cos(i1) and sin(i1)
    cosi1 = (im_d%costo - ray%cost * cosmu) / (ray%sint * sinmu)
    sini1 = sin(ray%phi - im_d%phio - pi) * im_d%sinto / sinmu

    ! Use tan(i1) = sin(i1)/cos(i1) identity,
    ! add pi to place between 0 and 2 pi
    ri1 = atan2(sini1,cosi1) + pi

    ! If i1 is > pi, then we are in wrong quadrant
    ! adjust by 2 pi
    if (ri1 > pi) then
      ri3 = twopi - ri1
      cosi3 = cos(ri3)
      sini3 = sin(ri3)
      ! Use trig identities
      sin2i3 = 2.0_dp * sini3 * cosi3
      cos2i3 = 2.0_dp * cosi3 * cosi3 - 1.0_dp

      ! Form matrix elements of R(theta) * L(-i1)
      a11 = p1
      a12 = p2 * cos2i3
      a13 = p2 * sin2i3

      ! if "direct" (extreamly rare) forward/backward scattering return
      if (cosmu == 1.0_dp) then
          return
      else
        if (cosmu == -1.0_dp) then
          ray%fu = -ray%fu
          return
        end if
      end if

      ! if observation angle is within bounds, proceed to calculate i2
      ! else, catch special cases for directly scattering into observation direction
      ! (extreamly rare)
      if (abs(im_d%costo) < 1.0_dp) then
        sini2 = sini3 * ray%sint / im_d%sinto
        bott = im_d%sinto * sinmu
        cosi2 = ray%cost / bott - im_d%costo * cosmu / bott
      else
        sini2 = 0.0_dp
        if (im_d%costo >= 1.0_dp) then
          cosi2 = -1.0_dp
        end if
        if (im_d%costo <= -1.0_dp) then
          cosi2 = 1.0_dp
        end if
      end if

      ! trig identities for i2
      sin2i2 = 2.0_dp * sini2 * cosi2
      cos2i2 = 2.0_dp * cosi2 * cosi2 - 1.0_dp
      sin2 = sin2i2 * sin2i3
      cos2 = cos2i2 * cos2i3
      sin2cos1 = sin2i2 * cos2i3
      cos2sin1 = cos2i2 * sin2i3

      ! Matrix elements
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

    !! else if(ri1.le.pi) then
    else
      ! Correct quadrant
      cosi1 = cos(ri1)
      sini1 = sin(ri1)
      ! Use trig identities
      sin2i1 = 2.0_dp * sini1 * cosi1
      cos2i1 = 2.0_dp * cosi1 * cosi1 - 1.0_dp
      ! Matrix elements
      a11 = p1
      a12 = p2 * cos2i1
      a13 = -p2 * sin2i1

      ! if "direct" (extreamly rare) forward/backward scattering return
      if (cosmu == 1.0_dp) then
        return
      else
        if (cosmu == -1.0_dp) then
          ray%fu = -ray%fu
          return
        end if
      end if

      ! if observation angle is within bounds, proceed to calculate i2
      ! else, catch special cases for directly scattering into observation direction
      ! (extreamly rare)
      if (abs(im_d%costo) < 1.0_dp) then
        sini2 = sini1 * ray%sint / im_d%sinto
        bott = im_d%sinto * sinmu
        cosi2 = ray%cost / bott - im_d%costo * cosmu / bott
      else
        sini2 = 0.0_dp
        if (im_d%costo >= 1.0_dp) then
          cosi2 = -1.0_dp
        end if
        if (im_d%costo <= -1.0_dp) then
          cosi2 = 1.0_dp
        end if
      end if

      ! trig identities for i2
      sin2i2 = 2.0_dp * sini2 * cosi2
      cos2i2 = 2.0_dp * cosi2 * cosi2 - 1.0_dp
      sin2 = sin2i2 * sin2i1
      cos2 = cos2i2 * cos2i1
      sin2cos1 = sin2i2 * cos2i1
      cos2sin1 = cos2i2 * sin2i1

      ! Matrix elements
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

    ! Normalisation factor is the new I stokes parameter equation
    a_norm = (a11 * ray%fi + a12 * ray%fq + a13 * ray%fu)
    !a_norm = a11

    ! New stokes vectors - normalised to I
    si = 1.0_dp!(a11*ray%fi+a12*ray%fq+a13*ray%fu)/a_norm
    sq = (a21*ray%fi + a22*ray%fq + a23*ray%fu + a24*ray%fv) / a_norm
    su = (a31*ray%fi + a32*ray%fq + a33*ray%fu + a34*ray%fv) / a_norm
    sv = (a42*ray%fq + a43 * ray%fu + a44 * ray%fv) / a_norm

    !if (abs(sq) > abs(si)) then
    !  print*,'scattp: Q > I',ray%id,ray%n_scat
    !  print*, si, sq, a_norm
    !  print*, ray%fi,ray%fq,ray%fu
    !  print*,a11,a12,a13
    !  print*,a21,a22,a23
    !  print*,calpha,ri1
    !  print*,im_d%costo,im_d%sinto,im_d%phio
    !  print*, '--'
    !endif

    ! New stokes vector is new Stokes matrix * old Stokes matrix
    ray%fi = si * f_norm
    ray%fq = sq * f_norm
    ray%fu = su * f_norm
    ray%fv = sv * f_norm

  end subroutine scat_peel

end module mc_k_peeloff_scatt
