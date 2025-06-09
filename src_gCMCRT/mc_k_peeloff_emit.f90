module mc_k_peeloff_emit
  use mc_precision
  use mc_data_mod
  use mc_class_pac
  use mc_class_grid
  use mc_class_imag
  use mc_k_raytrace
  use ieee_arithmetic
  use cudafor
  implicit none

contains

  attributes(device) subroutine peeloff_emit(ph)
    implicit none

    type(pac), intent(in) :: ph

    type(pac) :: ray
    integer :: xl, yl, istat
    real(dp) :: wfac, peel_fac
    real(dp) :: phot, photq, photu, ximage, yimage

    !! We copy all data from the packet to a 'virtual' ray to trace
    !! leaving the origonal packet untouched
    ray = ph

    ! weighting factor for isotropic emission
    wfac = 1.0_dp/fourpi
    !wfac = 1.0_dp/pi

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
      !print *,'peeloff_emit: tflag', ray%p_flag
      istat = atomicadd(im_d%fail_pemit, 1)
      return
    endif

    ! Zero or negative tau_peel
    if ((ray%tau < 0.0_dp) .or. (ieee_is_nan(ray%tau) .eqv. .True.)) then
      !print *,'peeloff_emit: tau_peel negative/NaN',ray%tau
      istat = atomicadd(im_d%fail_pemit, 1)
      return
    end if

    ! Peeloff weighting
    peel_fac = wfac * exp(-ray%tau) * ph%wght

    if (ieee_is_nan(peel_fac) .eqv. .True.) then
      !print *,'peeloff_emit: peel_fac NaN',peel_fac, wfac, exp(-ray%tau), ph%wght
      istat = atomicadd(im_d%fail_pemit, 1)
      return
    end if

    phot = peel_fac * ray%fi
    photq = peel_fac * ray%fq
    photu = peel_fac * ray%fu

    if (ieee_is_nan(phot) .eqv. .True.) then
      !print*, 'peeloff_emit: NaN phot'
      istat = atomicadd(im_d%fail_pemit, 1)
      return
    endif

    !! Add energy to total counters
    !istat = atomicadd(p_noise(ph%bin_idx,na), 1)
    istat = atomicadd(im_d%fsum, phot)
    istat = atomicadd(im_d%qsum, photq)
    istat = atomicadd(im_d%usum, photu)

    if (do_images_d .eqv. .True.) then

      !! Bin the photon into the image according to its position and
      ! direction of travel.
      ximage = im_d%rimage + ph%zp * im_d%sinto &
        & - ph%yp * im_d%costo * im_d%sinpo - ph%xp * im_d%costo * im_d%cospo
      yimage = im_d%rimage + ph%yp * im_d%cospo - ph%xp * im_d%sinpo
      xl = int(im_d%x_pix * ximage / (2.0_dp * im_d%rimage)) + 1
      yl = int(im_d%y_pix * yimage / (2.0_dp * im_d%rimage)) + 1

      if ((xl > im_d%x_pix) .or. (xl < 1)) then
        print*, 'peeloff_emit: xl out of bounds',xl
        print*, ph%xp,ph%yp,ph%zp,ray%tau
        return
      endif
      if ((yl > im_d%y_pix) .or. (yl < 1)) then
        print*, 'peeloff_emit: yl out of bounds',yl
        print*,  ph%xp,ph%yp,ph%zp,ray%tau
        return
      endif

      !! Add weighted peeloff energy to images
      !print*, 'p', f_d(xl,yl), xl, yl, phot
      istat = atomicadd(f_d(xl,yl), phot)
      istat = atomicadd(q_d(xl,yl), photq)
      istat = atomicadd(u_d(xl,yl), photu)
    end if

    if (do_cf_d .eqv. .True.) then
      istat = atomicadd(cf_d(ph%c(1),ph%c(2),ph%c(3)), phot)
    end if

  end subroutine peeloff_emit

end module mc_k_peeloff_emit
