module mc_views
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  implicit none


  public :: set_phase_views_em, set_transit_views, transit_lens_geom

contains


  !! Calculate the emission viewing angles for the peel-off scheme, and set up
  !! the eclipse (secondary occultation) simulator.
  !!
  !! Produces, as module arrays of length n_loop:
  !!   viewphi_n, viewthet_n : viewing longitude / polar angle per phase (deg)
  !!   phi_key               : orbital phase per entry (sort key, 0..1)
  !!   occ_state             : tri-state occultation routing flag
  !!                             0 = fully visible  -> fsum_occ = fsum
  !!                             1 = partial        -> per-packet geometric mask
  !!                             2 = total          -> fsum_occ = 0 (host zeroes it)
  !!   xstar, ystar, zstar   : stellar centre projected into the per-phase
  !!                            image basis; zstar is along the observer axis (cm)
  !!
  !! Conventions (circular, zero-obliquity, synchronous rotation):
  !!   viewphi = modulo(180 - f, 360) : 180 deg (night) at transit, 0 deg (day) at eclipse
  !!   viewthet = inc                 : holds ONLY for zero obliquity; sub-observer lat = 90 - inc
  subroutine set_phase_views_em(n_phase, n_ecl, n_loop)
    implicit none

    integer, intent(in)  :: n_phase, n_ecl
    integer, intent(out) :: n_loop

    integer :: i, j, k, m, t, n_keep, n_side, n_tot, je
    integer, allocatable, dimension(:) :: idx

    real(dp) :: phi_orb, f, R_p, R_s, a, P
    real(dp) :: phi_ecl, inc_r, dphi_out, dphi_in, arg_out, arg_in
    real(dp) :: phi_lo, phi_hi, key_i, delta, R_in, R_out, u, s
    real(dp) :: sinto, costo, sinpo, cospo
    real(dp) :: e1x, e2x, obsx

    real(dp), dimension(n_phase) :: a_phi, a_theta
    real(dp), dimension(n_ecl)   :: a_phi_ecl, a_theta_ecl, phi_ecl_arr
    real(dp), allocatable, dimension(:) :: tphi, tthet, tkey

    logical :: has_tot, use_bands

    !! ------------------------------------------------------------------
    !! Re-entry guard: free any module arrays from a previous call so the
    !! allocates below are always safe.
    !! ------------------------------------------------------------------
    if (allocated(viewphi_n))  deallocate(viewphi_n)
    if (allocated(viewthet_n)) deallocate(viewthet_n)
    if (allocated(phi_key))    deallocate(phi_key)
    if (allocated(occ_state))  deallocate(occ_state)
    if (allocated(xstar))      deallocate(xstar)
    if (allocated(ystar))      deallocate(ystar)
    if (allocated(zstar))      deallocate(zstar)

    !! ------------------------------------------------------------------
    !! Phase section: coarse, evenly spaced phase grid
    !! ------------------------------------------------------------------
    do i = 1, n_phase
      phi_orb    = real(i-1,dp)/real(n_phase,dp)     ! 0 .. <1, evenly spaced
      f          = 360.0_dp * phi_orb                ! true anomaly (deg), circular
      a_phi(i)   = modulo(180.0_dp - f, 360.0_dp)    ! viewing longitude
      a_theta(i) = inc                               ! viewthet; = inc for zero obliquity
    end do

    !! ------------------------------------------------------------------
    !! No-eclipse path (also the brown-dwarf / isolated-rotator case):
    !! plain phase curve, no occultation.
    !! ------------------------------------------------------------------
    if (do_eclipse .eqv. .False.) then
      n_loop = n_phase
      allocate(viewphi_n(n_loop), viewthet_n(n_loop), phi_key(n_loop), occ_state(n_loop))
      viewphi_n(:)  = a_phi(:)
      viewthet_n(:) = a_theta(:)
      do i = 1, n_phase
        phi_key(i) = real(i-1,dp)/real(n_phase,dp)
      end do
      occ_state(:) = 0
      return
    end if

    !! ------------------------------------------------------------------
    !! Eclipse section
    !! ------------------------------------------------------------------
    if (n_ecl < 2) then
      print*, 'set_phase_views_em: do_eclipse on but n_ecl < 2; increase n_ecl'
      stop
    end if

    ! Planet radius = outermost atmosphere level (top-of-atmosphere radius).
    ! NB: confirm H(grid%n_lev) is an ABSOLUTE radius, not height-above-reference.
    R_p = H(grid%n_lev)

    ! Convert system variables to cgs
    R_s      = Rs * Rsun
    R_s_sq_d = R_s**2                ! device-side R_s^2 for the peel-off mask
    a        = sm_ax * Au
    P        = orbital_period * daysec   ! (unused by geometry; kept for timing)

    inc_r   = inc * pi / 180.0_dp
    phi_ecl = 0.5_dp                 ! circular orbit -> secondary eclipse at phase 0.5

    ! --- outer contact (1st / 4th): sky-plane separation = R_s + R_p ---
    arg_out = (R_s + R_p)**2 - (a*cos(inc_r))**2
    if (arg_out <= 0.0_dp) then
      ! Planet never even grazes the stellar disk -> no occultation at all.
      ! Fall back to the plain phase curve.
      n_loop = n_phase
      allocate(viewphi_n(n_loop), viewthet_n(n_loop), phi_key(n_loop), occ_state(n_loop))
      viewphi_n(:)  = a_phi(:)
      viewthet_n(:) = a_theta(:)
      do i = 1, n_phase
        phi_key(i) = real(i-1,dp)/real(n_phase,dp)
      end do
      occ_state(:) = 0
      return
    end if
    dphi_out = asin(min(1.0_dp, sqrt(arg_out)/(a*sin(inc_r)))) / (2.0_dp*pi)

    ! --- inner contact (2nd / 3rd): sky-plane separation = R_s - R_p ---
    arg_in  = (R_s - R_p)**2 - (a*cos(inc_r))**2
    has_tot = (arg_in > 0.0_dp)
    if (has_tot) then
      dphi_in = asin(min(1.0_dp, sqrt(arg_in)/(a*sin(inc_r)))) / (2.0_dp*pi)
    else
      dphi_in = 0.0_dp               ! grazing eclipse: no flat bottom, bands merge
    end if

    ! --- distribute n_ecl points: cluster in the partial (ingress/egress) bands ---
    if (has_tot .and. n_ecl >= 6) then
      use_bands = .true.
      n_side = max(2, (n_ecl - 2)/2)  ! points per ingress / egress band (incl. contacts)
      n_tot  = n_ecl - 2*n_side       ! remainder -> sparse floor across totality
    else
      use_bands = .false.             ! grazing, or too few points to split into bands
      n_side = n_ecl
      n_tot  = 0
    end if

    ! --- build the eclipse phase list (single source of truth for keys AND angles) ---
    je = 0
    if (use_bands) then
      ! ingress: 1st contact -> 2nd contact (inclusive)
      do j = 1, n_side
        je = je + 1
        phi_ecl_arr(je) = (phi_ecl - dphi_out) &
          & + (dphi_out - dphi_in)*real(j-1,dp)/real(n_side-1,dp)
      end do
      ! totality floor: strictly interior, so the contacts aren't duplicated
      do j = 1, n_tot
        je = je + 1
        phi_ecl_arr(je) = (phi_ecl - dphi_in) &
          & + (2.0_dp*dphi_in)*real(j,dp)/real(n_tot+1,dp)
      end do
      ! egress: 3rd contact -> 4th contact (inclusive)
      do j = 1, n_side
        je = je + 1
        phi_ecl_arr(je) = (phi_ecl + dphi_in) &
          & + (dphi_out - dphi_in)*real(j-1,dp)/real(n_side-1,dp)
      end do
    else
      ! single window, cosine-clustered toward both contacts (Chebyshev-like)
      do j = 1, n_ecl
        u = real(j-1,dp)/real(n_ecl-1,dp)
        s = 0.5_dp*(1.0_dp - cos(pi*u))          ! dense near s=0 and s=1
        je = je + 1
        phi_ecl_arr(je) = (phi_ecl - dphi_out) + 2.0_dp*dphi_out*s
      end do
    end if
    ! je == n_ecl here
    if (je /= n_ecl) then
      print*, 'set_phase_views_em: eclipse phase count mismatch', je, n_ecl
      stop
    end if

    ! --- viewing angles at those eclipse phases ---
    do j = 1, n_ecl
      f = 360.0_dp * phi_ecl_arr(j)
      a_phi_ecl(j)   = modulo(180.0_dp - f, 360.0_dp)
      a_theta_ecl(j) = inc            ! = inc for zero obliquity (see header)
    end do

    !! ------------------------------------------------------------------
    !! Merge: drop coarse phases inside the (outer) eclipse window, keep the
    !! rest, append the dense eclipse set, then sort by orbital phase.
    !! ------------------------------------------------------------------
    phi_lo = phi_ecl - dphi_out
    phi_hi = phi_ecl + dphi_out

    ! count coarse phases OUTSIDE the eclipse window
    n_keep = 0
    do i = 1, n_phase
      phi_orb = real(i-1,dp)/real(n_phase,dp)
      if (phi_orb <= phi_lo .or. phi_orb >= phi_hi) n_keep = n_keep + 1
    end do
    n_loop = n_keep + n_ecl

    allocate(viewphi_n(n_loop), viewthet_n(n_loop), phi_key(n_loop))

    ! --- fill: surviving coarse phases, then the dense eclipse set ---
    k = 0
    do i = 1, n_phase
      phi_orb = real(i-1,dp)/real(n_phase,dp)
      if (phi_orb > phi_lo .and. phi_orb < phi_hi) cycle   ! strict: drop interior
      k = k + 1
      phi_key(k)    = phi_orb
      viewphi_n(k)  = a_phi(i)
      viewthet_n(k) = a_theta(i)
    end do
    do j = 1, n_ecl
      k = k + 1
      phi_key(k)    = phi_ecl_arr(j)      ! same phase that set a_phi_ecl(j)
      viewphi_n(k)  = a_phi_ecl(j)
      viewthet_n(k) = a_theta_ecl(j)
    end do
    ! k == n_loop here

    ! --- argsort by phi_key (insertion sort on an index array) ---
    allocate(idx(n_loop))
    do i = 1, n_loop
      idx(i) = i
    end do
    do i = 2, n_loop
      t = idx(i); key_i = phi_key(t); m = i - 1
      do while (m >= 1)
        if (phi_key(idx(m)) <= key_i) exit
        idx(m+1) = idx(m); m = m - 1
      end do
      idx(m+1) = t
    end do

    ! --- gather all three arrays into sorted order ---
    allocate(tphi(n_loop), tthet(n_loop), tkey(n_loop))
    do i = 1, n_loop
      tphi(i)  = viewphi_n(idx(i))
      tthet(i) = viewthet_n(idx(i))
      tkey(i)  = phi_key(idx(i))
    end do
    viewphi_n  = tphi
    viewthet_n = tthet
    phi_key    = tkey
    deallocate(tphi, tthet, tkey, idx)

    !! ------------------------------------------------------------------
    !! Stellar centre projected into the per-phase image basis.
    !! For a synchronously rotating planet the stellar direction is fixed in
    !! the planet frame.  Project that fixed vector into the image basis used
    !! by the peel-off image coordinates before applying the occultation mask.
    !! ------------------------------------------------------------------
    allocate(xstar(n_loop), ystar(n_loop), zstar(n_loop))

    do i = 1, n_loop
      sinto = sin(viewthet_n(i)*pi/180.0_dp)
      costo = cos(viewthet_n(i)*pi/180.0_dp)
      sinpo = sin(viewphi_n(i)*pi/180.0_dp)
      cospo = cos(viewphi_n(i)*pi/180.0_dp)

      e1x = -costo * cospo
      e2x = -sinpo
      obsx = sinto * cospo

      xstar(i) = a * e1x
      ystar(i) = a * e2x
      zstar(i) = a * obsx
    end do

    !! ------------------------------------------------------------------
    !! Per-phase occultation routing (tri-state).
    !! ------------------------------------------------------------------
    allocate(occ_state(n_loop))

    R_in  = R_s - R_p        ! separation below this -> total eclipse
    R_out = R_s + R_p        ! separation above this -> no occultation

    do i = 1, n_loop
      delta = sqrt(xstar(i)**2 + ystar(i)**2)   ! sky-plane star-planet separation

      if (zstar(i) <= 0.0_dp) then
        occ_state(i) = 0        ! star on far side (transit side) -> never mask
      else if (delta >= R_out) then
        occ_state(i) = 0        ! star clear of the planet disk
      else if (delta <= R_in) then
        occ_state(i) = 2        ! planet fully behind the star (total)
      else
        occ_state(i) = 1        ! ingress / egress annulus (partial)
      end if
    end do

  end subroutine set_phase_views_em


  !! Calculate viewing/orbital geometry for a primary-transit light curve.
  !!
  !! Produces module arrays of length n_loop:
  !!   viewphi_n, viewthet_n : observer direction per phase (deg)
  !!   phi_key               : wrapped orbital phase in [0,1)
  !!   trans_state           : 0 = out of transit, 1 = partial, 2 = full
  !!   xstar, ystar, zstar   : stellar centre projected into the per-phase
  !!                            image basis; zstar is along the observer axis (cm)
  subroutine set_transit_views(n_trans, n_ingress, n_loop)
    implicit none

    integer, intent(in)  :: n_trans, n_ingress
    integer, intent(out) :: n_loop

    integer :: i, j, n_side, n_tot, je
    real(dp) :: R_p, R_s, a, inc_r
    real(dp) :: dphi_out, dphi_in, arg_out, arg_in
    real(dp) :: delta, R_in, R_out, u, s, phi_unwrapped
    real(dp) :: sinto, costo, sinpo, cospo
    real(dp) :: e1x, e2x, obsx
    logical :: has_tot, use_bands
    real(dp), allocatable, dimension(:) :: phi_arr

    if (allocated(viewphi_n))   deallocate(viewphi_n)
    if (allocated(viewthet_n))  deallocate(viewthet_n)
    if (allocated(phi_key))     deallocate(phi_key)
    if (allocated(occ_state))   deallocate(occ_state)
    if (allocated(trans_state)) deallocate(trans_state)
    if (allocated(xstar))       deallocate(xstar)
    if (allocated(ystar))       deallocate(ystar)
    if (allocated(zstar))       deallocate(zstar)

    if (n_trans < 2) then
      print*, 'set_transit_views: n_trans < 2; increase n_trans'
      stop
    end if

    R_p = H(grid%n_lev)
    R_s = Rs * Rsun
    R_s_sq_d = R_s**2
    a = sm_ax * Au
    inc_r = inc * pi / 180.0_dp

    arg_out = (R_s + R_p)**2 - (a*cos(inc_r))**2
    if (arg_out <= 0.0_dp) then
      print*, 'set_transit_views: planet never transits stellar disk'
      stop
    end if
    dphi_out = asin(min(1.0_dp, sqrt(arg_out)/(a*sin(inc_r)))) / (2.0_dp*pi)

    arg_in = (R_s - R_p)**2 - (a*cos(inc_r))**2
    has_tot = (arg_in > 0.0_dp)
    if (has_tot) then
      dphi_in = asin(min(1.0_dp, sqrt(arg_in)/(a*sin(inc_r)))) / (2.0_dp*pi)
    else
      dphi_in = 0.0_dp
    end if

    if (has_tot .and. n_ingress >= 2 .and. n_trans >= 2*n_ingress + 1) then
      use_bands = .true.
      n_side = n_ingress
      n_tot = n_trans - 2*n_side
    else
      use_bands = .false.
      n_side = n_trans
      n_tot = 0
    end if

    n_loop = n_trans
    allocate(phi_arr(n_loop))

    je = 0
    if (use_bands) then
      do j = 1, n_side
        je = je + 1
        phi_arr(je) = -dphi_out + (dphi_out - dphi_in)*real(j-1,dp)/real(n_side-1,dp)
      end do
      do j = 1, n_tot
        je = je + 1
        phi_arr(je) = -dphi_in + (2.0_dp*dphi_in)*real(j,dp)/real(n_tot+1,dp)
      end do
      do j = 1, n_side
        je = je + 1
        phi_arr(je) = dphi_in + (dphi_out - dphi_in)*real(j-1,dp)/real(n_side-1,dp)
      end do
    else
      do j = 1, n_trans
        u = real(j-1,dp)/real(n_trans-1,dp)
        s = 0.5_dp*(1.0_dp - cos(pi*u))
        je = je + 1
        phi_arr(je) = -dphi_out + 2.0_dp*dphi_out*s
      end do
    end if

    if (je /= n_loop) then
      print*, 'set_transit_views: phase count mismatch', je, n_loop
      stop
    end if

    allocate(viewphi_n(n_loop), viewthet_n(n_loop), phi_key(n_loop))
    allocate(trans_state(n_loop))
    allocate(xstar(n_loop), ystar(n_loop), zstar(n_loop))

    do i = 1, n_loop
      phi_unwrapped = phi_arr(i)
      phi_key(i) = modulo(phi_unwrapped, 1.0_dp)

      viewphi_n(i) = modulo(180.0_dp - 360.0_dp*phi_unwrapped, 360.0_dp)
      viewthet_n(i) = inc

      ! For a synchronously rotating planet the stellar direction is fixed in
      ! the planet frame.  Project that vector onto the same image-plane basis
      ! used by source_pac_inc_3D_transit, so the finite stellar-disk mask and
      ! packet coordinates are in the same basis.
      sinto = sin(viewthet_n(i)*pi/180.0_dp)
      costo = cos(viewthet_n(i)*pi/180.0_dp)
      sinpo = sin(viewphi_n(i)*pi/180.0_dp)
      cospo = cos(viewphi_n(i)*pi/180.0_dp)

      e1x = -costo * cospo
      e2x = -sinpo
      obsx = sinto * cospo

      xstar(i) = a * e1x
      ystar(i) = a * e2x
      zstar(i) = a * obsx

      delta = sqrt(xstar(i)**2 + ystar(i)**2)
      R_in  = R_s - R_p
      R_out = R_s + R_p

      if (zstar(i) >= 0.0_dp) then
        trans_state(i) = 0
      else if (delta >= R_out) then
        trans_state(i) = 0
      else if (delta <= R_in) then
        trans_state(i) = 2
      else
        trans_state(i) = 1
      end if
    end do

    deallocate(phi_arr)

  end subroutine set_transit_views


  !! Star/planet overlap ("lens") geometry for the transit LC partial-phase
  !! sampler.  Works in a frame rotated so the star centre lies on +x at
  !! distance d = sqrt(xs^2 + ys^2); psi rotates that frame back to the (e1,e2)
  !! sky basis.  Returns the analytic circle-circle overlap area A_lens and a
  !! tight bounding box [x_lo,x_hi] x [-y_max,y_max] of the lens (all in cm).
  subroutine transit_lens_geom(R_p, R_s, xs, ys, A_lens, psi, x_lo, x_hi, y_max)
    implicit none

    real(dp), intent(in)  :: R_p, R_s, xs, ys
    real(dp), intent(out) :: A_lens, psi, x_lo, x_hi, y_max

    real(dp) :: d, x_star, a1, a2, tri

    d   = sqrt(xs**2 + ys**2)
    psi = atan2(ys, xs)

    if (d >= R_p + R_s) then
      ! No overlap (out of transit).
      A_lens = 0.0_dp
      x_lo = 0.0_dp; x_hi = 0.0_dp; y_max = 0.0_dp
      return
    else if (d <= abs(R_p - R_s)) then
      ! One disk fully inside the other (full transit / total block).
      A_lens = pi * min(R_p, R_s)**2
      x_lo   = max(-R_p, d - R_s)
      x_hi   = min( R_p, d + R_s)
      y_max  = min(R_p, R_s)
      return
    end if

    ! Partial lens: planet at origin, star centre at (d,0).
    x_star = (d**2 + R_p**2 - R_s**2) / (2.0_dp*d)
    y_max  = sqrt(max(R_p**2 - x_star**2, 0.0_dp))
    x_lo   = max(-R_p, d - R_s)
    x_hi   = min( R_p, d + R_s)

    ! Analytic circle-circle intersection area.
    a1  = R_p**2 * acos(min(1.0_dp, max(-1.0_dp, (d**2 + R_p**2 - R_s**2)/(2.0_dp*d*R_p))))
    a2  = R_s**2 * acos(min(1.0_dp, max(-1.0_dp, (d**2 + R_s**2 - R_p**2)/(2.0_dp*d*R_s))))
    tri = 0.5_dp * sqrt(max(0.0_dp, &
      & (-d + R_p + R_s)*(d + R_p - R_s)*(d - R_p + R_s)*(d + R_p + R_s)))
    A_lens = a1 + a2 - tri

  end subroutine transit_lens_geom

end module mc_views
