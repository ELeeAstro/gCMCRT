module CK_table_RO
  use optools_data_mod
  use optools_aux, only : sort2, locate, linear_log_interp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none


  private
  public :: RO_CK, RO_CK_2, RO_CK_RORR

contains

  ! TODO: Send each table's g weights and convolve them onto a new Gw grid.

  !! Random overlap with resorting and conservative weighted-bin rebinning.
  !! The working distribution contains absolute VMR-weighted cross sections,
  !! sum_s(q_s*k_s), so skipped species cannot affect a later normalisation.
  subroutine RO_CK_RORR(z,nG,Gw,CK_work,CK_RO)
    implicit none

    integer, intent(in) :: z,nG
    real(kind=dp), dimension(nCK,nG), intent(in) :: CK_work
    real(kind=dp), dimension(nG), intent(in) :: Gw
    real(kind=dp), dimension(nG), intent(out) :: CK_RO

    real(kind=dp), parameter :: VMR_skip = 1.0e-30_dp
    real(kind=dp), parameter :: rebin_eps = 1.0e-8_dp

    integer :: i, j, g, m, s, nG2, n_active
    real(kind=dp) :: q, weight_sum
    real(kind=dp) :: target_lo, target_hi, source_lo, source_hi
    real(kind=dp) :: overlap, bin_width, bin_sum, boundary_tol
    real(kind=dp), dimension(nG) :: Gw_norm
    real(kind=dp), dimension(nG*nG) :: k_mix, wt_mix
    real(kind=dp), dimension(0:nG) :: target_cum
    real(kind=dp), dimension(0:nG*nG) :: source_cum

    if (nG < 1) then
      print*, 'ERROR - RORR requires at least one g ordinate - STOPPING'
      print*, 'Layer, nG: ', z, nG
      stop
    end if

    if (any(.not. ieee_is_finite(Gw)) .or. any(Gw <= 0.0_dp)) then
      print*, 'ERROR - RORR g weights must be finite and positive - STOPPING'
      print*, 'Layer, min/max weight: ', z, minval(Gw), maxval(Gw)
      stop
    end if

    weight_sum = sum(Gw)
    if (.not. ieee_is_finite(weight_sum) .or. weight_sum <= 0.0_dp) then
      print*, 'ERROR - Invalid total RORR g weight - STOPPING'
      print*, 'Layer, weight sum: ', z, weight_sum
      stop
    end if

    if (.not. ieee_is_finite(N_lay(z)) .or. N_lay(z) < 0.0_dp) then
      print*, 'ERROR - Invalid layer number density in RORR - STOPPING'
      print*, 'Layer, number density: ', z, N_lay(z)
      stop
    end if

    ! RORR operates on probability weights, irrespective of whether the input
    ! quadrature weights sum to one or two.
    Gw_norm(:) = Gw(:)/weight_sum
    target_cum(0) = 0.0_dp
    do g = 1, nG
      target_cum(g) = target_cum(g-1) + Gw_norm(g)
    end do
    target_cum(nG) = 1.0_dp

    nG2 = nG*nG
    n_active = 0
    CK_RO(:) = 0.0_dp

    do s = 1, nCK
      q = VMR_lay(CK_tab(s)%iVMR,z)

      if (.not. ieee_is_finite(q) .or. q < 0.0_dp) then
        print*, 'ERROR - CK VMR must be finite and non-negative - STOPPING'
        print*, 'Layer, species, VMR: ', z, CK_tab(s)%sp, q
        stop
      end if

      ! Apply the requested abundance-only cut-off consistently, including to
      ! the first or only species in a mixture.
      if (q < VMR_skip) cycle

      if (any(.not. ieee_is_finite(CK_work(s,:))) .or. any(CK_work(s,:) < 0.0_dp)) then
        print*, 'ERROR - Invalid CK opacity supplied to RORR - STOPPING'
        print*, 'Layer, species: ', z, CK_tab(s)%sp
        print*, 'Minimum/maximum opacity: ', minval(CK_work(s,:)), maxval(CK_work(s,:))
        stop
      end if

      if (n_active == 0) then
        CK_RO(:) = q*CK_work(s,:)
        n_active = 1
        cycle
      end if

      ! Random-overlap convolution of the current absolute mixture with the
      ! next VMR-weighted species distribution.
      do i = 1, nG
        do j = 1, nG
          m = (i-1)*nG + j
          k_mix(m) = CK_RO(i) + q*CK_work(s,j)
          wt_mix(m) = Gw_norm(i)*Gw_norm(j)
        end do
      end do

      call sort2(nG2,k_mix,wt_mix)

      source_cum(0) = 0.0_dp
      do m = 1, nG2
        source_cum(m) = source_cum(m-1) + wt_mix(m)
      end do
      source_cum(nG2) = 1.0_dp

      ! Conservatively reblock the sorted distribution into bins whose widths
      ! are the original quadrature weights. As in SOCRATES, cumulative source
      ! boundaries within a relative tolerance of 1e-8 above a target boundary
      ! are treated as coincident, avoiding spurious microscopic bin splits.
      ! Other source terms crossing a target boundary contribute to both
      ! adjacent bins in exact proportion to their cumulative-weight overlap.
      m = 1
      do g = 1, nG
        target_lo = target_cum(g-1)
        target_hi = target_cum(g)
        bin_width = target_hi - target_lo
        boundary_tol = rebin_eps*target_hi
        bin_sum = 0.0_dp

        do while (m <= nG2)
          source_lo = source_cum(m-1)
          source_hi = source_cum(m)

          ! Snap a source boundary that exceeds the target boundary only by
          ! the SOCRATES tolerance. Moving the shared source boundary preserves
          ! the total probability and avoids creating a negligible split.
          if (source_hi >= target_hi .and. &
              & source_hi - target_hi <= boundary_tol .and. &
              & target_hi > source_lo) then
            if (m == nG2) then
              source_cum(m) = target_hi
              source_hi = target_hi
            else if (target_hi < source_cum(m+1)) then
              source_cum(m) = target_hi
              source_hi = target_hi
            end if
          end if

          if (source_hi <= target_lo) then
            m = m + 1
            cycle
          end if
          if (source_lo >= target_hi) exit

          overlap = min(target_hi,source_hi) - max(target_lo,source_lo)
          if (overlap > 0.0_dp) bin_sum = bin_sum + overlap*k_mix(m)

          ! Keep m unchanged when this source interval crosses the target
          ! boundary; its unused fraction belongs to the next target bin.
          if (source_hi >= target_hi) exit
          m = m + 1
        end do

        CK_RO(g) = bin_sum/bin_width
      end do

      n_active = n_active + 1
    end do

    if (n_active == 0) then
      CK_RO(:) = 0.0_dp
      return
    end if

    ! Convert the VMR-weighted cross-section distribution to cm^-1.
    CK_RO(:) = N_lay(z)*CK_RO(:)

  end subroutine RO_CK_RORR

  !! RO_CK_2 is based on the CHIMERA method and is retained for comparison.
  subroutine RO_CK_2(z,nG,Gw,Gx,CK_work,CK_RO)
    implicit none

    integer, intent(in) :: z,nG
    real(kind=dp), dimension(nCK,nG), intent(in) :: CK_work
    real(kind=dp), dimension(nG), intent(in) :: Gw, Gx
    real(kind=dp), dimension(nG), intent(out) :: CK_RO

    integer :: i, j, g, s, nG2
    real(kind=dp) :: VMR_tot, VMR_cum
    real(kind=dp), dimension(nG*nG) :: k_mix
    real(kind=dp), dimension(nG*nG) :: wt_mix
    real(kind=dp), dimension(0:nG*nG) :: intg, x

    integer :: ix, ix1
    real(kind=dp) :: xval, x0, x1, y0, y1, yval

    ! Handle the case with only one CK table.
    if (nCK == 1) then
      CK_RO(:) = CK_work(1,:) * VMR_lay(CK_tab(1)%iVMR,z) * N_lay(z)
      return
    end if

    ! Proceed to mix the k-tables using random overlap
    nG2 = nG * nG

    !! Start RO procedure
    ! Initialise the mixture with the first k table and its VMR.
    VMR_tot = VMR_lay(CK_tab(1)%iVMR,z)
    CK_RO(:) = CK_work(1,:)
    ! Loop over the remaining tables.
    do s = 2, nCK
      ! A zero-abundance species has no contribution to the mixture.
      if (VMR_lay(CK_tab(s)%iVMR,z) == 0.0_dp) then
        cycle
      end if

      ! Track current cumulative VMR (_cum) and next VMR (_tot)
      VMR_cum = VMR_tot
      VMR_tot = VMR_tot + VMR_lay(CK_tab(s)%iVMR,z)
      ! Construct the nG-by-nG random-overlap opacity and weight combinations
      ! following Amundsen et al. (2017).
      do i = 1, nG
        do j = 1, nG
          k_mix((i-1)*nG+j) = (VMR_cum*CK_RO(i) + VMR_lay(CK_tab(s)%iVMR,z)*CK_work(s,j)) &
            & / VMR_tot
          wt_mix((i-1)*nG+j) = Gw(i) * Gw(j)
        end do
      end do

      ! Sort the mixed k coefficients and their associated weights.
      call sort2(nG2, k_mix, wt_mix)

      ! Reconstruct the cumulative x (g) coordinate.
      ! Find the cumulative sum of the mixed weights.
      intg(0) = 0.0_dp
      intg(1) = wt_mix(1)
      do g = 2, nG2
        intg(g) = intg(g-1) + wt_mix(g)
      end do
      ! Normalise the cumulative weights. Mapping to [-1,1] is unnecessary
      ! because the input g coordinate spans [0,1].
      x(:) = intg(:)/maxval(intg)

      ! A larger weight occupies a wider interval in cumulative probability,
      ! making its associated opacity more likely to be sampled.

      ! Interpolate the mixed k table onto the original x grid.
      do g = 1, nG
        xval = Gx(g)
        if (xval <= x(1)) then
          CK_RO(g) = k_mix(1)
        else if (xval >= x(nG2)) then
          CK_RO(g) = k_mix(nG2)
        else
          call locate(x(1:nG2),xval,ix)
          ix1 = ix + 1
          call linear_log_interp(xval, x(ix), x(ix1), k_mix(ix), k_mix(ix1), CK_RO(g))
        end if
      end do

    end do

    ! Scale the randomly overlapped opacities by the total VMR and the layer
    ! number density.
    CK_RO(:) = VMR_tot * N_lay(z) * CK_RO(:)

  end subroutine RO_CK_2

  !! This older RO version is based on the NEMESIS method.
  subroutine RO_CK(z,nG,Gw,CK_work,CK_RO)
    implicit none

    integer, intent(in) :: z,nG
    integer :: s
    real(kind=dp), dimension(nCK,nG), intent(inout) :: CK_work
    real(kind=dp), dimension(nG) :: Gw
    real(kind=dp), dimension(nG), intent(out) :: CK_RO

    integer :: g, g1, g2, nloop, ig
    real(kind=dp) :: q1, q2, sumr, frac
    real(kind=dp) :: g_work(nG*nG+1), g_dist(0:nG*nG), weight(nG*nG), contri(nG*nG)

    ! Handle the case with only one CK table.
    if (nCK == 1) then
      CK_RO(:) = CK_work(1,:) * VMR_lay(CK_tab(1)%iVMR,z) * N_lay(z)
      return
    end if

    ! Mix the k tables using random overlap.

    ! VMR trackers.
    q1 = 0.0_dp
    q2 = 0.0_dp

    ! Loop over the CK tables.
    do s = 1, nCK-1

      ! q1 = cumulative VMR, q2 = VMR of next species
      q1 =  q1 + VMR_lay(CK_tab(s)%iVMR,z)
      q2 =  VMR_lay(CK_tab(s+1)%iVMR,z)


      ! Construct the convolved weights and opacity contributions.
      nloop = 0
      do g1 = 1, nG
        do g2 = 1, nG
          nloop = nloop+1
          weight(nloop) = Gw(g1)*Gw(g2)
          contri(nloop) = (CK_work(s,g1)*q1 + CK_work(s+1,g2)*q2)/(q1+q2)
          !print*, 'contri', k_work_1D(j,n,g1,z), k_work_1D(j+1,n,g2,z), contri(nloop)
        end do
      end do

      ! Construct the target cumulative g-weight boundaries.
      g_work(1) = 0.0_dp
      do g = 1, nG
        g_work(g+1) = g_work(g) + Gw(g)
        !print*, g, g_work(g), del_g(j,g)
      end do

      ! Sort the contributions, keeping their weights in the same order.
      call sort2(nloop, contri, weight)

      ! Construct the new cumulative weight distribution.
      g_dist(0) = 0.0_dp
      g_dist(1) = weight(1)
      !print*,g_dist(0)
      !print*,g_dist(1)
      do g = 2, nloop
        g_dist(g) = weight(g) + g_dist(g-1)
        !print*, g, g_dist(g)
      enddo

      ! Initialise the work array to zero.
      CK_RO(:) = 0.0_dp

      ig = 1
      sumr = 0.0
      do g = 1, nloop
        !print*,'C',g,ig,g_dist(g),g_work(ig),g_work(ig+1)
        if ((g_dist(g) < g_work(ig+1)) .and. (ig <= ng)) then
          CK_RO(ig) = CK_RO(ig) + contri(g) * weight(g)
          sumr = sumr + weight(g)
          !print*,'A',g,ig,k_work_1D_ro(n,ig,z), sum, contri(g) * weight(g)
        else
          frac = (g_work(ig+1)-g_dist(g-1)) / (g_dist(g)-g_dist(g-1))
          CK_RO(ig) = CK_RO(ig) + frac * contri(g)*weight(g)
          sumr = sumr + frac * weight(g)
          CK_RO(ig) = CK_RO(ig) / sumr
          ig = ig + 1
          if (ig <= ng)then
            sumr = (1.0_dp-frac) * weight(g)
            CK_RO(ig) = CK_RO(ig) + (1.0_dp-frac)*contri(g)*weight(g)
            !print*,'B',g,ig,k_work_1D_ro(n,ig,z),sum, frac, contri(g) * weight(g)
          endif
        endif
      end do
      if (ig == ng) then
        CK_RO(ig) = CK_RO(ig)/sumr
      endif
      ! Replace the k-table work array with the random-overlap result.
      CK_work(s+1,:) = CK_RO(:)

    end do ! s loop

    ! The final CK_work index contains the overlap of all species.
    ! Multiply it by the cumulative VMR and layer number density.
    ! Units are now cm^-1.
    CK_RO(:) = CK_work(nCK,:) * N_lay(z) * (q1 + q2)

  end subroutine RO_CK

end module CK_table_RO
