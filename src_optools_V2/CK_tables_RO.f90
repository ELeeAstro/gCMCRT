module CK_table_RO
  use optools_data_mod
  use optools_aux, only : sort2, locate, linear_log_interp
  implicit none


  private
  public :: RO_CK, RO_CK_2

contains

  ! TODO: Send each tables g-weights and convolve to a new Gw grid

  !! _2 is the new RO version - based on the CHIMERA method (more accurate when less nG)
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

    ! Check if only 1 Ck table
    if (nCK == 1) then
      CK_RO(:) = CK_work(1,:) * VMR_lay(CK_tab(1)%iVMR,z) * N_lay(z)
      return
    end if

    ! Proceed to mix the k-tables using random overlap
    nG2 = nG * nG

    !! Start RO procedure
    ! Initial k-table is 1st table * VMR
    VMR_tot = VMR_lay(CK_tab(1)%iVMR,z)
    CK_RO(:) = CK_work(1,:)
    ! Loop over all other tables
    do s = 2, nCK
      ! Track current cumulative VMR (_cum) and next VMR (_tot)
      VMR_cum = VMR_tot
      VMR_tot = VMR_tot + VMR_lay(CK_tab(s)%iVMR,z)
      ! Skip this species if very low contribution to the band
      if (VMR_lay(CK_tab(s)%iVMR,z)*sum(CK_work(s,:)*Gw(:)) < 1.0e-60_dp) then
        cycle
      end if
      ! Loop nG by nG, perform random k and weight mixing following Amundsen et al. (2017)
      do i = 1, nG
        do j = 1, nG
          k_mix((i-1)*nG+j) = (VMR_cum*CK_RO(i) + VMR_lay(CK_tab(s)%iVMR,z)*CK_work(s,j)) &
            & / VMR_tot
          wt_mix((i-1)*nG+j) = Gw(i) * Gw(j)
        end do
      end do

      ! Sort the mixed k tables and assosiated weights
      call sort2(nG2, k_mix, wt_mix)

      ! Now reconstruct the x (g) coordinate
      ! Find cumulative sum of the mixed weights
      intg(0) = 0.0_dp
      intg(1) = wt_mix(1)
      do g = 2, nG2
        intg(g) = intg(g-1) + wt_mix(g)
      end do
      ! Normalised cumulative sum of weights
      x(:) = intg(:)/maxval(intg) !*2.0_dp - 1.0_dp !(*2 - 1 not needed here, as here weights go 0-1)

      ! Note:, I belive this works in this case, as due to the larger the weight the more
      ! likely the opacity of that x coordinate is to be sampled in a probabilistic sense
      ! (i.e. takes up more range in the x coordinate), so the cumulative weights normalised gives the fraction of the
      ! importance of that g-ordinate to the total opacity distribution

      ! Now interpolate to the origional x grid, this is the mixed k-table.
      do g = 1, nG
        xval = Gx(g)
        call locate(x(:),xval,ix)
        ix1 = ix + 1
        call linear_log_interp(xval, x(ix), x(ix1), k_mix(ix), k_mix(ix1), CK_RO(g))
      end do

    end do

    ! Now scale the randomly overlapped opacities with the total VMR and number
    ! density of the layer
    CK_RO(:) = VMR_tot * N_lay(z) * CK_RO(:)

  end subroutine RO_CK_2

  !! _2 is the older RO version - based on the NEMESIS method
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

    ! Check if only 1 Ck table
    if (nCK == 1) then
      CK_RO(:) = CK_work(1,:) * VMR_lay(CK_tab(1)%iVMR,z) * N_lay(z)
      return
    end if

    ! Proceed to mix the k-tables using random overlap

    ! q value trackers
    q1 = 0.0_dp
    q2 = 0.0_dp

    ! Perform a loop over number of CK tables
    do s = 1, nCK-1

      ! q1 = cumulative VMR, q2 = VMR of next species
      q1 =  q1 + VMR_lay(CK_tab(s)%iVMR,z)
      q2 =  VMR_lay(CK_tab(s+1)%iVMR,z)


      ! Find convolved weight and contributions
      nloop = 0
      do g1 = 1, nG
        do g2 = 1, nG
          nloop = nloop+1
          weight(nloop) = Gw(g1)*Gw(g2)
          contri(nloop) = (CK_work(s,g1)*q1 + CK_work(s+1,g2)*q2)/(q1+q2)
          !print*, 'contri', k_work_1D(j,n,g1,z), k_work_1D(j+1,n,g2,z), contri(nloop)
        end do
      end do

      ! Find the cumulative g-ordinance weight
      g_work(1) = 0.0_dp
      do g = 1, nG
        g_work(g+1) = g_work(g) + Gw(g)
        !print*, g, g_work(g), del_g(j,g)
      end do

      ! Sort the contributions, keeping their weights in the same order.
      call sort2(nloop, contri, weight)

      ! The new culmulative distribution of weights
      g_dist(0) = 0.0_dp
      g_dist(1) = weight(1)
      !print*,g_dist(0)
      !print*,g_dist(1)
      do g = 2, nloop
        g_dist(g) = weight(g) + g_dist(g-1)
        !print*, g, g_dist(g)
      enddo

      ! 0 the work array
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
      ! Replace k work array with ro work array
      CK_work(s+1,:) = CK_RO(:)

    end do ! s loop

    ! last CK_work index is fully overlapped over all species
    ! Multiply by the cumulative VMR of all species * layer number density
    ! Units are now [cm-1]
    CK_RO(:) = CK_work(nCK,:) * N_lay(z) * (q1 + q2)

  end subroutine RO_CK

end module CK_table_RO
