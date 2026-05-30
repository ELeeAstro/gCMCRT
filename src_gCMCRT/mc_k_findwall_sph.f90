module mc_k_findwall_sph
  use mc_data_mod
  use mc_class_pac
  use mc_class_grid
  use cudafor
  implicit none

contains

  attributes(device) subroutine findwall_1D_sph(ph,ioffset,dcell)
    implicit none

    type(pac), intent(in) :: ph
    integer, intent(out) :: ioffset
    real(dp), intent(out) :: dcell

    integer :: iface
    real(dp) :: ds

    call radface(r_d(ph%c(3)),r_d(ph%c(3)+1),ph%xp,ph%yp,ph%zp,ph%nxp,ph%nyp,ph%nzp,ds,iface)

    ioffset = -1 + 2*iface
    dcell = ds

  end subroutine findwall_1D_sph

  attributes(device) subroutine findwall_sph(ph,ioffset,dcell)
    implicit none

    type(pac), intent(in) :: ph
    integer, dimension(3) :: cell
    integer, dimension(3), intent(out) :: ioffset
    integer, dimension(3) :: idist
    integer :: i, iface
    real(dp) :: xp,yp,zp,nxp,nyp,nzp
    real(dp), intent(out) :: dcell
    real(dp), dimension(3) :: dist
    real(dp) :: ds, R1, R2, A1, A2, B1, B2, C1, C2, D1, D2, T1, T2
    real(dp) :: eps_face

    cell(1) = ph%c(1)
    cell(2) = ph%c(2)
    cell(3) = ph%c(3)

    R1 = r_d(cell(1))
    R2 = r_d(cell(1)+1)
    A1 = aarr_d(cell(2))
    A2 = aarr_d(cell(2)+1)
    B1 = barr_d(cell(2))
    B2 = barr_d(cell(2)+1)
    C1 = 0.0_dp
    C2 = 0.0_dp
    D1 = 0.0_dp
    D2 = 0.0_dp
    T1 = tan2thet_d(cell(3))
    T2 = tan2thet_d(cell(3)+1)

    xp = ph%xp
    yp = ph%yp
    zp = ph%zp
    nxp = ph%nxp
    nyp = ph%nyp
    nzp = ph%nzp

    call radface(R1,R2,xp,yp,zp,nxp,nyp,nzp,ds,iface)
    dist(1) = ds
    idist(1) = iface
    !print*, ds, iface, '1'

    call phiface(A1,B1,C1,D1,A2,B2,C2,D2,xp,yp,zp,nxp,nyp,nzp,ds,iface)
    dist(2) = ds
    idist(2) = iface
    !print*, ds, iface, '2'

    call thetaface(T1,T2,xp,yp,zp,nxp,nyp,nzp,ds,iface)
    dist(3) = ds
    idist(3) = iface
    !print*, ds, iface, '3'

    ! Find max distance to cell face
    dcell = max(dist(1),dist(2),dist(3))

    ! Minimum distance to cell face
    do i = 1, 3
      if ((dist(i) > 0.0_dp) .and. (dist(i) < dcell)) then
        dcell = dist(i)
      end if
    end do

    eps_face = max(1.0e-12_dp * max(1.0_dp, abs(dcell)), 1.0e-14_dp)

    ioffset(1) = 0
    ioffset(2) = 0
    ioffset(3) = 0
    do i = 1, 3
      if ((dist(i) > 0.0_dp) .and. (abs(dist(i) - dcell) <= eps_face)) then
        ! "shortest distance" offsets are either -1 or 1
        ioffset(i) = -1 + 2 * idist(i)
      end if
    end do

  end subroutine findwall_sph

  attributes(device) subroutine radface(R1,R2,xp,yp,zp,nxp,nyp,nzp,dr,iface)
    implicit none

    integer, intent(out) :: iface
    integer, dimension(4) :: ind
    integer :: i, npos
    real(dp), intent(in) :: R1,R2,xp,yp,zp,nxp,nyp,nzp
    real(dp), intent(out) :: dr
    real(dp), dimension(4) :: root, posroot
    real(dp) :: bb, cc, det

    ! IND array indicates which "face" is associated with which roots
    ind(1) = 0
    ind(2) = 0
    ind(3) = 1
    ind(4) = 1

    bb = 2.0_dp*(xp*nxp+yp*nyp+zp*nzp)
    cc = xp*xp + yp*yp + zp*zp - R1*R1

    det = bb*bb - 4.0_dp*cc

    if (det < 0.0_dp) then ! no roots
      root(1) = -999.0_dp
      root(2) = -999.0_dp
    else
      det = sqrt(det)
      root(1) = (-bb+det)/2.0_dp
      root(2) = (-bb-det)/2.0_dp
    end if

    cc = xp*xp + yp*yp + zp*zp - R2*R2

    det = bb*bb - 4.0_dp*cc

    if (det < 0.0_dp) then ! no roots
      root(3) = -999.0_dp
      root(4) = -999.0_dp
    else
      det = sqrt(det)
      root(3) = (-bb+det)/2.0_dp
      root(4) = (-bb-det)/2.0_dp
    end if

    ! find only positive roots
    npos = 0
    do i = 1, 4
      if (root(i) > 0.0_dp) then
        npos = npos + 1
        posroot(npos) = root(i)
      end if
    end do

    ! all roots < 0
    if (npos == 0) then
      dr = -1.0_dp
    else
      if (npos == 1) then
        dr = posroot(1)
      else if (npos == 2) then
        dr = min(posroot(1),posroot(2))
      else if (npos == 3) then
        dr = min(posroot(1),posroot(2),posroot(3))
      else
        dr = min(posroot(1),posroot(2),posroot(3),posroot(4))
      end if

      do i = 1, 4
        if (root(i) == dr) then
          iface = ind(i)
        end if
      end do

    end if

  end subroutine radface

  attributes(device) subroutine phiface(A1,B1,C1,D1,A2,B2,C2,D2, &
                                        xp,yp,zp,nxp,nyp,nzp,ds,iface)

    implicit none

    integer, intent(out) :: iface
    integer, dimension(2) :: ind
    integer :: i

    real(kind=dp), intent(in) :: xp,yp,zp,nxp,nyp,nzp
    real(kind=dp), intent(in) :: A1,B1,C1,D1,A2,B2,C2,D2
    real(kind=dp), intent(out) :: ds

    real(kind=dp), dimension(2) :: root
    real(kind=dp) :: denom1, denom2
    real(kind=dp) :: numer1, numer2
    real(kind=dp) :: rootmax, rootmin
    real(kind=dp) :: eps1, eps2

    ! IND array indicates which face is associated with which root.
    ind(1) = 0
    ind(2) = 1

    denom1 = A1*nxp + B1*nyp + C1*nzp
    denom2 = A2*nxp + B2*nyp + C2*nzp

    numer1 = A1*xp + B1*yp + C1*zp + D1
    numer2 = A2*xp + B2*yp + C2*zp + D2

    ! Scale-aware near-parallel tests.
    eps1 = 1.0e-14_dp * max(1.0_dp, abs(A1*nxp) + abs(B1*nyp) + abs(C1*nzp))
    eps2 = 1.0e-14_dp * max(1.0_dp, abs(A2*nxp) + abs(B2*nyp) + abs(C2*nzp))

    if (abs(denom1) <= eps1) then
       root(1) = -2.0_dp
    else
       root(1) = -numer1 / denom1
    end if

    if (abs(denom2) <= eps2) then
       root(2) = -2.0_dp
    else
       root(2) = -numer2 / denom2
    end if

    rootmax = max(root(1), root(2))
    rootmin = min(root(1), root(2))

    if (rootmax <= 0.0_dp) then
       ds = -1.0_dp
    else if (rootmin < 0.0_dp) then
       ds = rootmax
    else
       ds = rootmin
    end if

    iface = -1
    do i = 1, 2
       if (root(i) == ds) then
          iface = ind(i)
       end if
    end do

  end subroutine phiface

  attributes(device) subroutine thetaface(tan2th1,tan2th2, &
                                          xp,yp,zp,nxp,nyp,nzp,ds,iface)

    implicit none

    integer, intent(out) :: iface
    integer, dimension(4) :: ind
    integer :: i, npos

    real(kind=dp), intent(in) :: xp,yp,zp,nxp,nyp,nzp
    real(kind=dp), intent(in) :: tan2th1, tan2th2
    real(kind=dp), intent(out) :: ds

    real(kind=dp), dimension(4) :: root, posroot
    real(kind=dp) :: aa, bb, cc, det
    real(kind=dp) :: tan2tha, tan2thb
    real(kind=dp) :: eps_quad, eps_lin

    ! Initialise roots defensively.
    root(1) = -999.0_dp
    root(2) = -999.0_dp
    root(3) = -999.0_dp
    root(4) = -999.0_dp

    if (tan2th2 < tan2th1) then
       tan2tha = tan2th2
       tan2thb = tan2th1

       ind(1) = 1
       ind(2) = 1
       ind(3) = 0
       ind(4) = 0
    else
       tan2tha = tan2th1
       tan2thb = tan2th2

       ind(1) = 0
       ind(2) = 0
       ind(3) = 1
       ind(4) = 1
    end if

    ! First theta face.
    if (tan2tha < 0.0_dp) then
       ! Polar/special face represented as z = 0 plane.
       if (abs(nzp) <= 1.0e-14_dp) then
          root(1) = -999.0_dp
          root(2) = -999.0_dp
       else
          root(1) = -zp / nzp
          root(2) = -999.0_dp
       end if
    else
       aa = nxp*nxp + nyp*nyp - nzp*nzp*tan2tha
       bb = 2.0_dp * (xp*nxp + yp*nyp - zp*nzp*tan2tha)
       cc = xp*xp + yp*yp - zp*zp*tan2tha

       eps_quad = 1.0e-14_dp * max(1.0_dp, abs(nxp*nxp) + abs(nyp*nyp) + &
                                   abs(nzp*nzp*tan2tha))
       eps_lin  = 1.0e-14_dp * max(1.0_dp, abs(bb))

       if (abs(aa) <= eps_quad) then
          ! Degenerate quadratic: solve bb*s + cc = 0.
          if (abs(bb) <= eps_lin) then
             root(1) = -999.0_dp
             root(2) = -999.0_dp
          else
             root(1) = -cc / bb
             root(2) = -999.0_dp
          end if
       else
          det = bb*bb - 4.0_dp*aa*cc

          if (det < 0.0_dp) then
             root(1) = -999.0_dp
             root(2) = -999.0_dp
          else
             det = sqrt(det)
             root(1) = (-bb + det) / (2.0_dp*aa)
             root(2) = (-bb - det) / (2.0_dp*aa)
          end if
       end if
    end if

    ! Second theta face.
    if (tan2thb < 0.0_dp) then
       if (abs(nzp) <= 1.0e-14_dp) then
          root(3) = -999.0_dp
          root(4) = -999.0_dp
       else
          root(3) = -zp / nzp
          root(4) = -999.0_dp
       end if
    else
       aa = nxp*nxp + nyp*nyp - nzp*nzp*tan2thb
       bb = 2.0_dp * (xp*nxp + yp*nyp - zp*nzp*tan2thb)
       cc = xp*xp + yp*yp - zp*zp*tan2thb

       eps_quad = 1.0e-14_dp * max(1.0_dp, abs(nxp*nxp) + abs(nyp*nyp) + &
                                   abs(nzp*nzp*tan2thb))
       eps_lin  = 1.0e-14_dp * max(1.0_dp, abs(bb))

       if (abs(aa) <= eps_quad) then
          if (abs(bb) <= eps_lin) then
             root(3) = -999.0_dp
             root(4) = -999.0_dp
          else
             root(3) = -cc / bb
             root(4) = -999.0_dp
          end if
       else
          det = bb*bb - 4.0_dp*aa*cc

          if (det < 0.0_dp) then
             root(3) = -999.0_dp
             root(4) = -999.0_dp
          else
             det = sqrt(det)
             root(3) = (-bb + det) / (2.0_dp*aa)
             root(4) = (-bb - det) / (2.0_dp*aa)
          end if
       end if
    end if

    ! Select smallest positive root.
    npos = 0

    do i = 1, 4
       if (root(i) > 0.0_dp) then
          npos = npos + 1
          posroot(npos) = root(i)
       end if
    end do

    if (npos == 0) then
       ds = -1.0_dp
       iface = -1
    else if (npos == 1) then
       ds = posroot(1)
    else if (npos == 2) then
       ds = min(posroot(1), posroot(2))
    else if (npos == 3) then
       ds = min(posroot(1), posroot(2), posroot(3))
    else
       ds = min(posroot(1), posroot(2), posroot(3), posroot(4))
    end if

    if (npos > 0) then
       iface = -1
       do i = 1, 4
          if (root(i) == ds) then
             iface = ind(i)
          end if
       end do
    end if

  end subroutine thetaface

end module mc_k_findwall_sph
