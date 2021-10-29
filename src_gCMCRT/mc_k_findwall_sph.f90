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

    do i = 1, 3
      ioffset(i) = 0
      if (dist(i) == dcell) then
        !"shortest distance" offsets are either -1 or 1
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

  attributes(device) subroutine phiface(A1,B1,C1,D1,A2,B2,C2,D2,xp,yp,zp,nxp,nyp,nzp,ds,iface)
    implicit none

    integer, intent(out) :: iface
    integer, dimension(2) :: ind
    integer :: i
    real(kind=dp), intent(in) ::  xp,yp,zp,nxp,nyp,nzp
    real(kind=dp), intent(in) :: A1,B1,C1,D1,A2,B2,C2,D2
    real(kind=dp), intent(out) :: ds
    real(kind=dp), dimension(2) :: root
    real(kind=dp) ::  denom1,denom2,rootmax,rootmin

    ! IND arry, indicates which "face" is associated with which roots
    ind(1) = 0
    ind(2) = 1

    denom1 = (A1*nxp + B1*nyp + C1*nzp)
    if (denom1 == 0.0_dp) then
      root(1) = -2.0_dp
    else
      root(1) = -(A1*xp + B1*yp +C1*zp + D1) / denom1
    end if

    denom2 = (A2*nxp + B2*nyp + C2*nzp)
    if (denom2 == 0.0_dp) then
       root(2) = -2.0_dp
    else
       root(2) = -(A2*xp + B2*yp +C2*zp +D2) / denom2
    end if

    rootmax = max(root(1),root(2))
    rootmin = min(root(1),root(2))


    if (rootmax <= 0.0_dp) then        !both roots .le. 0
       ds = -1.0_dp
    else if (rootmin < 0.0_dp) then   !one root .lt. 0
       ds = rootmax
    else                            !both roots ge 0
       ds = rootmin
    end if

    !find which index "t" originated
    iface = -1
    do i = 1, 2
       if (root(i) == ds) then
         iface=ind(i)
       end if
    end do

  end subroutine phiface

  attributes(device) subroutine thetaface(tan2th1,tan2th2,xp,yp,zp,nxp,nyp,nzp,ds,iface)
    implicit none

    integer, intent(out) :: iface
    integer, dimension(4) :: ind
    integer :: i, npos
    real(kind=dp), intent(in) :: xp,yp,zp,nxp,nyp,nzp
    real(kind=dp), intent(in) :: tan2th1, tan2th2
    real(kind=dp), intent(out) :: ds
    real(kind=dp), dimension(4) :: root, posroot
    real(kind=dp) :: aa,bb,cc,tan2tha,tan2thb,det

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
    endif

    if (tan2tha < 0.0_dp) then
       call phiface(0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,&
       & xp,yp,zp,nxp,nyp,nzp,ds,iface)
       root(1)=ds
       root(2)=-999.0_dp
    else
       aa = nxp*nxp + nyp*nyp - nzp*nzp*tan2tha
       bb = 2.0_dp*(xp*nxp + yp*nyp - zp*nzp*tan2tha)
       cc = xp*xp + yp*yp - zp*zp * tan2tha
       det = bb*bb - 4.0_dp*aa*cc
       if (det < 0.0_dp) then
         root(1) = -999.0_dp
         root(2) = -999.0_dp
       else
         det = sqrt(det)
         root(1) = (-bb+det)/(2.0_dp*aa)
         root(2) = (-bb-det)/(2.0_dp*aa)
       endif
    endif

    aa = nxp*nxp + nyp*nyp - nzp*nzp*tan2thb
    bb = 2.0_dp*(xp*nxp+yp*nyp-zp*nzp*tan2thb)
    cc = xp*xp + yp*yp - zp*zp*tan2thb
    det = bb*bb - 4.0_dp*aa*cc
    if (det < 0.0_dp) then
      root(3) = -999.0_dp
      root(4) = -999.0_dp
    else
      det = sqrt(det)
      root(3) = (-bb+det)/(2.0_dp*aa)
      root(4) = (-bb-det)/(2.0_dp*aa)
    endif

    npos = 0
    do i = 1, 4
       if (root(i) > 0.0_dp) then
         npos = npos + 1
         posroot(npos) = root(i)
       endif
    end do

    if (npos == 0) then
      ds = -1.0_dp
      iface = -1
    else
      if (npos == 1) then
        ds=posroot(1)
      else if (npos == 2) then
        ds=min(posroot(1),posroot(2))
      else if (npos == 3) then
        ds=min(posroot(1),posroot(2),posroot(3))
      else
        ds=min(posroot(1),posroot(2),posroot(3),posroot(4))
      end if

      do i = 1, 4
         if (root(i) == ds) then
            iface = ind(i)
          end if
      end do

    end if

  end subroutine thetaface

end module mc_k_findwall_sph
