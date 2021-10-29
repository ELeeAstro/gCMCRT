module mc_k_scatt_mat
  use mc_precision
  use mc_data_mod
  use cudafor
  implicit none

contains

  attributes(device) subroutine isomat(p1,p2,p3,p4)
    implicit none

    real(dp), intent(out) :: p1, p2, p3, p4

    ! Isotropic scattering matrix elements
    p1 = 1.0_dp
    p2 = 0.0_dp
    p3 = 1.0_dp
    p4 = 0.0_dp

  end subroutine isomat

  attributes(device) subroutine Raymat(p1,p2,p3,p4,cosmu,cosmu2)
    implicit none

    real(dp), intent(in) :: cosmu, cosmu2
    real(dp), intent(out) :: p1, p2, p3, p4

    ! Non-depolarising (isotropic particle) Rayleigh scattering
    p1 = 1.0_dp + cosmu2
    p2 = (-1.0_dp + cosmu2)
    p3 = 2.0_dp * cosmu
    p4 = 0.0_dp

  end subroutine Raymat

  attributes(device) subroutine Raymat_depol(p1,p2,p3,p4,cosmu,cosmu2,sinmu2)
    implicit none

     real(dp), intent(in) :: cosmu, cosmu2, sinmu2
     real(dp), intent(out) :: p1, p2, p3, p4
     real(dp) :: dpol

     ! Assume air depol factor for now
     dpol = 0.0279_dp

     ! Depolarising (anisotropic particle) Rayleigh scattering
     p1 = 1.0_dp + cosmu2
     p2 = -sinmu2
     p3 = 2.0_dp * cosmu
     p4 = 0.0_dp

   end subroutine Raymat_depol

 !*****************************************************************************
 attributes(device) subroutine dustmat_HG(p1,p2,p3,p4,cosmu,cosmu2,hgg,g2)
   implicit none

 !                             a.d. code october 25, 1989
 !                             revised baw apr 10, 1990
 !       f90 update Mar 23, 2016 - E. Lee
 !     **********************************************************************
 !
 !     this program calculates the elements of the phase matrix for a
 !     simple representation of the mrn dust mixture using the algorithms
 !     for the ultraviolet region due to richard l. white ap.j. 229, 954,
 !     1979.
 !
 !     ***********************************************************************
 !
 !         bmu = cos(angle) of scattering (i.e. angle between incident
 !             photon and scattered photon)
 !         g = mean value of cosine of scattering angle (henyey-greenstein)
 !         pl = peak linear polarization
 !         pc = peak value of linear to circular conversion
 !         sc = asymmetry of the circular polarization.
 !         p1 = intensity phase function
 !         p2 = polarization function
 !         p3 = skew polarization
 !         p4 = circular polarization
 !
 !         the scattering matrix for (i,q,u,v) is of the form
 !
 !              p1    p2    0    0
 !              p2    p1    0    0
 !              0     0     p3  -p4
 !              0     0     p4   p3
 !
 !     **********************************************************************

 ! Careful, I think this routine may rely on real and dp floating point differences.
    real(dp), intent(in) ::  hgg, g2
    real(dp), intent(inout) :: cosmu
    real(dp), intent(out) :: p1, p2, p3, p4
    real(dp) :: rphi, f, f2, c, cosmu2

    if (abs(cosmu) > 1.0_dp) then
      print *, 'in HG dustmat, cosmu > 1', cosmu
      if (cosmu > 1.0_dp) then
        cosmu = 1.0_dp
        cosmu2 = 1.0_dp
      else
        cosmu = -1.0_dp
        cosmu2 = 1.0_dp
      end if
    end if

    ! HG function
    p1 = (1.0_dp - g2) / (1.0_dp + g2 - 2.0_dp * hgg * cosmu)**1.5_dp

    p2 = p1 * (-pl_d * (1.0_dp - cosmu2) / (1.0_dp + cosmu2))
    p3 = p1 * (2.0_dp * cosmu) / (1.0_dp + cosmu2)

    ! angle in degrees!
    rphi = acos(cosmu) * 180.0_dp / pi
    f = 3.13_dp * rphi * exp(-7.0_dp * rphi / 180.0_dp)
    ! now convert to radians
    f2 = (rphi + sc_d * f) * pi / 180.0_dp
    ! fphi= (1+3.13*sc*exp(-7.0*phi/pi))*phi
    c = (cos(f2))**2
    p4 = -pc_d * p1 * (1.0_dp - c) / (1.0_dp + c)

  end subroutine dustmat_HG

  !*****************************************************************************
  attributes(device) subroutine dustmat_TTHG(p1,p2,p3,p4,cosmu,cosmu2,g1,g2,gb1,gb2,alph)
    implicit none

     real(dp), intent(in) ::  g1, g2, gb1, gb2, alph
     real(dp), intent(inout) :: cosmu
     real(dp), intent(out) :: p1, p2, p3, p4
     real(dp) :: rphi, f, f2, c, cosmu2

     if (abs(cosmu) > 1.0_dp) then
       print *, 'in TTHG dustmat, cosmu > 1', cosmu
       if (cosmu > 1.0_dp) then
         cosmu = 1.0_dp
         cosmu2 = 1.0_dp
       else
         cosmu = -1.0_dp
         cosmu2 = 1.0_dp
       end if
     end if

     !! TTHG function
     p1 = (alph * ((1.0_dp - g2)/(1.0_dp + g2 - 2.0_dp*g1*cosmu)**1.5_dp)) &
       & + ((1.0_dp - alph) * ((1.0_dp - gb2)/(1.0_dp + gb2 - 2.0_dp*gb1*cosmu)**1.5_dp))

     p2 = p1 * (-pl_d * (1.0_dp - cosmu2) / (1.0_dp + cosmu2))
     p3 = p1 * (2.0_dp * cosmu) / (1.0_dp + cosmu2)

     ! angle in degrees!
     rphi = acos(cosmu) * 180.0_dp / pi
     f = 3.13_dp * rphi * exp(-7.0_dp * rphi / 180.0_dp)
     ! now convert to radians
     f2 = (rphi + sc_d * f) * pi / 180.0_dp
     ! fphi= (1+3.13*sc*exp(-7.0*phi/pi))*phi
     c = (cos(f2))**2
     p4 = -pc_d * p1 * (1.0_dp - c) / (1.0_dp + c)

   end subroutine dustmat_TTHG

   !*****************************************************************************
   attributes(device) subroutine dustmat_Draine(p1,p2,p3,p4,cosmu,cosmu2,G1,G2,alph)
     implicit none

      real(dp), intent(in) ::  G1, G2, alph
      real(dp), intent(inout) :: cosmu
      real(dp), intent(out) :: p1, p2, p3, p4
      real(dp) :: rphi, f, f2, c, cosmu2

      if (abs(cosmu) > 1.0_dp) then
        print *, 'in Draine dustmat, cosmu > 1', cosmu
        if (cosmu > 1.0_dp) then
          cosmu = 1.0_dp
          cosmu2 = 1.0_dp
        else
          cosmu = -1.0_dp
          cosmu2 = 1.0_dp
        end if
      end if

      ! Draine (2003) function
      p1 = ((1.0_dp - G2)/(1.0_dp + alph*(1.0_dp + 2.0_dp*G2)/3.0_dp)) &
      & * ((1.0_dp + alph*cosmu2)/(1.0_dp + G2 - 2.0_dp*G1*cosmu)**1.5_dp)

      p2 = p1 * (-pl_d * (1.0_dp - cosmu2) / (1.0_dp + cosmu2))
      p3 = p1 * (2.0_dp * cosmu) / (1.0_dp + cosmu2)

      ! angle in degrees!
      rphi = acos(cosmu) * 180.0_dp / pi
      f = 3.13_dp * rphi * exp(-7.0_dp * rphi / 180.0_dp)
      ! now convert to radians
      f2 = (rphi + sc_d * f) * pi / 180.0_dp
      ! fphi= (1+3.13*sc*exp(-7.0*phi/pi))*phi
      c = (cos(f2))**2
      p4 = -pc_d * p1 * (1.0_dp - c) / (1.0_dp + c)

    end subroutine dustmat_Draine

  attributes(device) subroutine dustmat_2(p1,p2,p3,p4,cosmu,cosmu2,alph,g1,g2)
    implicit none

    real(dp), intent(in) :: cosmu, cosmu2, alph, g1, g2
    real(dp), intent(out) :: p1, p2, p3, p4

    real(dp) :: bee, beecir, rtheta, B, Bcir, f, f2, c

    p1 = alph * (1.0_dp - g1**2) / (1.0_dp + g1**2 - 2.0_dp * g1 * cosmu)**1.5_dp + &
      & (1.0_dp - alph) * (1.0_dp - g2**2) / (1.0_dp + g2**2 - 2.0_dp * g2 * cosmu)**1.5_dp

    bee = 1.0_dp
    beecir = 1.0_dp

    rtheta = acos(cosmu) * 180.0_dp / pi
    if (rtheta < 35.0_dp) then
      B = bee*(1.0_dp - abs(rtheta - 145.0_dp)/35.0_dp)
    else
      B = 0.0_dp
    end if

    p2 = p1 * (-pl_d * (1.0_dp - cosmu2) / (1.0_dp + cosmu2) + B)
    p3 = p1 * (2.0_dp * cosmu) / (1.0_dp + cosmu2)

    if (rtheta >= 140.0_dp) then
      Bcir = -beecir * (1.0_dp - abs(rtheta - 160.0_dp)/20.0_dp)
    else if(rtheta >= 110.0_dp) then
      Bcir = 0.6_dp*beecir * (1.0_dp - ((rtheta - 110.0_dp)/30.0_dp)**2)
    else if(rtheta >= 55.0_dp) then
      Bcir = 0.6_dp*beecir * (rtheta - 55.0_dp)/55.0_dp
    else
      Bcir = 0.0_dp
    end if
    f = 3.13_dp * rtheta * exp(-7.0_dp * rtheta / 180.0_dp)
    ! now convert to radians
    f2 = (rtheta + sc_d * f) * pi / 180.0_dp
    ! fphi= (1+3.13*sc*exp(-7.0*phi/pi))*phi
    c = (cos(f2))**2
    p4 = p1 * (-pc_d * (1.0_dp - c) / (1.0_dp + c) + Bcir)

  end subroutine dustmat_2



end module mc_k_scatt_mat
