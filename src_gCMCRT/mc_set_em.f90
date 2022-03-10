module mc_set_em
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_opacset, only : k_tot_abs
  implicit none


  real(dp), parameter :: tau_cuttoff = 30.0_dp
  logical :: first_call = .True.
  real(dp) :: n_cells
  integer, allocatable, dimension(:,:) :: itau3

  private :: F

contains

  subroutine set_grid_em(l,n)
    implicit none

    integer, intent(in) :: l
    integer, intent(in), optional :: n
    integer :: i, j, k, g,  itau, itaul

    real(dp) :: tau_sum, densav, wl_eff, BBf

    if (first_call .eqv. .True.) then
      !! Allocate all arrays
      allocate(itau3(grid%n_phi-1,grid%n_theta-1))
      allocate(l_cell(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      if (ck .eqv. .True.) then
        allocate(l_cell_g(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
        allocate(cell_gord_cdf(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
        allocate(cell_gord_cdf_d(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
        allocate(cell_gord_wght(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
        allocate(cell_gord_wght_d(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      end if
      first_call = .False.
    end if

    if (oneD .eqv. .True.) then
      !! If 1D profile input then only 1 vertical profile required, use index j = 1, k = 1

      ! Calculate
      itau3(1,1) = grid%n_lev - 1
      tau_sum = 0.0_dp
       do i = grid%n_lev, 2, -1
         itaul = i - 1
          densav = 0.0_dp
          !if (ck .eqv. .True.) then
            !do g = 1, ng
            !  densav = densav + gord_w(g) * rhokap(g,i-1,1,1)/grid%r_del
            !end do
          !else if (lbl .eqv. .True.) then
            densav = rhokap(1,i-1,1,1)/grid%r_del
          !end if
          tau_sum = tau_sum + densav * abs((H(i)-H(i-1)))
           if (tau_sum > tau_cuttoff) then
             if (itaul < itau3(1,1)) then
               itau3(1,1) = itaul
             end if
             exit
           end if
           if (i == 2) then
             itau3(1,1) = itaul
           end if
       end do

       itau3(:,:) = itau3(1,1)
       print*, 'itau: ', tau_sum, itau3(1,1), PG(itau3(1,1),1,1)/bar

     else if (threeD .eqv. .True.) then

       ! Calculate
       itau3(:,:) = grid%n_lev - 1
       do k = 1, grid%n_theta-1
         do j = 1, grid%n_phi-1
           tau_sum = 0.0_dp
            do i = grid%n_lev, 2, -1
              itaul = i - 1
               densav = 0.0_dp
               !if (ck .eqv. .True.) then
                 !do g = 1, ng
                   densav = densav + gord_w(g) * rhokap(g,i-1,j,k)/grid%r_del
                 !end do
               !else if (lbl .eqv. .True.) then
                 densav = rhokap(1,i-1,j,k)/grid%r_del
               !end if
               tau_sum = tau_sum + densav * abs((H(i)-H(i-1)))
                if (tau_sum > tau_cuttoff) then
                  if (itaul < itau3(j,k)) then
                    itau3(j,k) = itaul
                  end if
                  exit
                end if
            end do
          end do
        end do

        print*, 'min, max itau: ', minval(itau3(:,:)),maxval(itau3(:,:))


     end if


       n_cells = 0.0_dp
       grid%lumtot = 0.0_dp

         do k = 1, grid%n_theta-1
           do j = 1, grid%n_phi-1

              !if (phiarr(j) > pi/2.0_dp .and. phiarr(j) < (pi + pi/2.0_dp)) then
               !  l_cell(:,j,k) = 0.0_dp
               !  cell_gord_cdf(:,:,j,k) = 0.0_dp
               !  cell_gord_wght(:,:,j,k) = 0.0_dp
               !  cycle
               !end if

             do i = 1, grid%n_lay

               if (i <= itau3(j,k) .or. i == grid%n_lay) then
                 l_cell(i,j,k) = 0.0_dp
                 if (ck .eqv. .True.) then
                   cell_gord_cdf(:,i,j,k) = 0.0_dp
                end if
                 cycle
               end if

               if (lbl .eqv. .True.) then
                 if (doppler_on .eqv. .True.) then
                   wl_eff = wl(l)*(1.0_dp - v_los(n,i,j,k)/c_s)
                   l_cell(i,j,k) = fourpi * RH(i,j,k) * v_cell(i,j,k) * k_tot_abs(1,i,j,k) * BB(wl_eff,TG(i,j,k))
                   !print*, i,j,k,wl_eff, BB(wl_eff,TG(i,j,k)), k_tot_abs(1,i,j,k) , l_cell(i,j,k)
                 else
                   l_cell(i,j,k) = fourpi * RH(i,j,k) * v_cell(i,j,k) * k_tot_abs(1,i,j,k) * BB(wl(l),TG(i,j,k))
                 end if
                 grid%lumtot = grid%lumtot  + l_cell(i,j,k)
                 n_cells = n_cells + 1.0_dp
                 cycle
               end if

               if (do_BB_band .eqv. .True.) then
                 !print*, TG(i,j,k), wl_e(l),wl_e(l+1)
                 !BBf = (BB_band(wl_e(l),wl_e(l+1),TG(i,j,k)) * 1000.0_dp) / (wl(l) * 1.0e-4_dp)
                 ! print*,'band', i, wl_e(l),wl_e(l+1), BBf
                 !BBf = (tplkavg(wl_e(l),wl_e(l+1),TG(i,j,k)) * 1000.0_dp) !/ (wl(l) * 1.0e-4_dp)
                 ! print*,'tplkavg', i, wl_e(l),wl_e(l+1), BBf
                 !BBf = planck1(wl(l),TG(i,j,k),wl_e(l),wl_e(l+1)) * (1.0_dp / (wl(l) * 1.0e-4_dp)) !/ (wl(l) * 1.0e-4_dp) !* 10.0_dp
                 ! print*,'planck1', i, wl_e(l),wl_e(l+1), BBf
                 BBf = BB(wl(l),TG(i,j,k))
                 ! print*,'single', i, wl(l), BBf
               else
                 BBf = BB(wl(l),TG(i,j,k))
                 !print*,'single', i, wl(l), BBf * (wl(l) * 1e-4_dp)
               end if


              do g = 1, ng
                l_cell_g(g,i,j,k) = gord_w(g) * fourpi * RH(i,j,k) * v_cell(i,j,k) * k_tot_abs(g,i,j,k) * BBf
              end do
              l_cell(i,j,k) = sum(l_cell_g(:,i,j,k))
              grid%lumtot = grid%lumtot + l_cell(i,j,k)
              n_cells = n_cells + 1.0_dp

              !print*, i, l_cell(i,j,k), k_tot_abs(1,i,j,k), k_tot_abs(16,i,j,k)

         cell_gord_cdf(1,i,j,k) = 0.0_dp
          do g = 2, ng
              cell_gord_cdf(g,i,j,k) = cell_gord_cdf(g-1,i,j,k) + l_cell_g(g,i,j,k)/l_cell(i,j,k)
          end do


          do g = 1, ng
            !cell_gord_wght(g,i,j,k) = (real(ng,dp) * gord_w(g) *  l_cell_g(g,i,j,k))/l_cell(i,j,k)
            cell_gord_wght(g,i,j,k) = (1.0_dp / ((1.0_dp - xi_g) + &
            & xi_g * (l_cell(i,j,k) / real(ng,dp) &
            & / l_cell_g(g,i,j,k))))
          end do

      end do
      !stop
    end do
  end do

    if (ck .eqv. .True.) then
      ! Send relevent GPU data
      cell_gord_cdf_d(:,:,:,:) = cell_gord_cdf(:,:,:,:)
      cell_gord_wght_d(:,:,:,:) = cell_gord_wght(:,:,:,:)
    end if

  end subroutine set_grid_em


  real(dp) function BB(wl_in, T_in)
    implicit none

    ! Planck function in wavelength units

    real(kind=dp), intent(in) :: wl_in, T_in
    real(kind=dp) :: left, right
    real(kind=dp) :: wl_cm

    wl_cm = wl_in * 1.0e-4_dp

    left = (2.0_dp * hpl * c_s**2)/(wl_cm)**5
    right = 1.0_dp / (exp((hpl * c_s) / (wl_cm * kb * T_in)) - 1.0_dp)
    BB = left * right

  end function BB


  real(dp) function BB_band(wl1_in, wl2_in, T_in)
    implicit none

    real(dp), intent(in) :: wl1_in, wl2_in, T_in

    integer :: i, b, j, intitera
    real(dp) :: iB1, iB2, wn1, wn2
    real(dp) :: x, x2, x3, itera, summ, dn

      !! Constants for BB_band function
    real(dp), parameter :: kb_mks = 1.380649e-23_dp ! J K-1
    real(dp), parameter :: h_mks = 6.62607015e-34_dp ! J s
    real(dp), parameter :: cs_mks = 2.99792458e8_dp ! m s-1
    real(dp), parameter :: c1 = (h_mks * cs_mks) / kb_mks
    real(dp), parameter :: c2 = cs_mks**2
    real(dp), parameter :: n2 = 2.0_dp * h_mks * c2

    !! Code for integrating the blckbody function between two wavenumbers
    !! This is a method that uses a sum convergence
    !! Taken from: spectralcalc.com/blackbody/inband_radiance.html

    if (T_in < 1.0e-6_dp) then

      BB_band = 1.0e-99_dp

    else

      ! First wavenumber

      wn1 = 1.0_dp/(wl1_in * 1e-4_dp)

      x = c1 * 100.0_dp * wn1/T_in
      x2 = x**2
      x3 = x**3

      itera = 2.0_dp + 20.0_dp/x
      if (itera > 512) then
        itera = 512
      end if
      intitera = int(itera)

      summ = 0.0_dp
      do j = 1, intitera + 1
        dn = 1.0_dp/real(j,dp)
        summ = summ + exp(-min(real(j,dp)*x,300.0_dp)) * &
        & (x3 + (3.0_dp * x2 + 6.0_dp*(x+dn)*dn)*dn)*dn
      end do

      iB1 = n2 * (T_in/c1)**4 * summ

      ! Second wavenumber

      wn2 = 1.0_dp/(wl2_in * 1e-4_dp)

      x = c1 * 100.0_dp * wn2/T_in
      x2 = x**2
      x3 = x**3

      itera = 2.0_dp + 20.0_dp/x
      if (itera > 512) then
        itera = 512
      end if
      intitera = int(itera)

      summ = 0.0_dp
      do j = 1, intitera + 1
        dn = 1.0_dp/real(j,dp)
        summ = summ + exp(-min(real(j,dp)*x,300.0_dp)) * &
        & (x3 + (3.0_dp * x2 + 6.0_dp*(x+dn)*dn)*dn)*dn
      end do

      iB2 = n2 * (T_in/c1)**4 * summ

      BB_band = max((iB2 - iB1), 1.0e-99_dp) ! Units W m-2 str-1

    end if

  end function BB_band

  real(dp) function F(x)
    implicit none

    real(dp), intent(in) :: x
    F = x**3 / (exp(x) - 1.0_dp)
  end function F

  real(dp) function tplkavg (wl1_in, wl2_in, T_in)
    implicit none

    real(dp), intent(in) :: wl1_in, wl2_in, T_in

    real(dp) :: WNUMLO, WNUMHI, T
    real(dp), parameter :: A1 = 1.0_dp/3.0_dp
    real(dp), parameter :: A2 = -1.0_dp/8.0_dp
    real(dp), parameter :: A3 = 1.0_dp/60.0_dp
    real(dp), parameter :: A4 = -1.0_dp/5040.0_dp
    real(dp), parameter :: A5 = 1.0_dp/272160.0_dp
    real(dp), parameter :: A6 = -1.0_dp/13305600.0_dp
    real(dp) :: EXM,X,arg,HH,DEL,c1,OLDVAL,VAL,VAL0
    real(dp) ::  wvn
    integer ::  I,N,K,M,MMAX

    integer ::  SMALLV
    real(dp), parameter :: C2 = 1.438786_dp
    real(dp), parameter :: SIGMA = 5.670374419e-8_dp
    real(dp), parameter :: VCUT = 1.5_dp
    real(dp), dimension(7), parameter :: VCP = (/ 10.25_dp, 5.7_dp, 3.9_dp, 2.9_dp, 2.3_dp, 1.9_dp, 0.0_dp /)
    real(dp), parameter :: VMAX = log(HUGE(1.0_dp))
    real(dp), parameter :: EPSIL = EPSILON(1.0_dp)
    real(dp), parameter :: SIGDPI = SIGMA/pi
    real(dp), parameter :: CONC = 15.0_dp/pi**4

    real(dp), dimension(2) :: D, P, V
    real(dp) :: EX, MV, VSQ

!      F(X) = X**3 / ( DEXP(X) - 1 )

    T = T_in

    ! VMAX = log(HUGE(1.0_dp))
    ! EPSIL = EPSILON(1.0_dp)
    ! SIGDPI = SIGMA/pi
    ! CONC = 15.0_dp/pi**4

    WNUMLO = 1.0_dp/(wl2_in*1e-4_dp)
    WNUMHI = 1.0_dp/(wl1_in*1e-4_dp)

    if ( T < 0.0_dp .or. WNUMHI < WNUMLO .OR. WNUMLO < 0.0_dp) then
      print*, 'PLKAVG--TEMPERATURE OR WAVENUMS. WRONG'
    end if

    IF ( T < 1.0E-4_dp )  THEN
       TPLKAVG = 0.0_dp
       RETURN
    ENDIF

    IF ( wnumhi == wnumlo ) THEN
       wvn  =  wnumhi
       arg  = exp( - C2 * wvn / T)
       tplkavg = c1 * (wvn**3) * arg / ( 1.0_dp - arg )
       RETURN
    ENDIF

    V(1) = C2 * WNUMLO / T
    V(2) = C2 * WNUMHI / T

    IF ( V(1).GT.EPSIL .AND. V(2).LT.VMAX .AND. &
 &     (WNUMHI-WNUMLO)/WNUMHI .LT. 1.0e-2_dp )  THEN

!                          ** WAVENUMBERS ARE VERY CLOSE.  GET INTEGRAL
!                          ** BY ITERATING SIMPSON RULE TO CONVERGENCE.
     HH = V(2) - V(1)
     OLDVAL = 0.0_dp
     VAL0 = F( V(1) ) + F( V(2) )

     DO N = 1, 10
        DEL = HH / (2.0_dp*real(N,dp))
        VAL = VAL0
        DO K = 1, 2*N-1
           VAL = VAL + 2.0_dp*(1.0_dp+MOD(K,2)) * F( V(1) + real(K,dp)*DEL )
        end do
        VAL = DEL/3.0_dp * VAL
        IF ( ABS( (VAL-OLDVAL)/VAL ) .LE. 1.0E-6_dp) then
          TPLKAVG = SIGDPI * T**4 * CONC * VAL
          return
        end if
        OLDVAL = VAL
     end do

     TPLKAVG = SIGDPI * T**4 * CONC * VAL
     return

    END IF

    SMALLV = 0
    DO I = 1, 2

      IF( V(I).LT.VCUT )  THEN
!                                   ** USE POWER SERIES
        SMALLV = SMALLV + 1
        VSQ = V(I)**2
        P(I) =  CONC * VSQ * V(I) * ( A1 + V(I) * ( A2 + V(I) * &
 &                ( A3 + VSQ * ( A4 + VSQ * ( A5 + VSQ*A6 ) ) ) ) )
      ELSE

!                    ** USE EXPONENTIAL SERIES
        MMAX = 0
!                                ** FIND UPPER LIMIT OF SERIES
        do M = 1, 7
          MMAX = MMAX + 1
          IF ( V(I).LT.VCP( MMAX ) ) then
            continue
          else
            exit
          end if
        end do

        EX = exp( - V(I) )
        EXM = 1.0_dp
        D(I) = 0.0_dp

        DO M = 1, MMAX
           MV = M * V(I)
           EXM = EX * EXM
           D(I) = D(I) + &
 &                EXM * ( 6.0_dp + MV*( 6.0_dp + MV*( 3.0_dp + MV ) ) ) / M**4
        end do

        D(I) = CONC * D(I)
     END IF
    end do

  IF ( SMALLV .EQ. 2 ) THEN
!                                    ** WNUMLO AND WNUMHI BOTH SMALL
     TPLKAVG = P(2) - P(1)

  ELSE IF ( SMALLV .EQ. 1 ) THEN
!                                    ** WNUMLO SMALL, WNUMHI LARGE
     TPLKAVG = 1.0_dp - P(1) - D(2)

  ELSE
!                                    ** WNUMLO AND WNUMHI BOTH LARGE
     TPLKAVG = D(1) - D(2)

  END IF

    TPLKAVG = SIGDPI * T**4 * TPLKAVG
    IF( TPLKAVG.EQ.0.0_dp ) then
     print*, 'PLKAVG--RETURNS ZERO; POSSIBLE UNDERFLOW'
    end if

  end function tplkavg

  real(dp) function planck1(wl_in, T_in, wl1_e_in, wl2_e_in)
    implicit none

    real(dp), intent(in) :: wl_in, T_in, wl1_e_in, wl2_e_in

    integer, parameter :: NBB = 4

    integer :: i
    real(dp) :: wav, dw
    real(dp) :: planck, wavnum, expterm

    wav = 1.0_dp/(wl_in*1.0e-4_dp)
    dw = 1.0_dp/(wl1_e_in*1.0e-4_dp) - 1.0_dp/(wl2_e_in*1.0e-4_dp)

    planck = 0.0_dp
    do i = -NBB, NBB, 1
      wavnum = wav + real(i,dp)*dw/(2.0_dp*real(NBB,dp))
      expterm = exp(-1.438769_dp * wavnum/T_in)
      planck = planck + 1.191e-5_dp * wavnum**3 * expterm / &
      & (1.0_dp - expterm)
    end do
    planck1 = planck/(2.0_dp * real(NBB,dp) + 1.0_dp)
    if (planck1 < 1.0e-300_dp) then
      planck1 = 1.0e-300_dp
    end if

  end function planck1


end module mc_set_em
