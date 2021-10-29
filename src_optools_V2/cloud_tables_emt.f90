module cloud_tables_emt
  use optools_data_mod
  use optools_aux, only : gauss
  implicit none

  private
  public :: emt_cl

contains

  subroutine emt_cl(z,n_work,k_work,N_eff)
    implicit none

    integer, intent(in) :: z
    integer :: s
    real(kind=dp), dimension(ncl), intent(in) :: n_work, k_work
    complex(kind=dp), intent(out) :: N_eff
    complex(kind=dp), dimension(ncl) :: N_inc, e_inc
    complex(kind=dp) :: e_eff0, N_eff0, e_eff
    real(kind=dp), dimension(ncl) :: Vol
    logical :: errflag

    Vol(:) = VMR_cl_lay(cl_tab(:)%iVMR,z)

    N_inc(:) = cmplx(n_work(:),k_work(:),kind=dp)

    N_eff = (0.0_dp, 0.0_dp)
    e_eff = (0.0_dp, 0.0_dp)

    errflag = .False.

    N_eff0 = (0.0_dp, 0.0_dp)
    do s = 1, ncl
      N_eff0 = N_eff0 + Vol(s) * N_inc(s)
      e_inc(s) = m2e(N_inc(s))
    end do

    ! print*, sum(nd_cl_lay(:,z))
    ! print*, N_inc
    ! print*, N_eff0
    ! print*, e_inc

    select case(imix)
    case(1)
      ! Try Bruggeman approach
      call NR(N_inc(:),Vol(:),N_eff0,N_eff,errflag)
      ! If NR fails use LLL method
      if (errflag .eqv. .True.) then
        e_eff0 = (0.0_dp,0.0_dp)
        do s = 1, ncl
          e_eff0 = e_eff0 + Vol(s)*e_inc(s)**(third)
        end do
        e_eff = e_eff0**3
        N_eff = e2m(e_eff)
        errflag = .False.
      end if

    case(2)
      ! Use LLL method
      e_eff0 = (0.0_dp,0.0_dp)
      do s = 1, ncl
        e_eff0 = e_eff0 + Vol(s)*e_inc(s)**(third)
      end do
      e_eff = e_eff0**3
      N_eff = e2m(e_eff)
      errflag = .False.

    case default
      print*, 'ERROR - EMT mixing method not valid - STOPPING'
      print*, 'imix: ', imix
      stop
    end select

! print*, N_eff
!
! stop

  end subroutine emt_cl

  ! ------------- Functions for LLL theory ------------- !!
  pure complex(kind=dp) function e2m(e)
    implicit none

    real(kind=dp) :: ereal, eimag, n, k
    real(kind=dp) :: sqrte2
    complex(kind=dp), intent(in) :: e

    ereal = real(e,kind=dp)
    eimag = aimag(e)
    sqrte2 = sqrt(ereal*ereal + eimag*eimag)
    n = sqrt(0.5_dp * ( ereal + sqrte2))
    k = sqrt(0.5_dp * (-ereal + sqrte2))
    e2m = cmplx(n,k,kind=dp)
  end function e2m

  pure complex(kind=dp) function m2e(m)
    implicit none

    real(kind=dp) :: ereal, eimag, n, k
    complex(kind=dp), intent(in) :: m

    n = real(m,kind=dp)
    k = aimag(m)
    ereal = n*n - k*k
    eimag = 2.0_dp * n * k
    m2e = cmplx(ereal,eimag,kind=dp)
  end function m2e

  !!! Combine using Bruggeman formula
  ! define function to be minimized
  pure subroutine Bruggeman(M_eff,M_inc,V_inc,FF)
    implicit none

    integer :: i
    real(kind=dp), dimension(:), intent(in) :: V_inc
    real(kind=dp), intent(out) :: FF(2)
    complex(kind=dp), dimension(:), intent(in) ::  M_inc
    complex(kind=dp), intent(in) :: M_eff
    complex(kind=dp) :: fun, mm2,mmi2

    mm2 = M_eff**2
    fun = cmplx(0.0_dp,0.0_dp,kind=dp)
    do i = 1, size(V_inc)
      mmi2 = M_inc(i)**2
      fun = fun + V_inc(i)*(mmi2 - mm2)/(mmi2 + 2.0_dp*mm2)
    end do

    FF(1) = real(fun,kind=dp)
    FF(2) = aimag(fun)

  end subroutine Bruggeman

  ! -------------------------------------------------------------------------------------------------
  ! A program for function minimization using the Newton-Raphson method.
  ! -------------------------------------------------------------------------------------------------
  pure subroutine NR(M_inc,V_inc,M_eff0,M_eff,unphysical)
    implicit none

    integer, parameter :: itmax = 30
    integer :: it
    real(kind=dp),dimension(:), intent(in)  :: V_inc
    real(kind=dp), dimension(2,2) :: DF
    real(kind=dp), dimension(2) ::  Fold, FF, FF1, FF2, FF3, FF4
    real(kind=dp), dimension(2) :: corr, xx, xnew
    real(kind=dp) :: qual, de1, de2
    complex(kind=dp), dimension(:), intent(in) :: M_inc
    complex(kind=dp), intent(in)  :: M_eff0
    complex(kind=dp), intent(out) :: M_eff
    logical, intent(inout) :: unphysical

    M_eff = M_eff0

    do it = 1, itmax
      xx(1) = real(M_eff,kind=dp)
      xx(2) = aimag(M_eff)
      call Bruggeman(M_eff,M_inc,V_inc,FF)
      qual = FF(1)*FF(1) + FF(2)*FF(2)
      de1 = xx(1)*1.0e-5_dp
      de2 = xx(2)*1.0e-5_dp
      call Bruggeman(M_eff+cmplx(de1,0.0_dp,kind=dp),M_inc,V_inc,FF1)
      call Bruggeman(M_eff-cmplx(de1,0.0_dp,kind=dp),M_inc,V_inc,FF2)
      call Bruggeman(M_eff+cmplx(0.0_dp,de2,kind=dp),M_inc,V_inc,FF3)
      call Bruggeman(M_eff-cmplx(0.0_dp,de2,kind=dp),M_inc,V_inc,FF4)
      DF(1,1) = (FF1(1)-FF2(1)) / (2.0_dp*de1)
      DF(1,2) = (FF3(1)-FF4(1)) / (2.0_dp*de2)
      DF(2,1) = (FF1(2)-FF2(2)) / (2.0_dp*de1)
      DF(2,2) = (FF3(2)-FF4(2)) / (2.0_dp*de2)
      Fold = FF
      call gauss(2,2,DF,corr,FF)
      corr = -corr
      call eff_pullback(2,xx,corr,Fold,xnew,unphysical,M_inc,V_inc)
      if (unphysical) then
  !        print*, qual, 'qual'
        exit
      end if
      M_eff = cmplx(xnew(1),xnew(2),kind=dp)
      if (abs(qual) < 1.0e-12_dp) then
        exit
      end if
    end do

  end subroutine NR

  ! Pull Back Eff
  pure subroutine eff_pullback(N,xx,dx,Fold,xnew,unphysical,M_inc,V_inc)
    implicit none

    integer, intent(in) :: N
    integer, parameter :: itmax=20
    integer :: it
    real(kind=dp), dimension(:),intent(in) :: V_inc
    real(kind=dp), dimension(N), intent(in) :: xx,Fold
    real(kind=dp), dimension(N), intent(inout) :: dx
    real(kind=dp), dimension(N), intent(out) :: xnew
    real(kind=dp), dimension(2) :: Fnew(2)
    real(kind=dp) :: fac,qold,qnew
    complex(kind=dp), dimension(:),intent(in) :: M_inc
    complex(kind=dp) :: mm
    logical, intent(out) :: unphysical

    qold = Fold(1)*Fold(1)+Fold(2)*Fold(2)
    fac = 1.0_dp

    do it = 1, itmax
      xnew = xx + fac*dx

      if ((xnew(1) > 0.0_dp).and.(xnew(2) > 0.0_dp)) then
        mm = cmplx(xnew(1),xnew(2),kind=dp)
        call Bruggeman(mm,M_inc,V_inc,Fnew)
        qnew = Fnew(1)*Fnew(1)+Fnew(2)*Fnew(2)
          !write(*,*) it,qold,qnew
        unphysical = .False.
        if (qnew < qold) then
          exit
        end if
      else
          !write(*,*) it,"negative (n,k)",xnew
        unphysical = .True.
      endif

        fac = fac*0.7_dp
    enddo

  end subroutine eff_pullback

end module cloud_tables_emt
