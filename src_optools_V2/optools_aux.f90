module optools_aux
  use optools_data_mod, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none

  integer, parameter :: locate_below = -1
  integer, parameter :: locate_inside = 0
  integer, parameter :: locate_above = 1

  private
  public :: locate, locate_triplet, locate_below, locate_inside, locate_above, &
    & sort2, gauss, &
    & linear_interp, linear_log_interp, bilinear_interp, bilinear_log_interp, &
    & trapz, Bezier_interp

contains

  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(kind=dp), intent(in) :: xval, y1, y2, x1, x2
    real(kind=dp), intent(out) :: yval
    real(kind=dp) :: norm

    if (x1 >= x2) then
      print*, 'Error in linear_interp: x1 >= x2 - STOPPING', x1, x2
      stop
    end if

    norm = 1.0_dp / (x2 - x1)

    yval = (y1 * (x2 - xval) + y2 * (xval - x1)) * norm

  end subroutine linear_interp

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / log10(x2/x1)

    yval = 10.0_dp**((ly1 * log10(x2/xval) + ly2 * log10(xval/x1)) * norm)

  end subroutine linear_log_interp

  subroutine bilinear_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
    implicit none

    real(kind=dp), intent(in) :: xval, yval, x1, x2, y1, y2, a11, a21, a12, a22
    real(kind=dp), intent(out) :: aval
    real(kind=dp) :: norm

    if (x1 >= x2) then
      print*, 'Error in bi_linear_interp: x1 >= x2 - STOPPING', x1, x2
      stop
    else if (y1 >= y2) then
      print*, 'Error in bi_linear_interp: y1 >= y2 - STOPPING', y1, y2
      stop
    end if

    norm = 1.0_dp / (x2 - x1) / (y2 - y1)

    aval = a11 * (x2 - xval) * (y2 - yval) * norm &
      & + a21 * (xval - x1) * (y2 - yval) * norm &
      & + a12 * (x2 - xval) * (yval - y1) * norm &
      & + a22 * (xval - x1) * (yval - y1) * norm

  end subroutine bilinear_interp

  subroutine bilinear_log_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
    implicit none

    real(kind=dp), intent(in) :: xval, yval, x1, x2, y1, y2, a11, a21, a12, a22
    real(kind=dp) :: lxval, lyval, lx1, lx2, ly1, ly2, la11, la21, la12, la22
    real(kind=dp), intent(out) :: aval
    real(kind=dp) :: norm

    lxval = log10(xval) ; lyval = log10(yval)
    lx1 = log10(x1) ; lx2 = log10(x2) ; ly1 = log10(y1) ; ly2 = log10(y2)
    la11 = log10(a11) ; la21 = log10(a21) ; la12 = log10(a12) ; la22 = log10(a22)

    if (lx1 >= lx2) then
      print*, 'Error in bilinear_log_interp: lx1 >= lx2 - STOPPING', lx1, lx2, x1, x2
      stop
    else if (ly1 >= ly2) then
      print*, 'Error in bilinear_log_interp: ly1 >= ly2 - STOPPING', ly1, ly2, y1, y2
      stop
    end if

    norm = 1.0_dp / (lx2 - lx1) / (ly2 - ly1)

    aval = la11 * (lx2 - lxval) * (ly2 - lyval) * norm &
      & + la21 * (lxval - lx1) * (ly2 - lyval) * norm &
      & + la12 * (lx2 - lxval) * (lyval - ly1) * norm &
      & + la22 * (lxval - lx1) * (lyval - ly1) * norm

    aval = 10.0_dp**(aval)

  end subroutine bilinear_log_interp

  ! Perform Bezier interpolation
  subroutine Bezier_interp(xi, yi, x, y)
    implicit none

    real(dp), dimension(3), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: dx, dx1, dy, dy1, w, yc, t, wlim, wlim1
    real(dp) :: denom, denom1, grad_scale, grad_tol, y_linear

    if (any(.not. ieee_is_finite(xi)) .or. any(.not. ieee_is_finite(yi)) .or. &
      & .not. ieee_is_finite(x)) then
      print*, 'Error in Bezier_interp: non-finite input - STOPPING', x, xi, yi
      stop
    end if

    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)

    if (dx <= 0.0_dp .or. dx1 <= 0.0_dp) then
      print*, 'Error in Bezier_interp: stencil grid must be strictly increasing - STOPPING', xi
      stop
    end if

    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    ! Return exact endpoints and clamp any round-off excursion outside the stencil.
    if (x <= xi(1)) then
      y = yi(1)
      return
    else if (x >= xi(3)) then
      y = yi(3)
      return
    else if (x == xi(2)) then
      y = yi(2)
      return
    end if

    ! Calculate the ordinary linear result for the active interval. This is
    ! the safe fallback for flat or ill-conditioned Bezier stencils.
    if (x < xi(2)) then
      t = (x - xi(1))/dx
      y_linear = (1.0_dp - t)*yi(1) + t*yi(2)
    else
      t = (x - xi(2))/dx1
      y_linear = (1.0_dp - t)*yi(2) + t*yi(3)
    end if

    grad_scale = max(1.0_dp, abs(yi(1)), abs(yi(2)), abs(yi(3)))
    grad_tol = 100.0_dp*epsilon(1.0_dp)*grad_scale

    ! The limiter ratios divide by both gradients. A flat interval or a
    ! change in gradient sign is better represented by the linear fallback.
    if (abs(dy) <= grad_tol .or. abs(dy1) <= grad_tol .or. &
      & (dy > 0.0_dp .and. dy1 < 0.0_dp) .or. &
      & (dy < 0.0_dp .and. dy1 > 0.0_dp)) then
      y = y_linear
      return
    end if

    denom = 1.0_dp - (dy1/dy)*(dx/dx1)
    denom1 = 1.0_dp - (dy/dy1)*(dx1/dx)

    if (.not. ieee_is_finite(denom) .or. .not. ieee_is_finite(denom1) .or. &
      & abs(denom) <= 100.0_dp*epsilon(1.0_dp) .or. &
      & abs(denom1) <= 100.0_dp*epsilon(1.0_dp)) then
      y = y_linear
      return
    end if

    if (x < xi(2)) then
      ! left hand side interpolation
      w = dx1/(dx + dx1)
      wlim = 1.0_dp + 1.0_dp/denom
      wlim1 = 1.0_dp/denom1
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) - dx/2.0_dp * (w*dy/dx + (1.0_dp - w)*dy1/dx1)
      t = (x - xi(1))/dx
      if (.not. ieee_is_finite(yc) .or. yc < min(yi(1),yi(2)) .or. &
        & yc > max(yi(1),yi(2))) then
        y = y_linear
        return
      end if
      y = (1.0_dp - t)**2*yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else
      ! right hand side interpolation
      w = dx/(dx + dx1)
      wlim = 1.0_dp/denom
      wlim1 = 1.0_dp + 1.0_dp/denom1
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (w*dy1/dx1 + (1.0_dp - w)*dy/dx)
      t = (x - xi(2))/dx1
      if (.not. ieee_is_finite(yc) .or. yc < min(yi(2),yi(3)) .or. &
        & yc > max(yi(2),yi(3))) then
        y = y_linear
        return
      end if
      y = (1.0_dp - t)**2*yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

    if (.not. ieee_is_finite(y)) then
      y = y_linear
    end if

  end subroutine Bezier_interp

  pure subroutine locate(arr, var, idx)
    implicit none

    integer, intent(out) :: idx
    real(kind=dp), dimension(:), intent(in) :: arr
    real(kind=dp),intent(in) ::  var
    integer :: jl, jm, ju

    ! Search an array using bi-section/binary search (numerical methods)
    ! Then return array index that is lower than var in arr

    jl = 0
    ju = size(arr)+1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if (var > arr(jm)) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate

  ! Locate a consecutive three-point stencil on a strictly increasing grid.
  ! Region reports whether var is below, inside, or above the full grid, while
  ! exact_idx identifies exact grid matches so callers can bypass interpolation.
  subroutine locate_triplet(arr, var, idx, region, exact_idx)
    implicit none

    real(kind=dp), dimension(:), intent(in) :: arr
    real(kind=dp), intent(in) :: var
    integer, dimension(3), intent(out) :: idx
    integer, intent(out) :: region, exact_idx
    integer :: lower, start

    if (size(arr) < 3) then
      print*, 'Error in locate_triplet: at least 3 grid points are required - STOPPING', size(arr)
      stop
    end if

    if (.not. ieee_is_finite(var)) then
      print*, 'Error in locate_triplet: target is not finite - STOPPING', var
      stop
    end if

    call locate(arr,var,lower)
    exact_idx = 0

    if (var < arr(1)) then
      region = locate_below
    else if (var > arr(size(arr))) then
      region = locate_above
    else
      region = locate_inside
      if (lower < size(arr)) then
        if (var == arr(lower+1)) exact_idx = lower + 1
      end if
    end if

    start = max(1,min(lower-1,size(arr)-2))
    idx = [start,start+1,start+2]

  end subroutine locate_triplet

  pure subroutine sort2(N,RA,RB)
    integer, intent(in) :: N
    integer :: L, IR, I, J
    real(kind=dp), dimension(N), intent(inout) :: RA, RB
    real(kind=dp) :: RRA, RRB
    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
      L=L-1
      RRA=RA(L)
      RRB=RB(L)
    ELSE
      RRA=RA(IR)
      RRB=RB(IR)
      RA(IR)=RA(1)
      RB(IR)=RB(1)
      IR=IR-1
      IF(IR.EQ.1)THEN
        RA(1)=RRA
        RB(1)=RRB
        RETURN
      ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
      IF(J.LT.IR)THEN
        IF(RA(J).LT.RA(J+1)) THEN
          J=J+1
        ENDIF
      ENDIF
      IF(RRA.LT.RA(J))THEN
        RA(I)=RA(J)
        RB(I)=RB(J)
        I=J
        J=J+J
      ELSE
        J=IR+1
      ENDIF
      GO TO 20
    ENDIF
    RA(I)=RRA
    RB(I)=RRB
    GO TO 10
  end subroutine sort2

  !**********************************************************************
    pure subroutine gauss(Nd,N,a,x,b)
      implicit none
  !**********************************************************************
  !*****                                                            *****
  !*****   This routine solves the linear system                     *****
  !*****   (( a )) * ( x ) = ( b ) for x.                            *****
  !*****   The algorithm reduces matrix a to triangular form.        *****
  !*****                                                            *****
  !*****   INPUT:   Nd = dimensions of the vectors and matrix        *****
  !*****             N = dimension of the linear system (N<=Nd)      *****
  !*****              a = (N x N) matrix                            *****
  !*****             b = vector of length N                          *****
  !*****   OUTPUT:    x = vector of length N                          *****
  !*****                                                            *****
  !**********************************************************************
  !*
      integer, intent(in) :: Nd, N
      integer :: i, j, k, kmax
      real(kind=dp), dimension(Nd,Nd), intent(inout) :: a
      real(kind=dp), dimension(Nd), intent(inout) :: b
      real(kind=dp), dimension(Nd), intent(out) :: x
      real(kind=dp) :: c, amax

      do i = 1, N-1
  !*       ------------------------------------------
  !*       ***  Maximum-pivot row exchange for row i  ***
  !*       ------------------------------------------
        kmax = i
        amax = abs(a(i,i))
        do k = i+1, N
          if (abs(a(k,i)) > amax) then
            amax = abs(a(k,i))
            kmax = k
          endif
        end do

        if (kmax /= i) then
          do j = 1, N
            c = a(i,j)
            a(i,j) = a(kmax,j)
            a(kmax,j) = c
          end do
          c = b(i)
          b(i) = b(kmax)
          b(kmax) = c
        end if
  !*
  !*       ---------------------------------
  !*       ***  Reduce to triangular form  ***
  !*       ---------------------------------
        do k = i+1, N
          c = a(k,i) / a(i,i)
          a(k,i) = 0.0_dp
          do j = i+1, N
            a(k,j) = a(k,j) - c * a(i,j)
          end do
          b(k) = b(k) - c * b(i)
        end do
  !*
      end do
  !*
  !*     --------------------------
  !*     ***  Solve for x  ***
  !*     --------------------------
      do i = N, 1, -1
        c = 0.0_dp
        if (i < N) then
          do j = i+1, N
            c = c + a(i,j) * x(j)
          end do
        end if
        x(i) = (b(i) - c) / a(i,i)
      end do

    end subroutine gauss

    !! Integration functions from fortranwiki

    pure function trapz(x, y) result(r)
      !! Calculates the integral of an array y with respect to x using the trapezoid
      !! approximation. Note that the mesh spacing of x does not have to be uniform.
      real(kind=dp), intent(in)  :: x(:)         !! Variable x
      real(kind=dp), intent(in)  :: y(size(x))   !! Function y(x)
      real(kind=dp)              :: r            !! Integral ∫y(x)·dx

      ! Integrate using the trapezoidal rule
      associate(n => size(x))
        r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2.0_dp
      end associate
    end function trapz

    ! function PCHI(x, y, a, b) result(r)
    !   !! This function constructs a piecewise cubic Hermitian interpolation of an array y(x) based on
    !   !! discrete numerical data, and subsequently evaluates the integral of the interpolation in the
    !   !! range (a,b). Note that the mesh spacing of x does not necessarily have to be uniform.
    !   real(kind=dp), intent(in)  :: x(:)        !! Variable x
    !   real(kind=dp), intent(in)  :: y(size(x))  !! Function y(x)
    !   real(kind=dp), intent(in)  :: a           !! Left endpoint
    !   real(kind=dp), intent(in)  :: b           !! Right endpoint
    !   real(kind=dp)              :: r           !! Integral ∫y(x)·dx
    !
    !   real(kind=dp), external    :: dpchqa
    !   real(kind=dp)              :: d(size(x))
    !   integer               :: err
    !
    !   ! Create a PCHIP interpolation of the input data
    !   !call dpchez(size(x), x, y, d, .false., 0, 0, err)
    !
    !   ! Integrate the interpolation in the provided range
    !   r = dpchqa(size(x), x, y, d, a, b, err)
    ! end function PCHI


end module optools_aux
