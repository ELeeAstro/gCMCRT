module mc_k_aux
  use mc_precision
  use mc_class_pac
  implicit none


contains

  attributes(device) subroutine locate(arr, var, idx)
    implicit none

    real(dp), dimension(:), intent(in) :: arr
    real(dp),intent(in) ::  var

    integer, intent(out) :: idx

    integer :: n, jl, jm, ju

    n = size(arr)

    ! Search an array using bi-section (numerical methods)
    ! Then return array index that is lower than var in arr
    jl = 0
    ju = n + 1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if ((arr(n) > arr(1)).eqv.(var > arr(jm))) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate

  subroutine locate_host(arr, var, idx)
    implicit none

    real(dp), dimension(:), intent(in) :: arr
    real(dp),intent(in) ::  var

    integer, intent(out) :: idx

    integer :: n, jl, jm, ju

    n = size(arr)

    ! Search an array using bi-section (numerical methods)
    ! Then return array index that is lower than var in arr
    jl = 0
    ju = n + 1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if ((arr(n) > arr(1)).eqv.(var > arr(jm))) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate_host

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    if (x1 >= x2) then
      print*, 'Error in linear_log_interp: x1 >= x2 - STOPPING', x1, x2
      stop
    end if

    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / log10(x2/x1)

    yval = 10.0_dp**((ly1 * log10(x2/xval) + ly2 * log10(xval/x1)) * norm)

  end subroutine linear_log_interp

  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp), intent(out) :: yval

    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp) :: norm

    if (x1 >= x2) then
      print*, 'Error in linear_interp: x1 >= x2 - STOPPING', x1, x2
      stop
    end if

    lxval = xval
    ly1 = y1
    ly2 = y2
    lx1 = x1
    lx2 = x2

    norm = 1.0_dp / (lx2 - lx1)

    yval = (ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm

  end subroutine linear_interp

  ! attributes(device) subroutine set_seed(ph)
  !   ! Note: this should be called with a different seed in each thread or
  !   ! process.
  !   implicit none
  !   type(pac),intent(inout) :: ph
  !
  !   ph%i = 97
  !   ph%j = 33
  !   ph%c_ran = 0.0_dp
  !
  !   ph%iseed = -abs(ph%seed)
  !   call set_seed_64(abs(ph%seed), 987654321, ph)
  ! end subroutine set_seed
  !
  ! attributes(device) subroutine set_seed_64(seed1,seed2,ph)
  !   implicit none
  !   integer,intent(in) :: seed1,seed2
  !   type(pac), intent(inout) :: ph
  !   integer i,j,x,y
  !   real(dp) :: s,t
  !   x=seed1
  !   y=seed2
  !   do i=1,97
  !      s=0.0_dp
  !      t=0.5_dp
  !      do j=1,53
  !         x=mod(6969*x,65543)
  !         y=mod(8888*x,65579)
  !         if (iand(ieor(x,y),32) > 0) then
  !           s=s+t
  !         end if
  !         t=0.5_dp*t
  !      end do
  !      ph%u(i)=s
  !   end do
  ! end subroutine set_seed_64
  !
  ! attributes(device) function ran2(ph)
  !   ! Random number between 0 and 1
  !   ! Based on "The 64-bit universal RNG", Marsaglia & Tsang (2004)
  !   ! Taken from Thomas P. Robitaille
  !   implicit none
  !   real(dp) :: xi
  !   type(pac), intent(inout) :: ph
  !   real(dp) :: x
  !   real(dp), parameter :: r = 9007199254740881.0_dp / 9007199254740992.0_dp
  !   real(dp), parameter :: d = 362436069876.0_dp / 9007199254740992.0_dp
  !   real(dp) :: ran2
  !
  !   x = ph%u(ph%i) - ph%u(ph%j)
  !   if (x < 0.0) then
  !      x = x + 1.0_dp
  !   end if
  !   ph%u(ph%i) = x
  !   ph%i = ph%i - 1
  !   if (ph%i == 0) then
  !     ph%i = 97
  !   end if
  !   ph%j = ph%j - 1
  !   if (ph%j == 0) then
  !     ph%j = 97
  !   end if
  !   ph%c_ran = ph%c_ran - d
  !   if (ph%c_ran < 0.0_dp) then
  !     ph%c_ran = ph%c_ran + r
  !   end if
  !   x = x - ph%c_ran
  !   xi = x
  !   if (x <= 0.0_dp) then
  !     xi = x + 1.0_dp
  !   end if
  !
  !   ran2 = xi

  ! end function ran2




end module mc_k_aux
