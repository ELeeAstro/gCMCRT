module LHS_sampling_mod
  use mc_precision
  implicit none

  public :: LHS_sample_2D
  !private ::

contains

  subroutine LHS_sample_2D(N, x, y, method, optimise, max_iter, opt_method)
    implicit none
    integer, intent(in) :: N, max_iter
    character(len=*), intent(in) :: method, opt_method
    logical, intent(in) :: optimise
    real(dp), intent(out) :: x(N), y(N)

    integer :: i, j, k, ia, ib, d_iter, col
    real(dp) :: r(N), d0, d1, tmp
    integer :: perm(N)
    real(dp), allocatable :: sample(:,:)

    call random_seed()

    select case (trim(method))
    case ('lhs')
      ! Generate base LHS sample in 2D
      allocate(sample(N,2))

      ! LHS for x
      call random_number(r)
      do i = 1, N
        sample(i,1) = (i - 1 + r(i)) / real(N,dp)
      end do
      call random_permutation(N, perm)
      sample(:,1) = sample(perm,1)

      ! LHS for y
      call random_number(r)
      do i = 1, N
        sample(i,2) = (i - 1 + r(i)) / real(N,dp)
      end do
      call random_permutation(N, perm)
      sample(:,2) = sample(perm,2)

      if (optimise) then
        d0 = compute_min_dist(sample(:,1), sample(:,2), N)

        select case (trim(opt_method))
        case ('maximin')
          ! Pairwise maximin full sweep
          do k = 1, max_iter
            d1 = d0
            ia = -1; ib = -1
            do i = 1, N-1
              do j = i+1, N
                call swap(sample(:,2), i, j)
                tmp = compute_min_dist(sample(:,1), sample(:,2), N)
                if (tmp > d1) then
                  d1 = tmp
                  ia = i
                  ib = j
                end if
                call swap(sample(:,2), i, j) ! revert
              end do
            end do
            if (ia > 0) then
              call swap(sample(:,2), ia, ib)
              d0 = d1
            else
              exit
            end if
          end do

        case ('random_cd')
          ! Random coordinate descent on y
          do d_iter = 1, max_iter
            call random_number(r)
            ia = 1 + int(r(1) * N)
            call random_number(r)
            ib = 1 + int(r(1) * N)
            if (ia == ib) cycle

            call swap(sample(:,2), ia, ib)
            d1 = compute_min_dist(sample(:,1), sample(:,2), N)

            if (d1 >= d0) then
              d0 = d1
            else
              call swap(sample(:,2), ia, ib)
            end if
          end do

        case ('cpe')
          ! Columnwise pairwise exchange
          do k = 1, max_iter
            d1 = d0
            ia = -1; ib = -1
            do i = 1, N-1
              do j = i+1, N
                call swap(sample(:,2), i, j)
                tmp = compute_min_dist(sample(:,1), sample(:,2), N)
                if (tmp > d1) then
                  d1 = tmp
                  ia = i
                  ib = j
                end if
                call swap(sample(:,2), i, j)
              end do
            end do
            if (ia > 0) then
              call swap(sample(:,2), ia, ib)
              d0 = d1
            else
              exit
            end if
          end do

        case default
          print *, 'Unknown optimisation method: ', trim(opt_method)
          stop 1
        end select
      end if

      x = sample(:,1)
      y = sample(:,2)
      deallocate(sample)

    case ('random')
      call random_number(x)
      call random_number(y)

    case default
      print *, "Unknown sampling method: ", trim(method)
      stop 1
    end select
  end subroutine LHS_sample_2D

  function compute_min_dist(x, y, N) result(dmin)
    implicit none
    integer, intent(in) :: N
    real(dp), intent(in) :: x(N), y(N)
    real(dp) :: dmin, dx, dy, d2
    integer :: i, j

    dmin = HUGE(0.0_dp)
    do i = 1, N - 1
      do j = i + 1, N
        dx = x(i) - x(j)
        dy = y(i) - y(j)
        d2 = dx*dx + dy*dy
        if (d2 < dmin) dmin = d2
      end do
    end do
    dmin = sqrt(dmin)
  end function compute_min_dist

  subroutine swap(a, i, j)
    implicit none
    real(dp), intent(inout) :: a(:)
    integer, intent(in) :: i, j
    real(dp) :: tmp
    tmp = a(i)
    a(i) = a(j)
    a(j) = tmp
  end subroutine swap

  subroutine random_permutation(N, perm)
    implicit none
    integer, intent(in) :: N
    integer, intent(out) :: perm(N)
    integer :: i, j, tmp
    real(dp) :: r

    do i = 1, N
      perm(i) = i
    end do

    do i = N, 2, -1
      call random_number(r)
      j = int(r * i) + 1
      tmp = perm(i)
      perm(i) = perm(j)
      perm(j) = tmp
    end do
  end subroutine random_permutation

end module LHS_sampling_mod
