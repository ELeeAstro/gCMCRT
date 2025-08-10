module LHS_sampling_mod
  use mc_precision
  implicit none

  public :: LHS_sample
  private :: random_permutation

contains

  !--------------------------------------------------------------------
  ! Generate a Latin Hypercube sample sequence (no optimisation).
  !
  ! Inputs:
  !   N         : number of points
  !   centered  : (optional) if .true., place points at bin centres
  !               (i-0.5)/N; otherwise jitter within bins (default)
  !
  ! Outputs:
  !   x, y, z   : arrays of size N with values in [0,1)
  !
  ! Notes:
  !   - This routine does NOT call random_seed(). Seed once in your
  !     main program if you want non-reproducible runs each time.
  !   - Each axis is stratified into N bins and independently permuted.
  !--------------------------------------------------------------------
  subroutine LHS_sample(N, Ndim, x, y, z, center)
    implicit none

    integer, intent(in) :: N, Ndim
    logical, intent(in)  :: center

    real(dp), dimension(N), intent(out) :: x, y, z

    integer :: i, k
    real(dp), dimension(N) :: r
    integer, dimension(N) :: perm
    logical :: use_center

    real(dp), dimension(N,Ndim) :: samp

    ! Loop over the number of dimensions
    do k = 1, Ndim
      if (center .eqv. .True.) then
        do i = 1, N
          samp(i,k) = (i - 0.5_dp) / real(N,dp)
        end do
      else
        call random_number(r)
        do i = 1, N
          samp(i,k) = (i - 1 + r(i)) / real(N,dp)
        end do
      end if
      call random_permutation(N, perm)
      samp(:,k) = samp(perm,k)
    end do

    if (Ndim == 1) then
      x(:) = samp(:,1)
      y(:) = 0.0_dp
      z(:) = 0.0_dp
    else if  (Ndim == 2) then
      x(:) = samp(:,1)
      y(:) = samp(:,2)
      z(:) = 0.0_dp
    else
      x(:) = samp(:,1)
      y(:) = samp(:,2)
      z(:) = samp(:,3)
    end if


  end subroutine LHS_sample

  ! Fisher–Yates shuffle for 1..N
  subroutine random_permutation(N, perm)
    integer, intent(in)  :: N
    integer, intent(out) :: perm(N)
    integer :: i, j, tmp
    real(dp) :: r

    do i = 1, N
      perm(i) = i
    end do

    do i = N, 2, -1
      call random_number(r)
      j = int(r * real(i,dp)) + 1   ! j in [1, i]
      tmp     = perm(i)
      perm(i) = perm(j)
      perm(j) = tmp
    end do
  end subroutine random_permutation

end module LHS_sampling_mod
