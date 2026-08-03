module LHS_sampling_mod
  use mc_precision
  use random_cpu
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
  !   - This routine does not seed the generator. Seed once in the calling
  !     program; reuse a seed for reproducibility or change it between runs.
  !   - Each axis is stratified into N bins and independently permuted.
  !--------------------------------------------------------------------
  subroutine LHS_sample(N, Ndim, x, y, z, center)
    implicit none

    integer, intent(in) :: N, Ndim
    logical, intent(in)  :: center

    real(dp), dimension(N), intent(out) :: x, y, z

    if (N < 1) then
      print*, 'ERROR: LHS sample count must be positive.'
      stop
    end if
    if ((Ndim < 1) .or. (Ndim > 3)) then
      print*, 'ERROR: LHS dimension must be between 1 and 3.'
      print*, 'Ndim: ', Ndim
      stop
    end if

    ! Generate one axis at a time.  The former N-by-Ndim automatic array and
    ! automatic permutation array could exhaust the stack at production Nph.
    call sample_axis(N, center, x)

    if (Ndim == 1) then
      y(:) = 0.0_dp
      z(:) = 0.0_dp
    else if  (Ndim == 2) then
      call sample_axis(N, center, y)
      z(:) = 0.0_dp
    else
      call sample_axis(N, center, y)
      call sample_axis(N, center, z)
    end if


  end subroutine LHS_sample

  subroutine sample_axis(N, center, axis)
    implicit none

    integer, intent(in) :: N
    logical, intent(in) :: center
    real(dp), dimension(N), intent(out) :: axis

    integer :: i
    integer, allocatable, dimension(:) :: perm
    real(dp), allocatable, dimension(:) :: strata
    real(dp) :: r

    allocate(perm(N), strata(N))

    if (center .eqv. .True.) then
      do i = 1, N
        strata(i) = (real(i,dp) - 0.5_dp) / real(N,dp)
      end do
    else
      do i = 1, N
        call rng_uniform(r)
        strata(i) = (real(i-1,dp) + r) / real(N,dp)
      end do
    end if

    call random_permutation(N, perm)
    do i = 1, N
      axis(i) = strata(perm(i))
    end do

    deallocate(perm, strata)

  end subroutine sample_axis

  ! Fisher–Yates shuffle for 1..N
  subroutine random_permutation(N, perm)
    implicit none

    integer, intent(in)  :: N

    integer, dimension(N), intent(out) :: perm

    integer :: i, j, tmp
    real(dp) :: r

    do i = 1, N
      perm(i) = i
    end do

    do i = N, 2, -1
      !call random_number(r)
      call rng_uniform(r)
      j = int(r * real(i,dp)) + 1   ! j in [1, i]
      tmp     = perm(i)
      perm(i) = perm(j)
      perm(j) = tmp
    end do
  end subroutine random_permutation

end module LHS_sampling_mod
