module random_cpu
  use iso_fortran_env, only: int64
  use mc_precision
  implicit none
  private
  public :: rng_seed, rng_uniform, rng_next_uint64

  integer(int64) :: s(4) = [0_int64,0_int64,0_int64,0_int64]  ! xoshiro256+ state

contains

  pure function rotl(x,k) result(r)
    integer(int64), intent(in) :: x
    integer,        intent(in) :: k
    integer(int64)             :: r
    ! rotate-left using logical shifts and OR
    r = ior( ishft(x,  k), ishft(x, k-64) )
  end function rotl

  pure function splitmix64_step(z) result(znext)
    ! One step of SplitMix64 (returns next state)
    integer(int64), intent(in) :: z
    integer(int64)             :: znext
    integer(int64) :: t

    t = z + int(Z'9E3779B97F4A7C15', int64)
    t = ieor(t, ishft(t, 13));  t = t * int(Z'BF58476D1CE4E5B9', int64)
    t = ieor(t, ishft(t,  7));  t = t * int(Z'94D049BB133111EB', int64)
    t = ieor(t, ishft(t, 17))
    znext = t
  end function splitmix64_step

  subroutine rng_seed(iseed)
    ! Seed xoshiro256+ via SplitMix64 stream
    integer, intent(in) :: iseed
    integer(int64) :: z
    integer :: i

    z = int(iseed,int64)
    do i = 1, 4
      z   = splitmix64_step(z + int(i, int64))
      s(i) = z
    end do
    ! Avoid all-zero state
    if (all(s == 0_int64)) s(1) = int(Z'1', int64)
  end subroutine rng_seed

  subroutine rng_next_uint64(x)
    ! xoshiro256+ next 64-bit unsigned integer (as signed int64 container)
    integer(int64), intent(out) :: x
    integer(int64) :: t, res, s0, s1, s2, s3

    s0 = s(1); s1 = s(2); s2 = s(3); s3 = s(4)

    res = s0 + s3
    t   = ishft(s1, 17)

    s2 = ieor(s2, s0)
    s3 = ieor(s3, s1)
    s1 = ieor(s1, s2)
    s0 = ieor(s0, s3)
    s2 = ieor(s2, t)
    s3 = rotl(s3, 45)

    s(1) = s0; s(2) = s1; s(3) = s2; s(4) = s3
    x = res

  end subroutine rng_next_uint64

  subroutine rng_uniform(u)
    ! Uniform double in [0,1)
    real(dp), intent(out) :: u
    integer(int64) :: x
    integer(int64) :: hi53
    call rng_next_uint64(x)
    hi53 = ishft(x, -11)                     ! keep top 53 bits
    u = real(hi53, dp) * 2.0_dp**(-53)
  end subroutine rng_uniform

end module random_cpu

