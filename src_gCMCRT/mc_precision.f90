module mc_precision
  use, intrinsic :: iso_fortran_env
  implicit none

  !! Single and Double precision kinds
  integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64

  !! Single and Double precision kinds
  ! integer, parameter :: sp = kind(0.0e0)
  ! integer, parameter :: dp = kind(0.0d0)

end module mc_precision
