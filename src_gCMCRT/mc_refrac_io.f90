module mc_refrac_io
  use mc_precision
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none

  private

  integer, parameter, public :: REFRAC_IO_SUCCESS = 0
  integer, parameter, public :: REFRAC_IO_BAD_DIMENSIONS = 1
  integer, parameter, public :: REFRAC_IO_BAD_GEOMETRY = 2
  integer, parameter, public :: REFRAC_IO_BAD_VALUE = 3

  public :: unpack_refrac_record

contains

  subroutine unpack_refrac_record(record, oneD, threeD, refrac_grid, status, &
      bad_index)
    implicit none

    real(sp), intent(in) :: record(:)
    logical, intent(in) :: oneD, threeD
    real(dp), intent(out) :: refrac_grid(:,:,:)
    integer, intent(out) :: status, bad_index

    integer :: n_lay, n_lon, n_lat
    integer :: z, j, k, n, expected_size

    status = REFRAC_IO_SUCCESS
    bad_index = 0
    refrac_grid(:,:,:) = 0.0_dp

    if (oneD .eqv. threeD) then
      status = REFRAC_IO_BAD_GEOMETRY
      return
    end if

    n_lay = size(refrac_grid,1)
    n_lon = size(refrac_grid,2)
    n_lat = size(refrac_grid,3)
    if (n_lay < 1 .or. n_lon < 1 .or. n_lat < 1) then
      status = REFRAC_IO_BAD_DIMENSIONS
      return
    end if

    if (oneD) then
      expected_size = n_lay
    else
      expected_size = n_lay*n_lon*n_lat
    end if
    if (size(record) /= expected_size) then
      status = REFRAC_IO_BAD_DIMENSIONS
      return
    end if

    do n = 1, size(record)
      if (.not. ieee_is_finite(real(record(n),dp)) .or. &
          real(record(n),dp) < 0.0_dp .or. real(record(n),dp) >= 1.0_dp) then
        status = REFRAC_IO_BAD_VALUE
        bad_index = n
        return
      end if
    end do

    if (oneD) then
      do z = 1, n_lay
        refrac_grid(z,:,:) = real(record(z),dp)
      end do
    else
      n = 1
      do k = 1, n_lat
        do j = 1, n_lon
          do z = 1, n_lay
            refrac_grid(z,j,k) = real(record(n),dp)
            n = n + 1
          end do
        end do
      end do
    end if

  end subroutine unpack_refrac_record

end module mc_refrac_io
