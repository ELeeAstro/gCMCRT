! Module to calculate the capital G for the Draine scattering phase function, for parameterised alpha and g
! i.e. alpha and g are known.
! Here we take alpha from the namelist, Draine_alp = 1 corresponds to the Cornette & Shanks (1992) phase function

module mc_Draine_G
  use mc_precision
  use mc_data_mod, only : Draine_alp
  use mc_class_grid
  use cudafor
  use libm
  implicit none

  logical :: first_call = .True.

contains

  subroutine Draine_G()
    implicit none

    integer :: i, j, k
    real(dp) :: a1, b1, g, g2, g3

    if (first_call .eqv. .True.) then
      ! Need to limit alpha to non-zero small value to avoid numerical issues
      Draine_alp = max(1e-4_dp,Draine_alp)
      Draine_alp_d = Draine_alp
      allocate(Dgg(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(Dgg_d(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      first_call = .False.
    end if

    do k = 1, grid%n_theta-1
      do j = 1, grid%n_phi-1
        do i = 1, grid%n_lay
          if (Draine_alp == 0.0_dp) then
            Dgg(i,j,k) = cld_g(i,j,k)
            cycle
          end if
          if (cld_g(i,j,k) == 0.0_dp) then
            Dgg(i,j,k) = 0.0_dp
            cycle
          end if
          g = cld_g(i,j,k)
          g2 = cld_g(i,j,k)**2
          g3 = cld_g(i,j,k)**3
          a1 = 1.0_dp/2.0_dp + 5.0_dp/(6.0_dp*Draine_alp) - (25.0_dp/81.0_dp)*g2
          b1 = (125.0_dp/729.0_dp)*g3 + 5.0_dp/(9.0_dp*Draine_alp)*g
          Dgg(i,j,k) = cbrt(sqrt(a1**3 + b1**2) + b1) - cbrt(sqrt(a1**3 + b1**2) - b1) + (5.0_dp/9.0_dp)*g
        end do
      end do
    end do

    !! Draine G function
    Dgg_d(:,:,:) = Dgg(:,:,:)

  end subroutine Draine_G


end module mc_Draine_G
