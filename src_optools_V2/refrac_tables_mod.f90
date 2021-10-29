module refrac_tables_mod
  use optools_data_mod
  use Ray_tables_mod, only : refrac_index_calc
  implicit none


  subroutine calc_refrac_table()
    implicit none



    !! Begin openMP loops
    !$omp parallel default (none), &
    !$omp& private (l,z,s), &
    !$omp& shared (nwl,nlay,nRay,Ray_out,VMR_lay,iVMR,N_lay,RH_lay,Ray_work,wl)

    ! Perform Rayleigh scattering cross section calculation
    ! Species loops are inside subroutines
    do l = 1, nwl
      !$omp single
      if (mod(l,nwl/10) == 0) then
        print*, l, wl(l), nwl
      end if
      !$omp end single
      ! Find the refractive index and King factor for this wavelength
      call refrac_index_calc(l)
      ! Find the cross section for each species for this wavelength
      call Ray_xsec_calc(l)

      !$omp do schedule (static)
      do z = 1, nlay
        ! Find the Rayleigh cross section opacity for this layer
        A_out(z) = 0.0_dp
        do s = 1, nRay
          ! Sum species cross sections and *VMR*N_lay to convert from cm2 molecule-1 to cm-1 of atmosphere
          Ray_out(z) = Ray_out(z) + VMR_lay(iVMR(s),z) * N_lay(z) * Ray_work(s)
        end do
        ! Convert to cm2 g-1 of atmosphere
        Ray_out(z) = Ray_out(z)/RH_lay(z)
      end do
      !$omp end do

      !$omp single
      ! Output CMCRT formatted Rayleigh scattering table for layers
      call output_Ray_table()
      !$omp end single

    end do
    !$omp end parallel


  end subroutine calc_refrac_table




end module refrac_tables_mod
