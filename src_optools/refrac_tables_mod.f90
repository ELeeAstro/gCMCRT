module refrac_tables_mod
  use optools_data_mod
  use gas_refractivity_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none

  private
  public :: calc_refrac_table

contains

  subroutine calc_refrac_table()
    implicit none

    integer :: i, j, l, s, z, reclen
    integer, allocatable :: iVMR(:), model(:)
    real(dp), allocatable :: n_reference(:), density_reference(:)
    real(dp), allocatable :: alpha_volume(:), refrac_out(:)
    real(sp), allocatable :: refrac_write(:)
    real(dp) :: species_density, species_term, mixture_term
    real(dp) :: n_squared, refractivity
    logical :: found, water_valid

    if (.not. refrac) return

    if (nrefrac < 1) then
      print*, 'ERROR - Refractivity requested without any refractive species'
      stop
    end if

    allocate(iVMR(nrefrac), model(nrefrac))
    allocate(n_reference(nrefrac), density_reference(nrefrac))
    allocate(alpha_volume(nrefrac))
    allocate(refrac_out(nlay), refrac_write(nlay))

    ! Resolve the requested species once and reject duplicates explicitly.
    do s = 1, nrefrac
      do i = 1, s-1
        if (trim(refrac_name(s)) == trim(refrac_name(i))) then
          print*, 'ERROR - Duplicate refractive species: ', trim(refrac_name(s))
          stop
        end if
      end do

      found = .False.
      do j = 1, ngas
        if (trim(refrac_name(s)) == trim(g_name(j))) then
          iVMR(s) = j
          found = .True.
          exit
        end if
      end do
      if (.not. found) then
        print*, 'ERROR - Refractive species is absent from the profile: ', &
          trim(refrac_name(s))
        stop
      end if
    end do

    print*, ' ~~ Performing gas refractivity calculation and output ~~ '
    print*, ' ~~ Output quantity is nu = n - 1 ~~ '

    inquire(iolength=reclen) refrac_write
    open(newunit=urefrac, file='refrac.cmcrt', action='readwrite', &
      form='unformatted', status='replace', access='direct', recl=reclen)

    do l = 1, nwl
      if (wl(l) <= 0.0_dp .or. .not. ieee_is_finite(wl(l))) then
        print*, 'ERROR - Invalid wavelength for refractivity at index: ', l, wl(l)
        stop
      end if

      ! Gather wavelength-dependent species data outside the layer loop.
      do s = 1, nrefrac
        call gas_refractivity_data(refrac_name(s), wl(l), wn(l), model(s), &
          n_reference(s), density_reference(s), alpha_volume(s))
        if (model(s) == REFRAC_MODEL_UNAVAILABLE) then
          print*, 'ERROR - No neutral-gas refractivity model for species: ', &
            trim(refrac_name(s))
          print*, 'Atomic H and free electrons require separate dispersion models.'
          stop
        end if
      end do

      do z = 1, nlay
        if (N_lay(z) < 0.0_dp .or. .not. ieee_is_finite(N_lay(z))) then
          print*, 'ERROR - Invalid total number density in layer: ', z, N_lay(z)
          stop
        end if
        if (TG_lay(z) <= 0.0_dp .or. .not. ieee_is_finite(TG_lay(z))) then
          print*, 'ERROR - Invalid temperature in layer: ', z, TG_lay(z)
          stop
        end if

        mixture_term = 0.0_dp
        do s = 1, nrefrac
          if (VMR_lay(iVMR(s),z) < 0.0_dp .or. &
              .not. ieee_is_finite(VMR_lay(iVMR(s),z))) then
            print*, 'ERROR - Invalid VMR for refractive species/layer: ', &
              trim(refrac_name(s)), z, VMR_lay(iVMR(s),z)
            stop
          end if

          species_density = VMR_lay(iVMR(s),z)*N_lay(z)

          select case(model(s))
          case(REFRAC_MODEL_REFERENCE_INDEX)
            if (n_reference(s) <= 0.0_dp .or. &
                density_reference(s) <= 0.0_dp .or. &
                .not. ieee_is_finite(n_reference(s))) then
              print*, 'ERROR - Invalid reference refractivity data: ', &
                trim(refrac_name(s)), l
              stop
            end if
            species_term = species_density/density_reference(s) * &
              ((n_reference(s)**2-1.0_dp)/(n_reference(s)**2+2.0_dp))

          case(REFRAC_MODEL_POLARIZABILITY)
            species_term = (4.0_dp*pi/3.0_dp)*species_density*alpha_volume(s)

          case(REFRAC_MODEL_WATER)
            call water_lorentz_lorenz(wl(l), TG_lay(z), species_density, &
              species_term, water_valid)
            if (.not. water_valid) then
              print*, 'ERROR - H2O refractivity model is invalid at wavelength/layer: ', &
                wl(l), z
              print*, 'The current H2O model requires 0.2 < wavelength < 2.5 micron.'
              stop
            end if

          case default
            print*, 'ERROR - Invalid refractivity model for species: ', &
              trim(refrac_name(s))
            stop
          end select

          if (species_term < 0.0_dp .or. .not. ieee_is_finite(species_term)) then
            print*, 'ERROR - Invalid Lorentz-Lorenz contribution: ', &
              trim(refrac_name(s)), l, z, species_term
            stop
          end if
          mixture_term = mixture_term + species_term
        end do

        if (mixture_term < 0.0_dp .or. mixture_term >= 1.0_dp .or. &
            .not. ieee_is_finite(mixture_term)) then
          print*, 'ERROR - Invalid mixture Lorentz-Lorenz term: ', &
            l, z, mixture_term
          stop
        end if

        n_squared = (1.0_dp+2.0_dp*mixture_term)/(1.0_dp-mixture_term)
        ! Stable evaluation of sqrt(n_squared)-1 for weak upper atmospheres.
        refractivity = 3.0_dp*mixture_term / &
          ((1.0_dp-mixture_term)*(sqrt(n_squared)+1.0_dp))

        if (refractivity < 0.0_dp .or. .not. ieee_is_finite(refractivity)) then
          print*, 'ERROR - Invalid mixture refractivity: ', l, z, refractivity
          stop
        end if
        refrac_out(z) = refractivity
      end do

      refrac_write(:) = real(max(refrac_out(:),0.0_dp),sp)
      write(urefrac,rec=l) refrac_write

      if (mod(l,max(1,nwl/10)) == 0 .or. l == 1 .or. l == nwl) then
        print*, l, wl(l), minval(refrac_out), maxval(refrac_out)
      end if
    end do

    close(urefrac)
    deallocate(iVMR, model, n_reference, density_reference, alpha_volume)
    deallocate(refrac_out, refrac_write)

    print*, ' ~~ Refractivity table complete ~~ '

  end subroutine calc_refrac_table

end module refrac_tables_mod
