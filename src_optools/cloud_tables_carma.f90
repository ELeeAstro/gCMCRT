module cloud_tables_carma
  use optools_data_mod
  use mie_routines, only : shexqnn2
  use ieee_arithmetic
  implicit none

contains

  subroutine calc_carma(z,l,n_work,k_work, cl_out_k, cl_out_a, cl_out_g)
    implicit none

    integer, intent(in) :: z, l
    real(dp), dimension(ncl),  intent(in) :: n_work, k_work

    real(kind=dp), intent(out) :: cl_out_k, cl_out_a, cl_out_g

    integer :: s, b 
    real(dp) :: a, nd, kext, ksca, kscag

    complex(kind=dp) :: eps
    real(kind=dp) :: rQext, rQsca, rQabs, rQbk, rQpr, ralbedo, rg
    real(kind=dp) :: x, xsec, xsec_ext, xsec_sca

    !! Parameters for MIEX; keep them fixed here for now.
    integer, parameter :: rnang = 1
    integer :: rier
    complex(kind=dp), dimension(rnang) :: rSA1, rSA2
    logical, parameter :: rdoSA = .False.


    kext = 0.0_dp
    ksca = 0.0_dp
    kscag = 0.0_dp

    ! Now to loop over species
    do s = 1, ncl
      eps = cmplx(n_work(s),k_work(s),kind=dp)
      ! Do loop over bin
      do b = 1, nmode
        a = a_C_cl_lay(cl_tab(s)%iVMR,b)
        nd = nd_C_cl_lay(cl_tab(s)%iVMR,b,z)

        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(nd) .or. &
          & .not. ieee_is_finite(real(eps,kind=dp)) .or. &
          & .not. ieee_is_finite(aimag(eps)) .or. a <= 0.0_dp .or. nd < 0.0_dp) then
          print*, 'ERROR - Invalid CARMA cloud input passed to Mie theory - STOPPING'
          print*, 'Layer, wavelength index, wavelength [um]: ', z, l, wl(l)
          print*, 'Species, size bin, radius [cm], number density: ', &
            & cl_tab(s)%sp, b, a, nd
          print*, 'Refractive index: ', eps
          stop
        end if

        if (nd < 1e-99_dp) then
          cycle
        end if

        xsec = pi * a**2

        x = (twopi * a*1e4_dp)/wl(l)
        if (x <= 1.00001e-6_dp) then
          x = 1.00001e-6_dp
        else if (x >= 100000.0_dp) then
          x = 100000.0_dp
        end if

        !! Use the MIEX Mie theory routine
        ! Take care with per-thread memory use when x is large.
        call shexqnn2(eps, x, rQext, rQsca, rQabs, rQbk, rQpr, ralbedo, rg, &
          & rier, rSA1, rSA2, rdoSA, rnang)

        if (rier /= 0 .or. .not. ieee_is_finite(rQext) .or. &
          & .not. ieee_is_finite(rQsca) .or. rQext < 0.0_dp .or. rQsca < 0.0_dp .or. &
          & (rQsca > 0.0_dp .and. .not. ieee_is_finite(rg))) then
          print*, 'ERROR - Corrupted CARMA Mie output - STOPPING'
          print*, 'Layer, wavelength index, wavelength [um]: ', z, l, wl(l)
          print*, 'Species, size bin, radius [cm], number density: ', &
            & cl_tab(s)%sp, b, a, nd
          print*, 'Refractive index, size parameter, status: ', eps, x, rier
          print*, 'Qext, Qsca, g: ', rQext, rQsca, rg
          stop
        end if

        if (rQext <= 0.0_dp .or. rQsca <= 0.0_dp) then
          ralbedo = 0.0_dp
          rg = 0.0_dp
        else
          ralbedo = min(1.0_dp,max(0.0_dp,rQsca/rQext))
          rg = min(1.0_dp,max(-1.0_dp,rg))
        end if

        kext = kext +  rQext * xsec * nd
        ksca = ksca + rQext * xsec * nd * ralbedo
        kscag = kscag + rQext * xsec * nd * ralbedo * rg
      end do
    end do

    if (.not. ieee_is_finite(kext) .or. .not. ieee_is_finite(ksca) .or. &
      & .not. ieee_is_finite(kscag) .or. kext < 0.0_dp .or. ksca < 0.0_dp) then
      print*, 'ERROR - Corrupted integrated CARMA cloud output - STOPPING'
      print*, 'Layer, wavelength index, wavelength [um]: ', z, l, wl(l)
      print*, 'Extinction, scattering, scattering-g: ', kext, ksca, kscag
      stop
    end if

    cl_out_k = max(kext,1.0e-199_dp)
    if (kext <= 0.0_dp .or. ksca <= 0.0_dp) then
      cl_out_a = 0.0_dp
      cl_out_g = 0.0_dp
    else
      cl_out_a = min(1.0_dp,max(0.0_dp,ksca/kext))
      cl_out_g = min(1.0_dp,max(-1.0_dp,kscag/ksca))
    end if

    !print*, kext, ksca, kscag

  end subroutine calc_carma
 
end module cloud_tables_carma
