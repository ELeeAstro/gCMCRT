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

    !! Parameters for MIEX - set as paramaters to avoid messing up for now
    integer, parameter :: rnang = 1
    integer :: rier
    complex(kind=dp), dimension(rnang) :: rSA1, rSA2
    logical, parameter :: rdoSA = .False.


    kext = 0.0_dp
    ksca = 0.0_dp
    kscag = 0.0_dp

    ! Now to loop over species
    do s = 1, ncl
      eps = cmplx(n_work(cl_tab(s)%iVMR),k_work(cl_tab(s)%iVMR),kind=dp)
      ! Do loop over bin
      do b = 1, nmode
        a = a_C_cl_lay(cl_tab(s)%iVMR,b)
        nd = nd_C_cl_lay(cl_tab(s)%iVMR,b,z)

        if (nd < 1e-99_dp) then
          cycle
        end if

        xsec = pi * (a*1.0e-4_dp)**2

        x = (twopi * a)/wl(l)
        if (x <= 1.00001e-6_dp) then
          x = 1.00001e-6_dp
        else if (x >= 100000.0_dp) then
          x = 100000.0_dp
        end if

        !! Use the MIEX Mie theory routine
        ! - careful with memory in parallel for large x
        call shexqnn2(eps, x, rQext, rQsca, rQabs, rQbk, rQpr, ralbedo, rg, &
          & rier, rSA1, rSA2, rdoSA, rnang)

        kext = kext +  rQext * xsec * nd
        ksca = ksca + rQext * xsec * nd * ralbedo
        kscag = kscag + rQext * xsec * nd * ralbedo * rg

        if ((ieee_is_normal(kext) .eqv. .False.) .or. (kext == 0.0_dp) .or.  (rier /= 0)) then
          print*, 'cl: NaN in carma mie theory: ', l, s, b, a, nd, x, rQext, ralbedo, rg
        end if
      end do
    end do

    kext = max(kext,1e-99_dp)
    ksca = max(ksca,1e-99_dp)

    cl_out_k = kext
    cl_out_a = ksca/kext
    cl_out_g = kscag/ksca

    !print*, kext, ksca, kscag

  end subroutine calc_carma
 
end module cloud_tables_carma
