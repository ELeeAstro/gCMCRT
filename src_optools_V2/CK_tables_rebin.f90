module CK_tables_rebin
  use optools_data_mod
  implicit none

  public :: rebin_CK_tables

contains

  subroutine rebin_CK_tables(nrebin,nG)
    implicit none

    integer, intent(in) :: nrebin, nG

    integer :: u, j, i, s, l, nnbin
    real(dp), allocatable, dimension(:) :: wn_e
    real(dp), allocatable, dimension(:) :: wl_nc, wl_ne

    real(dp), allocatable, dimension(:) :: k_work
    real(dp), allocatable, dimension(:,:) :: kg_work, g_work
    real(dp) :: k_min, k_max

    nnbin = nwl/nrebin

    print*, ' ~~ Performing CK rebinning ~~ ', nrebin
    print*, ' number of new bins', nnbin

    allocate(wn_e(nwl-1))
    allocate(wn_ne(nnbin+1),wl_nc(nnbin))

    ! We need to calculate the wavenumber grid edges of the k-table, but we are given the wavelength bin center
    ! therefore ignore the first and last bin, as they cannot be accuratly rebinned
    ! We use the first table as the wavenumber grid
    ! this fundamentally assumes that the wavelength grid is equal for all k-tables

    do l = 1, nwl-1
      wn_e(l) = (CK_tab(1)%wn(l+1) + CK_tab(1)%wn(l))/2.0_dp
    end do


    print*, ' ~~ Please wait... ~~ '

    allocate(k_work(100),kg_work(nrebin,nG))

    do s = 1, nCK

      ! Perform rebinning for each pressure and temperature for this k-table
      do j = 1, CK_tab(s)%nP
        do i = 1, CK_tab(s)%nT
          do n = 1, nnbin
            kg_work(n,:) = CK_tab(s)%k_abs(l,j,i,:) 
          end do

        end do
      end do



    end do

    print*, 'Outputting rebinned wavelength scale and changing optools wavelength grid'
    nwl = nwl - 2
    deallocate(iwl,wl,wl_cm,wl_A,wn,freq)
    allocate(iwl(nwl),wl(nwl),wl_cm(nwl),wl_A(nwl),wn(nwl),freq(nwl))


    open(newunit=u,file='wavelengths_rebin.wl',action='readwrite',status='unknown',form='formatted')
    write(u,*) nwl
    do l = 1, nwl
      write(u,*) l, wl_c(l)
    end do

    print*, ' ~~ Quest completed ~~ '


  end subroutine rebin_CK_tables

end module CK_tables_rebin
