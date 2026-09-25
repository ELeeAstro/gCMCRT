module lbl_tables_mod
  use optools_data_mod
  use lbl_tables_read
  use lbl_tables_interp
  use lbl_tables_combine
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
#ifdef GPU_OFFLOAD
  use omp_lib, only : omp_get_num_devices, omp_is_initial_device
  use, intrinsic :: iso_fortran_env, only : int64
#endif
  implicit none

#ifdef GPU_OFFLOAD
#ifndef LBL_GPU_BLOCK
#define LBL_GPU_BLOCK 64
#endif
#endif

  logical :: first_call = .True.

  real(kind=dp), allocatable, dimension(:) :: lbl_out
  real(kind=sp), allocatable, dimension(:) :: lbl_write

  ! Namelist options
  integer :: iopts
  integer, allocatable, dimension(:) :: form
  character(len=150), allocatable, dimension(:) :: paths
  logical :: interp_wl

  namelist /lbl_nml/ iopts, form, paths, interp_wl

  private :: output_lbl_table, validate_lbl_tables, calc_lbl_table_cpu
#ifdef GPU_OFFLOAD
  private :: calc_lbl_table_gpu, locate_triplet_gpu, Bezier_interp_gpu
  !$omp declare target(locate_triplet_gpu,Bezier_interp_gpu)
#endif
  public :: calc_lbl_table


contains

  subroutine calc_lbl_table()
    implicit none

    integer :: s, j
    logical :: exists

    ! The array-allocation order is important because it follows the namelist order.

    ! Allocate number of lbl tables
    allocate(lbl_tab(nlbl))

    ! Allocate required number of arrays from namelist options
    allocate(form(nlbl),paths(nlbl))

    ! Read lbl namelist parameters
    read(u_nml, nml=lbl_nml)

    if (interp_wl .eqv. .True.) then
      print*, 'ERROR - LBL wavelength interpolation is not supported - STOPPING'
      print*, 'Set interp_wl = .False. and use the calculation wavelength grid in the LBL tables.'
      stop
    end if

    ! Allocate work arrays
    allocate(lbl_out(nlay),lbl_write(nlay))

    ! Give the classes some global data from par and namelists
    lbl_tab(:)%sp = lbl_name(:)
    lbl_tab(:)%form = form(:)
    lbl_tab(:)%path = paths(:)

    ! Find the PRF VMR indices of the LBL species.
    do s = 1, nlbl
      exists = .False.
      do j = 1, ngas
        if (lbl_tab(s)%sp == g_name(j)) then
          lbl_tab(s)%iVMR = j
          exists = .True.
          exit
        end if
      end do
      if (exists .eqv. .False.) then
        print*, 'ERROR - Specified lbl species not found in prf VMR list - STOPPING'
        print*, 'Species: ', lbl_tab(s)%sp
        stop
      end if
    end do

    ! Read the lbl tables
    call read_lbl_tables()

    call validate_lbl_tables()

#ifdef GPU_OFFLOAD
    call calc_lbl_table_gpu()
#else
    call calc_lbl_table_cpu()
#endif

    print*, ' ~~ Quest completed ~~ '

    !deallocate all allocated arrays
    deallocate(lbl_tab)
    deallocate(form,paths)
    deallocate(lbl_out,lbl_write)
    ! Close the LBL I/O unit.
    close(ulbl)

  end subroutine calc_lbl_table

  subroutine calc_lbl_table_cpu()
    implicit none

    integer :: l, z
    real(kind=dp), allocatable, dimension(:) :: lbl_work
    real(kind=dp) :: lbl_comb

    allocate(lbl_work(nlbl))
    lbl_work(:) = 0.0_dp
    lbl_comb = 0.0_dp

    print*, ' ~~ Performing lbl interpolation, combining and output ~~ '
    print*, ' ~~ Please wait... ~~ '

    !$omp parallel default (none), &
    !$omp& private (l,z), &
    !$omp& shared (nwl,wl,nlay,lbl_out,RH_lay,interp_wl), &
    !$omp& firstprivate(lbl_work, lbl_comb)

    do l = 1, nwl
      !$omp single
      if (mod(l,max(1,nwl/10)) == 0) then
        print*, l, wl(l), nwl
      end if
      !$omp end single

      !$omp do schedule (dynamic)
      do z = 1, nlay
        call interp_lbl_tables_Bezier(l,z,lbl_work(:))
        call combine_lbl_opacity(z,lbl_work(:),lbl_comb)
        lbl_out(z) = lbl_comb/RH_lay(z)
      end do
      !$omp end do

      !$omp single
      call output_lbl_table(l)
      !$omp end single
    end do
    !$omp end parallel

    deallocate(lbl_work)

  end subroutine calc_lbl_table_cpu

#ifdef GPU_OFFLOAD
  pure subroutine locate_triplet_gpu(arr,offset,n,var,idx,region,exact_idx)
    implicit none

    real(kind=dp), dimension(*), intent(in) :: arr
    integer(kind=int64), intent(in) :: offset
    integer, intent(in) :: n
    real(kind=dp), intent(in) :: var
    integer, dimension(3), intent(out) :: idx
    integer, intent(out) :: region, exact_idx
    integer :: lower, upper, middle, start

    lower = 0
    upper = n + 1
    do while (upper-lower > 1)
      middle = (upper+lower)/2
      if (var > arr(offset+int(middle-1,kind=int64))) then
        lower = middle
      else
        upper = middle
      end if
    end do

    exact_idx = 0
    if (var < arr(offset)) then
      region = -1
    else if (var > arr(offset+int(n-1,kind=int64))) then
      region = 1
    else
      region = 0
      if (lower < n) then
        if (var == arr(offset+int(lower,kind=int64))) exact_idx = lower + 1
      end if
    end if

    start = max(1,min(lower-1,n-2))
    idx(1) = start
    idx(2) = start + 1
    idx(3) = start + 2

  end subroutine locate_triplet_gpu

  pure subroutine Bezier_interp_gpu(xi,yi,x,y)
    implicit none

    real(kind=dp), dimension(3), intent(in) :: xi, yi
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(out) :: y
    real(kind=dp) :: dx, dx1, dy, dy1, w, yc, t, wlim, wlim1
    real(kind=dp) :: denom, denom1, grad_scale, grad_tol, y_linear

    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if (x <= xi(1)) then
      y = yi(1)
      return
    else if (x >= xi(3)) then
      y = yi(3)
      return
    else if (x == xi(2)) then
      y = yi(2)
      return
    end if

    if (x < xi(2)) then
      t = (x-xi(1))/dx
      y_linear = (1.0_dp-t)*yi(1) + t*yi(2)
    else
      t = (x-xi(2))/dx1
      y_linear = (1.0_dp-t)*yi(2) + t*yi(3)
    end if

    grad_scale = max(1.0_dp,abs(yi(1)),abs(yi(2)),abs(yi(3)))
    grad_tol = 100.0_dp*epsilon(1.0_dp)*grad_scale
    if (abs(dy) <= grad_tol .or. abs(dy1) <= grad_tol .or. &
      & (dy > 0.0_dp .and. dy1 < 0.0_dp) .or. &
      & (dy < 0.0_dp .and. dy1 > 0.0_dp)) then
      y = y_linear
      return
    end if

    denom = 1.0_dp - (dy1/dy)*(dx/dx1)
    denom1 = 1.0_dp - (dy/dy1)*(dx1/dx)
    if (denom /= denom .or. denom1 /= denom1 .or. &
      & abs(denom) > huge(denom) .or. abs(denom1) > huge(denom1) .or. &
      & abs(denom) <= 100.0_dp*epsilon(1.0_dp) .or. &
      & abs(denom1) <= 100.0_dp*epsilon(1.0_dp)) then
      y = y_linear
      return
    end if

    if (x < xi(2)) then
      w = dx1/(dx+dx1)
      wlim = 1.0_dp + 1.0_dp/denom
      wlim1 = 1.0_dp/denom1
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) w = 1.0_dp
      yc = yi(2) - dx/2.0_dp*(w*dy/dx+(1.0_dp-w)*dy1/dx1)
      t = (x-xi(1))/dx
      if (yc /= yc .or. abs(yc) > huge(yc) .or. yc < min(yi(1),yi(2)) .or. &
        & yc > max(yi(1),yi(2))) then
        y = y_linear
        return
      end if
      y = (1.0_dp-t)**2*yi(1) + 2.0_dp*t*(1.0_dp-t)*yc + t**2*yi(2)
    else
      w = dx/(dx+dx1)
      wlim = 1.0_dp/denom
      wlim1 = 1.0_dp + 1.0_dp/denom1
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) w = 1.0_dp
      yc = yi(2) + dx1/2.0_dp*(w*dy1/dx1+(1.0_dp-w)*dy/dx)
      t = (x-xi(2))/dx1
      if (yc /= yc .or. abs(yc) > huge(yc) .or. yc < min(yi(2),yi(3)) .or. &
        & yc > max(yi(2),yi(3))) then
        y = y_linear
        return
      end if
      y = (1.0_dp-t)**2*yi(2) + 2.0_dp*t*(1.0_dp-t)*yc + t**2*yi(3)
    end if

    if (y /= y .or. abs(y) > huge(y)) y = y_linear

  end subroutine Bezier_interp_gpu

  subroutine calc_lbl_table_gpu()
    implicit none

    integer, parameter :: block_size = LBL_GPU_BLOCK
    integer :: s, itab, ptab, l, l_start, l_end, n_block, lb, z, j
    integer :: T_region, P_region, T_exact, P_exact, T_fixed, P_fixed
    integer :: iT_idx(3), iP_idx(3)
    integer :: iT1, iT2, iT3, iP1, iP2, iP3
    integer :: offload_active
    integer, allocatable :: lbl_nT(:), lbl_nP(:), lbl_iVMR(:)
    integer(kind=int64) :: total_T, total_P, total_cells, T_pos, P_pos
    integer(kind=int64) :: block_values, block_pos, k_idx
    integer(kind=int64), allocatable :: lbl_T_offset(:), lbl_P_offset(:)
    integer(kind=int64), allocatable :: lbl_K_offset(:)
    real(kind=dp), allocatable :: lbl_T_flat(:), lbl_lT_flat(:)
    real(kind=dp), allocatable :: lbl_P_flat(:), lbl_lP_flat(:)
    real(kind=dp), allocatable :: lbl_K_block(:), lbl_out_block(:,:)
    real(kind=dp) :: T_layer, P_layer, lT_layer, lP_layer
    real(kind=dp) :: log_k, k_value, lbl_comb
    real(kind=dp) :: lTa(3), lPa(3), lka(3), lka_lbl(3)

    if (block_size < 1) then
      print*, 'ERROR - GPU LBL wavelength block must be positive - STOPPING'
      stop 1
    end if

    if (any(TG_lay <= 0.0_dp) .or. any(.not. ieee_is_finite(TG_lay))) then
      print*, 'ERROR - Invalid atmospheric temperature for GPU LBL interpolation - STOPPING'
      stop 1
    end if
    if (any(PG_lay <= 0.0_dp) .or. any(.not. ieee_is_finite(PG_lay))) then
      print*, 'ERROR - Invalid atmospheric pressure for GPU LBL interpolation - STOPPING'
      stop 1
    end if
    if (any(RH_lay <= 0.0_dp) .or. any(.not. ieee_is_finite(RH_lay)) .or. &
      & any(N_lay < 0.0_dp) .or. any(.not. ieee_is_finite(N_lay))) then
      print*, 'ERROR - Invalid atmospheric density for GPU LBL interpolation - STOPPING'
      stop 1
    end if
    do s = 1, nlbl
      if (any(VMR_lay(lbl_tab(s)%iVMR,:) < 0.0_dp) .or. &
        & any(.not. ieee_is_finite(VMR_lay(lbl_tab(s)%iVMR,:)))) then
        print*, 'ERROR - Invalid LBL species VMR for GPU interpolation: ', lbl_tab(s)%sp
        stop 1
      end if
    end do

    allocate(lbl_nT(nlbl),lbl_nP(nlbl),lbl_iVMR(nlbl))
    allocate(lbl_T_offset(nlbl),lbl_P_offset(nlbl),lbl_K_offset(nlbl))

    total_T = 0_int64
    total_P = 0_int64
    total_cells = 0_int64
    do s = 1, nlbl
      lbl_nT(s) = lbl_tab(s)%nT
      lbl_nP(s) = lbl_tab(s)%nP
      lbl_iVMR(s) = lbl_tab(s)%iVMR
      lbl_T_offset(s) = total_T + 1_int64
      lbl_P_offset(s) = total_P + 1_int64
      lbl_K_offset(s) = total_cells*int(block_size,kind=int64) + 1_int64
      total_T = total_T + int(lbl_nT(s),kind=int64)
      total_P = total_P + int(lbl_nP(s),kind=int64)
      total_cells = total_cells + &
        & int(lbl_nT(s),kind=int64)*int(lbl_nP(s),kind=int64)
    end do

    block_values = total_cells*int(block_size,kind=int64)
    allocate(lbl_T_flat(total_T),lbl_lT_flat(total_T))
    allocate(lbl_P_flat(total_P),lbl_lP_flat(total_P))
    allocate(lbl_K_block(block_values),lbl_out_block(nlay,block_size))
    lbl_K_block(:) = 0.0_dp
    lbl_out_block(:,:) = 0.0_dp

    do s = 1, nlbl
      T_pos = lbl_T_offset(s)
      P_pos = lbl_P_offset(s)
      lbl_T_flat(T_pos:T_pos+int(lbl_nT(s)-1,kind=int64)) = lbl_tab(s)%T(:)
      lbl_lT_flat(T_pos:T_pos+int(lbl_nT(s)-1,kind=int64)) = lbl_tab(s)%lT(:)
      lbl_P_flat(P_pos:P_pos+int(lbl_nP(s)-1,kind=int64)) = lbl_tab(s)%P(:)
      lbl_lP_flat(P_pos:P_pos+int(lbl_nP(s)-1,kind=int64)) = lbl_tab(s)%lP(:)
    end do

    offload_active = 0
    !$omp target map(from: offload_active)
    if (.not. omp_is_initial_device()) offload_active = 1
    !$omp end target
    if (offload_active /= 1) then
      print*, 'ERROR - LBL OpenMP target region did not execute on a GPU - STOPPING'
      print*, 'OpenMP target devices visible: ', omp_get_num_devices()
      stop 1
    end if

    print*, ' ~~ LBL OpenMP GPU offload active; visible devices: ', omp_get_num_devices()
    print*, ' ~~ LBL GPU wavelength block size: ', block_size
    print*, ' ~~ LBL GPU block storage [MiB]: ', &
      & real(block_values,kind=dp)*real(storage_size(1.0_dp)/8,kind=dp)/(1024.0_dp**2)
    print*, ' ~~ Performing lbl interpolation, combining and output ~~ '
    print*, ' ~~ Please wait... ~~ '

    !$omp target data map(to: TG_lay,PG_lay,N_lay,RH_lay,VMR_lay,nlay,nlbl) &
    !$omp& map(to: lbl_nT,lbl_nP,lbl_iVMR,lbl_T_offset,lbl_P_offset,lbl_K_offset) &
    !$omp& map(to: lbl_T_flat,lbl_lT_flat,lbl_P_flat,lbl_lP_flat) &
    !$omp& map(alloc: lbl_K_block,lbl_out_block)

    do l_start = 1, nwl, block_size
      l_end = min(nwl,l_start+block_size-1)
      n_block = l_end-l_start+1

      do s = 1, nlbl
        do itab = 1, lbl_nT(s)
          do ptab = 1, lbl_nP(s)
            block_pos = lbl_K_offset(s) &
              & + int(itab-1,kind=int64)*int(lbl_nP(s),kind=int64) &
              & * int(block_size,kind=int64) &
              & + int(ptab-1,kind=int64)*int(block_size,kind=int64)
            lbl_K_block(block_pos:block_pos+int(n_block-1,kind=int64)) = &
              & lbl_tab(s)%lk_abs(l_start:l_end,ptab,itab)
          end do
        end do
      end do
      !$omp target update to(lbl_K_block)

      !$omp target teams loop collapse(2) &
      !$omp& private(s,j,T_region,P_region,T_exact,P_exact,T_fixed,P_fixed) &
      !$omp& private(iT_idx,iP_idx,iT1,iT2,iT3,iP1,iP2,iP3,k_idx) &
      !$omp& private(T_layer,P_layer,lT_layer,lP_layer,log_k,k_value,lbl_comb) &
      !$omp& private(lTa,lPa,lka,lka_lbl)
      do lb = 1, n_block
        do z = 1, nlay
          T_layer = TG_lay(z)
          P_layer = PG_lay(z)
          lT_layer = log10(T_layer)
          lP_layer = log10(P_layer)
          lbl_comb = 0.0_dp

          do s = 1, nlbl
            call locate_triplet_gpu(lbl_T_flat,lbl_T_offset(s),lbl_nT(s), &
              & T_layer,iT_idx,T_region,T_exact)
            call locate_triplet_gpu(lbl_P_flat,lbl_P_offset(s),lbl_nP(s), &
              & P_layer,iP_idx,P_region,P_exact)

            iT1 = iT_idx(1)
            iT2 = iT_idx(2)
            iT3 = iT_idx(3)
            iP1 = iP_idx(1)
            iP2 = iP_idx(2)
            iP3 = iP_idx(3)

            lTa(1) = lbl_lT_flat(lbl_T_offset(s)+int(iT1-1,kind=int64))
            lTa(2) = lbl_lT_flat(lbl_T_offset(s)+int(iT2-1,kind=int64))
            lTa(3) = lbl_lT_flat(lbl_T_offset(s)+int(iT3-1,kind=int64))
            lPa(1) = lbl_lP_flat(lbl_P_offset(s)+int(iP1-1,kind=int64))
            lPa(2) = lbl_lP_flat(lbl_P_offset(s)+int(iP2-1,kind=int64))
            lPa(3) = lbl_lP_flat(lbl_P_offset(s)+int(iP3-1,kind=int64))

            T_fixed = T_exact
            if (T_region == -1) T_fixed = 1
            if (T_region == 1) T_fixed = lbl_nT(s)
            P_fixed = P_exact
            if (P_region == -1) P_fixed = 1
            if (P_region == 1) P_fixed = lbl_nP(s)

            if (T_fixed > 0 .and. P_fixed > 0) then
              k_idx = lbl_K_offset(s) &
                & + int(T_fixed-1,kind=int64)*int(lbl_nP(s),kind=int64) &
                & * int(block_size,kind=int64) &
                & + int(P_fixed-1,kind=int64)*int(block_size,kind=int64) &
                & + int(lb-1,kind=int64)
              log_k = lbl_K_block(k_idx)

            else if (T_fixed > 0) then
              do j = 1, 3
                k_idx = lbl_K_offset(s) &
                  & + int(T_fixed-1,kind=int64)*int(lbl_nP(s),kind=int64) &
                  & * int(block_size,kind=int64) &
                  & + int(iP_idx(j)-1,kind=int64)*int(block_size,kind=int64) &
                  & + int(lb-1,kind=int64)
                lka(j) = lbl_K_block(k_idx)
              end do
              call Bezier_interp_gpu(lPa,lka,lP_layer,log_k)

            else if (P_fixed > 0) then
              do j = 1, 3
                k_idx = lbl_K_offset(s) &
                  & + int(iT_idx(j)-1,kind=int64)*int(lbl_nP(s),kind=int64) &
                  & * int(block_size,kind=int64) &
                  & + int(P_fixed-1,kind=int64)*int(block_size,kind=int64) &
                  & + int(lb-1,kind=int64)
                lka(j) = lbl_K_block(k_idx)
              end do
              call Bezier_interp_gpu(lTa,lka,lT_layer,log_k)

            else
              do j = 1, 3
                k_idx = lbl_K_offset(s) &
                  & + int(iT_idx(j)-1,kind=int64)*int(lbl_nP(s),kind=int64) &
                  & * int(block_size,kind=int64) &
                  & + int(iP1-1,kind=int64)*int(block_size,kind=int64) &
                  & + int(lb-1,kind=int64)
                lka(1) = lbl_K_block(k_idx)
                k_idx = lbl_K_offset(s) &
                  & + int(iT_idx(j)-1,kind=int64)*int(lbl_nP(s),kind=int64) &
                  & * int(block_size,kind=int64) &
                  & + int(iP2-1,kind=int64)*int(block_size,kind=int64) &
                  & + int(lb-1,kind=int64)
                lka(2) = lbl_K_block(k_idx)
                k_idx = lbl_K_offset(s) &
                  & + int(iT_idx(j)-1,kind=int64)*int(lbl_nP(s),kind=int64) &
                  & * int(block_size,kind=int64) &
                  & + int(iP3-1,kind=int64)*int(block_size,kind=int64) &
                  & + int(lb-1,kind=int64)
                lka(3) = lbl_K_block(k_idx)
                call Bezier_interp_gpu(lPa,lka,lP_layer,lka_lbl(j))
              end do
              call Bezier_interp_gpu(lTa,lka_lbl,lT_layer,log_k)
            end if

            k_value = 10.0_dp**log_k
            lbl_comb = lbl_comb + k_value*N_lay(z)*VMR_lay(lbl_iVMR(s),z)
          end do

          lbl_out_block(z,lb) = lbl_comb/RH_lay(z)
        end do
      end do
      !$omp end target teams loop

      !$omp target update from(lbl_out_block)
      if (any(.not. ieee_is_finite(lbl_out_block(:,1:n_block)))) then
        print*, 'ERROR - Non-finite LBL opacity returned by GPU in wavelength block: ', &
          & l_start, l_end
        stop 1
      end if

      do lb = 1, n_block
        l = l_start+lb-1
        lbl_out(:) = lbl_out_block(:,lb)
        call output_lbl_table(l)
        if (mod(l,max(1,nwl/10)) == 0) print*, l, wl(l), nwl
      end do
    end do

    !$omp end target data

    deallocate(lbl_nT,lbl_nP,lbl_iVMR)
    deallocate(lbl_T_offset,lbl_P_offset,lbl_K_offset)
    deallocate(lbl_T_flat,lbl_lT_flat,lbl_P_flat,lbl_lP_flat)
    deallocate(lbl_K_block,lbl_out_block)

  end subroutine calc_lbl_table_gpu
#endif

  subroutine validate_lbl_tables()
    implicit none

    integer :: s, i
    real(kind=dp), parameter :: range_tol = 1.0e-2_dp
    real(kind=dp) :: scale

    do s = 1, nlbl
      if (.not. allocated(lbl_tab(s)%T) .or. .not. allocated(lbl_tab(s)%lT) .or. &
        & .not. allocated(lbl_tab(s)%P) .or. .not. allocated(lbl_tab(s)%lP) .or. &
        & .not. allocated(lbl_tab(s)%wl) .or. .not. allocated(lbl_tab(s)%lk_abs)) then
        print*, 'ERROR - LBL table did not provide all required interpolation data - STOPPING'
        print*, 'Species, path: ', lbl_tab(s)%sp, trim(lbl_tab(s)%path)
        stop
      end if

      if (lbl_tab(s)%nT < 3 .or. lbl_tab(s)%nP < 3) then
        print*, 'ERROR - LBL Bezier interpolation requires at least 3 T and P points - STOPPING'
        print*, 'Species, nT, nP: ', lbl_tab(s)%sp, lbl_tab(s)%nT, lbl_tab(s)%nP
        stop
      end if

      if (lbl_tab(s)%nwl /= nwl) then
        print*, 'ERROR - LBL table wavelength count does not match wavelengths.wl - STOPPING'
        print*, 'Species, table nwl, calculation nwl: ', lbl_tab(s)%sp, lbl_tab(s)%nwl, nwl
        stop
      end if

      if (any(.not. ieee_is_finite(lbl_tab(s)%T)) .or. &
        & any(.not. ieee_is_finite(lbl_tab(s)%P)) .or. &
        & any(.not. ieee_is_finite(lbl_tab(s)%wl)) .or. &
        & any(.not. ieee_is_finite(lbl_tab(s)%lk_abs))) then
        print*, 'ERROR - LBL interpolation data must be finite - STOPPING'
        print*, 'Species, path: ', lbl_tab(s)%sp, trim(lbl_tab(s)%path)
        stop
      end if

      if (any(lbl_tab(s)%T <= 0.0_dp) .or. any(lbl_tab(s)%P <= 0.0_dp) .or. &
        & any(lbl_tab(s)%wl <= 0.0_dp)) then
        print*, 'ERROR - LBL temperature, pressure, and wavelength grids must be positive - STOPPING'
        print*, 'Species, path: ', lbl_tab(s)%sp, trim(lbl_tab(s)%path)
        stop
      end if

      do i = 2, lbl_tab(s)%nT
        if (lbl_tab(s)%T(i) <= lbl_tab(s)%T(i-1)) then
          print*, 'ERROR - LBL temperature grid must be strictly increasing - STOPPING'
          print*, 'Species, index, values: ', lbl_tab(s)%sp, i, &
            & lbl_tab(s)%T(i-1), lbl_tab(s)%T(i)
          stop
        end if
      end do

      do i = 2, lbl_tab(s)%nP
        if (lbl_tab(s)%P(i) <= lbl_tab(s)%P(i-1)) then
          print*, 'ERROR - LBL pressure grid must be strictly increasing - STOPPING'
          print*, 'Species, index, values: ', lbl_tab(s)%sp, i, &
            & lbl_tab(s)%P(i-1), lbl_tab(s)%P(i)
          stop
        end if
      end do

      ! Bin centres can legitimately differ slightly between the LBL table
      ! and wavelengths.wl (e.g. different centring conventions), so only
      ! check the overall range rather than an exact per-point match.
      scale = max(1.0_dp, abs(lbl_tab(s)%wl(1)), abs(wl(1)))
      if (abs(lbl_tab(s)%wl(1) - wl(1)) > range_tol*scale) then
        print*, 'ERROR - LBL wavelength grid range does not match wavelengths.wl - STOPPING'
        print*, 'Species, table wl(1), calculation wl(1): ', lbl_tab(s)%sp, &
          & lbl_tab(s)%wl(1), wl(1)
        stop
      end if

      scale = max(1.0_dp, abs(lbl_tab(s)%wl(nwl)), abs(wl(nwl)))
      if (abs(lbl_tab(s)%wl(nwl) - wl(nwl)) > range_tol*scale) then
        print*, 'ERROR - LBL wavelength grid range does not match wavelengths.wl - STOPPING'
        print*, 'Species, table wl(nwl), calculation wl(nwl): ', lbl_tab(s)%sp, &
          & lbl_tab(s)%wl(nwl), wl(nwl)
        stop
      end if
    end do

  end subroutine validate_lbl_tables


  subroutine output_lbl_table(l)
    implicit none

    integer, intent(in) :: l
    integer :: reclen

    if (first_call .eqv. .True.) then
      inquire(iolength=reclen) lbl_write
      ! Output lbl-table in 1D flattened 3D CMCRT format lbl.cmcrt (single precision)
      open(newunit=ulbl, file='lbl.cmcrt', action='readwrite', &
      & form='unformatted',status='replace',access='direct',recl=reclen)
      first_call = .False.
    end if

    ! Convert to single precision on output and protect against underflow.
    lbl_write(:) = real(max(lbl_out(:),1.0e-30_dp),kind=sp)
    write(ulbl,rec=l) lbl_write

  end subroutine output_lbl_table

end module lbl_tables_mod
