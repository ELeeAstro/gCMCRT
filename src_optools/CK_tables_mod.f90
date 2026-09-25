module CK_tables_mod
  use optools_data_mod
  use CK_tables_read
  use CK_tables_interp
  use CK_table_RO
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
#ifdef GPU_OFFLOAD
  use omp_lib, only : omp_get_num_devices, omp_is_initial_device
  use, intrinsic :: iso_fortran_env, only : int64
#endif
  !use CK_table_rebin
  implicit none

#ifdef GPU_OFFLOAD
#ifndef CK_GPU_BLOCK
#define CK_GPU_BLOCK 8
#endif
#endif

  logical :: first_call = .True.

  real(kind=dp), allocatable, dimension(:,:) :: CK_out
  real(kind=sp), allocatable, dimension(:,:) :: CK_write

  real(kind=dp), allocatable, dimension(:) :: Gw, Gx

  ! Namelist variables
  integer :: iopts, nG, gdist
  real(kind=dp) :: gmin1, gmax1, gmin2, gmax2
  integer, allocatable, dimension(:) :: form
  character(len=150), allocatable, dimension(:) :: paths
  logical :: pre_mixed, rebin
  integer :: nrebin
  logical :: interp_wl = .False.

  namelist /CK_nml/ iopts, form, paths, nG, gmin1, gmax1, gmin2, gmax2, &
    & pre_mixed, rebin, nrebin, interp_wl

  private :: output_CK_table, output_CK_gord, validate_CK_tables, calc_CK_table_cpu
#ifdef GPU_OFFLOAD
  private :: calc_CK_table_gpu, locate_triplet_gpu, Bezier_interp_gpu, sort2_gpu
  !$omp declare target(locate_triplet_gpu,Bezier_interp_gpu,sort2_gpu)
#endif
  public :: calc_CK_table


contains

  subroutine calc_CK_table()
    implicit none

    integer :: s, j
    logical :: exists

    ! The array-allocation order is important because it follows the namelist order.

    ! Allocate number of CK tables
    allocate(CK_tab(nCK))

    ! Allocate required number of arrays from namelist options
    allocate(form(nCK),paths(nCK))

    ! Read CK namelist parameters. Keep interp_wl as a compatibility-only
    ! input so older namelists containing .False. continue to parse.
    interp_wl = .False.
    read(u_nml, nml=CK_nml)

    if (interp_wl .eqv. .True.) then
      print*, 'ERROR - CK wavelength interpolation has been removed - STOPPING'
      print*, 'Use index-aligned wavelength bins from the correlated-k table.'
      stop 1
    end if

    if (rebin .eqv. .True.) then
      print*, 'ERROR - CK rebinning is not implemented - STOPPING'
      print*, 'Set rebin = .False. in CK_nml.'
      stop
    end if

    if (nG < 1) then
      print*, 'ERROR - CK_nml nG must be positive - STOPPING'
      stop
    end if

    ! Allocate work arrays
    allocate(CK_out(nG,nlay),CK_write(nG,nlay))

    allocate(Gx(nG),Gw(nG))

    ! Give the classes some global data from par and namelists
    CK_tab(:)%sp = CK_name(:)
    CK_tab(:)%form = form(:)
    CK_tab(:)%path = paths(:)

    ! Find the PRF VMR indices of the CK species.
    if (pre_mixed .eqv. .False.) then
      do s = 1, nCK
        exists = .False.
        do j = 1, ngas
          if (CK_tab(s)%sp == g_name(j)) then
            CK_tab(s)%iVMR = j
            exists = .True.
            exit
          end if
        end do
        if (exists .eqv. .False.) then
          print*, 'ERROR - Specified CK species not found in prf VMR list - STOPPING'
          print*, 'Species: ', CK_tab(s)%sp
          stop
        end if
      end do
    end if

    ! Read the CK tables
    call read_CK_tables(pre_mixed)

    ! Validate table dimensions and grids before interpolation or random overlap
    call validate_CK_tables()

    ! Rebin each k table if requested.
    !if (rebin .eqv. .True.) then
      !call rebin_CK_tables(nrebin,nG)
    !end if

    ! Use the first table's g grid and weights.
    Gx(:) = CK_tab(1)%Gx(:)
    Gw(:) = CK_tab(1)%Gw(:)

#ifdef GPU_OFFLOAD
    call calc_CK_table_gpu()
#else
    call calc_CK_table_cpu()
#endif

    print*, ' ~~ Quest completed ~~ '

    ! Write the g ordinates and weights.
    call output_CK_gord()

    ! Deallocate all arrays.
    deallocate(CK_out,CK_write)
    deallocate(Gx,Gw)
    deallocate(CK_tab)
    deallocate(form,paths)
    ! Close the CK output unit.
    close(uCK)

  end subroutine calc_CK_table

  subroutine calc_CK_table_cpu()
    implicit none

    integer :: l, z
    real(kind=dp), allocatable, dimension(:) :: CK_RO
    real(kind=dp), allocatable, dimension(:,:) :: CK_work

    allocate(CK_work(nCK,nG),CK_RO(nG))
    CK_work(:,:) = 0.0_dp
    CK_RO(:) = 0.0_dp

    if (pre_mixed .eqv. .False.) then
      print*, ' ~~ Performing CK interpolation, RO and output ~~ '
    else
      print*, ' ~~ Performing CK premixed interpolation and output ~~ '
    end if
    print*, ' ~~ Please wait... ~~ '

    !! Begin OpenMP loops.
    !$omp parallel default (none), &
    !$omp& private (l,z), &
    !$omp& shared (nwl,wl,nlay,CK_out,RH_lay,nG,Gw,pre_mixed,N_lay), &
    !$omp& firstprivate(CK_work,CK_RO)


    ! Interpolate the CK tables to the layer temperature and pressure.
    ! Species loops are contained within the interpolation subroutines.
    do l = 1, nwl
      !$omp single
      if (mod(l,max(1,nwl/10)) == 0) then
        print*, l, wl(l), nwl
      end if
      !$omp end single
      !$omp do schedule (dynamic)
      do z = 1, nlay

        ! Correlated-k wavelength bins are index-aligned. Only temperature and
        ! pressure are interpolated.
        call interp_CK_tables_Bezier(l,z,nG,CK_work(:,:))

        if (pre_mixed .eqv. .True.) then
          ! Return the interpolated pre-mixed CK table to the output array.
          CK_out(:,z) = (CK_work(1,:)*N_lay(z))/RH_lay(z)
          !CK_out(:,z) = (CK_work(1,:))/RH_lay(z)
        else
          ! Perform random overlap with resorting and rebinning for all species.
          call RO_CK_RORR(z,nG,Gw(:),CK_work(:,:),CK_RO(:))
          !call RO_CK_2(z,nG,Gw(:),Gx(:),CK_work(:,:),CK_RO(:))
          !call RO_CK(z,nG,Gw(:),CK_work(:,:),CK_RO(:))
          ! Convert the overlapped result to cm^2 g^-1 of atmosphere.
          CK_out(:,z) = CK_RO(:)/RH_lay(z)
        end if


      end do
      !$omp end do

      !$omp single
      ! Write the CMCRT-formatted CK table for this wavelength bin.
      call output_CK_table(l)
      !$omp end single

    end do
    !$omp end parallel

    deallocate(CK_work,CK_RO)

  end subroutine calc_CK_table_cpu

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

  pure subroutine sort2_gpu(n,offset,ra,rb)
    implicit none

    integer, intent(in) :: n
    integer(kind=int64), intent(in) :: offset
    real(kind=dp), dimension(*), intent(inout) :: ra, rb
    integer :: l, ir, i, j
    real(kind=dp) :: rra, rrb

    if (n <= 1) return

    l = n/2 + 1
    ir = n
    do
      if (l > 1) then
        l = l - 1
        rra = ra(offset+int(l,kind=int64))
        rrb = rb(offset+int(l,kind=int64))
      else
        rra = ra(offset+int(ir,kind=int64))
        rrb = rb(offset+int(ir,kind=int64))
        ra(offset+int(ir,kind=int64)) = ra(offset+1_int64)
        rb(offset+int(ir,kind=int64)) = rb(offset+1_int64)
        ir = ir - 1
        if (ir == 1) then
          ra(offset+1_int64) = rra
          rb(offset+1_int64) = rrb
          return
        end if
      end if

      i = l
      j = l + l
      do while (j <= ir)
        if (j < ir) then
          if (ra(offset+int(j,kind=int64)) < &
              & ra(offset+int(j+1,kind=int64))) j = j + 1
        end if
        if (rra < ra(offset+int(j,kind=int64))) then
          ra(offset+int(i,kind=int64)) = ra(offset+int(j,kind=int64))
          rb(offset+int(i,kind=int64)) = rb(offset+int(j,kind=int64))
          i = j
          j = j + j
        else
          exit
        end if
      end do
      ra(offset+int(i,kind=int64)) = rra
      rb(offset+int(i,kind=int64)) = rrb
    end do

  end subroutine sort2_gpu

  subroutine calc_CK_table_gpu()
    implicit none

    integer, parameter :: block_size = CK_GPU_BLOCK
    real(kind=dp), parameter :: VMR_skip = 1.0e-30_dp
    real(kind=dp), parameter :: rebin_eps = 1.0e-8_dp
    integer :: s, g, i, j, m, itab, ptab
    integer :: l, l_start, l_end, n_block, lb, z
    integer :: T_region, P_region, T_exact, P_exact, T_fixed, P_fixed
    integer :: iT_idx(3), iP_idx(3), iT1, iT2, iT3, iP1, iP2, iP3
    integer :: nG2, n_active, offload_active
    integer, allocatable :: CK_nT(:), CK_nP(:), CK_iVMR(:)
    integer(kind=int64) :: total_T, total_P, total_cells
    integer(kind=int64) :: T_pos, P_pos, block_values, block_pos, k_idx
    integer(kind=int64) :: cell_count, cell_idx, work_stride, work_values
    integer(kind=int64) :: nG2_64, mix_stride, mix_values
    integer(kind=int64) :: cumulative_stride, cumulative_values
    integer(kind=int64) :: work_base, work_idx, mix_base, cumulative_base
    integer(kind=int64), allocatable :: CK_T_offset(:), CK_P_offset(:), CK_K_offset(:)
    real(kind=dp), allocatable :: CK_T_flat(:), CK_lT_flat(:)
    real(kind=dp), allocatable :: CK_P_flat(:), CK_lP_flat(:)
    real(kind=dp), allocatable :: CK_K_block(:), CK_work_block(:)
    real(kind=dp), allocatable :: k_mix_block(:), wt_mix_block(:)
    real(kind=dp), allocatable :: source_cumulative_block(:)
    real(kind=dp), allocatable :: CK_out_block(:,:,:)
    real(kind=dp), allocatable :: Gw_norm(:), target_cumulative(:)
    real(kind=dp) :: weight_sum, T_layer, P_layer, lT_layer, lP_layer
    real(kind=dp) :: log_k, q, target_lo, target_hi, source_lo, source_hi
    real(kind=dp) :: overlap, bin_width, bin_sum, boundary_tol
    real(kind=dp) :: lTa(3), lPa(3), lka(3), lka_ck(3)
    real(kind=dp) :: table_mib, scratch_mib

    if (block_size < 1) then
      print*, 'ERROR - GPU CK wavelength block must be positive - STOPPING'
      stop 1
    end if

    if (any(TG_lay <= 0.0_dp) .or. any(.not. ieee_is_finite(TG_lay))) then
      print*, 'ERROR - Invalid atmospheric temperature for GPU CK interpolation - STOPPING'
      stop 1
    end if
    if (any(PG_lay <= 0.0_dp) .or. any(.not. ieee_is_finite(PG_lay))) then
      print*, 'ERROR - Invalid atmospheric pressure for GPU CK interpolation - STOPPING'
      stop 1
    end if
    if (any(RH_lay <= 0.0_dp) .or. any(.not. ieee_is_finite(RH_lay)) .or. &
      & any(N_lay < 0.0_dp) .or. any(.not. ieee_is_finite(N_lay))) then
      print*, 'ERROR - Invalid atmospheric density for GPU CK interpolation - STOPPING'
      stop 1
    end if
    if (.not. pre_mixed) then
      do s = 1, nCK
        if (any(VMR_lay(CK_tab(s)%iVMR,:) < 0.0_dp) .or. &
          & any(.not. ieee_is_finite(VMR_lay(CK_tab(s)%iVMR,:)))) then
          print*, 'ERROR - Invalid CK species VMR for GPU interpolation: ', CK_tab(s)%sp
          stop 1
        end if
      end do
    end if

    if (any(Gw <= 0.0_dp) .or. any(.not. ieee_is_finite(Gw))) then
      print*, 'ERROR - GPU CK g weights must be finite and positive - STOPPING'
      stop 1
    end if
    weight_sum = sum(Gw)
    if (.not. ieee_is_finite(weight_sum) .or. weight_sum <= 0.0_dp) then
      print*, 'ERROR - Invalid total GPU CK g weight - STOPPING'
      stop 1
    end if

    allocate(CK_nT(nCK),CK_nP(nCK),CK_iVMR(nCK))
    allocate(CK_T_offset(nCK),CK_P_offset(nCK),CK_K_offset(nCK))

    total_T = 0_int64
    total_P = 0_int64
    total_cells = 0_int64
    do s = 1, nCK
      CK_nT(s) = CK_tab(s)%nT
      CK_nP(s) = CK_tab(s)%nP
      CK_iVMR(s) = 1
      if (.not. pre_mixed) CK_iVMR(s) = CK_tab(s)%iVMR
      CK_T_offset(s) = total_T + 1_int64
      CK_P_offset(s) = total_P + 1_int64
      CK_K_offset(s) = total_cells*int(block_size,kind=int64) + 1_int64
      total_T = total_T + int(CK_nT(s),kind=int64)
      total_P = total_P + int(CK_nP(s),kind=int64)
      total_cells = total_cells + int(CK_nT(s),kind=int64) &
        & * int(CK_nP(s),kind=int64)*int(nG,kind=int64)
    end do

    if (total_cells > huge(0_int64)/int(block_size,kind=int64)) then
      print*, 'ERROR - GPU CK packed table size overflow - STOPPING'
      stop 1
    end if
    block_values = total_cells*int(block_size,kind=int64)

    cell_count = int(block_size,kind=int64)*int(nlay,kind=int64)
    work_stride = int(nCK,kind=int64)*int(nG,kind=int64)
    if (cell_count > huge(0_int64)/max(1_int64,work_stride)) then
      print*, 'ERROR - GPU CK interpolation workspace size overflow - STOPPING'
      stop 1
    end if
    work_values = max(1_int64,cell_count*work_stride)

    nG2_64 = int(nG,kind=int64)*int(nG,kind=int64)
    if (nG2_64 > int(huge(nG2),kind=int64)) then
      print*, 'ERROR - GPU CK nG squared exceeds the supported integer range - STOPPING'
      stop 1
    end if
    nG2 = int(nG2_64)
    mix_stride = nG2_64
    cumulative_stride = mix_stride + 1_int64
    if (pre_mixed) then
      mix_values = 1_int64
      cumulative_values = 1_int64
    else
      if (cell_count > huge(0_int64)/max(1_int64,mix_stride) .or. &
        & cell_count > huge(0_int64)/cumulative_stride) then
        print*, 'ERROR - GPU CK random-overlap workspace size overflow - STOPPING'
        stop 1
      end if
      mix_values = max(1_int64,cell_count*mix_stride)
      cumulative_values = max(1_int64,cell_count*cumulative_stride)
    end if

    allocate(CK_T_flat(total_T),CK_lT_flat(total_T))
    allocate(CK_P_flat(total_P),CK_lP_flat(total_P))
    allocate(CK_K_block(block_values),CK_work_block(work_values))
    allocate(k_mix_block(mix_values),wt_mix_block(mix_values))
    allocate(source_cumulative_block(cumulative_values))
    allocate(CK_out_block(nG,nlay,block_size))
    allocate(Gw_norm(nG),target_cumulative(nG+1))

    CK_K_block(:) = 0.0_dp
    CK_out_block(:,:,:) = 0.0_dp
    Gw_norm(:) = Gw(:)/weight_sum
    target_cumulative(1) = 0.0_dp
    do g = 1, nG
      target_cumulative(g+1) = target_cumulative(g) + Gw_norm(g)
    end do
    target_cumulative(nG+1) = 1.0_dp

    do s = 1, nCK
      T_pos = CK_T_offset(s)
      P_pos = CK_P_offset(s)
      CK_T_flat(T_pos:T_pos+int(CK_nT(s)-1,kind=int64)) = CK_tab(s)%T(:)
      CK_lT_flat(T_pos:T_pos+int(CK_nT(s)-1,kind=int64)) = CK_tab(s)%lT(:)
      CK_P_flat(P_pos:P_pos+int(CK_nP(s)-1,kind=int64)) = CK_tab(s)%P(:)
      CK_lP_flat(P_pos:P_pos+int(CK_nP(s)-1,kind=int64)) = CK_tab(s)%lP(:)
    end do

    offload_active = 0
    !$omp target map(from: offload_active)
    if (.not. omp_is_initial_device()) offload_active = 1
    !$omp end target
    if (offload_active /= 1) then
      print*, 'ERROR - CK OpenMP target region did not execute on a GPU - STOPPING'
      print*, 'OpenMP target devices visible: ', omp_get_num_devices()
      stop 1
    end if

    table_mib = real(block_values,kind=dp)*real(storage_size(1.0_dp)/8,kind=dp) &
      & /(1024.0_dp**2)
    scratch_mib = real(work_values+2_int64*mix_values+cumulative_values,kind=dp) &
      & * real(storage_size(1.0_dp)/8,kind=dp)/(1024.0_dp**2)
    print*, ' ~~ CK OpenMP GPU offload active; visible devices: ', omp_get_num_devices()
    print*, ' ~~ CK GPU wavelength block size: ', block_size
    print*, ' ~~ CK GPU packed table storage [MiB]: ', table_mib
    print*, ' ~~ CK GPU scratch storage [MiB]: ', scratch_mib
    if (pre_mixed) then
      print*, ' ~~ Performing CK premixed interpolation and output ~~ '
    else
      print*, ' ~~ Performing CK interpolation, RO and output ~~ '
    end if
    print*, ' ~~ Please wait... ~~ '

    !$omp target data map(to: TG_lay,PG_lay,N_lay,RH_lay,VMR_lay,nlay,nCK,nG,nG2,pre_mixed) &
    !$omp& map(to: work_stride,mix_stride,cumulative_stride) &
    !$omp& map(to: CK_nT,CK_nP,CK_iVMR,CK_T_offset,CK_P_offset,CK_K_offset) &
    !$omp& map(to: CK_T_flat,CK_lT_flat,CK_P_flat,CK_lP_flat,Gw_norm,target_cumulative) &
    !$omp& map(alloc: CK_K_block,CK_work_block,k_mix_block,wt_mix_block) &
    !$omp& map(alloc: source_cumulative_block,CK_out_block)

    do l_start = 1, nwl, block_size
      l_end = min(nwl,l_start+block_size-1)
      n_block = l_end-l_start+1

      do s = 1, nCK
        do g = 1, nG
          do itab = 1, CK_nT(s)
            do ptab = 1, CK_nP(s)
              block_pos = CK_K_offset(s) &
                & + ((int(g-1,kind=int64)*int(CK_nT(s),kind=int64) &
                & + int(itab-1,kind=int64))*int(CK_nP(s),kind=int64) &
                & + int(ptab-1,kind=int64))*int(block_size,kind=int64)
              CK_K_block(block_pos:block_pos+int(n_block-1,kind=int64)) = &
                & CK_tab(s)%lk_abs(l_start:l_end,ptab,itab,g)
            end do
          end do
        end do
      end do
      !$omp target update to(CK_K_block)

      !$omp target teams loop collapse(2) &
      !$omp& private(s,g,i,j,m,T_region,P_region,T_exact,P_exact,T_fixed,P_fixed) &
      !$omp& private(iT_idx,iP_idx,iT1,iT2,iT3,iP1,iP2,iP3,n_active) &
      !$omp& private(cell_idx,work_base,work_idx,mix_base,cumulative_base,k_idx) &
      !$omp& private(T_layer,P_layer,lT_layer,lP_layer,log_k,q) &
      !$omp& private(target_lo,target_hi,source_lo,source_hi,overlap,bin_width) &
      !$omp& private(bin_sum,boundary_tol,lTa,lPa,lka,lka_ck)
      do lb = 1, n_block
        do z = 1, nlay
          cell_idx = int(lb-1,kind=int64)*int(nlay,kind=int64) &
            & + int(z-1,kind=int64)
          work_base = cell_idx*work_stride
          mix_base = cell_idx*mix_stride
          cumulative_base = cell_idx*cumulative_stride
          T_layer = TG_lay(z)
          P_layer = PG_lay(z)
          lT_layer = log10(T_layer)
          lP_layer = log10(P_layer)

          do s = 1, nCK
            call locate_triplet_gpu(CK_T_flat,CK_T_offset(s),CK_nT(s), &
              & T_layer,iT_idx,T_region,T_exact)
            call locate_triplet_gpu(CK_P_flat,CK_P_offset(s),CK_nP(s), &
              & P_layer,iP_idx,P_region,P_exact)

            iT1 = iT_idx(1)
            iT2 = iT_idx(2)
            iT3 = iT_idx(3)
            iP1 = iP_idx(1)
            iP2 = iP_idx(2)
            iP3 = iP_idx(3)
            lTa(1) = CK_lT_flat(CK_T_offset(s)+int(iT1-1,kind=int64))
            lTa(2) = CK_lT_flat(CK_T_offset(s)+int(iT2-1,kind=int64))
            lTa(3) = CK_lT_flat(CK_T_offset(s)+int(iT3-1,kind=int64))
            lPa(1) = CK_lP_flat(CK_P_offset(s)+int(iP1-1,kind=int64))
            lPa(2) = CK_lP_flat(CK_P_offset(s)+int(iP2-1,kind=int64))
            lPa(3) = CK_lP_flat(CK_P_offset(s)+int(iP3-1,kind=int64))

            T_fixed = T_exact
            if (T_region == -1) T_fixed = 1
            if (T_region == 1) T_fixed = CK_nT(s)
            P_fixed = P_exact
            if (P_region == -1) P_fixed = 1
            if (P_region == 1) P_fixed = CK_nP(s)

            do g = 1, nG
              if (T_fixed > 0 .and. P_fixed > 0) then
                k_idx = CK_K_offset(s) &
                  & + ((int(g-1,kind=int64)*int(CK_nT(s),kind=int64) &
                  & + int(T_fixed-1,kind=int64))*int(CK_nP(s),kind=int64) &
                  & + int(P_fixed-1,kind=int64))*int(block_size,kind=int64) &
                  & + int(lb-1,kind=int64)
                log_k = CK_K_block(k_idx)
              else if (T_fixed > 0) then
                do j = 1, 3
                  k_idx = CK_K_offset(s) &
                    & + ((int(g-1,kind=int64)*int(CK_nT(s),kind=int64) &
                    & + int(T_fixed-1,kind=int64))*int(CK_nP(s),kind=int64) &
                    & + int(iP_idx(j)-1,kind=int64))*int(block_size,kind=int64) &
                    & + int(lb-1,kind=int64)
                  lka(j) = CK_K_block(k_idx)
                end do
                call Bezier_interp_gpu(lPa,lka,lP_layer,log_k)
              else if (P_fixed > 0) then
                do j = 1, 3
                  k_idx = CK_K_offset(s) &
                    & + ((int(g-1,kind=int64)*int(CK_nT(s),kind=int64) &
                    & + int(iT_idx(j)-1,kind=int64))*int(CK_nP(s),kind=int64) &
                    & + int(P_fixed-1,kind=int64))*int(block_size,kind=int64) &
                    & + int(lb-1,kind=int64)
                  lka(j) = CK_K_block(k_idx)
                end do
                call Bezier_interp_gpu(lTa,lka,lT_layer,log_k)
              else
                do j = 1, 3
                  k_idx = CK_K_offset(s) &
                    & + ((int(g-1,kind=int64)*int(CK_nT(s),kind=int64) &
                    & + int(iT_idx(j)-1,kind=int64))*int(CK_nP(s),kind=int64) &
                    & + int(iP1-1,kind=int64))*int(block_size,kind=int64) &
                    & + int(lb-1,kind=int64)
                  lka(1) = CK_K_block(k_idx)
                  k_idx = CK_K_offset(s) &
                    & + ((int(g-1,kind=int64)*int(CK_nT(s),kind=int64) &
                    & + int(iT_idx(j)-1,kind=int64))*int(CK_nP(s),kind=int64) &
                    & + int(iP2-1,kind=int64))*int(block_size,kind=int64) &
                    & + int(lb-1,kind=int64)
                  lka(2) = CK_K_block(k_idx)
                  k_idx = CK_K_offset(s) &
                    & + ((int(g-1,kind=int64)*int(CK_nT(s),kind=int64) &
                    & + int(iT_idx(j)-1,kind=int64))*int(CK_nP(s),kind=int64) &
                    & + int(iP3-1,kind=int64))*int(block_size,kind=int64) &
                    & + int(lb-1,kind=int64)
                  lka(3) = CK_K_block(k_idx)
                  call Bezier_interp_gpu(lPa,lka,lP_layer,lka_ck(j))
                end do
                call Bezier_interp_gpu(lTa,lka_ck,lT_layer,log_k)
              end if
              work_idx = work_base + int(s-1,kind=int64)*int(nG,kind=int64) &
                & + int(g,kind=int64)
              CK_work_block(work_idx) = 10.0_dp**log_k
            end do
          end do

          if (pre_mixed) then
            do g = 1, nG
              CK_out_block(g,z,lb) = CK_work_block(work_base+int(g,kind=int64)) &
                & * N_lay(z)/RH_lay(z)
            end do
          else
            do g = 1, nG
              CK_out_block(g,z,lb) = 0.0_dp
            end do
            n_active = 0

            do s = 1, nCK
              q = VMR_lay(CK_iVMR(s),z)
              if (q < VMR_skip) cycle

              if (n_active == 0) then
                do g = 1, nG
                  work_idx = work_base + int(s-1,kind=int64)*int(nG,kind=int64) &
                    & + int(g,kind=int64)
                  CK_out_block(g,z,lb) = q*CK_work_block(work_idx)
                end do
                n_active = 1
                cycle
              end if

              do i = 1, nG
                do j = 1, nG
                  m = (i-1)*nG+j
                  work_idx = work_base + int(s-1,kind=int64)*int(nG,kind=int64) &
                    & + int(j,kind=int64)
                  k_mix_block(mix_base+int(m,kind=int64)) = &
                    & CK_out_block(i,z,lb) + q*CK_work_block(work_idx)
                  wt_mix_block(mix_base+int(m,kind=int64)) = Gw_norm(i)*Gw_norm(j)
                end do
              end do

              call sort2_gpu(nG2,mix_base,k_mix_block,wt_mix_block)

              source_cumulative_block(cumulative_base+1_int64) = 0.0_dp
              do m = 1, nG2
                source_cumulative_block(cumulative_base+int(m+1,kind=int64)) = &
                  & source_cumulative_block(cumulative_base+int(m,kind=int64)) &
                  & + wt_mix_block(mix_base+int(m,kind=int64))
              end do
              source_cumulative_block(cumulative_base+int(nG2+1,kind=int64)) = 1.0_dp

              m = 1
              do g = 1, nG
                target_lo = target_cumulative(g)
                target_hi = target_cumulative(g+1)
                bin_width = target_hi-target_lo
                boundary_tol = rebin_eps*target_hi
                bin_sum = 0.0_dp

                do while (m <= nG2)
                  source_lo = source_cumulative_block(cumulative_base+int(m,kind=int64))
                  source_hi = source_cumulative_block(cumulative_base+int(m+1,kind=int64))

                  if (source_hi >= target_hi .and. &
                    & source_hi-target_hi <= boundary_tol .and. target_hi > source_lo) then
                    if (m == nG2) then
                      source_cumulative_block(cumulative_base+int(m+1,kind=int64)) = target_hi
                      source_hi = target_hi
                    else if (target_hi < source_cumulative_block( &
                      & cumulative_base+int(m+2,kind=int64))) then
                      source_cumulative_block(cumulative_base+int(m+1,kind=int64)) = target_hi
                      source_hi = target_hi
                    end if
                  end if

                  if (source_hi <= target_lo) then
                    m = m + 1
                    cycle
                  end if
                  if (source_lo >= target_hi) exit

                  overlap = min(target_hi,source_hi)-max(target_lo,source_lo)
                  if (overlap > 0.0_dp) then
                    bin_sum = bin_sum + overlap*k_mix_block(mix_base+int(m,kind=int64))
                  end if
                  if (source_hi >= target_hi) exit
                  m = m + 1
                end do
                CK_out_block(g,z,lb) = bin_sum/bin_width
              end do
              n_active = n_active + 1
            end do

            if (n_active > 0) then
              do g = 1, nG
                CK_out_block(g,z,lb) = CK_out_block(g,z,lb)*N_lay(z)/RH_lay(z)
              end do
            end if
          end if
        end do
      end do
      !$omp end target teams loop

      !$omp target update from(CK_out_block)
      if (any(.not. ieee_is_finite(CK_out_block(:,:,1:n_block)))) then
        print*, 'ERROR - Non-finite CK opacity returned by GPU in wavelength block: ', &
          & l_start, l_end
        stop 1
      end if

      do lb = 1, n_block
        l = l_start+lb-1
        CK_out(:,:) = CK_out_block(:,:,lb)
        call output_CK_table(l)
        if (mod(l,max(1,nwl/10)) == 0) print*, l, wl(l), nwl
      end do
    end do

    !$omp end target data

    deallocate(CK_nT,CK_nP,CK_iVMR)
    deallocate(CK_T_offset,CK_P_offset,CK_K_offset)
    deallocate(CK_T_flat,CK_lT_flat,CK_P_flat,CK_lP_flat)
    deallocate(CK_K_block,CK_work_block,k_mix_block,wt_mix_block)
    deallocate(source_cumulative_block,CK_out_block,Gw_norm,target_cumulative)

  end subroutine calc_CK_table_gpu
#endif

  subroutine validate_CK_tables()
    implicit none

    integer :: s, i, wavelength_direction
    real(kind=dp), parameter :: grid_tol = 1.0e-10_dp
    real(kind=dp), parameter :: range_tol = 1.0e-2_dp
    real(kind=dp) :: scale

    if (nCK < 1) then
      print*, 'ERROR - Correlated-k opacity is enabled but no CK tables are configured - STOPPING'
      stop
    end if

    if (nwl < 1) then
      print*, 'ERROR - Correlated-k calculation requires at least one wavelength - STOPPING'
      stop
    end if

    wavelength_direction = 1
    if (nwl > 1) then
      if (wl(nwl) < wl(1)) wavelength_direction = -1
      do i = 2, nwl
        if (real(wavelength_direction,kind=dp)*(wl(i)-wl(i-1)) <= 0.0_dp) then
          print*, 'ERROR - Calculation wavelength grid must be strictly ordered - STOPPING'
          print*, 'Index, values: ', i, wl(i-1), wl(i)
          stop
        end if
      end do
    end if

    do s = 1, nCK
      if (.not. allocated(CK_tab(s)%Gx) .or. .not. allocated(CK_tab(s)%Gw) .or. &
        & .not. allocated(CK_tab(s)%wl) .or. .not. allocated(CK_tab(s)%T) .or. &
        & .not. allocated(CK_tab(s)%lT) .or. .not. allocated(CK_tab(s)%P) .or. &
        & .not. allocated(CK_tab(s)%lP) .or. .not. allocated(CK_tab(s)%lk_abs)) then
        print*, 'ERROR - CK table did not provide all required grids - STOPPING'
        print*, 'Species, path: ', CK_tab(s)%sp, trim(CK_tab(s)%path)
        stop
      end if

      if (CK_tab(s)%nG /= nG) then
        print*, 'ERROR - CK table nG does not match the CK namelist - STOPPING'
        print*, 'Species, file nG, namelist nG: ', CK_tab(s)%sp, CK_tab(s)%nG, nG
        stop
      end if

      if (CK_tab(s)%nwl /= nwl) then
        print*, 'ERROR - CK table wavelength count does not match wavelengths.wl - STOPPING'
        print*, 'Species, table nwl, calculation nwl: ', CK_tab(s)%sp, CK_tab(s)%nwl, nwl
        stop
      end if

      if (CK_tab(s)%nT < 3 .or. CK_tab(s)%nP < 3) then
        print*, 'ERROR - CK Bezier interpolation requires at least 3 T and P points - STOPPING'
        print*, 'Species, nT, nP: ', CK_tab(s)%sp, CK_tab(s)%nT, CK_tab(s)%nP
        stop
      end if

      if (any(CK_tab(s)%T <= 0.0_dp) .or. any(CK_tab(s)%P <= 0.0_dp)) then
        print*, 'ERROR - CK temperature and pressure grids must be positive - STOPPING'
        print*, 'Species, path: ', CK_tab(s)%sp, trim(CK_tab(s)%path)
        stop
      end if

      if (any(.not. ieee_is_finite(CK_tab(s)%T)) .or. &
        & any(.not. ieee_is_finite(CK_tab(s)%P)) .or. &
        & any(.not. ieee_is_finite(CK_tab(s)%wl)) .or. &
        & any(.not. ieee_is_finite(CK_tab(s)%Gx)) .or. &
        & any(.not. ieee_is_finite(CK_tab(s)%Gw)) .or. &
        & any(.not. ieee_is_finite(CK_tab(s)%lk_abs))) then
        print*, 'ERROR - CK interpolation data must be finite - STOPPING'
        print*, 'Species, path: ', CK_tab(s)%sp, trim(CK_tab(s)%path)
        stop
      end if

      do i = 2, CK_tab(s)%nT
        if (CK_tab(s)%T(i) <= CK_tab(s)%T(i-1)) then
          print*, 'ERROR - CK temperature grid must be strictly increasing - STOPPING'
          print*, 'Species, index, values: ', CK_tab(s)%sp, i, &
            & CK_tab(s)%T(i-1), CK_tab(s)%T(i)
          stop
        end if
      end do

      do i = 2, CK_tab(s)%nP
        if (CK_tab(s)%P(i) <= CK_tab(s)%P(i-1)) then
          print*, 'ERROR - CK pressure grid must be strictly increasing - STOPPING'
          print*, 'Species, index, values: ', CK_tab(s)%sp, i, &
            & CK_tab(s)%P(i-1), CK_tab(s)%P(i)
          stop
        end if
      end do

      if (any(CK_tab(s)%Gw <= 0.0_dp)) then
        print*, 'ERROR - CK g weights must be positive - STOPPING'
        print*, 'Species, path: ', CK_tab(s)%sp, trim(CK_tab(s)%path)
        stop
      end if

      do i = 2, nG
        if (CK_tab(s)%Gx(i) <= CK_tab(s)%Gx(i-1)) then
          print*, 'ERROR - CK g ordinates must be strictly increasing - STOPPING'
          print*, 'Species, index, values: ', CK_tab(s)%sp, i, &
            & CK_tab(s)%Gx(i-1), CK_tab(s)%Gx(i)
          stop
        end if
      end do

      do i = 2, nwl
        if (real(wavelength_direction,kind=dp) &
          & * (CK_tab(s)%wl(i)-CK_tab(s)%wl(i-1)) <= 0.0_dp) then
          print*, 'ERROR - CK wavelength grid ordering does not match wavelengths.wl - STOPPING'
          print*, 'Species, index, values: ', CK_tab(s)%sp, i, &
            & CK_tab(s)%wl(i-1), CK_tab(s)%wl(i)
          stop
        end if
      end do

      ! Bin centres can legitimately differ slightly between the CK table
      ! (arithmetic mean of edges) and wavelengths.wl (which may use a
      ! different centring convention). Check each positional bin with a
      ! modest tolerance, but never interpolate between wavelength bins.
      do i = 1, nwl
        scale = max(1.0_dp, abs(CK_tab(s)%wl(i)), abs(wl(i)))
        if (abs(CK_tab(s)%wl(i)-wl(i)) > range_tol*scale) then
          print*, 'ERROR - CK wavelength bins do not match wavelengths.wl - STOPPING'
          print*, 'Species, index, table wavelength, calculation wavelength: ', &
            & CK_tab(s)%sp, i, CK_tab(s)%wl(i), wl(i)
          stop
        end if
      end do
    end do

    do s = 2, nCK
      do i = 1, nG
        scale = max(1.0_dp, abs(CK_tab(s)%Gx(i)), abs(CK_tab(1)%Gx(i)))
        if (abs(CK_tab(s)%Gx(i) - CK_tab(1)%Gx(i)) > grid_tol*scale) then
          print*, 'ERROR - CK Gx grids differ between species - STOPPING'
          print*, 'Species, index, value, reference: ', CK_tab(s)%sp, i, &
            & CK_tab(s)%Gx(i), CK_tab(1)%Gx(i)
          stop
        end if

        scale = max(1.0_dp, abs(CK_tab(s)%Gw(i)), abs(CK_tab(1)%Gw(i)))
        if (abs(CK_tab(s)%Gw(i) - CK_tab(1)%Gw(i)) > grid_tol*scale) then
          print*, 'ERROR - CK Gw grids differ between species - STOPPING'
          print*, 'Species, index, value, reference: ', CK_tab(s)%sp, i, &
            & CK_tab(s)%Gw(i), CK_tab(1)%Gw(i)
          stop
        end if
      end do
    end do

  end subroutine validate_CK_tables

  subroutine output_CK_table(l)
    implicit none

    integer, intent(in) :: l
    integer :: reclen

    if (first_call .eqv. .True.) then
      ! Output k-table in 1D flattened 3D CMCRT format CK.cmcrt (single precision)
      inquire(iolength=reclen) CK_write
      open(newunit=uCK, file='CK.cmcrt', action='readwrite', &
              & form='unformatted', status='replace', access='direct',recl=reclen)
      first_call = .False.
    end if

    ! Convert to single precision on output and protect against underflow.
    CK_write = real(max(CK_out,1.0e-30_dp),kind=sp)
    write(uCK,rec=l) CK_write

  end subroutine output_CK_table

  subroutine output_CK_gord()
    implicit none
    integer :: g, u_g
    real(kind=dp) :: sum1

    print*, ' ~~ Outputting gord.cmcrt ~~ '

    ! Write the g ordinates and their weights to gord.cmcrt.
    open(newunit=u_g, file='gord.cmcrt', action='readwrite',form='formatted')
    write(u_g,*) nG

    sum1 = 0.0_dp
    do g = 1, nG
      write(u_g,*) Gx(g), Gw(g)
      sum1 = sum1 + Gw(g)
      !print*, g,  Gx(g), Gw(g), sum1
    end do

    print*,'G-ordinance sums:', sum(Gx(:)), sum(Gw(:))

    close(u_g)

    print*, ' ~~ Quest completed  ~~'

end subroutine output_CK_gord


end module CK_tables_mod
