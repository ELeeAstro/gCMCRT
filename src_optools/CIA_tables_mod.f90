module CIA_tables_mod
  use optools_data_mod
  use optools_table_class
  use CIA_tables_read, only : read_CIA_tables
  use CIA_tables_interp, only : interp_CIA_tables, interp_CIA_tables_Bezier
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
#ifdef GPU_OFFLOAD
  use omp_lib, only : omp_get_num_devices, omp_is_initial_device
  use CIA_tables_Hminus, only : An_ff1, Bn_ff1, Cn_ff1, Dn_ff1, En_ff1, Fn_ff1, &
    & An_ff2, Bn_ff2, Cn_ff2, Dn_ff2, En_ff2, Fn_ff2, Cn_bf
#endif
  implicit none

  logical :: first_call = .True.

  real(kind=dp), allocatable, dimension(:) :: CIA_out
  real(kind=sp), allocatable, dimension(:) :: CIA_write

  ! Namelist variables
  integer :: iopts
  integer, allocatable, dimension(:) :: form
  character(len=150), allocatable, dimension(:) :: paths

  namelist /CIA_nml/ iopts, form, paths

  private :: find_CIA_consituents, output_CIA_table, validate_CIA_tables, &
    & calc_CIA_table_cpu
#ifdef GPU_OFFLOAD
  private :: calc_CIA_table_gpu
#endif
  public :: calc_CIA_table

contains

  !! Driver routine for CIA table calculation

  subroutine calc_CIA_table()
    implicit none

    integer :: s, i, j, ni
    logical :: exists

    ! The array-allocation order is important because it follows the namelist order.

    ! Allocate number of CIA tables
    allocate(CIA_tab(nCIA))

    ! Allocate required number of arrays from namelist options
    allocate(form(nCIA),paths(nCIA))

    ! Read CIA namelist parameters
    read(u_nml, nml=CIA_nml)

    ! Allocate work arrays
    allocate(CIA_out(nlay),CIA_write(nlay))

    ! Give the classes some global data from par and namelists
    CIA_tab(:)%sp = CIA_name(:)
    CIA_tab(:)%form = form(:)
    CIA_tab(:)%path = paths(:)

    ! Find the CIA constituents from lookup table
    call find_CIA_consituents()

    ! Find the PRF VMR indices of the CIA constituent species.
    do s = 1, nCIA

      ! Check for 3 species special
      if (CIA_tab(s)%i3 .eqv. .True.) then
        ni = 3
      else
        ni = 2
      end if

      do i = 1, ni
        exists = .False.

        do j = 1, ngas
          if (ni == 2) then
            if (CIA_tab(s)%sp_con(i) == g_name(j)) then
              CIA_tab(s)%iVMR(i) = j
              exists = .True.
              exit
            end if
          else if (ni == 3) then
            if (CIA_tab(s)%sp_con_3(i) == g_name(j)) then
              CIA_tab(s)%iVMR_3(i) = j
              exists = .True.
              exit
            end if
          end if
        end do

        if (exists .eqv. .False.) then
          print*, 'ERROR - Specified CIA species component not found in prf VMR list - STOPPING'
          if (ni == 2) then
            print*, 'Species 2 part: ', CIA_tab(s)%sp, CIA_tab(s)%sp_con(i)
          else if (ni == 3) then
            print*, 'Species 3 part: ', CIA_tab(s)%sp, CIA_tab(s)%sp_con_3(i)
          end if
          stop
        end if

      end do
    end do

    ! Read the CIA tables
    call read_CIA_tables()

    call validate_CIA_tables()

#ifdef GPU_OFFLOAD
    call calc_CIA_table_gpu()
#else
    call calc_CIA_table_cpu()
#endif

    !deallocate all allocated arrays
    deallocate(CIA_out,CIA_write)
    deallocate(CIA_tab)
    deallocate(form,paths)
    ! Close the CIA I/O unit.
    close(uCIA)

    print*, ' ~~ Quest completed  ~~ '

  end subroutine calc_CIA_table

  subroutine calc_CIA_table_cpu()
    implicit none

    integer :: l, z
    real(kind=dp) :: CIA_work

    CIA_work = 0.0_dp

    print*, ' ~~ Performing CIA interpolation and output ~~ '
    print*, ' ~~ Please wait... ~~ '

    !$omp parallel default (none), &
    !$omp& private (l,z), &
    !$omp& shared (nwl,nlay,CIA_out,RH_lay,wl), &
    !$omp& firstprivate(CIA_work)

    do l = 1, nwl
      !$omp single
      if (mod(l,max(1,nwl/10)) == 0) then
        print*, l, wl(l), nwl
      end if
      !$omp end single

      !$omp do schedule (dynamic)
      do z = 1, nlay
        call interp_CIA_tables(l,z,CIA_work)
        CIA_out(z) = CIA_work/RH_lay(z)
      end do
      !$omp end do

      !$omp single
      call output_CIA_table(l)
      !$omp end single
    end do
    !$omp end parallel

  end subroutine calc_CIA_table_cpu

#ifdef GPU_OFFLOAD
  subroutine calc_CIA_table_gpu()
    implicit none

    integer, parameter :: CIA_MODEL_HMINUS = 1
    integer, parameter :: CIA_MODEL_HEMINUS = 2
    integer, parameter :: CIA_MODEL_HEMINUS_BELL = 3
    integer, parameter :: CIA_MODEL_H2MINUS_BELL = 4
    integer, parameter :: CIA_MODEL_H2O = 5
    integer :: s, sn, l, z, j, n, nTs, nrec
    integer :: iwn, iwn1, iT, iT1, jl, jm, ju
    integer :: max_nset, total_T, total_wn, total_tab
    integer :: T_pos, wn_pos, tab_pos, offload_active
    integer, allocatable :: CIA_form(:), CIA_model(:), CIA_nset(:)
    integer, allocatable :: CIA_iVMR(:,:), CIA_iVMR_3(:,:)
    integer, allocatable :: CIA_nT(:,:), CIA_irec(:,:)
    integer, allocatable :: CIA_T_offset(:,:), CIA_wn_offset(:,:), CIA_tab_offset(:,:)
    real(kind=dp), allocatable :: CIA_T_flat(:), CIA_wn_flat(:), CIA_wl_flat(:)
    real(kind=dp), allocatable :: CIA_tab_flat(:)
    real(kind=dp), allocatable :: CIA_wn_s(:,:), CIA_wn_e(:,:)
    real(kind=dp), allocatable :: CIA_Tmin(:,:), CIA_Tmax(:,:)
    real(kind=dp) :: CIA_work, CIA_value
    real(kind=dp) :: xval, yval, x0, x1, y0, y1
    real(kind=dp) :: a00, a10, a01, a11, norm
    real(kind=dp) :: lxval, lyval, lx0, lx1, ly0, ly1
    real(kind=dp) :: la00, la10, la01, la11
    real(kind=dp) :: T, T5040, a, b, c, kff, kbf, fbf, xbf, sff
    real(kind=dp), dimension(6,6) :: Hminus_ff1, Hminus_ff2
    real(kind=dp), dimension(6) :: Hminus_bf

    ! OpenMP does not portably deep-copy allocatable members of derived-type
    ! arrays. Pack the HITRAN tables into simple contiguous arrays once.
    max_nset = 1
    total_T = 0
    total_wn = 0
    total_tab = 0
    do s = 1, nCIA
      if (CIA_tab(s)%form == 4) then
        max_nset = max(max_nset,CIA_tab(s)%nset)
        do sn = 1, CIA_tab(s)%nset
          total_T = total_T + CIA_tab(s)%nT(sn)
          total_wn = total_wn + CIA_tab(s)%irec(sn)
          total_tab = total_tab + CIA_tab(s)%nT(sn)*CIA_tab(s)%irec(sn)
        end do
      else if (CIA_tab(s)%form == 2) then
        total_T = total_T + CIA_tab(s)%nT(1)
        total_wn = total_wn + CIA_tab(s)%nwl
        total_tab = total_tab + CIA_tab(s)%nT(1)*CIA_tab(s)%nwl
      end if
    end do

    allocate(CIA_form(nCIA),CIA_model(nCIA),CIA_nset(nCIA))
    allocate(CIA_iVMR(2,nCIA),CIA_iVMR_3(3,nCIA))
    allocate(CIA_nT(max_nset,nCIA),CIA_irec(max_nset,nCIA))
    allocate(CIA_T_offset(max_nset,nCIA),CIA_wn_offset(max_nset,nCIA))
    allocate(CIA_tab_offset(max_nset,nCIA))
    allocate(CIA_wn_s(max_nset,nCIA),CIA_wn_e(max_nset,nCIA))
    allocate(CIA_Tmin(max_nset,nCIA),CIA_Tmax(max_nset,nCIA))
    allocate(CIA_T_flat(max(1,total_T)),CIA_wn_flat(max(1,total_wn)))
    allocate(CIA_wl_flat(max(1,total_wn)))
    allocate(CIA_tab_flat(max(1,total_tab)))

    CIA_form(:) = 0
    CIA_model(:) = 0
    CIA_nset(:) = 0
    CIA_iVMR(:,:) = 1
    CIA_iVMR_3(:,:) = 1
    CIA_nT(:,:) = 0
    CIA_irec(:,:) = 0
    CIA_T_offset(:,:) = 1
    CIA_wn_offset(:,:) = 1
    CIA_tab_offset(:,:) = 1
    CIA_wn_s(:,:) = 0.0_dp
    CIA_wn_e(:,:) = 0.0_dp
    CIA_Tmin(:,:) = 0.0_dp
    CIA_Tmax(:,:) = 0.0_dp
    CIA_T_flat(:) = 0.0_dp
    CIA_wn_flat(:) = 0.0_dp
    CIA_wl_flat(:) = 0.0_dp
    CIA_tab_flat(:) = 0.0_dp

    Hminus_ff1(:,1) = An_ff1
    Hminus_ff1(:,2) = Bn_ff1
    Hminus_ff1(:,3) = Cn_ff1
    Hminus_ff1(:,4) = Dn_ff1
    Hminus_ff1(:,5) = En_ff1
    Hminus_ff1(:,6) = Fn_ff1
    Hminus_ff2(:,1) = An_ff2
    Hminus_ff2(:,2) = Bn_ff2
    Hminus_ff2(:,3) = Cn_ff2
    Hminus_ff2(:,4) = Dn_ff2
    Hminus_ff2(:,5) = En_ff2
    Hminus_ff2(:,6) = Fn_ff2
    Hminus_bf(:) = Cn_bf

    T_pos = 1
    wn_pos = 1
    tab_pos = 1
    do s = 1, nCIA
      CIA_form(s) = CIA_tab(s)%form
      if (CIA_form(s) == 4) then
        CIA_iVMR(:,s) = CIA_tab(s)%iVMR(:)
        CIA_nset(s) = CIA_tab(s)%nset
        do sn = 1, CIA_nset(s)
          nTs = CIA_tab(s)%nT(sn)
          nrec = CIA_tab(s)%irec(sn)
          CIA_nT(sn,s) = nTs
          CIA_irec(sn,s) = nrec
          CIA_T_offset(sn,s) = T_pos
          CIA_wn_offset(sn,s) = wn_pos
          CIA_tab_offset(sn,s) = tab_pos
          CIA_wn_s(sn,s) = CIA_tab(s)%wn_s(sn)
          CIA_wn_e(sn,s) = CIA_tab(s)%wn_e(sn)
          CIA_Tmin(sn,s) = CIA_tab(s)%Tmin(sn)
          CIA_Tmax(sn,s) = CIA_tab(s)%Tmax(sn)

          CIA_T_flat(T_pos:T_pos+nTs-1) = CIA_tab(s)%T(sn,1:nTs)
          CIA_wn_flat(wn_pos:wn_pos+nrec-1) = CIA_tab(s)%wn(sn,1:nrec)
          do j = 1, nTs
            CIA_tab_flat(tab_pos+(j-1)*nrec:tab_pos+j*nrec-1) = &
              & CIA_tab(s)%tab(sn,1:nrec,j)
          end do

          T_pos = T_pos + nTs
          wn_pos = wn_pos + nrec
          tab_pos = tab_pos + nTs*nrec
        end do
        cycle
      end if

      select case(trim(CIA_tab(s)%sp))
      case('H-')
        CIA_model(s) = CIA_MODEL_HMINUS
        CIA_iVMR_3(:,s) = CIA_tab(s)%iVMR_3(:)
      case('He-')
        CIA_iVMR(:,s) = CIA_tab(s)%iVMR(:)
        if (CIA_form(s) == 2) then
          CIA_model(s) = CIA_MODEL_HEMINUS_BELL
        else
          CIA_model(s) = CIA_MODEL_HEMINUS
        end if
      case('H2-')
        CIA_model(s) = CIA_MODEL_H2MINUS_BELL
        CIA_iVMR(:,s) = CIA_tab(s)%iVMR(:)
      case('H2O')
        CIA_model(s) = CIA_MODEL_H2O
        CIA_iVMR(:,s) = CIA_tab(s)%iVMR(:)
      case default
        print*, 'ERROR - CIA GPU special model not found - STOPPING'
        print*, 'Species: ', CIA_tab(s)%sp, CIA_form(s)
        stop 1
      end select

      if (CIA_form(s) == 2) then
        nTs = CIA_tab(s)%nT(1)
        nrec = CIA_tab(s)%nwl
        CIA_nset(s) = 1
        CIA_nT(1,s) = nTs
        CIA_irec(1,s) = nrec
        CIA_T_offset(1,s) = T_pos
        CIA_wn_offset(1,s) = wn_pos
        CIA_tab_offset(1,s) = tab_pos
        CIA_T_flat(T_pos:T_pos+nTs-1) = CIA_tab(s)%T(1,1:nTs)
        CIA_wl_flat(wn_pos:wn_pos+nrec-1) = CIA_tab(s)%wl(1,1:nrec)
        do j = 1, nTs
          CIA_tab_flat(tab_pos+(j-1)*nrec:tab_pos+j*nrec-1) = &
            & CIA_tab(s)%tab(1,1:nrec,j)
        end do
        T_pos = T_pos + nTs
        wn_pos = wn_pos + nrec
        tab_pos = tab_pos + nTs*nrec
      end if
    end do

    offload_active = 0
    !$omp target map(from: offload_active)
    if (.not. omp_is_initial_device()) offload_active = 1
    !$omp end target
    if (offload_active /= 1) then
      print*, 'ERROR - CIA OpenMP target region did not execute on a GPU - STOPPING'
      print*, 'OpenMP target devices visible: ', omp_get_num_devices()
      stop 1
    end if
    print*, ' ~~ CIA OpenMP GPU offload active; visible devices: ', omp_get_num_devices()
    print*, ' ~~ Performing CIA interpolation and output ~~ '
    print*, ' ~~ Please wait... ~~ '

    !$omp target data map(to: VMR_lay,N_lay,TG_lay,RH_lay,wl,wl_A,wn,freq,nCIA,nlay) &
    !$omp& map(to: CIA_form,CIA_model,CIA_nset,CIA_iVMR,CIA_iVMR_3,CIA_nT,CIA_irec) &
    !$omp& map(to: CIA_T_offset,CIA_wn_offset,CIA_tab_offset) &
    !$omp& map(to: CIA_wn_s,CIA_wn_e,CIA_Tmin,CIA_Tmax) &
    !$omp& map(to: CIA_T_flat,CIA_wn_flat,CIA_wl_flat,CIA_tab_flat) &
    !$omp& map(to: Hminus_ff1,Hminus_ff2,Hminus_bf) map(alloc: CIA_out)

    do l = 1, nwl
      if (mod(l,max(1,nwl/10)) == 0) then
        print*, l, wl(l), nwl
      end if

      !$omp target teams loop &
      !$omp& private(s,sn,j,n,nTs,nrec,iwn,iwn1,iT,iT1,jl,jm,ju) &
      !$omp& private(CIA_work,CIA_value,xval,yval,x0,x1,y0,y1) &
      !$omp& private(a00,a10,a01,a11,norm,lxval,lyval,lx0,lx1,ly0,ly1) &
      !$omp& private(la00,la10,la01,la11,T,T5040,a,b,c,kff,kbf,fbf,xbf,sff)
      do z = 1, nlay
        CIA_work = 0.0_dp

        do s = 1, nCIA
          if (CIA_form(s) /= 4) then
            CIA_value = 0.0_dp

            select case(CIA_model(s))
            case(CIA_MODEL_HMINUS)
              ! H- bound-free and free-free opacity from John (1988).
              T = TG_lay(z)
              T5040 = 5040.0_dp/T

              if (wl(l) > 1.6419_dp .or. wl(l) < 0.125_dp) then
                xbf = 0.0_dp
              else
                fbf = 0.0_dp
                do n = 1, 6
                  fbf = fbf + Hminus_bf(n) &
                    & * (1.0_dp/wl(l)-1.0_dp/1.6419_dp) &
                    & **((real(n,kind=dp)-1.0_dp)/2.0_dp)
                end do
                xbf = 1.0e-18_dp*wl(l)**3 &
                  & * (1.0_dp/wl(l)-1.0_dp/1.6419_dp)**(3.0_dp/2.0_dp)*fbf
              end if

              sff = 0.0_dp
              if (wl(l) >= 0.3645_dp) then
                do n = 1, 6
                  sff = sff + T5040**((real(n,kind=dp)+1.0_dp)/2.0_dp) &
                    & * (wl(l)**2*Hminus_ff2(n,1)+Hminus_ff2(n,2) &
                    & + Hminus_ff2(n,3)/wl(l)+Hminus_ff2(n,4)/wl(l)**2 &
                    & + Hminus_ff2(n,5)/wl(l)**3+Hminus_ff2(n,6)/wl(l)**4)
                end do
                kff = 1.0e-29_dp*sff
              else if (wl(l) > 0.1823_dp) then
                do n = 1, 6
                  sff = sff + T5040**((real(n,kind=dp)+1.0_dp)/2.0_dp) &
                    & * (wl(l)**2*Hminus_ff1(n,1)+Hminus_ff1(n,2) &
                    & + Hminus_ff1(n,3)/wl(l)+Hminus_ff1(n,4)/wl(l)**2 &
                    & + Hminus_ff1(n,5)/wl(l)**3+Hminus_ff1(n,6)/wl(l)**4)
                end do
                kff = 1.0e-29_dp*sff
              else
                kff = 0.0_dp
              end if

              kbf = xbf*VMR_lay(CIA_iVMR_3(1,s),z)*N_lay(z)
              kff = kff*(VMR_lay(CIA_iVMR_3(2,s),z)*N_lay(z) &
                & * VMR_lay(CIA_iVMR_3(3,s),z)*N_lay(z))*kb*T
              CIA_value = kbf + kff

            case(CIA_MODEL_HEMINUS)
              ! Analytic He- free-free opacity.
              T = TG_lay(z)
              a = 3.397e-46_dp + (-5.216e-31_dp+7.039e-15_dp/freq(l))/freq(l)
              b = -4.116e-42_dp + (1.067e-26_dp+8.135e-11_dp/freq(l))/freq(l)
              c = 5.081e-37_dp + (-8.724e-23_dp-5.659e-8_dp/freq(l))/freq(l)
              kff = a*T + b + c/T
              CIA_value = kff*VMR_lay(CIA_iVMR(1,s),z)*N_lay(z) &
                & * VMR_lay(CIA_iVMR(2,s),z)*N_lay(z)

            case(CIA_MODEL_HEMINUS_BELL,CIA_MODEL_H2MINUS_BELL)
              ! Bell He- or H2- table interpolation.
              nrec = CIA_irec(1,s)
              nTs = CIA_nT(1,s)
              xval = wl_A(l)
              yval = 5040.0_dp/TG_lay(z)

              jl = 0
              ju = nrec + 1
              do while (ju-jl > 1)
                jm = (ju+jl)/2
                if (xval > CIA_wl_flat(CIA_wn_offset(1,s)+jm-1)) then
                  jl = jm
                else
                  ju = jm
                end if
              end do
              iwn = jl
              iwn1 = iwn + 1

              jl = 0
              ju = nTs + 1
              do while (ju-jl > 1)
                jm = (ju+jl)/2
                if (yval > CIA_T_flat(CIA_T_offset(1,s)+jm-1)) then
                  jl = jm
                else
                  ju = jm
                end if
              end do
              iT = jl
              iT1 = iT + 1

              if (iwn >= 1 .and. iwn1 <= nrec .and. &
                  & iT >= 1 .and. iT1 <= nTs) then
                x0 = CIA_wl_flat(CIA_wn_offset(1,s)+iwn-1)
                x1 = CIA_wl_flat(CIA_wn_offset(1,s)+iwn1-1)
                y0 = CIA_T_flat(CIA_T_offset(1,s)+iT-1)
                y1 = CIA_T_flat(CIA_T_offset(1,s)+iT1-1)
                a00 = CIA_tab_flat(CIA_tab_offset(1,s)+(iT-1)*nrec+iwn-1)
                a10 = CIA_tab_flat(CIA_tab_offset(1,s)+(iT-1)*nrec+iwn1-1)
                a01 = CIA_tab_flat(CIA_tab_offset(1,s)+iT*nrec+iwn-1)
                a11 = CIA_tab_flat(CIA_tab_offset(1,s)+iT*nrec+iwn1-1)
                norm = 1.0_dp/(x1-x0)/(y1-y0)
                CIA_value = a00*(x1-xval)*(y1-yval)*norm &
                  & + a10*(xval-x0)*(y1-yval)*norm &
                  & + a01*(x1-xval)*(yval-y0)*norm &
                  & + a11*(xval-x0)*(yval-y0)*norm
                CIA_value = 1.0e-26_dp*CIA_value &
                  & * VMR_lay(CIA_iVMR(1,s),z)*N_lay(z) &
                  & * (VMR_lay(CIA_iVMR(2,s),z)*N_lay(z)*kb*TG_lay(z))
              end if

            case(CIA_MODEL_H2O)
              CIA_value = 1.0e-30_dp &
                & * (1.0_dp+4000.0_dp*(wl(l)-0.3_dp)/(10.0_dp-0.3_dp))*1.0e4_dp
              CIA_value = CIA_value*VMR_lay(CIA_iVMR(1,s),z)*N_lay(z)
            end select

            CIA_work = CIA_work + CIA_value
            cycle
          end if

          sn = 0
          if (CIA_nset(s) > 1) then
            do j = 1, CIA_nset(s)
              if (wn(l) >= CIA_wn_s(j,s) .and. wn(l) <= CIA_wn_e(j,s)) then
                if (TG_lay(z) >= CIA_Tmin(j,s) .and. &
                    & TG_lay(z) <= CIA_Tmax(j,s)) then
                  sn = j
                  exit
                end if
              end if
            end do

            if (sn == 0) then
              do j = 1, CIA_nset(s)
                if (wn(l) >= CIA_wn_s(j,s) .and. wn(l) <= CIA_wn_e(j,s)) then
                  if (TG_lay(z) <= CIA_Tmin(j,s) .or. &
                      & TG_lay(z) >= CIA_Tmax(j,s)) then
                    sn = j
                    exit
                  end if
                end if
              end do
            end if
            if (sn == 0) cycle
          else
            sn = 1
          end if

          nrec = CIA_irec(sn,s)
          nTs = CIA_nT(sn,s)
          xval = wn(l)

          ! Binary search for the lower wavenumber table index.
          jl = 0
          ju = nrec + 1
          do while (ju-jl > 1)
            jm = (ju+jl)/2
            if (xval > CIA_wn_flat(CIA_wn_offset(sn,s)+jm-1)) then
              jl = jm
            else
              ju = jm
            end if
          end do
          iwn = jl
          iwn1 = iwn + 1
          if (iwn < 1 .or. iwn1 > nrec) cycle

          ! Binary search for the lower temperature table index.
          yval = TG_lay(z)
          jl = 0
          ju = nTs + 1
          do while (ju-jl > 1)
            jm = (ju+jl)/2
            if (yval > CIA_T_flat(CIA_T_offset(sn,s)+jm-1)) then
              jl = jm
            else
              ju = jm
            end if
          end do
          iT = jl
          iT1 = iT + 1

          x0 = CIA_wn_flat(CIA_wn_offset(sn,s)+iwn-1)
          x1 = CIA_wn_flat(CIA_wn_offset(sn,s)+iwn1-1)

          if (iT < 1) then
            a00 = CIA_tab_flat(CIA_tab_offset(sn,s)+iwn-1)
            a10 = CIA_tab_flat(CIA_tab_offset(sn,s)+iwn1-1)
            norm = 1.0_dp/log10(x1/x0)
            CIA_value = 10.0_dp**((log10(a00)*log10(x1/xval) + &
              & log10(a10)*log10(xval/x0))*norm)
          else if (iT1 > nTs) then
            a00 = CIA_tab_flat(CIA_tab_offset(sn,s)+(nTs-1)*nrec+iwn-1)
            a10 = CIA_tab_flat(CIA_tab_offset(sn,s)+(nTs-1)*nrec+iwn1-1)
            norm = 1.0_dp/log10(x1/x0)
            CIA_value = 10.0_dp**((log10(a00)*log10(x1/xval) + &
              & log10(a10)*log10(xval/x0))*norm)
          else
            y0 = CIA_T_flat(CIA_T_offset(sn,s)+iT-1)
            y1 = CIA_T_flat(CIA_T_offset(sn,s)+iT1-1)
            a00 = CIA_tab_flat(CIA_tab_offset(sn,s)+(iT-1)*nrec+iwn-1)
            a10 = CIA_tab_flat(CIA_tab_offset(sn,s)+(iT-1)*nrec+iwn1-1)
            a01 = CIA_tab_flat(CIA_tab_offset(sn,s)+iT*nrec+iwn-1)
            a11 = CIA_tab_flat(CIA_tab_offset(sn,s)+iT*nrec+iwn1-1)

            lxval = log10(xval)
            lyval = log10(yval)
            lx0 = log10(x0)
            lx1 = log10(x1)
            ly0 = log10(y0)
            ly1 = log10(y1)
            la00 = log10(a00)
            la10 = log10(a10)
            la01 = log10(a01)
            la11 = log10(a11)
            norm = 1.0_dp/(lx1-lx0)/(ly1-ly0)
            CIA_value = la00*(lx1-lxval)*(ly1-lyval)*norm &
              & + la10*(lxval-lx0)*(ly1-lyval)*norm &
              & + la01*(lx1-lxval)*(lyval-ly0)*norm &
              & + la11*(lxval-lx0)*(lyval-ly0)*norm
            CIA_value = 10.0_dp**CIA_value
          end if

          CIA_work = CIA_work + CIA_value &
            & * VMR_lay(CIA_iVMR(1,s),z)*N_lay(z) &
            & * VMR_lay(CIA_iVMR(2,s),z)*N_lay(z)
        end do

        CIA_out(z) = CIA_work/RH_lay(z)
      end do
      !$omp end target teams loop

      !$omp target update from(CIA_out)
      if (any(.not. ieee_is_finite(CIA_out))) then
        print*, 'ERROR - Non-finite CIA opacity returned by GPU at wavelength: ', l, wl(l)
        stop 1
      end if
      call output_CIA_table(l)
    end do

    !$omp end target data

    deallocate(CIA_form,CIA_model,CIA_nset,CIA_iVMR,CIA_iVMR_3,CIA_nT,CIA_irec)
    deallocate(CIA_T_offset,CIA_wn_offset,CIA_tab_offset)
    deallocate(CIA_wn_s,CIA_wn_e,CIA_Tmin,CIA_Tmax)
    deallocate(CIA_T_flat,CIA_wn_flat,CIA_wl_flat,CIA_tab_flat)

  end subroutine calc_CIA_table_gpu
#endif

  subroutine validate_CIA_tables()
    implicit none

    integer :: s, sn, i

    do s = 1, nCIA
      if (CIA_tab(s)%form /= 4) cycle

      if (.not. allocated(CIA_tab(s)%nT) .or. .not. allocated(CIA_tab(s)%irec) .or. &
        & .not. allocated(CIA_tab(s)%T) .or. .not. allocated(CIA_tab(s)%wn) .or. &
        & .not. allocated(CIA_tab(s)%ltab) .or. .not. allocated(CIA_tab(s)%Tmin) .or. &
        & .not. allocated(CIA_tab(s)%Tmax) .or. .not. allocated(CIA_tab(s)%wn_s) .or. &
        & .not. allocated(CIA_tab(s)%wn_e)) then
        print*, 'ERROR - CIA table did not provide all required interpolation data - STOPPING'
        print*, 'Species, path: ', CIA_tab(s)%sp, trim(CIA_tab(s)%path)
        stop
      end if

      do sn = 1, CIA_tab(s)%nset
        if (CIA_tab(s)%nT(sn) < 1) then
          print*, 'ERROR - CIA table set must contain at least one temperature - STOPPING'
          print*, 'Species, set, nT: ', CIA_tab(s)%sp, sn, CIA_tab(s)%nT(sn)
          stop
        end if

        if (CIA_tab(s)%irec(sn) < 3) then
          print*, 'ERROR - CIA Bezier interpolation requires at least 3 wavenumber points - STOPPING'
          print*, 'Species, set, nwn: ', CIA_tab(s)%sp, sn, CIA_tab(s)%irec(sn)
          stop
        end if

        if (any(.not. ieee_is_finite(CIA_tab(s)%T(sn,1:CIA_tab(s)%nT(sn)))) .or. &
          & any(.not. ieee_is_finite(CIA_tab(s)%wn(sn,1:CIA_tab(s)%irec(sn)))) .or. &
          & any(.not. ieee_is_finite(CIA_tab(s)%ltab(sn,1:CIA_tab(s)%irec(sn), &
          & 1:CIA_tab(s)%nT(sn))))) then
          print*, 'ERROR - CIA interpolation data must be finite - STOPPING'
          print*, 'Species, set, path: ', CIA_tab(s)%sp, sn, trim(CIA_tab(s)%path)
          stop
        end if

        if (any(CIA_tab(s)%T(sn,1:CIA_tab(s)%nT(sn)) <= 0.0_dp)) then
          print*, 'ERROR - CIA temperature grid must be positive - STOPPING'
          print*, 'Species, set: ', CIA_tab(s)%sp, sn
          stop
        end if

        do i = 2, CIA_tab(s)%nT(sn)
          if (CIA_tab(s)%T(sn,i) <= CIA_tab(s)%T(sn,i-1)) then
            print*, 'ERROR - CIA temperature grid must be strictly increasing - STOPPING'
            print*, 'Species, set, index, values: ', CIA_tab(s)%sp, sn, i, &
              & CIA_tab(s)%T(sn,i-1), CIA_tab(s)%T(sn,i)
            stop
          end if
        end do

        do i = 2, CIA_tab(s)%irec(sn)
          if (CIA_tab(s)%wn(sn,i) <= CIA_tab(s)%wn(sn,i-1)) then
            print*, 'ERROR - CIA wavenumber grid must be strictly increasing - STOPPING'
            print*, 'Species, set, index, values: ', CIA_tab(s)%sp, sn, i, &
              & CIA_tab(s)%wn(sn,i-1), CIA_tab(s)%wn(sn,i)
            stop
          end if
        end do
      end do
    end do

  end subroutine validate_CIA_tables

  subroutine find_CIA_consituents()
    implicit none

    integer :: s

    do s = 1, nCIA

      select case(CIA_tab(s)%sp)

      case('H2-H2')
        CIA_tab(s)%sp_con(1) = 'H2'
        CIA_tab(s)%sp_con(2) = 'H2'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 113

      case('H2-He','He-H2')
        CIA_tab(s)%sp_con(1) = 'H2'
        CIA_tab(s)%sp_con(2) = 'He'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 334

      case('H-')
        CIA_tab(s)%sp_con_3(1) = 'H-'
        CIA_tab(s)%sp_con_3(2) = 'H'
        CIA_tab(s)%sp_con_3(3) = 'e-'

        CIA_tab(s)%i3 = .True.

      case('He-')
        CIA_tab(s)%sp_con(1) = 'He'
        CIA_tab(s)%sp_con(2) = 'e-'

      case('H2-')
        CIA_tab(s)%sp_con(1) = 'H2'
        CIA_tab(s)%sp_con(2) = 'e-'

      case('H2-H','H-H2')
        CIA_tab(s)%sp_con(1) = 'H2'
        CIA_tab(s)%sp_con(2) = 'H'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 4

      case('H-He','He-H')
        CIA_tab(s)%sp_con(1) = 'He'
        CIA_tab(s)%sp_con(2) = 'H'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 10

      case('CO2-CO2')
        CIA_tab(s)%sp_con(1) = 'CO2'
        CIA_tab(s)%sp_con(2) = 'CO2'

        CIA_tab(s)%nset = 3
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 10
        CIA_tab(s)%nT(2) = 6
        CIA_tab(s)%nT(3) = 3
        !CIA_tab(s)%nT(4) = 1

      case('CO2-He','He-CO2')
        CIA_tab(s)%sp_con(1) = 'CO2'
        CIA_tab(s)%sp_con(2) = 'He'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 1

      case('CO2-H2','H2-CO2')
        CIA_tab(s)%sp_con(1) = 'CO2'
        CIA_tab(s)%sp_con(2) = 'H2'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 4

      case('CO2-Ar','Ar-CO2')
        CIA_tab(s)%sp_con(1) = 'CO2'
        CIA_tab(s)%sp_con(2) = 'Ar'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 21

      case('N2-H2O','H2O-N2')
        CIA_tab(s)%sp_con(1) = 'N2'
        CIA_tab(s)%sp_con(2) = 'H2O'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 21

      case('N2-H2','H2-N2')
        CIA_tab(s)%sp_con(1) = 'N2'
        CIA_tab(s)%sp_con(2) = 'H2'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 10

      case('N2-He','He-N2')
        CIA_tab(s)%sp_con(1) = 'N2'
        CIA_tab(s)%sp_con(2) = 'He'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 1
   
      case('N2-N2')
        CIA_tab(s)%sp_con(1) = 'N2'
        CIA_tab(s)%sp_con(2) = 'N2'

        CIA_tab(s)%nset = 6
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 14
        CIA_tab(s)%nT(2) = 10
        CIA_tab(s)%nT(3) = 10
        CIA_tab(s)%nT(4) = 5
        CIA_tab(s)%nT(5) = 5
        CIA_tab(s)%nT(6) = 14

      case('N2-CH4','CH4-N2')
        CIA_tab(s)%sp_con(1) = 'N2'
        CIA_tab(s)%sp_con(2) = 'CH4'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 10

      case('H2O-H2O')
        CIA_tab(s)%sp_con(1) = 'H2O'
        CIA_tab(s)%sp_con(2) = 'H2O'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 21

      case('CH4-CH4')
        CIA_tab(s)%sp_con(1) = 'CH4'
        CIA_tab(s)%sp_con(2) = 'CH4'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 10

      case('CH4-He','He-CH4')
        CIA_tab(s)%sp_con(1) = 'CH4'
        CIA_tab(s)%sp_con(2) = 'He'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 10

      case('H2-CH4','CH4-H2')
        CIA_tab(s)%sp_con(1) = 'H2'
        CIA_tab(s)%sp_con(2) = 'CH4'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 10

      case('CH4-Ar','Ar-CH4')
        CIA_tab(s)%sp_con(1) = 'CH4'
        CIA_tab(s)%sp_con(2) = 'Ar'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 5

      case('O2-CO2','CO2-O2')
        CIA_tab(s)%sp_con(1) = 'O2'
        CIA_tab(s)%sp_con(2) = 'CO2'

        CIA_tab(s)%nset = 1
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 1

      case('O2-N2','N2-O2')
        CIA_tab(s)%sp_con(1) = 'O2'
        CIA_tab(s)%sp_con(2) = 'N2'

        CIA_tab(s)%nset = 5
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 7
        CIA_tab(s)%nT(2) = 5
        CIA_tab(s)%nT(3) = 5
        CIA_tab(s)%nT(4) = 1
        CIA_tab(s)%nT(5) = 1

      case('O2-O2')
        CIA_tab(s)%sp_con(1) = 'O2'
        CIA_tab(s)%sp_con(2) = 'O2'

        CIA_tab(s)%nset = 8
        allocate(CIA_tab(s)%nT(CIA_tab(s)%nset))
        CIA_tab(s)%nT(1) = 15
        CIA_tab(s)%nT(2) = 1
        CIA_tab(s)%nT(3) = 1
        CIA_tab(s)%nT(4) = 1
        CIA_tab(s)%nT(5) = 1
        CIA_tab(s)%nT(6) = 1
        CIA_tab(s)%nT(7) = 4
        CIA_tab(s)%nT(8) = 5

      case('H2O')
        CIA_tab(s)%sp_con(1) = 'H2O'
        CIA_tab(s)%sp_con(2) = 'H2O'

      case default
        print*, 'ERROR - CIA species constituents could not be found - STOPPING'
        print*, 'Species: ', CIA_tab(s)%sp
        stop
      end select

    end do

  end subroutine find_CIA_consituents

  subroutine output_CIA_table(l)
    implicit none

    integer, intent(in) :: l
    integer :: reclen

    if (first_call .eqv. .True.) then
      !print*, 'Outputting CIA.cmcrt'
      inquire(iolength=reclen) CIA_write
      ! Output k-table in 1D or flattened 3D CMCRT format k_CMCRT.ktb (single precision)
      open(newunit=uCIA, file='CIA.cmcrt', action='readwrite', &
      & form='unformatted',status='replace',access='direct',recl=reclen)
      first_call = .False.
    end if

    ! Convert to single precision on output and protect against underflow.
    CIA_write(:) = real(max(CIA_out(:),1.0e-30_dp),kind=sp)
    write(uCIA,rec=l) CIA_write

  end subroutine output_CIA_table

end module CIA_tables_mod
