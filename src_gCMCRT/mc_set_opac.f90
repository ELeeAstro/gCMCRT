module mc_opacset
  use mc_precision
  use mc_data_mod
  use mc_class_grid
  use mc_k_aux
  use ieee_arithmetic
  use mc_Draine_G, only : Draine_G
  implicit none

  logical :: first_call = .True.

  integer :: u_k, u_conti, u_Ray, u_cld_k, u_cld_a, u_cld_g
  integer :: id_u_k, id_u_conti, id_u_Ray, id_u_cld_k, id_u_cld_a, id_u_cld_g
  integer :: dop_pad

  !! Grid opacities
  real(dp), dimension(:,:,:,:), allocatable :: k_gas_abs, k_tot_abs
  real(dp), dimension(:,:,:), allocatable ::  k_gas_Ray

  integer :: iwl_rest
  real(dp), dimension(:), allocatable :: wl_pad
  real(dp), dimension(:,:,:,:), allocatable :: k_gas_abs_dop, k_tot_abs_dop
  real(dp), dimension(:,:,:,:), allocatable ::  k_gas_Ray_dop

  ! Dummy variable for reading in k-tables
  real(sp), allocatable, dimension(:) :: k_dum
  real(sp), allocatable, dimension(:) :: conti_dum_arr, Ray_dum_arr, lbl_dum_arr
  real(sp), allocatable, dimension(:) :: cld_k_dum_arr, cld_a_dum_arr, cld_g_dum_arr
  real(sp), allocatable, dimension(:,:) :: ck_dum_arr


contains

  subroutine set_grid_opac()
    implicit none

    integer :: k, j, i, g
    real(dp) :: k_tot_ext(ng), k_tot_scat, cld_sca


    if (oneD .eqv. .True.) then
      ! If 1D, only need to do 1 vertical column for whole atmosphere and g loop for indexing
      do i = 1, grid%n_lay
          cld_sca = cld_ssa(i,1,1)*cld_ext(i,1,1)
          k_tot_abs(:,i,1,1) = k_gas_abs(:,i,1,1) + (1.0_dp - cld_ssa(i,1,1))*cld_ext(i,1,1)
          k_tot_scat = k_gas_Ray(i,1,1) + cld_sca
          k_tot_ext(:) = k_tot_abs(:,i,1,1) +  k_gas_Ray(i,1,1) + cld_ext(i,1,1)
          rhokap(:,i,1,1) =  RH(i,1,1) * k_tot_ext(:) * grid%r_del
          ssa(:,i,1,1) = min(k_tot_scat/k_tot_ext(:), 0.99_dp)
          gg(i,1,1) = cld_g(i,1,1)
          dorg(i,1,1) = k_gas_Ray(i,1,1)/k_tot_scat
      end do

      ! Copy relvent array data to 3D
      do k = 1, grid%n_theta-1
        do j = 1, grid%n_phi-1
          k_tot_abs(:,:,j,k) = k_tot_abs(:,:,1,1)
          rhokap(:,:,j,k) = rhokap(:,:,1,1)
          ssa(:,:,j,k) = ssa(:,:,1,1)
          gg(:,j,k) = gg(:,1,1)
          dorg(:,j,k) = dorg(:,1,1)
        end do
      end do

    else if (threeD .eqv. .True.) then
      do k = 1, grid%n_theta-1
        do j = 1, grid%n_phi-1
           do i = 1, grid%n_lay
             cld_sca = cld_ssa(i,j,k)*cld_ext(i,j,k)
             k_tot_abs(:,i,j,k) = k_gas_abs(:,i,j,k) + (1.0_dp - cld_ssa(i,j,k))*cld_ext(i,j,k)
             k_tot_scat = k_gas_Ray(i,j,k) +  cld_sca
             k_tot_ext(:) = k_tot_abs(:,i,j,k) +  k_gas_Ray(i,j,k) + cld_ext(i,j,k)
             rhokap(:,i,j,k) =  RH(i,j,k) * k_tot_ext(:) * grid%r_del
             ssa(:,i,j,k) = min(k_tot_scat/k_tot_ext(:), 0.99_dp)
             gg(i,j,k) = cld_g(i,j,k)!(cld_g(i,j,k)*cld_sca) / (cld_sca + k_gas_Ray(i,j,k))
             dorg(i,j,k) = k_gas_Ray(i,j,k)/k_tot_scat
           end do
        end do
      end do
  end if

    !! Send relvant arrays to device memory
    rhokap_d(:,:,:,:) = rhokap(:,:,:,:)
    ssa_d(:,:,:,:) = ssa(:,:,:,:)
    g_d(:,:,:) = gg(:,:,:)
    dorg_d(:,:,:) = dorg(:,:,:)

  end subroutine set_grid_opac

  subroutine read_next_opac(l)
    implicit none

    integer, intent(in) :: l
    integer :: NX_dum, n_bins_dum, ng_dum, z, g, j, k, n
    integer :: reclen
    real(sp) :: Ray_dum, k_lbl_dum, conti_dum

    if (l == n_wl+1) then
      return
    end if

    if (first_call .eqv. .True.) then

      if (ck .eqv. .True.) then

        ! Allocate dummy variable for reading in k-table
        allocate(k_dum(ng))
        allocate(ck_dum_arr(ng,grid%n_cell))
        inquire(iolength=reclen) ck_dum_arr

        ! Read k-table in 1D or 3D CMCRT format (single precision)
        open(newunit=u_k, file='CK.cmcrt', status='old', action='read',&
         & form='unformatted', asynchronous='yes', access='direct',recl=reclen)

        print*, 'unit _k :', u_k
        print*, '- Complete -'

      else if (lbl .eqv. .True.) then

        allocate(lbl_dum_arr(grid%n_cell))
        inquire(iolength=reclen) lbl_dum_arr

        ! Read lbl table in 1D or 3D CMCRT format (single precision)
        open(newunit=u_k, file='lbl.cmcrt', status='old', action='read',&
        & form='unformatted', asynchronous='yes',access='direct',recl=reclen)

        print*, 'unit _lbl :', u_k
        print*, '- Complete -'

      end if

      if (inc_CIA .eqv. .True.) then

        allocate(conti_dum_arr(grid%n_cell))
        inquire(iolength=reclen) conti_dum_arr

        open(newunit=u_conti, file='CIA.cmcrt', status='old', action='read', &
          & form='unformatted', asynchronous='yes', access='direct',recl=reclen)

        print*, 'unit _conti :', u_conti
        print*, '- Complete -'

      end if

      if (inc_Ray .eqv. .True.) then

        allocate(Ray_dum_arr(grid%n_cell))
        inquire(iolength=reclen) Ray_dum_arr

        open(newunit=u_Ray, file='Rayleigh.cmcrt', status='old', action='read', &
         & form='unformatted', asynchronous='yes', access='direct',recl=reclen)

        print*, 'unit _Ray :', u_Ray
        print*, '- Complete -'

      end if

      if (inc_cld .eqv. .True.) then

        allocate(cld_k_dum_arr(grid%n_cell))
        allocate(cld_g_dum_arr(grid%n_cell))
        allocate(cld_a_dum_arr(grid%n_cell))
        inquire(iolength=reclen) cld_k_dum_arr

        open(newunit=u_cld_k, file='cl_k.cmcrt', status='old', action='read', &
          & form='unformatted', asynchronous='yes', access='direct',recl=reclen)

        print*, 'unit _cld_k :', u_cld_k

        open(newunit=u_cld_a, file='cl_a.cmcrt', status='old', action='read', &
          & form='unformatted', asynchronous='yes', access='direct',recl=reclen)

        print*, 'unit _cld_a :', u_cld_a

        open(newunit=u_cld_g, file='cl_g.cmcrt', status='old', action='read', &
          & form='unformatted', asynchronous='yes', access='direct',recl=reclen)

        print*, 'unit _cld_g :', u_cld_g
        print*, '- Complete -'

      end if

      if (inc_lbl .eqv. .True.) then
        ng = 1
        ng_d = ng
      end if

      allocate(k_gas_abs(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(k_gas_Ray(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(k_tot_abs(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))

      allocate(cld_ext(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(cld_ssa(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(cld_g(grid%n_lay,grid%n_phi-1,grid%n_theta-1))

      allocate(dorg(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(dorg_d(grid%n_lay,grid%n_phi-1,grid%n_theta-1))

      allocate(rhokap(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(rhokap_d(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))

      allocate(ssa(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(ssa_d(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))

      allocate(gg(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(g_d(grid%n_lay,grid%n_phi-1,grid%n_theta-1))

      first_call = .False.


    end if


    k_gas_abs(:,:,:,:) = 0.0_dp
    k_gas_Ray(:,:,:) = 0.0_dp
    cld_ext(:,:,:) = 0.0_dp
    cld_ssa(:,:,:) = 0.0_dp
    cld_g(:,:,:) = 0.0_dp

    if (oneD .eqv. .True.) then

      if (inc_ck .eqv. .True.) then
        read(u_k,rec=l) ck_dum_arr
      else if (lbl .eqv. .True.) then
        read(u_k,rec=l) lbl_dum_arr
      end if
      if (inc_cld .eqv. .True.) then
        read(u_cld_k,rec=l) cld_k_dum_arr
        read(u_cld_a,rec=l) cld_a_dum_arr
        read(u_cld_g,rec=l) cld_g_dum_arr
      end if
      if (inc_CIA .eqv. .True.) then
        read(u_conti,rec=l) conti_dum_arr
      end if
      if (inc_Ray .eqv. .True.) then
       read(u_Ray,rec=l) Ray_dum_arr
      end if

      if (inc_ck .eqv. .True.) then
        wait(u_k)
        do z = 1, grid%n_lay
            do g = 1, ng
              k_gas_abs(g,z,:,:) = k_gas_abs(g,z,:,:) + real(ck_dum_arr(g,z),dp)
            end do
        end do
      else if (lbl .eqv. .True.) then
        wait(u_k)
        do z = 1, grid%n_lay
            k_gas_abs(1,z,:,:) = k_gas_abs(1,z,:,:) + real(lbl_dum_arr(z),dp)
        end do
      end if

      if (inc_cld .eqv. .True.) then
        wait(u_cld_k)
        wait(u_cld_a)
        wait(u_cld_g)

        do z = 1, grid%n_lay
          cld_ext(z,:,:) = real(cld_k_dum_arr(z),dp)
          cld_ssa(z,:,:) = real(cld_a_dum_arr(z),dp)
          cld_g(z,:,:) = real(cld_g_dum_arr(z),dp)
        end do
      end if

      if (inc_CIA .eqv. .True.) then
        wait(u_conti)
        do z = 1, grid%n_lay
          k_gas_abs(:,z,:,:) = k_gas_abs(:,z,:,:) + real(conti_dum_arr(z),dp)
        end do
      end if

      if (inc_Ray .eqv. .True.) then
        wait(u_Ray)
        do z = 1, grid%n_lay
          k_gas_Ray(z,:,:) = real(Ray_dum_arr(z),dp)
        end do
      end if

   else if (threeD .eqv. .True.) then

     if (inc_ck .eqv. .True.) then
       read(u_k,rec=l) ck_dum_arr
     else if (lbl .eqv. .True.) then
       read(u_k,rec=l) lbl_dum_arr
     end if
     if (inc_cld .eqv. .True.) then
       read(u_cld_k,rec=l) cld_k_dum_arr
       read(u_cld_a,rec=l) cld_a_dum_arr
       read(u_cld_g,rec=l) cld_g_dum_arr
     end if
     if (inc_CIA .eqv. .True.) then
       read(u_conti,rec=l) conti_dum_arr
     end if
     if (inc_Ray .eqv. .True.) then
       read(u_Ray,rec=l) Ray_dum_arr
     end if

      !print*, 'ck reading 1 : ', ck_dum_arr(:,1)!, ck_dum_arr(ng,1)
      !print*, 'ck reading n_cell: ', ck_dum_arr(:,grid%n_cell)

     if (inc_ck .eqv. .True.) then
       wait(u_k)
       n = 1
       do k = 1, grid%n_theta-1
         do j = 1, grid%n_phi-1
           do z = 1, grid%n_lay
             k_gas_abs(:,z,j,k) = k_gas_abs(:,z,j,k) + real(ck_dum_arr(:,n),dp)
             n = n + 1
          end do
         end do
       end do

     else if (lbl .eqv. .True.) then
       wait(u_k)
       !print*, 'lbl reading: ', lbl_dum_arr(1), lbl_dum_arr(grid%n_cell)
       n = 1
       do k = 1, grid%n_theta-1
         do j = 1, grid%n_phi-1
           do z = 1, grid%n_lay
             k_gas_abs(1,z,j,k) = k_gas_abs(1,z,j,k) + real(lbl_dum_arr(n),dp)
             n = n + 1
           end do
         end do
       end do

     end if

     if (inc_cld .eqv. .True.) then
       wait(u_cld_k)
       wait(u_cld_a)
       wait(u_cld_g)
      !
      ! print*, 'cld_k reading: ', cld_k_dum_arr(1), cld_k_dum_arr(grid%n_cell)
      ! print*, 'cld_a reading: ', cld_a_dum_arr(1), cld_a_dum_arr(grid%n_cell)
      ! print*, 'cld_g reading: ', cld_g_dum_arr(1), cld_g_dum_arr(grid%n_cell)

      n = 1
      do k = 1, grid%n_theta-1
        do j = 1, grid%n_phi-1
          do z = 1, grid%n_lay
            cld_ext(z,j,k) = real(cld_k_dum_arr(n),dp)
            cld_ssa(z,j,k) = real(cld_a_dum_arr(n),dp)
            cld_g(z,j,k) = real(cld_g_dum_arr(n),dp)
            !print*, k,j,z, k_gas_abs(ng,z,j,k), cld_k_dum_arr(n), cld_a_dum_arr(n), cld_g_dum_arr(n)

            n = n + 1
         end do
        end do
      end do
     end if

     if (inc_CIA .eqv. .True.) then
       wait(u_conti)
       !print*, 'Conti reading: ', conti_dum_arr(1), conti_dum_arr(grid%n_cell)
       n = 1
        do k = 1, grid%n_theta-1
          do j = 1, grid%n_phi-1
            do z = 1, grid%n_lay
              k_gas_abs(:,z,j,k) = k_gas_abs(:,z,j,k) + real(conti_dum_arr(n),dp)
              n = n + 1
           end do
          end do
        end do
    end if

    if (inc_Ray .eqv. .True.) then
      wait(u_Ray)
       !print*, 'Ray reading: ', Ray_dum_arr(1), Ray_dum_arr(grid%n_cell)
       n = 1
       do k = 1, grid%n_theta-1
         do j = 1, grid%n_phi-1
           do z = 1, grid%n_lay
             k_gas_Ray(z,j,k) = real(Ray_dum_arr(n),dp)
             n = n + 1
          end do
         end do
       end do
    end if

   end if

   if (do_Draine .eqv. .True.) then
     call Draine_G()
   end if

  end subroutine read_next_opac


  subroutine read_next_opac_doppler(ll)
    implicit none

    integer, intent(in) :: ll
    integer :: NX_dum, n_bins_dum, ng_dum, z, g, j, k, n, l
    real(sp) :: Ray_dum, k_lbl_dum, conti_dum
    integer :: pad_idx1, pad_idx2, reclen
    real(dp) :: wl_test1, wl_test2

    if (ll == n_wl+1) then
      return
    end if


    if (first_call .eqv. .True.) then

      allocate(lbl_dum_arr(grid%n_cell))
      inquire(iolength=reclen) lbl_dum_arr

      ! Read lbl table in 1D or 3D CMCRT format (single precision)
      open(newunit=u_k, file='lbl.cmcrt', status='old', action='read', &
        & form='unformatted', asynchronous='yes', access='direct',recl=reclen)

      print*, 'unit _lbl :', u_k
      print*, '- Complete -'

      allocate(conti_dum_arr(grid%n_cell))
      inquire(iolength=reclen) conti_dum_arr

      open(newunit=u_conti, file='CIA.cmcrt', status='old', action='read', &
        & form='unformatted', asynchronous='yes', access='direct',recl=reclen)

      print*, 'unit _conti :', u_conti
      print*, '- Complete -'

      if (inc_Ray .eqv. .True.) then

        allocate(Ray_dum_arr(grid%n_cell))
        inquire(iolength=reclen) Ray_dum_arr

        open(newunit=u_Ray, file='Rayleigh.cmcrt', status='old', action='read', &
          & form='unformatted', asynchronous='yes', access='direct',recl=reclen)

        print*, 'unit _Ray :', u_Ray
        print*, '- Complete -'

      end if

      ng = 1
      ng_d = ng

      if (allocated(k_gas_abs) .eqv. .False.) then
        wl_test1 = wl(1) * (1.0_dp - minval(v_los(:,:,:,:))/c_s)
        wl_test2 = wl(n_wl) * (1.0_dp - maxval(v_los(:,:,:,:))/c_s)
        call locate_host(wl,wl_test1,pad_idx1)
        call locate_host(wl,wl_test2,pad_idx2)
        dop_pad = (max(pad_idx1, n_wl-pad_idx2) + 1) * 2 ! Ensure dop_pad is an even number
        !! Allocate doppler padded arrays
        allocate(wl_pad(dop_pad))
        allocate(k_gas_abs_dop(dop_pad,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
        allocate(k_gas_Ray_dop(dop_pad,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
        allocate(k_tot_abs_dop(dop_pad,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
        !! Allocate normal lbl arrays
        allocate(k_gas_abs(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
        allocate(k_gas_Ray(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
        allocate(k_tot_abs(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
        allocate(cld_ext(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
        allocate(cld_ssa(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
        allocate(cld_g(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      end if

      allocate(dorg(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(dorg_d(grid%n_lay,grid%n_phi-1,grid%n_theta-1))

      allocate(rhokap(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(rhokap_d(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))

      allocate(ssa(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(ssa_d(ng,grid%n_lay,grid%n_phi-1,grid%n_theta-1))

      allocate(gg(grid%n_lay,grid%n_phi-1,grid%n_theta-1))
      allocate(g_d(grid%n_lay,grid%n_phi-1,grid%n_theta-1))

      !! Now read the first dop_pad opacities into memory

      k_gas_abs(:,:,:,:) = 0.0_dp
      k_gas_Ray(:,:,:) = 0.0_dp
      cld_ext(:,:,:) = 0.0_dp
      cld_ssa(:,:,:) = 0.0_dp
      cld_g(:,:,:) = 0.0_dp

      print*, 'reading doppler padding'

      k_gas_abs_dop(:,:,:,:) = 0.0_dp
      k_gas_Ray_dop(:,:,:,:) = 0.0_dp

      do l = 1, dop_pad
        print*, l, wl(l), dop_pad

        wl_pad(l) = wl(l)

        if (oneD .eqv. .True.) then

          do z = 1, grid%n_lay
            read(u_k,rec=ll+(l-1)) k_lbl_dum
            k_gas_abs_dop(l,z,:,:) = k_gas_abs_dop(l,z,:,:) + real(k_lbl_dum,dp)
          end do

           do z = 1, grid%n_lay
              read(u_conti,rec=ll+(l-1)) conti_dum
              k_gas_abs_dop(l,z,:,:) = k_gas_abs_dop(l,z,:,:) + real(conti_dum,dp)
           end do

           do z = 1, grid%n_lay
              read(u_Ray,rec=ll+(l-1)) Ray_dum
              k_gas_Ray_dop(l,z,:,:) = k_gas_Ray_dop(l,z,:,:) + real(Ray_dum,dp)
           end do

           else if (threeD .eqv. .True.) then

             read(u_k,rec=ll+(l-1)) lbl_dum_arr(:)
             n = 1
             do k = 1, grid%n_theta-1
               do j = 1, grid%n_phi-1
                 do z = 1, grid%n_lay
                   k_gas_abs_dop(l,z,j,k) = k_gas_abs_dop(l,z,j,k) + real(lbl_dum_arr(n),dp)
                   n = n + 1
                 end do
               end do
             end do


           read(u_conti,rec=ll+(l-1)) conti_dum_arr(:)
           n = 1
           do k = 1, grid%n_theta-1
             do j = 1, grid%n_phi-1
               do z = 1, grid%n_lay
                 k_gas_abs_dop(l,z,j,k) = k_gas_abs_dop(l,z,j,k) + real(conti_dum_arr(n),dp)
                 n = n + 1
               end do
             end do
           end do

           if (inc_Ray .eqv. .True.) then
             read(u_Ray,rec=ll+(l-1)) Ray_dum_arr(:)
             n = 1
             do k = 1, grid%n_theta-1
               do j = 1, grid%n_phi-1
                 do z = 1, grid%n_lay
                   k_gas_Ray_dop(l,z,j,k) = k_gas_Ray_dop(l,z,j,k) + real(Ray_dum_arr(n),dp)
                   n = n + 1
                 end do
              end do
            end do
          end if

        end if

      end do

      iwl_rest = 1

      first_call = .False.

    end if

   !! Shift opacities, unless dop_pad/2 within the wavelength grid edges

   if (ll <= dop_pad/2 .or. ll >= (n_wl - dop_pad/2)) then
     ! We don't need to read any more data in - only increase ticker
   else
     ! We are at the center of padded wavelength region - shift everything in preparation for next loop
     ! EKH-Lee note: This next loop causes a lot of overhead - try find an improvement
     !do l = 1, dop_pad-1
       !wl_pad(l) = wl_pad(l+1)
       !k_gas_abs_dop(l,:,:,:) = k_gas_abs_dop(l+1,:,:,:)
       !k_gas_Ray_dop(l,:,:,:) = k_gas_Ray_dop(l+1,:,:,:)
     !end do


     wl_pad(:) = eoshift(wl_pad(:), shift = 1, boundary = wl(ll + dop_pad/2), dim = 1)
     !wl_pad(1:dop_pad-1) = wl_pad(2:dop_pad)
     !wl_pad(dop_pad) = wl(ll + dop_pad/2) ! wl(ll + dop_pad/2 + 1) ! New end wavelength is + dop_pad/2 +1 indexes ahead

     k_gas_abs_dop(:,:,:,:) = eoshift(k_gas_abs_dop(:,:,:,:), shift = 1, boundary = 0.0_dp, dim = 1)
     !k_gas_abs_dop(1:dop_pad-1,:,:,:) = k_gas_abs_dop(2:dop_pad,:,:,:)
     !k_gas_abs_dop(dop_pad,:,:,:) = 0.0_dp

     if (inc_Ray .eqv. .True.) then
       k_gas_Ray_dop(:,:,:,:) = eoshift(k_gas_Ray_dop(:,:,:,:), shift = 1, boundary = 0.0_dp, dim = 1)
       !k_gas_Ray_dop(1:dop_pad-1,:,:,:) = k_gas_Ray_dop(2:dop_pad,:,:,:)
       !k_gas_Ray_dop(dop_pad,:,:,:) = 0.0_dp
     end if

     if (oneD .eqv. .True.) then

       do z = 1, grid%n_lay
         read(u_k,rec=ll+dop_pad/2) k_lbl_dum
         k_gas_abs_dop(dop_pad,:,:,:) = k_gas_abs_dop(dop_pad,:,:,:) + real(k_lbl_dum,dp)
       end do

       do z = 1, grid%n_lay
         read(u_conti,rec=ll+dop_pad/2) conti_dum
         k_gas_abs_dop(dop_pad,:,:,:) = k_gas_abs_dop(dop_pad,:,:,:) + real(conti_dum,dp)
       end do

       do z = 1, grid%n_lay
         read(u_Ray,rec=ll+dop_pad/2) Ray_dum
         k_gas_Ray_dop(dop_pad,:,:,:) = k_gas_Ray_dop(dop_pad,:,:,:) + real(Ray_dum,dp)
       end do

     else if (threeD .eqv. .True.) then

       if (lbl .eqv. .True.) then
         read(u_k,rec=ll+dop_pad/2) lbl_dum_arr
       end if
       if (inc_CIA .eqv. .True.) then
         read(u_conti,rec=ll+dop_pad/2) conti_dum_arr
       end if
       if (inc_Ray .eqv. .True.) then
         read(u_Ray,rec=ll+dop_pad/2) Ray_dum_arr
       end if

     if (lbl .eqv. .True.) then
       wait(u_k)
         n = 1
         do k = 1, grid%n_theta-1
           do j = 1, grid%n_phi-1
             do z = 1, grid%n_lay
               k_gas_abs_dop(dop_pad,z,j,k) = k_gas_abs_dop(dop_pad,z,j,k) + real(lbl_dum_arr(n),dp)
               n = n + 1
             end do
           end do
         end do
      end if

      if (inc_CIA .eqv. .True.) then
        wait(u_conti)
       n = 1
       do k = 1, grid%n_theta-1
         do j = 1, grid%n_phi-1
           do z = 1, grid%n_lay
             k_gas_abs_dop(dop_pad,z,j,k) = k_gas_abs_dop(dop_pad,z,j,k) + real(conti_dum_arr(n),dp)
             n = n + 1
           end do
         end do
       end do
     end if

       if (inc_Ray .eqv. .True.) then
         wait(u_Ray)
         n = 1
         do k = 1, grid%n_theta-1
           do j = 1, grid%n_phi-1
             do z = 1, grid%n_lay
               k_gas_Ray_dop(dop_pad,z,j,k) = k_gas_Ray_dop(dop_pad,z,j,k) + real(Ray_dum_arr(n),dp)
               n = n + 1
             end do
          end do
        end do
       end if

     end if

   end if

   ! Increase the iwl_rest ticker
   iwl_rest = iwl_rest + 1

 end subroutine read_next_opac_doppler

  subroutine shift_opac(n,ll)
    implicit none

    integer, intent(in) :: n, ll
    integer :: k, j, z
    integer :: wl_idx, wl_idx1
    real(dp) :: y1, y2, yval, wl_eff
    logical :: outr

    outr = .False.

     !! Now find the shifted opacity by interpolating from the dop array block
     do k = 1, grid%n_theta-1
       do j = 1, grid%n_phi-1
         do z = 1, grid%n_lay

           !! Find effective wavelength for the line of sight velocity
           wl_eff = wl(ll)*(1.0_dp - v_los(n,z,j,k)/c_s)
           call locate_host(wl_pad(:),wl_eff,wl_idx)
           wl_idx1 = wl_idx + 1

           !! Do some error checking
           if (wl_idx == 0) then
             ! blue shift is out of wavelength bounds, use rest frame opacity
             if (outr .eqv. .False.) then 
               print*, 'blueshift out of range: ', z, j, k, wl_eff, wl_pad(1),ll,dop_pad 
               outr = .True.
             end if
             k_gas_abs(1,z,j,k) = k_gas_abs_dop(1,z,j,k)
             k_gas_Ray(z,j,k) = k_gas_Ray_dop(1,z,j,k)
             cycle
           else if (wl_idx == dop_pad) then
             if (outr .eqv. .False.) then
               print*, 'redshift out of range: ', z, j, k, wl_eff, wl_pad(dop_pad),ll,dop_pad
               outr = .True.
             end if
             k_gas_abs(1,z,j,k) = k_gas_abs_dop(dop_pad,z,j,k)
             k_gas_Ray(z,j,k) = k_gas_Ray_dop(dop_pad,z,j,k)
             cycle
           end if

           !! Interpolate dop arrays to find abs
           y1 = k_gas_abs_dop(wl_idx,z,j,k)
           y2 = k_gas_abs_dop(wl_idx1,z,j,k)
           call linear_log_interp(wl_eff, wl_pad(wl_idx), wl_pad(wl_idx1), y1, y2, yval)
           k_gas_abs(1,z,j,k) = yval

           if (inc_Ray .eqv. .True.) then
             !! Interpolate dop arrays to find Ray
             y1 = k_gas_Ray_dop(wl_idx,z,j,k)
             y2 = k_gas_Ray_dop(wl_idx1,z,j,k)
             call linear_log_interp(wl_eff, wl_pad(wl_idx), wl_pad(wl_idx1), y1, y2, yval)
             k_gas_Ray(z,j,k) = yval
           else
             k_gas_Ray(z,j,k) = 0.0_dp
           end if

           if (k_gas_abs(1,z,j,k) < 0.0_dp .or. ieee_is_nan(k_gas_abs(1,z,j,k)) .or. .not.ieee_is_finite(k_gas_abs(1,z,j,k))) then
             print*, 'Encountered negative k_gas value --> check opacities'
             print*, z, j, k, k_gas_abs(1,z,j,k), wl_eff, wl_pad(wl_idx), wl_pad(wl_idx1), k_gas_abs_dop(wl_idx,z,j,k), k_gas_abs_dop(wl_idx1,z,j,k)
             stop
           end if

         end do
      end do
    end do
    
  end subroutine shift_opac

end module mc_opacset
