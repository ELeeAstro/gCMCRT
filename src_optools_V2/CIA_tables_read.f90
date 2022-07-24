module CIA_tables_read
  use optools_data_mod
  implicit none


  ! NEMESIS format parameters and arrays
  ! integer, parameter :: NUMPAIRS = 9, NUMT = 25, NUMWN = 1501
  ! real(kind=dp), parameter :: dnu1 = 20.0_dp ! = 20.0 for Jean-Loup's .tab and 10.0 for the rest of NEMESIS
  ! real(kind=dp), dimension(NUMT) :: temps
  ! real(kind=dp), allocatable, dimension(:) :: temps2
  ! integer :: itemps2
  ! real(kind=sp), dimension(NUMT) :: tempk1
  ! real(kind=sp), dimension(NUMPAIRS,NUMT,NUMWN) :: kcia
  ! real(kind=dp), dimension(NUMWN) :: vvk1, wlk1

  private :: read_CIA_HITRAN, read_CIA_NEMESIS, read_CIA_Bell
  public :: read_CIA_tables

contains

  subroutine read_CIA_tables()
    implicit none

    integer :: s

    print*, ' ~~ Reading in CIA tables ~~ '

    do s = 1, nCIA


      select case (CIA_tab(s)%form)
      case(0)
        print*, ' - Skipping special species: ', s,  CIA_tab(s)%sp, CIA_tab(s)%form

      case(1)
        print*, ' - Reading NEMESIS CIA table: ', s,  CIA_tab(s)%sp, CIA_tab(s)%form
        call read_CIA_NEMESIS(s)

      case(2)
        print*, ' - Reading Bell CIA table: ', s,  CIA_tab(s)%sp, CIA_tab(s)%form
        call read_CIA_Bell(s)

      case(4)
        print*, ' - Reading HITRAN CIA table for species: ', s, CIA_tab(s)%sp, CIA_tab(s)%form
        call read_CIA_HITRAN(s)
      case default
        print*, 'ERROR - CIA table format integer not valid - STOPPING'
        print*, 'Species: ', s, CIA_tab(s)%sp, CIA_tab(s)%form
        stop
      end select
    end do

    print*, ' ~~ Quest completed  ~~ '

  end subroutine read_CIA_tables

  subroutine read_CIA_HITRAN(s)
    use :: iso_fortran_env
    implicit none

    integer, intent(in) :: s
    integer :: u, i, n, j
    integer :: stat

    character(len=20) :: name
    integer :: nrec, mnT, mnrec
    real(kind=dp) :: wn_s, wn_e, temp_r, kmax, dum

    open(newunit=u,file=trim(CIA_tab(s)%path),form='formatted', status='old',action='read')

    ! Allocate CIA table temperature arrays
    mnT = maxval(CIA_tab(s)%nT(:))
    allocate(CIA_tab(s)%T(CIA_tab(s)%nset,mnT),CIA_tab(s)%lT(CIA_tab(s)%nset,mnT))
    allocate(CIA_tab(s)%irec(CIA_tab(s)%nset))
    allocate(CIA_tab(s)%Tmin(CIA_tab(s)%nset))
    allocate(CIA_tab(s)%Tmax(CIA_tab(s)%nset))
    allocate(CIA_tab(s)%wn_s(CIA_tab(s)%nset))
    allocate(CIA_tab(s)%wn_e(CIA_tab(s)%nset))

    ! Read to find meta data first

    do j = 1, CIA_tab(s)%nset
      do n = 1, CIA_tab(s)%nT(j)
        read(u,*,iostat=stat) name, wn_s, wn_e, nrec, temp_r, kmax, dum
        if (n == 1) then
          CIA_tab(s)%irec(j) = nrec
          CIA_tab(s)%wn_s(j) = wn_s
          CIA_tab(s)%wn_e(j) = wn_e
          CIA_tab(s)%Tmin(j) = temp_r
        else if (n == CIA_tab(s)%nT(j)) then
          CIA_tab(s)%Tmax(j) = temp_r
        end if
        CIA_tab(s)%T(j,n) = temp_r
        CIA_tab(s)%lT(j,n) = log10(CIA_tab(s)%T(j,n))
        ! Check if end of file reached
        if (stat == iostat_end) then
          print*,'Reached end of HITRAN CIA file: ', j, n, CIA_tab(s)%sp, CIA_tab(s)%path
        end if
        do i = 1, nrec
          read(u,*)
        end do
      end do
    end do

    rewind(u)

    mnrec = maxval(CIA_tab(s)%irec(:))
    allocate(CIA_tab(s)%wn(CIA_tab(s)%nset,mnrec))
    allocate(CIA_tab(s)%tab(CIA_tab(s)%nset,mnrec,mnT),CIA_tab(s)%ltab(CIA_tab(s)%nset,mnrec,mnT))

    ! Read and allocate data until error (end of file)
    do j = 1, CIA_tab(s)%nset
      do n = 1, CIA_tab(s)%nT(j)
        read(u,*)
        ! Read the record data
        do i = 1, CIA_tab(s)%irec(j)
          read(u,*) CIA_tab(s)%wn(j,i), CIA_tab(s)%tab(j,i,n)
          CIA_tab(s)%tab(j,i,n) = max(abs(CIA_tab(s)%tab(j,i,n)),1e-99_dp) 
          CIA_tab(s)%ltab(j,i,n) = log10(CIA_tab(s)%tab(j,i,n))
          !print*, i, CIA_tab(s)%wn(j,i), CIA_tab(s)%tab(j,i,n)
        end do
      end do
    end do

    close(u)

  end subroutine read_CIA_HITRAN

  subroutine read_CIA_Bell(s)
    implicit none

    integer, intent(in) :: s

    integer :: u, n

    allocate(CIA_tab(s)%nT(1))
    CIA_tab(s)%nT = 8

    if (trim(CIA_tab(s)%sp) == 'He-') then
      CIA_tab(s)%nwl = 16
    else if (trim(CIA_tab(s)%sp) == 'H2-') then
      CIA_tab(s)%nwl = 18
    end if

    open(newunit=u,file=trim(CIA_tab(s)%path),action='read',status='old')
    read(u,*)

    allocate(CIA_tab(s)%wl(1,CIA_tab(s)%nwl))
    allocate(CIA_tab(s)%wn(1,CIA_tab(s)%nwl))
    allocate(CIA_tab(s)%T(1,CIA_tab(s)%nT(1)))
    allocate(CIA_tab(s)%tab(1,CIA_tab(s)%nwl,CIA_tab(s)%nT(1)))

    read(u,*) CIA_tab(s)%T(1,:)

    do n = 1, CIA_tab(s)%nwl
      read(u,*) CIA_tab(s)%wl(1,n), CIA_tab(s)%tab(1,n,:)
      !print*, s, n, CIA(s)%wl(n), CIA(s)%kap(n,:)
    end do
    !CIA(s)%wl(:) = CIA(s)%wl(:) * 1.0e-4_dp
    CIA_tab(s)%wn(1,:) = 1.0_dp/(CIA_tab(s)%wl(1,:) * 1.0e-8_dp)

    close(u)

  end subroutine read_CIA_Bell

  subroutine read_CIA_NEMESIS(s)
    implicit none
  !
    integer, intent(in) :: s
  !
  !   open(newunit=u,file=trim(CIA_tab(s)%path),action='read',status='old')
  !
  !   do i = 1, NUMT
  !     tempk1(i) = real(temps(i),kind=sp)
  !   end do
  !
  !   read(u) kcia
  !
  !   do i = 1, NUMWN
  !    vvk1(i) = (i-1)*dnu1
  !   end do
  !
  !  close(u)
  !
  !  ! Remove 0 temperatures from temperature grid
  !  itemps2 = 0
  !  do i = 1, NUMT
  !    if (temps(i) == 0.0_dp) then
  !      itemps2 = itemps2 + 1
  !    end if
  !  end do
  !
  !  itemps2 = NUMT - itemps2
  !  allocate(temps2(itemps2))
  !  temps2(:) = temps(1:itemps2)
  !
  end subroutine read_CIA_NEMESIS


end module CIA_tables_read
