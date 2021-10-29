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

  private :: read_CIA_HITRAN, read_CIA_NEMESIS
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
        exit

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

  end subroutine

  subroutine read_CIA_HITRAN(s)
    use :: iso_fortran_env
    implicit none

    integer, intent(in) :: s
    integer :: u, i, n
    integer :: stat

    character(len=20) :: name
    integer :: nrec
    real(kind=dp) :: wn_s, wn_e, temp_r, kmax, dum

    open(newunit=u,file=trim(CIA_tab(s)%path),form='formatted', status='old',action='read')

    ! Allocate CIA table temperature arrays
    allocate(CIA_tab(s)%T(CIA_tab(s)%nT))

    ! Read and allocate data until error (end of file)
    do n = 1, CIA_tab(s)%nT
      read(u,*,iostat=stat) name, wn_s, wn_e, nrec, temp_r, kmax, dum
      !print*, n, name, wn_s, wn_e, nrec, temp_r, kmax, dum
      if (n == 1) then
        ! Allocate CIA table wn and table value array
        allocate(CIA_tab(s)%wn(nrec))
        allocate(CIA_tab(s)%tab(nrec,CIA_tab(s)%nT))
        CIA_tab(s)%nwl = nrec
      end if
      ! Check if end of file reached
      if (stat == iostat_end) then
        print*,'Reached end of HITRAN CIA file: ', CIA_tab(s)%sp, CIA_tab(s)%path
        exit
      else
        ! Temperature point of table
        CIA_tab(s)%T(n) = temp_r
        ! Read the record data
        do i = 1, nrec
          read(u,*) CIA_tab(s)%wn(i), CIA_tab(s)%tab(i,n)
          !print*, i, CIA_tab(s)%wn(i), CIA_tab(s)%tab(i,n)
        end do
      end if

    end do

    close(u)

  end subroutine read_CIA_HITRAN

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
