program gpuCMCRT
  use mc_data_mod
  use mc_class_imag
  implicit none

  integer, parameter :: len_max = 10
  integer :: i , nargs
  character (len_max) :: arg, arg2

  print*, 'Command line arguments: '
  nargs = command_argument_count()
  do i = 0, nargs
    call get_command_argument (i, arg)
    if (arg == '-vphi') then
      call get_command_argument (i+1, arg2)
      vphi_arg = arg2
      print*, trim(arg)//' '//trim(vphi_arg)
      read(vphi_arg , *) im%vphi
      cmd_vphi = .True.
      if (im%vphi > 360.0_dp .or. im%vphi < 0.0_dp) then
        print*, 'Invalid vphi detected: ', im%vphi
        stop
      end if

    end if
  end do

  print*, 'Reading namelist'
  call read_namelist()

  print*, 'Starting Experiment'

  select case(trim(xper))
  case('1D_pp')
    !call exp_1D_pp_atm()
  case('1D_sph')
    !call exp_1D_sph_atm()
  case('2D_car_gal')
    !call exp_3D_cart_galaxy()
  case('3D_sph_tests')
    call exp_3D_sph_atm()
  case('3D_sph_pol')
    call exp_3D_sph_atm_pol()
  case('3D_sph_alb')
    call exp_3D_sph_atm_albedo()
  case('3D_sph_trans')
    call exp_3D_sph_atm_transmission()
  case('3D_sph_em')
    call exp_3D_sph_atm_em()
  case('3D_sph_em_hi')
    call exp_3D_sph_atm_em_hires()
  case('3D_sph_trans_hi')
    call exp_3D_sph_atm_trans_hires()
  case default
    print*, 'Invalid experiment selected: ', trim(xper)
  end select

end program gpuCMCRT
