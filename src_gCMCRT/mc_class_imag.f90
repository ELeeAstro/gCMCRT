module mc_class_imag
  use mc_precision
  use mc_data_mod
  use cudafor
  implicit none

  type imag
    integer :: im_id
    integer :: x_pix, y_pix
    real(dp) :: sinto,costo,sinpo,cospo,phio
    real(dp) :: obsx,obsy,obsz
    real(dp) :: vtheta, vphi
    real(dp) :: rimage
    real(dp) :: fsum, qsum, usum
    integer :: fail_pscat
    integer :: fail_pemit
  end type imag

  type(imag) :: im             ! Host image data
  type(imag), device :: im_d   ! Device image data

  real(dp), allocatable, dimension(:,:) :: f, q, u, im_err
  real(dp), allocatable, dimension(:,:), device :: f_d, q_d, u_d, im_err_d

  integer, allocatable, dimension(:) :: f_im, q_im, u_im

contains

  subroutine set_image()
    implicit none
    logical, save :: first_call = .True.

    !print*, "setting image"


    !! Observation direction and image set up
    im%sinto = sin(im%vtheta*pi/180.0_dp)
    im%costo = cos(im%vtheta*pi/180.0_dp)
    im%sinpo = sin(im%vphi*pi/180.0_dp)
    im%cospo = cos(im%vphi*pi/180.0_dp)
    im%phio = atan2(im%sinpo,im%cospo)
    im%obsx = im%sinto*im%cospo
    im%obsy = im%sinto*im%sinpo
    im%obsz = im%costo


    if (do_images .eqv. .True. .and. first_call .eqv. .True.) then
      im%x_pix = xpix
      im%y_pix = ypix
      im%rimage = rimage
      ! Calculate in CPU
      allocate(f(im%x_pix,im%y_pix),q(im%x_pix,im%y_pix),u(im%x_pix,im%y_pix))
      allocate(im_err(im%x_pix,im%y_pix))

      !! Give to the GPU
      do_images_d = do_images
      allocate(f_d(im%x_pix,im%y_pix),q_d(im%x_pix,im%y_pix),u_d(im%x_pix,im%y_pix))
      allocate(im_err_d(im%x_pix,im%y_pix))

      first_call = .False.
    end if

    !print*, ' - Complete - '


  end subroutine set_image

  subroutine output_im(n,l)
    implicit none

    integer, intent(in) :: l, n
    logical, save :: first_call = .True.
    integer :: nn
    character (len=8) :: fmt
    character (len=3) :: n_str

    if (first_call .eqv. .True.) then
      allocate(f_im(n_phase))
      allocate(q_im(n_phase))
      allocate(u_im(n_phase))
      fmt = '(I3.3)'
      do nn = 1, n_phase
        write(n_str,fmt) nn 
        open(newunit=f_im(nn), file='f_im_'//trim(n_str)//'.txt', action='readwrite',form='unformatted')
      end do
      !open(newunit=q_im, file='q_im.txt', action='readwrite',form='unformatted')
      !open(newunit=u_im, file='u_im.txt', action='readwrite',form='unformatted')
      first_call = .False.
    end if

    write(f_im(n)) real(f(:,:))
    !write(q_im) real(q(:,:))
    !write(u_im) real(u(:,:))

    !flush(f_im)
    !flush(q_im)
    !flush(u_im)

  end subroutine output_im


end module mc_class_imag
