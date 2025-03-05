!!! E.K.H. Lee -  Translation into fortran 90 of the c++ code LX-MIE from Kitzmann & Heng (2018)
!!! original comments have been retained
!!! - input : ri - refractive index (with negative k, m = n - k), x size parameter
!!! - output : q_ext, qsca, q_abs efficencies and asymmetry parameter (g)
!!! - 
module lxmie_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  integer, parameter :: cf_max_terms = 10000001
  real(dp), parameter :: cf_epsilon = 1.0e-10_dp

  public :: lxmie
  private :: calcMieCoefficients, calcMieEfficiencies, calcAsymmetryParameter, &
    & an, anReal, startingANcontinuedFractions, startingANcontinuedFractionsReal

contains

  subroutine lxmie(ri, x, q_ext, q_sca, q_abs, g)
    implicit none

    complex(dp), intent(in) :: ri
    real(dp), intent(in) :: x

    real(dp), intent(out) :: q_ext, q_sca, q_abs, g

    integer :: nb
    complex(dp), allocatable, dimension(:) :: mie_coeff_a, mie_coeff_b

    !! Calculate the maximum number of terms in the Mie series
    !! See Eq. 22
    nb = int(x + 4.3_dp * x**(1.0_dp/3.0_dp)) + 2

    !! Allocate the mie coffiecent arrays
    allocate(mie_coeff_a(nb),mie_coeff_b(nb))

    call calcMieCoefficients(nb, ri, x, mie_coeff_a, mie_coeff_b)
    call calcMieEfficiencies(nb, x, mie_coeff_a, mie_coeff_b, q_ext, q_sca, q_abs)
    call calcAsymmetryParameter(nb, q_sca, x, mie_coeff_a, mie_coeff_b, g)

    deallocate(mie_coeff_a,mie_coeff_b)

  end subroutine lxmie


  !! Calculates the Mie coefficients a and b that are used for the construction of the Mie series, see Eq. 17
  !! The required coefficients A_n are calculated via downward recursion, B_n and C_n by upward recursion
  !! a and b are evaluated up to a number of "nb_mie_terms"
  subroutine calcMieCoefficients(nb, ri, x, mie_coeff_a, mie_coeff_b)
    implicit none

    integer, intent(in) :: nb
    complex(dp), intent(in) :: ri
    real(dp), intent(in) :: x

    complex(dp), dimension(nb), intent(out) :: mie_coeff_a, mie_coeff_b

    integer :: n
    real(dp) :: dn
    complex(dp) :: m, mx
    complex(dp) :: C_n, D_n, A_n
    real(dp) :: A_n_r

    m = ri
    mx = m * x

    !! First, calculate A_n via backward recursion
    !! Note that we need A_N(mx), a complex number, and A_N(x), a real number (see Eq. 17)
    !! We use the Mie a-coefficient to store the complex values and the Mie b-coefficients for the real numbers
    mie_coeff_a(nb) = startingANcontinuedFractions(nb, mx)
    mie_coeff_b(nb)%re = startingANcontinuedFractionsReal(nb, x)

    !! backward recursion
    do n = nb, 2, -1
      dn = real(n,dp)
      mie_coeff_a(n-1) = dn/mx - 1.0_dp/(dn/mx + mie_coeff_a(n))
      mie_coeff_b(n-1)%re = dn/x - 1.0_dp/(dn/x + real(mie_coeff_b(n),dp))
    end do

    !! Now we do a forward recursion to calculate B_n, C_n, and the Mie coefficients a_n and b_n
    C_n = cmplx(0.0_dp,0.0_dp,dp)
    D_n = cmplx(0.0_dp,-1.0_dp,dp)

    !! n = 1
    C_n = cmplx(1.0_dp, (cos(x) + x*sin(x))/(sin(x) - x*cos(x)), dp)
    C_n = 1.0_dp/C_n
    D_n = -1.0_dp/x + 1.0_dp/(1.0_dp/x - D_n)

    A_n = mie_coeff_a(1)
    A_n_r = real(mie_coeff_b(1),dp)

    mie_coeff_a(1) = C_n * (A_n/m - A_n_r) / (A_n/m - D_n)
    mie_coeff_b(1) = C_n * (A_n*m - A_n_r) / (A_n*m - D_n)

    !! n > 1
    do n = 2, nb

      dn = real(n,dp)
      A_n = mie_coeff_a(n)
      A_n_r = real(mie_coeff_b(n),dp)

      D_n = -dn/x + 1.0_dp/(dn/x - D_n)
      C_n = C_n * (D_n + dn/x)/(A_n_r + dn/x)

      mie_coeff_a(n) = C_n * (A_n/m - A_n_r) / (A_n/m - D_n)
      mie_coeff_b(n) = C_n * (A_n*m - A_n_r) / (A_n*m - D_n)

    end do

  end subroutine calcMieCoefficients

  !! Calculates the Mie efficiencies, see Eq. 1
  !! The absorption efficiency is calculated as the difference of the extinction and scattering efficiencies
  subroutine calcMieEfficiencies(nb, x, mie_coeff_a, mie_coeff_b, q_ext, q_sca, q_abs)
    implicit none

    integer, intent(in) :: nb
    real(dp), intent(in) :: x
    complex(dp), dimension(nb+1), intent(in) :: mie_coeff_a, mie_coeff_b

    real(dp), intent(out) :: q_ext, q_sca, q_abs

    integer :: n
    real(dp) :: dn

    q_ext = 0.0_dp
    q_sca = 0.0_dp

    do n = 1, nb
      dn = real(n,dp)
      q_sca = q_sca + (2.0_dp*dn + 1.0_dp) * (abs(mie_coeff_a(n)) * abs(mie_coeff_a(n)) + abs(mie_coeff_b(n))*abs(mie_coeff_b(n)))
      q_ext = q_ext + (2.0_dp*dn + 1.0_dp) * real(mie_coeff_a(n) + mie_coeff_b(n),dp)
    end do

    q_sca = q_sca*(2.0_dp/x**2)
    q_ext = q_ext*(2.0_dp/x**2)

    q_abs = q_ext - q_sca

  end subroutine calcMieEfficiencies

  !! Calculate and return the asymmetry parameter
  !! See Bohren&Huffman, page 120, for details on the equation
  subroutine calcAsymmetryParameter(nb, q_sca, x, mie_coeff_a, mie_coeff_b, g)
    implicit none

    integer, intent(in) :: nb
    real(dp), intent(in) :: q_sca, x
    complex(dp), dimension(nb+1), intent(in) :: mie_coeff_a, mie_coeff_b

    real(dp), intent(out) :: g

    integer :: n
    real(dp) :: dn

    g = 0.0_dp

    do n = 1, nb-1
      dn = real(n,dp)
      g = g + dn*(dn + 2.0_dp)/(dn + 1.0_dp) * real(mie_coeff_a(n) * conjg(mie_coeff_a(n+1)) + &
        & mie_coeff_b(n)*conjg(mie_coeff_b(n+1)),dp) + &
        & (2.0_dp*dn + 1.0_dp)/(dn * (dn + 1.0_dp)) * real(mie_coeff_a(n)*conjg(mie_coeff_b(n)),dp)
    end do 

    g = g * (4.0_dp/(x**2 * q_sca))

  end subroutine calcAsymmetryParameter

  !! Calculate the starting value of A_N via the method of continued fractions by Lentz (1976)
  !! Convergence is reached if two consecutive terms differ by less than "continued_fraction_epsilon"
  !! Returns A_N
  !! This is the version for a complex A_N 
  complex(dp) function startingANcontinuedFractions(nb, mx)
    implicit none

    integer, intent(in) :: nb
    complex(dp), intent(in) :: mx

    integer :: i
    real(dp) :: nu, con
    complex(dp) :: a_i
    complex(dp) :: function_numerator, function_denominator
    complex(dp) :: a_numerator, a_denominator

    nu = real(nb,dp) + 0.5_dp

    !! starting values
    function_numerator = cmplx(1.0_dp,1.0_dp,dp)
    function_denominator = cmplx(1.0_dp,1.0_dp,dp)

    !! n = 1
    a_numerator = an(1, nu, mx)
    a_denominator = 1.0_dp

    function_numerator = function_numerator * a_numerator
    function_denominator = function_denominator * a_denominator

    !! n = 2
    a_numerator = an(2, nu, mx) + 1.0_dp/a_numerator
    a_denominator = an(2, nu, mx)

    function_numerator = function_numerator * a_numerator
    function_denominator = function_denominator * a_denominator    

    do i = 3, cf_max_terms

      a_i = an(i, nu, mx)

      a_numerator = a_i + 1.0_dp/a_numerator
      a_denominator = a_i + 1.0_dp/a_denominator

      function_numerator = function_numerator * a_numerator
      function_denominator = function_denominator * a_denominator    

      con = abs((abs(a_numerator) - abs(a_denominator))/abs(a_numerator))
      if (con < cf_epsilon) then
        exit
      end if
    end do

    startingANcontinuedFractions = function_numerator/function_denominator &
      & - real(nb,dp)/mx

  end function startingANcontinuedFractions

  !! Calculate the starting value of A_N via the method of continued fractions by Lentz (1976)
  !! Convergence is reached if two consecutive terms differ by less than "continued_fraction_epsilon"
  !! Returns A_N
  !! This is the version for a real A_N
  real(dp) function startingANcontinuedFractionsReal(nb, x)
    implicit none

    integer, intent(in) :: nb
    real(dp), intent(in) :: x

    integer :: i
    real(dp) :: nu, con, a_i
    real(dp) :: function_numerator, function_denominator
    real(dp) :: a_numerator, a_denominator

    nu = real(nb,dp) + 0.5_dp

    !! starting values
    function_numerator = 1.0_dp
    function_denominator = 1.0_dp

    !! n = 1
    a_numerator = anReal(1, nu, x)
    a_denominator = 1.0_dp

    function_numerator = function_numerator * a_numerator
    function_denominator = function_denominator * a_denominator

    !! n = 2
    a_numerator = anReal(2, nu, x) + 1.0_dp/a_numerator
    a_denominator = anReal(2, nu, x)

    function_numerator = function_numerator * a_numerator
    function_denominator = function_denominator * a_denominator    

    do i = 3, cf_max_terms

      a_i = anReal(i, nu, x)

      a_numerator = a_i + 1.0_dp/a_numerator
      a_denominator = a_i + 1.0_dp/a_denominator

      function_numerator = function_numerator * a_numerator
      function_denominator = function_denominator * a_denominator

      con = abs((a_numerator - a_denominator)/a_numerator)
      if (con < cf_epsilon) then
        exit
      end if
    end do

    startingANcontinuedFractionsReal = function_numerator/function_denominator &
      & - real(nb,dp)/x

  end function startingANcontinuedFractionsReal

  !! Calculates single terms a_n used in the evaluation of the continued fractions by Lentz (1976)
  !! See Eq. 9 in Lentz (1976)
  !! This is the version for a complex a_n 
  complex(dp) function an(n, nu, z)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(in) :: nu
    complex(dp), intent(in) :: z

    if (mod(n,2) == 0) then
      an = -1.0_dp * 2.0_dp * (nu + real(n,dp) - 1.0_dp) * 1.0_dp/z
    else
      an = 2.0_dp * (nu + real(n,dp) - 1.0_dp) * 1.0_dp/z
    end if

  end function an

  !! Calculates single terms a_n used in the evaluation of the continued fractions by Lentz (1976)
  !! See Eq. 9 in Lentz (1976)
  !! This is the version for a real-numbered a_n
  real(dp) function anReal(n, nu, z)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(in) :: nu
    real(dp), intent(in) :: z

    if (mod(n,2) == 0) then
      anReal = -1.0_dp * 2.0_dp * (nu + real(n,dp) - 1.0_dp) * 1.0_dp/z
    else
      anReal = 2.0_dp * (nu + real(n,dp) - 1.0_dp) * 1.0_dp/z
    end if

  end function anReal

end module lxmie_mod
