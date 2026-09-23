module gas_refractivity_mod
  use optools_data_mod, only : dp, pi, amu, N_stp
  implicit none

  private

  integer, parameter, public :: REFRAC_MODEL_UNAVAILABLE = 0
  integer, parameter, public :: REFRAC_MODEL_REFERENCE_INDEX = 1
  integer, parameter, public :: REFRAC_MODEL_POLARIZABILITY = 2
  integer, parameter, public :: REFRAC_MODEL_WATER = 3

  public :: gas_refractivity_data, water_lorentz_lorenz

contains

  subroutine gas_refractivity_data(name, wavelength_um, wavenumber_cm, model, &
      n_reference, number_density_reference, alpha_volume)
    implicit none

    character(len=*), intent(in) :: name
    real(dp), intent(in) :: wavelength_um, wavenumber_cm
    integer, intent(out) :: model
    real(dp), intent(out) :: n_reference, number_density_reference
    real(dp), intent(out) :: alpha_volume

    real(dp) :: A, B, C

    model = REFRAC_MODEL_UNAVAILABLE
    n_reference = 1.0_dp
    number_density_reference = 0.0_dp
    alpha_volume = 0.0_dp

    select case(trim(name))
    case('H2')
      A = 13.58e-5_dp
      B = 7.52e-3_dp
      n_reference = 1.0_dp + A*(1.0_dp + B/wavelength_um**2)
      number_density_reference = 2.65163e19_dp
      model = REFRAC_MODEL_REFERENCE_INDEX

    case('He')
      A = 2283.0_dp
      B = 1.8102e13_dp
      C = 1.5342e10_dp
      n_reference = refractive_index_wavenumber(wavenumber_cm,A,B,C)
      number_density_reference = 2.546899e19_dp
      model = REFRAC_MODEL_REFERENCE_INDEX

    case('CO')
      A = 22851.0_dp
      B = 0.456e14_dp
      C = 71427.0_dp**2
      n_reference = refractive_index_wavenumber(wavenumber_cm,A,B,C)
      number_density_reference = 2.546899e19_dp
      model = REFRAC_MODEL_REFERENCE_INDEX

    case('CO2')
      n_reference = 1.1427e3_dp * &
        ((5799.25_dp/(128908.9_dp**2-wavenumber_cm**2)) + &
         (120.05_dp/(89223.8_dp**2-wavenumber_cm**2)) + &
         (5.3334_dp/(75037.5_dp**2-wavenumber_cm**2)) + &
         (4.3244_dp/(67837.7_dp**2-wavenumber_cm**2)) + &
         (0.1218145e-6_dp/(2418.136_dp**2-wavenumber_cm**2)))
      n_reference = n_reference + 1.0_dp
      number_density_reference = 2.546899e19_dp
      model = REFRAC_MODEL_REFERENCE_INDEX

    case('CH4')
      n_reference = 1.0_dp + &
        (46662.0_dp + 4.02e-6_dp*wavenumber_cm**2)/1.0e8_dp
      number_density_reference = 2.546899e19_dp
      model = REFRAC_MODEL_REFERENCE_INDEX

    case('H2O')
      ! The water model consumes the local partial density and temperature.
      model = REFRAC_MODEL_WATER

    case('O2')
      A = 20564.8_dp
      B = 2.480899e13_dp
      C = 4.09e9_dp
      n_reference = refractive_index_wavenumber(wavenumber_cm,A,B,C)
      number_density_reference = 2.68678e19_dp
      model = REFRAC_MODEL_REFERENCE_INDEX

    case('N2')
      if (wavenumber_cm > 21360.0_dp) then
        A = 5677.465_dp
        B = 318.81874e12_dp
      else
        A = 6498.2_dp
        B = 307.4335e12_dp
      end if
      C = 14.4e9_dp
      n_reference = refractive_index_wavenumber(wavenumber_cm,A,B,C)
      number_density_reference = 2.546899e19_dp
      model = REFRAC_MODEL_REFERENCE_INDEX

    case('NH3')
      A = 37.0e-5_dp
      B = 12.0e-3_dp
      n_reference = 1.0_dp + A*(1.0_dp + B/wavelength_um**2)
      number_density_reference = N_stp
      model = REFRAC_MODEL_REFERENCE_INDEX

    case('Ar')
      A = 6432.135_dp
      B = 286.06021e12_dp
      C = 14.4e9_dp
      n_reference = refractive_index_wavenumber(wavenumber_cm,A,B,C)
      number_density_reference = 2.546899e19_dp
      model = REFRAC_MODEL_REFERENCE_INDEX

    case('N2O')
      n_reference = 1.0_dp + &
        (46890.0_dp + 4.12e-6_dp*wavenumber_cm**2)/1.0e8_dp
      number_density_reference = 2.546899e19_dp
      model = REFRAC_MODEL_REFERENCE_INDEX

    case('SF6')
      n_reference = 1.0_dp + &
        (71517.0_dp + 4.996e-6_dp*wavenumber_cm**2)/1.0e8_dp
      number_density_reference = 2.546899e19_dp
      model = REFRAC_MODEL_REFERENCE_INDEX

    case('HCl')
      alpha_volume = 2.515_dp*1.0e-24_dp
      model = REFRAC_MODEL_POLARIZABILITY
    case('HCN')
      alpha_volume = 2.593_dp*1.0e-24_dp
      model = REFRAC_MODEL_POLARIZABILITY
    case('H2S')
      alpha_volume = 3.631_dp*1.0e-24_dp
      model = REFRAC_MODEL_POLARIZABILITY
    case('OCS')
      alpha_volume = 5.090_dp*1.0e-24_dp
      model = REFRAC_MODEL_POLARIZABILITY
    case('SO2')
      alpha_volume = 3.882_dp*1.0e-24_dp
      model = REFRAC_MODEL_POLARIZABILITY
    case('C2H2')
      alpha_volume = 3.487_dp*1.0e-24_dp
      model = REFRAC_MODEL_POLARIZABILITY
    case('PH3')
      alpha_volume = 4.237_dp*1.0e-24_dp
      model = REFRAC_MODEL_POLARIZABILITY
    case('SO3')
      alpha_volume = 4.297_dp*1.0e-24_dp
      model = REFRAC_MODEL_POLARIZABILITY
    end select

  end subroutine gas_refractivity_data


  subroutine water_lorentz_lorenz(wavelength_um, temperature, number_density, &
      lorentz_term, valid)
    implicit none

    real(dp), intent(in) :: wavelength_um, temperature, number_density
    real(dp), intent(out) :: lorentz_term
    logical, intent(out) :: valid

    real(dp), parameter :: molecular_weight_h2o = 18.01528_dp
    real(dp), parameter :: rho_s = 1000.0_dp
    real(dp), parameter :: lambda_s = 0.589_dp
    real(dp), parameter :: temperature_s = 273.15_dp
    real(dp), parameter :: lambda_ir = 5.432937_dp
    real(dp), parameter :: lambda_uv = 0.229202_dp
    real(dp), parameter :: coeff(8) = (/ &
      0.244257733_dp, 0.974634476e-2_dp, -0.373234996e-2_dp, &
      0.268678472e-3_dp, 0.158920570e-2_dp, 0.245934259e-2_dp, &
      0.900704920_dp, -0.166626219e-1_dp /)

    real(dp) :: density_h2o, sigma, lambda_ratio, theta

    lorentz_term = 0.0_dp
    valid = .False.

    if (number_density < 0.0_dp .or. temperature <= 0.0_dp) return
    if (wavelength_um <= 0.2_dp .or. wavelength_um >= 2.5_dp) return
    if (number_density <= 0.0_dp) then
      valid = .True.
      return
    end if

    density_h2o = number_density*molecular_weight_h2o*amu
    sigma = (density_h2o*1.0e3_dp)/rho_s
    lambda_ratio = wavelength_um/lambda_s
    theta = temperature/temperature_s

    lorentz_term = sigma * &
      (coeff(1) + coeff(2)*sigma + coeff(3)*theta + &
       coeff(4)*lambda_ratio**2*theta + coeff(5)/lambda_ratio**2 + &
       coeff(6)/(lambda_ratio**2-lambda_uv**2) + &
       coeff(7)/(lambda_ratio**2-lambda_ir**2) + coeff(8)*sigma**2)

    valid = lorentz_term >= 0.0_dp .and. lorentz_term < 1.0_dp

  end subroutine water_lorentz_lorenz


  real(dp) function refractive_index_wavenumber(wavenumber_cm,A,B,C)
    implicit none

    real(dp), intent(in) :: wavenumber_cm, A, B, C

    refractive_index_wavenumber = &
      1.0_dp + (A + B/(C-wavenumber_cm**2))/1.0e8_dp

  end function refractive_index_wavenumber

end module gas_refractivity_mod
