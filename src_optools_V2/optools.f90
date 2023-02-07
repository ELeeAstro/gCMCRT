program optools
  use optools_data_mod
  use read_io, only : read_optools_par, read_prf, read_wavelengths, read_clprf
  use CK_tables_mod, only : calc_CK_table
  use lbl_tables_mod, only : calc_lbl_table
  use Ray_tables_mod, only : calc_Ray_table
  use CIA_tables_mod, only : calc_CIA_table
  use cloud_tables_mod, only: calc_cloud_table
  use xsec_tables_mod, only: calc_xsec_table
  implicit none

  !! Idea is run this opacity toolkit separatly, and output opacity variables for each cell
  !! Then run MC code to read in the opacity files and perform the radiative transfer calcultions

  ! Read basic parameter file
  call read_optools_par()

  ! Read the 1D or 3D flattened 3D T-p and gas Volume Mixing Ratio (VMR) file
  call read_prf()

  ! If cloud opacity wanted, read cloud prf
  if (cloud_opc .eqv. .True.) then
      call read_clprf()
  end if

  ! Read wavelength array file wavelength.wl at bin centers
  call read_wavelengths()

  ! Read input k-tables and perform interpolation (+ random overlap)
  ! Output CMCRT k formatted table file
  if (corr_k .eqv. .True.) then
    call calc_CK_table()
  end if

  ! Read input lbl tables and perform interpolation
  ! Output CMCRT lbl formatted table file
  if (lbl .eqv. .True.) then
    call calc_lbl_table()
  end if

  ! Read continuum data tables and perform interpolation
  ! Output CMCRT contiuum formatted file
  if (conti .eqv. .True.) then
    call calc_CIA_table()
  end if

  ! Calculate gas phase Rayleigh scattering parameters
  ! Output CMCRT Rayleigh formatted file
  if (Ray_scat .eqv. .True.) then
    call calc_Ray_table()
  end if

  ! Read cloud particle inputs and parameters and perform Mie calculation
  ! Output CMCRT Cloud Properties formatted file
  if (cloud_opc .eqv. .True.) then
    call calc_cloud_table()
  end if

  if (xsec_opc .eqv. .True.) then
    call calc_xsec_table()
  end if

  ! Read gas phase optical constants and find refraction parameters
  ! Output CMCRT refrac file
  !if (refrac .eqv. .True.) then
    !call calc_refrac()
  !end if

  ! Read photochemistry xsec inputs and parameters
  ! Output CMCRT phchem file
  !if (ph_chem .eqv. .True.) then
    !call read_phchem()
    !call calc_phchem()
  !end if

end program optools
