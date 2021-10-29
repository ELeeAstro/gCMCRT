module lbl_tables_combine
  use optools_data_mod
  implicit none

contains

  subroutine combine_lbl_opacity(z,lbl_work,lbl_comb)
    implicit none

    integer, intent(in) :: z
    real(kind=dp), dimension(nlbl), intent(in) :: lbl_work
    real(kind=dp), intent(out) :: lbl_comb

    integer :: s

    ! Combine all the lbl species weighted by the layer VMR
    ! Convert to units of [cm-1] by * layer number density

    lbl_comb = 0.0_dp
    do s = 1, nlbl
      lbl_comb = lbl_comb + lbl_work(s) * N_lay(z) * VMR_lay(lbl_tab(s)%iVMR,z)
    end do

  end subroutine combine_lbl_opacity

end module lbl_tables_combine
