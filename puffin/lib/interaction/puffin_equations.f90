! Copyright 2012-2018, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause

module Equations


use puffin_kinds, only: WP, IP
use Globals, only: iUndStart_G, iUndEnd_G, iUndMain_G
use rhs_vars, only: sInv2rho, sp2, sField4ElecReal, sField4ElecImag, bxu, byu, bzu
use GlobalTypes, only: tUndulator

implicit none (type, external)
private

public :: adjundplace, alct_e_srtcts, bxu, byu, bzu, dalct_e_srtcts, &
           sField4ElecImag, sField4ElecReal, sInv2rho, sp2


contains





  subroutine alct_e_srtcts(ar_sz)

    implicit none (type, external)

! Allocate the arrays used in the calculation of
! the electron eqns

    integer(kind=ip), intent(in) :: ar_sz

    allocate(sp2(ar_sz), sField4ElecReal(ar_sz), &
             sField4ElecImag(ar_sz))! , Lj(ar_sz))

    allocate(bxu(ar_sz), byu(ar_sz), bzu(ar_sz))

  end subroutine alct_e_srtcts



  subroutine dalct_e_srtcts()

    implicit none (type, external)

! Allocate the arrays used in the calculation of
! the electron eqns

    deallocate(sp2, sField4ElecReal, &
             sField4ElecImag)! , Lj(ar_sz))

    deallocate(bxu, byu, bzu)

  end subroutine dalct_e_srtcts



  subroutine adjUndPlace(szl, und)

! Sets und%undulator_position based on current z position szl.
! Uses und%model_undulator_ends, und%z_start_undulator, und%z_end_undulator.
! Replaces the global iUndPlace_G; no longer touches any globals.

    real(kind=wp), intent(in) :: szl
    type(tUndulator), intent(inout) :: und

      if (und%model_undulator_ends) then

        if (szl < 0) then

          print*, "undulator section not recognised, sz < 0!!"
          stop

        else if (sZl <= und%z_start_undulator) then

          und%undulator_position = iUndStart_G

        else if (sZl >= und%z_end_undulator) then

          und%undulator_position = iUndEnd_G

        else if ((sZl > und%z_start_undulator) .and. (sZl < und%z_end_undulator)) then

          und%undulator_position = iUndMain_G

        else

          print*, "undulator section not recognised, sz > z_end_undulator!!"
          stop

        end if

      else

        und%undulator_position = iUndMain_G

      end if

  end subroutine adjUndPlace




end module equations
