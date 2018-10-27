! Copyright 2012-2018, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause

!> @author
!> Lawrence Campbell,
!> University of Strathclyde, 
!> Glasgow, UK
!> @brief
!> Module which calculates d/dz of the local electron and field variables, and 
!> then sums the global regions together.

module Derivative

use paratype, only: wp, ip
implicit none

contains

!> Subroutine which calculate d/dz of the local electron and field 
!> variables, and then sums the global regions together.
!> @param sz position in undulator module.
!> @param sAr real field 
!> @param sAi imaginary field 
!> @param sx electron macroparticles' x position

  subroutine derivs(sz, sAr, sAi, sx, sy, sz2, spr, spi, sp2, &
                    sdx, sdy, sdz2, sdpr, sdpi, sdp2, sdAr, sdAi, &
                    tScale, tMPI)
    use typeCalcParams, only: fcalcParams
    use typeScale, only: fScale
    use rhs, only: getrhs
    use parafield, only: upd8da, fz2, bz2
    use globals !, only: fieldMesh, nz2_g, fx_G, fy_G, kx_und_G, ky_und_G, n2col, &
                !       sZFS, sZFE, iOutInfo_G, sLengthOfElmX_G, sLengthOfElmY_G, &
                !       sLengthOfElmZ2_G, sKBetaXSF_G, sKBetaYSF_G, &
                !       qFocussing_G, qParrOK_G, qInnerXYOK_G
    use ParallelInfoType, only: cParallelInfoType
    use TransformInfoType
    implicit none

    real(kind=wp), intent(in)  :: sz
    real(kind=wp), contiguous, intent(in)  :: sAr(:), sAi(:)
    real(kind=wp), contiguous, intent(in)  :: sx(:), sy(:), sz2(:), &
                                  spr(:), spi(:), sp2(:)

    real(kind=wp), contiguous, intent(inout)  :: sdx(:), sdy(:), sdz2(:), &
                                  sdpr(:), sdpi(:), sdp2(:)

    real(kind=wp), contiguous, intent(inout) :: sdAr(:), sdAi(:)
    type(fScale), intent(in) :: tScale
    type(cParallelInfoType), intent(in) :: tMPI

    type(fcalcParams) :: tVars
    integer(kind=ip) :: error, iArEr, maxEl
    logical :: qOKL

    maxEl = maxval(procelectrons_G)

    call tVars%init(sLengthOfElmX_G, sLengthOfElmY_G, sLengthOfElmZ2_G, &
                    nspinDX, nspinDY, nz2_g, fz2, bz2, fieldMesh, &
                    fx_G, fy_G, kx_und_G, ky_und_G, n2col, sZFS, sZFE, &
                    sKBetaXSF_G, sKBetaYSF_G, zUndType_G, procelectrons_G(1), maxEl, &
                    tTransInfo_G%qOneD, qFocussing_G, qParrOK_G, qInnerXYOK_G, &
                    qElectronFieldCoupling_G, qElectronsEvolve_G, qFieldEvolve_G, &
                    iOutInfo_G, tScale)

!     Get RHS of field eqn and d/dz of electron variables

    call getrhs(sz, &
                sAr, sAi, &
                sx, sy, sz2, &
                spr, spi, sp2, &
                sdx, sdy, sdz2, &
                sdpr, sdpi, sdp2, &
                sdAr, sdAi, &
                tVars, tScale)

!    update fields in buffers

    call upd8da(sdAr, sdAi)

!     Check to see if parallel field setup OK...

    if (.not. tVars%qPArrOK) then
      iArEr = 1_ip
    else
      iArEr = 0_ip
    end if

    call mpi_allreduce(mpi_in_place, iArEr, 1, &
          MPI_INTEGER, &
          MPI_SUM, MPI_COMM_WORLD, error)

    if (iArEr > 0_ip) then
      tVars%qPArrOK = .false.
      if ((tMPI%qRoot) .and. (tVars%ioutInfo > 2)) then
        print*, 'electron outside parallel bounds!'
        print*, 'Emergency redistribute!!!'
        print*, 'If this happens often, then &
              & it is possible the parallel &
              & tuning parameters are inefficient, &
              & and not suitable...'
      end if
    end if

    if (.not. tVars%qInnerXYOK) then
      iArEr = 1_ip
    else
      iArEr = 0_ip
    end if

    call mpi_allreduce(mpi_in_place, iArEr, 1, &
          MPI_INTEGER, &
          MPI_SUM, MPI_COMM_WORLD, error)

    if (iArEr > 0_ip) then
      tVars%qInnerXYOK = .false.
      if ((tMPI%qRoot) .and. (tVars%ioutInfo > 2) ) then
        print*, 'electron outside transverse bounds!'
        print*, 'Emergency redistribute!!!'
      end if
    end if

    qPArrOK_G = tVars%qParrOK
    qInnerXYOK_G = tVars%qInnerXYOK
    call tVars%destroy()

  end subroutine derivs

end module Derivative
