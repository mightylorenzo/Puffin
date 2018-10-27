! Copyright 2012-2018, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Module to calculate d/dz of the field values and electron macroparticle
!> coordinates.

module rhs

use paratype, only: wp, ip
use MPEquations
use wigglerVar
use FiElec1D
use FiElec

implicit none

contains

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Calculate d/dz of radiation field and electron macroparticle coordinates in
!> 6D phase space
!> @param[in] sz zbar
!> @param[in] sAr Real (x) component of A_perp
!> @param[in] sAi Imaginary (-y) component of A_perp
!> @param[in] sx scaled electron x coordinates
!> @param[in] sy scaled electron y coordinates
!> @param[in] sz2 scaled electron z2 coordinates
!> @param[in] spr scaled real (px) p_perp electron coordinates
!> @param[in] spi scaled real (-py) p_perp electron coordinates
!> @param[in] sgam scaled energy (gamma) electron coordinates
!> @param[out] sdx d/dz of scaled electron x coordinates
!> @param[out] sdy d/dz of scaled electron y coordinates
!> @param[out] sdz2 d/dz of scaled electron z2 coordinates
!> @param[out] sdpr d/dz of scaled real (px) p_perp electron coordinates
!> @param[out] sdpi d/dz of scaled real (-py) p_perp electron coordinates
!> @param[out] sdgam d/dz of scaled energy (gamma) electron coordinates
!> @param[out] sDADzr d/dz of real (x) component of A_perp
!> @param[out] sDADzi d/dz of real (-y) component of A_perp

  subroutine getrhs(sz, &
                    sAr, sAi, &
                    sx, sy, sz2, &
                    spr, spi, sgam, &
                    sdx, sdy, sdz2, &
                    sdpr, sdpi, sdgam, &
                    sDADzr, sDADzi, &
                    tVars, tBFields, tScale)

  use typeCalcParams, only: fcalcParams
  use typeScale, only: fScale
  use bfields, only: fbfields
  implicit none

!> Inputs %%%
!
! sZ - Propagation distance
! sA - current radiation field vals
! sy - current electron coordinates in all dimensions
!
! Output
! sb  - d/dz of electron phase space positions
! sDADz - RHS of field source term

  real(kind=wp), intent(in) :: sz
  real(kind=wp), contiguous, intent(in) :: sAr(:), sAi(:)
  real(kind=wp), contiguous, intent(in)  :: sx(:), sy(:), sz2(:), &
                                            spr(:), spi(:), sgam(:)

  real(kind=wp), contiguous, intent(inout)  :: sdx(:), sdy(:), sdz2(:), &
                                               sdpr(:), sdpi(:), sdgam(:)

  real(kind=wp), contiguous,  intent(inout) :: sDADzr(:), sDADzi(:)
  type(fcalcParams), intent(inout) :: tVars
  type(fScale), intent(in) :: tScale
  type(fbfields), intent(in) :: tBFields

  integer(kind=ipl) :: i, z2node
  integer :: error

!     Initialise right hand side to zero

  tVars%sField4ElecReal(:) = 0.0_WP
  tVars%sField4ElecImag(:) = 0.0_WP

!  call rhs_tmsavers(sz)  ! This can be moved later...

!$OMP PARALLEL

  call tScale%getP2(tVars%sp2, sgam, spr, spi)

  if (tVars%q1D) then

    tVars%p_nodes = int(sz2 / tVars%dz2, kind=ip) + 1_IP - (tVars%fz2-1)

  else

!$OMP WORKSHARE

    tVars%p_nodes = (int( (sx+tVars%halfx)  / tVars%dx, kind=ip)  + 1_IP) + &
              (int( (sy+tVars%halfy)  / tVars%dy, kind=ip) * tVars%inNX )  + &   !  y 'slices' before primary node
              (tVars%inNX * tVars%inNY * &
                              int(sz2  / tVars%dz2, kind=ip) ) - &
                              (tVars%fz2-1)*tVars%inNN  ! transverse slices before primary node

!$OMP END WORKSHARE

  end if

  if (tVars%q1D) then

!    call getInterps_1D(sz2, tVars)
!    if (tVars%qPArrOK) then
!      call getFFelecs_1D(sAr, sAi, tVars)
!      call getSource_1D(sDADzr, sDADzi,  spr, spi, sgam, tVars)
!    end if

  else

    call getInterps_3D(sx, sy, sz2, tVars)
    if ((tVars%qPArrOK) .and. (tVars%qInnerXYOK)) then
      call getFFelecs_3D(sAr, sAi, tVars)
      call getSource_3D(sDADzr, sDADzi, spr, spi, sgam, s_chi_bar_G, tVars)
    end if

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Calculate electron d/dz of electron equations - if needed

  if (.not. tVars%qElectronFieldCoupling) then
    tVars%sField4ElecReal = 0.0_WP
    tVars%sField4ElecImag = 0.0_WP
  end if


    if (tVars%qElectronsEvolve) then

        call tBFields%getBFields(sx, sy, sz, tScale, &
                        tVars%bxu, tVars%byu, tVars%bzu)

        call dz2dz_f(sx, sy, sz2, spr, spi, sgam, tVars, &
                     sdz2)

        call dxdz_f(sx, sy, sz2, spr, spi, sgam, tVars,  &
                    sdx)

        call dydz_f(sx, sy, sz2, spr, spi, sgam, tVars, &
                    sdy)

        call dppdz_r_f(sx, sy, sz2, spr, spi, sgam, tVars, &
                       sdpr)

        call dppdz_i_f(sx, sy, sz2, spr, spi, sgam, tVars, &
                       sdpi)

        call dgamdz_f(sx, sy, sz2, spr, spi, sgam, tVars, &
                     sdgam)

    end if

!$OMP END PARALLEL

!     Switch field off

    if (.not. tVars%qFieldEvolve) then
       sDADzr = 0.0_WP
       sDADzi = 0.0_WP
    end if

!     if electrons not allowed to evolve then

    if (.not. tVars%qElectronsEvolve) then
       sdpr = 0.0_wp
       sdpi = 0.0_wp
       sdgam = 0.0_wp
       sdx   = 0.0_wp
       sdy   = 0.0_wp
       sdz2 = 0.0_wp
    end if

  end subroutine getrhs

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Initialize data used in the calculation of d/dz of electron beam + radiation
!> field quantities
!> @param[in] sz zbar

! subroutine rhs_tmsavers(sz)
! 
! use rhs_vars
! 
! real(kind=wp), intent(in) :: sz
! 
!   ioutside=0
! 
! 
! !     Define the size of each element
! 
!   dx = sLengthOfElmX_G
!   dy = sLengthOfElmY_G
!   dz2 = sLengthOfElmZ2_G
! 
!   dV3 = sLengthOfElmX_G*sLengthOfElmY_G*sLengthOfElmZ2_G
! 
! 
! !     Time savers
! 
!   sInv2rho    = 1.0_WP/(2.0_WP * sRho_G)
! 
!   ZOver2rho   = sz * sInv2rho
!   salphaSq    = (2.0_WP * sGammaR_G * sRho_G / sAw_G)**2
! 
!   un = sqrt(fx_G**2.0_WP + fy_G**2.0_WP)
! 
! 
! !     number of transverse nodes
! 
!   ntrans = NX_G * NY_G
! 
! !     Diff between real and imaginary nodes in the reduced system
! 
!   retim = nspinDX*nspinDY*nZ2_G
! 
!   econst = sAw_G/(sRho_G*sqrt(2.0_WP*(fx_G**2.0_WP+fy_G**2.0_WP)))
! 
!   nc = 2.0_WP*saw_G**2/(fx_G**2.0_WP + fy_G**2.0_WP)
! 
!   nb = 2.0_WP * sRho_G / ((fx_G**2.0_WP+fy_G**2.0_WP)*sEta_G)
! 
!   maxEl = maxval(procelectrons_G)
!   qoutside=.FALSE.
!   iOutside=0_IP
! 
! !  halfx = ((ReducedNX_G-1) / 2.0_WP) * sLengthOfElmX_G
! !  halfy = ((ReducedNY_G-1) / 2.0_WP) * sLengthOfElmY_G
! 
!   halfx = ((nspinDX-1) / 2.0_WP) * sLengthOfElmX_G
!   halfy = ((nspinDY-1) / 2.0_WP) * sLengthOfElmY_G
! 
! 
! end subroutine rhs_tmsavers

end module rhs
