! ################################################
! Copyright 2012-2018, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause
! ################################################
!
!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Module containing the routines which calculate the scaled magnetic b-fields
!> for the Lorentz force in Puffin


module bfields

use paratype, only: wp, ip
implicit none
private

integer(kind=ip), private :: iUndPlace_G
integer(kind=ip), parameter, private :: iUndStart_G = 1_ip, &
                                        iUndEnd_G = 2_ip, &
                                        iUndMain_G = 0_ip

public :: getBFields
private :: adjUndPlace, getBXfield, getBYfield, getBZfield

contains

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Calculate the b-fields. Calls bx, by and bz subroutines
!> @param[in] sx Macroparticle coordinates in xbar
!> @param[in] sy Macroparticle coordinates in ybar
!> @param[in] sz current scaled distance though undulator, zbar
!> @param[out] bxj scaled b-field in x direction for macroparticles
!> @param[out] byj scaled b-field in y direction for macroparticles
!> @param[out] bzj scaled b-field in z direction for macroparticles

  subroutine getBFields(sx, sy, sZ, tVars, &
                        bxj, byj, bzj)

  use typecalcParams, only: fcalcParams
  implicit none

!   subroutine to calculate the scaled magnetic fields
!   at a given zbar

  real(kind=wp), contiguous, intent(in) :: sx(:), sy(:)
  real(kind=wp), intent(in) :: sz
  type(fcalcParams), intent(in) :: tVars
  real(kind=wp), contiguous, intent(out) :: bxj(:), byj(:), bzj(:)
  real(kind=wp) :: alph

  call adjUndPlace(sz, tVars)
!  call getAlpha(sz, tVars, alph)

  call getBXfield(sx, sy, sz, tVars, bxj)
  call getBYfield(sx, sy, sz, tVars, byj)
  call getBZfield(sx, sy, sz, tVars, bzj)
  
!  bxj = alph * bxj
!  byj = alph * byj
!  bzj = alph * bzj

  end subroutine getBFields

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Determine whether we are in the ends or the active section of undulator.

subroutine adjUndPlace(szl, tVars)
  use typecalcParams, only: fcalcParams
  use globals, only: qUndEnds_G
  implicit none
  real(kind=wp) :: szl
  type(fcalcParams), intent(in) :: tVars

  if (qUndEnds_G) then
    if (szl < 0) then
      print*, 'undulator section not recognised, sz < 0!!'
      stop
    else if (sZl <= tVars%sZFS) then
      iUndPlace_G = iUndStart_G
    else if (sZl >= tVars%sZFE) then
      iUndPlace_G = iUndEnd_G
    else if ((sZl > tVars%sZFS) .and. (sZl < tVars%sZFE)) then
      iUndPlace_G = iUndMain_G
    else
      print*, 'undulator section not recognised, sz > sZFE!!'
      stop
    end if
  else
    iUndPlace_G = iUndMain_G
  end if

end subroutine adjUndPlace

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Calculates tuning of wiggler \f$ \alpha \f$

! subroutine getAlpha(tVars, sZ, n2col)
!   use typecalcParams, only: fcalcParams
!   implicit none
!   real(kind=wp), intent(in) :: sZ
!   type(fcalcParams), intent(in) :: tVars
!   real(kind=wp), intent(out) :: n2col
! 
!   if ((sZ >= tVars%sZFS) .and. (sZ <= tVars%sZFE)) then
!     n2col = tVars%n2col0  + tVars%undgrad*(sz - tVars%sZFS)  ! linear taper
!   else if (sZ < tVars%sZFS) then
!     n2col = tVars%n2col0
!   else if (sZ > tVars%sZFE) then
!     tVars%n2col = tVars%n2col0  + tVars%undgrad*(tVars%sZFE - tVars%sZFS)
!   end if
! 
! end subroutine getAlpha


subroutine getBXfield(sx, sy, sz, tVars, bxj)

  use typecalcParams, only: fcalcParams
  use typesAndConstants, only: pi
  implicit none

  real(kind=wp), contiguous, intent(in) :: sx(:), sy(:)
  real(kind=wp), intent(in) :: sz
  type(fcalcParams), intent(in) :: tVars
  real(kind=wp), contiguous, intent(out) :: bxj(:)

!    Local vars:-

  real(kind=wp) :: szt

!      cc1 = sqrt(tVars%eta) / (2_wp * sRho_G * tVars%kyu)




  szt = sZ
  szt = szt / 2_wp / tVars%rho


!  ####################################################
!    Curved pole case - planar wiggler with focusing
!    in both x and y (electron wiggles in x)

  if (tVars%zUndType == 'curved') then

    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
      bxj = tVars%kxu / tVars%kyu * sinh(tVars%kxu * sx) &
            * sinh(tVars%kyu * sy) &
            * szt / 4_wp / pi * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - tVars%sZFE
      szt = szt / 2_wp / tVars%rho

!$OMP WORKSHARE
      bxj = tVars%kxu / tVars%kyu * sinh(tVars%kxu * sx) &
            * sinh(tVars%kyu * sy) &
            * (-szt / 4_wp / pi + 1_wp) * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      bxj = tVars%kxu / tVars%kyu * sinh(tVars%kxu * sx) &
            * sinh(tVars%kyu * sy) &
            * sin(szt)
!$OMP END WORKSHARE

    end if

!    END curved pole field description
!  ####################################################








!  ####################################################
!    Plane-pole case - planar wiggler with focusing
!    only in y (and electron will wiggle in x)

  else if (tVars%zUndType == 'planepole')  then

    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
      bxj = 0_wp
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - tVars%sZFE
      szt = szt / 2_wp / tVars%rho

!$OMP WORKSHARE
      bxj = 0_wp
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      bxj = 0_wp
!$OMP END WORKSHARE

    end if

!    END plane pole undulator field description
!  ####################################################







!  ####################################################
!    Helical case - helical wiggler with focusing
!    in x and y (and electron will wiggle in x and y)

  else if (tVars%zUndType == 'helical')  then

    if (iUndPlace_G == iUndStart_G) then


!      bxj = sin(szt / 8_wp) * &
!               cos(szt / 8_wp) * sin(szt) / 4_wp   &
!             +  sin(szt/8_wp)**2_wp  * cos(szt)
!$OMP WORKSHARE
      bxj = (sZ - pi * tVars%rho) / (6_wp * pi*tVars%rho) * cos(szt)
!$OMP END WORKSHARE

      if (sZ < pi * tVars%rho) then 
!$OMP WORKSHARE
        bxj = 0.0_wp
!$OMP END WORKSHARE
else if (sZ > 7_wp * pi * tVars%rho) then
!$OMP WORKSHARE
        bxj = cos(szt)
!$OMP END WORKSHARE
      end if

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - tVars%sZFE
!      szt = szt / 2_wp / tVars%rho


!      bxj = - cos(szt / 8_wp) * &
!              sin(szt / 8_wp) * sin(szt)  / 4_wp  &
!            +  cos(szt/8_wp)**2_wp  * cos(szt)

!$OMP WORKSHARE
      bxj =  - (szt - 7.0_wp * pi * tVars%rho) / (6_wp * pi * tVars%rho) * &
               cos(szt / 2_wp / tVars%rho)
!$OMP END WORKSHARE

      if (sZt < pi * tVars%rho) then 
!$OMP WORKSHARE
        bxj = cos(szt / 2_wp / tVars%rho)
!$OMP END WORKSHARE
else if (sZt > 7_wp * pi * tVars%rho) then
!$OMP WORKSHARE
        bxj = 0.0_wp
!$OMP END WORKSHARE
      end if

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      bxj = cos(szt)
!$OMP END WORKSHARE

    end if

!    END helical undulator field description
!  ####################################################






  else






!  ####################################################
!    'puffin' elliptical undulator...
!    with variable x and y polarization...

    if (iUndPlace_G == iUndStart_G) then


!      bxj = tVars%ux * sin(szt / 8_wp) * &
!               cos(szt / 8_wp) * sin(szt) / 4_wp   &
!             +  sin(szt/8_wp)**2_wp  * cos(szt)
!$OMP WORKSHARE
      bxj = tVars%ux * (sZ - pi * tVars%rho) / (6_wp * pi*tVars%rho) * cos(szt)
!$OMP END WORKSHARE

      if (sZ < pi * tVars%rho) then 
!$OMP WORKSHARE
        bxj = 0.0_wp
!$OMP END WORKSHARE
else if (sZ > 7_wp * pi * tVars%rho) then
!$OMP WORKSHARE
        bxj = tVars%ux * cos(szt)
!$OMP END WORKSHARE
      end if

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - tVars%sZFE
      !szt = szt / 2_wp / tVars%rho

!$OMP WORKSHARE
      bxj =  -tVars%ux * (szt - 7.0_wp * pi * tVars%rho) / (6_wp * pi * tVars%rho) * &
               cos(szt / 2_wp / tVars%rho)
!$OMP END WORKSHARE

      if (sZt < pi * tVars%rho) then 
!$OMP WORKSHARE
        bxj = tVars%ux * cos(szt / 2_wp / tVars%rho)
!$OMP END WORKSHARE
else if (sZt > 7_wp * pi * tVars%rho) then
!$OMP WORKSHARE
        bxj = 0.0_wp
!$OMP END WORKSHARE
      end if
      

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      bxj = tVars%ux*cos(szt)
!$OMP END WORKSHARE

    end if

!    END elliptical undulator description
!  ####################################################




  end if


!   Focusing component (non-physical)

    if (tVars%qSF) then

!$OMP WORKSHARE
        bxj = sqrt(tVars%eta) * tVars%kbySF**2.0_wp / tVars%kappa &
              * sy + bxj
!$OMP END WORKSHARE

    end if


  end subroutine getBXfield




!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Calculate the by-fields for given coordinates
!> @param[in] sx Macroparticle coordinates in xbar
!> @param[in] sy Macroparticle coordinates in ybar
!> @param[in] sz current scaled distance though undulator, zbar
!> @param[out] byj scaled b-field in y direction for macroparticles

  subroutine getBYfield(sx, sy, sz, tVars, byj)

  use typecalcParams, only: fcalcParams
  use typesAndConstants, only: pi
  implicit none

  real(kind=wp), contiguous, intent(in) :: sx(:), sy(:)
  real(kind=wp), intent(in) :: sz
  type(fcalcParams), intent(in) :: tVars
  real(kind=wp), contiguous, intent(out) :: byj(:)

!    Local vars:-

  real(kind=wp) :: szt

  szt = sZ
  szt = szt / 2_wp / tVars%rho


!  ####################################################
!    Curved pole case - planar wiggler with focusing
!    in both x and y (electron wiggles in x)


  if (tVars%zUndType == 'curved') then

    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
      byj = cosh(tVars%kxu * sx) &
            * cosh(tVars%kyu * sy) &
            *  szt / 4_wp / pi * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - tVars%sZFE
      szt = szt / 2_wp / tVars%rho

!$OMP WORKSHARE
      byj = cosh(tVars%kxu * sx) &
            * cosh(tVars%kyu * sy) &
            * (-szt / 4_wp / pi + 1_wp) * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      byj = cosh(tVars%kxu * sx) &
            * cosh(tVars%kyu * sy) &
            * sin(szt)
!$OMP END WORKSHARE

    end if

!    END curved pole field description
!  ####################################################







!  ####################################################
!    Plane-pole case - planar wiggler with focusing
!    only in y (and electron will wiggle in x)



  else if (tVars%zUndType == 'planepole')  then

    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
!      byj = cosh( sqrt(tVars%eta) / 2_wp / tVars%rho * sy) * &
!            (  (- sin(szt / 8_wp) * &
!               cos(szt / 8_wp) * cos(szt) / 4_wp   &
!             +  sin(szt/8_wp)**2_wp  * sin(szt) )  )


      !byj = sin(szt/8_wp)**2_wp * sin(szt)
      byj = szt / 4_wp / pi * sin(szt)
!$OMP END WORKSHARE

!    print*, "hehehe"

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - tVars%sZFE
      szt = szt / 2_wp / tVars%rho

!$OMP WORKSHARE
!      byj = cosh( sqrt(tVars%eta) / 2_wp / tVars%rho * sy) * &
!            (  cos(szt / 8_wp) * &
!              sin(szt / 8_wp) * cos(szt)  / 4_wp  &
!            +  cos(szt/8_wp)**2_wp  * sin(szt)  )
            
      byj = (-szt / 4_wp / pi + 1_wp) * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
!      byj = cosh( sqrt(tVars%eta) / 2_wp / tVars%rho * sy) &
!            * sin(szt)
      byj = sin(szt)
!$OMP END WORKSHARE

    end if

!    END plane pole undulator field description
!  ####################################################









!  ####################################################
!    Helical case - helical wiggler with focusing
!    in x and y (and electron will wiggle in x and y)

  else if (tVars%zUndType == 'helical')  then

    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
      byj = szt / 4_wp / pi * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - tVars%sZFE
      szt = szt / 2_wp / tVars%rho

!$OMP WORKSHARE
      byj = (-szt / 4_wp / pi + 1_wp) * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      byj = sin(szt)
!$OMP END WORKSHARE

    end if

!    END helical undulator field description
!  ####################################################








  else






!  ####################################################
!    'puffin' elliptical undulator...
!    with variable x and y polarization...


    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
      byj = tVars%uy * szt / 4_wp / pi * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - tVars%sZFE
      szt = szt / 2_wp / tVars%rho

!$OMP WORKSHARE
      byj = tVars%uy * (-szt / 4_wp / pi + 1_wp) * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      byj = tVars%uy * sin(szt)
!$OMP END WORKSHARE

    end if

!    END elliptical undulator description
!  ####################################################




  end if

!   Focusing component (non-physical)

    if (tVars%qSF) then

!$OMP WORKSHARE
      byj = -sqrt(tVars%eta) * tVars%kbxSF**2.0_wp / tVars%kappa &
            * sx + byj
!$OMP END WORKSHARE

    end if

  end subroutine getBYfield


!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Calculate the bz-fields for given coordinates
!> @param[in] sx Macroparticle coordinates in xbar
!> @param[in] sy Macroparticle coordinates in ybar
!> @param[in] sz current scaled distance though undulator, zbar
!> @param[out] bzj scaled b-field in z direction for macroparticles

subroutine getBZfield(sx, sy, sz, tVars, bzj)
  use typecalcParams, only: fcalcParams
  use typesAndConstants, only: pi
  implicit none
  real(kind=wp), contiguous, intent(in) :: sx(:), sy(:)
  real(kind=wp), intent(in) :: sz
  type(fcalcParams), intent(in) :: tVars
  real(kind=wp), contiguous, intent(out) :: bzj(:)

!    Local vars:-

  real(kind=wp) :: szt

  szt = sZ
  szt = szt / 2.0_wp / tVars%rho

  if (tVars%q1D) then

!  ####################################################
!    1D case - no z-component of magnetic field

!$OMP WORKSHARE
    bzj = 0.0_wp
!$OMP END WORKSHARE

  else

!  ####################################################
!    Curved pole case - planar wiggler with focusing
!    in both x and y (electron wiggles in x)


  if (tVars%zUndType == 'curved') then

    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
      bzj = sqrt(tVars%eta) / 2 / tVars%rho / tVars%kyu &
                * cosh(tVars%kxu * sx) &
            * sinh(tVars%kyu * sy) &
            * (szt / 4_wp / pi) * cos(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - tVars%sZFE
      szt = szt / 2_wp / tVars%rho

!$OMP WORKSHARE
      bzj = sqrt(tVars%eta) / 2_wp / tVars%rho / tVars%kxu &
            * cosh(tVars%kxu * sx) * sinh(tVars%kyu * sy) &
            * (-szt / 4_wp / pi + 1_wp) * cos(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then


!$OMP WORKSHARE
      bzj = sqrt(tVars%eta) / 2_wp / tVars%rho / tVars%kxu * &
           cosh(tVars%kxu * sx) * sinh(tVars%kyu * sy) &
            * cos(szt)
!$OMP END WORKSHARE

    end if

!    END curved pole field description
!  ####################################################






!  ####################################################
!    Plane-pole case - planar wiggler with focusing
!    only in y (and electron will wiggle in x)

  else if (tVars%zUndType == 'planepole')  then

    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
      bzj = sinh( sqrt(tVars%eta) / 2_wp / tVars%rho * sy) * &
            (szt / 4_wp / pi) * cos(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - tVars%sZFE
      szt = szt / 2_wp / tVars%rho

!$OMP WORKSHARE
      bzj = sinh( sqrt(tVars%eta) / 2_wp / tVars%rho * sy) * &
            (-szt / 4_wp / pi + 1_wp) * cos(szt) 
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      bzj = sinh( sqrt(tVars%eta) / 2_wp / tVars%rho * sy) &
            * cos(szt)
!$OMP END WORKSHARE

    end if

!    END plane pole undulator field description
!  ####################################################






!  ####################################################
!    Helical case - helical wiggler with focusing
!    in x and y (and electron will wiggle in x and y)

  else if (tVars%zUndType == 'helical')  then

    if (iUndPlace_G == iUndStart_G) then

! ...from x-comp:

!$OMP WORKSHARE
      bzj = - sqrt(tVars%eta) / 2 / tVars%rho * (sZ - pi * tVars%rho) / (6_wp * pi*tVars%rho) &
            * sx * sin(szt)            
!$OMP END WORKSHARE

      if (sZ < pi * tVars%rho) then 
!$OMP WORKSHARE
        bzj = 0.0_wp
!$OMP END WORKSHARE
else if (sZ > 7_wp * pi * tVars%rho) then
!$OMP WORKSHARE
        bzj = - sqrt(tVars%eta) / 2 / tVars%rho * sx * sin(szt)
!$OMP END WORKSHARE
      end if
      

! ...and from y-comp:

!$OMP WORKSHARE
      bzj = bzj + szt / 4_wp / pi * sqrt(tVars%eta) / 2 / tVars%rho * sy * cos(szt)
!$OMP END WORKSHARE

!  !$OMP WORKSHARE
!        bzj = sqrt(tVars%eta) / 2 / tVars%rho * (     &
!              sx *  sin(szt) )    + &
!              sy * ( -1/32_wp * cos(szt/4_wp) * cos(szt) + &
!                      1/4_wp * sin(szt/4_wp) * sin(szt) + &
!                      sin(szt/8_wp)**2 * cos(szt) ) )
!  !$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - tVars%sZFE
      !szt = szt / 2_wp / tVars%rho

! ...from x-comp:

!$OMP WORKSHARE
      bzj = - sqrt(tVars%eta) / 2 / tVars%rho * (szt - 7.0_wp * pi * tVars%rho) / &
              (6_wp * pi * tVars%rho) * sx * sin(szt / 2_wp / tVars%rho)            
!$OMP END WORKSHARE

      if (szt < pi * tVars%rho) then 
!$OMP WORKSHARE
        bzj = - sqrt(tVars%eta) / 2 / tVars%rho * sx * sin(szt / 2_wp / tVars%rho)
!$OMP END WORKSHARE
else if (szt > 7_wp * pi * tVars%rho) then
!$OMP WORKSHARE
        bzj = 0.0_wp
!$OMP END WORKSHARE
      end if
      

! ...and from y-comp:

!$OMP WORKSHARE
      bzj = bzj + (-szt / 8_wp / pi / tVars%rho + 1_wp) * &
            sqrt(tVars%eta) / 2 / tVars%rho * sy * cos(szt / 2_wp / tVars%rho)
!$OMP END WORKSHARE












! !$OMP WORKSHARE
!       bzj = sqrt(tVars%eta) / 2 / tVars%rho * (     &
!             sx * ( -1/32_wp * cos(szt/4_wp) * sin(szt) - &
!                     1/4_wp * sin(szt/4_wp) * cos(szt) - &
!                     cos(szt/8_wp)**2 * sin(szt) )    + &
!             sy * ( 1/32_wp * cos(szt/4_wp) * cos(szt) - &
!                     1/4_wp * sin(szt/4_wp) * sin(szt) + &
!                     cos(szt/8_wp)**2 * cos(szt) ) )
! !$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      bzj = sqrt(tVars%eta) / 2 / tVars%rho * &
            ( -sx * sin(szt)  + sy * cos(szt) )
!$OMP END WORKSHARE

    end if

!    END helical undulator field description
!  ####################################################





  else




!  ####################################################
!    'puffin' elliptical undulator...
!    with variable x and y polarization...

    if (iUndPlace_G == iUndStart_G) then

! ...from x-comp:

!$OMP WORKSHARE
      bzj = - sqrt(tVars%eta) / 2 / tVars%rho * (sZ - pi * tVars%rho) / (6_wp * pi*tVars%rho) &
            * tVars%ux * sx * sin(szt)            
!$OMP END WORKSHARE

      if (sZ < pi * tVars%rho) then 
!$OMP WORKSHARE
        bzj = 0.0_wp
!$OMP END WORKSHARE
else if (sZ > 7_wp * pi * tVars%rho) then
!$OMP WORKSHARE
        bzj = - sqrt(tVars%eta) / 2 / tVars%rho * tVars%ux * sx * sin(szt)
!$OMP END WORKSHARE
      end if
      

! ...and from y-comp:

!$OMP WORKSHARE
      bzj = bzj + szt / 4_wp / pi * sqrt(tVars%eta) / 2 / tVars%rho * tVars%uy * sy * cos(szt)
!$OMP END WORKSHARE





!!$OMP WORKSHARE
!      bzj = sqrt(tVars%eta) / 2 / tVars%rho * (     &
!        tVars%ux*sx * ( 1/32_wp * cos(szt/4_wp) * sin(szt) + &
!                    1/4_wp * sin(szt/4_wp) * cos(szt) - &
!                    sin(szt/8_wp)**2 * sin(szt) )    + &
!        tVars%uy*sy * ( -1/32_wp * cos(szt/4_wp) * cos(szt) + &
!                    1/4_wp * sin(szt/4_wp) * sin(szt) + &
!                    sin(szt/8_wp)**2 * cos(szt) ) )
!!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - tVars%sZFE
!      szt = szt / 2_wp / tVars%rho








! ...from x-comp:

!$OMP WORKSHARE
      bzj = - sqrt(tVars%eta) / 2 / tVars%rho * (szt - 7.0_wp * pi * tVars%rho) / &
              (6_wp * pi * tVars%rho) * tVars%ux * sx * sin(szt / 2_wp / tVars%rho)            
!$OMP END WORKSHARE

      if (szt < pi * tVars%rho) then 
!$OMP WORKSHARE
        bzj = - sqrt(tVars%eta) / 2 / tVars%rho * tVars%ux * sx * sin(szt / 2_wp / tVars%rho)
!$OMP END WORKSHARE
else if (szt > 7_wp * pi * tVars%rho) then
!$OMP WORKSHARE
        bzj = 0.0_wp
!$OMP END WORKSHARE
      end if
      

! ...and from y-comp:

!$OMP WORKSHARE
      bzj = bzj + (-szt / 8_wp / pi / tVars%rho + 1_wp) * &
            sqrt(tVars%eta) / 2 / tVars%rho * tVars%uy * sy * cos(szt / 2_wp / tVars%rho)
!$OMP END WORKSHARE













!!$OMP WORKSHARE
!      bzj = sqrt(tVars%eta) / 2 / tVars%rho * (     &
!        tVars%ux*sx * ( -1/32_wp * cos(szt/4_wp) * sin(szt) - &
!                    1/4_wp * sin(szt/4_wp) * cos(szt) - &
!                    cos(szt/8_wp)**2 * sin(szt) )    + &
!        tVars%uy*sy * ( 1/32_wp * cos(szt/4_wp) * cos(szt) - &
!                    1/4_wp * sin(szt/4_wp) * sin(szt) + &
!                    cos(szt/8_wp)**2 * cos(szt) ) )
!!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      bzj = sqrt(tVars%eta) / 2 / tVars%rho * &
            ( -tVars%ux*sx * sin(szt)  + tVars%uy*sy * cos(szt) )
!$OMP END WORKSHARE

    end if

!    END elliptical undulator description
!  ####################################################




  end if

  end if

end subroutine getBZfield

end module bfields
