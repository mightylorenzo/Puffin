! ################################################
! Copyright 2012-2018, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause
! ################################################

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

  type, public :: fbfield
    real(kind=wp) :: taper, alpha0
    real(kind=wp) :: sZFS, sZFE
    real(kind=wp) :: kbxSF, kbySF, kxu, kyu, ux, uy
    character(32_ip) :: zUndType
    integer(kind=ip), private :: iUndPlace
    logical :: qUndEnds
    logical :: qSF
  contains
    procedure :: getBFields
    procedure, private :: getBXfield
    procedure, private :: getBYfield
    procedure, private :: getBZfield
    procedure, private :: getActive
    procedure :: getAlpha
  end type fbfield

contains

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Calculate the b-fields. Calls bx, by and bz subroutines.
!> @param[in] sx Macroparticle coordinates in xbar
!> @param[in] sy Macroparticle coordinates in ybar
!> @param[in] sz current scaled distance though this undulator module, \f$ \bar{z} \f$
!> @param[out] bxj scaled b-field in x direction for macroparticles
!> @param[out] byj scaled b-field in y direction for macroparticles
!> @param[out] bzj scaled b-field in z direction for macroparticles

  subroutine getBFields(self, sx, sy, sZ, tScale, &
                        bxj, byj, bzj)
    use typeScale, only: fScale
    implicit none
    class(fbfield), intent(in) :: self
    real(kind=wp), contiguous, intent(in) :: sx(:), sy(:)
    real(kind=wp), intent(in) :: sz
    type(fScale), intent(in) :: tScale
    real(kind=wp), contiguous, intent(out) :: bxj(:), byj(:), bzj(:)
    real(kind=wp) :: alph

    call self%adjUndPlace(sz)
    alph = self%getAlpha(sz)

    call self%getBXfield(sx, sy, sz, tScale, bxj)
    call self%getBYfield(sx, sy, sz, tScale, byj)
    call self%getBZfield(sx, sy, sz, tScale, bzj)
  
    bxj = alph * bxj
    byj = alph * byj
    bzj = alph * bzj

  end subroutine getBFields

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Determine whether we are in the ends or the active section of undulator.

  subroutine adjUndPlace(self, szl)
    implicit none
    class(fbfields), intent(inout) :: self
    real(kind=wp), intent(in) :: szl

    if (self%qUndEnds) then
      if (szl < 0) then
        print*, 'undulator section not recognised, sz < 0!!'
        stop
      else if (sZl <= self%sZFS) then
        iUndPlace_G = iUndStart_G
      else if (sZl >= self%sZFE) then
        iUndPlace_G = iUndEnd_G
      else if ((sZl > self%sZFS) .and. (sZl < self%sZFE)) then
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

  function getAlpha(self, sZ) result(alpha)
    implicit none
    real(kind=wp), intent(in) :: sZ
    real(kind=wp) :: alpha
    class(fbfields), intent(in) :: self
    real(kind=wp), intent(out) :: n2col

    if ((sZ >= self%sZFS) .and. (sZ <= self%sZFE)) then
      alpha = self%alpha0  + self%undgrad*(sz - self%sZFS)  ! linear taper
    else if (sZ < self%sZFS) then
      alpha = self%alpha0
    else if (sZ > self%sZFE) then
      alpha = self%alpha0  + self%undgrad*(self%sZFE - self%sZFS)
    end if
    return
  end function getAlpha

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Get scaled bx field \f$ b_x \f$ at input coordinates.

  subroutine getBXfield(self, sx, sy, sz, tScale, bxj)
    use typeScale, only: fScale
    use typesAndConstants, only: pi
    implicit none
    class(fbfields), intent(in) :: self
    real(kind=wp), contiguous, intent(in) :: sx(:), sy(:)
    real(kind=wp), intent(in) :: sz
    type(fScale), intent(in) :: tScale
    real(kind=wp), contiguous, intent(out) :: bxj(:)
    real(kind=wp) :: szt

    szt = sZ
    szt = szt / 2_wp / tScale%rho

!  ####################################################
!    Curved pole case - planar wiggler with focusing
!    in both x and y (electron wiggles in x)

    if (self%zUndType == 'curved') then
      if (iUndPlace_G == iUndStart_G) then
!$OMP WORKSHARE
        bxj = self%kxu / self%kyu * sinh(self%kxu * sx) &
              * sinh(self%kyu * sy) &
              * szt / 4_wp / pi * sin(szt)
!$OMP END WORKSHARE
      else if (iUndPlace_G == iUndEnd_G) then
        szt = sZ - self%sZFE
        szt = szt / 2_wp / tScale%rho
!$OMP WORKSHARE
        bxj = self%kxu / self%kyu * sinh(self%kxu * sx) &
              * sinh(self%kyu * sy) &
              * (-szt / 4_wp / pi + 1_wp) * sin(szt)
!$OMP END WORKSHARE
      else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
        bxj = self%kxu / self%kyu * sinh(self%kxu * sx) &
              * sinh(self%kyu * sy) &
              * sin(szt)
!$OMP END WORKSHARE

      end if

!    END curved pole field description
!  ####################################################
!  ####################################################
!    Plane-pole case - planar wiggler with focusing
!    only in y (and electron will wiggle in x)

else if (self%zUndType == 'planepole')  then
      if (iUndPlace_G == iUndStart_G) then
!$OMP WORKSHARE
        bxj = 0_wp
!$OMP END WORKSHARE
      else if (iUndPlace_G == iUndEnd_G) then
        szt = sZ - self%sZFE
        szt = szt / 2_wp / tScale%rho
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

else if (self%zUndType == 'helical')  then
      if (iUndPlace_G == iUndStart_G) then
!$OMP WORKSHARE
        bxj = (sZ - pi * tScale%rho) / (6_wp * pi*tScale%rho) * cos(szt)
!$OMP END WORKSHARE
        if (sZ < pi * tScale%rho) then 
!$OMP WORKSHARE
          bxj = 0.0_wp
!$OMP END WORKSHARE
else if (sZ > 7_wp * pi * tScale%rho) then
!$OMP WORKSHARE
          bxj = cos(szt)
!$OMP END WORKSHARE
        end if
      else if (iUndPlace_G == iUndEnd_G) then
        szt = sZ - self%sZFE
!$OMP WORKSHARE
        bxj =  - (szt - 7.0_wp * pi * tScale%rho) / (6_wp * pi * tScale%rho) * &
                 cos(szt / 2_wp / tScale%rho)
!$OMP END WORKSHARE
        if (sZt < pi * tScale%rho) then 
!$OMP WORKSHARE
          bxj = cos(szt / 2_wp / tScale%rho)
!$OMP END WORKSHARE
else if (sZt > 7_wp * pi * tScale%rho) then
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
!$OMP WORKSHARE
        bxj = self%ux * (sZ - pi * tScale%rho) / (6_wp * pi*tScale%rho) * cos(szt)
!$OMP END WORKSHARE
        if (sZ < pi * tScale%rho) then 
!$OMP WORKSHARE
          bxj = 0.0_wp
!$OMP END WORKSHARE
else if (sZ > 7_wp * pi * tScale%rho) then
!$OMP WORKSHARE
          bxj = self%ux * cos(szt)
!$OMP END WORKSHARE
        end if
      else if (iUndPlace_G == iUndEnd_G) then
        szt = sZ - self%sZFE
!$OMP WORKSHARE
        bxj =  -self%ux * (szt - 7.0_wp * pi * tScale%rho) / (6_wp * pi * tScale%rho) * &
                 cos(szt / 2_wp / tScale%rho)
!$OMP END WORKSHARE
        if (sZt < pi * tScale%rho) then 
!$OMP WORKSHARE
          bxj = self%ux * cos(szt / 2_wp / tScale%rho)
!$OMP END WORKSHARE
else if (sZt > 7_wp * pi * tScale%rho) then
!$OMP WORKSHARE
          bxj = 0.0_wp
!$OMP END WORKSHARE
        end if
      else if (iUndPlace_G == iUndMain_G) then
!$OMP WORKSHARE
        bxj = self%ux*cos(szt)
!$OMP END WORKSHARE
      end if
!    END elliptical undulator description
!  ####################################################
    end if
!   Strong, in-undulator focusing component (non-physical)
    if (self%qSF) then
!$OMP WORKSHARE
        bxj = sqrt(tScale%eta) * self%kbySF**2.0_wp / tScale%kappa &
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

  subroutine getBYfield(self, sx, sy, sz, tScale, byj)

  use typeScale, only: fScale
  use typesAndConstants, only: pi
  implicit none
  class(fbfields), intent(in) :: self
  real(kind=wp), contiguous, intent(in) :: sx(:), sy(:)
  real(kind=wp), intent(in) :: sz
  type(fScale), intent(in) :: tScale
  real(kind=wp), contiguous, intent(out) :: byj(:)

!    Local vars:-

  real(kind=wp) :: szt

  szt = sZ
  szt = szt / 2_wp / tScale%rho


!  ####################################################
!    Curved pole case - planar wiggler with focusing
!    in both x and y (electron wiggles in x)


  if (self%zUndType == 'curved') then

    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
      byj = cosh(self%kxu * sx) &
            * cosh(self%kyu * sy) &
            *  szt / 4_wp / pi * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - self%sZFE
      szt = szt / 2_wp / tScale%rho

!$OMP WORKSHARE
      byj = cosh(self%kxu * sx) &
            * cosh(self%kyu * sy) &
            * (-szt / 4_wp / pi + 1_wp) * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      byj = cosh(self%kxu * sx) &
            * cosh(self%kyu * sy) &
            * sin(szt)
!$OMP END WORKSHARE

    end if

!    END curved pole field description
!  ####################################################







!  ####################################################
!    Plane-pole case - planar wiggler with focusing
!    only in y (and electron will wiggle in x)



else if (self%zUndType == 'planepole')  then

    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
!      byj = cosh( sqrt(tScale%eta) / 2_wp / tScale%rho * sy) * &
!            (  (- sin(szt / 8_wp) * &
!               cos(szt / 8_wp) * cos(szt) / 4_wp   &
!             +  sin(szt/8_wp)**2_wp  * sin(szt) )  )


      !byj = sin(szt/8_wp)**2_wp * sin(szt)
      byj = szt / 4_wp / pi * sin(szt)
!$OMP END WORKSHARE

!    print*, "hehehe"

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - self%sZFE
      szt = szt / 2_wp / tScale%rho

!$OMP WORKSHARE
!      byj = cosh( sqrt(tScale%eta) / 2_wp / tScale%rho * sy) * &
!            (  cos(szt / 8_wp) * &
!              sin(szt / 8_wp) * cos(szt)  / 4_wp  &
!            +  cos(szt/8_wp)**2_wp  * sin(szt)  )
            
      byj = (-szt / 4_wp / pi + 1_wp) * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
!      byj = cosh( sqrt(tScale%eta) / 2_wp / tScale%rho * sy) &
!            * sin(szt)
      byj = sin(szt)
!$OMP END WORKSHARE

    end if

!    END plane pole undulator field description
!  ####################################################









!  ####################################################
!    Helical case - helical wiggler with focusing
!    in x and y (and electron will wiggle in x and y)

else if (self%zUndType == 'helical')  then

    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
      byj = szt / 4_wp / pi * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - self%sZFE
      szt = szt / 2_wp / tScale%rho

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
      byj = self%uy * szt / 4_wp / pi * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - self%sZFE
      szt = szt / 2_wp / tScale%rho

!$OMP WORKSHARE
      byj = self%uy * (-szt / 4_wp / pi + 1_wp) * sin(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      byj = self%uy * sin(szt)
!$OMP END WORKSHARE

    end if

!    END elliptical undulator description
!  ####################################################




  end if

!   Focusing component (non-physical)

    if (self%qSF) then

!$OMP WORKSHARE
      byj = -sqrt(tScale%eta) * self%kbxSF**2.0_wp / tScale%kappa &
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

subroutine getBZfield(self, sx, sy, sz, tScale, bzj)
  use typeScale, only: fScale
  use typesAndConstants, only: pi
  implicit none
  class(fScale), intent(in) :: self
  real(kind=wp), contiguous, intent(in) :: sx(:), sy(:)
  real(kind=wp), intent(in) :: sz
  type(fScale), intent(in) :: tScale
  real(kind=wp), contiguous, intent(out) :: bzj(:)

!    Local vars:-

  real(kind=wp) :: szt

  szt = sZ
  szt = szt / 2.0_wp / tScale%rho

  if (tScale%q1D) then

!  ####################################################
!    1D case - no z-component of magnetic field

!$OMP WORKSHARE
    bzj = 0.0_wp
!$OMP END WORKSHARE

  else

!  ####################################################
!    Curved pole case - planar wiggler with focusing
!    in both x and y (electron wiggles in x)


  if (self%zUndType == 'curved') then

    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
      bzj = sqrt(tScale%eta) / 2 / tScale%rho / self%kyu &
                * cosh(self%kxu * sx) &
            * sinh(self%kyu * sy) &
            * (szt / 4_wp / pi) * cos(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - self%sZFE
      szt = szt / 2_wp / tScale%rho

!$OMP WORKSHARE
      bzj = sqrt(tScale%eta) / 2_wp / tScale%rho / self%kxu &
            * cosh(self%kxu * sx) * sinh(self%kyu * sy) &
            * (-szt / 4_wp / pi + 1_wp) * cos(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then


!$OMP WORKSHARE
      bzj = sqrt(tScale%eta) / 2_wp / tScale%rho / self%kxu * &
           cosh(self%kxu * sx) * sinh(self%kyu * sy) &
            * cos(szt)
!$OMP END WORKSHARE

    end if

!    END curved pole field description
!  ####################################################






!  ####################################################
!    Plane-pole case - planar wiggler with focusing
!    only in y (and electron will wiggle in x)

else if (self%zUndType == 'planepole')  then

    if (iUndPlace_G == iUndStart_G) then

!$OMP WORKSHARE
      bzj = sinh( sqrt(tScale%eta) / 2_wp / tScale%rho * sy) * &
            (szt / 4_wp / pi) * cos(szt)
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - self%sZFE
      szt = szt / 2_wp / tScale%rho

!$OMP WORKSHARE
      bzj = sinh( sqrt(tScale%eta) / 2_wp / tScale%rho * sy) * &
            (-szt / 4_wp / pi + 1_wp) * cos(szt) 
!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      bzj = sinh( sqrt(tScale%eta) / 2_wp / tScale%rho * sy) &
            * cos(szt)
!$OMP END WORKSHARE

    end if

!    END plane pole undulator field description
!  ####################################################






!  ####################################################
!    Helical case - helical wiggler with focusing
!    in x and y (and electron will wiggle in x and y)

else if (self%zUndType == 'helical')  then

    if (iUndPlace_G == iUndStart_G) then

! ...from x-comp:

!$OMP WORKSHARE
      bzj = - sqrt(tScale%eta) / 2 / tScale%rho * (sZ - pi * tScale%rho) / (6_wp * pi*tScale%rho) &
            * sx * sin(szt)            
!$OMP END WORKSHARE

      if (sZ < pi * tScale%rho) then 
!$OMP WORKSHARE
        bzj = 0.0_wp
!$OMP END WORKSHARE
else if (sZ > 7_wp * pi * tScale%rho) then
!$OMP WORKSHARE
        bzj = - sqrt(tScale%eta) / 2 / tScale%rho * sx * sin(szt)
!$OMP END WORKSHARE
      end if
      

! ...and from y-comp:

!$OMP WORKSHARE
      bzj = bzj + szt / 4_wp / pi * sqrt(tScale%eta) / 2 / tScale%rho * sy * cos(szt)
!$OMP END WORKSHARE

!  !$OMP WORKSHARE
!        bzj = sqrt(tScale%eta) / 2 / tScale%rho * (     &
!              sx *  sin(szt) )    + &
!              sy * ( -1/32_wp * cos(szt/4_wp) * cos(szt) + &
!                      1/4_wp * sin(szt/4_wp) * sin(szt) + &
!                      sin(szt/8_wp)**2 * cos(szt) ) )
!  !$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - self%sZFE
      !szt = szt / 2_wp / tScale%rho

! ...from x-comp:

!$OMP WORKSHARE
      bzj = - sqrt(tScale%eta) / 2 / tScale%rho * (szt - 7.0_wp * pi * tScale%rho) / &
              (6_wp * pi * tScale%rho) * sx * sin(szt / 2_wp / tScale%rho)            
!$OMP END WORKSHARE

      if (szt < pi * tScale%rho) then 
!$OMP WORKSHARE
        bzj = - sqrt(tScale%eta) / 2 / tScale%rho * sx * sin(szt / 2_wp / tScale%rho)
!$OMP END WORKSHARE
else if (szt > 7_wp * pi * tScale%rho) then
!$OMP WORKSHARE
        bzj = 0.0_wp
!$OMP END WORKSHARE
      end if
      

! ...and from y-comp:

!$OMP WORKSHARE
      bzj = bzj + (-szt / 8_wp / pi / tScale%rho + 1_wp) * &
            sqrt(tScale%eta) / 2 / tScale%rho * sy * cos(szt / 2_wp / tScale%rho)
!$OMP END WORKSHARE












! !$OMP WORKSHARE
!       bzj = sqrt(tScale%eta) / 2 / tScale%rho * (     &
!             sx * ( -1/32_wp * cos(szt/4_wp) * sin(szt) - &
!                     1/4_wp * sin(szt/4_wp) * cos(szt) - &
!                     cos(szt/8_wp)**2 * sin(szt) )    + &
!             sy * ( 1/32_wp * cos(szt/4_wp) * cos(szt) - &
!                     1/4_wp * sin(szt/4_wp) * sin(szt) + &
!                     cos(szt/8_wp)**2 * cos(szt) ) )
! !$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      bzj = sqrt(tScale%eta) / 2 / tScale%rho * &
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
      bzj = - sqrt(tScale%eta) / 2 / tScale%rho * (sZ - pi * tScale%rho) / (6_wp * pi*tScale%rho) &
            * self%ux * sx * sin(szt)            
!$OMP END WORKSHARE

      if (sZ < pi * tScale%rho) then 
!$OMP WORKSHARE
        bzj = 0.0_wp
!$OMP END WORKSHARE
else if (sZ > 7_wp * pi * tScale%rho) then
!$OMP WORKSHARE
        bzj = - sqrt(tScale%eta) / 2 / tScale%rho * self%ux * sx * sin(szt)
!$OMP END WORKSHARE
      end if
      

! ...and from y-comp:

!$OMP WORKSHARE
      bzj = bzj + szt / 4_wp / pi * sqrt(tScale%eta) / 2 / tScale%rho * self%uy * sy * cos(szt)
!$OMP END WORKSHARE





!!$OMP WORKSHARE
!      bzj = sqrt(tScale%eta) / 2 / tScale%rho * (     &
!        self%ux*sx * ( 1/32_wp * cos(szt/4_wp) * sin(szt) + &
!                    1/4_wp * sin(szt/4_wp) * cos(szt) - &
!                    sin(szt/8_wp)**2 * sin(szt) )    + &
!        self%uy*sy * ( -1/32_wp * cos(szt/4_wp) * cos(szt) + &
!                    1/4_wp * sin(szt/4_wp) * sin(szt) + &
!                    sin(szt/8_wp)**2 * cos(szt) ) )
!!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndEnd_G) then

      szt = sZ - self%sZFE
!      szt = szt / 2_wp / tScale%rho








! ...from x-comp:

!$OMP WORKSHARE
      bzj = - sqrt(tScale%eta) / 2 / tScale%rho * (szt - 7.0_wp * pi * tScale%rho) / &
              (6_wp * pi * tScale%rho) * self%ux * sx * sin(szt / 2_wp / tScale%rho)            
!$OMP END WORKSHARE

      if (szt < pi * tScale%rho) then 
!$OMP WORKSHARE
        bzj = - sqrt(tScale%eta) / 2 / tScale%rho * self%ux * sx * sin(szt / 2_wp / tScale%rho)
!$OMP END WORKSHARE
else if (szt > 7_wp * pi * tScale%rho) then
!$OMP WORKSHARE
        bzj = 0.0_wp
!$OMP END WORKSHARE
      end if
      

! ...and from y-comp:

!$OMP WORKSHARE
      bzj = bzj + (-szt / 8_wp / pi / tScale%rho + 1_wp) * &
            sqrt(tScale%eta) / 2 / tScale%rho * self%uy * sy * cos(szt / 2_wp / tScale%rho)
!$OMP END WORKSHARE













!!$OMP WORKSHARE
!      bzj = sqrt(tScale%eta) / 2 / tScale%rho * (     &
!        self%ux*sx * ( -1/32_wp * cos(szt/4_wp) * sin(szt) - &
!                    1/4_wp * sin(szt/4_wp) * cos(szt) - &
!                    cos(szt/8_wp)**2 * sin(szt) )    + &
!        self%uy*sy * ( 1/32_wp * cos(szt/4_wp) * cos(szt) - &
!                    1/4_wp * sin(szt/4_wp) * sin(szt) + &
!                    cos(szt/8_wp)**2 * cos(szt) ) )
!!$OMP END WORKSHARE

    else if (iUndPlace_G == iUndMain_G) then

!$OMP WORKSHARE
      bzj = sqrt(tScale%eta) / 2 / tScale%rho * &
            ( -self%ux*sx * sin(szt)  + self%uy*sy * cos(szt) )
!$OMP END WORKSHARE

    end if

!    END elliptical undulator description
!  ####################################################




  end if

  end if

end subroutine getBZfield

end module bfields
