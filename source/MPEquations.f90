! Copyright 2012-2018, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause

module MPequations

  use paratype, only: wp, ip
  implicit none

contains

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Calculate d/dz of real pperp

  subroutine dppdz_r_f(sx, sy, sz2, spr, spi, sgam, tVars, sdpr)
    use typecalcParams, only: fcalcParams
  	implicit none
    real(kind=wp), contiguous, intent(in) :: sx(:), sy(:), sz2(:), spr(:), &
                                             spi(:), sgam(:)
    type(fcalcParams), intent(in) :: tVars
    real(kind=wp), contiguous, intent(out) :: sdpr(:)

!$OMP WORKSHARE
    sdpr = tVars%sInv2rho * ( tVars%n2col * tVars%byu  & 
                        - tVars%eta * tVars%sp2 / tVars%kappa**2 *    &
                        tVars%sField4ElecReal ) & 
           + tVars%kappa * spi / sgam * (1 + tVars%eta * tVars%sp2) &
               * tVars%n2col * tVars%bzu
!$OMP END WORKSHARE

  end subroutine dppdz_r_f

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Calculate d/dz of imag pperp

  subroutine dppdz_i_f(sx, sy, sz2, spr, spi, sgam, tVars, sdpi)
    use typecalcParams, only: fcalcParams
    implicit none
    real(kind=wp), contiguous, intent(in) :: sx(:), sy(:), sz2(:), spr(:), &
                                             spi(:), sgam(:) 
    type(fcalcParams), intent(in) :: tVars
    real(kind=wp), contiguous, intent(out) :: sdpi(:)

!$OMP WORKSHARE
    sdpi = tVars%sInv2rho * (  tVars%n2col * tVars%bxu  & 
           - tVars%eta * tVars%sp2 / tVars%kappa**2.0_wp * &
                        tVars%sField4ElecImag ) & 
           - tVars%kappa * spr / sgam * (1.0_wp + tVars%eta * tVars%sp2) &
               * tVars%n2col * tVars%bzu
!$OMP END WORKSHARE

  end subroutine dppdz_i_f

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Calculate d/dz of scaled energy Gamma

  subroutine dgamdz_f(sx, sy, sz2, spr, spi, sgam, tVars, sdgam)
    use typecalcParams, only: fcalcParams
    implicit none
    real(kind=wp), contiguous, intent(in) :: sx(:), sy(:), sz2(:), spr(:), &
                                             spi(:), sgam(:)
    type(fcalcParams), intent(in) :: tVars
    real(kind=wp), contiguous, intent(out) :: sdgam(:)

!$OMP WORKSHARE
    sdgam = -tVars%rho * ( 1.0_wp + tVars%eta * tVars%sp2 ) / sgam * 2.0_wp *   &
           ( spr * tVars%sField4ElecReal + spi * tVars%sField4ElecImag ) 
!$OMP END WORKSHARE

  end subroutine dgamdz_f

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Calculate d/dz of x-coord

  subroutine dxdz_f(sx, sy, sz2, spr, spi, sgam, tVars, sdx)
    use typecalcParams, only: fcalcParams
    implicit none
    real(kind=wp), contiguous,  intent(in) :: sx(:), sy(:), sz2(:), spr(:), &
                                              spi(:), sgam(:)
    type(fcalcParams), intent(in) :: tVars
    real(kind=wp), contiguous, intent(out) :: sdx(:)

!$OMP WORKSHARE
    sdx = 2.0_wp * tVars%rho * tVars%kappa / sqrt(tVars%eta) * &
          (1.0_wp + tVars%eta * tVars%sp2) / sgam *  &
          spr
!$OMP END WORKSHARE

  end subroutine dxdz_f

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Calculate d/dz of y-coord

  subroutine dydz_f(sx, sy, sz2, spr, spi, sgam, tVars, sdy)
    use typecalcParams, only: fcalcParams
    implicit none
    real(kind=wp), contiguous, intent(in) :: sx(:), sy(:), sz2(:), spr(:), &
                                             spi(:), sgam(:)
    type(fcalcParams), intent(in) :: tVars
    real(kind=wp), contiguous, intent(out) :: sdy(:)

!$OMP WORKSHARE
    sdy = - 2.0_wp * tVars%rho * tVars%kappa / sqrt(tVars%eta) * &
          (1.0_wp + tVars%eta * tVars%sp2) / sgam *  &
          spi
!$OMP END WORKSHARE

  end subroutine dydz_f

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Calculate d/dz of z2-coord

  subroutine dz2dz_f(sx, sy, sz2, spr, spi, sgam, tVars, sdz2)
    use typecalcParams, only: fcalcParams
    implicit none
    real(kind=wp), contiguous, intent(in) :: sx(:), sy(:), sz2(:), spr(:), &
                                             spi(:), sgam(:)
    type(fcalcParams), intent(in) :: tVars
    real(kind=wp), contiguous, intent(out) :: sdz2(:)

!$OMP WORKSHARE
    sdz2 = tVars%sp2
!$OMP END WORKSHARE

  end subroutine dz2dz_f

end module MPequations

