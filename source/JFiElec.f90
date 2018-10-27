! Copyright 2012-2018, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause

!> @author
!> Lawrence Campbell,
!> University of Strathclyde, 
!> Glasgow, UK
!> @brief
!> Module containing routines dealing with the interpolation of the 
!> macroparticles to the 3D field mesh.


module FiElec

use paratype, only: wp, ip
implicit none

contains

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Get interpolants for interpolating between macroparticles and mesh.
!> @param[in] sx X-coords of macroparticles
!> @param[in] sx Y-coords of macroparticles
!> @param[in] sx Z2-coords of macroparticles
!> @param[inout] tVars Custom Fortran type

  subroutine getInterps_3D(sx, sy, sz2, tVars)
    use typeCalcParams, only: fcalcParams
    use globals, only: iTemporal
    implicit none
    real(kind=wp), intent(in) :: sx(:), sy(:), sz2(:)
    type(fcalcParams), intent(inout) :: tVars

    integer(kind=ip) :: xnode, ynode, z2node
    integer(kind=ip) :: i

    real(kind=wp) :: locx, locy, locz2, &
                     x_in1, x_in2, y_in1, y_in2, z2_in1, z2_in2

!$OMP DO PRIVATE(xnode, ynode, z2node, locx, locy, locz2, &
!$OMP x_in1, x_in2, y_in1, y_in2, z2_in1, z2_in2)
    do i = 1, tVars%maxEl
      if (i<=tVars%iNMPsLoc) then 

!                  Get surrounding nodes 

        xnode = floor( (sx(i) + tVars%halfx ) / tVars%dx)  + 1_IP
        locx = sx(i) + tVars%halfx - real(xnode  - 1_IP, kind=wp) * tVars%dx
        x_in2 = locx / tVars%dx
        x_in1 = (1.0_wp - x_in2)

        ynode = floor( (sy(i) + tVars%halfy )  / tVars%dy)  + 1_IP
        locy = sy(i) + tVars%halfy - real(ynode  - 1_IP, kind=wp) * tVars%dy
        y_in2 = locy / tVars%dy
        y_in1 = (1.0_wp - y_in2)

        z2node = floor(sz2(i)  / tVars%dz2)  + 1_IP
        locz2 = sz2(i) - real(z2node  - 1_IP, kind=wp) * tVars%dz2
        z2_in2 = locz2 / tVars%dz2
        z2_in1 = (1.0_wp - z2_in2)

        if ((xnode >= tVars%inNX) .or. (xnode < 1))then
          tVars%qInnerXYOK = .false.
          tVars%qPArrOK = .false.
        end if

        if ((ynode >= tVars%inNY) .or. (ynode < 1)) then
          tVars%qInnerXYOK = .false.
          tVars%qPArrOK = .false.
        end if

        if (tVars%fieldMesh == itemporal) then
          if (z2node >= tVars%NZ2) then
            print*, 'Z2 coord is too large!! with node:', z2node, &
                    ' and pos ', sz2(i)
          STOP
          end if
        end if

        if (z2node >= tVars%bz2) then
          tVars%qPArrOK = .false.
        end if

!                  Get weights for interpolant

        tVars%lis_GR(1,i) = x_in1 * y_in1 * z2_in1
        tVars%lis_GR(2,i) = x_in2 * y_in1 * z2_in1
        tVars%lis_GR(3,i) = x_in1 * y_in2 * z2_in1
        tVars%lis_GR(4,i) = x_in2 * y_in2 * z2_in1
        tVars%lis_GR(5,i) = x_in1 * y_in1 * z2_in2
        tVars%lis_GR(6,i) = x_in2 * y_in1 * z2_in2
        tVars%lis_GR(7,i) = x_in1 * y_in2 * z2_in2
        tVars%lis_GR(8,i) = x_in2 * y_in2 * z2_in2

      end if
    end do
!$OMP END DO

  end subroutine getInterps_3D

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Get field for electron macroparticles.
!> @param[in] sAr Real field at macroparticle positions
!> @param[in] sAi Imaginary field at macroparticle positions
!> @param[inout] tVars Custom Fortran type with data

  subroutine getFFelecs_3D(sAr, sAi, tVars)
    use typeCalcParams, only: fcalcParams
    implicit none
    real(kind=wp), contiguous, intent(in) :: sAr(:), sAi(:)
    type(fcalcParams), intent(inout) :: tVars
    integer(kind=ip) :: i

!$OMP DO
    do i = 1, tVars%iNMPsLoc

      tVars%sField4ElecReal(i) = tVars%lis_GR(1,i) * sAr(tVars%p_nodes(i)) + tVars%sField4ElecReal(i)
      tVars%sField4ElecReal(i) = tVars%lis_GR(2,i) * sAr(tVars%p_nodes(i) + 1_ip) + tVars%sField4ElecReal(i)
      tVars%sField4ElecReal(i) = tVars%lis_GR(3,i) * sAr(tVars%p_nodes(i) + tVars%inNX) + tVars%sField4ElecReal(i)
      tVars%sField4ElecReal(i) = tVars%lis_GR(4,i) * sAr(tVars%p_nodes(i) + tVars%inNX + 1_ip) + tVars%sField4ElecReal(i)
      tVars%sField4ElecReal(i) = tVars%lis_GR(5,i) * sAr(tVars%p_nodes(i) + tVars%inNN) + tVars%sField4ElecReal(i)
      tVars%sField4ElecReal(i) = tVars%lis_GR(6,i) * sAr(tVars%p_nodes(i) + tVars%inNN + 1_ip) + tVars%sField4ElecReal(i)
      tVars%sField4ElecReal(i) = tVars%lis_GR(7,i) * sAr(tVars%p_nodes(i) + tVars%inNN + tVars%inNX) + tVars%sField4ElecReal(i)
      tVars%sField4ElecReal(i) = tVars%lis_GR(8,i) * sAr(tVars%p_nodes(i) + tVars%inNN + tVars%inNX + 1) + tVars%sField4ElecReal(i)
  
      tVars%sField4ElecImag(i) = tVars%lis_GR(1,i) * sAi(tVars%p_nodes(i)) + tVars%sField4ElecImag(i)
      tVars%sField4ElecImag(i) = tVars%lis_GR(2,i) * sAi(tVars%p_nodes(i)  + 1_ip) + tVars%sField4ElecImag(i)
      tVars%sField4ElecImag(i) = tVars%lis_GR(3,i) * sAi(tVars%p_nodes(i)  + tVars%inNX) + tVars%sField4ElecImag(i)
      tVars%sField4ElecImag(i) = tVars%lis_GR(4,i) * sAi(tVars%p_nodes(i)  + tVars%inNX + 1_ip) + tVars%sField4ElecImag(i)
      tVars%sField4ElecImag(i) = tVars%lis_GR(5,i) * sAi(tVars%p_nodes(i)  + tVars%inNN) + tVars%sField4ElecImag(i)
      tVars%sField4ElecImag(i) = tVars%lis_GR(6,i) * sAi(tVars%p_nodes(i)  + tVars%inNN + 1_ip) + tVars%sField4ElecImag(i)
      tVars%sField4ElecImag(i) = tVars%lis_GR(7,i) * sAi(tVars%p_nodes(i)  + tVars%inNN + tVars%inNX) + tVars%sField4ElecImag(i)
      tVars%sField4ElecImag(i) = tVars%lis_GR(8,i) * sAi(tVars%p_nodes(i)  + tVars%inNN + tVars%inNX + 1) + tVars%sField4ElecImag(i)
  
    end do 
!$OMP END DO

  end subroutine getFFelecs_3D

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Get da/dz
!> @param[in] sDADzr d/dz of real field
!> @param[in] sDADzi d/dz of imaginary field
!> @param[in] spr Real component of pperp
!> @param[in] spi Imaginary component of pperp
!> @param[in] sgam Scaled electron energy \f% \Gamma \f$
!> @param[in] wts Scaled macroparticle weights
!> @param[inout] tVars Custom Fortran type with data

  subroutine getSource_3D(sDADzr, sDADzi, spr, spi, sgam, wts, tVars)
    use typeCalcParams, only: fcalcParams
    implicit none
    real(kind=wp), intent(inout) :: sDADzr(:), sDADzi(:)
    real(kind=wp), intent(in) :: spr(:), spi(:)
    real(kind=wp), intent(in) :: sgam(:)
    real(kind=wp), intent(in) :: wts(:)
    type(fcalcParams), intent(in) :: tVars

    integer(kind=ip) :: i
    real(kind=wp) :: dadzRInst, dadzIInst

!$OMP DO PRIVATE(dadzRInst, dadzIInst)
    do i = 1, tVars%maxEl
  
      if (i<=tVars%iNMPsLoc) then

!                  Get 'instantaneous' dAdz

        dadzRInst = ((wts(i)/tVars%dV3) * (1 + tVars%eta * tVars%sp2(i) ) &
                          * spr(i) / sgam(i) )
    
        !$OMP ATOMIC
        sDADzr(tVars%p_nodes(i)) =                         &
          tVars%lis_GR(1,i) * dadzRInst + sDADzr(tVars%p_nodes(i))
      
        !$OMP ATOMIC
        sDADzr(tVars%p_nodes(i) + 1_ip) =                  &
          tVars%lis_GR(2,i) * dadzRInst + sDADzr(tVars%p_nodes(i) + 1_ip)                

        !$OMP ATOMIC
        sDADzr(tVars%p_nodes(i) + tVars%inNX) =           &
          tVars%lis_GR(3,i) * dadzRInst + sDADzr(tVars%p_nodes(i) + tVars%inNX)          

        !$OMP ATOMIC
        sDADzr(tVars%p_nodes(i) + tVars%inNX + 1_ip) =    &
          tVars%lis_GR(4,i) * dadzRInst + sDADzr(tVars%p_nodes(i) + tVars%inNX + 1_ip)   

        !$OMP ATOMIC
        sDADzr(tVars%p_nodes(i) + tVars%inNN) =                &
          tVars%lis_GR(5,i) * dadzRInst + sDADzr(tVars%p_nodes(i) + tVars%inNN)               

        !$OMP ATOMIC
        sDADzr(tVars%p_nodes(i) + tVars%inNN + 1_ip) =         &
          tVars%lis_GR(6,i) * dadzRInst + sDADzr(tVars%p_nodes(i) + tVars%inNN + 1_ip)         

        !$OMP ATOMIC
        sDADzr(tVars%p_nodes(i) + tVars%inNN + tVars%inNX) =  &
          tVars%lis_GR(7,i) * dadzRInst + sDADzr(tVars%p_nodes(i) + tVars%inNN + tVars%inNX)   

        !$OMP ATOMIC
        sDADzr(tVars%p_nodes(i) + tVars%inNN + tVars%inNX + 1_ip) = &
          tVars%lis_GR(8,i) * dadzRInst + sDADzr(tVars%p_nodes(i) + tVars%inNN + tVars%inNX + 1_ip)

!                   Imaginary part

        dadzIInst = ((wts(i) / tVars%dV3) * (1 + tVars%eta * tVars%sp2(i) ) &
                          * spi(i) / sgam(i) ) 

        !$OMP ATOMIC
        sDADzi(tVars%p_nodes(i)) =                             & 
          tVars%lis_GR(1,i) * dadzIInst + sDADzi(tVars%p_nodes(i))                        

        !$OMP ATOMIC
        sDADzi(tVars%p_nodes(i) + 1_ip) =                      & 
          tVars%lis_GR(2,i) * dadzIInst + sDADzi(tVars%p_nodes(i) + 1_ip)           

        !$OMP ATOMIC
        sDADzi(tVars%p_nodes(i) + tVars%inNX) =               & 
          tVars%lis_GR(3,i) * dadzIInst + sDADzi(tVars%p_nodes(i) + tVars%inNX)           

        !$OMP ATOMIC
        sDADzi(tVars%p_nodes(i) + tVars%inNX + 1_ip) =        & 
          tVars%lis_GR(4,i) * dadzIInst + sDADzi(tVars%p_nodes(i) + tVars%inNX + 1_ip)    

        !$OMP ATOMIC
        sDADzi(tVars%p_nodes(i) + tVars%inNN) =                    & 
          tVars%lis_GR(5,i) * dadzIInst + sDADzi(tVars%p_nodes(i) + tVars%inNN)               

        !$OMP ATOMIC
        sDADzi(tVars%p_nodes(i) + tVars%inNN + 1_ip) =             & 
          tVars%lis_GR(6,i) * dadzIInst + sDADzi(tVars%p_nodes(i) + tVars%inNN + 1_ip)       

        !$OMP ATOMIC
        sDADzi(tVars%p_nodes(i) + tVars%inNN + tVars%inNX) =      & 
          tVars%lis_GR(7,i) * dadzIInst + sDADzi(tVars%p_nodes(i) + tVars%inNN + tVars%inNX)  

        !$OMP ATOMIC
        sDADzi(tVars%p_nodes(i) + tVars%inNN + tVars%inNX + 1_ip) =  & 
          tVars%lis_GR(8,i) * dadzIInst + sDADzi(tVars%p_nodes(i) + tVars%inNN + tVars%inNX + 1_ip)

      end if
  
    end do 
!$OMP END DO

  end subroutine getSource_3D

end module FiElec
