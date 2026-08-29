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

use puffin_kinds, only: WP, IP
use ArrayFunctions, only: tErrorLog_G, log_error
use Globals, only: NX_G, NY_G, NZ2_G, ntrndsi_G, nspinDX, nspinDY, sLengthOfElmX_G, &
  sLengthOfElmY_G, sLengthOfElmZ2_G, procelectrons_G, iNumberElectrons_G, qElectronsEvolve_G, &
  qFieldEvolve_G, qElectronFieldCoupling_G, s_chi_bar_G, fieldMesh, iTemporal
use Equations, only: alct_e_srtcts, &
  dalct_e_srtcts, adjundplace, sInv2rho, sp2, sField4ElecReal, sField4ElecImag, bxu, byu, bzu
use wigglerVar, only: getalpha
use rhs_kernels, only: prep_kernel_1D, prep_kernel_3D, source_kernel_1D, &
  source_kernel_3D, gather_kernel_1D, gather_kernel_3D, zero_efield_kernel, &
  equation_kernel
use ParaField, only: fz2, bz2, tTransInfo_G
use bfields, only: getbfields
use GlobalTypes, only: tUndulator, tFELFrame, tSimulationContext


use Functions, only: tProcInfo_G
implicit none (type, external)
private

public :: getrhs, IP, log_error, sp2, tErrorLog_G, tProcInfo_G, WP


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
!> @param[out] qOK Error flag
!> @param[in] und Undulator parameters (replaces undulator globals in chain)
!> @param[in] frame FEL frame parameters (replaces physics globals in chain)

  subroutine getrhs(sz, &
                    sAr, sAi, &
                    sx, sy, sz2, &
                    spr, spi, sgam, &
                    sdx, sdy, sdz2, &
                    sdpr, sdpi, sdgam, &
                    sDADzr, sDADzi, &
                    qOK, ctx)

  use rhs_vars, only: p_nodes, halfx, halfy, lis_GR, dx, dy, dz2, dV3, sInv2rho, sp2, &
    sField4ElecReal, sField4ElecImag, bxu, byu, bzu, WP, IP

  implicit none (type, external)

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

  real(kind=wp), contiguous,  intent(inout) :: sDADzr(:), sDADzi(:) !!!!!!!
  logical, intent(inout) :: qOK
  type(tSimulationContext), intent(inout) :: ctx

  logical :: qOKL
  logical :: qCoupleOK

!     Begin

  qOK = .false.
  qOKL = .false.

!     SETUP AND INITIALISE THE PARTICLE'S POSITION
!     ALLOCATE THE ARRAYS

!  allocate(Lj(iNumberElectrons_G))
  allocate(p_nodes(iNumberElectrons_G))
!  allocate(p_nodes2(ispt))
!  allocate(tmp1(500000))
!  allocate(tmp2(ispt))

  call alct_e_srtcts(iNumberElectrons_G)

  if (tTransInfo_G%qOneD) then
    allocate(lis_GR(2,iNumberElectrons_G))
  else
    allocate(lis_GR(8,iNumberElectrons_G))
  end if

  call rhs_tmsavers(sz, ctx%und, ctx%frame)  ! This can be moved later...

!     Adjust undulator tuning

  call getAlpha(sZ, ctx%und)
  call adjUndPlace(sZ, ctx%und)



!$OMP PARALLEL

!     Pass 1: p2, mesh node index, interpolation weights, bounds flags.
!     One worksharing loop for the lot -- see rhs_kernels.

  if (tTransInfo_G%qOneD) then

    call prep_kernel_1D(sz2, spr, spi, sgam, sp2, p_nodes, lis_GR,          &
                        ctx%frame%eta, ctx%frame%gamma_ref, ctx%frame%aw,   &
                        dz2, fz2, bz2, NZ2_G,                               &
                        (fieldMesh == iTemporal), ctx%flags,                &
                        procelectrons_G(1))

  else

    call prep_kernel_3D(sx, sy, sz2, spr, spi, sgam, sp2, p_nodes, lis_GR,  &
                        ctx%frame%eta, ctx%frame%gamma_ref, ctx%frame%aw,   &
                        dx, dy, dz2, halfx, halfy,                          &
                        nspinDX, nspinDY, ntrndsi_G, fz2, bz2, NZ2_G,       &
                        (fieldMesh == iTemporal), ctx%flags,                &
                        procelectrons_G(1))

  end if

!     The prep loop's closing barrier has run, so every thread now sees the
!     settled bounds flags. If any macroparticle fell outside the parallel
!     field bounds we skip the gather and the scatter for all of them, exactly
!     as the unfused chain did.

  if (tTransInfo_G%qOneD) then
    qCoupleOK = ctx%flags%parallel_arrays_ok
  else
    qCoupleOK = (ctx%flags%parallel_arrays_ok) .and. (ctx%flags%inner_xy_ok)
  end if

!     b-fields still cost three worksharing loops of their own. They only read
!     sx/sy, so hoisting them above the field coupling changes nothing.

  if (qElectronsEvolve_G) then

    call getBFields(sx, sy, sz, &
                    bxu, byu, bzu, ctx%und, ctx%frame)

  end if

!     Pass 2 is three worksharing loops that cost one barrier between them.
!     They are kept apart so each stays vectorisable: the scatter carries
!     !$OMP ATOMIC, the gather/zero choice is a branch, and the equations must
!     be branchless. See rhs_kernels for why fusing them was slower.
!
!     The NOWAIT chain relies on every loop here having the same trip count and
!     an explicit SCHEDULE(STATIC), so a thread reads back only its own writes.

!     Pass 2a: scatter this rank's source term onto the mesh.

  if (qCoupleOK) then

    if (tTransInfo_G%qOneD) then

      call source_kernel_1D(sDADzr, sDADzi, spr, spi, sgam, sp2,            &
                            s_chi_bar_G, p_nodes, lis_GR, dV3,              &
                            ctx%frame%eta, procelectrons_G(1))

    else

      call source_kernel_3D(sDADzr, sDADzi, spr, spi, sgam, sp2,            &
                            s_chi_bar_G, p_nodes, lis_GR, dV3,              &
                            ctx%frame%eta, nspinDX, ntrndsi_G,              &
                            procelectrons_G(1))

    end if

  end if

!     Pass 2b: gather the field onto the macroparticles. Exactly one of these
!     branches runs and each writes every element, so the arrays need no
!     separate pre-zeroing.

  if (qCoupleOK .and. qElectronFieldCoupling_G) then

    if (tTransInfo_G%qOneD) then

      call gather_kernel_1D(sAr, sAi, sField4ElecReal, sField4ElecImag,     &
                            p_nodes, lis_GR, procelectrons_G(1))

    else

      call gather_kernel_3D(sAr, sAi, sField4ElecReal, sField4ElecImag,     &
                            p_nodes, lis_GR, nspinDX, ntrndsi_G,            &
                            procelectrons_G(1))

    end if

  else

    call zero_efield_kernel(sField4ElecReal, sField4ElecImag, &
                            procelectrons_G(1))

  end if

!     Pass 2c: the six electron equations.

  if (qElectronsEvolve_G) then

    call equation_kernel(spr, spi, sgam, sp2,                               &
                         sField4ElecReal, sField4ElecImag,                  &
                         bxu, byu, bzu,                                     &
                         sdx, sdy, sdz2, sdpr, sdpi, sdgam,                 &
                         sInv2rho, ctx%frame%eta, ctx%frame%kappa,          &
                         ctx%frame%rho, ctx%und%n2col, procelectrons_G(1))

  end if

!$OMP END PARALLEL



!    if (qFieldEvolve_G) then

!     Sum dadz from different MPI processes together

!        call sum2RootArr(sDADz,ReducedNX_G*ReducedNY_G*NZ2_G*2,0)

!     Boundary condition dadz = 0 at head of field

!        if (tProcInfo_G%qRoot) sDADz(1:ReducedNX_G*ReducedNY_G) = 0.0_WP
!        if (tProcInfo_G%qRoot) sDADz(ReducedNX_G*ReducedNY_G*NZ2_G + 1: &
!                                     ReducedNX_G*ReducedNY_G*NZ2_G + &
!                                     ReducedNX_G*ReducedNY_G) = 0.0_WP

        !if (tTransInfo_G%qOneD) then
        !  if (tProcInfo_G%qRoot) sDADz=sDADz !sDADz=6.0_WP*sDADz
        !else
        !   if (tProcInfo_G%qRoot) sDADz=sDADz !216.0_WP/8.0_WP*sDADz
        !end if

!    end if

!     Switch field off

    if (.not. qFieldEvolve_G) then
       sDADzr = 0.0_WP
       sDADzi = 0.0_WP
    end if

!     if electrons not allowed to evolve then

    if (.not. qElectronsEvolve_G) then
       sdpr = 0.0_wp
       sdpi = 0.0_wp
       sdgam = 0.0_wp
       sdx   = 0.0_wp
       sdy   = 0.0_wp
       sdz2 = 0.0_wp
    end if

!     Deallocate arrays

!    deallocate(i_n4e,N,iNodeList_Re,iNodeList_Im,i_n4ered)
!    deallocate(sField4ElecReal,sField4ElecImag,Lj,dp2f)
    !deallocate(Lj)
    deallocate(lis_GR)
    deallocate(p_nodes)
    call dalct_e_srtcts()


    ! Set the error flag and exit

    qOK = .true.

    goto 2000

    call log_error("Error in rhs:getrhs",tErrorLog_G)
    print*,"Error in rhs:getrhs"
2000 continue

  end subroutine getrhs



!        #########################################


!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Initialize data used in the calculation of d/dz of electron beam + radiation
!> field quantities
!> @param[in] sz zbar
!> @param[in] und Undulator parameters
!> @param[in] frame FEL frame parameters

subroutine rhs_tmsavers(sz, und, frame)

use rhs_vars, only: iOutside, retim, ntrans, halfx, halfy, nc, nb, ZOver2rho, salphaSq, &
  sInv2rho, econst, un, dV3, dx, dy, dz2, qoutside, WP

real(kind=wp), intent(in) :: sz
type(tUndulator), intent(in) :: und
type(tFELFrame), intent(in) :: frame

  ioutside=0


!     Define the size of each element

  dx = sLengthOfElmX_G
  dy = sLengthOfElmY_G
  dz2 = sLengthOfElmZ2_G

  dV3 = sLengthOfElmX_G*sLengthOfElmY_G*sLengthOfElmZ2_G


!     Time savers

  sInv2rho    = 1.0_WP/(2.0_WP * frame%rho)

  ZOver2rho   = sz * sInv2rho
  salphaSq    = (2.0_WP * frame%gamma_ref * frame%rho / frame%aw)**2

  un = sqrt(und%fx**2.0_WP + und%fy**2.0_WP)


!     number of transverse nodes

  ntrans = NX_G * NY_G

!     Diff between real and imaginary nodes in the reduced system

  retim = nspinDX*nspinDY*nZ2_G

  econst = frame%aw/(frame%rho*sqrt(2.0_WP*(und%fx**2.0_WP+und%fy**2.0_WP)))

  nc = 2.0_WP*frame%aw**2/(und%fx**2.0_WP + und%fy**2.0_WP)

  nb = 2.0_WP * frame%rho / ((und%fx**2.0_WP+und%fy**2.0_WP)*frame%eta)

  qoutside=.FALSE.
  iOutside=0_IP

!  halfx = ((ReducedNX_G-1) / 2.0_WP) * sLengthOfElmX_G
!  halfy = ((ReducedNY_G-1) / 2.0_WP) * sLengthOfElmY_G

  halfx = ((nspinDX-1) / 2.0_WP) * sLengthOfElmX_G
  halfy = ((nspinDY-1) / 2.0_WP) * sLengthOfElmY_G


end subroutine rhs_tmsavers

end module rhs
