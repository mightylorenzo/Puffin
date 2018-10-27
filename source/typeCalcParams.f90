! Copyright 2012-2018, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Module to store the Fortran type for calculations in Jrhs.f90

module typeCalcParams

  use paratype, only: wp, ip
  
  implicit none
  private

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Essentially just a lazy collection of vars needed throughout the source and electron
!> beam calculations in Puffin.

  type, public :: fcalcParams
  
    real(kind=wp) :: eta, aw, gamma0, kappa, rho, lw
    real(kind=wp) :: ux, uy, n2col, kxu, kyu, sZFS, sZFE, kbxSF, kbySF
    character(32_ip) :: zUndType
    integer(kind=ip) :: iOutside   ! Not *currently* used
    integer(kind=ip) :: iOutInfo

  !     Loop counters

    integer(kind=ip) :: maxEl, iNMPsLoc

  !     For interpolation

    real(kind=wp) :: halfx, halfy

  !    Shortcuts

    real(kind=wp) :: sInv2rho
    real(kind=wp) :: dV3, dx, dy, dz2

    integer(kind=ip) :: nz2, inNX, inNY, inNN
    integer(kind=ip) :: fz2, bz2, fieldMesh

    real(kind=wp), allocatable :: sp2(:), sField4ElecReal(:), &
                                  sField4ElecImag(:)

    real(kind=wp), allocatable :: bxu(:), byu(:), bzu(:)
    real(kind=wp), allocatable :: lis_GR(:,:)

  !    For index referencing

    integer(kind=ip), allocatable :: p_nodes(:)

    logical :: q1D, qSF, qParrOK, qInnerXYOK, qElectronFieldCoupling, &
               qElectronsEvolve, qFieldEvolve

  contains
    procedure :: init
    procedure :: destroy
  end type fcalcParams

contains

  subroutine init(self, dx, dy, dz2, inNX, inNY, nz2, fz2, bz2, fieldMesh, ux, uy, kxu, kyu, n2col, sZFS, sZFE, &
                  kbxSF, kbySF, zUndType, iNMPsLoc, maxEl, q1D, qSF, qParrOK, qInnerXYOK, &
                  qElectronFieldCoupling, qElectronsEvolve, qFieldEvolve, &
                  iOutInfo, tScale)

    use typeScale, only: fScale
    implicit none
    class(fcalcParams), intent(inout) :: self
    type(fScale), intent(in) :: tScale
    real(kind=wp), intent(in) :: dx, dy, dz2
    real(kind=wp), intent(in) :: ux, uy, n2col, kxu, kyu, sZFS, sZFE, kbxSF, kbySF
    character(32_ip) :: zUndType
    integer(kind=ip), intent(in) :: iNMPsLoc, maxEl, inNX, inNY, nz2, fz2, bz2, fieldMesh, iOutInfo
    logical, intent(in) :: q1D, qSF, qParrOK, qInnerXYOK, qElectronFieldCoupling, &
                           qElectronsEvolve, qFieldEvolve

    self%dx = dx
    self%dy = dy
    self%dz2 = dz2
    self%dV3 = dx*dy*dz2
    self%inNX = inNX
    self%inNY = inNY
    self%inNN = inNX * inNY
    self%nz2 = nz2
    self%fz2 = fz2
    self%bz2 = bz2
    self%fieldMesh = fieldMesh

    self%rho = tScale%rho
    self%aw = tScale%aw
    self%gamma0 = tScale%gamma0
    self%eta = tScale%eta
    self%lw = tScale%lambda_w
    self%kappa = tScale%kappa

    self%ux = ux
    self%uy = uy
    self%n2col = n2col
    self%kxu = kxu
    self%kyu = kyu
    self%sZFS = sZFS
    self%sZFE = sZFE
    self%kbxSF = kbxSF
    self%kbySF = kbySF
    self%zUndType = zUndType
    
    self%iNMPsLoc = iNMPsLoc
    self%maxEl = maxEl

    self%sInv2rho = 1.0_wp / 2.0_wp / tScale%rho

    self%halfx = ((inNX-1) / 2.0_WP) * dx
    self%halfy = ((inNY-1) / 2.0_WP) * dy
    
    self%q1D = q1D
    self%qSF = qSF
    self%qParrOK = qParrOK
    self%qInnerXYOK = qInnerXYOK
    self%qElectronFieldCoupling = qElectronFieldCoupling
    self%qElectronsEvolve = qElectronsEvolve
    self%qFieldEvolve = qFieldEvolve
    
    self%iOutInfo = iOutInfo

    allocate(self%p_nodes(iNMPsLoc))
    allocate(self%sp2(iNMPsLoc), self%sField4ElecReal(iNMPsLoc), &
             self%sField4ElecImag(iNMPsLoc))
    allocate(self%bxu(iNMPsLoc), self%byu(iNMPsLoc), self%bzu(iNMPsLoc))

    if (q1D) then
      allocate(self%lis_GR(2,iNMPsLoc))
    else
      allocate(self%lis_GR(8,iNMPsLoc))
    end if

  end subroutine init

  subroutine destroy(self)
    use typeScale, only: fScale
    implicit none
    class(fcalcParams), intent(inout) :: self

    deallocate(self%p_nodes)
    deallocate(self%sp2, self%sField4ElecReal, &
             self%sField4ElecImag)
    deallocate(self%bxu, self%byu, self%bzu)
    deallocate(self%lis_GR)

  end subroutine destroy

end module typeCalcParams


