! ###############################################
! Copyright 2012-2017, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause
! ###############################################

!> @author
!> Lawrence Campbell,
!> University of Strathclyde, 
!> Glasgow, UK
!> @brief
!> This module contains physical constants, and identifiers for array indices
!> for size arrays in Puffin.

module defs_physConstants

  use paratype
  implicit none

  real(kind=wp),    parameter :: pi  = 3.14159265358979323846
  real(kind=wp),    parameter :: c   = 2.99792458e8
  real(kind=wp),    parameter :: e_0 = 8.854187817e-12
  real(kind=wp),    parameter :: m_e = 9.1093826e-31
  real(kind=wp),    parameter :: q_e = 1.60217653e-19

  complex(kind=wp), parameter :: ci  = (0.0_WP,1.0_WP)

end module defs_physConstants



