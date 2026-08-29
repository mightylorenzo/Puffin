! Copyright 2012-2018, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> Fused per-macroparticle kernels for the d/dz evaluation in getrhs.
!>
!> Every step of the getrhs chain -- p2, the mesh node index, the interpolation
!> weights, the field gather, the source scatter and the six electron equations
!> -- is element-local in the macroparticle index.  Only the scatter touches
!> another particle's data, and that is already serialised with !$OMP ATOMIC.
!> The chain therefore needs far fewer barriers than the fourteen the unfused
!> version used, but it cannot collapse to a single loop, for two reasons.
!>
!> First, the bounds check: getrhs skips the gather and the scatter entirely,
!> for every particle, if ANY particle landed outside the parallel field
!> bounds.  That is a global decision and needs a real barrier before it.
!>
!> Second, vectorisation.  !$OMP ATOMIC lowers to a compare-exchange retry
!> loop, and gfortran will not vectorise an outer loop containing one; nor will
!> it reliably vectorise a loop whose body sits under loop-invariant IFs.  A
!> single fully fused loop measured 2.5x SLOWER than the unfused chain on the
!> 1D benchmarks, because the six equations lost their SIMD -- far more than
!> the barriers were worth.  So the branchy and atomic parts are peeled off
!> into their own loops and the arithmetic loops are left branchless.
!>
!> The resulting pass 2 is three loops that cost one barrier between them:
!>
!>   source_kernel_*  atomic scatter        SCHEDULE(STATIC), NOWAIT
!>   gather_kernel_*  field -> particle     SCHEDULE(STATIC), NOWAIT
!>   equation_kernel  the six equations     SCHEDULE(STATIC), barrier
!>
!> The NOWAITs are safe.  The scatter shares nothing with the other two.  The
!> equations do read what the gather wrote, at the same index -- and OpenMP
!> guarantees that two worksharing loops in the same parallel region with the
!> same iteration count and an explicit SCHEDULE(STATIC) give each thread the
!> same iterations, so every read is of that thread's own write.  The explicit
!> SCHEDULE(STATIC) is load-bearing; do not drop it.
!>
!> All module data is taken as dummy arguments rather than by use-association.
!> That is the same reasoning as the getSource_*_kernel routines this replaces:
!> the !$OMP ATOMIC updates act as memory barriers, so a loop reading module
!> arrays directly must reload each array descriptor every iteration, putting
!> the linker's symbol placement on the critical path.  See
!> benchmark/RESULTS-2026-08-07-globals-refactor.md.
!>
!> NUMERICS: the arithmetic here is a transcription of the routines that used to
!> do this work, kept operation-for-operation and association-for-association so
!> that results stay bit-identical.  Three places where that matters and the
!> obvious tidy-up would be wrong:
!>
!>   * p_nodes uses int(), the interpolation node uses floor().  These differ
!>     for negative arguments, which is reachable for a particle outside the
!>     transverse mesh.  They are computed separately on purpose.
!>   * the gather accumulates term by term (acc = w*A + acc), not as one
!>     eight-term sum.  Floating-point addition is commutative but not
!>     associative, so a single expression would re-associate and change the
!>     last bits.
!>   * the 1D and 3D source weights associate differently -- 1D forms
!>     (chi/dV)*(1+eta*p2)/gam and then multiplies by p_perp, 3D forms
!>     (chi/dV)*(1+eta*p2)*p_perp and then divides by gam.  Preserved as-is.

module rhs_kernels

use puffin_kinds, only: WP, IP, IPL
use GlobalTypes, only: tSimulationFlags

implicit none (type, external)
private

public :: prep_kernel_1D, prep_kernel_3D, source_kernel_1D, source_kernel_3D, &
          gather_kernel_1D, gather_kernel_3D, zero_efield_kernel, equation_kernel

contains


!> 1D prep: p2, longitudinal node index, linear interpolation weights.

subroutine prep_kernel_1D(sz2, spr, spi, sgam, sp2, p_nodes, lis_GR,       &
                          seta, sgamma0, saw, dz2, fz2, bz2, nz2,          &
                          qTemporal, flags, nLocalElecs)

real(kind=wp), contiguous, intent(in) :: sz2(:), spr(:), spi(:), sgam(:)
real(kind=wp), contiguous, intent(inout) :: sp2(:)
integer(kind=ip), contiguous, intent(inout) :: p_nodes(:)
real(kind=wp), contiguous, intent(inout) :: lis_GR(:,:)
real(kind=wp), intent(in) :: seta, sgamma0, saw, dz2
integer(kind=ip), intent(in) :: fz2, bz2, nz2
logical, intent(in) :: qTemporal
type(tSimulationFlags), intent(inout) :: flags
integer(kind=ipl), intent(in) :: nLocalElecs

integer(kind=ipl) :: i
integer(kind=ip) :: z2node
real(kind=wp) :: locz2, u, rt

!$OMP DO PRIVATE(z2node, locz2, u, rt)
  do i = 1, nLocalElecs

!     p2 from gamma (see gtop2:getP2 -- written this way to avoid the
!     cancellation in 1/sqrt(1-u) - 1 as u -> 0)

      u = ( 1.0_wp + saw**2*(spr(i)**2 + spi(i)**2) ) / (sgamma0**2 * sgam(i)**2)
      rt = sqrt(1.0_wp - u)
      sp2(i) = u / (seta * rt * (1.0_wp + rt))

!     Primary mesh node for this particle (int(), deliberately not floor())

      p_nodes(i) = int(sz2(i) / dz2, kind=ip) + 1_IP - (fz2-1)

!     Surrounding nodes and interpolation weights (floor(), deliberately)

      z2node = floor(sz2(i)  / dz2)  + 1_IP
      locz2 = sz2(i) - REAL(z2node  - 1_IP, kind=wp) * dz2

      if (qTemporal) then
        if (z2node >= nz2) then
          print*, "Z2 coord is too large!! with node:", z2node, &
                  " and pos ", sz2(i)
          STOP
        end if
      end if

      if (z2node >= bz2) then
        flags%parallel_arrays_ok = .false.
      end if

      lis_GR(1,i) = (1.0_wp - locz2/dz2)
      lis_GR(2,i) = 1 - lis_GR(1,i)

  end do
!$OMP END DO

end subroutine prep_kernel_1D



!> 3D prep: p2, primary node index, trilinear interpolation weights.

subroutine prep_kernel_3D(sx, sy, sz2, spr, spi, sgam, sp2, p_nodes, lis_GR, &
                          seta, sgamma0, saw, dx, dy, dz2, halfx, halfy,     &
                          nspinDX, nspinDY, ntrndsi, fz2, bz2, nz2,          &
                          qTemporal, flags, nLocalElecs)

real(kind=wp), contiguous, intent(in) :: sx(:), sy(:), sz2(:), spr(:), spi(:), sgam(:)
real(kind=wp), contiguous, intent(inout) :: sp2(:)
integer(kind=ip), contiguous, intent(inout) :: p_nodes(:)
real(kind=wp), contiguous, intent(inout) :: lis_GR(:,:)
real(kind=wp), intent(in) :: seta, sgamma0, saw, dx, dy, dz2, halfx, halfy
integer(kind=ip), intent(in) :: nspinDX, nspinDY, ntrndsi, fz2, bz2, nz2
logical, intent(in) :: qTemporal
type(tSimulationFlags), intent(inout) :: flags
integer(kind=ipl), intent(in) :: nLocalElecs

integer(kind=ipl) :: i
integer(kind=ip) :: xnode, ynode, z2node
real(kind=wp) :: locx, locy, locz2, u, rt, &
                 x_in1, x_in2, y_in1, y_in2, z2_in1, z2_in2

!$OMP DO PRIVATE(xnode, ynode, z2node, locx, locy, locz2, u, rt, &
!$OMP x_in1, x_in2, y_in1, y_in2, z2_in1, z2_in2)
  do i = 1, nLocalElecs

!     p2 from gamma

      u = ( 1.0_wp + saw**2*(spr(i)**2 + spi(i)**2) ) / (sgamma0**2 * sgam(i)**2)
      rt = sqrt(1.0_wp - u)
      sp2(i) = u / (seta * rt * (1.0_wp + rt))

!     Primary mesh node for this particle (int(), deliberately not floor())

      p_nodes(i) = (int( (sx(i)+halfx)  / dx, kind=ip)  + 1_IP) + &
                   (int( (sy(i)+halfy)  / dy, kind=ip) * nspinDX )  + &
                   (nspinDX * nspinDY * &
                                   int(sz2(i)  / dz2, kind=ip) ) - &
                                   (fz2-1)*ntrndsi

!     Surrounding nodes and interpolation weights (floor(), deliberately)

      xnode = floor( (sx(i) + halfx ) / dx)  + 1_IP
      locx = sx(i) + halfx - real(xnode  - 1_IP, kind=wp) * dx
      x_in2 = locx / dx
      x_in1 = (1.0_wp - x_in2)

      ynode = floor( (sy(i) + halfy )  / dy)  + 1_IP
      locy = sy(i) + halfy - real(ynode  - 1_IP, kind=wp) * dy
      y_in2 = locy / dy
      y_in1 = (1.0_wp - y_in2)

      z2node = floor(sz2(i)  / dz2)  + 1_IP
      locz2 = sz2(i) - real(z2node  - 1_IP, kind=wp) * dz2
      z2_in2 = locz2 / dz2
      z2_in1 = (1.0_wp - z2_in2)

      if ((xnode >= nspinDX) .or. (xnode < 1)) then
        flags%inner_xy_ok = .false.
        flags%parallel_arrays_ok = .false.
      end if

      if ((ynode >= nspinDY) .or. (ynode < 1)) then
        flags%inner_xy_ok = .false.
        flags%parallel_arrays_ok = .false.
      end if

      if (qTemporal) then
        if (z2node >= nz2) then
          print*, "Z2 coord is too large!! with node:", z2node, &
                  " and pos ", sz2(i)
          STOP
        end if
      end if

      if (z2node >= bz2) then
        flags%parallel_arrays_ok = .false.
      end if

      lis_GR(1,i) = x_in1 * y_in1 * z2_in1
      lis_GR(2,i) = x_in2 * y_in1 * z2_in1
      lis_GR(3,i) = x_in1 * y_in2 * z2_in1
      lis_GR(4,i) = x_in2 * y_in2 * z2_in1
      lis_GR(5,i) = x_in1 * y_in1 * z2_in2
      lis_GR(6,i) = x_in2 * y_in1 * z2_in2
      lis_GR(7,i) = x_in1 * y_in2 * z2_in2
      lis_GR(8,i) = x_in2 * y_in2 * z2_in2

  end do
!$OMP END DO

end subroutine prep_kernel_3D



!> 1D source scatter.
!>
!> Kept in its own loop rather than fused into elec_kernel_1D. The !$OMP ATOMIC
!> updates lower to compare-exchange retry loops, and gfortran will not
!> vectorise an outer loop that contains them ("loop nest containing two or more
!> consecutive inner loops"). Fusing the scatter in therefore cost the six
!> electron equations their SIMD, which measured 2.5x slower on the 1D
!> benchmarks than the unfused chain -- far more than the barriers saved.
!>
!> It costs no barrier: this loop writes only sDADzr/sDADzi and reads only what
!> prep_kernel wrote, while elec_kernel writes only sField4Elec* and the six
!> derivative arrays. There is no dependency between them, so this one ends
!> NOWAIT and the barrier at the end of elec_kernel covers both.

subroutine source_kernel_1D(sDADzr, sDADzi, spr, spi, sgam, sp2,           &
                            s_chi_bar, p_nodes, lis_GR, dV3, seta,         &
                            nLocalElecs)

real(kind=wp), contiguous, intent(inout) :: sDADzr(:), sDADzi(:)
real(kind=wp), contiguous, intent(in) :: spr(:), spi(:), sgam(:), sp2(:)
real(kind=wp), contiguous, intent(in) :: s_chi_bar(:)
integer(kind=ip), contiguous, intent(in) :: p_nodes(:)
real(kind=wp), contiguous, intent(in) :: lis_GR(:,:)
real(kind=wp), intent(in) :: dV3, seta
integer(kind=ipl), intent(in) :: nLocalElecs

integer(kind=ipl) :: i
integer(kind=ip) :: pn
real(kind=wp) :: w, dadzRInst, dadzIInst

!$OMP DO SCHEDULE(STATIC) PRIVATE(pn, w, dadzRInst, dadzIInst)
  do i = 1, nLocalElecs

      pn = p_nodes(i)

      w = (s_chi_bar(i)/dV3) * (1 + seta * sp2(i) ) / sgam(i)

      dadzRInst = w * spr(i)

      !$OMP ATOMIC
      sDADzr(pn) =                                 &
        lis_GR(1,i) * dadzRInst + sDADzr(pn)

      !$OMP ATOMIC
      sDADzr(pn + 1_ip) =                          &
        lis_GR(2,i) * dadzRInst + sDADzr(pn + 1_ip)

      dadzIInst = w * spi(i)

      !$OMP ATOMIC
      sDADzi(pn) =                                 &
        lis_GR(1,i) * dadzIInst + sDADzi(pn)

      !$OMP ATOMIC
      sDADzi(pn + 1_ip) =                          &
        lis_GR(2,i) * dadzIInst + sDADzi(pn + 1_ip)

  end do
!$OMP END DO NOWAIT

end subroutine source_kernel_1D



!> 3D source scatter. Separate loop for the same reason as source_kernel_1D.

subroutine source_kernel_3D(sDADzr, sDADzi, spr, spi, sgam, sp2,           &
                            s_chi_bar, p_nodes, lis_GR, dV3, seta,         &
                            nspinDX, ntrndsi, nLocalElecs)

real(kind=wp), contiguous, intent(inout) :: sDADzr(:), sDADzi(:)
real(kind=wp), contiguous, intent(in) :: spr(:), spi(:), sgam(:), sp2(:)
real(kind=wp), contiguous, intent(in) :: s_chi_bar(:)
integer(kind=ip), contiguous, intent(in) :: p_nodes(:)
real(kind=wp), contiguous, intent(in) :: lis_GR(:,:)
real(kind=wp), intent(in) :: dV3, seta
integer(kind=ip), intent(in) :: nspinDX, ntrndsi
integer(kind=ipl), intent(in) :: nLocalElecs

integer(kind=ipl) :: i
integer(kind=ip) :: pn
real(kind=wp) :: dadzRInst, dadzIInst

!$OMP DO SCHEDULE(STATIC) PRIVATE(pn, dadzRInst, dadzIInst)
  do i = 1, nLocalElecs

      pn = p_nodes(i)

      dadzRInst = ((s_chi_bar(i)/dV3) * (1 + seta * sp2(i) ) &
                        * spr(i) / sgam(i) )

      !$OMP ATOMIC
      sDADzr(pn) =                                 &
        lis_GR(1,i) * dadzRInst + sDADzr(pn)

      !$OMP ATOMIC
      sDADzr(pn + 1_ip) =                          &
        lis_GR(2,i) * dadzRInst + sDADzr(pn + 1_ip)

      !$OMP ATOMIC
      sDADzr(pn + nspinDX) =                       &
        lis_GR(3,i) * dadzRInst + sDADzr(pn + nspinDX)

      !$OMP ATOMIC
      sDADzr(pn + nspinDX + 1_ip) =                &
        lis_GR(4,i) * dadzRInst + sDADzr(pn + nspinDX + 1_ip)

      !$OMP ATOMIC
      sDADzr(pn + ntrndsi) =                       &
        lis_GR(5,i) * dadzRInst + sDADzr(pn + ntrndsi)

      !$OMP ATOMIC
      sDADzr(pn + ntrndsi + 1_ip) =                &
        lis_GR(6,i) * dadzRInst + sDADzr(pn + ntrndsi + 1_ip)

      !$OMP ATOMIC
      sDADzr(pn + ntrndsi + nspinDX) =             &
        lis_GR(7,i) * dadzRInst + sDADzr(pn + ntrndsi + nspinDX)

      !$OMP ATOMIC
      sDADzr(pn + ntrndsi + nspinDX + 1) =         &
        lis_GR(8,i) * dadzRInst + sDADzr(pn + ntrndsi + nspinDX + 1)

      dadzIInst = ((s_chi_bar(i)/dV3) * (1 + seta * sp2(i) ) &
                        * spi(i) / sgam(i) )

      !$OMP ATOMIC
      sDADzi(pn) =                                 &
        lis_GR(1,i) * dadzIInst + sDADzi(pn)

      !$OMP ATOMIC
      sDADzi(pn + 1_ip) =                          &
        lis_GR(2,i) * dadzIInst + sDADzi(pn + 1_ip)

      !$OMP ATOMIC
      sDADzi(pn + nspinDX) =                       &
        lis_GR(3,i) * dadzIInst + sDADzi(pn + nspinDX)

      !$OMP ATOMIC
      sDADzi(pn + nspinDX + 1_ip) =                &
        lis_GR(4,i) * dadzIInst + sDADzi(pn + nspinDX + 1_ip)

      !$OMP ATOMIC
      sDADzi(pn + ntrndsi) =                       &
        lis_GR(5,i) * dadzIInst + sDADzi(pn + ntrndsi)

      !$OMP ATOMIC
      sDADzi(pn + ntrndsi + 1_ip) =                &
        lis_GR(6,i) * dadzIInst + sDADzi(pn + ntrndsi + 1_ip)

      !$OMP ATOMIC
      sDADzi(pn + ntrndsi + nspinDX) =             &
        lis_GR(7,i) * dadzIInst + sDADzi(pn + ntrndsi + nspinDX)

      !$OMP ATOMIC
      sDADzi(pn + ntrndsi + nspinDX + 1) =         &
        lis_GR(8,i) * dadzIInst + sDADzi(pn + ntrndsi + nspinDX + 1)

  end do
!$OMP END DO NOWAIT

end subroutine source_kernel_3D


!> 1D field gather: interpolate the radiation field onto each macroparticle.
!>
!> Branchless and free of atomics so it vectorises. Accumulated term by term to
!> match the old getFFelecs_1D bit for bit -- floating-point addition is
!> commutative but not associative, so folding these into one expression would
!> change the last bits.

subroutine gather_kernel_1D(sAr, sAi, sField4ElecReal, sField4ElecImag,     &
                            p_nodes, lis_GR, nLocalElecs)

real(kind=wp), contiguous, intent(in) :: sAr(:), sAi(:)
real(kind=wp), contiguous, intent(inout) :: sField4ElecReal(:), sField4ElecImag(:)
integer(kind=ip), contiguous, intent(in) :: p_nodes(:)
real(kind=wp), contiguous, intent(in) :: lis_GR(:,:)
integer(kind=ipl), intent(in) :: nLocalElecs

integer(kind=ipl) :: i
integer(kind=ip) :: pn
real(kind=wp) :: sfr, sfi

!$OMP DO SCHEDULE(STATIC) PRIVATE(pn, sfr, sfi)
  do i = 1, nLocalElecs

      pn = p_nodes(i)

      sfr = 0.0_wp
      sfr = lis_GR(1,i) * sAr(pn) + sfr
      sfr = lis_GR(2,i) * sAr(pn + 1_ip) + sfr

      sfi = 0.0_wp
      sfi = lis_GR(1,i) * sAi(pn) + sfi
      sfi = lis_GR(2,i) * sAi(pn + 1_ip) + sfi

      sField4ElecReal(i) = sfr
      sField4ElecImag(i) = sfi

  end do
!$OMP END DO NOWAIT

end subroutine gather_kernel_1D



!> 3D field gather. Trilinear, otherwise as gather_kernel_1D.

subroutine gather_kernel_3D(sAr, sAi, sField4ElecReal, sField4ElecImag,     &
                            p_nodes, lis_GR, nspinDX, ntrndsi, nLocalElecs)

real(kind=wp), contiguous, intent(in) :: sAr(:), sAi(:)
real(kind=wp), contiguous, intent(inout) :: sField4ElecReal(:), sField4ElecImag(:)
integer(kind=ip), contiguous, intent(in) :: p_nodes(:)
real(kind=wp), contiguous, intent(in) :: lis_GR(:,:)
integer(kind=ip), intent(in) :: nspinDX, ntrndsi
integer(kind=ipl), intent(in) :: nLocalElecs

integer(kind=ipl) :: i
integer(kind=ip) :: pn
real(kind=wp) :: sfr, sfi

!$OMP DO SCHEDULE(STATIC) PRIVATE(pn, sfr, sfi)
  do i = 1, nLocalElecs

      pn = p_nodes(i)

      sfr = 0.0_wp
      sfr = lis_GR(1,i) * sAr(pn) + sfr
      sfr = lis_GR(2,i) * sAr(pn + 1_ip) + sfr
      sfr = lis_GR(3,i) * sAr(pn + nspinDX) + sfr
      sfr = lis_GR(4,i) * sAr(pn + nspinDX + 1_ip) + sfr
      sfr = lis_GR(5,i) * sAr(pn + ntrndsi) + sfr
      sfr = lis_GR(6,i) * sAr(pn + ntrndsi + 1_ip) + sfr
      sfr = lis_GR(7,i) * sAr(pn + ntrndsi + nspinDX) + sfr
      sfr = lis_GR(8,i) * sAr(pn + ntrndsi + nspinDX + 1) + sfr

      sfi = 0.0_wp
      sfi = lis_GR(1,i) * sAi(pn) + sfi
      sfi = lis_GR(2,i) * sAi(pn + 1_ip) + sfi
      sfi = lis_GR(3,i) * sAi(pn + nspinDX) + sfi
      sfi = lis_GR(4,i) * sAi(pn + nspinDX + 1_ip) + sfi
      sfi = lis_GR(5,i) * sAi(pn + ntrndsi) + sfi
      sfi = lis_GR(6,i) * sAi(pn + ntrndsi + 1_ip) + sfi
      sfi = lis_GR(7,i) * sAi(pn + ntrndsi + nspinDX) + sfi
      sfi = lis_GR(8,i) * sAi(pn + ntrndsi + nspinDX + 1) + sfi

      sField4ElecReal(i) = sfr
      sField4ElecImag(i) = sfi

  end do
!$OMP END DO NOWAIT

end subroutine gather_kernel_3D



!> Zero the per-particle field. Used instead of the gather when the particles
!> are not coupled to the field, or when the parallel bounds check failed.

subroutine zero_efield_kernel(sField4ElecReal, sField4ElecImag, nLocalElecs)

real(kind=wp), contiguous, intent(inout) :: sField4ElecReal(:), sField4ElecImag(:)
integer(kind=ipl), intent(in) :: nLocalElecs

integer(kind=ipl) :: i

!$OMP DO SCHEDULE(STATIC)
  do i = 1, nLocalElecs
      sField4ElecReal(i) = 0.0_wp
      sField4ElecImag(i) = 0.0_wp
  end do
!$OMP END DO NOWAIT

end subroutine zero_efield_kernel



!> The six electron equations. Identical in 1D and 3D -- the dimensionality
!> only ever entered through the gather and the scatter.
!>
!> Branchless, no atomics, all element-local: this is the loop that has to
!> vectorise, and it is the reason the gather and the scatter are peeled off
!> into loops of their own.

subroutine equation_kernel(spr, spi, sgam, sp2,                            &
                           sField4ElecReal, sField4ElecImag,               &
                           bxu, byu, bzu,                                  &
                           sdx, sdy, sdz2, sdpr, sdpi, sdgam,              &
                           sInv2rho, seta, skappa, srho, sn2col,           &
                           nLocalElecs)

real(kind=wp), contiguous, intent(in) :: spr(:), spi(:), sgam(:), sp2(:)
real(kind=wp), contiguous, intent(in) :: sField4ElecReal(:), sField4ElecImag(:)
real(kind=wp), contiguous, intent(in) :: bxu(:), byu(:), bzu(:)
real(kind=wp), contiguous, intent(inout) :: sdx(:), sdy(:), sdz2(:), &
                                            sdpr(:), sdpi(:), sdgam(:)
real(kind=wp), intent(in) :: sInv2rho, seta, skappa, srho, sn2col
integer(kind=ipl), intent(in) :: nLocalElecs

integer(kind=ipl) :: i
real(kind=wp) :: p2, gam, pr, pri, sfr, sfi, xk

!     Loop-invariant: this is exactly the leading factor the old dxdz_f/dydz_f
!     array expressions evaluated per element, hoisted. Same association, so
!     the result is unchanged.

  xk = 2 * srho * skappa / sqrt(seta)

!$OMP DO SCHEDULE(STATIC) PRIVATE(p2, gam, pr, pri, sfr, sfi)
  do i = 1, nLocalElecs

      p2 = sp2(i)
      gam = sgam(i)
      pr = spr(i)
      pri = spi(i)
      sfr = sField4ElecReal(i)
      sfi = sField4ElecImag(i)

      sdz2(i) = p2

      sdx(i) = xk * (1 + seta * p2) / gam * pr

      sdy(i) = - xk * (1 + seta * p2) / gam * pri

      sdpr(i) = sInv2rho * ( sn2col * byu(i)                    &
                           - seta * p2 / skappa**2 * sfr )      &
                + skappa * pri / gam * (1 + seta * p2)          &
                    * sn2col * bzu(i)

      sdpi(i) = sInv2rho * ( sn2col * bxu(i)                    &
                           - seta * p2 / skappa**2 * sfi )      &
                - skappa * pr / gam * (1 + seta * p2)           &
                    * sn2col * bzu(i)

      sdgam(i) = -srho * ( 1 + seta * p2 ) / gam * 2_wp *       &
                 ( pr * sfr + pri * sfi )

  end do
!$OMP END DO

end subroutine equation_kernel

end module rhs_kernels
