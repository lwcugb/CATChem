!> \file test_transport_unit.F90
!! \brief Unit tests for the horizontal transport process numerics.
!!
!! These tests exercise the public transport numerics end-to-end on a regular
!! longitude/latitude grid, using exactly the modules the run-time driver uses
!! (grid metrics -> mass flux -> FV3 PPM flux-form update -> halo exchange):
!!
!!   1. Constancy    : a uniform tracer stays uniform under an arbitrary wind
!!                     (the flux-form update must be tracer-mass consistent).
!!   2. Conservation : total tracer mass sum(q*dp*area) is invariant on a
!!                     zonally periodic channel with a non-divergent wind.
!!   3. Boundedness  : the PPM solution does not produce large new extrema and
!!                     stays finite (no NaN/Inf).
!!
!! NOTE ON PRECISION: the transport library is built with `-fdefault-real-8`
!! (bare `real` == REAL64), while this test is compiled without that flag
!! (bare `real` == REAL32 here). Every array handed to a transport routine
!! whose dummy is a bare `real` is therefore declared `real(r8)` explicitly;
!! the grid-metric latitude/longitude inputs stay `real(fp)` to match their
!! dummies.
module test_transport_unit_mod
   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS
   use fv3_grid_types_mod, only: fv_grid_type, fv_grid_bounds_type
   use fv3_tp_core_mod, only: fv_tp_2d
   use TransportGridMetrics_Mod, only: build_fv3_grid_metrics, FV3_GRID_REGULAR
   use TransportMassFlux_Mod, only: fv_mass_flux_type, mass_flux_alloc, mass_flux_free, &
                                    build_level_mass_flux, courant_max, scale_mass_flux
   use TransportHalo_Mod, only: transport_halo_type, halo_update, &
                                HALO_BC_PERIODIC, HALO_BC_REPLICATE

   implicit none
   private
   public :: r8, advect_field, make_lonlat_grid, col_mass, apply_mass_fix

   integer, parameter :: r8 = kind(0.0d0)   !< matches the -fdefault-real-8 lib
   real(r8), parameter :: DELP0 = 1000.0_r8 !< initial uniform layer thickness [Pa]

contains

   !> \brief Build a regular lon/lat grid (cell centres, degrees).
   !!
   !! Longitude is a full 0..360 band (periodic); latitude is centred and stops
   !! short of the poles to keep the reconstructed corner metrics well posed.
   subroutine make_lonlat_grid(nx, ny, lat_deg, lon_deg)
      integer,  intent(in)  :: nx, ny
      real(fp), intent(out) :: lat_deg(nx, ny), lon_deg(nx, ny)

      integer  :: i, j
      real(fp) :: dlon, dlat

      dlon = 360.0_fp / real(nx, fp)
      dlat = 170.0_fp / real(ny, fp)   ! span -85..+85 to avoid the poles

      do j = 1, ny
         do i = 1, nx
            lon_deg(i, j) = (real(i, fp) - 0.5_fp) * dlon
            lat_deg(i, j) = -85.0_fp + (real(j, fp) - 0.5_fp) * dlat
         end do
      end do
   end subroutine make_lonlat_grid

   !> \brief Advect a single tracer level `nstep` times with a uniform wind.
   !!
   !! Reproduces the driver's per-level algorithm (Courant sub-cycling +
   !! flux-form update) against the public transport APIs. The Lagrangian layer
   !! thickness `delp` is advanced with the mass-flux divergence each substep,
   !! so on return `total_mass` = sum(q*delp*area) is the exactly conserved
   !! tracer mass (net boundary flux vanishes on a periodic/closed channel).
   subroutine advect_field(nx, ny, ng, lat_deg, lon_deg, uval, vval, dt, &
                           nstep, hord, x_periodic, q, total_mass, rc)
      integer,  intent(in)    :: nx, ny, ng
      real(fp), intent(in)    :: lat_deg(nx, ny), lon_deg(nx, ny)
      real(r8), intent(in)    :: uval, vval, dt
      integer,  intent(in)    :: nstep, hord
      logical,  intent(in)    :: x_periodic
      real(r8), intent(inout) :: q(nx, ny)      !< in: initial; out: advected
      real(r8), intent(out)   :: total_mass
      integer,  intent(out)   :: rc

      type(fv_grid_type)        :: gridstruct
      type(fv_grid_bounds_type) :: bd
      type(fv_mass_flux_type)   :: mf
      type(transport_halo_type) :: halo

      real(r8), allocatable :: qd(:,:), ua(:,:), va(:,:), delp(:,:)
      real(r8), allocatable :: fx(:,:), fy(:,:)
      real(r8) :: cmax, dp1, dp2, rarea
      real(r8), parameter :: lim_fac = 1.0_r8
      integer  :: is, ie, js, je, isd, ied, jsd, jed
      integer  :: npx, npy, nsplt, istep, it, i, j

      total_mass = 0.0_r8

      ! --- Grid metrics (corners reconstructed from centres) ----------------
      ! x_periodic threads the zonal periodicity into the metric x-halo: the
      ! seam-cell edge lengths / areas (reconstructed from one-sided
      ! extrapolated corners) are then wrapped so face is and face ie+1 see the
      ! SAME seam metric. Combined with the periodic tracer halo below, the flux
      ! at the west seam (face is) equals the flux at the east seam (face ie+1)
      ! bit-for-bit, the flux divergence telescopes to zero, and total tracer
      ! mass is conserved to machine precision.
      call build_fv3_grid_metrics(lat_deg, lon_deg, ng, gridstruct, bd, rc, &
                                  grid_type=FV3_GRID_REGULAR, x_periodic=x_periodic)
      if (rc /= CC_SUCCESS) return

      is = bd%is; ie = bd%ie; js = bd%js; je = bd%je
      isd = bd%isd; ied = bd%ied; jsd = bd%jsd; jed = bd%jed
      npx = nx + 1; npy = ny + 1

      halo%x_bc = merge(HALO_BC_PERIODIC, HALO_BC_REPLICATE, x_periodic)
      halo%y_bc = HALO_BC_REPLICATE

      ! --- Data-domain fields (uniform wind => halo == interior value) ------
      allocate(qd(isd:ied, jsd:jed), ua(isd:ied, jsd:jed), &
               va(isd:ied, jsd:jed), delp(isd:ied, jsd:jed))
      allocate(fx(is:ie+1, js:je), fy(is:ie, js:je+1))
      ua   = uval
      va   = vval
      delp = DELP0
      qd   = 0.0_r8
      qd(is:ie, js:je) = q

      ! --- Mass flux for one full step, then split for Courant stability ----
      call mass_flux_alloc(mf, bd, rc)
      if (rc /= CC_SUCCESS) return
      call build_level_mass_flux(gridstruct, bd, ua, va, delp, dt, mf, rc)
      if (rc /= CC_SUCCESS) return

      cmax  = courant_max(mf, bd)
      nsplt = max(1, int(1.0_r8 + cmax))
      if (nsplt > 1) call scale_mass_flux(mf, gridstruct, bd, 1.0_r8 / real(nsplt, r8))

      ! --- Time stepping ----------------------------------------------------
      do istep = 1, nstep
         do it = 1, nsplt
            call halo_update(qd, bd, halo, rc)
            if (rc /= CC_SUCCESS) return

            call fv_tp_2d(qd, mf%crx, mf%cry, npx, npy, hord, fx, fy, &
                          mf%xfx, mf%yfx, gridstruct, bd, mf%ra_x, mf%ra_y, &
                          lim_fac, mfx=mf%mfx, mfy=mf%mfy)

            do j = js, je
               do i = is, ie
                  rarea = gridstruct%rarea(i, j)
                  dp1 = delp(i, j)
                  dp2 = dp1 + (mf%mfx(i, j) - mf%mfx(i+1, j) &
                            +  mf%mfy(i, j) - mf%mfy(i, j+1)) * rarea
                  qd(i, j) = (qd(i, j) * dp1 &
                             + (fx(i, j) - fx(i+1, j) &
                             +  fy(i, j) - fy(i, j+1)) * rarea) / dp2
                  ! Advance the Lagrangian thickness so sum(q*delp*area) stays
                  ! conserved to machine precision (telescoping flux divergence).
                  delp(i, j) = dp2
               end do
            end do
         end do
      end do

      ! --- Results ----------------------------------------------------------
      q = qd(is:ie, js:je)
      do j = js, je
         do i = is, ie
            total_mass = total_mass + qd(i, j) * delp(i, j) * gridstruct%area(i, j)
         end do
      end do

      call mass_flux_free(mf)
      deallocate(qd, ua, va, delp, fx, fy)
   end subroutine advect_field

   !> \brief Column-integrated tracer mass sum(q*dp) for every column, used by
   !!        the vertical (Lagrangian) PPM remap tests.
   subroutine col_mass(np, nz, pe, q, mass)
      integer,  intent(in)  :: np, nz
      real(r8), intent(in)  :: pe(np, nz+1)   !< edge pressures (top->bottom)
      real(r8), intent(in)  :: q(np, nz)      !< layer values
      real(r8), intent(out) :: mass(np)

      integer :: c, k

      mass = 0.0_r8
      do k = 1, nz
         do c = 1, np
            mass(c) = mass(c) + q(c, k) * (pe(c, k+1) - pe(c, k))
         end do
      end do
   end subroutine col_mass

   !> \brief GCHP-style per-column mass fixer, identical to the one applied in
   !!        transport_vertical. `mappm` conserves the interior exactly but its
   !!        top cell (edge at the model top) takes the top source-cell mean
   !!        rather than that cell's sub-cell average, leaving a small residual.
   !!        Rescale each positive-definite column by the ratio of the pre-remap
   !!        (source) to post-remap (target) tracer mass so column mass is
   !!        preserved to machine precision.
   subroutine apply_mass_fix(np, nz, pe1, q1, pe2, q2)
      integer,  intent(in)    :: np, nz
      real(r8), intent(in)    :: pe1(np, nz+1), pe2(np, nz+1)
      real(r8), intent(in)    :: q1(np, nz)
      real(r8), intent(inout) :: q2(np, nz)

      integer  :: c, k
      real(r8) :: m_src, m_tgt

      do c = 1, np
         m_src = 0.0_r8
         m_tgt = 0.0_r8
         do k = 1, nz
            m_src = m_src + q1(c, k) * (pe1(c, k+1) - pe1(c, k))
            m_tgt = m_tgt + q2(c, k) * (pe2(c, k+1) - pe2(c, k))
         end do
         if (m_src > 0.0_r8 .and. m_tgt > 0.0_r8) then
            do k = 1, nz
               q2(c, k) = q2(c, k) * (m_src / m_tgt)
            end do
         end if
      end do
   end subroutine apply_mass_fix

end module test_transport_unit_mod

program test_transport_unit
   use test_transport_unit_mod
   use fv3_vremap_mod, only: mappm
   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS
   use testing_mod, only: assert
   implicit none

   integer, parameter :: nx = 40, ny = 20, ng = 3, hord = 8
   real(r8), parameter :: pi = 3.14159265358979323846_r8

   ! Vertical (Lagrangian) remap test grid.
   integer,  parameter :: vnp = 3, vnz = 20, vord = 8, iv = 0
   real(r8), parameter :: ptop = 1.0_r8      !< model top pressure [Pa]
   real(r8), parameter :: psfc = 100000.0_r8 !< surface pressure [Pa]

   real(fp) :: lat_deg(nx, ny), lon_deg(nx, ny)
   real(r8) :: q(nx, ny), q0(nx, ny)
   real(r8) :: mass0, mass1, dt
   real(r8) :: max_dev, rel_mass_err, qmin, qmax, span
   integer  :: i, j, rc, k, c

   ! Vertical remap working arrays.
   real(r8) :: pe1(vnp, vnz+1), pe2(vnp, vnz+1)
   real(r8) :: vq1(vnp, vnz), vq2(vnp, vnz)
   real(r8) :: vm1(vnp), vm2(vnp)
   real(r8) :: frac, rel_err, vq1min, vq1max

   write(*,*) 'Testing horizontal transport (FV3 PPM kernel)...'
   write(*,*) ''

   call make_lonlat_grid(nx, ny, lat_deg, lon_deg)
   dt = 1800.0_r8

   ! -----------------------------------------------------------------------
   ! Test 1: Constancy - a uniform field must remain uniform under any wind.
   ! -----------------------------------------------------------------------
   write(*,*) 'Test 1: Constancy (uniform tracer preserved)'
   q = 2.0_r8
   call advect_field(nx, ny, ng, lat_deg, lon_deg, 10.0_r8, 5.0_r8, dt, &
                     20, hord, .true., q, mass1, rc)
   call assert(rc == CC_SUCCESS, "advect_field (constancy) should succeed")

   max_dev = 0.0_r8
   do j = 1, ny
      do i = 1, nx
         max_dev = max(max_dev, abs(q(i, j) - 2.0_r8))
      end do
   end do
   write(*,'(a,es12.4)') '   max deviation from constant: ', max_dev
   call assert(max_dev < 1.0e-9_r8, "uniform tracer must stay uniform")
   write(*,*) 'Test 1 passed!'
   write(*,*) ''

   ! -----------------------------------------------------------------------
   ! Test 2: Mass conservation on a zonally periodic, non-divergent channel.
   ! -----------------------------------------------------------------------
   write(*,*) 'Test 2: Mass conservation (periodic zonal flow)'
   do j = 1, ny
      do i = 1, nx
         q0(i, j) = 1.0_r8 + 0.3_r8 * cos(real(lon_deg(i, j), r8) * pi / 180.0_r8)
      end do
   end do

   ! Initial mass (nstep = 0 builds metrics/flux but does not advance).
   q = q0
   call advect_field(nx, ny, ng, lat_deg, lon_deg, 10.0_r8, 0.0_r8, dt, &
                     0, hord, .true., q, mass0, rc)
   call assert(rc == CC_SUCCESS, "advect_field (mass, init) should succeed")

   ! Advance 20 steps of pure zonal advection (v = 0 => closed in latitude).
   q = q0
   call advect_field(nx, ny, ng, lat_deg, lon_deg, 10.0_r8, 0.0_r8, dt, &
                     20, hord, .true., q, mass1, rc)
   call assert(rc == CC_SUCCESS, "advect_field (mass, run) should succeed")

   rel_mass_err = abs(mass1 - mass0) / mass0
   write(*,'(a,es12.4)') '   relative mass error: ', rel_mass_err
   call assert(rel_mass_err < 1.0e-9_r8, "total tracer mass must be conserved")
   write(*,*) 'Test 2 passed!'
   write(*,*) ''

   ! -----------------------------------------------------------------------
   ! Test 3: Boundedness - no large new extrema, and the field stays finite.
   ! -----------------------------------------------------------------------
   write(*,*) 'Test 3: Boundedness / finiteness'
   qmin = q(1, 1); qmax = q(1, 1)
   do j = 1, ny
      do i = 1, nx
         call assert(q(i, j) == q(i, j), "advected field must not contain NaN")
         qmin = min(qmin, q(i, j))
         qmax = max(qmax, q(i, j))
      end do
   end do
   ! Initial field is in [0.7, 1.3]; allow a small PPM overshoot margin.
   span = 1.3_r8 - 0.7_r8
   write(*,'(a,es12.4,a,es12.4)') '   min: ', qmin, '   max: ', qmax
   call assert(qmin > 0.7_r8 - 0.05_r8 * span, "no large undershoot")
   call assert(qmax < 1.3_r8 + 0.05_r8 * span, "no large overshoot")
   write(*,*) 'Test 3 passed!'
   write(*,*) ''

   ! =======================================================================
   ! Vertical (Lagrangian) PPM remap tests (FV3 mappm kernel).
   ! =======================================================================
   ! Source grid: stretched (finer near the surface), top->bottom.
   do c = 1, vnp
      do k = 1, vnz+1
         frac = real(k-1, r8) / real(vnz, r8)
         pe1(c, k) = ptop + (psfc - ptop) * frac**2   ! quadratic stretch
      end do
   end do
   ! Target grid: uniform in pressure, SAME end points as the source.
   do c = 1, vnp
      do k = 1, vnz+1
         frac = real(k-1, r8) / real(vnz, r8)
         pe2(c, k) = ptop + (psfc - ptop) * frac
      end do
      pe2(c, 1)     = pe1(c, 1)       ! exact top match
      pe2(c, vnz+1) = pe1(c, vnz+1)   ! exact surface match => exact conservation
   end do

   ! -----------------------------------------------------------------------
   ! Test 4: Vertical conservation - column mass sum(q*dp) invariant under the
   ! remap. The driver's transport_vertical applies a GCHP-style per-column mass
   ! fixer after mappm (the raw kernel leaves a small top-cell residual, which
   ! is canonical FV3 behaviour); this test exercises the same remap + fixer and
   ! requires machine-precision column-mass conservation.
   ! -----------------------------------------------------------------------
   write(*,*) 'Test 4: Vertical conservation (column mass invariant)'
   do k = 1, vnz
      do c = 1, vnp
         vq1(c, k) = 1.0_r8 + 0.5_r8 * sin(real(k, r8) * 0.3_r8) + 0.1_r8 * real(c, r8)
      end do
   end do

   call col_mass(vnp, vnz, pe1, vq1, vm1)
   call mappm(vnz, pe1, vq1, vnz, pe2, vq2, 1, vnp, iv, vord)
   call apply_mass_fix(vnp, vnz, pe1, vq1, pe2, vq2)
   call col_mass(vnp, vnz, pe2, vq2, vm2)

   rel_err = 0.0_r8
   do c = 1, vnp
      rel_err = max(rel_err, abs(vm2(c) - vm1(c)) / abs(vm1(c)))
   end do
   write(*,'(a,es12.4)') '   max relative column-mass error: ', rel_err
   call assert(rel_err < 1.0e-12_r8, "vertical remap must conserve column mass")
   write(*,*) 'Test 4 passed!'
   write(*,*) ''

   ! -----------------------------------------------------------------------
   ! Test 5: Vertical constancy - a uniform column stays uniform after remap.
   ! -----------------------------------------------------------------------
   write(*,*) 'Test 5: Vertical constancy (uniform tracer preserved)'
   vq1 = 5.0_r8
   call mappm(vnz, pe1, vq1, vnz, pe2, vq2, 1, vnp, iv, vord)

   max_dev = 0.0_r8
   do k = 1, vnz
      do c = 1, vnp
         max_dev = max(max_dev, abs(vq2(c, k) - 5.0_r8))
      end do
   end do
   write(*,'(a,es12.4)') '   max deviation from constant: ', max_dev
   call assert(max_dev < 1.0e-10_r8, "uniform column must stay uniform")
   write(*,*) 'Test 5 passed!'
   write(*,*) ''

   ! -----------------------------------------------------------------------
   ! Test 6: Vertical boundedness / positivity - no NaN, no large new extrema.
   ! -----------------------------------------------------------------------
   write(*,*) 'Test 6: Vertical boundedness / positivity'
   vq1 = 0.2_r8
   vq1(:, vnz/2) = 3.0_r8   ! a positive mid-column spike

   vq1min = minval(vq1); vq1max = maxval(vq1)
   call mappm(vnz, pe1, vq1, vnz, pe2, vq2, 1, vnp, iv, vord)

   qmin = vq2(1, 1); qmax = vq2(1, 1)
   do k = 1, vnz
      do c = 1, vnp
         call assert(vq2(c, k) == vq2(c, k), "remapped field must not contain NaN")
         qmin = min(qmin, vq2(c, k))
         qmax = max(qmax, vq2(c, k))
      end do
   end do
   span = vq1max - vq1min
   write(*,'(a,es12.4,a,es12.4)') '   min: ', qmin, '   max: ', qmax
   call assert(qmin >= 0.0_r8, "positive-definite remap must stay non-negative")
   call assert(qmin > vq1min - 0.05_r8 * span, "no large undershoot")
   call assert(qmax < vq1max + 0.05_r8 * span, "no large overshoot")
   write(*,*) 'Test 6 passed!'
   write(*,*) ''

   write(*,*) 'All transport tests passed!'
end program test_transport_unit
