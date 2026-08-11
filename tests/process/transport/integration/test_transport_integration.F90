!> \file test_transport_integration.F90
!! \brief End-to-end integration test for the FV3 transport process.
!!
!! Unlike the transport UNIT test (which drives the vendored kernels directly on
!! a hand-built grid), this test connects the process to the full model: it
!! builds a CATChemCore instance from a configuration file, registers the
!! transport process, sets up a global periodic lon-lat grid with winds and
!! surface pressure, adds the process to the run pipeline, and advances the
!! model with `core%run_timestep`. That path exercises the seam the unit test
!! deliberately bypasses: config parsing of the `transport:` block, the
!! ProcessManager dispatching a full-field (non-column) process, grid-metric
!! construction from the real GridManager geometry, and state plumbing through
!! the shared container.
!!
!! Two scenarios are run through the same core:
!!
!!   Phase A - pressure fixer OFF (no PS_NEXT supplied):
!!             a non-uniform (zonal-cosine) tracer advected by a uniform zonal
!!             wind on a periodic grid. The wind is non-divergent, so the layer
!!             thickness is unchanged (the vertical Lagrangian remap is a no-op)
!!             and total tracer mass sum(q*dp*area) must be conserved; the
!!             monotone PPM solution must stay bounded and finite.
!!
!!   Phase B - vertical reconciliation onto PS_NEXT (PS_NEXT supplied and marked
!!             populated): a UNIFORM tracer with a latitude-structured end-of-step
!!             surface pressure target. The vertical Lagrangian remap maps the
!!             advected column onto the PS_NEXT grid and the unconditional
!!             per-column mass rescale conserves each column's tracer mass. A
!!             uniform mixing ratio therefore does NOT stay globally uniform:
!!             conserving mass while the target air-column mass varies with
!!             latitude forces the mixing ratio to vary with latitude. Each
!!             column stays vertically uniform at the mass-conserving value and
!!             finite (see check_phase_b for the exact expected value).
!!
!! The run reuses the shared standalone configuration
!! (tests/Configs/Default/CATChem_new_config_standalone.yml), which already
!! declares a global periodic 72-level transport block (x_periodic: true,
!! horizontal + vertical, hord/vord = 8). The grid dimensions used here are
!! supplied programmatically via `with_grid`, overriding the file, and only the
!! transport process is added to the pipeline, so the emission/chemistry blocks
!! in that config are never exercised.
!!
!! The numeric accuracy of the vertical remap itself is validated to machine
!! precision by the unit test (Tests 4-6); here we verify the end-to-end wiring.
!!
!! PRECISION: this test is compiled WITHOUT `-fdefault-real-8`, so bare `real`
!! is REAL32 and `real(fp)` matches the core state (fp = kind(0.0)). It only
!! calls core APIs (all `real(fp)`) and the argument-free transport registration
!! routine, so there is no bare-`real` seam with the -r8 transport library.
!! Tracer mass is accumulated in double precision (kind=8) to avoid REAL32
!! summation error over the grid.
program test_transport_integration
   use precision_mod, only: fp
   use iso_fortran_env, only: output_unit, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use error_mod, only: CC_SUCCESS
   use CATChemCore_Mod, only: CATChemCoreType, CATChemBuilderType
   use StateManager_Mod, only: StateManagerType
   use ProcessManager_Mod, only: ProcessManagerType
   use GridManager_Mod, only: GridManagerType
   use MetState_Mod, only: MetStateType
   use ChemState_Mod, only: ChemStateType
   use TimeState_Mod, only: TimeStateType
   use met_utilities_mod, only: get_hybrid_ab
   use TransportProcessCreator_Mod, only: register_transport_process

   implicit none

   integer,  parameter :: dp = kind(0.0d0)     !< mass accumulation precision
   integer,  parameter :: nx = 24, ny = 12, nz = 72
   integer,  parameter :: nsteps = 3
   real(fp), parameter :: dt = 600.0_fp        !< transport timestep [s]
   real(fp), parameter :: PS0 = 1.0e5_fp       !< reference surface pressure [Pa]
   real(fp), parameter :: U0  = 10.0_fp        !< uniform zonal wind [m s-1]
   real(fp), parameter :: U0C = 40.0_fp        !< Phase C baroclinic wind amplitude [m s-1]
   real(fp), parameter :: Q0  = 1.0e-6_fp      !< reference tracer mixing ratio
   real(fp), parameter :: TOL = 1.0e-4_fp      !< REAL32-appropriate rel. tol.
   real(fp), parameter :: PI  = 3.14159265358979323846_fp
   real(fp), parameter :: REARTH = 6.371e6_fp  !< Earth radius [m]
   character(len=*), parameter :: config_file = './CATChem_new_config_standalone.yml'

   type(CATChemCoreType)             :: core
   type(CATChemBuilderType)          :: builder
   type(ProcessManagerType), pointer :: process_mgr
   type(StateManagerType),   pointer :: state_mgr
   type(MetStateType),       pointer :: met
   type(ChemStateType),      pointer :: chem
   type(TimeStateType),      pointer :: time_state

   real(fp) :: area(nx, ny), latc(nx, ny), lonc(nx, ny)
   real(dp) :: m0_phaseA
   real(dp) :: m0_phaseC
   logical  :: all_ok
   integer  :: rc, it

   all_ok = .true.

   write(output_unit,'(A)') '===================================='
   write(output_unit,'(A)') '=== TRANSPORT INTEGRATION TESTS  ==='
   write(output_unit,'(A)') '===================================='
   write(output_unit,'(A)') 'Driving the FV3 transport process end-to-end through'
   write(output_unit,'(A)') 'CATChemCore (config -> register -> run loop).'
   write(output_unit,'(A)') ''

   ! ----------------------------------------------------------------------
   ! Step 1: build the core (config, managers, state) on a nx*ny*nz grid.
   ! ----------------------------------------------------------------------
   write(output_unit,'(A)') 'Step 1: Initializing CATChem Core...'
   call builder%init()
   builder = builder%with_name('TransportIntegrationTest')
   builder = builder%with_config(config_file)
   builder = builder%with_grid(nx, ny, nz)
   call builder%build(core, rc)
   if (rc /= CC_SUCCESS) then
      write(error_unit,'(A)') 'ERROR: CATChemCore initialization failed'
      all_ok = .false.
      go to 999
   end if
   write(output_unit,'(A,I0,A,I0,A,I0,A)') '  ok CATChemCore initialized: ', &
      nx, ' x ', ny, ' x ', nz, ' (lon x lat x lev)'

   ! Resolve state handles.
   state_mgr  => core%get_state_manager()
   met        => state_mgr%get_met_state_ptr()
   chem       => state_mgr%get_chem_state_ptr()
   time_state => state_mgr%get_time_state_ptr()
   if (.not. (associated(state_mgr) .and. associated(met) .and. &
      associated(chem) .and. associated(time_state))) then
      write(error_unit,'(A)') 'ERROR: could not resolve state handles'
      all_ok = .false.
      go to 999
   end if
   if (chem%nSpeciesAdvect <= 0) then
      write(error_unit,'(A)') 'ERROR: no advected species (check species config)'
      all_ok = .false.
      go to 999
   end if
   write(output_unit,'(A,I0)') '  ok Advected species: ', chem%nSpeciesAdvect

   ! Transport reads its timestep from the shared time state.
   time_state%timestep = dt

   ! ----------------------------------------------------------------------
   ! Step 2: geometry (global periodic lon-lat grid) + register/add process.
   ! ----------------------------------------------------------------------
   write(output_unit,'(A)') ''
   write(output_unit,'(A)') 'Step 2: Setting up geometry and registering transport...'
   call setup_geometry()

   process_mgr => core%get_process_manager()
   call register_transport_process(process_mgr, rc)
   if (rc /= CC_SUCCESS) then
      write(error_unit,'(A)') 'ERROR: failed to register transport process'
      all_ok = .false.
      go to 999
   end if
   call core%add_process('transport', rc)
   if (rc /= CC_SUCCESS) then
      write(error_unit,'(A)') 'ERROR: failed to add transport process'
      all_ok = .false.
      go to 999
   end if
   write(output_unit,'(A)') '  ok transport registered and added to the run pipeline'

   ! ----------------------------------------------------------------------
   ! Phase A: pressure fixer OFF (no PS_NEXT). Non-uniform tracer, uniform
   !          zonal (non-divergent) wind -> mass conservation + boundedness.
   ! ----------------------------------------------------------------------
   write(output_unit,'(A)') ''
   write(output_unit,'(A)') 'Phase A: transport WITHOUT pressure fixer...'
   call setup_winds()
   call set_tracer_cosine()          ! non-uniform zonal tracer in [0.5, 1.5]*Q0
   m0_phaseA = total_tracer_mass()   ! reference mass BEFORE advection
   ! (PS_NEXT intentionally left unmarked -> fixer stays off.)
   do it = 1, nsteps
      call core%run_timestep(it, dt, rc)
      if (rc /= CC_SUCCESS) then
         write(error_unit,'(A,I0)') 'ERROR: Phase A timestep failed at step ', it
         all_ok = .false.
         go to 999
      end if
   end do
   call check_phase_a()

   ! ----------------------------------------------------------------------
   ! Phase C: VERTICAL transport. A surface-loaded tracer under a baroclinic,
   !          purely-zonal convergent wind (V = 0; U strong near the surface and
   !          vanishing aloft) on the periodic grid. Low-level horizontal mass
   !          convergence thickens the surface Lagrangian layers, and the
   !          vertical Lagrangian remap must then lift tracer mass off the
   !          surface into the layers above. Because the wind is purely zonal on
   !          a periodic grid (V = 0), there is no meridional boundary flux and
   !          global mass is conserved exactly, so we can assert BOTH:
   !            (1) total 3-D mass conservation through horizontal + vertical;
   !            (2) above-surface mass becomes non-zero -- direct proof the
   !                vertical remap redistributes mass, since horizontal
   !                advection alone can never move mass between levels.
   !          The pressure fixer stays OFF here (PS_NEXT is not set until Phase
   !          B, which runs next), so this is pure advective transport.
   ! ----------------------------------------------------------------------
   write(output_unit,'(A)') ''
   write(output_unit,'(A)') 'Phase C: VERTICAL transport (surface tracer + convergent wind)...'
   call setup_winds_baroclinic()
   call set_tracer_surface()
   m0_phaseC = total_tracer_mass()   ! reference 3-D mass BEFORE transport
   do it = 1, nsteps
      call core%run_timestep(it, dt, rc)
      if (rc /= CC_SUCCESS) then
         write(error_unit,'(A,I0)') 'ERROR: Phase C timestep failed at step ', it
         all_ok = .false.
         go to 999
      end if
   end do
   call check_phase_c()

   ! ----------------------------------------------------------------------
   ! Phase B: pressure fixer ON. Uniform tracer, uniform zonal wind, and a
   !          latitude-structured PS_NEXT target -> constancy + finiteness.
   ! ----------------------------------------------------------------------
   write(output_unit,'(A)') ''
   write(output_unit,'(A)') 'Phase B: transport WITH pressure fixer...'
   call setup_winds()
   call set_tracer_uniform()
   call setup_ps_next()              ! allocate + populate + mark PS_NEXT
   do it = 1, nsteps
      call core%run_timestep(it, dt, rc)
      if (rc /= CC_SUCCESS) then
         write(error_unit,'(A,I0)') 'ERROR: Phase B timestep failed at step ', it
         all_ok = .false.
         go to 999
      end if
   end do
   call check_phase_b()

   ! ----------------------------------------------------------------------
   call core%finalize(rc)
   if (rc /= CC_SUCCESS) write(error_unit,'(A)') 'WARNING: core finalize had issues'

999 continue
   write(output_unit,'(A)') ''
   write(output_unit,'(A)') '===================================='
   if (all_ok) then
      write(output_unit,'(A)') '=== ALL TRANSPORT TESTS PASSED!  ==='
   else
      write(output_unit,'(A)') '=== SOME TRANSPORT TESTS FAILED  ==='
   end if
   write(output_unit,'(A)') '===================================='
   if (.not. all_ok) stop 1

contains

   !> Build a global periodic lon-lat grid (cell centres, degrees) and the true
   !! spherical cell areas, and store them into MetState (LAT/LON/AREA_M2).
   !! Longitude spans a full 0..360 band (periodic); latitude is centred and
   !! stops short of the poles, matching the unit test's make_lonlat_grid.
   subroutine setup_geometry()
      integer  :: i, j
      real(fp) :: dlon, dlat, lat_s, lat_n
      dlon = 360.0_fp / real(nx, fp)
      dlat = 170.0_fp / real(ny, fp)          ! span -85..+85
      do j = 1, ny
         do i = 1, nx
            lonc(i, j) = (real(i, fp) - 0.5_fp) * dlon
            latc(i, j) = -85.0_fp + (real(j, fp) - 0.5_fp) * dlat
            lat_s = (latc(i, j) - 0.5_fp * dlat) * PI / 180.0_fp
            lat_n = (latc(i, j) + 0.5_fp * dlat) * PI / 180.0_fp
            area(i, j) = REARTH * REARTH * (dlon * PI / 180.0_fp) * &
               (sin(lat_n) - sin(lat_s))
         end do
      end do
      met%LON(:, :)     = lonc(:, :)
      met%LAT(:, :)     = latc(:, :)
      met%AREA_M2(:, :) = area(:, :)
   end subroutine setup_geometry

   !> Uniform zonal wind (U0, V=0) and a hybrid-sigma layer-thickness profile
   !! consistent with the pressure fixer's Ap/Bp table. A constant zonal wind on
   !! a periodic grid is non-divergent, so the layer thickness is preserved.
   subroutine setup_winds()
      integer  :: i, j, k
      real(fp), allocatable :: ap(:), bp(:)
      real(fp) :: dpk
      logical  :: ok
      call get_hybrid_ab(nz, ap, bp, ok)
      met%PS(:, :) = PS0
      do k = 1, nz
         if (ok) then
            dpk = (ap(k) - ap(k+1)) + (bp(k) - bp(k+1)) * PS0
         else
            dpk = PS0 / real(nz, fp)   ! fallback (should not happen for nz=72)
         end if
         do j = 1, ny
            do i = 1, nx
               met%U(i, j, k)    = U0
               met%V(i, j, k)    = 0.0_fp
               met%DELP(i, j, k) = dpk
            end do
         end do
      end do
   end subroutine setup_winds

   !> Latitude-structured end-of-step surface pressure so the vertical remap
   !! reconciles onto a non-trivial PS_NEXT, then mark it populated so the
   !! transport gate (is_field_set) engages the reconciliation.
   subroutine setup_ps_next()
      integer  :: i, j
      real(fp) :: latr
      if (.not. allocated(met%PS_NEXT)) allocate(met%PS_NEXT(nx, ny))
      do j = 1, ny
         do i = 1, nx
            latr = latc(i, j) * PI / 180.0_fp
            met%PS_NEXT(i, j) = PS0 * (1.0_fp + 0.01_fp * sin(latr))
         end do
      end do
      call met%mark_field_set('PS_NEXT')
   end subroutine setup_ps_next

   !> Non-uniform (zonal-cosine) initial tracer, bounded in [0.5, 1.5]*Q0.
   subroutine set_tracer_cosine()
      integer  :: s, isp, i, j, k
      real(fp) :: val
      do s = 1, chem%nSpeciesAdvect
         isp = chem%AdvectIndex(s)
         if (.not. associated(chem%ChemSpecies(isp)%conc)) cycle
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  val = Q0 * (1.0_fp + 0.5_fp * cos(lonc(i, j) * PI / 180.0_fp))
                  chem%ChemSpecies(isp)%conc(i, j, k) = val
               end do
            end do
         end do
      end do
   end subroutine set_tracer_cosine

   !> Uniform initial tracer (== Q0 everywhere).
   subroutine set_tracer_uniform()
      integer :: s, isp
      do s = 1, chem%nSpeciesAdvect
         isp = chem%AdvectIndex(s)
         if (.not. associated(chem%ChemSpecies(isp)%conc)) cycle
         chem%ChemSpecies(isp)%conc(:, :, :) = Q0
      end do
   end subroutine set_tracer_uniform

   !> Baroclinic, purely-zonal convergent wind used to drive VERTICAL transport.
   !! V = 0 everywhere (no meridional boundary flux => exact global mass
   !! conservation on the periodic grid). U = U0C*cos(lon)*wprof(k), with
   !! wprof = 1 through the lowest third of the (surface-first) column, tapering
   !! linearly to 0 by two-thirds and 0 above. dU/dx = -U0C*sin(lon)*wprof(k)
   !! gives low-level horizontal convergence/divergence that vanishes with
   !! height, so the surface Lagrangian layers deform while the upper layers do
   !! not -- exactly the structure that forces genuine vertical (Lagrangian
   !! remap) transport of a surface-confined tracer.
   subroutine setup_winds_baroclinic()
      integer  :: i, j, k, kbot
      real(fp), allocatable :: ap(:), bp(:)
      real(fp) :: dpk, wprof
      logical  :: ok
      call get_hybrid_ab(nz, ap, bp, ok)
      met%PS(:, :) = PS0
      kbot = max(1, nz / 3)          ! lowest third of the surface-first column
      do k = 1, nz
         if (ok) then
            dpk = (ap(k) - ap(k+1)) + (bp(k) - bp(k+1)) * PS0
         else
            dpk = PS0 / real(nz, fp)
         end if
         if (k <= kbot) then
            wprof = 1.0_fp
         else if (k < 2 * kbot) then
            wprof = real(2 * kbot - k, fp) / real(kbot, fp)
         else
            wprof = 0.0_fp
         end if
         do j = 1, ny
            do i = 1, nx
               met%U(i, j, k)    = U0C * cos(lonc(i, j) * PI / 180.0_fp) * wprof
               met%V(i, j, k)    = 0.0_fp
               met%DELP(i, j, k) = dpk
            end do
         end do
      end do
   end subroutine setup_winds_baroclinic

   !> Surface-loaded initial tracer: Q0 in the surface layer (k=1), zero above.
   !! Mirrors a species with surface-only emissions, so any tracer that appears
   !! at k>=2 after the run must have been transported there vertically.
   subroutine set_tracer_surface()
      integer :: s, isp
      do s = 1, chem%nSpeciesAdvect
         isp = chem%AdvectIndex(s)
         if (.not. associated(chem%ChemSpecies(isp)%conc)) cycle
         chem%ChemSpecies(isp)%conc(:, :, :) = 0.0_fp
         chem%ChemSpecies(isp)%conc(:, :, 1) = Q0
      end do
   end subroutine set_tracer_surface

   !> Total advected-tracer mass sum(conc*DELP*area), summed in double precision.
   function total_tracer_mass() result(m)
      real(dp) :: m
      integer  :: s, isp, i, j, k
      m = 0.0_dp
      do s = 1, chem%nSpeciesAdvect
         isp = chem%AdvectIndex(s)
         if (.not. associated(chem%ChemSpecies(isp)%conc)) cycle
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  m = m + real(chem%ChemSpecies(isp)%conc(i, j, k), dp) * &
                     real(met%DELP(i, j, k), dp) * real(area(i, j), dp)
               end do
            end do
         end do
      end do
   end function total_tracer_mass

   !> Advected-tracer mass held ABOVE the surface layer (levels k>=2), summed in
   !! double precision. It is exactly zero unless the vertical remap has moved
   !! mass off the surface, so it is a sharp on/off signal for vertical transport.
   function above_surface_mass() result(m)
      real(dp) :: m
      integer  :: s, isp, i, j, k
      m = 0.0_dp
      do s = 1, chem%nSpeciesAdvect
         isp = chem%AdvectIndex(s)
         if (.not. associated(chem%ChemSpecies(isp)%conc)) cycle
         do k = 2, nz
            do j = 1, ny
               do i = 1, nx
                  m = m + real(chem%ChemSpecies(isp)%conc(i, j, k), dp) * &
                     real(met%DELP(i, j, k), dp) * real(area(i, j), dp)
               end do
            end do
         end do
      end do
   end function above_surface_mass

   !> .true. if every advected tracer value is finite.
   logical function tracer_is_finite() result(ok)
      integer :: s, isp, i, j, k
      ok = .true.
      do s = 1, chem%nSpeciesAdvect
         isp = chem%AdvectIndex(s)
         if (.not. associated(chem%ChemSpecies(isp)%conc)) cycle
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  if (.not. ieee_is_finite(chem%ChemSpecies(isp)%conc(i, j, k))) then
                     ok = .false.
                     return
                  end if
               end do
            end do
         end do
      end do
   end function tracer_is_finite

   !> Phase A checks: mass conservation (non-divergent wind), boundedness of the
   !! monotone PPM solution, and finiteness.
   subroutine check_phase_a()
      real(dp) :: m1, rel
      real(fp) :: qmin, qmax
      integer  :: s, isp
      logical  :: ok

      m1 = total_tracer_mass()          ! live post-run advected mass

      qmin = huge(1.0_fp); qmax = -huge(1.0_fp)
      do s = 1, chem%nSpeciesAdvect
         isp = chem%AdvectIndex(s)
         if (.not. associated(chem%ChemSpecies(isp)%conc)) cycle
         qmin = min(qmin, minval(chem%ChemSpecies(isp)%conc))
         qmax = max(qmax, maxval(chem%ChemSpecies(isp)%conc))
      end do

      rel = abs(m1 - m0_phaseA) / m0_phaseA
      ok = .true.
      call report('Phase A: mass conservation', rel <= real(TOL, dp), ok)
      call report('Phase A: lower bound preserved', &
         qmin >= 0.5_fp * Q0 * (1.0_fp - TOL), ok)
      call report('Phase A: upper bound preserved', &
         qmax <= 1.5_fp * Q0 * (1.0_fp + TOL), ok)
      call report('Phase A: solution finite', tracer_is_finite(), ok)
      if (.not. ok) all_ok = .false.
   end subroutine check_phase_a

   !> Phase B checks (mass-conserving vertical reconciliation onto PS_NEXT).
   !! With a uniform initial mixing ratio, a uniform non-divergent wind and a
   !! latitude-structured PS_NEXT, the unconditional per-column mass rescale
   !! leaves every column VERTICALLY uniform at the value that conserves the
   !! column tracer mass. This test holds met%DELP fixed on the PS0 grid across
   !! all steps (a real run would advance PS -> PS_NEXT each step), so the
   !! reconciliation ratio (PS0-ptop)/(PS_NEXT-ptop) re-applies every step and
   !! after `nsteps`:
   !!     conc = Q0 * ((PS0 - ptop)/(PS_NEXT - ptop))**nsteps,  ptop = model top.
   !! A uniform field does NOT stay globally uniform here: conserving tracer mass
   !! while the target air-column mass varies with latitude REQUIRES the mixing
   !! ratio to vary with latitude (the mass-conserving default, replacing the
   !! removed mixing-ratio-preserving pfix path).
   subroutine check_phase_b()
      real(fp), allocatable :: ap(:), bp(:)
      real(fp) :: ptop, expect, dev
      !! Tolerance for the reconciliation target. Looser than the exact-mass
      !! Phase A tol: the remap fills the enlarged/shrunk PS_NEXT surface region
      !! by extrapolation and the per-column rescale supplies the mass-conserving
      !! factor R=(PS0-ptop)/(PS_NEXT-ptop); over nsteps of this fixed-DELP setup
      !! that leaves an O(1e-3) REAL32 remap residual on the analytic value.
      real(fp), parameter :: TOL_B = 5.0e-3_fp
      logical  :: ok, okab
      integer  :: s, isp, i, j, k

      call get_hybrid_ab(nz, ap, bp, okab)
      ! Model-top pressure = ap at the edge where bp -> 0 (ordering-agnostic).
      if (bp(1) <= bp(nz + 1)) then
         ptop = ap(1)
      else
         ptop = ap(nz + 1)
      end if

      dev = 0.0_fp
      do s = 1, chem%nSpeciesAdvect
         isp = chem%AdvectIndex(s)
         if (.not. associated(chem%ChemSpecies(isp)%conc)) cycle
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  expect = Q0 * ((PS0 - ptop) / (met%PS_NEXT(i, j) - ptop))**nsteps
                  dev = max(dev, &
                     abs(chem%ChemSpecies(isp)%conc(i, j, k) - expect) / expect)
               end do
            end do
         end do
      end do
      write(output_unit,'(A,ES12.4)') &
         '     max rel. deviation from mass-conserving target = ', dev

      ok = .true.
      call report('Phase B: hybrid Ap/Bp available', okab, ok)
      call report('Phase B: per-column mass-conserving rescale', dev <= TOL_B, ok)
      call report('Phase B: solution finite', tracer_is_finite(), ok)
      if (.not. ok) all_ok = .false.
   end subroutine check_phase_b

   !> Phase C checks: total 3-D mass conservation through horizontal + vertical
   !! transport, proof that the vertical remap redistributes mass off the surface
   !! (above-surface mass > 0), positivity, and finiteness. If the vertical remap
   !! were inactive or a no-op, above_surface_mass would be exactly zero, so the
   !! second check is the direct test that vertical transport is working.
   subroutine check_phase_c()
      real(dp) :: m1, rel, m_above, frac
      real(fp) :: qmin
      integer  :: s, isp
      logical  :: ok

      m1      = total_tracer_mass()       ! live post-run 3-D mass
      m_above = above_surface_mass()      ! mass now residing at k >= 2
      rel     = abs(m1 - m0_phaseC) / m0_phaseC
      frac    = m_above / m1

      qmin = huge(1.0_fp)
      do s = 1, chem%nSpeciesAdvect
         isp = chem%AdvectIndex(s)
         if (.not. associated(chem%ChemSpecies(isp)%conc)) cycle
         qmin = min(qmin, minval(chem%ChemSpecies(isp)%conc))
      end do

      write(output_unit,'(A,ES12.4)') '     above-surface mass fraction = ', frac

      ok = .true.
      call report('Phase C: total mass conserved (horizontal + vertical)', &
         rel <= real(TOL, dp), ok)
      call report('Phase C: vertical remap lifted mass off the surface', &
         m_above > 1.0e-6_dp * m1, ok)
      call report('Phase C: positivity preserved', qmin >= -TOL * Q0, ok)
      call report('Phase C: solution finite', tracer_is_finite(), ok)
      if (.not. ok) all_ok = .false.
   end subroutine check_phase_c

   !> Print a single pass/fail line and AND it into the running result.
   subroutine report(label, passed, acc)
      character(len=*), intent(in)    :: label
      logical,          intent(in)    :: passed
      logical,          intent(inout) :: acc
      if (passed) then
         write(output_unit,'(A,A)') '  ok  ', label
      else
         write(output_unit,'(A,A)') '  XX  ', label
      end if
      acc = acc .and. passed
   end subroutine report

end program test_transport_integration
