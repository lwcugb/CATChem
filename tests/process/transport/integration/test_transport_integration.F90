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
!!   Phase B - pressure fixer ON (PS_NEXT supplied and marked populated):
!!             a UNIFORM tracer with a latitude-structured end-of-step surface
!!             pressure target. The PJC/LLNL correction engages (global periodic
!!             grid, 72 levels, PS_NEXT populated); a uniform field must remain
!!             uniform (the correction mass flux and tracer flux telescope
!!             identically) and stay finite. Constancy is the end-to-end proof
!!             that the fixer flux is mass-consistent.
!!
!! The run reuses the shared standalone configuration
!! (tests/Configs/Default/CATChem_new_config_standalone.yml), which already
!! declares a global periodic 72-level transport block (x_periodic: true,
!! horizontal + vertical, hord/vord = 8). The grid dimensions used here are
!! supplied programmatically via `with_grid`, overriding the file, and only the
!! transport process is added to the pipeline, so the emission/chemistry blocks
!! in that config are never exercised.
!!
!! The numeric accuracy of the fixer correction itself is validated to machine
!! precision by the unit test (pfix Tests 7-8); here we verify the wiring.
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

   !> Latitude-structured end-of-step surface pressure so the PJC/LLNL fixer
   !! produces a non-trivial correction, then mark it populated so the
   !! transport gate (is_field_set) engages the fixer.
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

   !> Phase B checks: a uniform field must stay uniform (the fixer flux
   !! telescopes) and finite.
   subroutine check_phase_b()
      real(fp) :: dev
      integer  :: s, isp
      logical  :: ok
      dev = 0.0_fp
      do s = 1, chem%nSpeciesAdvect
         isp = chem%AdvectIndex(s)
         if (.not. associated(chem%ChemSpecies(isp)%conc)) cycle
         dev = max(dev, maxval(abs(chem%ChemSpecies(isp)%conc - Q0)))
      end do
      ok = .true.
      call report('Phase B: uniform field stays uniform', &
                  dev <= TOL * Q0, ok)
      call report('Phase B: solution finite', tracer_is_finite(), ok)
      if (.not. ok) all_ok = .false.
   end subroutine check_phase_b

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
