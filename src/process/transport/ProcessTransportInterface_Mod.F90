!> \file ProcessTransportInterface_Mod.F90
!! \brief Horizontal (and, later, vertical) tracer transport process for CATChem.
!!
!! Unlike the column-local processes (emission, chemistry, deposition, ...),
!! transport is INHERENTLY GRID-COUPLED: updating a cell requires its horizontal
!! neighbours (a halo) and the grid metrics. It therefore extends the base
!! \ref ProcessInterface (full-field `run(container)`), NOT ColumnProcessInterface.
!!
!! Numerics are provided by the vendored FV3 finite-volume flux operator
!! (\ref fv3_tp_core_mod, i.e. GFDL tp_core.F90). This module owns the CATChem
!! side of the seam: pulling winds/`DELP` from MetState and species from
!! ChemState, building grid metrics, doing the ESMF halo exchange, the Courant
!! sub-cycling, and the flux-form update. See the design plan for the full task
!! breakdown.
!!
!! STATUS: horizontal advection driver active (metrics -> mass fluxes -> halo ->
!! fv_tp_2d -> flux-form update, Courant sub-cycled). Vertical advection is the
!! FV3/GCHP vertically-Lagrangian PPM remap (mappm): the horizontally-deformed
!! layers are conservatively remapped back to the reference hybrid-sigma grid.
!! It is gated by do_vertical and needs no vertical velocity (no OMEGA).
!!
!! \author CATChem Development Team
!! \version 0.2.0
module ProcessTransportInterface_Mod

   ! Core CATChem infrastructure
   use precision_mod, only: fp
   use constants, only: MAX_LEN_NAME, g0
   use iso_fortran_env, only: output_unit
   use ProcessInterface_Mod, only: ProcessInterface
   use StateManager_Mod, only: StateManagerType
   use GridManager_Mod, only: GridManagerType
   use ConfigManager_Mod, only: ConfigManagerType
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Error, CC_Warning, ErrorManagerType
   use ChemState_Mod, only: ChemStateType
   use MetState_Mod, only: MetStateType

   ! Diagnostic system (transport budget fields written to the output files)
   use DiagnosticManager_Mod,   only: DiagnosticManagerType
   use DiagnosticInterface_Mod, only: DiagnosticRegistryType, DiagnosticFieldType, &
      DiagnosticDataType, DIAG_REAL_3D

   ! Vendored FV3 transport kernels + trimmed grid metric types
   use fv3_grid_types_mod, only: fv_grid_type, fv_grid_bounds_type
   use fv3_tp_core_mod, only: fv_tp_2d
   use fv3_vremap_mod, only: mappm
   use met_utilities_mod, only: get_hybrid_ab
   use TransportGridMetrics_Mod, only: build_fv3_grid_metrics
   use TransportHalo_Mod, only: transport_halo_type, halo_update, halo_global_max, &
      halo_global_sum, halo_is_root, &
      HALO_BC_REPLICATE, HALO_BC_PERIODIC
   use TransportMassFlux_Mod, only: fv_mass_flux_type, mass_flux_alloc, mass_flux_free, &
      build_level_mass_flux, courant_max, scale_mass_flux
   use TimeState_Mod, only: TimeStateType

   implicit none
   private

   public :: ProcessTransportInterface

   !> \brief Horizontal tracer-transport process (FV3 flux-form / PPM).
   type, extends(ProcessInterface) :: ProcessTransportInterface
      private

      ! Non-owning handles resolved at init from the StateManager container
      type(ChemStateType), pointer :: chem_state => null()
      type(MetStateType),  pointer :: met_state  => null()
      type(GridManagerType), pointer :: grid_mgr => null()

      ! FV3 kernel grid metrics + bounds (populated in a later milestone)
      type(fv_grid_type)        :: gridstruct
      type(fv_grid_bounds_type) :: bd
      logical :: metrics_ready = .false.

      ! Halo (ghost-cell) policy for the transport data domain. Defaults to the
      ! safe regional / zero-gradient closure; the driver / config sets
      ! x_bc = HALO_BC_PERIODIC for a global-longitude grid.
      type(transport_halo_type) :: halo

      ! Runtime options (parsed from config in a later milestone)
      integer :: hord = 8            !< PPM scheme id passed to fv_tp_2d (8 = monotone)
      integer :: vord = 8            !< PPM scheme id for the vertical remap (mappm kord; 8 = monotone)
      integer :: ng   = 3            !< halo (ghost) width for the FV3 data domain
      logical :: do_horizontal = .true.  !< enable horizontal advection
      logical :: do_vertical   = .true.  !< enable vertical (Lagrangian PPM remap) advection

      ! Diagnostics: when enabled, register + write per-species, per-level
      ! transport mixing-ratio and mass tendency fields to the output.
      ! diag_species selects which advected species to write (empty => all).
      logical :: diagnostics    = .false. !< write transport budget diagnostics
      logical :: diag_registered = .false. !< diagnostic fields registered
      character(len=32), allocatable :: diag_species(:) !< species to write (empty=all)

   contains
      ! Required ProcessInterface implementations
      procedure :: init     => transport_init
      procedure :: run      => transport_run
      procedure :: finalize => transport_finalize

      ! Capability registration
      procedure :: get_required_met_fields => transport_get_required_met_fields

      ! Diagnostic registration (overrides the base no-op)
      procedure :: register_diagnostics => transport_register_diagnostics
   end type ProcessTransportInterface

contains

   !> \brief Initialize the transport process: resolve state/grid handles.
   subroutine transport_init(this, container, rc)
      class(ProcessTransportInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_manager

      rc = CC_SUCCESS
      this%name        = 'transport'
      this%version     = '0.1.0'
      this%description = 'Horizontal tracer transport (FV3 flux-form / PPM kernel)'

      error_manager => container%get_error_manager()

      ! Resolve non-owning handles to the shared state
      this%chem_state => container%get_chem_state_ptr()
      this%met_state  => container%get_met_state_ptr()
      this%grid_mgr   => container%get_grid_manager()

      if (.not. associated(this%chem_state)) then
         if (associated(error_manager)) call error_manager%report_error( &
            1014, 'transport_init: chem state not available', rc)
         rc = CC_FAILURE
         return
      end if
      if (.not. associated(this%met_state)) then
         if (associated(error_manager)) call error_manager%report_error( &
            1014, 'transport_init: met state not available', rc)
         rc = CC_FAILURE
         return
      end if

      ! Grid metrics depend on cell lat/lon/area, which the host populates
      ! AFTER model initialization (e.g. the NUOPC cap assigns MetState%LAT/LON
      ! once the ESMF grid is known). They are therefore built lazily on the
      ! first run() call; see transport_ensure_metrics.
      this%metrics_ready = .false.

      ! Runtime options from configuration; sets active status.
      call transport_load_config(this, container)

      ! Register transport budget diagnostics (per-species, per-level tendency).
      ! Advected-species metadata and the grid shape are populated on the shared
      ! ChemState/GridManager before processes initialize, so this is safe here.
      ! Every PET registers the identical field set, keeping the (collective)
      ! diagnostic write balanced across ranks.
      if (this%diagnostics) then
         call this%register_diagnostics(container, rc)
         if (rc /= CC_SUCCESS) return
      end if

   end subroutine transport_init

   !> \brief Load transport runtime options from the configuration.
   !!
   !! All keys live under `processes/transport/` and are optional (sensible
   !! defaults are applied when absent, and when no config is attached):
   !!   - `activate`   (logical, .true.) : enable/disable the process entirely.
   !!   - `horizontal` (logical, .true.) : run horizontal advection.
   !!   - `vertical`   (logical, .true.) : run the vertical Lagrangian PPM remap
   !!                                      (requires `horizontal`: the remap acts
   !!                                      on the horizontally-deformed layers and
   !!                                      reconciles the surface pressure onto
   !!                                      PS_NEXT).
   !!   - `hord`       (integer, 8)      : FV3 PPM scheme id passed to fv_tp_2d.
   !!   - `vord`       (integer, 8)      : FV3 PPM scheme id (mappm kord) for the
   !!                                      vertical remap.
   !!   - `halo_width` (integer, 3)      : ghost-cell width of the data domain.
   !!   - `x_periodic` (logical, .false.): wrap the west/east (longitude) halo.
   !!   - `y_periodic` (logical, .false.): wrap the south/north (latitude) halo.
   !! A global longitude grid should set `x_periodic: true`; a regional grid
   !! keeps the default zero-gradient (replicated) closure on both axes.
   subroutine transport_load_config(this, container)
      class(ProcessTransportInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container

      type(ConfigManagerType), pointer :: config_manager
      integer :: cfg_rc
      logical :: activate_flag, x_periodic, y_periodic

      activate_flag = .true.
      config_manager => container%get_config_ptr()

      if (associated(config_manager)) then
         call config_manager%get_logical('processes/transport/activate',   activate_flag,      cfg_rc, .true.)
         call config_manager%get_logical('processes/transport/horizontal', this%do_horizontal, cfg_rc, .true.)
         call config_manager%get_logical('processes/transport/vertical',   this%do_vertical,   cfg_rc, .true.)
         call config_manager%get_integer('processes/transport/hord',       this%hord,          cfg_rc, 8)
         call config_manager%get_integer('processes/transport/vord',       this%vord,          cfg_rc, 8)
         call config_manager%get_integer('processes/transport/halo_width', this%ng,            cfg_rc, 3)
         if (this%ng < 1) this%ng = 3
         call config_manager%get_logical('processes/transport/x_periodic', x_periodic,         cfg_rc, .false.)
         call config_manager%get_logical('processes/transport/y_periodic', y_periodic,         cfg_rc, .false.)
         call config_manager%get_logical('processes/transport/diagnostics', this%diagnostics,   cfg_rc, .false.)
         call config_manager%get_array('processes/transport/diag_species', this%diag_species,    cfg_rc)

         if (x_periodic) then
            this%halo%x_bc = HALO_BC_PERIODIC
         else
            this%halo%x_bc = HALO_BC_REPLICATE
         end if
         if (y_periodic) then
            this%halo%y_bc = HALO_BC_PERIODIC
         else
            this%halo%y_bc = HALO_BC_REPLICATE
         end if
      end if

      if (activate_flag) then
         call this%activate()
      else
         call this%deactivate()
      end if

   end subroutine transport_load_config

   !> \brief Build the FV3 grid metrics from MetState geometry (lazy, once).
   !!
   !! Sourced from the geometry CATChem actually populates at run time:
   !! per-cell centre latitude/longitude and (when available) true cell areas.
   !! Corner-based great-circle metrics make this valid for rectilinear,
   !! curvilinear and cubed-sphere grids alike.
   subroutine transport_ensure_metrics(this, rc)
      class(ProcessTransportInterface), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: metrics_rc
      logical :: x_periodic

      rc = CC_SUCCESS
      if (this%metrics_ready) return

      if (.not. associated(this%met_state)) then
         rc = CC_FAILURE
         return
      end if

      ! Geometry may not be populated on the very first step; try again later.
      if (.not. (allocated(this%met_state%LAT) .and. allocated(this%met_state%LON))) then
         rc = CC_FAILURE
         return
      end if

      ! A zonally periodic longitude halo must also wrap the metric halo, so the
      ! seam-cell edge lengths/areas match across the west/east seam and the
      ! mass-flux divergence telescopes to zero (exact tracer-mass conservation).
      x_periodic = (this%halo%x_bc == HALO_BC_PERIODIC)

      if (allocated(this%met_state%AREA_M2)) then
         call build_fv3_grid_metrics(this%met_state%LAT, this%met_state%LON, this%ng, &
            this%gridstruct, this%bd, metrics_rc, &
            area_m2=this%met_state%AREA_M2, x_periodic=x_periodic)
      else
         call build_fv3_grid_metrics(this%met_state%LAT, this%met_state%LON, this%ng, &
            this%gridstruct, this%bd, metrics_rc, &
            x_periodic=x_periodic)
      end if

      if (metrics_rc == CC_SUCCESS) then
         this%metrics_ready = .true.
      else
         rc = CC_FAILURE
      end if

   end subroutine transport_ensure_metrics

   !> \brief Advance tracers one process timestep.
   !!
   !! Horizontal flux-form advection: build the FV3 grid metrics (lazily, on
   !! first use), pull the transport timestep from the shared time state, and
   !! run the per-level / per-species flux operator via transport_horizontal.
   !! When `do_vertical` is set, the horizontally-deformed (Lagrangian) layer
   !! thicknesses are captured and the tracers are conservatively remapped back
   !! to the reference hybrid-sigma grid via transport_vertical (the FV3/GCHP
   !! vertically-Lagrangian PPM remap). Vertical transport is entirely a product
   !! of horizontal mass convergence, so it requires horizontal advection.
   subroutine transport_run(this, container, rc)
      class(ProcessTransportInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer, parameter :: dpk = kind(0.0d0)
      type(TimeStateType), pointer :: time_state
      integer :: metrics_rc, nz, n_adv
      real    :: dt
      real, allocatable :: dp_lag(:,:,:)
      logical :: debug, did_vertical, do_diag
      real(dpk), allocatable :: mass_pre(:), mass_h(:), mass_v(:)
      real(fp),  allocatable :: gmin(:)
      real(fp),  allocatable :: conc_pre(:,:,:,:)

      rc = CC_SUCCESS
      if (.not. this%do_horizontal) return

      ! Build grid metrics on first use (geometry is populated post-init).
      if (.not. this%metrics_ready) then
         call transport_ensure_metrics(this, metrics_rc)
         if (metrics_rc /= CC_SUCCESS) then
            ! Geometry not ready yet (or grid unsuitable): skip advection this
            ! step without failing the model.
            return
         end if
      end if

      if (.not. associated(this%chem_state) .or. .not. associated(this%met_state)) then
         rc = CC_FAILURE
         return
      end if

      ! Nothing to transport if no species are flagged advected.
      if (this%chem_state%nSpeciesAdvect <= 0) return

      ! Transport timestep [s] from the shared time state.
      time_state => container%get_time_state_ptr()
      if (.not. associated(time_state)) then
         rc = CC_FAILURE
         return
      end if
      dt = real(time_state%get_timestep())

      ! Vertical remap is active only when enabled AND the met layer thicknesses
      ! are available; captured once so the debug branch and the run branch agree.
      did_vertical = (this%do_vertical .and. allocated(this%met_state%DELP))

      ! Global-mass budget diagnostics. Reductions inside transport_collect_global
      ! are COLLECTIVE, so every collect is guarded by `debug` only (env-driven,
      ! identical on all ranks) and never by the per-rank `rc`, keeping the
      ! collective call count balanced across PETs.
      debug = transport_debug_enabled()
      if (debug) then
         n_adv = this%chem_state%nSpeciesAdvect
         allocate(mass_pre(n_adv), mass_h(n_adv), mass_v(n_adv), gmin(n_adv))
         mass_pre = 0.0_dpk; mass_h = 0.0_dpk; mass_v = 0.0_dpk; gmin = 0.0_fp
         call transport_collect_global(this, mass_pre, gmin)
      end if

      ! Snapshot the pre-transport concentrations so the net transport tendency
      ! can be written to the output diagnostics after the step. Purely local
      ! (no reductions); gated identically on all PETs.
      do_diag = (this%diagnostics .and. this%diag_registered .and. dt > 0.0)
      if (do_diag) call transport_snapshot_conc(this, conc_pre)

      if (did_vertical) then
         ! Capture the Lagrangian layer thicknesses produced by the horizontal
         ! step so they can be remapped back to the reference grid.
         nz = size(this%met_state%DELP, 3)
         allocate(dp_lag(this%bd%is:this%bd%ie, this%bd%js:this%bd%je, nz))
         call transport_horizontal(this, dt, rc, dp_lag=dp_lag)
         if (debug) call transport_collect_global(this, mass_h, gmin)
         if (rc == CC_SUCCESS) call transport_vertical(this, dp_lag, rc)
         if (debug) call transport_collect_global(this, mass_v, gmin)
         deallocate(dp_lag)
      else
         call transport_horizontal(this, dt, rc)
         if (debug) call transport_collect_global(this, mass_h, gmin)
      end if

      ! Write the per-species, per-level transport tendency (post - pre)/dt to
      ! the diagnostic output fields.
      if (do_diag) then
         if (allocated(conc_pre)) then
            call transport_write_diagnostics(this, container, conc_pre, dt)
            deallocate(conc_pre)
         end if
      end if

      if (debug) then
         call transport_print_budget(this, mass_pre, mass_h, mass_v, gmin, did_vertical)
         deallocate(mass_pre, mass_h, mass_v, gmin)
      end if

   end subroutine transport_run

   !> \brief Horizontal flux-form advection of all advected species.
   !!
   !! For each vertical level: restagger the A-grid winds to the C-grid and build
   !! the FV3 transport quantities (`crx`,`cry`,`xfx`,`yfx`,`ra_x`,`ra_y`,`mfx`,
   !! `mfy`), size the Courant sub-cycling, then for every advected species run
   !! the FV3 flux operator `fv_tp_2d` and apply the mass-conserving flux-form
   !! update
   !!   `q = (q*dp1 + (fx-fx(i+1)+fy-fy(j+1))*rarea) / dp2`,
   !!   `dp2 = dp1 + (mfx-mfx(i+1)+mfy-mfy(j+1))*rarea`.
   !! Winds/`DELP` and the tracer are carried on the FV3 DATA domain (compute
   !! domain + halo ring); the halo is filled by TransportHalo per the configured
   !! per-axis policy (periodic wrap or zero-gradient replication).
   !!
   !! All transport work arrays are bare `real` (promoted to r8 by the library
   !! build flags) to match the vendored kernel and grid metrics; only the
   !! CATChem state (`conc`, `U`/`V`/`DELP`) is `real(fp)`, converted at the seam.
   subroutine transport_horizontal(this, dt, rc, dp_lag)
      class(ProcessTransportInterface), intent(inout) :: this
      real,    intent(in)  :: dt
      integer, intent(out) :: rc
      !> Optional Lagrangian layer thickness [Pa] after horizontal convergence,
      !! shape (is:ie, js:je, nz). When present it is filled per level with the
      !! final evolved `dp2`; the vertical remap consumes it. Requesting it does
      !! not change the horizontal solution.
      real, optional, intent(out) :: dp_lag(:,:,:)

      integer :: is, ie, js, je, isd, ied, jsd, jed
      integer :: nx, ny, nz, npx, npy
      integer :: i, j, k, s, isp, nsplt, it, hrc

      real, parameter :: LIM_FAC = 1.0   !< FV3 PPM limiter factor (default)

      real :: cmax
      type(fv_mass_flux_type) :: mf
      real, allocatable :: q(:,:), ua(:,:), va(:,:), delp(:,:)
      real, allocatable :: dp1(:,:), dp2(:,:)
      real, allocatable :: fx(:,:), fy(:,:)

      rc = CC_SUCCESS

      is  = this%bd%is;  ie  = this%bd%ie;  js  = this%bd%js;  je  = this%bd%je
      isd = this%bd%isd; ied = this%bd%ied; jsd = this%bd%jsd; jed = this%bd%jed
      nx  = ie; ny = je
      npx = nx + 1; npy = ny + 1

      if (.not. (allocated(this%met_state%U)    .and. &
         allocated(this%met_state%V)    .and. &
         allocated(this%met_state%DELP))) then
         rc = CC_FAILURE
         return
      end if
      nz = size(this%met_state%DELP, 3)

      ! Data-domain work arrays (compute domain + halo ring).
      allocate(q   (isd:ied, jsd:jed))
      allocate(ua  (isd:ied, jsd:jed))
      allocate(va  (isd:ied, jsd:jed))
      allocate(delp(isd:ied, jsd:jed))
      allocate(dp1 (isd:ied, jsd:jed))
      allocate(dp2 (isd:ied, jsd:jed))
      ! Face fluxes returned by the kernel (compute-domain faces).
      allocate(fx(is:ie+1, js:je))
      allocate(fy(is:ie,   js:je+1))

      call mass_flux_alloc(mf, this%bd, rc)
      if (rc /= CC_SUCCESS) return

      ! -----------------------------------------------------------------------
      ! No horizontal pressure fixer.  Following UFS-ATM / GCHP, tracers are
      ! advected with the wind-derived C-grid mass fluxes and the resulting
      ! (Lagrangian) surface pressure is reconciled to the met end-of-step
      ! surface pressure (PS_NEXT) entirely by the vertical remap in
      ! transport_vertical (its target grid is the reference hybrid coordinate
      ! at PS_NEXT).  This is fully column-local -> grid-agnostic (lon-lat and
      ! cubed-sphere) and decomposes naturally across PEs/tiles.
      ! -----------------------------------------------------------------------
      level_loop: do k = 1, nz

         ! --- Winds + layer thickness on the data domain (interior + halo) ---
         ua = 0.0; va = 0.0; delp = 0.0
         ua  (is:ie, js:je) = this%met_state%U   (:,:,k)
         va  (is:ie, js:je) = this%met_state%V   (:,:,k)
         delp(is:ie, js:je) = this%met_state%DELP(:,:,k)
         call halo_update(ua,   this%bd, this%halo, hrc)
         call halo_update(va,   this%bd, this%halo, hrc)
         call halo_update(delp, this%bd, this%halo, hrc)

         ! --- C-grid mass fluxes + Courant sub-cycling for this level --------
         call build_level_mass_flux(this%gridstruct, this%bd, ua, va, delp, dt, mf, rc)
         if (rc /= CC_SUCCESS) exit level_loop

         cmax  = courant_max(mf, this%bd)
         ! nsplt drives how many times the paired halo exchange is called in the
         ! sub-cycle below, so it MUST be identical on every PET. courant_max is
         ! a LOCAL (per-PET) maximum; reduce it to the global maximum first, or a
         ! PET with faster local winds would sub-cycle (and halo-exchange) more
         ! times than its neighbours and the off-PET Sendrecv would deadlock.
         ! (Serial / single-PET: halo_global_max is the identity.)
         call halo_global_max(cmax, hrc)
         nsplt = int(1.0 + cmax)
         if (nsplt < 1) nsplt = 1
         if (nsplt > 1) call scale_mass_flux(mf, this%gridstruct, this%bd, 1.0 / real(nsplt))

         ! --- Lagrangian layer thickness (species-independent) ---------------
         ! The flux-form evolution is linear in the (per-substep) mass flux, so
         ! after nsplt substeps every species reaches the same evolved dp2:
         !   dp_lag = delp + nsplt * (mfx-mfx(i+1)+mfy-mfy(j+1)) * rarea.
         ! Its column sum is the Lagrangian surface pressure produced by the
         ! winds; the vertical remap later reconciles it onto PS_NEXT.
         if (present(dp_lag)) then
            do j = js, je
               do i = is, ie
                  dp_lag(i,j,k) = delp(i,j) + real(nsplt) * &
                     (mf%mfx(i,j) - mf%mfx(i+1,j) + mf%mfy(i,j) - mf%mfy(i,j+1)) &
                     * this%gridstruct%rarea(i,j)
               end do
            end do
         end if

         ! --- Advect every species on this level -----------------------------
         species_loop: do s = 1, this%chem_state%nSpeciesAdvect
            isp = this%chem_state%AdvectIndex(s)
            if (isp < 1 .or. isp > size(this%chem_state%ChemSpecies)) cycle species_loop
            if (.not. associated(this%chem_state%ChemSpecies(isp)%conc)) cycle species_loop

            q = 0.0
            q(is:ie, js:je) = this%chem_state%ChemSpecies(isp)%conc(:,:,k)
            call halo_update(q, this%bd, this%halo, hrc)

            ! Layer thickness integrated with the (scaled) mass flux per sub-step.
            dp1(is:ie, js:je) = delp(is:ie, js:je)

            subcycle: do it = 1, nsplt
               call fv_tp_2d(q, mf%crx, mf%cry, npx, npy, this%hord, fx, fy, &
                  mf%xfx, mf%yfx, this%gridstruct, this%bd, &
                  mf%ra_x, mf%ra_y, LIM_FAC, mfx=mf%mfx, mfy=mf%mfy)

               do j = js, je
                  do i = is, ie
                     dp2(i,j) = dp1(i,j) + &
                        (mf%mfx(i,j) - mf%mfx(i+1,j) + mf%mfy(i,j) - mf%mfy(i,j+1)) &
                        * this%gridstruct%rarea(i,j)
                     q(i,j) = (q(i,j) * dp1(i,j) + &
                        (fx(i,j) - fx(i+1,j) + fy(i,j) - fy(i,j+1)) &
                        * this%gridstruct%rarea(i,j)) / dp2(i,j)
                  end do
               end do

               if (it < nsplt) then
                  call halo_update(q, this%bd, this%halo, hrc)
                  dp1(is:ie, js:je) = dp2(is:ie, js:je)
               end if
            end do subcycle

            ! Clip any tiny negative concentrations (e.g. from a cold start or
            ! PPM undershoot) to zero when writing back the tracer field.
            this%chem_state%ChemSpecies(isp)%conc(:,:,k) = &
               max(real(q(is:ie, js:je), fp), 0.0_fp)
         end do species_loop

      end do level_loop

      call mass_flux_free(mf)
      if (allocated(q))    deallocate(q)
      if (allocated(ua))   deallocate(ua)
      if (allocated(va))   deallocate(va)
      if (allocated(delp)) deallocate(delp)
      if (allocated(dp1))  deallocate(dp1)
      if (allocated(dp2))  deallocate(dp2)
      if (allocated(fx))   deallocate(fx)
      if (allocated(fy))   deallocate(fy)

   end subroutine transport_horizontal

   !> \brief Vertical (Lagrangian) conservative PPM remap of all advected species.
   !!
   !! This is the FV3/GCHP vertically-Lagrangian step: after horizontal advection
   !! the layers have deformed (their thicknesses are `dp_lag`), carrying the
   !! tracers on floating Lagrangian surfaces. Here every column is remapped back
   !! to the reference hybrid-sigma grid with the vendored FV3 `mappm` monotone
   !! PPM operator. The vertical redistribution of tracer mass IS this remap;
   !! there is no separate vertical-velocity advection (OMEGA is not used), which
   !! is exactly how UFS-ATM and GCHP transport constituents in the vertical.
   !!
   !! Conservation: `mappm` preserves the mass-weighted integral \f$\sum_k q\,\Delta p\f$
   !! exactly in the interior, but its top cell (whose edge coincides with the
   !! model top) takes the top source-cell mean rather than that cell's sub-cell
   !! average, leaving a small residual. The target column pressures are the
   !! reference hybrid coordinate evaluated at the end-of-step met surface
   !! pressure `PS_NEXT` when the host has supplied it, so the remap reconciles
   !! the wind-driven (Lagrangian) surface pressure onto the met surface pressure
   !! -- this is the FV3/GCHP way to close pressure and replaces the (removed)
   !! horizontal PJC pressure fixer. The model top edge matches by construction
   !! (Bp(top)=0); the surface edge is left at `PS_NEXT` so the target column air
   !! mass equals the met air mass (`PS = PS_NEXT`, pressure-consistent).
   !!
   !! Mass conservation: each positive-definite column is always rescaled by the
   !! ratio of its pre- to post-remap tracer mass, so the column (hence global)
   !! tracer mass is conserved to machine precision while the surface still sits
   !! at `PS_NEXT` (pressure stays consistent). Because wind-derived fluxes and an
   !! independently specified `PS_NEXT` are not exactly consistent (closing that
   !! gap exactly needs an elliptic pressure fixer), the small residual is
   !! absorbed as a uniform per-column mixing-ratio scaling. When the met winds
   !! and `PS` come from the same source the residual is ~0 and the rescale is a
   !! no-op (ratio ~= 1), so it is a harmless safety guard; when they are not, it
   !! keeps mass exact. The fix is purely column-local, so it is grid-agnostic
   !! (identical on lat-lon and cubed-sphere) and needs no global solve.
   !!
   !! When `PS_NEXT` is not populated (e.g. the final met slice), the target
   !! falls back to the Lagrangian surface pressure \f$p_{s}^{lag}=\sum_k \Delta p^{lag}\f$
   !! with the surface edge forced to match the source; that self-consistent case
   !! is rescaled by the same code so each column's tracer mass is conserved to
   !! machine precision.
   !!
   !! Ordering: CATChem is surface-first (`k=1` is the surface); `mappm` expects
   !! edges from model top to surface (increasing pressure), so columns are flipped
   !! on the way in and out.
   subroutine transport_vertical(this, dp_lag, rc)
      class(ProcessTransportInterface), intent(inout) :: this
      !> Lagrangian layer thickness [Pa], shape (is:ie, js:je, nz), surface-first.
      real,    intent(in)  :: dp_lag(:,:,:)
      integer, intent(out) :: rc

      integer :: is, ie, js, je, nz, np, i, j, k, c, s, isp
      logical :: ok, use_psn
      real(fp), allocatable :: ap(:), bp(:)
      real,     allocatable :: pe1(:,:), pe2(:,:), q1(:,:), q2(:,:)
      real                  :: ps_lag, ps_tgt, m_src, m_tgt

      rc = CC_SUCCESS

      is = this%bd%is; ie = this%bd%ie; js = this%bd%js; je = this%bd%je
      nz = size(dp_lag, 3)
      np = (ie - is + 1) * (je - js + 1)
      if (np < 1 .or. nz < 1) return

      ! Reference hybrid-sigma Ap/Bp (surface-first, Pa). If this level count has
      ! no table, skip the remap without failing the model.
      call get_hybrid_ab(nz, ap, bp, ok)
      if (.not. ok) return

      allocate(pe1(np, nz+1), pe2(np, nz+1))
      allocate(q1(np, nz), q2(np, nz))

      ! Reconcile onto the met end-of-step surface pressure when the host has
      ! populated it this step (tracked by the per-timestep populated-field
      ! registry); otherwise fall back to the self-consistent Lagrangian PS.
      use_psn = this%met_state%is_field_set('PS_NEXT') .and. &
         allocated(this%met_state%PS_NEXT)
      if (use_psn) then
         if (size(this%met_state%PS_NEXT,1) /= size(dp_lag,1) .or. &
            size(this%met_state%PS_NEXT,2) /= size(dp_lag,2)) use_psn = .false.
      end if

      ! --- Column edge pressures (top -> bottom), species-independent ---------
      c = 0
      do j = js, je
         do i = is, ie
            c = c + 1
            ! Source edges from the Lagrangian thicknesses. Model top pressure is
            ! ap(nz+1) (Bp(top)=0); accumulate downward. FV3 layer k (top-first)
            ! corresponds to CATChem layer nz+1-k (surface-first).
            pe1(c,1) = real(ap(nz+1))
            do k = 1, nz
               pe1(c,k+1) = pe1(c,k) + dp_lag(i,j, nz+1-k)
            end do
            ps_lag = pe1(c,nz+1)

            ! Target surface pressure: the met end-of-step PS when supplied
            ! (reconcile), else the Lagrangian PS (self-consistent fallback).
            if (use_psn) then
               ps_tgt = real(this%met_state%PS_NEXT(i,j))
            else
               ps_tgt = ps_lag
            end if

            ! Target edges = reference hybrid coordinate at ps_tgt. FV3 edge k
            ! (top-first) is surface-first edge nz+2-k.
            do k = 1, nz+1
               pe2(c,k) = real(ap(nz+2-k)) + real(bp(nz+2-k)) * ps_tgt
            end do
            ! Model top matches by construction (Bp(top)=0). For the fallback we
            ! also force the surface edge so each column conserves mass exactly;
            ! for the PS_NEXT target the surface edge stays at PS_NEXT so the
            ! remap adds/removes the reconciling surface mass.
            pe2(c,1) = pe1(c,1)
            if (.not. use_psn) pe2(c,nz+1) = pe1(c,nz+1)
         end do
      end do

      ! --- Remap every advected species ---------------------------------------
      species_loop: do s = 1, this%chem_state%nSpeciesAdvect
         isp = this%chem_state%AdvectIndex(s)
         if (isp < 1 .or. isp > size(this%chem_state%ChemSpecies)) cycle species_loop
         if (.not. associated(this%chem_state%ChemSpecies(isp)%conc)) cycle species_loop

         c = 0
         do j = js, je
            do i = is, ie
               c = c + 1
               do k = 1, nz
                  q1(c,k) = this%chem_state%ChemSpecies(isp)%conc(i,j, nz+1-k)
               end do
            end do
         end do

         call mappm(nz, pe1, q1, nz, pe2, q2, 1, np, 0, this%vord)

         ! Clip any small negative concentrations from the PPM remap to zero
         ! BEFORE the mass fixer, so the per-column rescale below redistributes
         ! onto the non-negative profile and column (hence global) tracer mass
         ! stays both conserved AND non-negative.
         do c = 1, np
            do k = 1, nz
               if (q2(c,k) < 0.0) q2(c,k) = 0.0
            end do
         end do

         ! --- Column tracer-mass fixer ----------------------------------------
         ! Rescale each positive-definite column by the ratio of its pre- to
         ! post-remap tracer mass so the column (hence global) tracer mass is
         ! conserved to machine precision. On the fallback path (target =
         ! Lagrangian PS) this also absorbs mappm's small top-cell residual. On
         ! the PS_NEXT path the surface still sits at PS_NEXT (pressure stays
         ! consistent) but the tracer mass is held fixed, so the wind/PS
         ! inconsistency becomes a uniform per-column mixing-ratio scaling instead
         ! of a burden change. When the met winds and PS come from the same source
         ! the ratio is ~1 (a no-op, harmless); when they are not, it keeps mass
         ! exact. Purely column-local, hence grid-agnostic (lat-lon and
         ! cubed-sphere) with no global solve.
         do c = 1, np
            m_src = 0.0
            m_tgt = 0.0
            do k = 1, nz
               m_src = m_src + q1(c,k) * (pe1(c,k+1) - pe1(c,k))
               m_tgt = m_tgt + q2(c,k) * (pe2(c,k+1) - pe2(c,k))
            end do
            if (m_src > 0.0 .and. m_tgt > 0.0) then
               do k = 1, nz
                  q2(c,k) = q2(c,k) * (m_src / m_tgt)
               end do
            end if
         end do

         c = 0
         do j = js, je
            do i = is, ie
               c = c + 1
               do k = 1, nz
                  this%chem_state%ChemSpecies(isp)%conc(i,j, nz+1-k) = real(q2(c,k), fp)
               end do
            end do
         end do
      end do species_loop

      deallocate(pe1, pe2, q1, q2)
      if (allocated(ap)) deallocate(ap)
      if (allocated(bp)) deallocate(bp)

   end subroutine transport_vertical

   !> \brief Finalize: release handles.
   subroutine transport_finalize(this, rc)
      class(ProcessTransportInterface), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      nullify(this%chem_state)
      nullify(this%met_state)
      nullify(this%grid_mgr)
      this%metrics_ready = .false.

      call this%deactivate()

   end subroutine transport_finalize

   !> \brief Meteorological fields required by transport.
   subroutine transport_get_required_met_fields(this, field_names)
      class(ProcessTransportInterface), intent(in) :: this
      character(len=MAX_LEN_NAME), allocatable, intent(out) :: field_names(:)

      ! Horizontal advection needs the winds and the layer pressure thickness.
      ! Vertical transport is a Lagrangian PPM remap driven by horizontal mass
      ! convergence (no vertical velocity / OMEGA required).
      allocate(field_names(3))
      field_names(1) = 'U'
      field_names(2) = 'V'
      field_names(3) = 'DELP'

   end subroutine transport_get_required_met_fields

   !> \brief Whether opt-in transport diagnostics are enabled.
   !!
   !! Controlled by the environment variable CATCHEM_TRANSPORT_DEBUG: any value
   !! other than unset / empty / "0" / "false" / "no" / "off" turns per-call
   !! reporting on. The result is cached on first query so the environment is
   !! read only once per run.
   logical function transport_debug_enabled()
      logical, save      :: checked = .false.
      logical, save      :: enabled = .false.
      character(len=32)  :: val
      integer            :: length, status

      if (.not. checked) then
         call get_environment_variable('CATCHEM_TRANSPORT_DEBUG', val, length, status)
         enabled = (status == 0 .and. length > 0)
         if (enabled) then
            select case (trim(adjustl(val)))
             case ('0', 'false', 'FALSE', 'no', 'NO', 'off', 'OFF')
               enabled = .false.
            end select
         end if
         checked = .true.
      end if
      transport_debug_enabled = enabled
   end function transport_debug_enabled

   !> \brief Collect the global (all-PET) tracer mass and global minimum per
   !!        advected species.
   !!
   !! For each advected species this sums the local tracer mass
   !! `sum(conc*DELP*AREA_M2)` and reduces it across every PET with
   !! `halo_global_sum`, and reduces the global minimum concentration with
   !! `halo_global_max` (min carried as -max so a single collective covers it).
   !! The per-PET local totals cannot show HORIZONTAL conservation on a
   !! decomposed run because advection moves mass ACROSS PET boundaries -- only
   !! the sum over every rank is conserved. Both reductions are the serial
   !! identity when no MPI backend is registered, so this is safe in every
   !! configuration. The reductions are COLLECTIVE, so they run for every
   !! advected species on every PET (mass 0 / min +inf for any locally skipped
   !! slot) to stay balanced.
   !!
   !! \param[in]  this  transport process (needs chem_state, met_state)
   !! \param[out] gmass global tracer mass per advected species [conc*Pa*m2]
   !! \param[out] gmin  global minimum concentration per advected species [conc]
   subroutine transport_collect_global(this, gmass, gmin)
      class(ProcessTransportInterface), intent(in) :: this
      integer, parameter :: dpk = kind(0.0d0)
      real(dpk), intent(out) :: gmass(:)
      real(fp),  intent(out) :: gmin(:)
      integer   :: nxl, nyl, nzl, i, j, k, s, isp, grc
      real(dpk) :: mass, cell_area
      real      :: gm, gmn
      logical   :: have_area

      gmass = 0.0_dpk
      gmin  = 0.0_fp
      if (.not. associated(this%chem_state) .or. &
         .not. associated(this%met_state)) return
      if (.not. allocated(this%met_state%DELP)) return
      have_area = allocated(this%met_state%AREA_M2)

      do s = 1, this%chem_state%nSpeciesAdvect
         isp  = this%chem_state%AdvectIndex(s)
         mass = 0.0_dpk
         gmn  = -huge(1.0)   ! min carried as -max; -huge => "no local data"
         if (isp >= 1 .and. isp <= size(this%chem_state%ChemSpecies)) then
            if (associated(this%chem_state%ChemSpecies(isp)%conc)) then
               associate (conc => this%chem_state%ChemSpecies(isp)%conc)
                  nxl = size(conc, 1); nyl = size(conc, 2); nzl = size(conc, 3)
                  do k = 1, nzl
                     do j = 1, nyl
                        do i = 1, nxl
                           if (have_area) then
                              cell_area = real(this%met_state%AREA_M2(i,j), dpk)
                           else
                              cell_area = 1.0_dpk
                           end if
                           mass = mass + real(conc(i,j,k), dpk) * &
                              real(this%met_state%DELP(i,j,k), dpk) * cell_area
                           gmn  = max(gmn, -real(conc(i,j,k)))
                        end do
                     end do
                  end do
               end associate
            end if
         end if
         gm = real(mass)
         call halo_global_sum(gm, grc)    ! collective; identity in serial
         call halo_global_max(gmn, grc)   ! collective; identity in serial
         gmass(s) = real(gm, dpk)
         gmin(s)  = real(-gmn, fp)
      end do
   end subroutine transport_collect_global

   !> \brief Print one combined global-mass budget line per advected species.
   !!
   !! Reports, on the root PET only, the global tracer mass before transport
   !! (`pre`), after horizontal advection (`post_h`) and, when the vertical remap
   !! ran, after it (`post_v`), plus the relative changes `dH=(post_h-pre)/pre`
   !! and `dV=(post_v-post_h)/post_h` and the global minimum concentration
   !! (`min`, confirms the non-negativity clip). Interpretation:
   !!   * horizontal advection conserves tracer mass, so `dH` should be ~0
   !!     (machine/round-off level, ~1e-6..1e-8 in single precision);
   !!   * the vertical remap conserves column mass on the fallback path (`dV`~0);
   !!     on the PS_NEXT path a small `dV` is EXPECTED (pressure reconciliation to
   !!     the met air mass), not a leak.
   !! Tracers that are identically zero (no burden) are skipped to cut clutter.
   subroutine transport_print_budget(this, mass_pre, mass_h, mass_v, gmin, did_vertical)
      class(ProcessTransportInterface), intent(in) :: this
      integer, parameter :: dpk = kind(0.0d0)
      real(dpk), intent(in) :: mass_pre(:), mass_h(:), mass_v(:)
      real(fp),  intent(in) :: gmin(:)
      logical,   intent(in) :: did_vertical
      real(dpk), parameter  :: tiny_m = 1.0e-300_dpk
      integer   :: s, isp
      real(dpk) :: dh, dv

      if (.not. halo_is_root()) return

      if (did_vertical) then
         write(output_unit,'(A)') &
            '  [transport-budget] global tracer mass  (pre -> post-horizontal -> post-vertical)'
      else
         write(output_unit,'(A)') &
            '  [transport-budget] global tracer mass  (pre -> post-horizontal)'
      end if

      do s = 1, this%chem_state%nSpeciesAdvect
         isp = this%chem_state%AdvectIndex(s)
         if (isp < 1 .or. isp > size(this%chem_state%ChemSpecies)) cycle
         ! Skip tracers that carry no mass at all (keeps the report short).
         if (mass_pre(s) <= tiny_m .and. mass_h(s) <= tiny_m) cycle
         dh = (mass_h(s) - mass_pre(s)) / max(mass_pre(s), tiny_m)
         if (did_vertical) then
            dv = (mass_v(s) - mass_h(s)) / max(mass_h(s), tiny_m)
            write(output_unit, &
               '(A,A,A,ES20.12,A,ES20.12,A,ES20.12,A,ES9.2,A,ES9.2,A,ES9.2)') &
               '    ', trim(this%chem_state%ChemSpecies(isp)%short_name), &
               ': pre=', mass_pre(s), ' post_h=', mass_h(s), ' post_v=', mass_v(s), &
               ' dH=', dh, ' dV=', dv, ' min=', gmin(s)
         else
            write(output_unit, &
               '(A,A,A,ES20.12,A,ES20.12,A,ES9.2,A,ES9.2)') &
               '    ', trim(this%chem_state%ChemSpecies(isp)%short_name), &
               ': pre=', mass_pre(s), ' post_h=', mass_h(s), &
               ' dH=', dh, ' min=', gmin(s)
         end if
      end do
   end subroutine transport_print_budget

   !> \brief Register the per-species, per-level transport budget diagnostics.
   !!
   !! Registers one DIAG_REAL_3D field `transport_tend_<species>` for every
   !! advected species on the shared 'transport' diagnostic registry. Called from
   !! transport_init (advected-species metadata and the grid shape are populated
   !! before processes initialize). Every PET registers the identical field set,
   !! so the collective diagnostic write stays balanced across ranks. Best effort:
   !! if the diagnostic manager / grid is unavailable this quietly leaves
   !! diagnostics off rather than failing the model.
   subroutine transport_register_diagnostics(this, container, rc)
      class(ProcessTransportInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(DiagnosticManagerType),  pointer :: diag_mgr => null()
      type(DiagnosticRegistryType), pointer :: registry => null()
      type(GridManagerType),        pointer :: grid_mgr => null()
      character(len=MAX_LEN_NAME) :: field_name, sp_name
      integer :: s, isp, nx, ny, nz, dims_3d(3)

      rc = CC_SUCCESS
      if (.not. this%diagnostics) return
      if (.not. associated(this%chem_state)) return
      if (this%chem_state%nSpeciesAdvect <= 0) return

      diag_mgr => container%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) return
      grid_mgr => container%get_grid_manager()
      if (.not. associated(grid_mgr)) return

      call grid_mgr%get_shape(nx, ny, nz)
      dims_3d = [nx, ny, nz]

      call diag_mgr%register_process('transport', rc)
      if (rc /= CC_SUCCESS) return
      call diag_mgr%get_process_registry('transport', registry, rc)
      if (rc /= CC_SUCCESS) return
      if (.not. associated(registry)) then
         rc = CC_FAILURE
         return
      end if

      do s = 1, this%chem_state%nSpeciesAdvect
         isp = this%chem_state%AdvectIndex(s)
         if (isp < 1 .or. isp > size(this%chem_state%ChemSpecies)) cycle
         sp_name = this%chem_state%ChemSpecies(isp)%short_name
         if (.not. transport_species_selected(this, trim(sp_name))) cycle

         ! Net transport tendency of the tracer mixing ratio, per level.
         write(field_name, '(A,A)') 'transport_tend_', trim(sp_name)
         call this%register_diagnostic_field(registry, trim(field_name), &
            'Net transport tendency of '//trim(sp_name)//' per level', &
            'kg kg-1 s-1', DIAG_REAL_3D, 'transport', dims_3d, rc=rc)
         if (rc /= CC_SUCCESS) return

         ! Net transport tendency of the tracer MASS in each layer [kg/s].
         write(field_name, '(A,A)') 'transport_mass_', trim(sp_name)
         call this%register_diagnostic_field(registry, trim(field_name), &
            'Net transport mass tendency of '//trim(sp_name)//' per level', &
            'kg s-1', DIAG_REAL_3D, 'transport', dims_3d, rc=rc)
         if (rc /= CC_SUCCESS) return
      end do

      this%diag_registered = .true.
   end subroutine transport_register_diagnostics

   !> \brief True if `short_name` is selected for diagnostic output.
   !!
   !! An empty (or unallocated) `diag_species` list selects *all* advected
   !! species, matching the convention of the global diagnostic list. A
   !! non-empty list restricts output to a case-insensitive match.
   logical function transport_species_selected(this, short_name) result(sel)
      class(ProcessTransportInterface), intent(in) :: this
      character(len=*), intent(in) :: short_name
      integer :: i

      sel = .true.
      if (.not. allocated(this%diag_species)) return
      if (size(this%diag_species) == 0) return
      sel = .false.
      do i = 1, size(this%diag_species)
         if (len_trim(this%diag_species(i)) == 0) cycle
         if (transport_same_name(this%diag_species(i), short_name)) then
            sel = .true.
            return
         end if
      end do
   end function transport_species_selected

   !> \brief Case-insensitive, trimmed comparison of two species names.
   logical function transport_same_name(a, b) result(same)
      character(len=*), intent(in) :: a, b
      same = (trim(transport_lower(adjustl(a))) == trim(transport_lower(adjustl(b))))
   end function transport_same_name

   !> \brief Lowercase an ASCII string.
   pure function transport_lower(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: out
      integer :: i, ic
      do i = 1, len(s)
         ic = iachar(s(i:i))
         if (ic >= iachar('A') .and. ic <= iachar('Z')) then
            out(i:i) = achar(ic + 32)
         else
            out(i:i) = s(i:i)
         end if
      end do
   end function transport_lower

   !> \brief Snapshot the pre-transport concentrations of every advected species.
   !!
   !! Fills `conc_pre(nx,ny,nz,nSpeciesAdvect)` (allocated here) with the current
   !! tracer fields so the net transport tendency can be formed after the step.
   !! Purely local; left unallocated if no advected species has a valid field.
   subroutine transport_snapshot_conc(this, conc_pre)
      class(ProcessTransportInterface), intent(in) :: this
      real(fp), allocatable, intent(out) :: conc_pre(:,:,:,:)
      integer :: s, isp, nx, ny, nz, n_adv

      n_adv = this%chem_state%nSpeciesAdvect
      nx = 0; ny = 0; nz = 0
      do s = 1, n_adv
         isp = this%chem_state%AdvectIndex(s)
         if (isp < 1 .or. isp > size(this%chem_state%ChemSpecies)) cycle
         if (.not. associated(this%chem_state%ChemSpecies(isp)%conc)) cycle
         nx = size(this%chem_state%ChemSpecies(isp)%conc, 1)
         ny = size(this%chem_state%ChemSpecies(isp)%conc, 2)
         nz = size(this%chem_state%ChemSpecies(isp)%conc, 3)
         exit
      end do
      if (nx == 0) return

      allocate(conc_pre(nx, ny, nz, n_adv))
      conc_pre = 0.0_fp
      do s = 1, n_adv
         isp = this%chem_state%AdvectIndex(s)
         if (isp < 1 .or. isp > size(this%chem_state%ChemSpecies)) cycle
         if (.not. associated(this%chem_state%ChemSpecies(isp)%conc)) cycle
         conc_pre(:,:,:,s) = this%chem_state%ChemSpecies(isp)%conc
      end do
   end subroutine transport_snapshot_conc

   !> \brief Write the transport diagnostics (mixing-ratio + mass tendency) per species.
   !!
   !! For every selected advected species stores, per level:
   !!  * `transport_tend_<species>`     = (conc - conc_pre)/dt            [kg kg-1 s-1]
   !!  * `transport_mass_<species>`     = (conc - conc_pre)/dt * DELP*AREA/g0  [kg s-1]
   !!
   !! Both are genuine transport tendencies (post minus pre over the step), the
   !! second weighted by the layer air mass so it reads as a tracer-mass budget
   !! term. Purely local (no reductions); best effort: any field that is missing
   !! / not ready / shape-mismatched is skipped. The mass tendency is skipped
   !! when DELP or the cell area are unavailable.
   subroutine transport_write_diagnostics(this, container, conc_pre, dt)
      class(ProcessTransportInterface), intent(in) :: this
      type(StateManagerType), intent(inout) :: container
      real(fp), intent(in) :: conc_pre(:,:,:,:)
      real,     intent(in) :: dt

      type(DiagnosticManagerType),  pointer :: diag_mgr   => null()
      type(DiagnosticRegistryType), pointer :: registry   => null()
      type(DiagnosticFieldType),    pointer :: diag_field => null()
      type(DiagnosticDataType),     pointer :: diag_data  => null()
      real(fp), pointer :: fptr(:,:,:) => null()
      real(fp), pointer :: conc(:,:,:) => null()
      character(len=MAX_LEN_NAME) :: field_name, sp_name
      logical :: have_mass
      integer :: s, isp, drc, i, j, k, nx, ny, nz

      if (dt <= 0.0) return
      diag_mgr => container%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) return
      call diag_mgr%get_process_registry('transport', registry, drc)
      if (drc /= CC_SUCCESS .or. .not. associated(registry)) return

      have_mass = associated(this%met_state)
      if (have_mass) have_mass = allocated(this%met_state%DELP) .and. &
         allocated(this%met_state%AREA_M2)

      do s = 1, this%chem_state%nSpeciesAdvect
         isp = this%chem_state%AdvectIndex(s)
         if (isp < 1 .or. isp > size(this%chem_state%ChemSpecies)) cycle
         if (.not. associated(this%chem_state%ChemSpecies(isp)%conc)) cycle
         sp_name = this%chem_state%ChemSpecies(isp)%short_name
         if (.not. transport_species_selected(this, trim(sp_name))) cycle
         conc => this%chem_state%ChemSpecies(isp)%conc

         ! --- Net transport tendency of the mixing ratio [kg kg-1 s-1] ---
         write(field_name, '(A,A)') 'transport_tend_', trim(sp_name)
         diag_field => registry%get_field_ptr(trim(field_name))
         if (associated(diag_field)) then
            if (diag_field%is_ready()) then
               diag_data => diag_field%get_data_ptr()
               if (associated(diag_data)) then
                  fptr => diag_data%get_real_3d_ptr()
                  if (associated(fptr)) then
                     if (size(fptr,1) == size(conc_pre,1) .and. &
                        size(fptr,2) == size(conc_pre,2) .and. &
                        size(fptr,3) == size(conc_pre,3)) then
                        fptr = (conc - conc_pre(:,:,:,s)) / dt
                     end if
                  end if
               end if
            end if
         end if

         ! --- Net transport tendency of the tracer mass [kg s-1] ---
         if (.not. have_mass) cycle
         write(field_name, '(A,A)') 'transport_mass_', trim(sp_name)
         diag_field => registry%get_field_ptr(trim(field_name))
         if (.not. associated(diag_field)) cycle
         if (.not. diag_field%is_ready()) cycle
         diag_data => diag_field%get_data_ptr()
         if (.not. associated(diag_data)) cycle
         fptr => diag_data%get_real_3d_ptr()
         if (.not. associated(fptr)) cycle
         nx = size(conc, 1); ny = size(conc, 2); nz = size(conc, 3)
         if (size(fptr,1) /= nx .or. size(fptr,2) /= ny .or. size(fptr,3) /= nz) cycle
         if (size(conc_pre,1) /= nx .or. size(conc_pre,2) /= ny .or. size(conc_pre,3) /= nz) cycle
         if (size(this%met_state%DELP,1) /= nx .or. size(this%met_state%DELP,2) /= ny .or. &
            size(this%met_state%DELP,3) /= nz) cycle
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  fptr(i,j,k) = (conc(i,j,k) - conc_pre(i,j,k,s)) / dt * &
                     this%met_state%DELP(i,j,k) * &
                     this%met_state%AREA_M2(i,j) / g0
               end do
            end do
         end do
      end do
   end subroutine transport_write_diagnostics

end module ProcessTransportInterface_Mod
