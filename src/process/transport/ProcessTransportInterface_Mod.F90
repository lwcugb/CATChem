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
   use ProcessInterface_Mod, only: ProcessInterface
   use StateManager_Mod, only: StateManagerType
   use GridManager_Mod, only: GridManagerType
   use ConfigManager_Mod, only: ConfigManagerType
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Error, CC_Warning, ErrorManagerType
   use ChemState_Mod, only: ChemStateType
   use MetState_Mod, only: MetStateType

   ! Vendored FV3 transport kernels + trimmed grid metric types
   use fv3_grid_types_mod, only: fv_grid_type, fv_grid_bounds_type
   use fv3_tp_core_mod, only: fv_tp_2d
   use fv3_vremap_mod, only: mappm
   use fv3_pfix_mod, only: pfix_correction
   use met_utilities_mod, only: get_hybrid_ab
   use TransportGridMetrics_Mod, only: build_fv3_grid_metrics
   use TransportHalo_Mod, only: transport_halo_type, halo_update, &
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
      logical :: do_vertical   = .false. !< enable vertical (Lagrangian PPM remap) advection

   contains
      ! Required ProcessInterface implementations
      procedure :: init     => transport_init
      procedure :: run      => transport_run
      procedure :: finalize => transport_finalize

      ! Capability registration
      procedure :: get_required_met_fields => transport_get_required_met_fields
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

   end subroutine transport_init

   !> \brief Load transport runtime options from the configuration.
   !!
   !! All keys live under `processes/transport/` and are optional (sensible
   !! defaults are applied when absent, and when no config is attached):
   !!   - `activate`   (logical, .true.) : enable/disable the process entirely.
   !!   - `horizontal` (logical, .true.) : run horizontal advection.
   !!   - `vertical`   (logical, .false.): run the vertical Lagrangian PPM remap
   !!                                      (requires `horizontal`: the remap acts
   !!                                      on the horizontally-deformed layers).
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
         call config_manager%get_logical('processes/transport/vertical',   this%do_vertical,   cfg_rc, .false.)
         call config_manager%get_integer('processes/transport/hord',       this%hord,          cfg_rc, 8)
         call config_manager%get_integer('processes/transport/vord',       this%vord,          cfg_rc, 8)
         call config_manager%get_integer('processes/transport/halo_width', this%ng,            cfg_rc, 3)
         if (this%ng < 1) this%ng = 3
         call config_manager%get_logical('processes/transport/x_periodic', x_periodic,         cfg_rc, .false.)
         call config_manager%get_logical('processes/transport/y_periodic', y_periodic,         cfg_rc, .false.)

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

      type(TimeStateType), pointer :: time_state
      integer :: metrics_rc, nz
      real    :: dt
      real, allocatable :: dp_lag(:,:,:)

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

      if (this%do_vertical .and. allocated(this%met_state%DELP)) then
         ! Capture the Lagrangian layer thicknesses produced by the horizontal
         ! step so they can be remapped back to the reference grid.
         nz = size(this%met_state%DELP, 3)
         allocate(dp_lag(this%bd%is:this%bd%ie, this%bd%js:this%bd%je, nz))
         call transport_horizontal(this, dt, rc, dp_lag=dp_lag)
         if (rc == CC_SUCCESS) call transport_vertical(this, dp_lag, rc)
         deallocate(dp_lag)
      else
         call transport_horizontal(this, dt, rc)
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

      ! --- PJC/LLNL pressure fixer (active only when the host supplies
      !     PS_NEXT and the domain is a single global periodic lon-lat grid) ---
      logical :: do_fix, ok
      real    :: cn, qxw, qys
      real(fp), allocatable :: ap(:), bp(:)
      real, allocatable :: dbk(:)
      real, allocatable :: areaC(:,:), dps_ctm(:,:), xcf(:,:), mmf(:)
      real, allocatable :: cfx(:,:), cfy(:,:), cfxq(:,:), cfyq(:,:)

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
      ! PRESSURE FIXER SET-UP (optional).  When the host supplies the end-of-step
      ! surface pressure (met%PS_NEXT) on a single global periodic lon-lat grid,
      ! solve the PJC/LLNL barotropic correction so the advected (Lagrangian)
      ! surface pressure closes onto PS_NEXT.  The correction is the per-column
      ! (dbk-weighted) mass-flux adjustment that drives FV3's own vertically
      ! integrated divergence to (PS_NEXT - PS).  See fv3_pfix_mod.
      !
      ! NOTE: PS_NEXT is auto-allocated by the MetState field generator, so
      ! `allocated(PS_NEXT)` is not a reliable trigger.  The fixer engages only
      ! when the host has actually populated PS_NEXT this timestep, tracked by
      ! the per-timestep populated-field registry (is_field_set).  When PS_NEXT
      ! is not supplied (e.g. the final met slice has no t+dt pressure), the
      ! fixer is skipped and transport runs its normal conservative path.
      ! -----------------------------------------------------------------------
      do_fix = this%met_state%is_field_set('PS_NEXT') .and. &
               allocated(this%met_state%PS_NEXT) .and. &
               allocated(this%met_state%PS)      .and. &
               allocated(this%met_state%AREA_M2) .and. &
               (this%halo%x_bc == HALO_BC_PERIODIC) .and. &
               (is == 1 .and. js == 1)
      if (do_fix) then
         if (size(this%met_state%PS_NEXT,1) /= nx .or. &
             size(this%met_state%PS_NEXT,2) /= ny .or. &
             size(this%met_state%PS,1) /= nx .or. &
             size(this%met_state%PS,2) /= ny .or. &
             size(this%met_state%AREA_M2,1) /= nx .or. &
             size(this%met_state%AREA_M2,2) /= ny) do_fix = .false.
      end if
      if (do_fix) then
         allocate(dbk(nz))
         call get_hybrid_ab(nz, ap, bp, ok)
         if (.not. ok) then
            do_fix = .false.
            deallocate(dbk)
         else
            do k = 1, nz
               dbk(k) = bp(k) - bp(k+1)
            end do
         end if
      end if
      if (do_fix) then
         allocate(areaC(nx,ny), dps_ctm(nx,ny), xcf(nx,ny), mmf(ny))
         allocate(cfx(is:ie+1, js:je), cfy(is:ie, js:je+1))
         allocate(cfxq(is:ie+1, js:je), cfyq(is:ie, js:je+1))
         areaC(:,:)   = real(this%met_state%AREA_M2(:,:))
         dps_ctm(:,:) = 0.0

         ! Pass 1: accumulate FV3's uncorrected vertically integrated divergence.
         prepass_loop: do k = 1, nz
            ua = 0.0; va = 0.0; delp = 0.0
            ua  (is:ie, js:je) = this%met_state%U   (:,:,k)
            va  (is:ie, js:je) = this%met_state%V   (:,:,k)
            delp(is:ie, js:je) = this%met_state%DELP(:,:,k)
            call halo_update(ua,   this%bd, this%halo, hrc)
            call halo_update(va,   this%bd, this%halo, hrc)
            call halo_update(delp, this%bd, this%halo, hrc)

            call build_level_mass_flux(this%gridstruct, this%bd, ua, va, delp, dt, mf, rc)
            if (rc /= CC_SUCCESS) exit prepass_loop
            cmax  = courant_max(mf, this%bd)
            nsplt = int(1.0 + cmax)
            if (nsplt < 1) nsplt = 1
            if (nsplt > 1) call scale_mass_flux(mf, this%gridstruct, this%bd, 1.0 / real(nsplt))

            do j = js, je
               do i = is, ie
                  dps_ctm(i,j) = dps_ctm(i,j) + real(nsplt) * &
                     (mf%mfx(i,j) - mf%mfx(i+1,j) + mf%mfy(i,j) - mf%mfy(i,j+1)) &
                     * this%gridstruct%rarea(i,j)
               end do
            end do
         end do prepass_loop

         if (rc == CC_SUCCESS) then
            call pfix_correction(nx, ny, areaC, &
                 real(this%met_state%PS), real(this%met_state%PS_NEXT), &
                 dps_ctm, xcf, mmf, cn)
         else
            do_fix = .false.
         end if
      end if

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
         nsplt = int(1.0 + cmax)
         if (nsplt < 1) nsplt = 1
         if (nsplt > 1) call scale_mass_flux(mf, this%gridstruct, this%bd, 1.0 / real(nsplt))

         ! --- Pressure-fixer correction fluxes for this level ---------------
         ! FV3 face fluxes that reproduce the PJC/LLNL barotropic correction:
         !   cfx(i,j) = xcf(i,j) * dbk(k) * area(j)   (west face of cell i)
         !   cfy(i,j) = mmf(j)   * dbk(k) * cn        (south face of cell j)
         ! divided per sub-step so nsplt sub-steps deliver the full correction.
         ! (cfx-cfx(i+1)+cfy-cfy(j+1))*rarea then equals the pfix cell divergence
         ! correction; summed over levels it closes FV3's divergence to PS_NEXT.
         if (do_fix) then
            do j = js, je
               do i = is, ie
                  cfx(i,j) = xcf(i,j) * dbk(k) * areaC(i,j) / real(nsplt)
               end do
               ! East face of the last cell wraps periodically to column 1.
               cfx(ie+1,j) = xcf(1,j) * dbk(k) * areaC(1,j) / real(nsplt)
            end do
            do j = js, je
               do i = is, ie
                  cfy(i,j) = mmf(j) * dbk(k) * cn / real(nsplt)
               end do
            end do
            ! North face of the northern-most cell carries no correction flux.
            cfy(is:ie, je+1) = 0.0
         end if

         ! --- Lagrangian layer thickness (species-independent) ---------------
         ! The flux-form evolution is linear in the (per-substep) mass flux, so
         ! after nsplt substeps every species reaches the same evolved dp2:
         !   dp_lag = delp + nsplt * (mfx-mfx(i+1)+mfy-mfy(j+1)) * rarea,
         ! plus the (nsplt-summed) pressure-fixer correction when active. With
         ! the fixer, sum_k dp_lag = PS_NEXT to machine precision.
         if (present(dp_lag)) then
            do j = js, je
               do i = is, ie
                  dp_lag(i,j,k) = delp(i,j) + real(nsplt) * &
                     (mf%mfx(i,j) - mf%mfx(i+1,j) + mf%mfy(i,j) - mf%mfy(i,j+1)) &
                     * this%gridstruct%rarea(i,j)
               end do
            end do
            if (do_fix) then
               do j = js, je
                  do i = is, ie
                     dp_lag(i,j,k) = dp_lag(i,j,k) + real(nsplt) * &
                        (cfx(i,j) - cfx(i+1,j) + cfy(i,j) - cfy(i,j+1)) &
                        * this%gridstruct%rarea(i,j)
                  end do
               end do
            end if
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

               ! Pressure-fixer tracer flux (upwind on the correction flux).
               ! Added as an extra flux-form term so the corrected mass flux and
               ! tracer flux telescope identically: a uniform field stays
               ! uniform and sum(q*dp) is conserved to machine precision.
               if (do_fix) then
                  do j = js, je
                     do i = is, ie+1
                        qxw = merge(q(i-1,j), q(i,j), cfx(i,j) > 0.0)
                        cfxq(i,j) = qxw * cfx(i,j)
                     end do
                  end do
                  do j = js, je+1
                     do i = is, ie
                        qys = merge(q(i,j-1), q(i,j), cfy(i,j) > 0.0)
                        cfyq(i,j) = qys * cfy(i,j)
                     end do
                  end do

                  do j = js, je
                     do i = is, ie
                        dp2(i,j) = dp1(i,j) + &
                           (mf%mfx(i,j) - mf%mfx(i+1,j) + mf%mfy(i,j) - mf%mfy(i,j+1) &
                          +   cfx(i,j)  -   cfx(i+1,j)  +   cfy(i,j)  -   cfy(i,j+1)) &
                           * this%gridstruct%rarea(i,j)
                        q(i,j) = (q(i,j) * dp1(i,j) + &
                           (fx(i,j)   - fx(i+1,j)   + fy(i,j)   - fy(i,j+1)  &
                          + cfxq(i,j) - cfxq(i+1,j) + cfyq(i,j) - cfyq(i,j+1)) &
                           * this%gridstruct%rarea(i,j)) / dp2(i,j)
                     end do
                  end do
               else
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
               end if

               if (it < nsplt) then
                  call halo_update(q, this%bd, this%halo, hrc)
                  dp1(is:ie, js:je) = dp2(is:ie, js:je)
               end if
            end do subcycle

            this%chem_state%ChemSpecies(isp)%conc(:,:,k) = real(q(is:ie, js:je), fp)
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
      if (allocated(ap))      deallocate(ap)
      if (allocated(bp))      deallocate(bp)
      if (allocated(dbk))     deallocate(dbk)
      if (allocated(areaC))   deallocate(areaC)
      if (allocated(dps_ctm)) deallocate(dps_ctm)
      if (allocated(xcf))     deallocate(xcf)
      if (allocated(mmf))     deallocate(mmf)
      if (allocated(cfx))     deallocate(cfx)
      if (allocated(cfy))     deallocate(cfy)
      if (allocated(cfxq))    deallocate(cfxq)
      if (allocated(cfyq))    deallocate(cfyq)

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
   !! reference hybrid coordinate evaluated at the Lagrangian surface pressure
   !! \f$p_{s}^{lag}=\sum_k \Delta p^{lag}\f$, with the top/surface edges forced to
   !! match the source exactly; a GCHP-style per-column mass fixer then rescales
   !! each column by the ratio of pre- to post-remap tracer mass, so each column's
   !! tracer mass is conserved to machine precision.
   !!
   !! Ordering: CATChem is surface-first (`k=1` is the surface); `mappm` expects
   !! edges from model top to surface (increasing pressure), so columns are flipped
   !! on the way in and out.
   !!
   !! \note The target uses the Lagrangian surface pressure, not the prescribed met
   !!       surface pressure. For offline winds that are not perfectly mass
   !!       consistent with the met pressure field these differ slightly; a
   !!       GCHP-style pressure fixer (future milestone) reconciles them. Tracer
   !!       mass is conserved by construction regardless.
   subroutine transport_vertical(this, dp_lag, rc)
      class(ProcessTransportInterface), intent(inout) :: this
      !> Lagrangian layer thickness [Pa], shape (is:ie, js:je, nz), surface-first.
      real,    intent(in)  :: dp_lag(:,:,:)
      integer, intent(out) :: rc

      integer :: is, ie, js, je, nz, np, i, j, k, c, s, isp
      logical :: ok
      real(fp), allocatable :: ap(:), bp(:)
      real,     allocatable :: pe1(:,:), pe2(:,:), q1(:,:), q2(:,:)
      real                  :: ps_lag, m_src, m_tgt

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

            ! Target edges = reference hybrid coordinate at the Lagrangian surface
            ! pressure. FV3 edge k (top-first) is surface-first edge nz+2-k.
            do k = 1, nz+1
               pe2(c,k) = real(ap(nz+2-k)) + real(bp(nz+2-k)) * ps_lag
            end do
            ! Force exact endpoint match => exact per-column mass conservation.
            pe2(c,1)    = pe1(c,1)
            pe2(c,nz+1) = pe1(c,nz+1)
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

         ! --- Column mass fixer (GCHP-style) -----------------------------------
         ! mappm conserves the mass-weighted integral in the interior, but its
         ! top cell (edge at the model top) takes the top source-cell mean rather
         ! than that cell's sub-cell average, leaving a small residual. Rescale
         ! each positive-definite column by the ratio of pre- to post-remap tracer
         ! mass so the column tracer mass is preserved to machine precision (the
         ! same role as GCHP's pressure/mass fixer).
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
   function transport_get_required_met_fields(this) result(field_names)
      class(ProcessTransportInterface), intent(in) :: this
      character(len=32), allocatable :: field_names(:)

      ! Horizontal advection needs the winds and the layer pressure thickness.
      ! Vertical transport is a Lagrangian PPM remap driven by horizontal mass
      ! convergence (no vertical velocity / OMEGA required).
      allocate(field_names(3))
      field_names(1) = 'U'
      field_names(2) = 'V'
      field_names(3) = 'DELP'

   end function transport_get_required_met_fields

end module ProcessTransportInterface_Mod
