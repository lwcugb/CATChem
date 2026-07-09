!> \file TransportMassFlux_Mod.F90
!! \brief Build FV3 C-grid mass fluxes from CATChem A-grid winds + DELP.
!!
!! The vendored FV3 flux operator (\ref fv3_tp_core_mod, `fv_tp_2d`) does NOT
!! take winds. It takes, on the C-grid, the pre-computed transport quantities:
!!   - `crx`, `cry` : Courant numbers               (dimensionless)
!!   - `xfx`, `yfx` : swept AREA fluxes             [m^2]
!!   - `ra_x`, `ra_y`: 1-D advective-update areas   [m^2]
!!   - `mfx`, `mfy` : layer MASS fluxes             [Pa m^2]
!! This module owns the "prep" that turns the geometry CATChem carries -- cell
!! centre winds `U`,`V` [m s-1] and layer pressure thickness `DELP` [Pa] -- into
!! those quantities, for one horizontal level at a time (fv_tp_2d is called per
!! level, and `DELP`/winds vary with level).
!!
!! ## Algorithm (offline / CTM prep, FV3 `fv_tracer2d` conventions)
!! For a single level, given A-grid `ua`,`va`,`delp` on the DATA domain and the
!! process timestep `dt`:
!!  1. **A -> C restagger** (2-point average):
!!       `uc(i,j) = 0.5*(ua(i-1,j) + ua(i,j))`   at west  faces (i = is..ie+1)
!!       `vc(i,j) = 0.5*(va(i,j-1) + va(i,j))`   at south faces (j = js..je+1)
!!     (Simple a2c; a higher-order A2D2C restagger is a later refinement.)
!!  2. **Courant numbers** with the C-grid cell spacing
!!       `dxc = 0.5*(dxa(i-1,j)+dxa(i,j))`, `dyc = 0.5*(dya(i,j-1)+dya(i,j))`:
!!       `crx = uc*dt/dxc`, `cry = vc*dt/dyc`.
!!  3. **Area fluxes** (upwind-biased A-grid width), as in FV3 `fv_tracer2d`:
!!       `xfx = crx * dxa(upwind) * dy`,  `yfx = cry * dya(upwind) * dx`.
!!     (The FV3 `sin_sg` cube-edge factor is 1 on regular / lat-lon grids and is
!!     omitted here; it is reinstated with the cubed-sphere metrics later.)
!!  4. **Advective areas**   `ra_x = area + xfx(i) - xfx(i+1)`,
!!                           `ra_y = area + yfx(j) - yfx(j+1)`.
!!  5. **Mass fluxes** (upwind `delp`)  `mfx = xfx * delp(upwind)`,
!!                                       `mfy = yfx * delp(upwind)`.
!! These are exactly the arrays `fv_tp_2d` reads; the tracer flux it returns is
!! `0.5*(fx+fx2)*mfx`, and the SAME `mfx`/`mfy` integrate `delp`, so the
!! flux-form tracer update conserves mass by construction.
!!
!! ## Courant sub-cycling
!! `courant_max` returns `cmax = max(|crx|,|cry|)` over the compute domain; the
!! driver forms `nsplt = int(1 + cmax)` and, when `nsplt > 1`, calls
!! `scale_mass_flux` with `frac = 1/nsplt` to divide every flux (and rebuild
!! `ra_x`/`ra_y`). The scaled fluxes are then reused across the `nsplt`
!! sub-steps while `delp` is integrated between them (FV3 `fv_tracer2d`).
!!
!! \note **Dry-air pressure fixer.** A faithful offline scheme (GCHP
!!       `GCHPctmEnv`) additionally corrects the mass fluxes so that advecting
!!       `delp` reproduces the host surface-pressure tendency. That correction
!!       needs the target dp tendency and a halo-wide flux adjustment, so it is
!!       deferred to the halo-exchange / conservation milestone; without it the
!!       scheme is still locally conservative but not pressure-consistent with
!!       the host.
!!
!! \author CATChem Development Team
!! \version 0.1.0
module TransportMassFlux_Mod

   use error_mod, only: CC_SUCCESS, CC_FAILURE
   use fv3_grid_types_mod, only: fv_grid_type, fv_grid_bounds_type

   implicit none
   private

   public :: fv_mass_flux_type
   public :: mass_flux_alloc
   public :: mass_flux_free
   public :: build_level_mass_flux
   public :: courant_max
   public :: scale_mass_flux

   !> \brief C-grid transport quantities for ONE horizontal level.
   !!
   !! Array shapes mirror the `fv_tp_2d` dummy arguments exactly, so a filled
   !! `fv_mass_flux_type` can be passed straight through to the kernel.
   type :: fv_mass_flux_type
      real, allocatable :: crx(:,:)   !< Courant X   (is:ie+1, jsd:jed)
      real, allocatable :: cry(:,:)   !< Courant Y   (isd:ied, js:je+1)
      real, allocatable :: xfx(:,:)   !< area flux X (is:ie+1, jsd:jed) [m^2]
      real, allocatable :: yfx(:,:)   !< area flux Y (isd:ied, js:je+1) [m^2]
      real, allocatable :: ra_x(:,:)  !< adv area X  (is:ie,   jsd:jed) [m^2]
      real, allocatable :: ra_y(:,:)  !< adv area Y  (isd:ied, js:je)   [m^2]
      real, allocatable :: mfx(:,:)   !< mass flux X (is:ie+1, js:je)   [Pa m^2]
      real, allocatable :: mfy(:,:)   !< mass flux Y (is:ie,   js:je+1) [Pa m^2]
      logical :: is_alloc = .false.       !< arrays allocated?
   end type fv_mass_flux_type

contains

   !> \brief Allocate a mass-flux workspace on the FV3 index sub-ranges.
   subroutine mass_flux_alloc(mf, bd, rc)
      type(fv_mass_flux_type),   intent(inout) :: mf
      type(fv_grid_bounds_type), intent(in)    :: bd
      integer,                   intent(out)   :: rc

      integer :: is, ie, js, je, isd, ied, jsd, jed

      rc = CC_SUCCESS
      is = bd%is; ie = bd%ie; js = bd%js; je = bd%je
      isd = bd%isd; ied = bd%ied; jsd = bd%jsd; jed = bd%jed

      call mass_flux_free(mf)

      allocate(mf%crx (is:ie+1, jsd:jed))
      allocate(mf%xfx (is:ie+1, jsd:jed))
      allocate(mf%cry (isd:ied, js:je+1))
      allocate(mf%yfx (isd:ied, js:je+1))
      allocate(mf%ra_x(is:ie,   jsd:jed))
      allocate(mf%ra_y(isd:ied, js:je))
      allocate(mf%mfx (is:ie+1, js:je))
      allocate(mf%mfy (is:ie,   js:je+1))

      mf%crx = 0.0; mf%xfx = 0.0
      mf%cry = 0.0; mf%yfx = 0.0
      mf%ra_x = 0.0; mf%ra_y = 0.0
      mf%mfx = 0.0; mf%mfy = 0.0
      mf%is_alloc = .true.

   end subroutine mass_flux_alloc

   !> \brief Release a mass-flux workspace.
   subroutine mass_flux_free(mf)
      type(fv_mass_flux_type), intent(inout) :: mf

      if (allocated(mf%crx))  deallocate(mf%crx)
      if (allocated(mf%xfx))  deallocate(mf%xfx)
      if (allocated(mf%cry))  deallocate(mf%cry)
      if (allocated(mf%yfx))  deallocate(mf%yfx)
      if (allocated(mf%ra_x)) deallocate(mf%ra_x)
      if (allocated(mf%ra_y)) deallocate(mf%ra_y)
      if (allocated(mf%mfx))  deallocate(mf%mfx)
      if (allocated(mf%mfy))  deallocate(mf%mfy)
      mf%is_alloc = .false.

   end subroutine mass_flux_free

   !> \brief Build the C-grid transport quantities for one level.
   !!
   !! \param[in]    gridstruct  FV3 grid metrics (area, dx, dy, dxa, dya).
   !! \param[in]    bd          FV3 index bounds.
   !! \param[in]    ua,va       A-grid (cell-centre) winds [m s-1], DATA domain.
   !! \param[in]    delp        layer pressure thickness [Pa], DATA domain.
   !! \param[in]    dt          transport timestep [s].
   !! \param[inout] mf          workspace to fill (must be alloc'd for this bd).
   !! \param[out]   rc          CC_SUCCESS / CC_FAILURE.
   !!
   !! `ua`,`va`,`delp` are expected on the DATA domain (halo-filled); the halo is
   !! supplied by the exchange stage. `crx` needs `ua(is-1)`, `cry` needs
   !! `va(js-1)`, both inside the data domain for ng >= 1.
   subroutine build_level_mass_flux(gridstruct, bd, ua, va, delp, dt, mf, rc)
      type(fv_grid_type),        intent(in)    :: gridstruct
      type(fv_grid_bounds_type), intent(in)    :: bd
      real, intent(in) :: ua  (bd%isd:bd%ied, bd%jsd:bd%jed)
      real, intent(in) :: va  (bd%isd:bd%ied, bd%jsd:bd%jed)
      real, intent(in) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed)
      real,                      intent(in)    :: dt
      type(fv_mass_flux_type),   intent(inout) :: mf
      integer,                   intent(out)   :: rc

      integer  :: is, ie, js, je, isd, ied, jsd, jed
      integer  :: i, j
      real     :: uc, vc, dxc, dyc

      rc = CC_SUCCESS
      if (.not. mf%is_alloc) then
         rc = CC_FAILURE
         return
      end if

      is = bd%is; ie = bd%ie; js = bd%js; je = bd%je
      isd = bd%isd; ied = bd%ied; jsd = bd%jsd; jed = bd%jed

      ! --- X faces: A->C restagger, Courant number, swept area flux ----------
      do j = jsd, jed
         do i = is, ie+1
            uc  = 0.5 * (ua(i-1,j) + ua(i,j))
            dxc = 0.5 * (gridstruct%dxa(i-1,j) + gridstruct%dxa(i,j))
            mf%crx(i,j) = uc * dt / dxc
            if (mf%crx(i,j) > 0.0) then
               mf%xfx(i,j) = mf%crx(i,j) * gridstruct%dxa(i-1,j) * gridstruct%dy(i,j)
            else
               mf%xfx(i,j) = mf%crx(i,j) * gridstruct%dxa(i,  j) * gridstruct%dy(i,j)
            end if
         end do
      end do

      ! --- Y faces: A->C restagger, Courant number, swept area flux ----------
      do j = js, je+1
         do i = isd, ied
            vc  = 0.5 * (va(i,j-1) + va(i,j))
            dyc = 0.5 * (gridstruct%dya(i,j-1) + gridstruct%dya(i,j))
            mf%cry(i,j) = vc * dt / dyc
            if (mf%cry(i,j) > 0.0) then
               mf%yfx(i,j) = mf%cry(i,j) * gridstruct%dya(i,j-1) * gridstruct%dx(i,j)
            else
               mf%yfx(i,j) = mf%cry(i,j) * gridstruct%dya(i,j  ) * gridstruct%dx(i,j)
            end if
         end do
      end do

      ! --- Advective-update areas -------------------------------------------
      call compute_advective_areas(gridstruct, bd, mf)

      ! --- Layer mass fluxes (upwind delp) ----------------------------------
      do j = js, je
         do i = is, ie+1
            if (mf%crx(i,j) > 0.0) then
               mf%mfx(i,j) = mf%xfx(i,j) * delp(i-1,j)
            else
               mf%mfx(i,j) = mf%xfx(i,j) * delp(i,j)
            end if
         end do
      end do
      do j = js, je+1
         do i = is, ie
            if (mf%cry(i,j) > 0.0) then
               mf%mfy(i,j) = mf%yfx(i,j) * delp(i,j-1)
            else
               mf%mfy(i,j) = mf%yfx(i,j) * delp(i,j)
            end if
         end do
      end do

   end subroutine build_level_mass_flux

   !> \brief (Re)compute ra_x / ra_y from the current area fluxes.
   subroutine compute_advective_areas(gridstruct, bd, mf)
      type(fv_grid_type),        intent(in)    :: gridstruct
      type(fv_grid_bounds_type), intent(in)    :: bd
      type(fv_mass_flux_type),   intent(inout) :: mf

      integer :: is, ie, js, je, isd, ied, jsd, jed
      integer :: i, j

      is = bd%is; ie = bd%ie; js = bd%js; je = bd%je
      isd = bd%isd; ied = bd%ied; jsd = bd%jsd; jed = bd%jed

      do j = jsd, jed
         do i = is, ie
            mf%ra_x(i,j) = gridstruct%area(i,j) + mf%xfx(i,j) - mf%xfx(i+1,j)
         end do
      end do
      do j = js, je
         do i = isd, ied
            mf%ra_y(i,j) = gridstruct%area(i,j) + mf%yfx(i,j) - mf%yfx(i,j+1)
         end do
      end do

   end subroutine compute_advective_areas

   !> \brief Maximum Courant number over the compute domain.
   !!
   !! `cmax = max(|crx|, |cry|)` on the u/v faces of the compute domain, used to
   !! size the sub-cycling: `nsplt = int(1 + cmax)`.
   function courant_max(mf, bd) result(cmax)
      type(fv_mass_flux_type),   intent(in) :: mf
      type(fv_grid_bounds_type), intent(in) :: bd
      real :: cmax

      integer :: is, ie, js, je
      integer :: i, j

      is = bd%is; ie = bd%ie; js = bd%js; je = bd%je
      cmax = 0.0

      do j = js, je
         do i = is, ie+1
            cmax = max(cmax, abs(mf%crx(i,j)))
         end do
      end do
      do j = js, je+1
         do i = is, ie
            cmax = max(cmax, abs(mf%cry(i,j)))
         end do
      end do

   end function courant_max

   !> \brief Scale every flux by `frac` and rebuild the advective areas.
   !!
   !! Used for Courant sub-cycling: with `frac = 1/nsplt` the scaled `crx`,`cry`,
   !! `xfx`,`yfx`,`mfx`,`mfy` describe one sub-step, and `ra_x`,`ra_y` are made
   !! consistent with the scaled `xfx`,`yfx`.
   subroutine scale_mass_flux(mf, gridstruct, bd, frac)
      type(fv_mass_flux_type),   intent(inout) :: mf
      type(fv_grid_type),        intent(in)    :: gridstruct
      type(fv_grid_bounds_type), intent(in)    :: bd
      real,                      intent(in)    :: frac

      mf%crx = mf%crx * frac
      mf%cry = mf%cry * frac
      mf%xfx = mf%xfx * frac
      mf%yfx = mf%yfx * frac
      mf%mfx = mf%mfx * frac
      mf%mfy = mf%mfy * frac

      call compute_advective_areas(gridstruct, bd, mf)

   end subroutine scale_mass_flux

end module TransportMassFlux_Mod
