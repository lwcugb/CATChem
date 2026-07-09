!> \file TransportGridMetrics_Mod.F90
!! \brief Build FV3 tp_core grid metrics from CATChem grid geometry.
!!
!! The vendored FV3 flux operator (\ref fv3_tp_core_mod) consumes a
!! \ref fv_grid_type (cell areas + along-axis spacings) and a
!! \ref fv_grid_bounds_type (compute/data index bounds). This module builds
!! those structures from the geometry CATChem actually carries at run time:
!! per-cell centre latitude/longitude and (optionally) true cell areas.
!!
!! ## Grid-type flexibility
!! All lengths are computed as GREAT-CIRCLE distances between grid CELL CORNERS.
!! This single formulation is correct for every grid CATChem may run on:
!!   - **Rectilinear** regular lat/lon (uniform or stretched).
!!   - **Curvilinear** / rotated / conformal grids.
!!   - **Cubed-sphere** tiles (gnomonic, per face).
!! The only per-grid difference is how the corners and the edge/corner flags are
!! supplied:
!!   - If the host provides staggered CORNER coordinates (curvilinear / cubed
!!     sphere), pass them via `lat_cor_deg`/`lon_cor_deg` for exact metrics.
!!   - Otherwise corners are reconstructed from the cell centres in 3-D
!!     (unit-vector) space with one-cell edge extrapolation. This is exact for a
!!     regular lat/lon grid and accurate for smooth curvilinear grids.
!!
!! ## Metric definitions (FV3 conventions)
!! These reproduce GFDL `tools/fv_grid_tools.F90` exactly (great-circle edges +
!! edge-midpoint widths); only the branching differs (see below).
!!   - `area(i,j)`  : cell area [m^2] (from AREA_M2 if given, else spherical
!!                    excess of the four corners; cf. FV3 `get_area`).
!!   - `rarea`      : 1/area.
!!   - `dx(i,j)`    : cell south edge = gc(p(i,j), p(i+1,j)).
!!                    (FV3 fv_grid_tools.F90: dx = great_circle_dist(grid(i,j),grid(i+1,j)))
!!   - `dy(i,j)`    : cell west  edge = gc(p(i,j), p(i,j+1)).
!!   - `dxa(i,j)`   : A-grid E-W width = gc(mid(west edge), mid(east edge)).
!!                    (FV3: mid_pt_sphere of the two y-edges, then great_circle_dist)
!!   - `dya(i,j)`   : A-grid N-S width = gc(mid(south edge), mid(north edge)).
!! The helpers mirror FV3: `midpoint` == `mid_pt_sphere` (latlon2xyz -> Cartesian
!! midpoint -> normalise); `gc_len` (atan2 of |cross|,dot) == `great_circle_dist`.
!! Inside `fv_tp_2d` only `area`, `dxa`, `dya` and the flags are read; `dx`,
!! `dy`, `rarea` feed the downstream mass-flux / flux-form update stages.
!!
!! ## Relation to FV3 / UFS (fv3atm) and GCHP
!!   - FV3 hard-codes a separate branch per grid_type: analytic
!!     `dl*radius*cos(lat)` for a pure lat/lon grid, and geodesic-from-corners
!!     for the cubed sphere. We use the SINGLE geodesic-from-corners form for
!!     every grid (the cubed-sphere method), which is negligibly different
!!     (O(dlon^2)) from the analytic values on a regular lat/lon grid.
!!   - GCHP/MAPL do not recompute these; they take the cubed-sphere grid's AREA
!!     field and feed the same tp_core. We mirror that by preferring the host's
!!     true `AREA_M2` when available, falling back to the FV3 corner area.
!!
!! \note Halo cells of the (static) metric arrays are currently filled by edge
!!       replication (exact in the x-halo for same-latitude neighbours; only
!!       approximate in the y-halo near the poles). FV3 instead halo-fills the
!!       corner array before computing metrics; that refinement, together with a
!!       global `da_min` reduction, belongs to the halo-exchange stage.
!!
!! \author CATChem Development Team
!! \version 0.2.0
module TransportGridMetrics_Mod

   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE
   use Constants, only: Re, PI_180

   use fv3_grid_types_mod, only: fv_grid_type, fv_grid_bounds_type
   use TransportHalo_Mod, only: halo_fill_scalar, HALO_BC_REPLICATE, HALO_BC_PERIODIC

   implicit none
   private

   public :: build_fv3_grid_metrics

   !> FV3 grid_type id for a regular (uniform / doubly-periodic) grid.
   !! Any value >= 3 selects the simple-branch stencils in xppm/yppm and
   !! disables the cubed-sphere edge treatment. Use 0 for a gnomonic
   !! cubed-sphere face (together with the appropriate corner flags).
   integer, parameter, public :: FV3_GRID_REGULAR = 4

contains

   !> \brief Populate an FV3 gridstruct + bounds from CATChem grid geometry.
   !!
   !! \param[in]  lat_deg     cell-centre latitude  [deg], shape (nx,ny)
   !! \param[in]  lon_deg     cell-centre longitude [deg], shape (nx,ny)
   !! \param[in]  ng          halo (ghost) width (>= 3 for PPM / hord=8)
   !! \param[inout] gridstruct FV3 grid metrics to (re)allocate + fill
   !! \param[out] bd          FV3 index bounds (compute + data domains)
   !! \param[out] rc          return code (CC_SUCCESS / CC_FAILURE)
   !! \param[in]  area_m2     OPTIONAL true cell areas [m^2], shape (nx,ny).
   !!                         If absent, areas are computed from the corners.
   !! \param[in]  lat_cor_deg OPTIONAL corner latitude  [deg], shape (nx+1,ny+1)
   !! \param[in]  lon_cor_deg OPTIONAL corner longitude [deg], shape (nx+1,ny+1)
   !! \param[in]  grid_type   OPTIONAL FV3 grid_type id (default FV3_GRID_REGULAR)
   !! \param[in]  bounded_domain OPTIONAL (default .false.)
   !! \param[in]  x_periodic  OPTIONAL (default .false.). When .true. the x-halo
   !!                         of every metric is filled by a periodic WRAP rather
   !!                         than edge replication. This is required for exact
   !!                         mass conservation on a zonally periodic lon/lat
   !!                         grid: the seam-edge lengths/areas (reconstructed
   !!                         from one-sided extrapolated corners) then match
   !!                         bit-for-bit across the west/east seam, so the
   !!                         mass-flux divergence telescopes to zero.
   !! \param[in]  sw_corner,se_corner,nw_corner,ne_corner OPTIONAL cube-edge
   !!                         flags (default .false.); only meaningful when
   !!                         bounded_domain is .false. (cubed sphere).
   subroutine build_fv3_grid_metrics(lat_deg, lon_deg, ng, gridstruct, bd, rc, &
                                     area_m2, lat_cor_deg, lon_cor_deg, &
                                     grid_type, bounded_domain, x_periodic, &
                                     sw_corner, se_corner, nw_corner, ne_corner)
      real(fp),                  intent(in)    :: lat_deg(:,:)
      real(fp),                  intent(in)    :: lon_deg(:,:)
      integer,                   intent(in)    :: ng
      type(fv_grid_type),        intent(inout) :: gridstruct
      type(fv_grid_bounds_type), intent(out)   :: bd
      integer,                   intent(out)   :: rc
      real(fp), optional,        intent(in)    :: area_m2(:,:)
      real(fp), optional,        intent(in)    :: lat_cor_deg(:,:)
      real(fp), optional,        intent(in)    :: lon_cor_deg(:,:)
      integer,  optional,        intent(in)    :: grid_type
      logical,  optional,        intent(in)    :: bounded_domain
      logical,  optional,        intent(in)    :: x_periodic
      logical,  optional,        intent(in)    :: sw_corner, se_corner, nw_corner, ne_corner

      integer :: nx, ny
      integer :: is, ie, js, je, isd, ied, jsd, jed
      integer :: i, j
      integer :: metric_x_bc
      real :: amin
      logical  :: have_corners
      ! Corner unit vectors, shape (nx+1, ny+1)
      real, allocatable :: px(:,:), py(:,:), pz(:,:)
      real :: wmx, wmy, wmz, emx, emy, emz   ! west/east edge midpoints
      real :: smx, smy, smz, nmx, nmy, nmz   ! south/north edge midpoints

      rc = CC_SUCCESS

      nx = size(lat_deg, 1)
      ny = size(lat_deg, 2)

      ! Transport needs a genuine 2-D horizontal grid.
      if (nx < 2 .or. ny < 2) then
         rc = CC_FAILURE
         return
      end if
      if (size(lon_deg,1) /= nx .or. size(lon_deg,2) /= ny) then
         rc = CC_FAILURE
         return
      end if
      if (ng < 3) then   ! PPM (hord=8) needs a 3-cell stencil halo
         rc = CC_FAILURE
         return
      end if
      if (present(area_m2)) then
         if (size(area_m2,1) /= nx .or. size(area_m2,2) /= ny) then
            rc = CC_FAILURE
            return
         end if
      end if

      have_corners = present(lat_cor_deg) .and. present(lon_cor_deg)
      if (have_corners) then
         if (size(lat_cor_deg,1) /= nx+1 .or. size(lat_cor_deg,2) /= ny+1 .or. &
             size(lon_cor_deg,1) /= nx+1 .or. size(lon_cor_deg,2) /= ny+1) then
            rc = CC_FAILURE
            return
         end if
      end if

      ! --- index bounds ----------------------------------------------------
      is  = 1;      ie  = nx
      js  = 1;      je  = ny
      isd = 1 - ng; ied = nx + ng
      jsd = 1 - ng; jed = ny + ng

      bd%is  = is;  bd%ie  = ie
      bd%js  = js;  bd%je  = je
      bd%isd = isd; bd%ied = ied
      bd%jsd = jsd; bd%jed = jed
      bd%ng  = ng

      ! --- (re)allocate metric arrays over the data domain -----------------
      call realloc_2d(gridstruct%area,  isd, ied, jsd, jed)
      call realloc_2d(gridstruct%rarea, isd, ied, jsd, jed)
      call realloc_2d(gridstruct%dx,    isd, ied, jsd, jed)
      call realloc_2d(gridstruct%dy,    isd, ied, jsd, jed)
      call realloc_2d(gridstruct%dxa,   isd, ied, jsd, jed)
      call realloc_2d(gridstruct%dya,   isd, ied, jsd, jed)

      ! --- grid corner unit vectors ----------------------------------------
      allocate(px(nx+1, ny+1), py(nx+1, ny+1), pz(nx+1, ny+1))
      if (have_corners) then
         call corners_from_input(lat_cor_deg, lon_cor_deg, px, py, pz)
      else
         call corners_from_centers(lat_deg, lon_deg, px, py, pz)
      end if

      ! --- edge lengths, cell widths, area over the compute domain ---------
      do j = js, je
         do i = is, ie
            ! Edge lengths (south edge, west edge).
            gridstruct%dx(i,j) = gc_len(px(i,j),   py(i,j),   pz(i,j), &
                                        px(i+1,j),  py(i+1,j),  pz(i+1,j))
            gridstruct%dy(i,j) = gc_len(px(i,j),   py(i,j),   pz(i,j), &
                                        px(i,j+1),  py(i,j+1),  pz(i,j+1))

            ! A-grid widths from edge midpoints.
            call midpoint(px(i,j),   py(i,j),   pz(i,j), &
                          px(i,j+1), py(i,j+1), pz(i,j+1), wmx, wmy, wmz)
            call midpoint(px(i+1,j),   py(i+1,j),   pz(i+1,j), &
                          px(i+1,j+1), py(i+1,j+1), pz(i+1,j+1), emx, emy, emz)
            gridstruct%dxa(i,j) = gc_len(wmx, wmy, wmz, emx, emy, emz)

            call midpoint(px(i,j),   py(i,j),   pz(i,j), &
                          px(i+1,j), py(i+1,j), pz(i+1,j), smx, smy, smz)
            call midpoint(px(i,j+1),   py(i,j+1),   pz(i,j+1), &
                          px(i+1,j+1), py(i+1,j+1), pz(i+1,j+1), nmx, nmy, nmz)
            gridstruct%dya(i,j) = gc_len(smx, smy, smz, nmx, nmy, nmz)

            ! Cell area: prefer the host's true area, else spherical excess.
            if (present(area_m2)) then
               gridstruct%area(i,j) = area_m2(i,j)
            else
               gridstruct%area(i,j) = quad_area( &
                  px(i,j),     py(i,j),     pz(i,j),     &   ! SW
                  px(i+1,j),   py(i+1,j),   pz(i+1,j),   &   ! SE
                  px(i+1,j+1), py(i+1,j+1), pz(i+1,j+1), &   ! NE
                  px(i,j+1),   py(i,j+1),   pz(i,j+1))       ! NW
            end if
         end do
      end do

      ! --- halo fill -------------------------------------------------------
      ! X: periodic WRAP for a zonally periodic grid (so the extrapolated
      ! seam-cell metrics match bit-for-bit across the west/east seam and mass
      ! conserves); edge REPLICATE otherwise. Y: always REPLICATE (polar fold
      ! is deferred).
      metric_x_bc = HALO_BC_REPLICATE
      if (present(x_periodic)) then
         if (x_periodic) metric_x_bc = HALO_BC_PERIODIC
      end if
      call fill_halo(gridstruct%area, is, ie, js, je, isd, ied, jsd, jed, metric_x_bc)
      call fill_halo(gridstruct%dx,   is, ie, js, je, isd, ied, jsd, jed, metric_x_bc)
      call fill_halo(gridstruct%dy,   is, ie, js, je, isd, ied, jsd, jed, metric_x_bc)
      call fill_halo(gridstruct%dxa,  is, ie, js, je, isd, ied, jsd, jed, metric_x_bc)
      call fill_halo(gridstruct%dya,  is, ie, js, je, isd, ied, jsd, jed, metric_x_bc)

      ! --- reciprocal area + minimum area ----------------------------------
      do j = jsd, jed
         do i = isd, ied
            if (gridstruct%area(i,j) > 0.0_fp) then
               gridstruct%rarea(i,j) = 1.0_fp / gridstruct%area(i,j)
            else
               gridstruct%rarea(i,j) = 0.0_fp
            end if
         end do
      end do
      amin = huge(1.0_fp)
      do j = js, je
         do i = is, ie
            amin = min(amin, gridstruct%area(i,j))
         end do
      end do
      gridstruct%da_min = amin

      ! --- flags -----------------------------------------------------------
      if (present(grid_type)) then
         gridstruct%grid_type = grid_type
      else
         gridstruct%grid_type = FV3_GRID_REGULAR
      end if
      if (present(bounded_domain)) then
         gridstruct%bounded_domain = bounded_domain
      else
         gridstruct%bounded_domain = .false.
      end if
      gridstruct%sw_corner = merge_flag(sw_corner)
      gridstruct%se_corner = merge_flag(se_corner)
      gridstruct%nw_corner = merge_flag(nw_corner)
      gridstruct%ne_corner = merge_flag(ne_corner)

      ! del6_u/del6_v and the USE_SG metrics are intentionally left
      ! unallocated: deln_flux (hyperdiffusion) and super-grid support are not
      ! active for regular-grid tracer advection.

      deallocate(px, py, pz)

   end subroutine build_fv3_grid_metrics

   ! =====================================================================
   ! Corner construction
   ! =====================================================================

   !> \brief Convert supplied corner lat/lon [deg] to unit vectors.
   subroutine corners_from_input(lat_cor_deg, lon_cor_deg, px, py, pz)
      real(fp), intent(in)  :: lat_cor_deg(:,:), lon_cor_deg(:,:)
      real,     intent(out) :: px(:,:), py(:,:), pz(:,:)
      integer :: i, j, ncx, ncy

      ncx = size(px, 1)
      ncy = size(px, 2)
      do j = 1, ncy
         do i = 1, ncx
            call ll2vec(lat_cor_deg(i,j), lon_cor_deg(i,j), px(i,j), py(i,j), pz(i,j))
         end do
      end do
   end subroutine corners_from_input

   !> \brief Reconstruct corner unit vectors from cell-centre lat/lon.
   !!
   !! Centres are converted to unit vectors, extended by a one-cell ghost ring
   !! via linear extrapolation in 3-D (then renormalised), and each corner is
   !! the normalised average of its four surrounding (extended) centres. Exact
   !! for a regular lat/lon grid; accurate for smooth curvilinear grids.
   subroutine corners_from_centers(lat_deg, lon_deg, px, py, pz)
      real(fp), intent(in)  :: lat_deg(:,:), lon_deg(:,:)
      real,     intent(out) :: px(:,:), py(:,:), pz(:,:)

      integer :: nx, ny, i, j
      real, allocatable :: cx(:,:), cy(:,:), cz(:,:)   ! (0:nx+1, 0:ny+1)

      nx = size(lat_deg, 1)
      ny = size(lat_deg, 2)
      allocate(cx(0:nx+1, 0:ny+1), cy(0:nx+1, 0:ny+1), cz(0:nx+1, 0:ny+1))

      ! Interior centres.
      do j = 1, ny
         do i = 1, nx
            call ll2vec(lat_deg(i,j), lon_deg(i,j), cx(i,j), cy(i,j), cz(i,j))
         end do
      end do

      ! Extrapolate the west/east ghost columns (interior rows).
      do j = 1, ny
         call extrap(cx(1,j), cy(1,j), cz(1,j), cx(2,j), cy(2,j), cz(2,j), &
                     cx(0,j), cy(0,j), cz(0,j))
         call extrap(cx(nx,j), cy(nx,j), cz(nx,j), cx(nx-1,j), cy(nx-1,j), cz(nx-1,j), &
                     cx(nx+1,j), cy(nx+1,j), cz(nx+1,j))
      end do

      ! Extrapolate the south/north ghost rows (all columns, incl. i-ghosts).
      do i = 0, nx+1
         call extrap(cx(i,1), cy(i,1), cz(i,1), cx(i,2), cy(i,2), cz(i,2), &
                     cx(i,0), cy(i,0), cz(i,0))
         call extrap(cx(i,ny), cy(i,ny), cz(i,ny), cx(i,ny-1), cy(i,ny-1), cz(i,ny-1), &
                     cx(i,ny+1), cy(i,ny+1), cz(i,ny+1))
      end do

      ! Corner = normalised average of the four surrounding (extended) centres.
      do j = 1, ny+1
         do i = 1, nx+1
            call normalize(cx(i-1,j-1) + cx(i,j-1) + cx(i-1,j) + cx(i,j), &
                           cy(i-1,j-1) + cy(i,j-1) + cy(i-1,j) + cy(i,j), &
                           cz(i-1,j-1) + cz(i,j-1) + cz(i-1,j) + cz(i,j), &
                           px(i,j), py(i,j), pz(i,j))
         end do
      end do

      deallocate(cx, cy, cz)
   end subroutine corners_from_centers

   ! =====================================================================
   ! Small spherical-geometry helpers (all on the unit sphere)
   ! =====================================================================

   !> \brief lat/lon [deg] -> unit vector.
   pure subroutine ll2vec(lat_deg, lon_deg, x, y, z)
      real(fp), intent(in)  :: lat_deg, lon_deg
      real,     intent(out) :: x, y, z
      real :: la, lo

      la = lat_deg * PI_180
      lo = lon_deg * PI_180
      x = cos(la) * cos(lo)
      y = cos(la) * sin(lo)
      z = sin(la)
   end subroutine ll2vec

   !> \brief Normalise a 3-vector (returns the zero vector if degenerate).
   pure subroutine normalize(x, y, z, ux, uy, uz)
      real, intent(in)  :: x, y, z
      real, intent(out) :: ux, uy, uz
      real :: mag

      mag = sqrt(x*x + y*y + z*z)
      if (mag > 0.0_fp) then
         ux = x / mag; uy = y / mag; uz = z / mag
      else
         ux = 0.0_fp; uy = 0.0_fp; uz = 0.0_fp
      end if
   end subroutine normalize

   !> \brief Linear extrapolation of a ghost point: g = normalize(2a - b).
   pure subroutine extrap(ax, ay, az, bx, by, bz, gx, gy, gz)
      real, intent(in)  :: ax, ay, az, bx, by, bz
      real, intent(out) :: gx, gy, gz

      call normalize(2.0_fp*ax - bx, 2.0_fp*ay - by, 2.0_fp*az - bz, gx, gy, gz)
   end subroutine extrap

   !> \brief Normalised midpoint of two unit vectors (point on the arc).
   pure subroutine midpoint(ax, ay, az, bx, by, bz, mx, my, mz)
      real, intent(in)  :: ax, ay, az, bx, by, bz
      real, intent(out) :: mx, my, mz

      call normalize(ax + bx, ay + by, az + bz, mx, my, mz)
   end subroutine midpoint

   !> \brief Great-circle distance [m] between two unit vectors.
   pure function gc_len(ax, ay, az, bx, by, bz) result(d)
      real, intent(in) :: ax, ay, az, bx, by, bz
      real :: d
      real :: cx, cy, cz, crossmag, dotp

      cx = ay*bz - az*by
      cy = az*bx - ax*bz
      cz = ax*by - ay*bx
      crossmag = sqrt(cx*cx + cy*cy + cz*cz)
      dotp = ax*bx + ay*by + az*bz
      d = Re * atan2(crossmag, dotp)
   end function gc_len

   !> \brief Area [m^2] of a spherical triangle (unit-vector vertices).
   !! Uses the Van Oosterom & Strackee formulation of the spherical excess.
   pure function tri_area(ax, ay, az, bx, by, bz, cx, cy, cz) result(a)
      real, intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz
      real :: a
      real :: num, den, ab, bc, ca

      ! num = | a . (b x c) |
      num = abs(ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx))
      ab = ax*bx + ay*by + az*bz
      bc = bx*cx + by*cy + bz*cz
      ca = cx*ax + cy*ay + cz*az
      den = 1.0_fp + ab + bc + ca
      a = Re * Re * 2.0_fp * atan2(num, den)
   end function tri_area

   !> \brief Area [m^2] of a spherical quad (vertices in order, unit vectors).
   pure function quad_area(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz) result(a)
      real, intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz
      real :: a

      a = tri_area(ax, ay, az, bx, by, bz, cx, cy, cz) &
        + tri_area(ax, ay, az, cx, cy, cz, dx, dy, dz)
   end function quad_area

   ! =====================================================================
   ! Array utilities
   ! =====================================================================

   !> \brief Resolve an optional logical flag (default .false.).
   pure function merge_flag(flag) result(v)
      logical, optional, intent(in) :: flag
      logical :: v

      if (present(flag)) then
         v = flag
      else
         v = .false.
      end if
   end function merge_flag

   !> \brief (Re)allocate a 2-D array with the given lower/upper bounds.
   subroutine realloc_2d(arr, lo1, hi1, lo2, hi2)
      real, allocatable, intent(inout) :: arr(:,:)
      integer, intent(in) :: lo1, hi1, lo2, hi2

      if (allocated(arr)) deallocate(arr)
      allocate(arr(lo1:hi1, lo2:hi2))
      arr = 0.0
   end subroutine realloc_2d

   !> \brief Fill the halo ring of a static metric array.
   !!
   !! Delegates to the shared halo routine (\ref TransportHalo_Mod) so field and
   !! metric halos share one implementation. The x-axis condition is caller
   !! supplied: on a zonally periodic lon/lat grid it MUST be a periodic wrap so
   !! the seam-cell metrics (built from one-sided extrapolated corners) match
   !! bit-for-bit across the west/east seam; edge replication there would leave a
   !! seam mismatch that breaks tracer-mass conservation. The y-axis is always
   !! replicated (the polar fold is deferred).
   subroutine fill_halo(arr, is, ie, js, je, isd, ied, jsd, jed, x_bc)
      real, allocatable, intent(inout) :: arr(:,:)
      integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed
      integer, intent(in) :: x_bc

      call halo_fill_scalar(arr, is, ie, js, je, isd, ied, jsd, jed, &
                            x_bc, HALO_BC_REPLICATE)
   end subroutine fill_halo

end module TransportGridMetrics_Mod
