!> \file fv3_grid_types_mod.F90
!! \brief Trimmed grid-metric derived types for the vendored FV3 tp_core kernel.
!!
!! These types replace `fv_grid_type` and `fv_grid_bounds_type` from the GFDL
!! FV3 dynamical core (`fv_arrays_mod`). They expose ONLY the members that the
!! vendored `fv3_tp_core_mod` (tp_core.F90) actually references, so the kernel
!! can be built in-tree WITHOUT any FMS / FV3 framework dependency.
!!
!! The member names and ranks are kept byte-for-byte compatible with the FV3
!! originals so that `fv3_tp_core_mod` needs no edits to its body — only its
!! `use` line points here instead of `fv_arrays_mod`.
!!
!! Metric arrays are populated by TransportGridMetrics_Mod from the CATChem
!! GridManager / ESMF grid. Members guarded by `USE_SG` in the kernel
!! (`rdxc`, `rdyc`, `sin_sg`) are declared here too so that cubed-sphere /
!! super-grid support can be enabled later without touching this module.
!!
!! \note Bare `real` is used to match tp_core.F90 exactly; the transport library
!!       is compiled with `-r8` / `-fdefault-real-8`, so `real` == `real(fp)`.
!!
!! \author CATChem Development Team
module fv3_grid_types_mod

   implicit none
   private

   public :: fv_grid_type
   public :: fv_grid_bounds_type

   !> \brief Local index bounds for a decomposition tile (replaces FV3 fv_grid_bounds_type).
   !!
   !! Compute domain is [is:ie, js:je]; data (halo-inclusive) domain is
   !! [isd:ied, jsd:jed]; ng is the halo (ghost) width.
   type :: fv_grid_bounds_type
      integer :: is  = 1   !< first compute-domain index, i
      integer :: ie  = 1   !< last  compute-domain index, i
      integer :: js  = 1   !< first compute-domain index, j
      integer :: je  = 1   !< last  compute-domain index, j
      integer :: isd = 1   !< first data-domain (halo) index, i
      integer :: ied = 1   !< last  data-domain (halo) index, i
      integer :: jsd = 1   !< first data-domain (halo) index, j
      integer :: jed = 1   !< last  data-domain (halo) index, j
      integer :: ng  = 3   !< halo (ghost) width
   end type fv_grid_bounds_type

   !> \brief Grid metrics used by the FV3 tp_core transport kernel
   !!        (trimmed replacement for FV3 fv_grid_type).
   type :: fv_grid_type
      ! --- cell-area metrics (data-domain sized: isd:ied, jsd:jed) ---
      real, allocatable :: area(:,:)    !< cell area [m^2]
      real, allocatable :: rarea(:,:)   !< 1 / area

      ! --- along-axis grid spacings (data-domain sized) ---
      real, allocatable :: dxa(:,:)     !< A-grid cell width  in x [m]
      real, allocatable :: dya(:,:)     !< A-grid cell height in y [m]
      real, allocatable :: dx(:,:)      !< distance between corners along x [m]
      real, allocatable :: dy(:,:)      !< distance between corners along y [m]

      ! --- del-n damping coefficients (used by deln_flux, non-USE_SG path) ---
      real, allocatable :: del6_u(:,:)  !< del-6 damping metric, u faces
      real, allocatable :: del6_v(:,:)  !< del-6 damping metric, v faces

      ! --- USE_SG (super-grid) metrics: unused unless -DUSE_SG is set ---
      real, allocatable :: rdxc(:,:)      !< 1 / dxc  (USE_SG only)
      real, allocatable :: rdyc(:,:)      !< 1 / dyc  (USE_SG only)
      real, allocatable :: sin_sg(:,:,:)  !< super-grid sine terms (USE_SG only)

      ! --- scalars / flags ---
      real    :: da_min = 0.0            !< minimum cell area (for damping scaling)
      integer :: grid_type = 0           !< 0-3 cubed-sphere variants; >=3 => regular/doubly-periodic
      logical :: bounded_domain = .false.!< .true. for regional/doubly-periodic (skips cube corners)
      logical :: sw_corner = .false.     !< tile owns SW cube corner
      logical :: se_corner = .false.     !< tile owns SE cube corner
      logical :: nw_corner = .false.     !< tile owns NW cube corner
      logical :: ne_corner = .false.     !< tile owns NE cube corner
   end type fv_grid_type

end module fv3_grid_types_mod
