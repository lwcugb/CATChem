!> \file TransportHalo_Mod.F90
!! \brief Halo (ghost-cell) exchange seam for the FV3 transport data domain.
!!
!! The FV3 flux operator works on a DATA domain `(isd:ied, jsd:jed)` that extends
!! the compute domain `(is:ie, js:je)` by `ng` ghost cells on every side. CATChem
!! carries its fields (tracers, winds, `DELP`) only on the compute domain, so
!! before each transport step the ghost ring must be filled from the neighbouring
!! data. This module is that seam: a single, well-defined place that fills the
!! halo of a data-domain 2-D field.
!!
!! ## Boundary conditions
!! Two per-axis conditions are supported today:
!!   - `HALO_BC_PERIODIC`  : wrap-around, for a globally periodic axis (the
!!                           longitude of a global lat/lon grid).
!!   - `HALO_BC_REPLICATE` : zero-gradient edge replication, for a regional
!!                           boundary or as a safe default.
!! Fields are cell-centred scalars from the halo's point of view; because the
!! only conditions here are periodic / replicate, the wind COMPONENTS `ua`,`va`
!! are filled with the same routine (no sign change is involved).
!!
!! ## Distributed / parallel backend
!! In a single-PET (serial / single-tile) build the ghost ring is filled in
!! process by the on-PET periodic / replicate closure below. For a decomposed
!! run the off-PET neighbour exchange is supplied by a backend that the host
!! registers with `set_halo_exchange_hook` (e.g. the MPI Cartesian backend in
!! `TransportHaloMPI_Mod`, or an ESMF `ESMF_FieldHalo` RouteHandle wrapper). A
!! registered backend owns the COMPLETE data-domain fill -- it exchanges the
!! interior neighbour rings across PETs AND applies the periodic / replicate
!! closure on any physical (non-neighbour) boundary -- so `halo_update` simply
!! delegates to it. When no backend is registered the on-PET closure runs, which
!! is exact for the serial / single-tile case. Callers never change.
!!
!! \note **Polar fold.** A global lat/lon grid whose poles lie on the y-boundary
!!       needs a cross-pole fold (longitude + 180, with a sign flip for vector
!!       components) rather than replication. That treatment is grid-orientation
!!       specific and is deferred; `HALO_BC_REPLICATE` on the y-axis is the
!!       current, safe placeholder near the poles.
!!
!! \note **Cubed-sphere vectors.** Across a cube edge the wind COMPONENTS rotate;
!!       a scalar exchange (this seam and the MPI backend) is correct for scalars
!!       and for lat-lon winds, but a cubed-sphere run must register a vector-
!!       aware backend (ESMF / FMS) for `ua`,`va`.
!!
!! \author CATChem Development Team
!! \version 0.1.0
module TransportHalo_Mod

   use error_mod, only: CC_SUCCESS, CC_FAILURE
   use fv3_grid_types_mod, only: fv_grid_bounds_type

   implicit none
   private

   public :: transport_halo_type
   public :: HALO_BC_REPLICATE, HALO_BC_PERIODIC
   public :: halo_fill_scalar
   public :: halo_update
   public :: halo_exchange_iface
   public :: set_halo_exchange_hook, clear_halo_exchange_hook, halo_exchange_registered
   public :: halo_global_reduce_iface
   public :: halo_global_max
   public :: set_halo_global_max_hook, clear_halo_global_max_hook
   public :: halo_global_sum
   public :: set_halo_global_sum_hook, clear_halo_global_sum_hook
   public :: halo_is_root_iface
   public :: halo_is_root
   public :: set_halo_is_root_hook, clear_halo_is_root_hook

   !> Zero-gradient edge replication (regional boundary / safe default).
   integer, parameter :: HALO_BC_REPLICATE = 0
   !> Wrap-around for a globally periodic axis (global-longitude x-halo).
   integer, parameter :: HALO_BC_PERIODIC  = 1

   !> \brief Off-PET halo-exchange backend signature.
   !!
   !! A registered backend fills the ENTIRE ghost ring of `arr` on the data
   !! domain `(isd:ied, jsd:jed)`: it exchanges the `ng`-wide neighbour rings
   !! across PETs and applies the `x_bc`/`y_bc` closure on any physical
   !! (non-neighbour) boundary. The compute domain `(is:ie, js:je)` holds valid
   !! data on entry.
   abstract interface
      subroutine halo_exchange_iface(arr, isd, ied, jsd, jed, is, ie, js, je, x_bc, y_bc, rc)
         integer, intent(in)    :: isd, ied, jsd, jed, is, ie, js, je, x_bc, y_bc
         real,    intent(inout) :: arr(isd:ied, jsd:jed)
         integer, intent(out)   :: rc
      end subroutine halo_exchange_iface
   end interface

   !> Registered off-PET exchange backend (null => serial on-PET closure).
   procedure(halo_exchange_iface), pointer, save :: off_pet_hook => null()

   !> \brief Off-PET scalar global-reduction backend signature.
   !!
   !! Reduces the scalar `val` in place across every PET of the decomposition
   !! (e.g. an MPI `Allreduce`). Needed so that quantities which must be
   !! globally consistent -- above all the Courant sub-cycling count
   !! `nsplt = int(1 + cmax)`, which drives how many times the paired halo
   !! exchange is called -- take the SAME value on every PET. Without it a PET
   !! with a locally larger Courant number would sub-cycle (and halo-exchange)
   !! more times than its neighbours and the paired exchange would deadlock.
   abstract interface
      subroutine halo_global_reduce_iface(val, rc)
         real,    intent(inout) :: val
         integer, intent(out)   :: rc
      end subroutine halo_global_reduce_iface
   end interface

   !> Registered global-max backend (null => serial identity, val unchanged).
   procedure(halo_global_reduce_iface), pointer, save :: global_max_hook => null()

   !> Registered global-sum backend (null => serial identity, val unchanged).
   !! Uses the same signature as the global-max hook; a registered backend does
   !! an MPI `Allreduce(MPI_SUM)`. Used only for reporting (e.g. the global
   !! tracer-mass conservation check in the transport debug output).
   procedure(halo_global_reduce_iface), pointer, save :: global_sum_hook => null()

   !> \brief Root-PET query signature (for report-once-on-root output).
   abstract interface
      logical function halo_is_root_iface()
      end function halo_is_root_iface
   end interface

   !> Registered root-query backend (null => serial, always root == .true.).
   procedure(halo_is_root_iface), pointer, save :: is_root_hook => null()

   !> \brief Halo policy for the transport data domain.
   type :: transport_halo_type
      integer :: x_bc = HALO_BC_REPLICATE  !< west/east (longitude) condition
      integer :: y_bc = HALO_BC_REPLICATE  !< south/north (latitude) condition
   end type transport_halo_type

contains

   !> \brief Fill the halo ring of a data-domain field (convenience wrapper).
   !!
   !! \param[inout] arr   field on the DATA domain `(bd%isd:bd%ied, bd%jsd:bd%jed)`;
   !!                     the compute domain `(is:ie, js:je)` must already be set.
   !! \param[in]    bd    FV3 index bounds.
   !! \param[in]    halo  halo policy (x/y boundary conditions).
   !! \param[out]   rc    CC_SUCCESS / CC_FAILURE.
   subroutine halo_update(arr, bd, halo, rc)
      type(fv_grid_bounds_type), intent(in)    :: bd
      real, intent(inout) :: arr(bd%isd:bd%ied, bd%jsd:bd%jed)
      type(transport_halo_type), intent(in)    :: halo
      integer,                   intent(out)   :: rc

      rc = CC_SUCCESS

      ! A registered distributed backend owns the complete fill (off-PET
      ! neighbour rings + physical-boundary closure); delegate to it. Otherwise
      ! run the serial on-PET periodic / replicate closure.
      if (associated(off_pet_hook)) then
         call off_pet_hook(arr, bd%isd, bd%ied, bd%jsd, bd%jed, &
                           bd%is, bd%ie, bd%js, bd%je, halo%x_bc, halo%y_bc, rc)
         return
      end if

      call halo_fill_scalar(arr, bd%is, bd%ie, bd%js, bd%je, &
                            bd%isd, bd%ied, bd%jsd, bd%jed, halo%x_bc, halo%y_bc)

   end subroutine halo_update

   !> \brief Register an off-PET halo-exchange backend (host / driver only).
   !!
   !! After registration every `halo_update` call delegates the full data-domain
   !! fill to `proc`. Pass a backend whose exchange is consistent with the grid
   !! decomposition (see `TransportHaloMPI_Mod`). Registering is a global,
   !! process-wide action; call `clear_halo_exchange_hook` to restore the serial
   !! on-PET closure.
   subroutine set_halo_exchange_hook(proc)
      procedure(halo_exchange_iface) :: proc
      off_pet_hook => proc
   end subroutine set_halo_exchange_hook

   !> \brief Remove any registered off-PET backend (restore serial closure).
   subroutine clear_halo_exchange_hook()
      off_pet_hook => null()
   end subroutine clear_halo_exchange_hook

   !> \brief Whether an off-PET halo-exchange backend is registered.
   logical function halo_exchange_registered()
      halo_exchange_registered = associated(off_pet_hook)
   end function halo_exchange_registered

   !> \brief Register an off-PET scalar global-max backend (host / driver only).
   !!
   !! Registered together with the exchange backend; call
   !! `clear_halo_global_max_hook` to restore the serial identity.
   subroutine set_halo_global_max_hook(proc)
      procedure(halo_global_reduce_iface) :: proc
      global_max_hook => proc
   end subroutine set_halo_global_max_hook

   !> \brief Remove any registered global-max backend (restore serial identity).
   subroutine clear_halo_global_max_hook()
      global_max_hook => null()
   end subroutine clear_halo_global_max_hook

   !> \brief Reduce a scalar to its global maximum across all PETs.
   !!
   !! Delegates to the registered backend (e.g. MPI `Allreduce(MPI_MAX)`); when
   !! no backend is registered (serial / single-PET) it is the identity, leaving
   !! `val` unchanged. This is a COLLECTIVE when a backend is registered, so it
   !! must be called the same number of times on every PET.
   subroutine halo_global_max(val, rc)
      real,    intent(inout) :: val
      integer, intent(out)   :: rc

      rc = CC_SUCCESS
      if (associated(global_max_hook)) call global_max_hook(val, rc)
   end subroutine halo_global_max

   !> \brief Register an off-PET scalar global-sum backend (host / driver only).
   subroutine set_halo_global_sum_hook(proc)
      procedure(halo_global_reduce_iface) :: proc
      global_sum_hook => proc
   end subroutine set_halo_global_sum_hook

   !> \brief Remove any registered global-sum backend (restore serial identity).
   subroutine clear_halo_global_sum_hook()
      global_sum_hook => null()
   end subroutine clear_halo_global_sum_hook

   !> \brief Reduce a scalar to its global sum across all PETs.
   !!
   !! Delegates to the registered backend (MPI `Allreduce(MPI_SUM)`); the serial
   !! identity leaves `val` unchanged. COLLECTIVE when a backend is registered,
   !! so it must be called the same number of times on every PET.
   subroutine halo_global_sum(val, rc)
      real,    intent(inout) :: val
      integer, intent(out)   :: rc

      rc = CC_SUCCESS
      if (associated(global_sum_hook)) call global_sum_hook(val, rc)
   end subroutine halo_global_sum

   !> \brief Register a root-PET query backend (host / driver only).
   subroutine set_halo_is_root_hook(proc)
      procedure(halo_is_root_iface) :: proc
      is_root_hook => proc
   end subroutine set_halo_is_root_hook

   !> \brief Remove any registered root-query backend (serial => always root).
   subroutine clear_halo_is_root_hook()
      is_root_hook => null()
   end subroutine clear_halo_is_root_hook

   !> \brief Whether this PET is the reduction root (for report-once output).
   !!
   !! Returns `.true.` in a serial / single-PET run (no backend registered), so
   !! report-once-on-root code prints exactly once in every configuration.
   logical function halo_is_root()
      halo_is_root = .true.
      if (associated(is_root_hook)) halo_is_root = is_root_hook()
   end function halo_is_root

   !> \brief Fill the halo ring of a data-domain field (explicit bounds).
   !!
   !! X halos are filled first over the interior rows, then Y halos over the full
   !! column range (including the freshly-filled x-halos) so the four corners are
   !! populated. Each axis uses `HALO_BC_PERIODIC` (wrap) or `HALO_BC_REPLICATE`
   !! (zero gradient).
   subroutine halo_fill_scalar(arr, is, ie, js, je, isd, ied, jsd, jed, x_bc, y_bc)
      integer, intent(in) :: is, ie, js, je, isd, ied, jsd, jed
      real, intent(inout) :: arr(isd:ied, jsd:jed)
      integer, intent(in) :: x_bc, y_bc

      integer :: i, j, h, nx, ny

      nx = ie - is + 1
      ny = je - js + 1

      ! --- West / east halos over the interior rows -------------------------
      if (x_bc == HALO_BC_PERIODIC) then
         do j = js, je
            do h = 1, is - isd
               arr(is-h, j) = arr(ie-h+1, j)   ! wrap from the east
            end do
            do h = 1, ied - ie
               arr(ie+h, j) = arr(is+h-1, j)   ! wrap from the west
            end do
         end do
      else   ! HALO_BC_REPLICATE
         do j = js, je
            do i = isd, is - 1
               arr(i,j) = arr(is,j)
            end do
            do i = ie + 1, ied
               arr(i,j) = arr(ie,j)
            end do
         end do
      end if

      ! --- South / north halos over all columns (incl. the x-halos) ---------
      if (y_bc == HALO_BC_PERIODIC) then
         do i = isd, ied
            do h = 1, js - jsd
               arr(i, js-h) = arr(i, je-h+1)
            end do
            do h = 1, jed - je
               arr(i, je+h) = arr(i, js+h-1)
            end do
         end do
      else   ! HALO_BC_REPLICATE
         do i = isd, ied
            do j = jsd, js - 1
               arr(i,j) = arr(i,js)
            end do
            do j = je + 1, jed
               arr(i,j) = arr(i,je)
            end do
         end do
      end if

   end subroutine halo_fill_scalar

end module TransportHalo_Mod
