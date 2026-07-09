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
!! ## Distributed / parallel backend (future)
!! In the current CATChem build each PET holds one local tile and the ESMF
!! `RouteHandle` used by the NUOPC cap is not reachable from a process, so the
!! fill is done in-process for the serial / single-tile case. When a
!! decomposition handle (ESMF VM + `ESMF_FieldHalo` RouteHandle, or an MPI
!! Cartesian communicator) becomes reachable from a process, it plugs in HERE:
!! `halo_update` becomes a wrapper that first does the off-PET exchange and then
!! applies the on-PET periodic / replicate closure below. Callers do not change.
!!
!! \note **Polar fold.** A global lat/lon grid whose poles lie on the y-boundary
!!       needs a cross-pole fold (longitude + 180, with a sign flip for vector
!!       components) rather than replication. That treatment is grid-orientation
!!       specific and is deferred; `HALO_BC_REPLICATE` on the y-axis is the
!!       current, safe placeholder near the poles.
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

   !> Zero-gradient edge replication (regional boundary / safe default).
   integer, parameter :: HALO_BC_REPLICATE = 0
   !> Wrap-around for a globally periodic axis (global-longitude x-halo).
   integer, parameter :: HALO_BC_PERIODIC  = 1

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
      call halo_fill_scalar(arr, bd%is, bd%ie, bd%js, bd%je, &
                            bd%isd, bd%ied, bd%jsd, bd%jed, halo%x_bc, halo%y_bc)

   end subroutine halo_update

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
