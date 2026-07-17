!> \file TransportHaloMPI_Mod.F90
!! \brief MPI Cartesian off-PET halo-exchange backend for the transport data domain.
!!
!! This is a concrete backend for the halo seam in \ref TransportHalo_Mod. On a
!! decomposed grid the compute domain `(is:ie, js:je)` of each PET is a patch of
!! a larger global grid; the FV3 flux operator needs an `ng`-wide ghost ring
!! filled from the neighbouring PETs. This module performs that fill with a
!! two-pass (x then y) MPI `Sendrecv` exchange over a Cartesian PET layout, then
!! applies a zero-gradient (replicate) closure on any physical boundary that has
!! no neighbour. The telescoping order (x over the interior rows first, then y
!! over the full width including the freshly-filled x-halos) reproduces exactly
!! the corner treatment of the serial \ref TransportHalo_Mod `halo_fill_scalar`,
!! so a single-PET run and an N-PET run fill the halo identically.
!!
!! ## Scope
!!   - Correct for scalars (`delp`, tracer `q`) on any Cartesian decomposition.
!!   - Correct for lat-lon winds (`ua`,`va`): longitude periodicity is carried by
!!     the wrap-around neighbour ranks; latitude physical boundaries replicate.
!!   - **NOT** sufficient for cubed-sphere winds: across a cube edge the vector
!!     COMPONENTS rotate, which a scalar exchange does not do. A cubed-sphere run
!!     must register a vector-aware backend (ESMF `ESMF_FieldHalo` / FMS
!!     `mpp_update_domains`) through the same seam instead of this module.
!!
!! ## Build
!! The backend is always compiled and MPI is always linked (MPI is a standard
!! CATChem build dependency). On a single-PET run the exchange reduces to the
!! same on-PET periodic / replicate closure as the serial path, so it is safe in
!! every configuration.
!!
!! ## Use
!! The host/driver, which owns the communicator and the decomposition, calls
!! `transport_halo_mpi_init(comm, npx, npy, px, py, x_periodic, y_periodic, rc)`
!! once after the grid is known; on success it registers this backend with
!! \ref TransportHalo_Mod. `npx`/`npy` are the number of PETs in x/y, `px`/`py`
!! are this PET's 0-based Cartesian coordinates (row-major rank = py*npx + px).
!!
!! \author CATChem Development Team
!! \version 0.1.0
module TransportHaloMPI_Mod

   use mpi
   use error_mod, only: CC_SUCCESS, CC_FAILURE
   use TransportHalo_Mod, only: set_halo_exchange_hook, clear_halo_exchange_hook, &
                                set_halo_global_max_hook, clear_halo_global_max_hook, &
                                HALO_BC_PERIODIC, HALO_BC_REPLICATE

   implicit none
   private

   public :: transport_halo_mpi_init
   public :: transport_halo_mpi_finalize
   public :: transport_halo_mpi_available

   ! --- Cartesian exchange state (set by init) ------------------------------
   logical, save :: is_initialized = .false.
   integer, save :: comm_c = MPI_COMM_NULL   !< communicator for the exchange
   integer, save :: nbr_w  = MPI_PROC_NULL   !< west  neighbour rank (or none)
   integer, save :: nbr_e  = MPI_PROC_NULL   !< east  neighbour rank (or none)
   integer, save :: nbr_s  = MPI_PROC_NULL   !< south neighbour rank (or none)
   integer, save :: nbr_n  = MPI_PROC_NULL   !< north neighbour rank (or none)

contains

   !> \brief Whether the MPI backend is compiled in and initialized.
   logical function transport_halo_mpi_available()
      transport_halo_mpi_available = is_initialized
   end function transport_halo_mpi_available

   !> \brief Initialize the Cartesian neighbour topology and register the backend.
   !!
   !! \param[in]  comm        MPI communicator spanning the decomposition.
   !! \param[in]  npx,npy     number of PETs along x (lon) and y (lat).
   !! \param[in]  px,py       this PET's 0-based Cartesian coordinates.
   !! \param[in]  x_periodic  wrap the x (longitude) neighbours (global grid).
   !! \param[in]  y_periodic  wrap the y (latitude) neighbours (rare; usually .false.).
   !! \param[out] rc          CC_SUCCESS / CC_FAILURE.
   subroutine transport_halo_mpi_init(comm, npx, npy, px, py, x_periodic, y_periodic, rc)
      integer, intent(in)  :: comm, npx, npy, px, py
      logical, intent(in)  :: x_periodic, y_periodic
      integer, intent(out) :: rc

      integer :: wx, ex, sy, ny_

      rc = CC_SUCCESS
      if (npx < 1 .or. npy < 1 .or. px < 0 .or. py < 0 .or. &
          px >= npx .or. py >= npy) then
         rc = CC_FAILURE
         return
      end if

      comm_c = comm

      ! West / east neighbours (x). Periodic wrap uses modular arithmetic;
      ! a non-periodic boundary has no neighbour (MPI_PROC_NULL -> replicate).
      if (x_periodic) then
         wx = modulo(px - 1, npx)
         ex = modulo(px + 1, npx)
         nbr_w = py * npx + wx
         nbr_e = py * npx + ex
      else
         if (px > 0)       then; nbr_w = py * npx + (px - 1); else; nbr_w = MPI_PROC_NULL; end if
         if (px < npx - 1) then; nbr_e = py * npx + (px + 1); else; nbr_e = MPI_PROC_NULL; end if
      end if

      ! South / north neighbours (y).
      if (y_periodic) then
         sy  = modulo(py - 1, npy)
         ny_ = modulo(py + 1, npy)
         nbr_s = sy  * npx + px
         nbr_n = ny_ * npx + px
      else
         if (py > 0)       then; nbr_s = (py - 1) * npx + px; else; nbr_s = MPI_PROC_NULL; end if
         if (py < npy - 1) then; nbr_n = (py + 1) * npx + px; else; nbr_n = MPI_PROC_NULL; end if
      end if

      call set_halo_exchange_hook(halo_mpi_exchange)
      call set_halo_global_max_hook(halo_mpi_global_max)
      is_initialized = .true.
   end subroutine transport_halo_mpi_init

   !> \brief Deregister the backend (restore the serial on-PET closure).
   subroutine transport_halo_mpi_finalize()
      call clear_halo_exchange_hook()
      call clear_halo_global_max_hook()
      is_initialized = .false.
      comm_c = MPI_COMM_NULL
      nbr_w = MPI_PROC_NULL; nbr_e = MPI_PROC_NULL
      nbr_s = MPI_PROC_NULL; nbr_n = MPI_PROC_NULL
   end subroutine transport_halo_mpi_finalize

   !> \brief Fill the full ghost ring of a data-domain field across PETs.
   !!
   !! Matches the \ref TransportHalo_Mod `halo_exchange_iface` signature. Two-pass
   !! Cartesian exchange: x over the interior rows, then y over the full column
   !! range (so the corners are populated), with a replicate closure on any
   !! physical (MPI_PROC_NULL) boundary. `x_bc`/`y_bc` select the physical-boundary
   !! closure; periodicity is already carried by the neighbour ranks.
   subroutine halo_mpi_exchange(arr, isd, ied, jsd, jed, is, ie, js, je, x_bc, y_bc, rc)
      integer, intent(in)    :: isd, ied, jsd, jed, is, ie, js, je, x_bc, y_bc
      real,    intent(inout) :: arr(isd:ied, jsd:jed)
      integer, intent(out)   :: rc

      integer :: ng, i, j, ierr
      integer :: nrow, ncol, cnt
      real, allocatable :: sbuf(:), rbuf(:)

      rc = CC_SUCCESS
      ierr = 0
      ng = is - isd
      if (ng < 1) return

      ! ---------------------------------------------------------------------
      ! Pass 1: west / east exchange over the interior rows (js:je).
      ! ---------------------------------------------------------------------
      nrow = je - js + 1
      cnt  = ng * nrow
      allocate(sbuf(cnt), rbuf(cnt))

      ! Send west interior (is:is+ng-1) to west; recv east halo (ie+1:ie+ng) from east.
      call pack_cols(arr, isd, ied, jsd, jed, is,        js, je, ng, sbuf)
      call sendrecv_bytes(sbuf, rbuf, cnt, nbr_w, nbr_e, ierr)
      if (ierr /= 0) then; rc = CC_FAILURE; go to 100; end if
      if (nbr_e /= MPI_PROC_NULL) call unpack_cols(arr, isd, ied, jsd, jed, ie+1, js, je, ng, rbuf)

      ! Send east interior (ie-ng+1:ie) to east; recv west halo (isd:is-1) from west.
      call pack_cols(arr, isd, ied, jsd, jed, ie-ng+1, js, je, ng, sbuf)
      call sendrecv_bytes(sbuf, rbuf, cnt, nbr_e, nbr_w, ierr)
      if (ierr /= 0) then; rc = CC_FAILURE; go to 100; end if
      if (nbr_w /= MPI_PROC_NULL) call unpack_cols(arr, isd, ied, jsd, jed, isd, js, je, ng, rbuf)

      ! Physical-boundary closure (no neighbour): replicate the edge column.
      if (nbr_w == MPI_PROC_NULL) then
         do j = js, je
            do i = isd, is - 1
               arr(i,j) = arr(is,j)
            end do
         end do
      end if
      if (nbr_e == MPI_PROC_NULL) then
         do j = js, je
            do i = ie + 1, ied
               arr(i,j) = arr(ie,j)
            end do
         end do
      end if
      deallocate(sbuf, rbuf)

      ! ---------------------------------------------------------------------
      ! Pass 2: south / north exchange over the FULL width (isd:ied), so the
      ! four corners (already carrying valid x-halo data) are propagated.
      ! ---------------------------------------------------------------------
      ncol = ied - isd + 1
      cnt  = ng * ncol
      allocate(sbuf(cnt), rbuf(cnt))

      ! Send south interior (js:js+ng-1) to south; recv north halo (je+1:je+ng) from north.
      call pack_rows(arr, isd, ied, jsd, jed, js,        ng, sbuf)
      call sendrecv_bytes(sbuf, rbuf, cnt, nbr_s, nbr_n, ierr)
      if (ierr /= 0) then; rc = CC_FAILURE; go to 100; end if
      if (nbr_n /= MPI_PROC_NULL) call unpack_rows(arr, isd, ied, jsd, jed, je+1, ng, rbuf)

      ! Send north interior (je-ng+1:je) to north; recv south halo (jsd:js-1) from south.
      call pack_rows(arr, isd, ied, jsd, jed, je-ng+1, ng, sbuf)
      call sendrecv_bytes(sbuf, rbuf, cnt, nbr_n, nbr_s, ierr)
      if (ierr /= 0) then; rc = CC_FAILURE; go to 100; end if
      if (nbr_s /= MPI_PROC_NULL) call unpack_rows(arr, isd, ied, jsd, jed, jsd, ng, rbuf)

      if (nbr_s == MPI_PROC_NULL) then
         do j = jsd, js - 1
            do i = isd, ied
               arr(i,j) = arr(i,js)
            end do
         end do
      end if
      if (nbr_n == MPI_PROC_NULL) then
         do j = je + 1, jed
            do i = isd, ied
               arr(i,j) = arr(i,je)
            end do
         end do
      end if

100   continue
      if (allocated(sbuf)) deallocate(sbuf)
      if (allocated(rbuf)) deallocate(rbuf)
      if (x_bc /= HALO_BC_PERIODIC .and. x_bc /= HALO_BC_REPLICATE) rc = CC_FAILURE
      if (y_bc /= HALO_BC_PERIODIC .and. y_bc /= HALO_BC_REPLICATE) rc = CC_FAILURE
   end subroutine halo_mpi_exchange

   !> \brief Pack `ng` columns starting at `i0` over rows [j0:j1] into a buffer.
   subroutine pack_cols(arr, isd, ied, jsd, jed, i0, j0, j1, ng, buf)
      integer, intent(in) :: isd, ied, jsd, jed, i0, j0, j1, ng
      real, intent(in)    :: arr(isd:ied, jsd:jed)
      real, intent(out)   :: buf(:)
      integer :: i, j, n
      n = 0
      do j = j0, j1
         do i = i0, i0 + ng - 1
            n = n + 1
            buf(n) = arr(i,j)
         end do
      end do
   end subroutine pack_cols

   !> \brief Unpack `ng` columns starting at `i0` over rows [j0:j1] from a buffer.
   subroutine unpack_cols(arr, isd, ied, jsd, jed, i0, j0, j1, ng, buf)
      integer, intent(in)  :: isd, ied, jsd, jed, i0, j0, j1, ng
      real, intent(inout)  :: arr(isd:ied, jsd:jed)
      real, intent(in)     :: buf(:)
      integer :: i, j, n
      n = 0
      do j = j0, j1
         do i = i0, i0 + ng - 1
            n = n + 1
            arr(i,j) = buf(n)
         end do
      end do
   end subroutine unpack_cols

   !> \brief Pack `ng` rows starting at `j0` over the full width into a buffer.
   subroutine pack_rows(arr, isd, ied, jsd, jed, j0, ng, buf)
      integer, intent(in) :: isd, ied, jsd, jed, j0, ng
      real, intent(in)    :: arr(isd:ied, jsd:jed)
      real, intent(out)   :: buf(:)
      integer :: i, j, n
      n = 0
      do j = j0, j0 + ng - 1
         do i = isd, ied
            n = n + 1
            buf(n) = arr(i,j)
         end do
      end do
   end subroutine pack_rows

   !> \brief Unpack `ng` rows starting at `j0` over the full width from a buffer.
   subroutine unpack_rows(arr, isd, ied, jsd, jed, j0, ng, buf)
      integer, intent(in)  :: isd, ied, jsd, jed, j0, ng
      real, intent(inout)  :: arr(isd:ied, jsd:jed)
      real, intent(in)     :: buf(:)
      integer :: i, j, n
      n = 0
      do j = j0, j0 + ng - 1
         do i = isd, ied
            n = n + 1
            arr(i,j) = buf(n)
         end do
      end do
   end subroutine unpack_rows

   !> \brief One paired send/recv using MPI_BYTE (kind-agnostic for the seam's
   !!        bare `real`, which the transport library promotes to r8).
   subroutine sendrecv_bytes(sbuf, rbuf, cnt, dest, src, ierr)
      real,    intent(in)  :: sbuf(:)
      real,    intent(out) :: rbuf(:)
      integer, intent(in)  :: cnt, dest, src
      integer, intent(out) :: ierr
      integer :: nbytes, status(MPI_STATUS_SIZE)

      ierr = 0
      if (cnt < 1) return
      nbytes = cnt * (storage_size(sbuf(1)) / 8)
      call MPI_Sendrecv(sbuf, nbytes, MPI_BYTE, dest, 1001, &
                        rbuf, nbytes, MPI_BYTE, src,  1001, &
                        comm_c, status, ierr)
   end subroutine sendrecv_bytes

   !> \brief Reduce a scalar to its global maximum across the exchange comm.
   !!
   !! Registered through the halo seam as the global-max backend. It guarantees
   !! that the Courant sub-cycling count `nsplt = int(1 + cmax)` is IDENTICAL on
   !! every PET, so the paired halo `Sendrecv` is called the same number of
   !! times everywhere (otherwise the exchange deadlocks). The MPI datatype is
   !! chosen from the storage size of the seam's bare `real` (r4 or r8 depending
   !! on the transport library build flags).
   subroutine halo_mpi_global_max(val, rc)
      real,    intent(inout) :: val
      integer, intent(out)   :: rc
      integer :: ierr, dtype

      rc = CC_SUCCESS
      if (comm_c == MPI_COMM_NULL) return
      if (storage_size(val) / 8 == 8) then
         dtype = MPI_REAL8
      else
         dtype = MPI_REAL4
      end if
      call MPI_Allreduce(MPI_IN_PLACE, val, 1, dtype, MPI_MAX, comm_c, ierr)
      if (ierr /= 0) rc = CC_FAILURE
   end subroutine halo_mpi_global_max

end module TransportHaloMPI_Mod
