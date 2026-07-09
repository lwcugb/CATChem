!> \file TransportProcessCreator_Mod.F90
!! \brief Factory for creating and registering the transport process.
!!
!! Follows the CATChem Process Factory pattern (cf. SettlingProcessCreator_Mod).
!!
!! \author CATChem Development Team
!! \version 0.1.0
module TransportProcessCreator_Mod

   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, CC_Error, CC_Warning, ErrorManagerType
   use ProcessInterface_Mod
   use ProcessTransportInterface_Mod

   implicit none
   private

   public :: create_transport_process
   public :: register_transport_process

contains

   !> Create a new transport process instance (uninitialized).
   !!
   !! @param[out] process  Allocated polymorphic process instance
   !! @param[out] rc       Return code
   subroutine create_transport_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessTransportInterface), allocatable :: transport_process
      integer :: alloc_stat

      rc = CC_SUCCESS

      allocate(transport_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         rc = CC_FAILURE
         return
      end if

      call move_alloc(transport_process, process)

   end subroutine create_transport_process

   !> Register the transport process with a ProcessManager.
   !!
   !! @param[inout] process_mgr The ProcessManager to register with
   !! @param[out]   rc          Return code
   subroutine register_transport_process(process_mgr, rc)
      use ProcessManager_Mod, only: ProcessManagerType

      type(ProcessManagerType), intent(inout) :: process_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      call process_mgr%register_process( &
         name='transport', &
         category='transport', &
         description='Horizontal tracer transport (FV3 flux-form / PPM kernel)', &
         creator=create_transport_process, &
         rc=rc &
         )

   end subroutine register_transport_process

end module TransportProcessCreator_Mod
