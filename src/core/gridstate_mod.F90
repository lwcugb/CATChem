!> \file gridstate_mod.F90
!!
!! \brief Module for grid state variables
!!
!! This module contains subroutines and functions related to the grid state.
!!
!! \ingroup core_modules
!!!>
module GridState_Mod

   USE Error_Mod
   USE precision_mod

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: Grid_Init_State
   type, public :: GridStateType
      CHARACTER(LEN=4) :: State = 'Grid'  !< Name of this state

      ! Integers
      integer :: nx = 1
      integer :: ny = 1
      integer :: number_of_levels  !< The number of vertical levels
      integer :: number_of_soil_layers  !< The number of soil layers

      ! Reals
      real(fp) :: area  !< Grid cell horizontal area [m^2]

   end type GridStateType

contains

   !> \brief Initialize a GridState object
   !!
   !! This subroutine initializes a GridState object.
   !!
   !! \param Config The input config object.
   !! \param GridState The GridState object to be initialized.
   !! \param RC The return code.
   !!
   !! \ingroup core_modules
   !!!>
   subroutine Grid_Init_State(GridState, RC)
      use Error_Mod, only : CC_SUCCESS
      use Config_Opt_Mod, Only : ConfigType
      implicit none

      ! type(ConfigType),    intent(in)    :: Config     !< Input Options object
      type(GridStateType), intent(inout) :: GridState  !< Grid State object
      INTEGER,             INTENT(OUT)   :: RC         !< Success or failure

      ! Local variables
      CHARACTER(LEN=512) :: errMsg
      CHARACTER(LEN=255) :: thisLoc

      ! Set error handling defaults
      RC = CC_SUCCESS
      errMsg = ''
      thisLoc = 'Grid_Init_State() -> at initializing GridState'

      ! initialize GridState
      GridState%nx=1
      GridState%ny=1
      GridState%number_of_levels=1  ! FIXME: use Config?
      GridState%area = 1._fp

   end subroutine Grid_Init_State

end module GridState_Mod
