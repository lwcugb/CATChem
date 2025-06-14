!> \file gridstate_mod.F90
!! \brief Module for grid state variables
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2023
!!
!! This module contains the GridStateType derived type and related subroutines
!! for managing grid state information in the CATChem atmospheric chemistry model.
!! It handles grid dimensions, spatial configuration, and initialization.
!!
module GridState_Mod

   USE Error_Mod
   USE precision_mod

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: Grid_Init_State
   !> Derived type for grid state information
   !!
   !! This type contains all grid-related parameters including dimensions,
   !! vertical levels, soil layers, and horizontal area information.
   !!
   !! @param State Name identifier for this state object
   !! @param nx Number of grid points in x direction
   !! @param ny Number of grid points in y direction
   !! @param number_of_levels Number of vertical atmospheric levels
   !! @param number_of_soil_layers Number of soil layers
   !! @param area Grid cell horizontal area [m^2]
   type, public :: GridStateType
      CHARACTER(LEN=4) :: State = 'Grid'  !< Name of this state

      ! Integers
      integer :: nx = 1                        !< Number of grid points in x direction
      integer :: ny = 1                        !< Number of grid points in y direction
      integer :: number_of_levels              !< The number of vertical levels
      integer :: number_of_soil_layers         !< The number of soil layers

      ! Reals
      real(fp) :: area                         !< Grid cell horizontal area [m^2]

   end type GridStateType

contains

   !> Initialize a GridState object
   !!
   !! This subroutine initializes a GridState object with default values
   !! for grid dimensions, levels, and horizontal area.
   !!
   !! @param GridState The GridState object to be initialized
   !! @param RC The return code indicating success (CC_SUCCESS) or failure
   !!
   !! @note Currently sets default values for single-cell configuration
   subroutine Grid_Init_State(GridState, RC)
      use Error_Mod, only : CC_SUCCESS
      use Config_Opt_Mod, Only : ConfigType
      implicit none

      ! type(ConfigType),    intent(in)    :: Config     ! Input Options object
      type(GridStateType), intent(inout) :: GridState  ! Grid State object
      INTEGER,             INTENT(OUT)   :: RC         ! Success or failure

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
